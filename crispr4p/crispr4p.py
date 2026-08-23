#!/usr/bin/env python3

import argparse
import os
import sys
import time
import multiprocessing
import queue

if __package__:
    from .cache import (
        cache_path,
        cache_exists,
        ensure_cache_dir,
        get_cached,
    )
    from .design import TableSorting, run_design
    from .coordinates import hit_coordinates, slice_bounds
    from .genome import GenomePamIndex
    from .guides import (
        find_guides,
        guide_primers,
        is_match,
        match_guide,
    )
    from .resources import (
        AnnotationParser,
        ReferenceResources,
        chromosomeFasta,
        read_fasta,
    )
    from .primers import (
        build_hr_dna,
        checking_primers,
        design_primers as _design_primers,
    )
else:  # Support direct script execution.
    from cache import (
        cache_path,
        cache_exists,
        ensure_cache_dir,
        get_cached,
    )
    from design import TableSorting, run_design
    from coordinates import hit_coordinates, slice_bounds
    from genome import GenomePamIndex
    from guides import (
        find_guides,
        guide_primers,
        is_match,
        match_guide,
    )
    from resources import (
        AnnotationParser,
        ReferenceResources,
        chromosomeFasta,
        read_fasta,
    )
    from primers import (
        build_hr_dna,
        checking_primers,
        design_primers as _design_primers,
    )

datapath = os.path.join(os.path.dirname(__file__), "../data/")

FASTA = datapath + 'Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa'
COORDINATES = datapath + 'COORDINATES.txt'
SYNONIMS = datapath + 'SYNONIMS.txt'
PRECOMPUTED = 'precomputed_stand_alone'
PRECOMPUTED_VERSION = 4
############### CONFIGURATION VALUES ###################
SEED_LENGTH = 20
UNIQUE_INDEX_LENGTH = (-12,-3)   # range of values selected for uniqueness


def timeit(method):

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print('%r (%r, %r) %2.2f sec' % \
              (method.__name__, args, kw, te-ts))
        return result

    return timed


class CPU_RAM:
    def getNumProccess(self):
        #return the number of process to run
        return 1
        # return multiprocessing.cpu_count()*3/4


class NGG(object):
    __slots__ = ('chromosome', 'pos', 'strand', 'seed', 'pam', 'primer')
    def __init__(self, chro, pos, strand, seed, pam):
        self.chromosome = chro
        self.pos = pos
        self.strand = strand
        self.seed = seed
        self.pam = pam
        self.primer = None


class PrimerDesign:
    '''
    Primer design for CRISPR.
    '''

    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    def __init__(
        self,
        sequenceFile,
        coordinates,
        synomins,
        verbose=False,
        precomputed_folder=PRECOMPUTED,
        regression=False,
        genome_index=None,
        reference_resources=None,
    ):
        self.sequenceFile_ = sequenceFile
        self.reference_resources = (
            reference_resources
            if reference_resources is not None
            else ReferenceResources.from_files(
                sequenceFile,
                coordinates,
                synomins,
            )
        )
        self.chromosomesData = self.reference_resources.chromosomes
        self._numAlternativeCheckings = 2
        self.annotationParser_ = self.reference_resources.annotations
        self.userNGGs = []
        self.tableNGGs = {}
        self.genome_index = genome_index
        self.NGGs = (
            genome_index.by_suffix
            if genome_index is not None
            else None
        )
        self.verbose = verbose
        self.precomputed_folder = precomputed_folder
        self.regression = regression
        if self.regression:
            self.getNGGsFromGenome()
        else:
            # todo
            self.NGGs = None


    def argumentParser(self):
        self.argp_ = argparse.ArgumentParser(description='crispr4p description')
        self.argp_.add_argument('--name', action='store', type=str, help='Name')
        self.argp_.add_argument('-cr','--chromosome', action='store', type=str, help='Chromosome')
        self.argp_.add_argument('-co','--coords', action='store', type=str, help='Coordinates')
        self.argp_.add_argument('--mismatch', action='store', type=int, default=0, help='Allowed amount of mismatches.')
        self.argp_.add_argument('--oligo', action='store', type=str, help='Oligo/sgRNA sequence (20bp seed or 23bp seed+PAM) to analyze.')

    def parseArgs(self, localArgs):
        if not hasattr(self, 'argp_'):
            self.argumentParser()
        self.argsList_ = self.argp_.parse_args(localArgs)

        if self.argsList_.oligo:
            return None
        elif self.argsList_.name:
            return self.annotationParser_.getCoordsFromName(self.argsList_.name)
        elif self.argsList_.coords and self.argsList_.chromosome:
            assert '...' in self.argsList_.coords, 'Coordinates need 3 dots in the middle'
            start, end = [x.strip() for x in self.argsList_.coords.split('...')]
            return self.argsList_.chromosome, start, end, '1'

        #print help and exit
        self.argp_.print_help()
        sys.exit()

    def checkCoords_(self, chromosome, start, end):
        '''
        Checks coordinates exists in this chromosome
            :param chromosome: string
            :param start: string
            :param end: string
            :return: Boolean
        '''
        crFasta = self.chromosomesData.get(chromosome, None)
        assert chromosome, 'Bad chromosome specified.'
        for x in (start, end):
            assert x.isdigit() and int(x)>0 and int(x)<=len(crFasta.sequence), \
            'Bad chromosomes specified'

        assert int(start) < int(end), 'Start "%s" must be smaller than end "%s".' % (start, end)

        return True

    def _getUserNGGs(self, crFasta, start, end):
        # Do not leave guides from an earlier query if discovery fails.
        self.userNGGs = []
        start_index, end_index = slice_bounds(start, end)
        self.userNGGs = find_guides(
            crFasta.sequence,
            crFasta.name,
            start_index,
            end_index - 1,
            hit_factory=NGG,
            reverse_complement=self.reverseComplement,
            seed_length=SEED_LENGTH,
        )

    def getPrimerGRNA(self, crFasta, start, end, ngg):
        start_index, end_index = slice_bounds(start, end)
        return guide_primers(
            crFasta.sequence,
            start_index,
            end_index - 1,
            ngg,
            self.reverseComplement,
        )

    def _genPrecomputedName(self, name, nMismatch, cr, start, end):
        ensure_cache_dir(self.precomputed_folder)
        systematic_name = None
        if name:
            #use systematic name (SPAC)
            systematic_name = [
                x for x in self.annotationParser_.synonims_ if name in x
            ][0][0]
        return cache_path(
            self.precomputed_folder,
            PRECOMPUTED_VERSION,
            nMismatch,
            cr,
            start,
            end,
            systematic_name=systematic_name,
        )

    @staticmethod
    def _isPrecomputed(precomputedName):
        if cache_exists(precomputedName):
            return True

    def run_(self, chromosome, start, end, nMismatch, name):
        '''
        Runs Primer design for CRISPR. giving a tuple
            :param coords: tuple(int, int, int)
            :return: tuple(1,2,3)
        '''
        return run_design(
            chromosome,
            start,
            end,
            nMismatch,
            regression=self.regression,
            chromosomes=self.chromosomesData,
            guide_matches=self.tableNGGs,
            ensure_genome_index=self.getNGGsFromGenome,
            discover_guides=self._getUserNGGs,
            match_guides=self.gRNA_Table,
            build_guide_primer=self.getPrimerGRNA,
            build_hr_dna=self.HR_DNA,
            build_checking_primers=self.CheckingPrimers,
        )

    def getNGGsFromGenome(self):
        '''
        Run at initialization.
        :return:
        '''
        if getattr(self, 'genome_index', None) is None:
            chromosome_sequences = {
                name: chromosome.sequence
                for name, chromosome in self.chromosomesData.items()
            }
            self.genome_index = GenomePamIndex.build(
                chromosome_sequences,
                hit_factory=NGG,
                seed_length=SEED_LENGTH,
            )

        self.NGGs = self.genome_index.by_suffix
        return self.genome_index

    def getOligoHitCoordinates(self, ngg):
        """
        Return an indexed oligo hit's reference PAM and Cas9 cut coordinates.

        NGG.pos is a 0-based position in whichever strand sequence was
        searched by getNGGsFromGenome(). Keep that internal value unchanged
        for matching, and normalize only when reporting a genomic hit.

        Returns:
            ((pam_start, pam_end), (cut_left, cut_right)), with all values
            expressed as 1-based reference-chromosome coordinates. The PAM
            interval is inclusive and the cut lies between cut_left/right.
        """
        chromosome = self.chromosomesData.get(ngg.chromosome)
        if chromosome is None:
            raise ValueError(
                'Chromosome "%s" not found for oligo hit' % ngg.chromosome
            )

        return hit_coordinates(
            ngg.pos,
            ngg.strand,
            ngg.pam,
            chromosome.sequence,
            self.reverseComplement,
        )

    @staticmethod
    def genomeCompare(g1, g2, nmismatch):
        return is_match(g1, g2, nmismatch)

    def _gRNA_Table_Worker(self, readDataQueue, storeDataQueue, nMismatch):
        '''

        :param readDataQueue:
        :param storeDataQueue:
        :return:
        '''
        while True:
            try:
                userNGG = readDataQueue.get_nowait()
            except queue.Empty:
                break
            storeDataQueue.put(self._single_table_worker(userNGG, nMismatch))

    def _single_table_worker(self, userNGG, nMismatch):
        return match_guide(userNGG, self.NGGs, nMismatch)

    def gRNA_Table(self, nMismatch):
        '''
        Match user ngg with genome nggs in parallel
            :param nMismatch: int
            :return: Tuple
        '''
        num_processes = CPU_RAM().getNumProccess()
        if num_processes > 1:

            #prepare data to read
            readData = multiprocessing.Queue()
            for n in self.userNGGs:
                readData.put(n)

            #queue to store data
            storeData = multiprocessing.Queue()

            #prepare parallel workers
            processList = []
            for w in range(num_processes):
                p = multiprocessing.Process(target=self._gRNA_Table_Worker, args=(readData, storeData, nMismatch,))
                p.start()
                processList.append(p)

            #collect data
            for x in range(len(self.userNGGs)):
                if self.verbose:
                    print('Generating NGG table:', x*100//len(self.userNGGs), '%')
                key, value = storeData.get()
                self.tableNGGs[key] = value

            #flush and close process
            readData.close()
            storeData.close()
            for p in processList:
                p.terminate()
            del processList

        else:

            for x in range(len(self.userNGGs)):
                if self.verbose:
                    print('Generating NGG table:', x * 100 // len(self.userNGGs), '%')
                key, value = self._single_table_worker(self.userNGGs[x], nMismatch)
                self.tableNGGs[key] = value


    def HR_DNA(self, crFasta, start, end):
        '''

            :param crFasta: string
            :param start: int
            :param end: int
            :return: Tuple
        '''
        start_index, end_index = slice_bounds(start, end)
        return build_hr_dna(
            crFasta.sequence,
            start_index,
            end_index,
            self.sequenceComplement_,
        )

    def CheckingPrimers(self, crFasta, start, end):
        '''

            :param crFasta: string
            :param start: int
            :param end: int
            :return: Tuple
        '''
        return self.CheckingPrimersWidth_(crFasta, start, end, 300)

    def CheckingPrimersWidth_(self, crFasta, start, end, width):
        start_index, end_index = slice_bounds(start, end)
        return checking_primers(
            crFasta.sequence,
            start_index,
            end_index,
            width,
            self._numAlternativeCheckings,
            primer_designer=_design_primers,
        )

    @staticmethod
    def sequenceComplement_(sequence):
        '''
        Returns the complement of an DNA sequence
            :param sequence: string
            :return: string
        '''
        return ''.join([PrimerDesign.complements[x] for x in sequence])

    @staticmethod
    def reverseComplement(sequence):
        return PrimerDesign.sequenceComplement_(sequence)[::-1]

    def run(self, chromosome, start, end, nMismatch, name):
        '''
        Runs Primer design for CRISPR.
            :param chromosome: string
            :param start: integer
            :param end: integer
            :param nMismatch: integer
            :param name: string
        '''
        self.checkCoords_(chromosome, start, end)
        precomputedName = self._genPrecomputedName(name, nMismatch, chromosome, start, end)
        return get_cached(
            precomputedName,
            lambda: self.run_(
                chromosome,
                int(start),
                int(end),
                nMismatch,
                name,
            ),
            self._isPrecomputed,
        )

    def runOligoQuery(self, oligo_seq, nMismatch):
        # 1. Clean input
        oligo_seq = oligo_seq.upper().strip()
        
        # We accept either 20bp or 23bp
        if len(oligo_seq) == 20:
            seed = oligo_seq
            pam = "NGG" # default
        elif len(oligo_seq) == 23:
            seed = oligo_seq[:20]
            pam = oligo_seq[20:]
        else:
            print(f"Error: Oligo sequence must be 20 bp (seed only) or 23 bp (seed + PAM). Received length: {len(oligo_seq)}")
            sys.exit(1)
            
        print(f"Querying S. pombe genome for oligo seed: {seed} with PAM: {pam} (Allowed mismatches: {nMismatch})")
        
        # Load NGGs from genome if not loaded
        if not self.NGGs:
            print("Loading genome data and indexing PAM sites...")
            self.getNGGsFromGenome()
            print("Genome indexed successfully.")
            
        query_ngg = NGG(chro='query', pos=0, strand=1, seed=seed, pam=pam)
        
        # Run search
        _, tableDict = self._single_table_worker(query_ngg, nMismatch)
        
        # Print summary of occurrences at different seed lengths
        print("\nSummary of genome occurrences matching the seed sequence:")
        print("-" * 65)
        print(f"{'Seed Length':12s} | {'Matching Sites (containing NGG or NAG PAM)':40s}")
        print("-" * 65)
        for length in (8, 10, 12, 14, 16, 18, 20):
            count = len(tableDict.get(length, []))
            print(f"{length:12d} | {count}")
        print("-" * 65)
        
        # If there are matching sites at 20bp, let's print their details!
        matches_20 = tableDict.get(20, [])
        if matches_20:
            print(f"\nDetails of {len(matches_20)} genomic target/off-target sites (full 20bp matches):")
            for idx, match in enumerate(matches_20):
                strand_str = "+" if match.strand == 1 else "-"
                pam_coords, cut_coords = self.getOligoHitCoordinates(match)
                print(
                    f"  {idx+1:2d}. Chromosome: {match.chromosome:4s} | "
                    f"PAM coordinates: {pam_coords[0]} - {pam_coords[1]} | "
                    f"Cut: {cut_coords[0]} | {cut_coords[1]} | "
                    f"Strand: {strand_str} | Sequence: {match.seed} | "
                    f"PAM: {match.pam}"
                )
        else:
            print("\nNo full 20bp matches found in the genome.")

    def runCL(self, localArgs):
        '''
        Run from Command line
            :param localArgs: list of strings
        '''
        # Keep the old entry point; normal CLI parsing lives in crispr4p.cli.
        self.argumentParser()
        self.argsList_ = self.argp_.parse_args(localArgs)
        
        if self.argsList_.oligo:
            self.runOligoQuery(self.argsList_.oligo, self.argsList_.mismatch)
            return

        parsed = self.parseArgs(localArgs)
        if parsed is None:
            return
        chromosome, start, end, strand = parsed

        #get primer and grna table
        name = self.annotationParser_.normalize_name(self.argsList_.name) if self.argsList_.name else None
        tablePos_grna, hr_dna, primercheck, gRNAs_match = self.run(chromosome, start, end, self.argsList_.mismatch, name)

        if not self.verbose:
            return

        for ind, elem in enumerate(tablePos_grna):

            #prints the position of the table and occurrences tuple
            print(ind+1, '-', elem[0], tablePos_grna[ind][2:])

            #prints grna report
            self.gRNA_report(elem[1])

        self.HR_DNA_report(hr_dna)
        if primercheck:
            self.CheckingPrimers_report(primercheck)


    def gRNA_report(self, gRNA, ):
        print('gRNA: ', gRNA[0], 'PAM: %d - %d' % gRNA[3], gRNA[5], gRNA[4])
        print('gRNAfw: ', gRNA[1])
        print('gRNArv: ', gRNA[2], '\n')

    def HR_DNA_report(self, hr_dna):
        print('HRfw: ', hr_dna[0])
        print('HRrv: ', hr_dna[1])
        print('Deleted DNA: ', hr_dna[2], '\n')

    def CheckingPrimers_report(self, primerDesigns):
        pm = primerDesigns[0]
        print('Check primer left: ', pm['PRIMER_LEFT_0_SEQUENCE'], 'TM:', pm['PRIMER_LEFT_0_TM'])
        print('Check primer right: ', pm['PRIMER_RIGHT_0_SEQUENCE'], 'TM:', pm['PRIMER_RIGHT_0_TM'])
        print('Deleted DNA product size: ', pm['PRIMER_PAIR_0_PRODUCT_SIZE'])
        print('Negative result product size: ', pm['negative_result'], '\n')

    def runWeb(self, name=None, cr=None, 
            start=None, end=None, strand=None, nMismatch=0):
        '''
        Function ready to be called from other sources
            :param name:
            :param cr:
            :param start:
            :param end:
            :param strand:
            :param nMismatch:
            :return:
        '''
        if name==None:
            if cr==None: raise ValueError('chromosome value (cr) must be given.')
            if start==None: raise ValueError('coordinate start index (start) must be given.')
            if end==None: raise ValueError('coordinate end index (end) must be given.')
            tablePos_grna, hr_dna, primercheck, gRNAs_match = self.run(cr, start, end, nMismatch, name)
        else:
            name = self.annotationParser_.normalize_name(name)
            cr, start, end, _ = self.annotationParser_.getCoordsFromName(name)
            tablePos_grna, hr_dna, primercheck, gRNAs_match = self.run(cr, start, end, nMismatch, name)

        return tablePos_grna, hr_dna, primercheck, name, cr, start, end

    def readsequence(self, sequenceFile):
        '''
        Returns a dictionary with header and data for the given file
            :param sequenceFile: string
            :return: dict
        '''
        # Direct callers expect a mutable dict.
        return dict(read_fasta(sequenceFile))


if __name__ == "__main__":
    # Support direct ``python crispr4p/crispr4p.py`` execution.
    if not __package__:
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if project_root not in sys.path:
            sys.path.insert(0, project_root)

    from crispr4p.cli import main

    raise SystemExit(main(sys.argv[1:]))
