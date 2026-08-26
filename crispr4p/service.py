"""CRISPR4P application service."""

import os
from dataclasses import replace
from pathlib import Path

from .annotations import GenomeAnnotations
from .crispr4p import NGG, PrimerDesign
from .diagnostics import CassetteOption, pcr_rflp
from .disruption import (
    CASSETTE_FORMATS,
    build_donor,
    load_cassettes,
    recut_sites,
    target_gene,
    target_strand,
    valid_rescue_site,
)
from .guides import match_guide
from .models import DesignResult, OligoAnalysisResult, OligoMatch
from .primers import (
    PrimerNotFoundError,
    insertion_checks as design_insertion_checks,
    insertion_primers as design_insertion_primers,
)
from .resources import GeneNameNotFoundError
from .restoration import build_donor as build_restoration_donor
from .spedit import has_bsai, make_oligos


PROJECT_DATA_DIRECTORY = Path(__file__).resolve().parent.parent / "data"


class OligoLengthError(ValueError):
    """Raised for an invalid oligo length."""

    def __init__(self, sequence_length):
        self.sequence_length = sequence_length
        super().__init__(
            "Oligo sequence must be 20 bp (seed only) or 23 bp "
            f"(seed + PAM). Received length: {sequence_length}"
        )


def parse_oligo(oligo_sequence):
    """Parse a 20-nt seed or 23-nt seed and PAM."""
    normalized_sequence = oligo_sequence.upper().strip()
    if len(normalized_sequence) == 20:
        return normalized_sequence, normalized_sequence, "NGG"
    if len(normalized_sequence) == 23:
        return normalized_sequence, normalized_sequence[:20], normalized_sequence[20:]
    raise OligoLengthError(len(normalized_sequence))


class Crispr4pService:
    """Application entry point for CRISPR4P operations."""

    def __init__(
        self,
        sequence_file,
        coordinates_file,
        synonyms_file,
        precomputed_folder="precomputed",
        designer_factory=PrimerDesign,
        designer_verbose=False,
        genome_index=None,
        reference_resources=None,
        annotation_gff_file=None,
        gene_viability_file=None,
        genome_annotations=None,
        annotation_factory=GenomeAnnotations.from_files,
        cassette_file=None,
        cassettes=None,
    ):
        self.sequence_file = os.fspath(sequence_file)
        self.coordinates_file = os.fspath(coordinates_file)
        self.synonyms_file = os.fspath(synonyms_file)
        self.precomputed_folder = os.fspath(precomputed_folder)
        self._designer_factory = designer_factory
        self.designer_verbose = designer_verbose
        self._genome_index = genome_index
        self._reference_resources = reference_resources
        self.annotation_gff_file = (
            os.fspath(annotation_gff_file)
            if annotation_gff_file is not None
            else None
        )
        self.gene_viability_file = (
            os.fspath(gene_viability_file)
            if gene_viability_file is not None
            else None
        )
        self._genome_annotations = genome_annotations
        self._annotation_factory = annotation_factory
        self.cassette_file = (
            os.fspath(cassette_file) if cassette_file is not None else None
        )
        self._cassettes = tuple(cassettes) if cassettes is not None else None

    @classmethod
    def from_project_data(
        cls,
        precomputed_folder="precomputed",
        designer_factory=PrimerDesign,
        designer_verbose=False,
        genome_index=None,
        reference_resources=None,
        genome_annotations=None,
        annotation_factory=GenomeAnnotations.from_files,
        cassettes=None,
    ):
        """Use the reference files included with CRISPR4P."""
        return cls(
            PROJECT_DATA_DIRECTORY
            / "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa",
            PROJECT_DATA_DIRECTORY / "COORDINATES.txt",
            PROJECT_DATA_DIRECTORY / "SYNONIMS.txt",
            precomputed_folder=precomputed_folder,
            designer_factory=designer_factory,
            designer_verbose=designer_verbose,
            genome_index=genome_index,
            reference_resources=reference_resources,
            annotation_gff_file=(
                PROJECT_DATA_DIRECTORY
                / "Schizosaccharomyces_pombe_all_chromosomes.gff3"
            ),
            gene_viability_file=(
                PROJECT_DATA_DIRECTORY / "gene_viability.tsv"
            ),
            genome_annotations=genome_annotations,
            annotation_factory=annotation_factory,
            cassette_file=PROJECT_DATA_DIRECTORY / "stop_cassettes.json",
            cassettes=cassettes,
        )

    def design_gene(self, name, n_mismatch=0):
        """Design for a gene name or systematic ID."""
        try:
            return self._run(name=name, n_mismatch=n_mismatch)
        except GeneNameNotFoundError as legacy_error:
            if (
                self._genome_annotations is None
                and self.annotation_gff_file is None
            ):
                raise

            gene = self._load_annotations().find_gene(name)
            if gene is None:
                raise legacy_error from None

            try:
                result = self._run(
                    name=gene.gene_id,
                    n_mismatch=n_mismatch,
                )
            except GeneNameNotFoundError:
                if gene.chromosome is None:
                    raise legacy_error from None
                result = self._run(
                    chromosome=gene.chromosome,
                    start=str(gene.start),
                    end=str(gene.end),
                    n_mismatch=n_mismatch,
                )

            return self._with_name(
                result,
                gene.name or gene.gene_id,
            )

    def design_region(
        self,
        chromosome,
        start,
        end,
        strand=None,
        n_mismatch=0,
    ):
        """Design for a chromosome interval."""
        return self._run(
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand,
            n_mismatch=n_mismatch,
        )

    def analyze_oligo(self, oligo_sequence, n_mismatch=0):
        """Search the genome for a guide sequence."""
        normalized_sequence, seed, pam = parse_oligo(oligo_sequence)

        spedit_forward, spedit_reverse = make_oligos(seed)
        designer = self._designer()
        designer.getNGGsFromGenome()
        self._keep_index(designer)

        query = NGG(
            chro="query",
            pos=0,
            strand=1,
            seed=seed,
            pam=pam,
        )
        _, legacy_matches = match_guide(
            query,
            designer.NGGs,
            n_mismatch,
        )

        match_counts = {
            length: len(legacy_matches.get(length, []))
            for length in (8, 10, 12, 14, 16, 18, 20)
        }
        full_matches = []
        for match in legacy_matches.get(20, []):
            pam_coordinates, cut_coordinates = (
                designer.getOligoHitCoordinates(match)
            )
            full_matches.append(
                OligoMatch(
                    chromosome=match.chromosome,
                    pam_coordinates=pam_coordinates,
                    cut_coordinates=cut_coordinates,
                    strand=match.strand,
                    seed=match.seed,
                    pam=match.pam,
                )
            )

        return OligoAnalysisResult(
            oligo_sequence=normalized_sequence,
            seed=seed,
            pam=pam,
            n_mismatch=n_mismatch,
            spedit_forward=spedit_forward,
            spedit_reverse=spedit_reverse,
            has_internal_bsai=has_bsai(seed),
            match_counts=match_counts,
            full_matches=full_matches,
        )

    def annotate_cut(self, chromosome, cut_coordinates):
        """Annotate a Cas9 cut boundary."""
        annotations = self._load_annotations()
        return annotations.annotate_cut(chromosome, cut_coordinates)

    def annotate_guide(self, guide):
        """Annotate one guide."""
        return self.annotate_cut(guide.chromosome, guide.cut_coordinates)

    def annotate_guides(self, guides):
        """Annotate guides in display order."""
        return tuple(self.annotate_guide(guide) for guide in guides)

    @property
    def cassettes(self):
        """Return the packaged stop cassettes."""
        if self._cassettes is None:
            if self.cassette_file is None:
                return ()
            self._cassettes = load_cassettes(self.cassette_file)
        return self._cassettes

    def cassette_choices(self, guides, annotations, target_name=None):
        """Return cassettes without a similar guide target at the junction."""
        guides = tuple(guides)
        annotations = tuple(annotations)
        if len(guides) != len(annotations):
            raise ValueError("guide and annotation counts must match")

        resources = self._load_resources()
        choices = []
        for guide, annotation in zip(guides, annotations):
            gene = target_gene(annotation, target_name)
            if gene is not None and not gene.is_protein_coding:
                choices.append(())
                continue

            chromosome = resources.chromosomes[guide.chromosome]
            coding_strand = target_strand(annotation, target_name)
            choices.append(
                tuple(
                    cassette
                    for cassette in self.cassettes
                    if not recut_sites(
                        chromosome.sequence,
                        guide.cut_coordinates,
                        guide.seed,
                        cassette,
                        coding_strand,
                    )
                )
            )
        return tuple(choices)

    def disruption_donors(
        self,
        guides,
        annotations,
        cassettes,
        arm_length,
        target_name=None,
    ):
        """Build donors for aligned groups of safe cassette choices."""
        guides = tuple(guides)
        annotations = tuple(annotations)
        cassettes = tuple(cassettes)
        if not len(guides) == len(annotations) == len(cassettes):
            raise ValueError("guide, annotation, and cassette counts must match")

        resources = self._load_resources()
        donors = []
        for guide, annotation, choices in zip(
            guides,
            annotations,
            cassettes,
        ):
            coding_strand = target_strand(annotation, target_name)
            if coding_strand is None:
                donors.append(())
                continue

            reference = resources.chromosomes[guide.chromosome].sequence
            donors.append(
                tuple(
                    build_donor(
                        reference,
                        guide.cut_coordinates,
                        cassette,
                        coding_strand,
                        arm_length,
                    )
                    for cassette in choices
                )
            )
        return tuple(donors)

    def insertion_primers(
        self,
        chromosome,
        cut_coordinates,
        arm_length=80,
        insert_length=23,
        window=300,
    ):
        """Design checking primers for one insertion site."""
        resources = self._load_resources()
        chromosome = str(chromosome).strip()
        try:
            reference = resources.chromosomes[chromosome].sequence
        except KeyError:
            raise ValueError(f"unknown chromosome: {chromosome}") from None

        return design_insertion_primers(
            reference,
            cut_coordinates,
            arm_length=arm_length,
            window=window,
            insert_length=insert_length,
        )

    def insertion_checks(
        self,
        chromosome,
        cut_coordinates,
        cassette_id,
        coding_strand,
        arm_length=80,
        window=300,
    ):
        """Design edit-spanning and junction-checking primer pairs."""
        resources = self._load_resources()
        chromosome = str(chromosome).strip()
        try:
            reference = resources.chromosomes[chromosome].sequence
        except KeyError:
            raise ValueError(f"unknown chromosome: {chromosome}") from None

        cassette = next(
            (item for item in self.cassettes if item.id == cassette_id),
            None,
        )
        if cassette is None:
            raise ValueError(f"unknown cassette: {cassette_id}")
        if coding_strand not in ("+", "-"):
            raise ValueError("coding strand must be + or -")

        return design_insertion_checks(
            reference,
            cut_coordinates,
            cassette.orient(coding_strand),
            arm_length=arm_length,
            window=window,
        )

    def cassette_options(
        self,
        chromosome,
        cut_coordinates,
        guide,
        cassette_id,
        coding_strand,
        arm_length=80,
        window=300,
    ):
        """Evaluate the selectable cassette formats at one cut site."""
        resources = self._load_resources()
        chromosome = str(chromosome).strip()
        try:
            reference = resources.chromosomes[chromosome].sequence
        except KeyError:
            raise ValueError(f"unknown chromosome: {chromosome}") from None

        cassette = next(
            (item for item in self.cassettes if item.id == cassette_id),
            None,
        )
        if cassette is None:
            raise ValueError(f"unknown cassette: {cassette_id}")
        if coding_strand not in ("+", "-"):
            raise ValueError("coding strand must be + or -")

        try:
            base_pair = design_insertion_primers(
                reference,
                cut_coordinates,
                arm_length=arm_length,
                window=window,
                insert_length=len(cassette.sequence),
            )
        except PrimerNotFoundError:
            base_pair = None

        core = cassette.orient(coding_strand)
        options = []
        for item in CASSETTE_FORMATS:
            insert = item.orient(cassette, coding_strand)
            recut = recut_sites(
                reference,
                cut_coordinates,
                guide,
                cassette,
                coding_strand,
                cassette_format=item,
            )
            rescue_ok = valid_rescue_site(
                reference,
                cut_coordinates,
                cassette,
                coding_strand,
                cassette_format=item,
            )
            pair = (
                replace(base_pair, insert_length=len(insert))
                if base_pair is not None
                else None
            )
            checks = None
            if base_pair is not None:
                try:
                    checks = design_insertion_checks(
                        reference,
                        cut_coordinates,
                        insert,
                        arm_length=arm_length,
                        window=window,
                        core=core,
                        spanning=base_pair,
                    )
                except PrimerNotFoundError:
                    pass

            digest = None
            if item.enzyme is not None and pair is not None:
                digest = pcr_rflp(
                    reference,
                    cut_coordinates,
                    insert,
                    pair,
                    item,
                    window=window,
                )

            available = not recut and rescue_ok
            if item.enzyme is not None:
                available = (
                    available
                    and digest is not None
                    and digest.available
                )
            options.append(
                CassetteOption(
                    cassette_format=item,
                    coding_sequence=item.sequence(cassette),
                    insert=insert,
                    available=available,
                    spanning=pair,
                    checks=checks,
                    digest=digest,
                )
            )
        return tuple(options)

    def restoration_donors(
        self,
        guides,
        annotations,
        arm_length,
        target_name=None,
    ):
        guides = tuple(guides)
        annotations = tuple(annotations)
        if len(guides) != len(annotations):
            raise ValueError("guide and annotation counts must match")

        resources = self._load_resources()
        donors = []
        for guide, annotation in zip(guides, annotations):
            if target_strand(annotation, target_name) is None:
                donors.append(None)
                continue

            reference = resources.chromosomes[guide.chromosome].sequence
            donors.append(
                build_restoration_donor(
                    reference,
                    guide.cut_coordinates,
                    arm_length,
                )
            )
        return tuple(donors)

    def _run(
        self,
        name=None,
        chromosome=None,
        start=None,
        end=None,
        strand=None,
        n_mismatch=0,
    ):
        designer = self._designer()
        result = DesignResult.from_legacy(
            designer.runWeb(
                name=name,
                cr=chromosome,
                start=start,
                end=end,
                strand=strand,
                nMismatch=n_mismatch,
            )
        )
        self._keep_index(designer)
        return result

    def _designer(self):
        # PrimerDesign stores mutable state, so each request gets an instance.
        options = {
            "precomputed_folder": self.precomputed_folder,
            "verbose": self.designer_verbose,
        }
        if self._genome_index is not None:
            options["genome_index"] = self._genome_index
        if self._reference_resources is not None:
            options["reference_resources"] = self._reference_resources

        designer = self._designer_factory(
            self.sequence_file,
            self.coordinates_file,
            self.synonyms_file,
            **options,
        )
        self._keep_resources(designer)
        return designer

    @staticmethod
    def _with_name(result, name):
        if result.name == name:
            return result
        legacy_result = result.to_legacy()
        return DesignResult.from_legacy(
            legacy_result[:3] + (name,) + legacy_result[4:]
        )

    @property
    def genome_index(self):
        return self._genome_index

    @property
    def reference_resources(self):
        return self._reference_resources

    @property
    def genome_annotations(self):
        """Return the annotation index after it has been loaded."""
        return self._genome_annotations

    def _load_annotations(self):
        if self._genome_annotations is not None:
            return self._genome_annotations
        if self.annotation_gff_file is None:
            raise RuntimeError("Cut-site annotation GFF3 is not configured")
        self._genome_annotations = self._annotation_factory(
            self.annotation_gff_file,
            self.gene_viability_file,
        )
        return self._genome_annotations

    def _load_resources(self):
        if self._reference_resources is None:
            self._designer()
        return self._reference_resources

    def _keep_index(self, designer):
        genome_index = getattr(designer, "genome_index", None)
        if genome_index is not None:
            self._genome_index = genome_index

    def _keep_resources(self, designer):
        resources = getattr(designer, "reference_resources", None)
        if resources is not None:
            self._reference_resources = resources
