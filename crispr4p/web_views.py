"""HTML views for CRISPR4P."""

import html
import json

from .annotations import VIABILITY_LABELS
from .disruption import STOP_CODONS, target_strand
from .models import DesignResult, OligoAnalysisResult
from .spedit import has_bsai, make_oligos


WEB_REGION_LABELS = {
    "CDS": "coding sequence (CDS)",
    "five_prime_UTR": "5' UTR",
    "three_prime_UTR": "3' UTR",
    "intron": "intron",
    "exon": "non-coding exon",
    "gene": "gene (no child feature at cut)",
}


def _region(feature_type):
    return WEB_REGION_LABELS.get(
        feature_type,
        feature_type.replace("_", " "),
    )


def _block(block):
    if block is None:
        return None
    return {
        "feature_type": block.feature_type,
        "label": _region(block.feature_type),
        "start": block.start,
        "end": block.end,
        "length": block.length,
    }


def _neighbor(neighbor):
    if neighbor is None:
        return None
    return {
        "block": _block(neighbor.block),
        "distance": neighbor.distance,
    }


def _context(context):
    within_block = context.relation == "within"
    if within_block:
        region = _region(context.block.feature_type)
    else:
        left = (
            _region(context.left.block.feature_type)
            if context.left is not None
            else "intergenic"
        )
        right = (
            _region(context.right.block.feature_type)
            if context.right is not None
            else "intergenic"
        )
        region = f"{left} / {right} boundary"

    cds = None
    if (
        within_block
        and context.block.feature_type == "CDS"
        and context.cds_position is not None
    ):
        cds = {
            "base": context.cds_position.base,
            "total": context.cds_position.total,
            "percent": round(context.cds_position.percent, 1),
        }

    cut_position = None
    if within_block:
        if context.strand == "+":
            five_prime_bases = context.lower_bases
            three_prime_bases = context.higher_bases
        else:
            five_prime_bases = context.higher_bases
            three_prime_bases = context.lower_bases
        cut_position = {
            "five_prime_bases": five_prime_bases,
            "three_prime_bases": three_prime_bases,
        }

    return {
        "transcript_id": context.transcript_id,
        "transcript_type": context.transcript_type.replace("_", " "),
        "transcript_strand": context.strand,
        "region": region,
        "crosses_boundary": not within_block,
        "cds": cds,
        "block": _block(context.block),
        "left_block": (
            _block(context.left.block)
            if context.left is not None
            else None
        ),
        "right_block": (
            _block(context.right.block)
            if context.right is not None
            else None
        ),
        "cut_position": cut_position,
        "upstream": _neighbor(context.upstream),
        "downstream": _neighbor(context.downstream),
    }


def _is_target(target_name, gene):
    if target_name is None:
        return False
    target = str(target_name).strip().casefold()
    identifiers = {gene.gene_id.casefold()}
    if gene.name:
        identifiers.add(gene.name.casefold())
    return target in identifiers


def _gene(gene, contexts, role, is_primary):
    display_name = (
        f"{gene.name} ({gene.gene_id})" if gene.name else gene.gene_id
    )
    return {
        "gene_id": gene.gene_id,
        "name": gene.name,
        "display_name": display_name,
        "gene_type": (
            gene.gene_type.replace("_", " ")
            if gene.gene_type
            else "not specified"
        ),
        "is_protein_coding": gene.is_protein_coding,
        "viability_status": gene.viability,
        "viability": VIABILITY_LABELS[gene.viability],
        "role": role,
        "is_primary": is_primary,
        "contexts": [_context(context) for context in contexts],
    }


def _nearby_gene(direction, nearby_gene):
    if nearby_gene is None:
        return None
    gene = nearby_gene.gene
    display_name = (
        f"{gene.name} ({gene.gene_id})" if gene.name else gene.gene_id
    )
    return {
        "direction": direction,
        "gene_id": gene.gene_id,
        "display_name": display_name,
        "gene_type": (
            gene.gene_type.replace("_", " ")
            if gene.gene_type
            else "not specified"
        ),
        "viability_status": gene.viability,
        "viability": VIABILITY_LABELS[gene.viability],
        "start": gene.start,
        "end": gene.end,
        "distance": nearby_gene.distance,
    }


def annotation_rows(guides, annotations, target_name=None):
    """Build browser data for each guide."""
    guides = tuple(guides)
    annotations = tuple(annotations)
    if len(guides) != len(annotations):
        raise ValueError("guide and annotation counts must match")

    candidates = []
    for guide, annotation in zip(guides, annotations):
        grouped_contexts = {}
        for context in annotation.transcripts:
            grouped_contexts.setdefault(context.gene.gene_id, []).append(context)

        target_gene = next(
            (
                contexts[0].gene
                for contexts in grouped_contexts.values()
                if _is_target(target_name, contexts[0].gene)
            ),
            None,
        )
        target_gene_id = target_gene.gene_id if target_gene else None

        genes = []
        for contexts in grouped_contexts.values():
            gene = contexts[0].gene
            is_primary = gene.gene_id == target_gene_id
            if is_primary:
                role = "Primary target"
            elif target_gene_id is not None:
                role = "Additional overlap"
            else:
                role = "Overlapping gene"
            genes.append(
                _gene(
                    gene,
                    contexts,
                    role,
                    is_primary,
                )
            )
        genes.sort(key=lambda gene: (not gene["is_primary"], gene["gene_id"]))

        candidates.append(
            {
                "chromosome": guide.chromosome,
                "seed": guide.seed,
                "gc_percent": round(guide.gc_percent, 1),
                "pam": guide.pam,
                "pam_coordinates": list(guide.pam_coordinates),
                "cut_coordinates": list(guide.cut_coordinates),
                "strand": "+" if guide.strand == 1 else "-",
                "coding_strand": target_strand(annotation, target_name),
                "coding_target": (
                    target_gene.is_protein_coding if target_gene else None
                ),
                "is_intergenic": annotation.is_intergenic,
                "gene_count": len(annotation.genes),
                "genes": genes,
                "nearest_genes": [
                    nearby
                    for nearby in (
                        _nearby_gene("Lower-coordinate", annotation.lower_gene),
                        _nearby_gene("Higher-coordinate", annotation.higher_gene),
                    )
                    if nearby is not None
                ],
            }
        )
    return candidates


def spedit_rows(table_pos_grna) -> list[dict]:
    """Build SpEDIT data for each guide row."""
    candidates = []

    for row in table_pos_grna:
        try:
            # row[1][0] is the 20-nt guide in the old result layout.
            guide = row[1][0]
            forward, reverse = make_oligos(guide)

            candidates.append(
                {
                    "guide": guide,
                    "forward": forward,
                    "reverse": reverse,
                    "has_internal_bsai": has_bsai(guide),
                    "error": None,
                }
            )

        except (IndexError, TypeError, ValueError) as error:
            candidates.append(
                {
                    "guide": "",
                    "forward": "",
                    "reverse": "",
                    "has_internal_bsai": False,
                    "error": str(error),
                }
            )

    return candidates


def cassette_data(cassette_choices):
    """Build a shared cassette catalog and guide-specific choice lists."""
    catalog = {}
    choices_by_guide = []
    for choices in cassette_choices:
        ids = []
        for cassette in choices:
            forward, reverse = make_oligos(cassette.guide)
            row = {
                "id": cassette.id,
                "sequence": cassette.sequence,
                "length": len(cassette.sequence),
                "guide": cassette.guide,
                "pam": cassette.pam,
                "gc_percent": round(cassette.gc_percent, 1),
                "spedit_forward": forward,
                "spedit_reverse": reverse,
                "has_internal_bsai": has_bsai(cassette.guide),
                "frames": [
                    " ".join(
                        f"{codon}*" if codon in STOP_CODONS else codon
                        for codon in frame
                    )
                    for frame in cassette.frames
                ],
            }
            previous = catalog.setdefault(str(cassette.id), row)
            if previous != row:
                raise ValueError("cassette IDs must identify one sequence")
            ids.append(cassette.id)
        choices_by_guide.append(ids)
    return {"catalog": catalog, "choices": choices_by_guide}


def donor_data(donor_choices):
    """Build one shared-arm record per guide."""
    rows = []
    for choices in donor_choices:
        if not choices:
            rows.append(None)
            continue

        first = choices[0]
        oligos = first.oligos
        row = {
            "coding_strand": first.coding_strand,
            "arm_length": first.arm_length,
            "left_arm": first.left_arm,
            "right_arm": first.right_arm,
            "hrfw_length": len(oligos.forward),
            "hrrv_length": len(oligos.reverse),
            "overlap_length": len(oligos.overlap),
            "hr_product_length": oligos.product_length,
        }
        for donor in choices[1:]:
            donor_oligos = donor.oligos
            if (
                donor.coding_strand != first.coding_strand
                or donor.arm_length != first.arm_length
                or donor.left_arm != first.left_arm
                or donor.right_arm != first.right_arm
                or len(donor_oligos.forward) != row["hrfw_length"]
                or len(donor_oligos.reverse) != row["hrrv_length"]
                or len(donor_oligos.overlap) != row["overlap_length"]
                or donor_oligos.product_length != row["hr_product_length"]
            ):
                raise ValueError(
                    "donors for one guide must share arms and oligo lengths"
                )
        rows.append(row)
    return rows


def restoration_data(donors):
    """Build guide-aligned restoration donor data."""
    rows = []
    for donor in donors:
        if donor is None:
            rows.append(None)
            continue
        oligos = donor.oligos
        rows.append({
            "arm_length": donor.arm_length,
            "left_arm": donor.left_arm,
            "right_arm": donor.right_arm,
            "hrfw_length": len(oligos.forward),
            "hrrv_length": len(oligos.reverse),
            "overlap_length": len(oligos.overlap),
            "hr_product_length": oligos.product_length,
        })
    return rows


def _compact_json(value):
    return json.dumps(value, separators=(",", ":"))


def render_query_error():
    return (
        '<font color="red"><h3>Error: Please fill either Name, '
        'Coordinates, or Oligo Sequence</h3></font>'
    )


def render_error(error):
    return f'<font color="red"><h3>ERROR during execution: {str(error)}</h3></font>'


def render_gene_error(query):
    safe_query = html.escape(str(query), quote=True)
    return (
        '<font color="red"><h3>Gene &quot;'
        f"{safe_query}"
        '&quot; was not found. Check the spelling or try a PomBase '
        "systematic gene ID, such as SPBPB2B2.01.</h3></font>"
    )


def render_ambiguous_gene_error(error):
    message = html.escape(str(error), quote=True)
    return f'<font color="red"><h3>{message}</h3></font>'


def render_oligo_error(sequence_length):
    return (
        '<font color="red"><h3>Error: Oligo sequence must be 20 bp '
        f'(seed only) or 23 bp (seed + PAM). Current length: {sequence_length}'
        '</h3></font>'
    )


def render_oligo(result: OligoAnalysisResult):
    """Render an oligo result."""
    if result.has_internal_bsai:
        spedit_warning = (
            '<strong style="color: #b00020;">'
            "Warning: this guide contains an internal BsaI recognition site."
            "</strong>"
        )
    else:
        spedit_warning = "No internal BsaI site detected."

    details_html = ""
    if result.full_matches:
        details_html += '<table border="1" cellpadding="5" cellspacing="0" style="border-collapse: collapse; width: 100%; font-family: monospace; font-size: 12px;">'
        details_html += '<tr style="background-color: #D1F0A6;"><th>#</th><th>Chromosome</th><th>PAM coordinates (1-based, inclusive)</th><th>Cas9 cut</th><th>Strand</th><th>Genomic Target Sequence (Seed)</th><th>PAM</th></tr>'
        for index, match in enumerate(result.full_matches):
            strand = "+" if match.strand == 1 else "-"
            details_html += (
                f'<tr><td>{index+1}</td><td>{match.chromosome}</td>'
                f'<td>{match.pam_coordinates[0]} - {match.pam_coordinates[1]}</td>'
                f'<td>{match.cut_coordinates[0]} | {match.cut_coordinates[1]}</td>'
                f'<td>{strand}</td><td>{match.seed}</td>'
                f'<td>{match.pam}</td></tr>'
            )
        details_html += '</table>'
    else:
        details_html = "<p>No full 20bp target/off-target matches found in the genome.</p>"

    block = f"""
        <div id="search_content">
          <div id="search_summary">
              <h4>Oligo Search Results:</h4>
              <b>Oligo Sequence (Query)</b>: {result.oligo_sequence}<br>
              <b>Seed Segment (20bp)</b>: {result.seed}<br>
              <b>Mismatches Allowed</b>: {result.n_mismatch}<br>
              <hr>

              <h4>SpEDIT/pLSB BsaI Golden Gate oligos</h4>

              <b>Forward oligo, 52 nt, 5&#8242;&rarr;3&#8242;</b>:
              <code>{result.spedit_forward}</code><br>

              <b>Reverse oligo, 52 nt, 5&#8242;&rarr;3&#8242;</b>:
              <code>{result.spedit_reverse}</code><br>

              <b>Internal BsaI site check</b>: {spedit_warning}<br>
          </div>

          <h3 class="toggle_header">Genome Match Summary</h3>
          <div style="padding: 10px;">
            <table border="1" cellpadding="5" cellspacing="0" style="border-collapse: collapse; width: 50%;">
              <thead>
                <tr style="background-color: #D1F0A6;">
                  <th>Seed Match Length</th>
                  <th>Matching Sites (adjacent to NGG/NAG PAM)</th>
                </tr>
              </thead>
              <tbody>
                <tr><td>8 bp</td><td>{result.match_counts.get(8, 0)}</td></tr>
                <tr><td>10 bp</td><td>{result.match_counts.get(10, 0)}</td></tr>
                <tr><td>12 bp</td><td>{result.match_counts.get(12, 0)}</td></tr>
                <tr><td>14 bp</td><td>{result.match_counts.get(14, 0)}</td></tr>
                <tr><td>16 bp</td><td>{result.match_counts.get(16, 0)}</td></tr>
                <tr><td>18 bp</td><td>{result.match_counts.get(18, 0)}</td></tr>
                <tr><td>20 bp</td><td>{result.match_counts.get(20, 0)}</td></tr>
              </tbody>
            </table>
          </div>

          <h3 class="toggle_header">Details of Full 20bp Matches</h3>
          <div style="padding: 10px;">
            {details_html}
          </div>
        </div>
        """
    return block


def render_design(
    result: DesignResult,
    guide_annotations,
    template_text,
    cassette_choices=(),
    disruption_donors=(),
    restoration_donors=(),
):
    """Render a design result with the HTML template."""
    primer = result.checking_primers[0] if result.checking_primers else {}

    def tm(value):
        try:
            return "%d &deg;C" % int(round(float(value)))
        except (ValueError, TypeError):
            return "- &deg;C"

    def primer_sequence(value):
        return f"5'-{value}-3'" if value else "-"

    context = {
        'name': result.name or '-',
        'chromosome': result.chromosome,
        'start': result.start,
        'end': result.end,
        'hrfw': result.hr_dna[0],
        'hrrv': result.hr_dna[1],
        'deleted_dna': result.hr_dna[2],
        'primer_left': primer_sequence(
            primer.get('PRIMER_LEFT_0_SEQUENCE')
        ),
        'left_tm': tm(primer.get('PRIMER_LEFT_0_TM')),
        'primer_right': primer_sequence(
            primer.get('PRIMER_RIGHT_0_SEQUENCE')
        ),
        'right_tm': tm(primer.get('PRIMER_RIGHT_0_TM')),
        'deleted_dna_size': str(
            primer.get('PRIMER_PAIR_0_PRODUCT_SIZE', '-')
        ) + " (bp)",
        'negative_result_size': str(
            primer.get('negative_result', '-')
        ) + " (bp)",
    }
    context['json_table'] = _compact_json(result.guide_table)
    context['spedit_json'] = _compact_json(
        spedit_rows(result.guide_table)
    )
    context['annotation_json'] = _compact_json(
        annotation_rows(
            result.guides,
            guide_annotations,
            target_name=result.name,
        )
    )
    context['cassette_json'] = _compact_json(cassette_data(cassette_choices))
    context['donor_json'] = _compact_json(donor_data(disruption_donors))
    context['restoration_json'] = _compact_json(
        restoration_data(restoration_donors)
    )

    return template_text % context
