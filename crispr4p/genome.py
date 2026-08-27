"""Genome-wide PAM index."""

import re
from collections.abc import Mapping
from types import MappingProxyType


class GenomePamIndex(Mapping):
    """Read-only suffix lookup for PAM-adjacent hits."""

    PAM_SUFFIXES = ("GG", "AG")
    COMPLEMENTS = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
    }

    def __init__(self, buckets):
        frozen_buckets = {
            suffix: tuple(hits)
            for suffix, hits in buckets.items()
        }
        self._buckets = MappingProxyType(frozen_buckets)
        self.hit_count = sum(len(hits) for hits in frozen_buckets.values())

    @classmethod
    def build(cls, chromosome_sequences, hit_factory, seed_length=20):
        """Index NGG/NAG sites on both strands."""
        buckets = {}
        for chromosome, forward_sequence in chromosome_sequences.items():
            strands = (
                (1, forward_sequence),
                (-1, cls.reverse_complement(forward_sequence)),
            )
            for strand, sequence in strands:
                for pam_suffix in cls.PAM_SUFFIXES:
                    # Lookahead keeps both overlapping PAMs in GGG.
                    for match in re.finditer(f"(?={pam_suffix})", sequence):
                        position = match.start()
                        if position < seed_length + 1:
                            continue
                        pam = sequence[position - 1:position + 2]
                        seed = sequence[
                            position - seed_length - 1:position - 1
                        ]
                        suffix = seed[-8:]
                        buckets.setdefault(suffix, []).append(
                            hit_factory(
                                chromosome,
                                position,
                                strand,
                                seed,
                                pam,
                            )
                        )

        return cls(buckets)

    @classmethod
    def reverse_complement(cls, sequence):
        return "".join(cls.COMPLEMENTS[base] for base in sequence)[::-1]

    @property
    def by_suffix(self):
        """Return hits grouped by 8-nt suffix."""
        return self._buckets

    def __getitem__(self, suffix):
        return self._buckets[suffix]

    def __iter__(self):
        return iter(self._buckets)

    def __len__(self):
        return len(self._buckets)
