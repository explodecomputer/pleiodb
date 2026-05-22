"""
ALID (Allele-ID) utilities.

ALID format: ``CHROM:POS_A1_A2``
  - A1 and A2 are in canonical (alphabetically sorted) order: A1 ≤ A2.
  - Effect allele = A2 (alphabetically second).
  - Alleles longer than ALLELE_MAX_LEN characters are compressed to
    ``{allele[:8]}~{sha256(allele)[:4]}`` (13 chars, tilde is the
    truncation marker).  Collision resistance is only required within a
    single genomic position; 4 hex chars (16 bits) is sufficient.

Public API
----------
compress_allele(allele) -> str
is_compressed(allele) -> bool
canonical_alid(chrom, pos, a1, a2) -> (alid_str, was_flipped)
"""

from __future__ import annotations

import hashlib

ALLELE_MAX_LEN = 20   # alleles with length > this are compressed
_COMPRESSED_LEN = 13  # 8 + 1 (tilde) + 4


def compress_allele(allele: str) -> str:
    """Return a compact representation for alleles longer than ALLELE_MAX_LEN.

    Alleles of 20 characters or fewer are returned unchanged.
    Longer alleles become ``{allele[:8]}~{sha256(allele)[:4]}``.
    """
    if len(allele) <= ALLELE_MAX_LEN:
        return allele
    h = hashlib.sha256(allele.encode()).hexdigest()[:4]
    return f"{allele[:8]}~{h}"


def is_compressed(allele: str) -> bool:
    """Return True when *allele* is a compressed ALID allele (contains '~')."""
    return "~" in allele


def canonical_alid(
    chrom: str,
    pos: int,
    a1: str,
    a2: str,
) -> tuple[str, bool]:
    """Build a canonical ALID string from its components.

    Parameters
    ----------
    chrom, pos : genomic coordinates
    a1, a2     : raw allele strings (order does not matter)

    Returns
    -------
    alid_str   : ``CHROM:POS_CA1_CA2`` where CA1 ≤ CA2 (after compression)
    was_flipped: True when the input ``a1 / a2`` order was swapped to
                 enforce canonical ordering.  Callers should negate EAF
                 when this flag is True.
    """
    ca1 = compress_allele(a1)
    ca2 = compress_allele(a2)

    if ca1 <= ca2:
        return f"{chrom}:{pos}_{ca1}_{ca2}", False
    else:
        return f"{chrom}:{pos}_{ca2}_{ca1}", True
