"""
Coordinate liftover utilities for handling hg19/hg38 build mismatches.
"""

from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger(__name__)

_BUILD_ALIASES: dict[str, str] = {
    "hg19": "hg19",
    "grch37": "hg19",
    "b37": "hg19",
    "37": "hg19",
    "hg38": "hg38",
    "grch38": "hg38",
    "b38": "hg38",
    "38": "hg38",
}


def normalise_build(build: str) -> str:
    key = build.lower().strip()
    canonical = _BUILD_ALIASES.get(key)
    if canonical is None:
        raise ValueError(
            f"Unknown genome build: {build!r}. "
            "Use hg19, hg38 (or GRCh37/GRCh38, b37/b38, 37/38)."
        )
    return canonical


def builds_differ(a: str | None, b: str | None) -> bool:
    """Return True when both builds are specified and they are different."""
    if a is None or b is None:
        return False
    return normalise_build(a) != normalise_build(b)


def make_lifted_lookup(
    variants: np.ndarray,
    from_build: str,
    to_build: str,
) -> dict[str, list[tuple[str, str, int]]]:
    """
    Lift variant positions from *from_build* to *to_build* and return a
    positional lookup mapping ``'chrom:pos'`` keys (in *to_build* coordinates)
    to ``[(a1, a2, row_idx), ...]`` entries.

    Both ``chr``-prefixed and bare chromosome forms are inserted so the lookup
    works regardless of the VCF's CHROM convention.

    Variants whose position cannot be lifted are omitted (NaN for any trait
    whose VCF uses *to_build*).

    Requires ``pyliftover``: ``pip install pyliftover``.
    """
    try:
        from pyliftover import LiftOver  # type: ignore
    except ImportError:
        raise ImportError(
            "pyliftover is required for cross-build matching: pip install pyliftover"
        )

    from_build = normalise_build(from_build)
    to_build = normalise_build(to_build)

    log.info("Lifting %d variants %s → %s", len(variants), from_build, to_build)
    lo = LiftOver(from_build, to_build)

    lookup: dict[str, list[tuple[str, str, int]]] = {}
    n_fail = 0

    for i, row in enumerate(variants):
        chrom = str(row["chrom"])
        pos = int(row["pos"])
        a1 = str(row["a1"])
        a2 = str(row["a2"])

        if not chrom or pos == 0:
            n_fail += 1
            continue

        chrom_in = chrom if chrom.startswith("chr") else f"chr{chrom}"
        result = lo.convert_coordinate(chrom_in, pos - 1)

        if not result:
            n_fail += 1
            continue

        new_chrom_full: str = result[0][0]
        new_pos = int(result[0][1]) + 1
        new_chrom_bare = new_chrom_full.replace("chr", "")

        for chrom_form in (new_chrom_bare, new_chrom_full):
            key = f"{chrom_form}:{new_pos}"
            lookup.setdefault(key, []).append((a1, a2, i))

    if n_fail:
        log.warning(
            "Liftover %s→%s: %d/%d variants failed to map",
            from_build, to_build, n_fail, len(variants),
        )

    return lookup
