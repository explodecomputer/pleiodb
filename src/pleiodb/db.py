"""
GWASDatabase — top-level handle that owns all sub-matrices and metadata.

Directory layout:
  {name}.pleiodb/
    meta.json          — dimensions, chunk sizes, dtypes, thresholds, var_y per trait
    variants.tsv       — V rows + header: alid, eaf  (replaces variants.npy + eaf.f16)
    traits.tsv         — T rows + header: trait_id, trait_name, N, K, neff_study,
                         var_y, n_variants, n_variants_var_y, n_sig_*
                         (replaces traits.npy + neff_base.f32)
    zscore.bin/.cidx   — V×T  int16  (z * 100, NA = -32768)
    neff.bin/.cidx     — V×T  uint16 (log2-encoded, NA = 0xFFFF)
    lambda.bin/.cidx   — T×T  float16 (symmetric)
    masks/
      {thr}.coo.zst    — COO-format significance hits: sorted (v_idx u4, t_idx u4) pairs
"""

from __future__ import annotations

import json
import struct
from pathlib import Path
from typing import Sequence

import numpy as np
import zstandard as zstd

from .alid import parse_alid
from .store import ChunkedMatrix
from .quantize import (
    decode_z, decode_neff, decode_eaf,
    encode_z, encode_neff, encode_eaf,
    reconstruct_beta_se, Z_SCALE, Z_NA, NEFF_NA,
)

_DCTX = zstd.ZstdDecompressor()
_CCTX = zstd.ZstdCompressor(level=3, threads=-1)

_TRAITS_DT = np.dtype([("id", "U64"), ("name", "U256")])

_VARIANTS_DT = np.dtype([
    ("id",    "U64"),
    ("chrom", "U10"),
    ("pos",   np.uint32),
    ("a1",    "U64"),
    ("a2",    "U64"),
])


def _load_variants_tsv(
    path: "Path",
) -> tuple[np.ndarray, np.ndarray]:
    """Parse ``variants.tsv`` → (variants structured array, eaf float32 array).

    Shared by the ``variants`` and ``eaf`` properties; both are populated in a
    single pass so the file is only read once.
    """
    rows = []
    eaf_list = []
    with open(path) as fh:
        header_seen = False
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if not header_seen:
                # First non-blank line is the header
                header_seen = True
                continue
            parts = line.split("\t", 1)
            alid = parts[0]
            raw_eaf = parts[1].strip() if len(parts) > 1 else ""
            eaf_val = float(raw_eaf) if raw_eaf else np.nan

            chrom, pos, a1, a2 = parse_alid(alid)
            rows.append((alid, chrom, pos, a1, a2))
            eaf_list.append(eaf_val)

    variants = np.array(rows, dtype=_VARIANTS_DT)
    eaf = np.array(eaf_list, dtype=np.float32)
    return variants, eaf


def _load_traits_tsv(
    path: "Path",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse ``traits.tsv`` → (traits structured array, neff_study, var_y).

    Shared by the ``traits``, ``neff_base``, and ``var_y`` properties; all
    three are populated in a single pass so the file is only read once.

    Returns
    -------
    traits      : structured array with fields id(U64), name(U256)
    neff_study  : float32(T) — study-level Neff per trait (used as Neff fallback)
    var_y       : float32(T) — per-trait phenotypic variance (NaN if unknown)
    """
    rows: list[tuple[str, str]] = []
    neff_list: list[float] = []
    var_y_list: list[float] = []
    with open(path) as fh:
        header: list[str] | None = None
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if header is None:
                header = line.split("\t")
                continue
            cells = dict(zip(header, line.split("\t")))
            rows.append((cells.get("trait_id", ""), cells.get("trait_name", "")))
            raw_neff = cells.get("neff_study", "").strip()
            neff_list.append(float(raw_neff) if raw_neff else np.nan)
            raw_vy = cells.get("var_y", "").strip()
            var_y_list.append(float(raw_vy) if raw_vy else np.nan)

    traits = np.array(rows, dtype=_TRAITS_DT)
    neff_study = np.array(neff_list, dtype=np.float32)
    var_y = np.array(var_y_list, dtype=np.float32)
    return traits, neff_study, var_y


class GWASDatabase:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._meta: dict = {}
        self._variants: np.ndarray | None = None
        self._traits: np.ndarray | None = None
        self._neff_base_arr: np.ndarray | None = None   # loaded together with _traits
        self._var_y_arr: np.ndarray | None = None        # loaded together with _traits
        self._zscore: ChunkedMatrix | None = None
        self._neff: ChunkedMatrix | None = None
        self._eaf: np.ndarray | None = None
        self._lambda: ChunkedMatrix | None = None
        self._rho: ChunkedMatrix | None = None

    # ------------------------------------------------------------------
    # Open / close
    # ------------------------------------------------------------------

    @classmethod
    def open(cls, path: str | Path) -> "GWASDatabase":
        db = cls(path)
        meta_path = db.path / "meta.json"
        if not meta_path.exists():
            raise FileNotFoundError(f"Not a pleiodb directory: {path}")
        db._meta = json.loads(meta_path.read_text())
        V, T = db._meta["V"], db._meta["T"]
        cs = db._meta.get("chunk_shape", [512, 512])
        db._zscore = ChunkedMatrix(db.path / "zscore", (V, T), np.int16, tuple(cs))
        db._neff = ChunkedMatrix(db.path / "neff", (V, T), np.uint16, tuple(cs))
        return db

    def close(self) -> None:
        pass  # no persistent file handles in read path

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    # ------------------------------------------------------------------
    # Lazy-loaded accessors
    # ------------------------------------------------------------------

    @property
    def V(self) -> int:
        return self._meta["V"]

    @property
    def T(self) -> int:
        return self._meta["T"]

    def _ensure_variants_loaded(self) -> None:
        """Load variants.tsv once, populating both _variants and _eaf caches."""
        if self._variants is None:
            self._variants, self._eaf = _load_variants_tsv(
                self.path / "variants.tsv"
            )

    @property
    def variants(self) -> np.ndarray:
        self._ensure_variants_loaded()
        return self._variants

    def _ensure_traits_loaded(self) -> None:
        """Load traits.tsv once, populating _traits, _neff_base_arr, _var_y_arr."""
        if self._traits is None:
            self._traits, self._neff_base_arr, self._var_y_arr = _load_traits_tsv(
                self.path / "traits.tsv"
            )

    @property
    def traits(self) -> np.ndarray:
        self._ensure_traits_loaded()
        return self._traits

    @property
    def eaf(self) -> np.ndarray:
        self._ensure_variants_loaded()
        return self._eaf

    @property
    def neff_base(self) -> np.ndarray:
        self._ensure_traits_loaded()
        return self._neff_base_arr

    @property
    def var_y(self) -> np.ndarray:
        """Per-trait phenotypic variance from traits.tsv (float32, NaN if unknown)."""
        self._ensure_traits_loaded()
        return self._var_y_arr

    @property
    def lambda_matrix(self) -> ChunkedMatrix:
        if self._lambda is None:
            T = self.T
            cs = self._meta.get("lambda_chunk_shape", [512, 512])
            self._lambda = ChunkedMatrix(self.path / "lambda", (T, T), np.float16, tuple(cs))
        return self._lambda

    @property
    def rho_matrix(self) -> ChunkedMatrix:
        """T×T sample-overlap-weighted phenotypic correlation matrix (float16)."""
        if self._rho is None:
            T = self.T
            cs = self._meta.get("rho_chunk_shape", [512, 512])
            self._rho = ChunkedMatrix(self.path / "rho", (T, T), np.float16, tuple(cs))
        return self._rho

    def get_rho_block(
        self, t_indices_row: Sequence[int], t_indices_col: Sequence[int]
    ) -> np.ndarray:
        """Return a sub-block of the rho matrix as float32."""
        tr = np.asarray(t_indices_row, np.int64)
        tc = np.asarray(t_indices_col, np.int64)
        t_min = int(min(tr.min(), tc.min()))
        t_max = int(max(tr.max(), tc.max())) + 1
        raw = self.rho_matrix.get_block(t_min, t_max, t_min, t_max)
        block = raw.astype(np.float32)
        return block[np.ix_(tr - t_min, tc - t_min)]

    # ------------------------------------------------------------------
    # Variant / trait index resolution
    # ------------------------------------------------------------------

    def variant_index(self, ids: Sequence[str]) -> np.ndarray:
        """Return integer indices for variant IDs (raises KeyError if any missing)."""
        id_arr = self.variants["id"]
        lookup = {v: i for i, v in enumerate(id_arr)}
        return np.array([lookup[x] for x in ids], dtype=np.int64)

    def trait_index(self, ids: Sequence[str]) -> np.ndarray:
        id_arr = self.traits["id"]
        lookup = {v: i for i, v in enumerate(id_arr)}
        return np.array([lookup[x] for x in ids], dtype=np.int64)

    def region_to_variant_indices(self, chrom: str, start: int, end: int) -> np.ndarray:
        """Return variant indices whose position falls within [start, end]."""
        v = self.variants
        mask = (v["chrom"] == chrom) & (v["pos"] >= start) & (v["pos"] <= end)
        return np.where(mask)[0].astype(np.int64)

    # ------------------------------------------------------------------
    # Z-score queries
    # ------------------------------------------------------------------

    def zscore_variant(self, v_idx: int) -> np.ndarray:
        """All T z-scores for one variant → float32(T)."""
        raw = self._zscore.get_block(v_idx, v_idx + 1, 0, self.T)[0]
        return decode_z(raw)

    def zscore_trait(self, t_idx: int) -> np.ndarray:
        """All V z-scores for one trait → float32(V)."""
        raw = self._zscore.get_block(0, self.V, t_idx, t_idx + 1)[:, 0]
        return decode_z(raw)

    def zscore_region(self, chrom: str, start: int, end: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Return (v_indices, z_matrix) for all variants in genomic region.
        z_matrix shape: (n_variants, T), float32.
        """
        v_idx = self.region_to_variant_indices(chrom, start, end)
        if len(v_idx) == 0:
            return v_idx, np.empty((0, self.T), dtype=np.float32)
        v_min, v_max = int(v_idx.min()), int(v_idx.max()) + 1
        raw = self._zscore.get_block(v_min, v_max, 0, self.T)
        z = decode_z(raw)[v_idx - v_min]
        return v_idx, z

    def zscore_block(
        self, v_indices: Sequence[int], t_indices: Sequence[int]
    ) -> np.ndarray:
        """Arbitrary V×T block (indices need not be contiguous) → float32."""
        v_idx = np.asarray(v_indices, np.int64)
        t_idx = np.asarray(t_indices, np.int64)
        v_min, v_max = int(v_idx.min()), int(v_idx.max()) + 1
        t_min, t_max = int(t_idx.min()), int(t_idx.max()) + 1
        raw = self._zscore.get_block(v_min, v_max, t_min, t_max)
        return decode_z(raw)[np.ix_(v_idx - v_min, t_idx - t_min)]

    # ------------------------------------------------------------------
    # Beta / SE reconstruction
    # ------------------------------------------------------------------

    def beta_se_variant(self, v_idx: int) -> tuple[np.ndarray, np.ndarray]:
        z = self.zscore_variant(v_idx)
        neff_raw = self._neff.get_block(v_idx, v_idx + 1, 0, self.T)[0]
        neff = decode_neff(neff_raw)
        eaf_v = self.eaf[v_idx]
        return reconstruct_beta_se(z, neff, np.array([eaf_v], dtype=np.float32))

    def beta_se_block(
        self, v_indices: Sequence[int], t_indices: Sequence[int]
    ) -> tuple[np.ndarray, np.ndarray]:
        v_idx = np.asarray(v_indices, np.int64)
        t_idx = np.asarray(t_indices, np.int64)
        v_min, v_max = int(v_idx.min()), int(v_idx.max()) + 1
        t_min, t_max = int(t_idx.min()), int(t_idx.max()) + 1
        z_raw = self._zscore.get_block(v_min, v_max, t_min, t_max)
        neff_raw = self._neff.get_block(v_min, v_max, t_min, t_max)
        z = decode_z(z_raw)[np.ix_(v_idx - v_min, t_idx - t_min)]
        neff = decode_neff(neff_raw)[np.ix_(v_idx - v_min, t_idx - t_min)]
        eaf = self.eaf[v_idx]
        return reconstruct_beta_se(z, neff, eaf)

    # ------------------------------------------------------------------
    # P-value / significance queries
    # ------------------------------------------------------------------

    def _z_to_pval_threshold(self, pval: float) -> float:
        from scipy.stats import norm  # type: ignore
        return abs(norm.ppf(pval / 2))

    def query_significant(
        self, pval: float = 5e-8, mask_name: str | None = None
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Return (v_indices, t_indices, z_values) for all cells passing pval.

        Tries the pre-computed COO mask first; falls back to scanning the
        full z-score matrix if no mask is available.
        """
        if mask_name is None:
            mask_name = f"{pval:.0e}"

        mask_path = self.path / "masks" / f"{mask_name}.coo.zst"
        if mask_path.exists():
            return self._load_coo_mask(mask_path, pval)

        # Full-scan fallback (slow for large databases)
        return self._scan_significant(pval)

    def _load_coo_mask(
        self, path: Path, pval: float
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        blob = _DCTX.decompress(path.read_bytes())
        pairs = np.frombuffer(blob, dtype=np.uint32).reshape(-1, 2)
        v_idx = pairs[:, 0].astype(np.int64)
        t_idx = pairs[:, 1].astype(np.int64)
        z_threshold = self._z_to_pval_threshold(pval)
        z = self.zscore_block(v_idx, t_idx)
        z_diag = np.array([z[i, i] for i in range(len(v_idx))])
        return v_idx, t_idx, z_diag

    def _scan_significant(
        self, pval: float
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        z_thresh = self._z_to_pval_threshold(pval)
        vs, ts, zs = [], [], []
        chunk_v = self._zscore.CV
        chunk_t = self._zscore.CT
        for vi in range(self._zscore.n_v_chunks):
            for ti in range(self._zscore.n_t_chunks):
                raw = self._zscore.get_raw_chunk(vi, ti)
                z = decode_z(raw)
                hit_v, hit_t = np.where(np.abs(z) >= z_thresh)
                vs.append(hit_v + vi * chunk_v)
                ts.append(hit_t + ti * chunk_t)
                zs.append(z[hit_v, hit_t])
        return (
            np.concatenate(vs).astype(np.int64),
            np.concatenate(ts).astype(np.int64),
            np.concatenate(zs).astype(np.float32),
        )

    # ------------------------------------------------------------------
    # Lambda (sample overlap) access
    # ------------------------------------------------------------------

    def get_lambda_block(
        self, t_indices_row: Sequence[int], t_indices_col: Sequence[int]
    ) -> np.ndarray:
        tr = np.asarray(t_indices_row, np.int64)
        tc = np.asarray(t_indices_col, np.int64)
        t_min = int(min(tr.min(), tc.min()))
        t_max = int(max(tr.max(), tc.max())) + 1
        raw = self.lambda_matrix.get_block(t_min, t_max, t_min, t_max)
        block = raw.astype(np.float32)
        return block[np.ix_(tr - t_min, tc - t_min)]

    # ------------------------------------------------------------------
    # Utility
    # ------------------------------------------------------------------

    def info(self) -> dict:
        m = dict(self._meta)
        try:
            zbin = self.path / "zscore.bin"
            nbin = self.path / "neff.bin"
            rbin = self.path / "rho.bin"
            m["zscore_size_GB"] = round(zbin.stat().st_size / 1e9, 2) if zbin.exists() else None
            m["neff_size_GB"] = round(nbin.stat().st_size / 1e9, 2) if nbin.exists() else None
            m["rho_present"] = rbin.exists()
            if not rbin.exists():
                m.setdefault("warnings", []).append(
                    "rho matrix absent — run 'pleiodb rho <db>' to compute"
                )
        except Exception:
            pass
        return m
