"""
Low-level chunked compressed binary storage.

Layout per matrix (e.g. 'zscore'):
  {name}.bin   — concatenated zstd-compressed chunks, row-major chunk order
  {name}.cidx  — uint64 array of byte offsets, length = n_chunks + 1

Chunk id for chunk (vi, ti):  vi * n_t_chunks + ti
cidx[chunk_id]   = start byte in .bin
cidx[chunk_id+1] = end byte   in .bin
"""

from __future__ import annotations

import os
import struct
from pathlib import Path
from typing import Sequence, Tuple

import numpy as np
import zstandard as zstd

ZSTD_LEVEL = 3
_COMPRESSOR = zstd.ZstdCompressor(level=ZSTD_LEVEL, threads=-1)
_DECOMPRESSOR = zstd.ZstdDecompressor()


class ChunkedMatrix:
    """
    Random-access, compressed 2-D matrix backed by two files (.bin + .cidx).

    dtype must be a fixed-width numpy dtype.  NA values are communicated by
    the caller via sentinel constants (e.g. INT16_MIN for z-scores).
    """

    def __init__(
        self,
        base_path: Path,
        shape: Tuple[int, int],
        dtype: np.dtype | str,
        chunk_shape: Tuple[int, int] = (512, 512),
    ):
        self.base_path = Path(base_path)
        self.shape = shape
        self.dtype = np.dtype(dtype)
        self.chunk_shape = chunk_shape
        self.V, self.T = shape
        self.CV, self.CT = chunk_shape
        self.n_v_chunks = (self.V + self.CV - 1) // self.CV
        self.n_t_chunks = (self.T + self.CT - 1) // self.CT
        self.n_chunks = self.n_v_chunks * self.n_t_chunks

        self._bin_path = self.base_path.with_suffix(".bin")
        self._cidx_path = self.base_path.with_suffix(".cidx")
        self._cidx_cache: np.ndarray | None = None

    # ------------------------------------------------------------------
    # Build API
    # ------------------------------------------------------------------

    def open_write(self) -> None:
        self._wfp = open(self._bin_path, "wb")
        self._offsets: list[int] = [0]

    def write_chunk(self, vi: int, ti: int, data: np.ndarray) -> None:
        """Append one chunk.  data shape must match the actual chunk footprint."""
        blob = _COMPRESSOR.compress(np.ascontiguousarray(data, self.dtype).tobytes())
        self._wfp.write(blob)
        self._offsets.append(self._offsets[-1] + len(blob))

    def close_write(self) -> None:
        self._wfp.close()
        np.array(self._offsets, dtype=np.uint64).tofile(self._cidx_path)

    # ------------------------------------------------------------------
    # Query API
    # ------------------------------------------------------------------

    def _load_cidx(self) -> np.ndarray:
        if self._cidx_cache is None:
            self._cidx_cache = np.fromfile(self._cidx_path, dtype=np.uint64)
        return self._cidx_cache

    def _chunk_id(self, vi: int, ti: int) -> int:
        return vi * self.n_t_chunks + ti

    def _chunk_bounds(self, vi: int, ti: int) -> Tuple[int, int, int, int]:
        v0 = vi * self.CV
        v1 = min(v0 + self.CV, self.V)
        t0 = ti * self.CT
        t1 = min(t0 + self.CT, self.T)
        return v0, v1, t0, t1

    def get_raw_chunk(self, vi: int, ti: int) -> np.ndarray:
        """Return one decompressed chunk in its native dtype."""
        cidx = self._load_cidx()
        cid = self._chunk_id(vi, ti)
        start, end = int(cidx[cid]), int(cidx[cid + 1])
        v0, v1, t0, t1 = self._chunk_bounds(vi, ti)
        with open(self._bin_path, "rb") as fh:
            fh.seek(start)
            blob = fh.read(end - start)
        raw = _DECOMPRESSOR.decompress(blob)
        return np.frombuffer(raw, dtype=self.dtype).reshape(v1 - v0, t1 - t0)

    def get_block(self, v_start: int, v_end: int, t_start: int, t_end: int) -> np.ndarray:
        """Return a rectangular block in the native dtype."""
        v_end = min(v_end, self.V)
        t_end = min(t_end, self.T)
        result = np.empty((v_end - v_start, t_end - t_start), dtype=self.dtype)

        vi0, vi1 = v_start // self.CV, (v_end - 1) // self.CV + 1
        ti0, ti1 = t_start // self.CT, (t_end - 1) // self.CT + 1

        # Sort chunk reads by file offset to maximise sequential throughput
        cidx = self._load_cidx()
        order = sorted(
            [(vi, ti) for vi in range(vi0, vi1) for ti in range(ti0, ti1)],
            key=lambda p: int(cidx[self._chunk_id(*p)]),
        )

        with open(self._bin_path, "rb") as fh:
            for vi, ti in order:
                cid = self._chunk_id(vi, ti)
                start, end_b = int(cidx[cid]), int(cidx[cid + 1])
                cv0, cv1, ct0, ct1 = self._chunk_bounds(vi, ti)

                fh.seek(start)
                blob = fh.read(end_b - start)
                chunk = np.frombuffer(
                    _DECOMPRESSOR.decompress(blob), dtype=self.dtype
                ).reshape(cv1 - cv0, ct1 - ct0)

                ov0, ov1 = max(v_start, cv0), min(v_end, cv1)
                ot0, ot1 = max(t_start, ct0), min(t_end, ct1)

                result[ov0 - v_start : ov1 - v_start, ot0 - t_start : ot1 - t_start] = (
                    chunk[ov0 - cv0 : ov1 - cv0, ot0 - ct0 : ot1 - ct0]
                )

        return result

    def get_rows(self, v_indices: Sequence[int]) -> np.ndarray:
        """Retrieve arbitrary rows (may not be contiguous)."""
        v_indices = np.asarray(v_indices, dtype=np.int64)
        result = np.empty((len(v_indices), self.T), dtype=self.dtype)
        vi_groups: dict[int, list[int]] = {}
        for pos, vi in enumerate(v_indices // self.CV):
            vi_groups.setdefault(int(vi), []).append(pos)

        for vi, positions in vi_groups.items():
            block = self.get_block(vi * self.CV, (vi + 1) * self.CV, 0, self.T)
            for pos in positions:
                row_in_chunk = int(v_indices[pos]) - vi * self.CV
                result[pos] = block[row_in_chunk]
        return result

    def get_cols(self, t_indices: Sequence[int]) -> np.ndarray:
        """Retrieve arbitrary columns (may not be contiguous)."""
        t_indices = np.asarray(t_indices, dtype=np.int64)
        result = np.empty((self.V, len(t_indices)), dtype=self.dtype)
        ti_groups: dict[int, list[int]] = {}
        for pos, ti in enumerate(t_indices // self.CT):
            ti_groups.setdefault(int(ti), []).append(pos)

        for ti, positions in ti_groups.items():
            block = self.get_block(0, self.V, ti * self.CT, (ti + 1) * self.CT)
            for pos in positions:
                col_in_chunk = int(t_indices[pos]) - ti * self.CT
                result[:, pos] = block[:, col_in_chunk]
        return result

    # ------------------------------------------------------------------
    # Metadata persistence
    # ------------------------------------------------------------------

    def meta_dict(self) -> dict:
        return {
            "shape": list(self.shape),
            "dtype": self.dtype.str,
            "chunk_shape": list(self.chunk_shape),
        }

    @classmethod
    def from_meta(cls, base_path: Path, meta: dict) -> "ChunkedMatrix":
        return cls(
            base_path,
            shape=tuple(meta["shape"]),
            dtype=np.dtype(meta["dtype"]),
            chunk_shape=tuple(meta["chunk_shape"]),
        )
