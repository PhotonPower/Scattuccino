"""Counter-based random numbers: one independent stream per packet (stage S4).

u = mix(seed, packet id, draw number) with the SplitMix64 finaliser, mapped to (0, 1) with 53 bits.
The k-th number of a packet does not depend on how many other packets exist, in which order they
are processed or how many voxel faces they cross – results are invariant under batching and under
refining voxels of one material, and a GPU kernel can reproduce them bit for bit (one stream per
thread).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from scattuccino._types import FloatArray, IntArray

_M1 = np.uint64(0xBF58476D1CE4E5B9)
_M2 = np.uint64(0x94D049BB133111EB)
_G = np.uint64(0x9E3779B97F4A7C15)


def _mix(x: Any) -> Any:
    x = (x ^ (x >> np.uint64(30))) * _M1
    x = (x ^ (x >> np.uint64(27))) * _M2
    return x ^ (x >> np.uint64(31))


class CounterRandom:
    """uniform(ids) returns one number per id and advances that id's counter."""

    def __init__(self, seed: int, packets: int) -> None:
        self.seed = np.uint64(seed)
        self.counter = np.zeros(packets, dtype=np.uint64)

    def uniform(self, ids: IntArray) -> FloatArray:
        ids = np.asarray(ids, dtype=np.int64)
        with np.errstate(over="ignore"):
            key = _mix(self.seed * _G + np.uint64(1)) ^ _mix(
                ids.astype(np.uint64) * _G + np.uint64(2)
            )
            x = _mix(key ^ (self.counter[ids] * _G + np.uint64(3)))
        self.counter[ids] += np.uint64(1)
        return np.asarray(((x >> np.uint64(11)).astype(np.float64) + 0.5) * 2.0**-53)
