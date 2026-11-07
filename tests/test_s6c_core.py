"""S6c: polarised transport in the C++ core against the Python reference, packet by packet (T7e)."""

from __future__ import annotations

import os

import numpy as np
import pytest

from scattuccino import fast
from scattuccino.polarisation import mie_mueller, rayleigh_mueller, simulate_polarised

if not fast.CORE_AVAILABLE:
    if os.environ.get("SCATTUCCINO_REQUIRE_CORE") == "1":
        raise RuntimeError("SCATTUCCINO_REQUIRE_CORE=1 but scattuccino._core is not built")
    pytest.skip("scattuccino._core not built", allow_module_level=True)

MIE = mie_mueller(1.59 / 1.33 + 0j, 2 * np.pi * 0.5 * 1.33 / 0.633)
CASES = [
    # matched indices, normal incidence, linear
    dict(
        d=1.0,
        mua=0.1,
        mus=10.0,
        table=MIE,
        stokes=(1.0, 1.0, 0.0, 0.0),
        n_slab=1.0,
        n_below=1.0,
        ci=1.0,
    ),
    # index steps, oblique, elliptical; Fresnel with complex amplitudes and total reflection
    dict(
        d=1.0,
        mua=0.1,
        mus=10.0,
        table=MIE,
        stokes=(1.0, 0.3, 0.5, 0.6),
        n_slab=1.4,
        n_below=1.5,
        ci=float(np.cos(0.6)),
    ),
    # Rayleigh scatterers, half space with index step
    dict(
        d=np.inf,
        mua=0.05,
        mus=5.0,
        table=rayleigh_mueller(),
        stokes=(1.0, 0.0, 0.0, 1.0),
        n_slab=1.33,
        n_below=1.0,
        ci=1.0,
    ),
]


@pytest.mark.parametrize("case", CASES)
def test_core_equals_python_packet_by_packet(case: dict) -> None:  # type: ignore[type-arg]
    kw = dict(n_slab=case["n_slab"], n_below=case["n_below"], cos_incident=case["ci"])
    s = np.array(case["stokes"])
    a = simulate_polarised(
        case["d"], case["mua"], case["mus"], case["table"], s, 4_000, seed=7, **kw
    )
    b = fast.simulate_polarised_slab(
        case["d"], case["mua"], case["mus"], case["table"], s, 4_000, seed=7, threads=3, **kw
    )
    assert np.array_equal(a.side, b.side)
    assert np.array_equal(a.scatterings, b.scatterings)
    assert np.array_equal(a.weight, b.weight)
    for x, y in ((a.stokes, b.stokes), (a.u, b.u), (a.v, b.v)):
        assert np.max(np.abs(x - y)) < 1e-10
    assert a.specular is not None and b.specular is not None
    assert abs(a.specular[0] - b.specular[0]) < 1e-15
    assert np.max(np.abs(a.specular[1] - b.specular[1])) < 1e-14
