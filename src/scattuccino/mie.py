"""Mie series for homogeneous spheres (Bohren & Huffman 1983, ch. 4), time factor e^(−iωt).

m = n_sphere/n_medium (complex, Im ≥ 0 absorbing), x = 2π r n_medium/λ₀.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.special import spherical_jn, spherical_yn

from scattuccino._types import FloatArray

ComplexArray = NDArray[np.complex128]


def coefficients(m: complex, x: float) -> tuple[FloatArray, ComplexArray, ComplexArray]:
    nmax = int(x + 4 * x ** (1 / 3) + 2)
    n = np.arange(1, nmax + 1)
    mx = m * x
    jn, yn = spherical_jn(n, x), spherical_yn(n, x)
    jnp, ynp = spherical_jn(n, x, derivative=True), spherical_yn(n, x, derivative=True)
    hn, hnp = jn + 1j * yn, jnp + 1j * ynp
    jm, jmp = spherical_jn(n, mx), spherical_jn(n, mx, derivative=True)
    psi_x, dpsi_x = x * jn, jn + x * jnp
    xi_x, dxi_x = x * hn, hn + x * hnp
    psi_m, dpsi_m = mx * jm, jm + mx * jmp
    a = (m * psi_m * dpsi_x - psi_x * dpsi_m) / (m * psi_m * dxi_x - xi_x * dpsi_m)
    b = (psi_m * dpsi_x - m * psi_x * dpsi_m) / (psi_m * dxi_x - m * xi_x * dpsi_m)
    return n.astype(np.float64), a, b


def s12(m: complex, x: float, theta: FloatArray) -> tuple[ComplexArray, ComplexArray]:
    n, a, b = coefficients(m, x)
    mu = np.cos(np.asarray(theta, dtype=np.float64))
    nn = len(n)
    pi = np.zeros((nn, mu.size))
    tau = np.zeros_like(pi)
    pi[0], tau[0] = 1.0, mu
    if nn > 1:
        pi[1] = 3 * mu
        tau[1] = 2 * mu * pi[1] - 3 * pi[0]
    for k in range(2, nn):
        j = k + 1
        pi[k] = ((2 * j - 1) * mu * pi[k - 1] - j * pi[k - 2]) / (j - 1)
        tau[k] = j * mu * pi[k] - (j + 1) * pi[k - 1]
    f = ((2 * n + 1) / (n * (n + 1)))[:, None]
    return (
        np.sum(f * (a[:, None] * pi + b[:, None] * tau), axis=0),
        np.sum(f * (a[:, None] * tau + b[:, None] * pi), axis=0),
    )


def qeff(m: complex, x: float) -> tuple[float, float]:
    """(Q_ext, Q_sca)."""
    n, a, b = coefficients(m, x)
    q_ext = 2 / x**2 * np.sum((2 * n + 1) * np.real(a + b))
    q_sca = 2 / x**2 * np.sum((2 * n + 1) * (np.abs(a) ** 2 + np.abs(b) ** 2))
    return float(q_ext), float(q_sca)
