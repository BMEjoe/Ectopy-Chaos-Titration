"""
Phase 12 generators (PREREGISTRATION.md section 4): the BP96 Fig. 1 ecological model, the Table 1 systems with fixed
standard parameters and sampling, and the static-transform nulls.  Every generator takes an integer seed (Phase 12
seeds start at 12000) that sets the initial condition; the first `burn` iterations / time units are discarded.
"""
from __future__ import annotations

import numpy as np

N_BENCH = 1000
ENTROPY12 = 20261107


def _rng(code, seed, *more):
    return np.random.default_rng(np.random.SeedSequence([ENTROPY12, code, int(seed), *map(int, more)]))


def _check(y):
    return bool(np.all(np.isfinite(y)) and np.std(y) > 1e-6)


def _retry(make, code, seed):
    """Draw initial conditions until the orbit is finite and non-degenerate (at most 20 draws; registered rule)."""
    rng = _rng(code, seed)
    for attempt in range(20):
        y = make(rng)
        if _check(y):
            return y, attempt
    raise FloatingPointError(f"no admissible orbit for code {code} seed {seed}")


# ------------------------------------------------------------------------------------------- maps
def ecological(seed, n=N_BENCH, lam=118.0, burn=1000):
    """BP96 Fig. 1 legend: x_k = lam y_{k-1} exp(-0.001 (y_{k-1} + x_{k-1}));
    y_k = 0.2 x_{k-1} exp(-0.07 (y_{k-1} + x_{k-1})) + 0.8 y_{k-1} exp(-0.05 (y_{k-1} + 0.5 x_{k-1})).  Observed: y."""
    def make(rng):
        x, y = rng.uniform(10.0, 100.0), rng.uniform(0.1, 1.0)
        out = np.empty(n)
        for k in range(burn + n):
            x, y = (lam * y * np.exp(-0.001 * (y + x)),
                    0.2 * x * np.exp(-0.07 * (y + x)) + 0.8 * y * np.exp(-0.05 * (y + 0.5 * x)))
            if k >= burn:
                out[k - burn] = y
        return out
    return _retry(make, 1, seed)[0]


def logistic(seed, n=N_BENCH, r=4.0, burn=1000):
    def make(rng):
        x = rng.uniform(0.1, 0.9)
        out = np.empty(n)
        for k in range(burn + n):
            x = r * x * (1 - x)
            if k >= burn:
                out[k - burn] = x
        return out
    return _retry(make, 2, seed)[0]


def henon(seed, n=N_BENCH, a=1.4, b=0.3, burn=1000):
    def make(rng):
        u, v = rng.uniform(-0.1, 0.1, 2)
        out = np.empty(n)
        for k in range(burn + n):
            u, v = 1 - a * u * u + v, b * u
            if k >= burn:
                out[k - burn] = u
        return out
    return _retry(make, 3, seed)[0]


def ikeda(seed, n=N_BENCH, u=0.9, burn=1000):
    """z_{n+1} = 1 + u z_n exp(i (0.4 - 6 / (1 + |z_n|^2))); observed Re z."""
    def make(rng):
        z = complex(*rng.uniform(-0.5, 0.5, 2))
        out = np.empty(n)
        for k in range(burn + n):
            z = 1 + u * z * np.exp(1j * (0.4 - 6.0 / (1 + abs(z) ** 2)))
            if k >= burn:
                out[k - burn] = z.real
        return out
    return _retry(make, 4, seed)[0]


def gopy_sna(seed, n=N_BENCH, sigma=1.5, burn=1000):
    """Grebogi, Ott, Pelikan & Yorke (1984) strange non-chaotic attractor (BP96 Table 1 ref. 21):
    x_{n+1} = 2 sigma tanh(x_n) cos(2 pi theta_n), theta_{n+1} = theta_n + omega mod 1, omega = golden mean."""
    om = (np.sqrt(5.0) - 1) / 2

    def make(rng):
        x, th = rng.uniform(0.1, 1.0), rng.uniform(0.0, 1.0)
        out = np.empty(n)
        for k in range(burn + n):
            x, th = 2 * sigma * np.tanh(x) * np.cos(2 * np.pi * th), (th + om) % 1.0
            if k >= burn:
                out[k - burn] = x
        return out
    return _retry(make, 5, seed)[0]


# ------------------------------------------------------------------------------------------- flows (RK4)
def _rk4(f, s, t, dt):
    k1 = f(t, s)
    k2 = f(t + dt / 2, s + dt / 2 * k1)
    k3 = f(t + dt / 2, s + dt / 2 * k2)
    k4 = f(t + dt, s + dt * k3)
    return s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def _flow(f, s0, n, ts, dt, transient, comp=0):
    every = int(round(ts / dt))
    n0 = int(round(transient / dt))
    out = np.empty(n)
    s, t = np.asarray(s0, float), 0.0
    for i in range(n0 + n * every):
        s = _rk4(f, s, t, dt)
        t += dt
        j = i + 1 - n0
        if j > 0 and j % every == 0:
            out[j // every - 1] = s[comp]
    return out


def rossler(seed, n=N_BENCH, a=0.2, b=0.2, c=5.7, ts=0.5, dt=0.01, transient=500.0):
    f = lambda t, s: np.array([-s[1] - s[2], s[0] + a * s[1], b + s[2] * (s[0] - c)])  # noqa: E731
    return _retry(lambda rng: _flow(f, rng.uniform(-1, 1, 3) + [1, 1, 0], n, ts, dt, transient), 6, seed)[0]


def duffing(seed, n=N_BENCH, delta=0.25, gamma=0.3, omega=1.0, transient_periods=100):
    """Double-well forced Duffing x'' + delta x' - x + x^3 = gamma cos(omega t); 12 samples per forcing period."""
    T0 = 2 * np.pi / omega
    f = lambda t, s: np.array([s[1], -delta * s[1] + s[0] - s[0] ** 3 + gamma * np.cos(omega * t)])  # noqa: E731
    return _retry(lambda rng: _flow(f, rng.uniform(-1, 1, 2), n, T0 / 12, T0 / 600, transient_periods * T0), 7, seed)[0]


def lorenz(seed, n=N_BENCH, sigma=10.0, rho=28.0, beta=8.0 / 3.0, ts=0.075, dt=0.0025, transient=100.0):
    """Sampling as the Phase 10 verification generator (ts = 0.075, dt = 0.0025; PB01 settings)."""
    f = lambda t, s: np.array([sigma * (s[1] - s[0]), s[0] * (rho - s[2]) - s[1], s[0] * s[1] - beta * s[2]])  # noqa: E731
    return _retry(lambda rng: _flow(f, rng.uniform([-10, -10, 10], [10, 10, 40]), n, ts, dt, transient), 8, seed)[0]


def mackey_glass(seed, n=N_BENCH, tau=100.0):
    """Phase 10 verification generator (ts = 12.5, dt = 0.05, transient 5,000; initial history from the seed)."""
    from experiments.phase10_final import verify_titration as VT
    return np.asarray(VT.mackey_glass(tau, n=n, seed=int(seed)), float)


SYSTEMS = {"ecological": ecological, "logistic": logistic, "henon": henon, "ikeda": ikeda, "gopy_sna": gopy_sna,
           "rossler": rossler, "duffing": duffing, "lorenz": lorenz, "mackey_glass_100": mackey_glass,
           "mackey_glass_17": lambda seed, n=N_BENCH: mackey_glass(seed, n, tau=17.0)}
TABLE1 = {"logistic": 70, "henon": 70, "ikeda": 65, "ecological": 50, "gopy_sna": 25, "rossler": 75, "duffing": 40,
          "lorenz": 50, "mackey_glass_100": 15, "mackey_glass_17": 15}


# ------------------------------------------------------------------------------------------- static transforms
def ar1(seed, n, phi, burn=500):
    rng = _rng(20, seed, int(round(100 * phi)), n)
    e = rng.standard_normal(n + burn)
    x = np.zeros(n + burn)
    for k in range(1, n + burn):
        x[k] = phi * x[k - 1] + e[k]
    x = x[burn:]
    return (x - x.mean()) / x.std()


TRANSFORMS = {
    "identity": lambda z: z,
    "cube": lambda z: z ** 3,
    "exp": lambda z: np.exp(z),                       # exp(x / sigma) on the standardised series
    "sqrt": lambda z: np.sign(z) * np.abs(z) ** 0.5,
    "tanh": lambda z: np.tanh(2 * z),                 # tanh(2 x / sigma)
    "N5": lambda z: np.exp(0.5 * z),                  # phase5_rr N5 warp: exp(WARP_A z), WARP_A = 0.5, then rescaled
}


def correlated_noise(y, rng):
    """Unit-SD noise with the series' amplitude spectrum (phase-randomised surrogate of y)."""
    from experiments.phase12_original import titration12 as T12
    s = T12.ft_surrogate((y - y.mean()) / y.std(), rng)
    return (s - s.mean()) / s.std()
