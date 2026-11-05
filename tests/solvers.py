"""Registry of solvers covered by the property tests.

Each entry supplies what the generic tests cannot infer: physical parameters
(set through the module's Fortran common block), a generator of physically
admissible states, and, where known, the flux function.  Flux functions here
are the standard textbook forms of each system (LeVeque 2002, Finite Volume
Methods for Hyperbolic Problems, Ch. 2-3, 13-14, 18); they are written from
the equations, not taken from the Fortran.

Every compiled solver module must appear either in ``SOLVERS`` or in
``NOT_COVERED`` with a reason; ``test_registry_is_complete`` enforces this so
new solvers are noticed.
"""

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np


@dataclass
class Solver:
    name: str                          # clawpack.riemann module name
    dim: int                           # 1, 2 or 3
    meqn: int
    mwaves: int
    params: dict = field(default_factory=dict)   # {common_block: {name: value}}
    states: Callable = None            # (rng, n) -> q, shape (meqn, n)
    flux: Optional[Callable] = None    # (q, ixy) -> f, shape (meqn, n)
    fwave: bool = False                # waves are f-waves (sum = flux jump)
    maux: int = 0
    aux: Optional[Callable] = None     # (rng, n) -> aux, shape (maux, n)
    jump_components: Optional[tuple] = None   # restrict sum-of-waves check
    swap: Optional[tuple] = None       # (mu, mv) normal/tangential momenta (2D)
    well_balanced: Optional[Callable] = None  # (rng, n) -> (q, aux) at rest
    speed_bounds: Optional[Callable] = None   # (q_l, q_r, ixy) -> (s_min, s_max)
    fname: str = ""

    def __post_init__(self):
        if not self.fname:
            self.fname = {1: "rp1", 2: "rpn2", 3: "rpn3"}[self.dim]

    def set_params(self, module):
        for block, values in self.params.items():
            common = getattr(module, block)
            for key, value in values.items():
                setattr(common, key, value)


GAMMA = 1.4
GRAV = 1.0
RHO, BULK = 2.0, 3.0           # acoustics: deliberately not 1 so terms differ
CC, ZZ = np.sqrt(BULK / RHO), np.sqrt(BULK * RHO)
ACOUSTICS = {"cparam": {"rho": RHO, "bulk": BULK, "cc": CC, "zz": ZZ}}


# ---- state generators ---------------------------------------------------

def normal(meqn):
    return lambda rng, n: rng.standard_normal((meqn, n))


def shallow(dim):
    def gen(rng, n):
        h = 0.5 + 1.5 * rng.random(n)
        vel = rng.standard_normal((dim, n))
        return np.vstack([h, h * vel])
    return gen


def euler(dim, extra=0):
    def gen(rng, n):
        rho = 0.5 + 1.5 * rng.random(n)
        vel = rng.standard_normal((dim, n))
        p = 0.5 + 1.5 * rng.random(n)
        energy = p / (GAMMA - 1) + 0.5 * rho * np.sum(vel**2, axis=0)
        rows = [rho, *(rho * vel), energy]
        rows += [rho * rng.random(n) for _ in range(extra)]
        return np.vstack(rows)
    return gen


def shallow_tracer(rng, n):
    return np.vstack([shallow(1)(rng, n), rng.random(n)])


def flat_bottom(rng, n):
    return np.full((1, n), -10.0)


def lake_at_rest(dim):
    """Still water over varying bathymetry: h + b = 0, zero momentum."""
    def gen(rng, n):
        b = -(0.5 + 1.5 * rng.random(n))
        q = np.vstack([-b, np.zeros((dim, n))])
        return q, b[np.newaxis, :]
    return gen


# ---- flux functions -----------------------------------------------------

def f_shallow(q, ixy=1):
    h = q[0]
    vel = q[1:] / h
    mu = ixy                            # normal momentum component index
    f = q * vel[mu - 1]
    f[mu] += 0.5 * GRAV * h**2
    return f


def f_euler(q, ixy=1, dim=None):
    dim = q.shape[0] - 2 if dim is None else dim
    rho, energy = q[0], q[dim + 1]
    vel = q[1:dim + 1] / rho
    p = (GAMMA - 1) * (energy - 0.5 * rho * np.sum(vel**2, axis=0))
    un = vel[ixy - 1]
    f = q * un
    f[ixy] += p
    f[dim + 1] += p * un
    return f


def f_acoustics_1d(q, ixy=1):
    p, u = q
    return np.array([BULK * u, p / RHO])


def f_acoustics_2d(q, ixy=1):
    p, u, v = q
    un = (u, v)[ixy - 1]
    f = np.zeros_like(q)
    f[0] = BULK * un
    f[ixy] = p / RHO
    return f


# ---- HLLE wave-speed estimates ------------------------------------------
# Einfeldt (1988): the HLLE speeds are s_1 = min(lambda_min(q_l), lambda_min(Roe))
# and s_2 = max(lambda_max(q_r), lambda_max(Roe)), with the Roe-averaged
# state's eigenvalues.  Sound speeds here come from the equation of state,
# not from the Fortran.

def _roe_weights(a_l, a_r):
    w_l, w_r = np.sqrt(a_l), np.sqrt(a_r)
    return w_l / (w_l + w_r), w_r / (w_l + w_r)


def einfeldt_shallow(q_l, q_r, ixy=1):
    h_l, h_r = q_l[0], q_r[0]
    u_l, u_r = q_l[ixy] / h_l, q_r[ixy] / h_r
    c_l, c_r = np.sqrt(GRAV * h_l), np.sqrt(GRAV * h_r)
    a, b = _roe_weights(h_l, h_r)
    u_hat = a * u_l + b * u_r
    c_hat = np.sqrt(GRAV * 0.5 * (h_l + h_r))
    return (np.minimum(u_l - c_l, u_hat - c_hat),
            np.maximum(u_r + c_r, u_hat + c_hat))


def einfeldt_euler(q_l, q_r, ixy=1):
    dim = q_l.shape[0] - 2

    def primitive(q):
        rho, energy = q[0], q[dim + 1]
        vel = q[1:dim + 1] / rho
        p = (GAMMA - 1) * (energy - 0.5 * rho * np.sum(vel**2, axis=0))
        return rho, vel, (energy + p) / rho, np.sqrt(GAMMA * p / rho)

    rho_l, vel_l, H_l, c_l = primitive(q_l)
    rho_r, vel_r, H_r, c_r = primitive(q_r)
    a, b = _roe_weights(rho_l, rho_r)
    vel_hat = a * vel_l + b * vel_r
    H_hat = a * H_l + b * H_r
    c_hat = np.sqrt((GAMMA - 1) * (H_hat - 0.5 * np.sum(vel_hat**2, axis=0)))
    un_l, un_r, un_hat = vel_l[ixy - 1], vel_r[ixy - 1], vel_hat[ixy - 1]
    return (np.minimum(un_l - c_l, un_hat - c_hat),
            np.maximum(un_r + c_r, un_hat + c_hat))


# ---- the registry -------------------------------------------------------

SOLVERS = [
    # 1D scalar
    Solver("advection_1D", 1, 1, 1, {"cparam": {"u": 0.7}}, normal(1),
           lambda q, ixy=1: 0.7 * q),
    Solver("burgers_1D", 1, 1, 1, {}, normal(1), lambda q, ixy=1: 0.5 * q**2),
    Solver("cubic_1D", 1, 1, 1, {}, normal(1), lambda q, ixy=1: q**3),
    Solver("traffic_1D", 1, 1, 1, {"cparam": {"umax": 1.0}},
           lambda rng, n: rng.random((1, n)), lambda q, ixy=1: q * (1 - q)),
    # 1D systems
    Solver("acoustics_1D", 1, 2, 2, ACOUSTICS, normal(2), f_acoustics_1d),
    Solver("euler_with_efix_1D", 1, 3, 3, {"cparam": {"gamma": GAMMA}},
           euler(1), f_euler),
    Solver("euler_hlle_1D", 1, 3, 2, {"cparam": {"gamma": GAMMA}},
           euler(1), f_euler, speed_bounds=einfeldt_euler),
    Solver("shallow_roe_with_efix_1D", 1, 2, 2, {"cparam": {"grav": GRAV}},
           shallow(1), f_shallow),
    Solver("shallow_hlle_1D", 1, 2, 2, {"cparam": {"grav": GRAV}},
           shallow(1), f_shallow, speed_bounds=einfeldt_shallow),
    Solver("shallow_roe_tracer_1D", 1, 3, 3, {"cparam": {"grav": GRAV}},
           shallow_tracer),
    Solver("shallow_bathymetry_fwave_1D", 1, 2, 2,
           {"cparam": {"grav": GRAV, "dry_tolerance": 1e-3, "sea_level": 0.0}},
           shallow(1), f_shallow, fwave=True, maux=1, aux=flat_bottom,
           well_balanced=lake_at_rest(1)),
    Solver("sw_aug_1D", 1, 2, 2,      # 3 augmented waves folded into 2
           {"cparam": {"grav": GRAV, "dry_tolerance": 1e-3}},
           shallow(1), f_shallow, fwave=True, maux=1, aux=flat_bottom,
           well_balanced=lake_at_rest(1)),
    # 2D normal solvers
    Solver("acoustics_2D", 2, 3, 2, ACOUSTICS, normal(3), f_acoustics_2d,
           swap=(1, 2), jump_components=None),
    Solver("shallow_roe_with_efix_2D", 2, 3, 3, {"cparam": {"grav": GRAV}},
           shallow(2), f_shallow, swap=(1, 2)),
    Solver("shallow_hlle_2D", 2, 3, 2, {"cparam": {"grav": GRAV}},
           shallow(2), f_shallow, swap=(1, 2), speed_bounds=einfeldt_shallow),
    Solver("euler_4wave_2D", 2, 4, 4, {"cparam": {"gamma": GAMMA}},
           euler(2), f_euler, swap=(1, 2)),
    Solver("euler_hlle_2D", 2, 4, 2, {"cparam": {"gamma": GAMMA}},
           euler(2), f_euler, swap=(1, 2), speed_bounds=einfeldt_euler),
    # 5 waves (4 Euler + tracer contact): matches static.py and the source;
    # euler_5wave_2D_constants.py says 4 (see README, "Known issues").
    Solver("euler_5wave_2D", 2, 5, 5, {"cparam": {"gamma": GAMMA}},
           euler(2, extra=1), swap=(1, 2)),
    Solver("shallow_bathymetry_fwave_2D", 2, 3, 3,
           {"cparam": {"grav": GRAV, "dry_tolerance": 1e-3, "sea_level": 0.0}},
           shallow(2), f_shallow, fwave=True, maux=1, aux=flat_bottom,
           swap=(1, 2), well_balanced=lake_at_rest(2)),
    Solver("sw_aug_2D", 2, 3, 3, {"cparam": {"grav": GRAV}},
           shallow(2), f_shallow, fwave=True, maux=1, aux=flat_bottom,
           swap=(1, 2), well_balanced=lake_at_rest(2)),
]

# acoustics_2D carries only the two acoustic waves; the zero-speed wave for
# the tangential velocity is omitted, so the waves sum to the jump only in the
# pressure and normal-velocity components.  Checked per direction in the test.
ACOUSTICS_2D_JUMP = {1: (0, 1), 2: (0, 2)}

# 1D/2D pairs that solve the same physics: with zero tangential velocity the
# 2D normal solver must reproduce the 1D solver (components mapped below).
DIMENSION_PAIRS = [
    ("shallow_roe_with_efix_1D", "shallow_roe_with_efix_2D"),
    ("shallow_hlle_1D", "shallow_hlle_2D"),
    ("euler_hlle_1D", "euler_hlle_2D"),
    ("shallow_bathymetry_fwave_1D", "shallow_bathymetry_fwave_2D"),
    ("sw_aug_1D", "sw_aug_2D"),
]

NOT_COVERED = {
    # needs per-cell coefficients in aux and a matching flux
    "acoustics_variable_1D": "variable coefficients in aux",
    "advection_color_1D": "velocity in aux",
    "traffic_vc_1D": "variable speed limit in aux",
    "traffic_vc_fwave_1D": "variable speed limit in aux",
    "traffic_vc_tracer_1D": "variable speed limit in aux",
    "nonlinear_elasticity_fwave_1D": "material properties in aux",
    "psystem_fwave_1D": "material properties in aux",
    "vc_acoustics_2D": "variable coefficients in aux",
    "vc_advection_2D": "velocities in aux",
    "vc_elasticity_2D": "material properties in aux",
    "acoustics_mapped_2D": "mapped-grid metrics in aux",
    "euler_mapgrid_2D": "mapped-grid metrics in aux",
    "shallow_sphere_2D": "sphere metrics in aux and extra common blocks",
    "psystem_2D": "material properties in aux",
    # need an admissible-state model not yet written
    "reactive_euler_with_efix_1D": "needs reactive-Euler state model; wave "
                                   "count disagrees with static.py (README)",
    "mhd_roe_1D": "needs MHD state and flux model",
    "euler_hlle_with_walls_2D": "wall flags in aux",
    "advection_2D": "needs direction-dependent velocity handling",
    "burgers_2D": "flux form not yet confirmed",
    "kpp_2D": "flux form not yet confirmed",
    # 3D normal solvers: next phase
    "burgers_3D": "3D not yet covered",
    "euler_3D": "3D not yet covered",
    "euler_mapgrid_3D": "3D not yet covered",
    "vc_acoustics_3D": "3D not yet covered",
    "vc_advection_3D": "3D not yet covered",
    # pointwise wrappers (rp1_ptwise.f90 drivers)
    "acoustics_1D_ptwise": "pointwise driver; covered via acoustics_1D physics later",
    "advection_1D_ptwise": "pointwise driver",
    "acoustics_2D_ptwise": "pointwise driver",
    # not a solver
    "mappedGrid": "helper module, not a Riemann solver",
}

BY_NAME = {s.name: s for s in SOLVERS}
