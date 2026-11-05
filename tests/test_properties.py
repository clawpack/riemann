"""Property tests for the compiled Riemann solvers.

Each test checks an identity that a correct solver satisfies for *any*
admissible input, so no reference data is stored.  References:

* LeVeque (2002), Finite Volume Methods for Hyperbolic Problems:
  Sec. 4.12 (wave-propagation form of Godunov's method: waves, fluctuations,
  waves summing to the jump); Sec. 15.3 (approximate Riemann solvers: the
  conservation requirement amdq + apdq = f(q_r) - f(q_l), Roe linearization,
  entropy fixes); Sec. 15.5 (f-wave form of the wave-propagation method).
* Bale, LeVeque, Mitran & Rossmanith (2002), SIAM J. Sci. Comput. 24(3),
  955-978: f-wave decomposition (f-waves sum to the flux difference).
* George (2008), J. Comput. Phys. 227(6), 3089-3113: augmented solver and
  preservation of the lake-at-rest steady state.
"""

import numpy as np
import pytest

import rp_harness as rh
from solvers import (ACOUSTICS_2D_JUMP, BY_NAME, DIMENSION_PAIRS, NOT_COVERED,
                     SOLVERS)

TOL = 1e-12
SEED = 20261105


def _directions(solver):
    return [None] if solver.dim == 1 else [1, 2]


CASES = [pytest.param(s, ixy, id=f"{s.name}" + (f"-ixy{ixy}" if ixy else ""))
         for s in SOLVERS for ixy in _directions(s)]


def _run(solver, q, aux=None, ixy=None):
    module = rh.load(solver.name)
    solver.set_params(module)
    return rh.solve(module, solver.fname, q, q, solver.mwaves, aux=aux,
                    ixy=ixy)


def _states(solver, seed=SEED):
    rng = np.random.default_rng(seed)
    q = solver.states(rng, rh.NUM_TOTAL)
    aux = solver.aux(rng, rh.NUM_TOTAL) if solver.aux else None
    return q, aux


# ---- bookkeeping --------------------------------------------------------

def _compiled_modules():
    import pkgutil
    import clawpack.riemann as pkg
    skip = {"riemann_tools", "riemann_interactive", "setup", "static",
            "write_constants"}
    return sorted(m.name for m in pkgutil.iter_modules(pkg.__path__)
                  if m.ispkg is False and m.name not in skip
                  and not m.name.endswith(("_constants", "_py")))


def test_registry_is_complete():
    """Every compiled solver is either tested or listed with a reason."""
    known = set(BY_NAME) | set(NOT_COVERED)
    missing = [m for m in _compiled_modules() if m not in known]
    assert not missing, f"add to tests/solvers.py: {missing}"
    stale = [m for m in known if m not in _compiled_modules()]
    assert not stale, f"no longer compiled, remove from registry: {stale}"


# static.py fwave entries that are missing or not a bool (see README).
STATIC_FWAVE_BAD = {"sw_aug_1D": "missing from static.fwave",
                    "sw_aug_2D": "static.fwave value is 3, not True"}


@pytest.mark.parametrize(
    "solver",
    [pytest.param(s, id=s.name, marks=pytest.mark.xfail(
        strict=True, reason=STATIC_FWAVE_BAD[s.name])
        if s.name in STATIC_FWAVE_BAD else ()) for s in SOLVERS])
def test_registry_matches_static(solver):
    """Registry sizes agree with riemann/static.py, which PyClaw uses."""
    from clawpack.riemann import static
    assert static.num_eqn.get(solver.name) == solver.meqn
    assert static.num_waves.get(solver.name) == solver.mwaves
    assert static.fwave.get(solver.name) is solver.fwave


# <name>_constants.py modules whose sizes disagree with static.py.
CONSTANTS_MISMATCH = {"euler_5wave_2D", "euler_mapgrid_2D",
                      "reactive_euler_with_efix_1D"}


def _constants_modules():
    import pkgutil
    import clawpack.riemann as pkg
    for m in sorted(pkgutil.iter_modules(pkg.__path__), key=lambda m: m.name):
        if m.name.endswith("_constants") and m.name != "write_constants":
            base = m.name[:-len("_constants")]
            marks = ()
            if base in CONSTANTS_MISMATCH:
                marks = pytest.mark.xfail(
                    strict=True, reason="disagrees with static.py (see README)")
            yield pytest.param(base, id=base, marks=marks)


@pytest.mark.parametrize("base", list(_constants_modules()))
def test_constants_match_static(base):
    """<name>_constants.py num_eqn/num_waves agree with static.py."""
    from clawpack.riemann import static
    constants = rh.load(base + "_constants")
    if base not in static.num_eqn:
        pytest.skip(f"{base} is not listed in static.py")
    for attr in ("num_eqn", "num_waves"):
        if hasattr(constants, attr):
            assert getattr(constants, attr) == getattr(static, attr)[base], attr


# Wrappers whose auxl/auxr are not declared as arrays in the Fortran source,
# so f2py exposes them as scalars and makes maux a required argument.
AUX_NOT_ARRAY = {("acoustics_1D", "rp1"), ("kpp_2D", "rpn2"),
                 ("burgers_3D", "rpn3")}


def _solver_functions():
    for name in _compiled_modules():
        module = rh.load(name)
        for fname in ("rp1", "rpn2", "rpn3"):
            if hasattr(module, fname):
                marks = ()
                if (name, fname) in AUX_NOT_ARRAY:
                    marks = pytest.mark.xfail(
                        strict=True,
                        reason="auxl/auxr not declared as arrays (see README)")
                yield pytest.param(name, fname, id=f"{name}.{fname}",
                                   marks=marks)


@pytest.mark.parametrize("name, fname", list(_solver_functions()))
def test_aux_arguments_are_arrays(name, fname):
    """auxl/auxr must be declared arrays so the interface is uniform."""
    assert rh.aux_is_array(getattr(rh.load(name), fname))


# ---- properties ---------------------------------------------------------

@pytest.mark.parametrize("solver, ixy", CASES)
def test_constant_state_gives_zero_fluctuations(solver, ixy):
    """Consistency: identical left and right states produce no fluctuation."""
    q, aux = _states(solver)
    q = np.repeat(q[:, :1], q.shape[1], axis=1)
    if aux is not None:
        aux = np.repeat(aux[:, :1], aux.shape[1], axis=1)
    out = _run(solver, q, aux, ixy)
    for key in ("amdq", "apdq", "wave"):
        rh.assert_close(rh.at_interfaces(out[key]),
                        np.zeros_like(rh.at_interfaces(out[key])), TOL, key)


@pytest.mark.parametrize("solver, ixy",
                         [c for c in CASES if c.values[0].flux or not c.values[0].fwave])
def test_waves_sum_to_jump(solver, ixy):
    """sum_p W_p = q_r - q_l, or for f-waves sum_p Z_p = f(q_r) - f(q_l)."""
    q, aux = _states(solver)
    out = _run(solver, q, aux, ixy)
    total = rh.at_interfaces(out["wave"]).sum(axis=1)
    if solver.fwave:
        fl, fr = rh.left_right(solver.flux(q, ixy or 1))
        expected = fr - fl
    else:
        ql, qr = rh.left_right(q)
        expected = qr - ql
    rows = slice(None)
    if solver.name == "acoustics_2D":
        rows = list(ACOUSTICS_2D_JUMP[ixy])
    rh.assert_close(total[rows], expected[rows], TOL, "sum of waves")


@pytest.mark.parametrize("solver, ixy",
                         [c for c in CASES if c.values[0].flux])
def test_fluctuations_sum_to_flux_difference(solver, ixy):
    """Conservation: amdq + apdq = f(q_r) - f(q_l)."""
    q, aux = _states(solver)
    out = _run(solver, q, aux, ixy)
    fl, fr = rh.left_right(solver.flux(q, ixy or 1))
    rh.assert_close(rh.at_interfaces(out["amdq"] + out["apdq"]), fr - fl, TOL,
                    "amdq + apdq")


@pytest.mark.parametrize("solver", [s for s in SOLVERS if s.swap],
                         ids=lambda s: s.name)
def test_rotation_invariance(solver):
    """ixy=1 on (.., a, b, ..) equals ixy=2 on (.., b, a, ..), permuted back."""
    q, aux = _states(solver)
    mu, mv = solver.swap
    perm = list(range(solver.meqn))
    perm[mu], perm[mv] = perm[mv], perm[mu]
    out1 = _run(solver, q, aux, 1)
    out2 = _run(solver, q[perm], aux, 2)
    rh.assert_close(rh.at_interfaces(out2["s"]), rh.at_interfaces(out1["s"]),
                    TOL, "s")
    for key in ("amdq", "apdq"):
        rh.assert_close(rh.at_interfaces(out2[key][perm]),
                        rh.at_interfaces(out1[key]), TOL, key)
    rh.assert_close(rh.at_interfaces(out2["wave"][perm]),
                    rh.at_interfaces(out1["wave"]), TOL, "wave")


@pytest.mark.parametrize(
    "solver, ixy",
    [c for c in CASES if c.values[0].well_balanced])
def test_lake_at_rest_is_preserved(solver, ixy):
    """Well-balancing: still water over varying bathymetry stays at rest."""
    q, aux = solver.well_balanced(np.random.default_rng(SEED), rh.NUM_TOTAL)
    out = _run(solver, q, aux, ixy)
    for key in ("amdq", "apdq"):
        rh.assert_close(rh.at_interfaces(out[key]),
                        np.zeros_like(rh.at_interfaces(out[key])), TOL, key)


# rpn2_shallow_roe_with_efix skips the 3-wave entropy-fix check whenever the
# shear wave is right-going (s(2) > 0); rp1 always checks it.  The two differ
# at interfaces with a transonic u+c rarefaction and ubar > 0.  Which one is
# intended is a numerical question for the maintainers (see README).
_PAIR_XFAIL = {"shallow_roe_with_efix_2D": "1D/2D entropy-fix logic differs (README)"}


@pytest.mark.parametrize(
    "name1, name2",
    [pytest.param(a, b, id=b, marks=pytest.mark.xfail(strict=True, reason=_PAIR_XFAIL[b])
                  if b in _PAIR_XFAIL else ()) for a, b in DIMENSION_PAIRS])
def test_2d_reduces_to_1d(name1, name2):
    """With zero tangential velocity the 2D solver matches the 1D solver."""
    s1, s2 = BY_NAME[name1], BY_NAME[name2]
    q1, aux = _states(s1)
    # insert a zero tangential momentum after the normal momentum (index 2)
    q2 = np.insert(q1, 2, 0.0, axis=0)
    rows = [i for i in range(s2.meqn) if i != 2]
    out1 = _run(s1, q1, aux)
    out2 = _run(s2, q2, aux, 1)
    for key in ("amdq", "apdq"):
        rh.assert_close(rh.at_interfaces(out2[key][rows]),
                        rh.at_interfaces(out1[key]), TOL, key)
        rh.assert_close(rh.at_interfaces(out2[key][2]),
                        np.zeros(rh.NUM_TOTAL - 1), TOL, key + " tangential")
