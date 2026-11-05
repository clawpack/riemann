# Riemann solver property tests

These tests call the compiled solvers directly. The meta-repository's meson
build turns each solver in `src/` into an f2py module
(`clawpack.riemann.<name>`), and the tests check identities that a correct
solver satisfies for *any* admissible input. No reference output files are
stored. A failure points at one solver and one property, not at a whole
simulation.

```
pip install --no-build-isolation -e .   # from the clawpack meta-repo root
cd riemann && pytest tests
```

## What is checked

| Test | Property | Reference |
|---|---|---|
| `test_constant_state_gives_zero_fluctuations` | consistency: q_l = q_r gives no waves and no fluctuations | LeVeque (2002) §4.12 |
| `test_waves_sum_to_jump` | waves sum to q_r − q_l; f-waves sum to f(q_r) − f(q_l) | LeVeque (2002) §4.12, §15.5; Bale et al. (2002) |
| `test_fluctuations_sum_to_flux_difference` | conservation: amdq + apdq = f(q_r) − f(q_l) | LeVeque (2002) §15.3 |
| `test_rotation_invariance` | `rpn2` with `ixy=1` and `ixy=2` agree once the momenta are swapped | rotational invariance of the equations |
| `test_lake_at_rest_is_preserved` | still water over varying bathymetry gives zero fluctuations | George (2008) |
| `test_hlle_speeds_match_einfeldt` | HLLE wave speeds equal the Einfeldt estimates (true left/right sound speeds and Roe-averaged speeds) | Einfeldt (1988) |
| `test_2d_reduces_to_1d` | `rpn2` with zero tangential velocity reproduces the matching `rp1` | |
| `test_registry_matches_static` | the test registry matches `static.py`, which PyClaw uses to size arrays | |
| `test_constants_match_static` | `<name>_constants.py` sizes match `static.py` | |
| `test_aux_arguments_are_arrays` | wrappers take `auxl`/`auxr` as arrays (uniform interface) | |
| `test_registry_is_complete` | every compiled module is either tested or listed in `NOT_COVERED` with a reason | |

States are random but seeded (default 40 interior cells plus ghost cells),
and checks use an absolute tolerance of 1e-12 scaled by the size of the data.
Flux functions in `solvers.py` are the textbook forms of each system, written
from the equations rather than from the Fortran.

To confirm that the checks can fail, these deliberate errors were each
caught when the tests were written: gravity off by 1%, γ off by 1e-6, the
wrong momentum permutation, a 1e-6 bump on a lake at rest, and an
inconsistent acoustic impedance.

**References**
- R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*, Cambridge, 2002.
- D. S. Bale, R. J. LeVeque, S. Mitran, J. A. Rossmanith, "A wave propagation
  method for conservation laws and balance laws with spatially varying flux
  functions", *SIAM J. Sci. Comput.* 24(3), 955–978, 2002.
- B. Einfeldt, "On Godunov-type methods for gas dynamics", *SIAM J. Numer.
  Anal.* 25(2), 294–318, 1988.
- D. L. George, "Augmented Riemann solvers for the shallow water equations
  over variable topography with steady states and inundation", *J. Comput.
  Phys.* 227(6), 3089–3113, 2008.

## Coverage

`solvers.py` lists the covered solvers: 1D scalar and system solvers and the
2D normal solvers (`rpn2`) for shallow water, Euler and acoustics. It also
lists every compiled module that is not yet covered, with the reason. Not yet
covered:
- transverse solvers (`rpt2`, `rpt3`, `rptt3`)
- all 3D solvers
- variable-coefficient solvers and mapped-grid solvers, which need a model of
  their `aux` arrays
- the pointwise wrappers

## GeoClaw's solvers are NOT tested here

This suite does **not** test GeoClaw's Riemann solvers, and it should. The
solvers GeoClaw actually runs live in this repository:

- `src/rpn2_geoclaw.f`, `src/rpt2_geoclaw.f`, `src/geoclaw_riemann_utils.f`
- the `_vec` and adjoint variants
- `src/rpn2_layered_shallow_water.f90` and its transverse solver, used by
  GeoClaw's multilayer code

But they are compiled only by GeoClaw's Makefiles, never by the meson build
that produces the modules tested here.

`sw_aug_1D`/`sw_aug_2D` are tested, and they contain an augmented solver also
named `riemann_aug_JCP`, but it is a **separate, forked copy**. Its
signature differs from GeoClaw's (GeoClaw's adds `pL, pR, rho` for
atmospheric-pressure forcing), so passing tests here say nothing about
GeoClaw's version.

**What it would take to include them.** The GeoClaw solvers `use` GeoClaw
modules for their parameters:

| Solver | Modules used |
|---|---|
| `rpn2_geoclaw.f` | `geoclaw_module` (`grav`, `dry_tolerance`, `rho`, `earth_radius`, `deg2rad`), `amr_module` (`mcapa`), `met_forcing_module` (`pressure_forcing`, `pressure_index`) |
| `rpt2_geoclaw.f` | `geoclaw_module` (`grav`, `dry_tolerance`, `coordinate_system`, `earth_radius`, `deg2rad`) |
| `rpn2_layered_shallow_water.f90` | the above plus `multilayer_module` |

Steps:
1. **Build a test module.** Compile those solver files together with the
   GeoClaw module sources (or minimal test-only stand-ins that declare the
   same variables) into an f2py extension. For example, add a
   `geoclaw_2D_test` target to the meson build, or build it in a pytest
   fixture with `f2py -c`. The stand-in modules must match the real ones'
   interfaces exactly, or the test drifts from what GeoClaw compiles. That
   argues for building against GeoClaw's real module sources from a GeoClaw
   checkout.
2. **Make the parameters settable from Python** (module variables are
   reachable through f2py), and set `coordinate_system = 1` (Cartesian) so
   the spherical metric terms drop out.
3. **Register the solver in `solvers.py`**: q = (h, hu, hv), aux(1) =
   bathymetry, fwave = True. The existing consistency, f-wave sum,
   conservation (flat bottom), rotation and lake-at-rest tests then apply
   unchanged.
4. **Add GeoClaw-specific properties:**
   - dry states and wet/dry interfaces (no negative depths, zero
     fluctuations between two dry cells)
   - pressure forcing: a pressure gradient balanced by a surface slope, with
     zero flow, should give zero fluctuations
   - agreement with `sw_aug_2D` on wet, flat-bottom, unforced states. That
     last check shows whether the two `riemann_aug_JCP` copies have diverged.
5. **Decide where it runs:** in this repository's CI (it needs a GeoClaw
   checkout) or in GeoClaw's (where changes to `geoclaw_module` happen).

## Known issues found by these tests

Each issue is marked as a strict expected failure, so fixing it turns the
test red until the mark is removed. The one exception is noted under the HLLE
entry below.

- **`auxl`/`auxr` not declared as arrays** in `rp1_acoustics.f90` and
  `rpn2_kpp.f90` (not declared at all) and in `rpn3_burgers.f90` (declared as
  scalars). Unused arguments let gfortran accept this, but f2py then wraps
  them as scalar floats and makes `maux` a required positional argument, so
  these Python signatures differ from every other solver's. `rpn2_kpp.f90`
  also lacks `intent(out)` on its outputs. Calls from Fortran are unaffected.
- **`static.py` / `*_constants.py` disagreements.**
  - `euler_5wave_2D_constants.num_waves` is 4; the solver writes 5 waves, and
    `static.py` says 5.
  - `euler_mapgrid_2D_constants.num_waves` is 3; `static.py` says 4.
  - `reactive_euler_with_efix_1D_constants.num_waves` is 3, but `static.py`
    says 4. The solver computes 3 waves yet sums the fluctuations over
    `mwaves`, so a caller using 4 reads a 4th wave and speed that the solver
    never set.
  - `static.fwave` has no entry for `sw_aug_1D` (PyClaw would treat it as a
    regular wave solver), and its `sw_aug_2D` entry is `3`, not `True`.
- **`shallow_roe_with_efix` 1D vs 2D entropy fix.** `rpn2` skips the 3-wave
  (u+c) transonic check whenever the shear wave is right-going
  (`if (s(2,i) > 0) goto 200`), while `rp1` always performs it. The two
  solvers give different fluctuations at interfaces with a transonic u+c
  rarefaction and positive Roe-averaged velocity. Both still conserve.
  Which behavior is intended is a question for the maintainers, not for
  this test.
- **Euler HLLE sound speed.** `rp1_euler_hlle.f90`, `rpn2_euler_hlle.f90`
  and `rpn2_euler_hlle_with_walls.f90` compute the left and right sound
  speeds as `sqrt(gamma1*H - u**2/2)` instead of
  `sqrt(gamma1*(H - u**2/2))`; the 2D kinetic term also omits v². In
  addition, `H_r` uses `p_l`. The Roe-averaged `c_hat` line is correct. As a
  result:
  - the HLLE speeds are underestimated (`test_hlle_speeds_match_einfeldt`,
    strict, fails on every platform);
  - when the computed c² is negative, the square root is NaN. Whether that
    NaN reaches the output depends on how the processor evaluates
    `MIN`/`MAX`: on x86-64 (CI) it does, on arm64 it did not. The
    conservation, wave-sum, rotation and 1D/2D tests for these solvers
    therefore carry a **non-strict** mark, and report XPASS where NaN
    doesn't show through.
