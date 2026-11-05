"""Call compiled Clawpack Riemann solvers directly from Python.

The meta-repository's meson build compiles each solver in ``src/`` into an
f2py extension module (``clawpack.riemann.<name>``).  These wrappers do not all
share one Python signature: argument names differ (``maxm``/``maxmx``,
``mbc``/``num_ghost``, ...), some return ``wave, s, amdq, apdq`` while others
fill caller-supplied arrays in place, and a few accept ``auxl``/``auxr`` as
scalars.  :func:`solve` reads each wrapper's f2py docstring and adapts, so the
tests can treat every solver the same way.

Index convention (as in Clawpack's ``step1``/``flux2``): with ``n`` cells
including ghost cells, the solver computes interface ``j`` (``1 <= j < n``)
between the left state ``qr[:, j-1]`` and the right state ``ql[:, j]``, and
stores the result at index ``j``.
"""

import importlib
import re

import numpy as np

NUM_GHOST = 2
NUM_CELLS = 40
NUM_TOTAL = NUM_CELLS + 2 * NUM_GHOST

# f2py argument name -> canonical name used here.
_ALIASES = {
    "maxm": "maxm", "maxmx": "maxm",
    "mx": "mx", "num_cells": "mx",
    "meqn": "meqn", "num_eqn": "meqn",
    "mwaves": "mwaves", "num_waves": "mwaves",
    "maux": "maux", "num_aux": "maux",
    "mbc": "mbc", "num_ghost": "mbc",
    "ixy": "ixy", "ixyz": "ixy",
    "ql": "ql", "qr": "qr", "auxl": "auxl", "auxr": "auxr",
    "wave": "wave", "fwave": "wave", "s": "s",
    "amdq": "amdq", "apdq": "apdq",
}

_OUTPUTS = ("wave", "s", "amdq", "apdq")


def load(name):
    """Import ``clawpack.riemann.<name>``."""
    return importlib.import_module("clawpack.riemann." + name)


def signature(func):
    """Parse an f2py docstring.

    Returns ``(returned, args, scalar_args)``: the names of returned values
    (empty for in-place wrappers), all argument names in order (required then
    optional), and the set of arguments f2py wraps as scalars.
    """
    doc = func.__doc__
    first = doc.splitlines()[0]
    returned = []
    if "=" in first:
        lhs, first = first.split("=", 1)
        returned = [a.strip() for a in lhs.split(",")]
    inner = first[first.index("(") + 1:first.rindex(")")]
    required, _, optional = inner.partition("[")
    args = [a.strip() for a in (required + "," + optional.rstrip("]")).split(",")
            if a.strip()]
    scalar_args = set(re.findall(r"^(\w+) : input (?:float|int)\s*$", doc, re.M))
    return returned, args, scalar_args


def aux_is_array(func):
    """True if the wrapper takes ``auxl``/``auxr`` as arrays (see README)."""
    _, args, scalar_args = signature(func)
    return not ({"auxl", "auxr"} & scalar_args)


def solve(module, fname, ql, qr, mwaves, aux=None, ixy=None):
    """Call ``module.<fname>`` (rp1/rpn2/rpn3) on cell states ``ql``/``qr``.

    ``ql``, ``qr``: arrays of shape ``(meqn, NUM_TOTAL)``.  ``aux``: optional
    array of shape ``(maux, NUM_TOTAL)`` used for both ``auxl`` and ``auxr``.
    Returns a dict with ``wave`` (meqn, mwaves, n), ``s`` (mwaves, n),
    ``amdq`` and ``apdq`` (meqn, n).
    """
    func = getattr(module, fname)
    returned, args, scalar_args = signature(func)
    meqn, n = ql.shape
    if aux is None:
        aux = np.zeros((1, n))
    maux = aux.shape[0]
    values = {
        "maxm": n - 2 * NUM_GHOST, "mx": n - 2 * NUM_GHOST,
        "meqn": meqn, "mwaves": mwaves, "maux": maux, "mbc": NUM_GHOST,
        "ixy": ixy,
        "ql": np.asfortranarray(ql, dtype=float),
        "qr": np.asfortranarray(qr, dtype=float),
        "auxl": np.asfortranarray(aux, dtype=float),
        "auxr": np.asfortranarray(aux, dtype=float),
        "wave": np.zeros((meqn, mwaves, n), order="F"),
        "s": np.zeros((mwaves, n), order="F"),
        "amdq": np.zeros((meqn, n), order="F"),
        "apdq": np.zeros((meqn, n), order="F"),
    }
    kwargs = {}
    for arg in args:
        key = _ALIASES.get(arg)
        if key is None:
            raise TypeError(f"{module.__name__}.{fname}: unknown argument {arg!r}")
        if returned and key in _OUTPUTS:
            continue
        value = values[key]
        if arg in scalar_args and isinstance(value, np.ndarray):
            value = 0.0         # wrapper declares this array as a scalar
        kwargs[arg] = value
    result = func(**kwargs)
    if returned:
        out = dict(zip((_ALIASES[r] for r in returned), result))
    else:
        out = {k: values[k] for k in _OUTPUTS}
    return out


def left_right(q):
    """Left and right states at interfaces 1..n-1 for cell states ``q``."""
    return q[..., :-1], q[..., 1:]


def at_interfaces(a):
    """Restrict a solver output to interfaces 1..n-1 (last axis)."""
    return a[..., 1:]


def assert_close(actual, expected, tol, what):
    """Absolute check scaled by the size of the data being compared."""
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    scale = max(1.0, float(np.max(np.abs(expected))) if expected.size else 1.0)
    err = float(np.max(np.abs(actual - expected))) if actual.size else 0.0
    assert np.all(np.isfinite(actual)), f"{what}: non-finite values"
    assert err <= tol * scale, f"{what}: max error {err:.3e} > {tol:.1e} * {scale:.3g}"
