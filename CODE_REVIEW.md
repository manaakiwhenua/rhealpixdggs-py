# Code Review — rHEALPixDGGS-py

**Date:** 2026-08-07
**Reviewed at:** commit `e10f3c2` (`v0.6.0`, branch `master`)
**Scope:** whole repository — correctness/robustness, test coverage, packaging/licensing, docs, CI/tooling, repo hygiene
**Method:** static reading of all 9 modules in `rhealpixdggs/`, cross-referencing against the project's own test suite, `CHANGES.rst`, and (at the maintainer's request, since correctness is the top priority for what is understood to be the only implementation of this DGGS) the founding publication — R. G. Gibb, *"The rHEALPix Discrete Global Grid System,"* IOP Conf. Ser.: Earth Environ. Sci. **34** 012012 (2016) — and its companion conference slide deck. A handful of findings below were additionally confirmed by executing the code (marked "confirmed by execution"); everything else is confirmed by reading, cited to file:line.

## Executive summary

This is a reasonably-maintained, well-tested-for-its-size library (77 runnable unit tests + 392 doctest examples, all passing; a clean, acyclic module dependency graph) with no CI, no linting, and a handful of confirmed correctness bugs — several only visible once you step outside the exact inputs the test suite happens to exercise (a specific `n`, a specific `N_side`, a specific cell shape). None of what follows suggests the core geometry is fundamentally unsound; most of the confirmed bugs are narrow, mechanical, and cheap to fix. The one finding worth real attention is structural rather than a bug: **the founding paper's central technical contribution — that rHEALPix cell IDs support the OGC-mandated DE-9IM topological predicates directly and efficiently — isn't actually implemented or exposed anywhere in the library.** That's a gap between the published compliance claim and the shipped API, not a bug in existing code, but it's the one item here that bears on whether the library fully delivers what the paper says it delivers.

**If you only fix five things:**
1. `n == 2` → `n = 2` in `Cell.boundary()` (cell.py:813) — a one-character typo that currently crashes with `ZeroDivisionError` on `boundary(n=1, plane=True)` for any non-quad/cap cell (**confirmed by execution**).
2. `Cell.boundary(plane=False)`'s 0.6.0 short-circuit silently returns 4 points instead of the documented `4n-4` for quad/cap cells when `n>2` (**confirmed by execution**) — either fix the short-circuit to honor `n`, or fix the docstring/changelog to disclose the limitation.
3. `except NameError` → `except (AttributeError, ModuleNotFoundError)` in `Projection.__call__` (projection_wrapper.py:103) — currently any of the three declared-but-unimplemented "homemade" projections (`isea`, `csea`, `qsc`) crash uncaught instead of degrading gracefully (**confirmed by execution**).
4. `compact_cells`/`compress_order_cells` hardcode "9 siblings = compress" (conversion.py:151,196) — silently never compresses anything for the shipped `WGS84_002` (`N_side=2`) preset, where a complete sibling group has 4 members, not 9.
5. Decide what to do about the missing DE-9IM predicates (see SPEC-1) — even thin, well-tested wrappers for `touches`/`disjoint`/`within`/`covers` around the existing common-parent/relative-child logic would close the gap the paper claims is closed.

## Findings dashboard

| ID | Finding | Tier | Category | Confidence | Fix effort |
|---|---|---|---|---|---|
| BUG-4 | `n == 2` typo in `boundary()` — guard is a no-op, crashes on `n<2` | S1 | Correctness | Confirmed (executed) | S |
| BUG-2 | `boundary(plane=False)` short-circuit breaks documented `4n-4` contract for quad/cap, `n>2` | S1 | Correctness | Confirmed (executed) | S–M |
| BUG-1 | `except NameError` can never catch the real errors in `Projection.__call__` | S1 | Correctness | Confirmed (executed) | S |
| BUG-5 | `compact_cells`/`compress_order_cells` hardcode a 9-sibling assumption; silently inert for `N_side=2` | S1 | Correctness | Confirmed (read + preset check) | S |
| BUG-3 | `print()`+`inf` sentinel contradicts "throw error" comment; no exception a caller can catch | S1–S2 | Correctness/Robustness | High (read) | M |
| SPEC-1 | Paper's DE-9IM predicates (touches/disjoint/within/covers/…) not implemented anywhere; `overlaps()` mislabeled per the paper's own terms | S2/S4 | Spec-conformance | Confirmed (paper + grep) | M–L |
| LAT-1 | `neighbors()` mutates shared singleton `ellipsoid.lon_0`, no try/finally | S2 | Robustness/Latent | High | S–M |
| LAT-2 | Input validation via bare `assert` (5 sites) — stripped under `-O`, wrong exception type | S2/S4 | Robustness | High | S |
| TEST-5 | `tests/test_PJ_healpix.c.py` (4 tests) silently never runs via `run_unittests.sh` | S4 | Test-coverage | Confirmed (executed) | S |
| INC-1 | `_find_cell_ring_start` self-acknowledged wrong for cap faces | S3 | Known-incomplete | High (author-confirmed) | L |
| INC-2 | Cap-cell handling incomplete in `cells_from_line`; antimeridian TODOs in `polyfill`/`linetrace` | S3 | Known-incomplete | High (author-confirmed) | L |
| TEST-1 | `auth_lat()`'s large-flattening branch (f>1/150) has zero coverage, unreachable by any predefined `Ellipsoid` | S4 | Test-coverage | High | S |
| TEST-2 | `projection_wrapper.py` has zero unittest coverage | S4 | Test-coverage | High | M |
| TEST-3 | `Cell.area/color/overlaps/region_overlaps` untested; `area_error_budget()` doctest-only; `test_centroid` itself is a disabled stub | S4 | Test-coverage | High | M |
| TEST-4 | `boundary()`'s 0.6.0 short-circuit has no dedicated regression test | S4 | Test-coverage | High | S |
| DOC-1 | Sphinx toctree omits `cell.py`, `conversion.py`, `rhp_wrappers.py` | S4 | Documentation | Confirmed | S |
| DOC-2 | README references nonexistent `requirements*.txt`; states `numpy<2` vs. `pyproject.toml`'s `numpy>=2.0` | S4 | Documentation | Confirmed | S |
| DOC-3 | `docs/source/conf.py` stale (`release="0.5.4"`, `copyright="2013, Alexander Raichev"`) | S4 | Documentation | Confirmed | S |
| PKG-1 | `pyproject.toml` license field (`GPL-3.0-or-later`) contradicts actual dual LGPL-3.0/MIT | S4 (elevated) | Packaging | Confirmed | S |
| PKG-2 | `poetry.lock` untracked by git | S4 | Packaging | Confirmed | S |
| PKG-3 | `pypi.txt` describes an outdated `setup.py`/`twine` release flow | S4 | Packaging | Confirmed | S |
| PROC-1 | No CI anywhere — the 77-test/392-doctest suite guards nothing on incoming PRs | S4 (elevated) | Process | Confirmed | M |
| PROC-2 | No lint/type-check enforcement (`black` declared, never run) | S4 | Process | Confirmed | M |
| PROC-3 | Classifiers list only Python 3.11; no CI evidence 3.12/3.13 work | S4 | Process | Confirmed | S |
| DEAD-1 | `RhealPolygon` stub class, zero references/tests | S4/S0 | Dead-code | Confirmed | S |
| DEAD-2 | Commented-out "Option 2" block in `centroid()` | S0 | Dead-code | Confirmed | S |
| DEAD-3 | ~150 lines of commented-out tests for a `CellFamily` class that doesn't exist anywhere in current source | S0 | Dead-code | Confirmed | S |
| DEAD-4 | Tracked build cruft: `.coverage`, `.bak`, `docs/.doctrees/*`, orphaned `.txt` fixtures | S0 | Dead-code/Hygiene | Confirmed | S |

Severity tiers: **S0** fix now, mechanical, zero risk · **S1** correctness bug, live · **S2** latent/conditional bug · **S3** known limitation, author-acknowledged · **S4** hygiene/coverage/process gap · **S5** strength, no action.

---

## Detailed findings

### BUG-4 — `boundary()`'s low-`n` guard is a no-op (`==` instead of `=`)

- **Location:** `rhealpixdggs/cell.py:812-813`
- **Tier / Category / Confidence / Effort:** S1 / Correctness / Confirmed by execution / S

```python
if n < 2:
    n == 2
```

That's a comparison expression, not an assignment — its result is discarded, so `n` is never actually clamped to 2. Confirmed:

```pycon
>>> c = Cell(RHEALPixDGGS(), ['N', 2])   # a "dart" cell — doesn't hit the quad/cap short-circuit
>>> c.ellipsoidal_shape
'dart'
>>> c.boundary(n=1, plane=True)
ZeroDivisionError: float division by zero        # at cell.py:818: delta = (w - 2*eps) / (n - 1)
>>> c.boundary(n=0, plane=True)
[]                                                 # silently empty instead of clamped to n=2
```

**Why it matters:** any caller that passes `n<2` (a plausible mistake, or a caller iterating `n` from 0) gets a crash or a silently-empty boundary instead of the documented minimum. This is masked for quad/cap cells because they now short-circuit to `vertices()` before reaching this code at all (see BUG-2) — the bug only fires for dart/skew_quad cells, which is presumably why the existing test suite (which happens to always call `boundary()` with `n>=2`) never caught it.

**Suggested fix:** `n = 2`. One character. Add a regression test for `boundary(n=0)`/`boundary(n=1)` on a dart cell while you're there (feeds TEST-4).

---

### BUG-2 — `boundary(plane=False)`'s quad/cap short-circuit silently ignores `n`

- **Location:** `rhealpixdggs/cell.py:806-809` (short-circuit); docstring contract at `cell.py:756-759`
- **Tier / Category / Confidence / Effort:** S1 / Correctness / Confirmed by execution / S–M

The docstring's first line is unconditional: *"Return a list of `4*n - 4` boundary points of this cell, `n` on each edge."* The 0.6.0 performance change added:

```python
if not plane and self.ellipsoidal_shape in ("quad", "cap"):
    return self.vertices(plane=False)
```

`vertices()` always returns exactly 4 points, regardless of `n`. Confirmed by execution on `Cell(RHEALPixDGGS(), ['P', 2])` (shape `"quad"`):

| `n` | documented `4n-4` | actual |
|---|---|---|
| 2 | 4 | 4 |
| 3 | 8 | **4** |
| 5 | 16 | **4** |

The existing doctest for this method (cell.py:774-796) only exercises a `"dart"`-shaped cell, which doesn't take the short-circuit — so the discrepancy is invisible to both the doctest and the unittest suite. The docstring does have a later paragraph (cell.py:768-772) describing this exact optimization as something *"callers can"* opt into via `ellipsoidal_shape` — but the primary contract statement at the top was never updated to reflect that it's now applied unconditionally inside `boundary()` itself, and there's no runtime warning.

**Why it matters:** any caller relying on `n` points per edge for `plane=False` — e.g. building a smoother polyline for rendering/simplification at higher resolution — silently gets 4 points instead, with no signal that anything changed. This is exactly the kind of thing CHANGES.rst 0.6.0 called "performance improvements," which is true, but the contract change deserves equal billing.

**Suggested fix:** either honor `n` for quad/cap cells too (interpolate along the straight/rotationally-symmetric edges implied by the 4 vertices — cheap, since these shapes are geometrically simple), or update the docstring's opening line and CHANGES.rst to state the `n=4` fixed-count exception explicitly for `plane=False` on quad/cap cells, and add a regression test either way (TEST-4).

---

### BUG-1 — `except NameError` can never catch the errors it's meant to catch

- **Location:** `rhealpixdggs/projection_wrapper.py:98-105`
- **Tier / Category / Confidence / Effort:** S1 / Correctness / Confirmed by execution / S

```python
HOMEMADE_PROJECTIONS = {"healpix", "rhealpix", "isea", "csea", "qsc"}
...
try:
    module = importlib.import_module("rhealpixdggs.pj_" + proj)
    f = getattr(module, proj)(a=a, e=e, **kwargs)
except NameError:
    print("Oops! Projection %s is not implemented." % proj)
    return
```

Three of the five declared `HOMEMADE_PROJECTIONS` — `isea`, `csea`, `qsc` — have no backing `rhealpixdggs/pj_*.py` module in the package at all. Confirmed by execution:

```pycon
>>> Projection(ellipsoid=WGS84_ELLIPSOID, proj="isea")(0, 0)
ModuleNotFoundError: No module named 'rhealpixdggs.pj_isea'
```

`importlib.import_module` raises `ModuleNotFoundError`; a present-module-but-missing-function case would raise `AttributeError` from `getattr`. Neither is `NameError` — this `except` clause reads like a leftover from a pre-`importlib` implementation (e.g. one that used `eval(proj)`, which *would* raise `NameError`). The intended graceful-degradation message is dead code for its stated purpose.

**Suggested fix:** `except (AttributeError, ModuleNotFoundError):`. Also worth a one-line decision: either implement `isea`/`csea`/`qsc`, or drop them from `HOMEMADE_PROJECTIONS` until they exist, so the set doesn't advertise unimplemented projections.

---

### BUG-5 — `compact_cells`/`compress_order_cells` hardcode a 9-sibling assumption

- **Location:** `rhealpixdggs/conversion.py:145-166` (`compress_order_cells`), `:169-204` (`compact_cells`)
- **Tier / Category / Confidence / Effort:** S1 / Correctness / Confirmed by reading + preset check / S

Both functions group cells by parent prefix (`cell[:-1]`) and merge a group into its parent only `if len(v) == 9` / `if len(children) == 9`. That's correct only for an `N_side=3` DGGS (9 children per cell, digits `0`-`8`). But `N_side` is a general, user-facing parameter (`RHEALPixDGGS.__init__`'s own docstring: *"An integer of size at least 2"*, cell.py:182), and the package **ships a preset with `N_side=2`**:

```python
# rhealpixdggs/dggs.py:1539-1541
WGS84_002 = RHEALPixDGGS(ellipsoid=WGS84_ELLIPSOID, north_square=0, south_square=0, N_side=2)
```

For `WGS84_002`, a cell has only 4 children (digits `0`-`3`), so a complete sibling group has at most 4 members — `len(v) == 9` is never true, and compaction silently never happens. No error, no warning; `compact_cells`/`compress_order_cells` just quietly do nothing for any DGGS other than the default `N_side=3` ones. `polyfill()` (`rhp_wrappers.py`) calls `compact_cells` internally (per CHANGES.rst 0.5.14), so anyone calling `polyfill(..., dggs=WGS84_002)` gets uncompressed, far more verbose output than the function's own docstring promises, with no indication why.

**Suggested fix:** replace the hardcoded `9` with `cell.rdggs.N_side ** 2` (or accept `N_side`/a sample cell as a parameter, since these are currently free functions operating on plain ID strings with no DGGS reference at all — that's also why they can't detect this themselves). Add a test using `WGS84_002` or another non-default `N_side`.

---

### SPEC-1 — The paper's DE-9IM predicates aren't implemented or exposed

- **Location:** absence, across `rhealpixdggs/cell.py`; contrast with Gibb (2016) Section 4 and the accompanying slide deck (slides 21-23)
- **Tier / Category / Confidence / Effort:** S2/S4 / Spec-conformance / Confirmed (paper cross-check + repo-wide grep) / M–L

The founding paper's central technical contribution is showing that cell IDs alone (via a shared-prefix "`ID_common_parent`"/"`ID_relative_child`" decomposition) can compute the full OGC-mandated DE-9IM predicate set — `equals`, `within`, `contains`, `touches` (with explicit regex patterns distinguishing edge- vs. corner-touching relative-child suffixes), `disjoint`, `covers`, `coveredBy`, and the paper explicitly notes `intersects`/`crosses`/`overlaps` are structurally impossible between two grid cells (Section 4(vii)): *"We note that because of the way cells are defined, no single cell can intersect, cross, or overlap another cell."*

A repo-wide search confirms none of these are implemented as named methods anywhere:

```
$ grep -rniE "\bdef (touches|disjoint|covers|coveredby|crosses|intersects|equals|within)\b" rhealpixdggs/ tests/
(no matches)
$ grep -rniE "\btouches\b|\bcoveredby\b|\bDE-?9IM\b" rhealpixdggs/ README.rst CHANGES.rst
(no matches)
```

The only cell-to-cell relational method is `Cell.overlaps()` (cell.py:970-983), which does a `zip(self.suid, other_cell.suid)` prefix comparison — i.e., it computes exactly the paper's "shares an `ID_common_parent`" ancestor/descendant test (paper §4(ii): *within*/*contains*), not overlap. Per the paper's own §4(vii) quoted above, calling this `overlaps()` is arguably a direct naming clash with the paper's own vocabulary, where true cell-cell "overlap" is stated to be impossible by construction.

**Why it matters:** this isn't a bug in `overlaps()` — it does what it does correctly and is tested (`test_cell.py`). But a maintainer, downstream integrator, or OGC-conformance auditor reading the paper this library implements would reasonably expect to find `touches`/`disjoint`/`within`/`covers` as first-class, documented cell methods, since the paper presents their cell-ID-based computation as the whole point. They aren't there. This is the one finding in this review that speaks to whether the library fully delivers the compliance story its own paper tells, rather than to a bug in delivered code.

**Suggested fix:** the underlying primitive (common-parent / relative-child string decomposition) is already half-built into `overlaps()`. Adding `within`/`contains` (as thin, correctly-named wrappers, since that's what `overlaps()` already computes), `equals`, and `touches` (using the paper's own regex patterns on the relative-child suffix — `(0|1|2)+|(2|5|8)+|(0|3|6)+|(6|7|8)+` for edge, `(0)+|(2)+|(6)+|(8)+` for corner, for `N_side=3`; would need generalizing for other `N_side`) would close most of the gap cheaply. `disjoint`/`covers`/`coveredBy` follow directly per the paper's §4(v)-(vi) as boolean combinations of the above. Consider renaming or deprecating `overlaps()`/`region_overlaps()` in favor of correctly-named `within`/`contains` methods, with a changelog entry consistent with the project's existing breaking-change-disclosure practice.

---

### BUG-3 — `print()` + `inf` sentinel instead of raising

- **Location:** `rhealpixdggs/pj_healpix.py:79-82`; mirrored in `rhealpixdggs/pj_rhealpix.py` (~327, ~403)
- **Tier / Category / Confidence / Effort:** S1–S2 / Correctness/Robustness / High (read) / M

```python
# Throw error if input coordinates are out of bounds.
if not in_healpix_image(x, y):
    print("Error (hsi): input coordinates (%.20f,%.20f) are out of bounds" % (x, y))
    return float("inf"), float("inf")
```

The comment says "throw error"; the code does neither — it prints to stdout (which a library should generally never do on a caller's behalf) and returns an `inf` sentinel with no exception a caller can catch. `inf` can then silently propagate into downstream area/containment arithmetic without ever raising, in the innermost numeric layer of the whole projection stack.

**Suggested fix:** raise a `ValueError` (or a small custom exception) with the same message instead of printing + returning a sentinel. Consider bundling this with LAT-2 and BUG-1 as one design-level change: adopt raising typed exceptions consistently at API boundaries, rather than three different informal "something's wrong" signals (print+sentinel here, bare `assert` in `cell.py`, wrong exception type in `projection_wrapper.py`).

---

### LAT-1 — `neighbors()` mutates a shared, process-wide singleton with no exception safety

- **Location:** `rhealpixdggs/cell.py:1408-1457` (dart & skew_quad branches); singletons at `dggs.py:1539-1549`
- **Tier / Category / Confidence / Effort:** S2 / Robustness/Latent / High / S–M

```python
old_lon_0 = self.rdggs.ellipsoid.lon_0
self.rdggs.ellipsoid.lon_0 = -self.nucleus(plane=False)[0]
...
self.rdggs.ellipsoid.lon_0 = old_lon_0
```

`self.rdggs.ellipsoid` is not a private, per-call object — for anyone using the shipped presets (`WGS84_002`, `WGS84_003`, `UNIT_003`, etc., `dggs.py:1539-1549`), it's a single module-level instance shared for the lifetime of the process. This block temporarily mutates `lon_0` in place, with **no `try/finally`** around the restore. `nucleus(plane=False)` — called *inside* that exact mutated-state window — routes through the projection stack (`pyproj`/homemade projections), which is a real, non-exotic path for an exception to fire (domain errors, `pyproj` errors) between the mutate and the restore. If that happens, `lon_0` is left permanently wrong on the shared ellipsoid for every subsequent use of that DGGS preset, for the rest of the process — silent, hard-to-diagnose corruption, not a crash.

Secondary, currently-inactive concern: this pattern is also not thread-safe. A repo-wide search found no `threading`/`asyncio`/`multiprocessing` use anywhere in the codebase or its tests, so this isn't being exercised today — but it's worth knowing about before this library gets used from, say, a web service handling concurrent requests against a shared preset instance.

**Suggested fix:** wrap the mutate/restore in `try`/`finally`. Longer-term, consider whether `neighbors()` needs to mutate the shared ellipsoid at all, versus computing the rotated/re-centered coordinates directly without touching shared state.

---

### LAT-2 — Input validation via bare `assert`

- **Location:** `rhealpixdggs/cell.py:154,157,160,166` (constructor), `:979` (`overlaps()`)
- **Tier / Category / Confidence / Effort:** S2/S4 / Robustness / High / S

Five sites validate constructor/method arguments with `assert`, e.g.:

```python
assert isinstance(suid, list) or isinstance(suid, tuple), (
    "Cell suid must be a list or tuple. Got %s." % suid
)
```

`assert` statements are compiled out entirely under `python -O`/`PYTHONOPTIMIZE`, silently disabling all of this validation. They also raise `AssertionError` for what's really API misuse, which conventionally should be `ValueError`/`TypeError`. Separately, `cell.py:979`'s `assert self.suid is not tuple()` is an identity check against a freshly-constructed literal — it happens to work because CPython interns the empty tuple, but that's an implementation detail, not a language guarantee, which matters given the project's own classifiers claim PyPy support.

**Suggested fix:** replace with explicit `if ...: raise ValueError(...)`/`TypeError(...)`.

---

### TEST-5 — A whole test file silently never runs

- **Location:** `tests/test_PJ_healpix.c.py`
- **Tier / Category / Confidence / Effort:** S4 / Test-coverage / Confirmed by execution / S

`run_unittests.sh` runs `python -m unittest discover tests`. Verbose discovery output was diffed against a full accounting of every `def test_` in the tree:

```
$ python3 -m unittest discover tests -v 2>&1 | grep -c "^test_"
77
```

77 tests actually run: `test_cell`(26) + `test_conversion`(14) + `test_dggs`(9 — the other 9 `def test_` lines in that file are commented out, see DEAD-3) + `test_ellipsoids`(1) + `test_pj_healpix`(3) + `test_pj_rhealpix`(5) + `test_rhp_wrappers`(13) + `test_utils`(6). **`tests/test_PJ_healpix.c.py`'s 4 tests (which compare the homemade HEALPix implementation against `pyproj`'s own C-library reference — a genuinely valuable cross-check) never appear at all**, confirmed both by their absence from `-v` output and by `python -m unittest tests.test_PJ_healpix.c` failing with `ModuleNotFoundError: No module named 'tests.test_PJ_healpix'`. The embedded second "extension" in the filename (`.c.py`) breaks `unittest`'s dotted-module-name resolution.

Running the file directly (`python3 tests/test_PJ_healpix.c.py`) does work, but produces plain stdout comparison output rather than pass/fail assertions in the normal sense — it looks designed to be run manually and its output captured (matching the adjacent, apparently-orphaned `tests/test_PJ_healpix.c.txt`/`c_degrees.txt` reference-output files, see DEAD-4). Whether this was a deliberate manual-only test or an accidental casualty of the filename, `run_unittests.sh` reporting "OK" gives no indication these 4 tests were skipped.

**Suggested fix:** rename to something `unittest discover` can resolve (e.g. `test_pj_healpix_c_reference.py`), or explicitly document that it's a manual-only comparison script excluded from `run_unittests.sh` on purpose.

---

### INC-1 / INC-2 — Known, self-acknowledged incompleteness (aggregated, not newly discovered)

- **INC-1 location:** `rhealpixdggs/rhp_wrappers.py:753-822` (`_find_cell_ring_start`, TODO at line 814); `tests/test_rhp_wrappers.py:634-648` (matching `@unittest.expectedFailure`)
- **INC-2 location:** `rhealpixdggs/dggs.py:1173-1178` (`cells_from_line`); `rhp_wrappers.py:532` (`polyfill`), `:639` (`linetrace`)
- **Tier / Category / Confidence:** S3 / Known-incomplete / High (author-confirmed)

These are already flagged by the code itself — a live `# TODO: this is wrong` in the ring-walking helper, a matching `expectedFailure` test documenting the actual-vs-expected wrong output for a cap-face case, explicit TODOs about cap-cell handling in `cells_from_line`, and two "decide what to do about the antimeridian" TODOs in `polyfill`/`linetrace`. `CHANGES.rst`'s 0.5.6 entry already discloses the linetrace limitation publicly. This review's contribution is aggregating them in one place and assessing scope, not discovering them.

Worth noting as a **positive**, not a gap: `cell_ring`/`k_ring` already warn about this by default —

```python
# rhp_wrappers.py:31
CELL_RING_WARNING = "WARNING: Implementation of cell rings is incomplete. Requesting a {0} ring that involves more than two resolution 0 cube faces will return unexpected results."
# rhp_wrappers.py:395-396, 490-491
if verbose:
    warn(str.format(CELL_RING_WARNING, "cell"))
```

— consistently with `linetrace`'s own `LINETRACE_WARNING`. This is a good, already-applied pattern; no action needed here.

**Suggested fix (aggregate):** these are genuinely hard geometry problems (cap-face/cube-corner handling, antimeridian crossing) and not a quick fix. Given the paper's own row/column base-N arithmetic (Section 4) already handles ordinary face-to-face adjacency cleanly and is faithfully reflected in `Cell.neighbor()`'s "naive neighbor + rotation correction" implementation (cell.py:1290-1347), it may be worth revisiting whether `_find_cell_ring_start`'s hand-rolled state machine could be simplified by walking the ring one `neighbor()` call at a time (already correctness-tested, 26 tests) rather than trying to jump directly to a ring's start cell — trading some performance for correctness/simplicity, at least until the underlying geometry problem is solved properly.

---

## Test-coverage findings (TEST-1 – TEST-4)

- **TEST-1** — `auth_lat()`'s large-flattening direct-formula branch (`f > 1/150`, `rhealpixdggs/utils.py`, ~lines 116-343) has zero test coverage anywhere, and **no predefined `Ellipsoid` in the codebase can reach it**: WGS84 (f≈0.00335) and its radians variant fall well under the 1/150≈0.00667 threshold and exercise the power-series branch instead; the sphere presets (`UNIT_SPHERE`, etc.) short-circuit via `e==0` before either branch. Since `Ellipsoid` is public, user-constructible API (e.g. someone modeling a fast-rotating or non-Earth body), this is a real, if narrow, gap. The power-series branch itself is well-validated — its docstring cites specific source-paper equations and `test_utils.py` round-trips both directions in both radians and degrees against arbitrary-precision reference values — so this finding is scoped specifically to the untested branch, not the whole function.
- **TEST-2** — `rhealpixdggs/projection_wrapper.py` has no dedicated unittest file; its only coverage is 5 doctest examples.
- **TEST-3** — `Cell.area()`, `color()`, `overlaps()`, `region_overlaps()` have neither a unittest nor a doctest. `RHEALPixDGGS.area_error_budget()` (a prominently-documented 0.5.16 feature, featured in `README.rst`'s "Error Budget" table) has a doctest but no dedicated unittest. Found while fixing DEAD-3: `tests/test_cell.py::test_centroid` (lines 549-626) is itself a disabled stub — the method body is `pass` followed by ~75 lines of commented-out, Python-2-era Monte Carlo comparison code (`print` statement syntax, an unimported `uniform`) that hasn't run in years. `Cell.centroid()` — which uses `scipy.integrate.dblquad` for two of its four cell-shape branches — currently has **no active unittest coverage at all**, only the one doctest example in its docstring.
- **TEST-4** — `boundary()`'s 0.6.0 short-circuit (BUG-2) has no dedicated regression test; add one alongside the fix.

## Documentation findings (DOC-1 – DOC-3)

- **DOC-1** — `docs/source/index.rst`'s toctree lists only `introduction, utils, pj_healpix, pj_rhealpix, ellipsoids, projection_wrapper, dggs` — **`cell.py`, `conversion.py`, and `rhp_wrappers.py` (arguably the most-used modules, including the entire H3-style public API) have no Sphinx page and never reach the built manual.** The `index.rst` header comment (`.. scenzgrid-dggs documentation master file.`) is also a pre-rename relic from the project's original name.
- **DOC-2** — `README.rst:126,132` reference `requirements.txt`/`requirements.dev.txt`, which don't exist (the project moved to Poetry). `README.rst`'s stated `NumPy >=1.25.2,<2` also directly contradicts `pyproject.toml:16`'s `numpy>=2.0`.
- **DOC-3** — `docs/source/conf.py:54,63` hardcode `copyright = "2013, Alexander Raichev"` and `release = "0.5.4"` — both several years and several major versions stale. Worth contrasting with `CITATION.CFF`, which **is** accurate and current (0.6.0, correct maintainer, dated 2026-05-05) — the staleness is uneven, not universal.

## Packaging findings (PKG-1 – PKG-3)

- **PKG-1** — `pyproject.toml:13` declares `license = "GPL-3.0-or-later"`; the repo is actually dual LGPL-3.0-or-later/MIT (`README.rst:203-212`, `LICENSE`, `LICENSE-LGPL`, `LICENSE-MIT`) — a loose end from the dual-licensing PR (#39, `c38eae1`) that updated the license files and README but not this field (and left the old GPLv3 classifier commented out at line 34 rather than replaced). Worth elevating slightly above a typical hygiene item: an incorrect SPDX-style license field can cause automated corporate dependency scanners to wrongly block adoption of a project that's actually **more** permissively licensed than declared.
- **PKG-2** — `poetry.lock` exists locally but is untracked by git (your own `git status` shows `?? poetry.lock`) — it's not gitignored, it's simply never been committed. Scope note: this affects contributor/CI build reproducibility only, not end-user `pip install` (wheel/sdist builds driven by `pyproject.toml` don't consume the lock file).
- **PKG-3** — `pypi.txt` (a manual release checklist) describes a `python setup.py sdist bdist_wheel` / `twine upload` flow; no `setup.py` exists anywhere in the repo, the actual release flow is `poetry build`/`poetry publish` per `README.rst:215-233`.

## Process/tooling findings (PROC-1 – PROC-3)

- **PROC-1** — No CI configuration anywhere (`.github/` doesn't exist; no `.travis.yml`/`.gitlab-ci.yml`/etc.). Worth elevating above a typical hygiene item for the same reason as PKG-1: the existing 77-test/392-doctest suite is genuinely decent for a library this size, but provides **zero regression protection on incoming PRs** — a bigger structural risk than the absence of linting.
- **PROC-2** — `black` is declared as a dev dependency (`pyproject.toml:50`) with no `[tool.black]` config section, no CI step, and no pre-commit hook — effectively decorative. No `mypy`/`ruff`/`flake8`/`pylint` config exists at all.
- **PROC-3** — `requires-python = ">=3.11"` with no upper bound is fine per current packaging norms (not a defect in itself), but classifiers list only `Python :: 3.11` and there's no CI to show 3.12/3.13 actually work.

## Dead code / repo hygiene (DEAD-1 – DEAD-4)

- **DEAD-1** — `RhealPolygon` (`rhealpixdggs/dggs.py`, ~line 1554) is a bare stub — `__init__` only, no other methods, zero references or tests anywhere in the codebase. Safe to remove, but note it's a public, importable class, so treat its removal as a (harmless) breaking change worth a changelog line, consistent with the project's existing practice of flagging breaking changes explicitly.
- **DEAD-2** — A ~90-line commented-out "Option 2" alternative implementation left inside `Cell.centroid()` (`cell.py:1073-1159`).
- **DEAD-3** — ~150 lines of commented-out `CellFamily` tests in `tests/test_dggs.py` (lines ~451-630) — for a class that doesn't exist anywhere in current source (no definition, no import, confirmed by a repo-wide search). The `WGS84_123` DGGS instance the block uses is still defined locally in the same file (line 44), so that part is fine — it's specifically `CellFamily` that's gone.
- **DEAD-4** — Tracked build artifacts that shouldn't be in version control: `rhealpixdggs/.coverage` (an 8KB coverage.py database, tracked despite matching the `.gitignore` pattern), `tests/test_PJ_healpix.c.py.bak` (a stray Python-2-era backup file, tracked), `docs/.doctrees/*` (11 tracked Sphinx binary pickle caches, including doctrees for pages that no longer exist under the project's original "scenzgrid" name, e.g. `scenzgrid_grids.doctree`). Additionally, `tests/test_PJ_healpix.c.txt`/`c_degrees.txt` (71KB/58KB) appear to be orphaned manual-run output with no code reading them (see TEST-5 for the likely explanation).

---

## Done well

- **Clean architecture.** A strictly layered module dependency graph with no import cycles: `utils.py` → `ellipsoids.py` → `pj_healpix.py` → `pj_rhealpix.py` → `projection_wrapper.py`, with `cell.py`/`dggs.py`/`conversion.py`/`rhp_wrappers.py` layered cleanly on top.
- **Test suite genuinely proportionate to project size.** 77 runnable unittest cases + 392 doctest examples, all passing at HEAD (confirmed by running `run_unittests.sh`/`run_doctests.sh`: `OK (expected failures=1)`, doctest exit code 0). Doctests use a defensive `my_round(...)` + `.tolist()` pattern specifically to keep float output reproducible across NumPy versions.
- **Faithful, and in places more complete, implementation of the paper's core adjacency scheme.** The paper (§4) sketches face-to-face adjacency via base-3/base-4 row/column arithmetic and states outright that it isn't giving "the complete solution" for cube-corner transitions ("this logic can be readily coded"). `Cell.neighbor()`/`neighbors()` (cell.py:1236-1458) implement this with an explicit rotation-correction step for exactly those transitions, and it's well covered by the test suite (26 `test_cell.py` cases) — this is the part of the spec the paper explicitly left as an exercise, and it's been done.
- **Existing incompleteness is honestly disclosed where it matters most.** `cell_ring`/`k_ring` warn by default (`CELL_RING_WARNING`) about exactly the geometry cases where they're known to be wrong, mirroring `linetrace`'s own warning — a good, consistently-applied pattern.
- **Changelog-based breaking-change disclosure already exists.** 0.6.0's `ellipsoidal_shape` change is explicitly flagged "Breaking change" in `CHANGES.rst` — worth extending the same scrutiny to changes described purely as "performance improvements" (see BUG-2).
- **`CITATION.CFF` is accurate and current** (version 0.6.0, correct current maintainer, dated) — a useful contrast to `docs/source/conf.py`'s staleness (DOC-3), showing metadata neglect is uneven rather than universal.
- **No copyleft-compatibility risk.** Runtime dependencies (numpy, scipy, matplotlib, pyproj, shapely) are entirely permissively (BSD/MIT-family) licensed — no license-compatibility issue for a dual LGPL/MIT project.
- **No meaningful security surface.** No `eval`/`exec`/`pickle`/`subprocess`/`os.system` anywhere in `rhealpixdggs/`. The one dynamic-dispatch site (`importlib.import_module("rhealpixdggs.pj_" + proj)` in `projection_wrapper.py`) is scoped to a hardcoded internal namespace prefix, not attacker-influenceable.
- **The 0.6.0 projection cache doesn't introduce a regression.** `Projection.__call__` reads `ellipsoid.lon_0`/`lat_0` live on every call rather than baking in values at construction time, so `RHEALPixDGGS._projection_cache` composes safely with the mutate/restore pattern in `neighbors()` (LAT-1's problem is the missing `try/finally`, not the cache).
- **SUID encode/decode is correctly factored, not duplicated.** `Cell.suid_from_index`/`index` and `RHEALPixDGGS.num_cells`/`child_order` share the same underlying combinatorics via `num_cells()`, and are round-trip tested (`test_cell.py::test_index`). Initially looked like possible duplication during this review; checked and it isn't.

---

## Appendix A — findings by file

| File | Findings |
|---|---|
| `rhealpixdggs/cell.py` | BUG-4, BUG-2, LAT-1, LAT-2, DEAD-2, TEST-3, TEST-4, SPEC-1 |
| `rhealpixdggs/projection_wrapper.py` | BUG-1, TEST-2 |
| `rhealpixdggs/conversion.py` | BUG-5 |
| `rhealpixdggs/pj_healpix.py` / `pj_rhealpix.py` | BUG-3 |
| `rhealpixdggs/dggs.py` | LAT-1 (singletons), DEAD-1, INC-2 |
| `rhealpixdggs/rhp_wrappers.py` | INC-1, INC-2 |
| `rhealpixdggs/utils.py` | TEST-1 |
| `tests/test_PJ_healpix.c.py` | TEST-5 |
| `tests/test_dggs.py` | DEAD-3 |
| `docs/source/index.rst`, `docs/source/conf.py` | DOC-1, DOC-3 |
| `README.rst` | DOC-2 |
| `pyproject.toml` | PKG-1, PROC-2, PROC-3 |
| repo root / `.github/` (absence) | PKG-2, PKG-3, PROC-1 |
| tracked build artifacts | DEAD-4 |

## Appendix B — quick-win checklist (mechanical, zero design judgment)

- [ ] Remove `RhealPolygon` stub class (`dggs.py` ~1554) — note as a breaking change in `CHANGES.rst` per existing practice (DEAD-1).
- [ ] Delete the commented-out "Option 2" block in `Cell.centroid()` (DEAD-2).
- [ ] Delete the commented-out `CellFamily` test block in `tests/test_dggs.py` (DEAD-3).
- [ ] `git rm` tracked build artifacts: `rhealpixdggs/.coverage`, `tests/test_PJ_healpix.c.py.bak`, `docs/.doctrees/*`; investigate and likely delete orphaned `tests/test_PJ_healpix.c.txt`/`c_degrees.txt` (DEAD-4).
- [ ] Fix `pyproject.toml:13`'s `license` field to reflect the actual dual LGPL-3.0-or-later/MIT terms; remove or replace the stale commented-out classifier at line 34 (PKG-1).
- [ ] Update `docs/source/conf.py`'s `copyright`/`release` values (DOC-3).
- [ ] Commit `poetry.lock` (PKG-2).

*(Real bugs — BUG-1 through BUG-5, LAT-1, LAT-2 — are one-line/small fixes too, but are listed in the detailed findings above rather than here, since each needs a moment's understanding of why before applying, not just deletion/text replacement.)*
