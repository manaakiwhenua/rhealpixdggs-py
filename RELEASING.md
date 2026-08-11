# Releasing

How to publish a new release of rHEALPixDGGS. Requires
[Poetry](https://python-poetry.org/) **version 2.0 or later** — this
project's `pyproject.toml` uses the PEP 621 `[project]` table, which older
Poetry versions cannot read (they fail with a cryptic `'name'` error).

## 1. Prepare the release

- [ ] Update the version number in `pyproject.toml` (the single source of
      truth: the docs and package metadata read it from there).
- [ ] Update `CHANGES.rst` with the release notes, including explicit
      **Breaking change:** callouts where applicable.
- [ ] Update `CITATION.CFF`: `version` and `date-released`.
- [ ] Check the copyright year range in `LICENSE` and `LICENSE-MIT` covers
      the current year (the docs' copyright line computes its end year at
      build time, but the license files are updated by hand).
- [ ] Run the tests and doctests:

  ```sh
  ./run_unittests.sh
  ./run_doctests.sh
  ```

## 2. Tag

```sh
git commit -am "Release <VERSION>"
git tag <VERSION>
git push
git push --tags
```

## 3. Build

Stale artifacts in `dist/` are easy to upload by accident, so clear it
first:

```sh
rm -rf dist/
poetry build
```

This produces the sdist (`.tar.gz`) and wheel (`.whl`) in `dist/`. Sanity
check: the wheel filename carries the new version number.

## 4. Publish

Optionally rehearse against TestPyPI:

```sh
poetry publish --repository testpypi
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple rhealpixdggs
```

Then publish for real:

```sh
poetry publish
```

## 5. After publishing

- The **conda package** is maintained at
  [conda-forge/rhealpixdggs-feedstock](https://github.com/conda-forge/rhealpixdggs-feedstock);
  its bot normally opens a version-bump PR automatically after the PyPI
  release appears.
- Consider creating a GitHub release from the tag, pasting the
  `CHANGES.rst` entry.
