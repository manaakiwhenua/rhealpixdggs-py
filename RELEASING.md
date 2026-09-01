# Releasing

How to publish a new release of rHEALPixDGGS. Requires
[Poetry](https://python-poetry.org/) **version 2.2 or later** — this
project's `pyproject.toml` uses the PEP 621 `[project]` table, which Poetry
before 2.0 cannot read (it fails with a cryptic `'name'` error), and the
built wheel must carry the PEP 639 `License-Expression` metadata, which only
poetry-core 2.2+ writes. Upgrade with `poetry self update`.

`scripts/release.py` drives the whole process. It is stdlib-only, so there is
nothing to install beyond Poetry itself.

## The short version

Write the release notes first (see [Before you start](#before-you-start)),
then, from a clean `master`:

```sh
python scripts/release.py check    0.7.0   # preflight only, changes nothing
python scripts/release.py prepare  0.7.0   # bump versions, test, build, verify
python scripts/release.py tag      0.7.0   # commit, tag, push
python scripts/release.py publish  0.7.0   # upload to PyPI
python scripts/release.py announce 0.7.0   # draft GitHub release: changelog + generated notes
python scripts/release.py finalize 0.7.0   # publish the draft GitHub release
```

Each stage re-runs the checks it depends on, so stopping to fix something
and starting again is safe. Nothing irreversible happens without you asking
for it by name, and every stage visible from outside your machine — `tag`,
which pushes, `publish`, which uploads, `announce`, which drafts the GitHub
release, and `finalize`, which makes it public — also prompts before acting. Add `--dry-run` to any stage to
see what it would do, or `--yes` to skip the prompts; either flag works
before or after the stage name.

## Before you start

The script checks everything it can, but two things are yours to write:

- [ ] **Release notes in `CHANGES.rst`**, under a section header for the new
      version, with explicit **Breaking change:** callouts where applicable.
      `check` fails if the section is missing or empty, but it cannot tell
      you whether the notes are any good.
- [ ] **The copyright year range in `LICENSE` and `LICENSE-MIT`**, if the
      year has rolled over. `check` fails if the range ends before the
      current year. (The docs' copyright line computes its end year at build
      time; the license files are maintained by hand.)

## What each stage does

### `check`

Read-only. Verifies that Poetry is present and new enough; that you are on a
clean `master` in step with `origin`; that the version is well formed and
newer than the current one; that `v<VERSION>` is not already tagged locally
or on `origin`; that `CHANGES.rst` has notes for it; that the copyright
years are current; and, via `gh`, that CI is green on the exact commit being
released. Run it as often as you like.

### `prepare`

Sets `version` in `pyproject.toml` (the single source of truth — the docs
and package metadata both read it from there) and `version` plus
`date-released` in `CITATION.CFF`. Then runs `scripts/run_unittests.sh` and
`scripts/run_doctests.sh`, clears `dist/` (stale artifacts are easy to
upload by accident), builds, and verifies the results:

- the wheel and sdist filenames carry the new version;
- the wheel's `METADATA` declares that version, the
  `LGPL-3.0-or-later OR MIT` license expression, all three license files,
  and a `text/markdown` README (without which PyPI will not render it);
- the sdist ships `tests/` and `docs/source/` and has **not** swept up build
  leftovers. Poetry's `include` patterns override gitignore, so anything
  transient under `docs/` stows away in the sdist — this is the check for
  that.

Stop here and `git diff` before going on.

### `tag`

Commits the version bump as `Release <VERSION>`, tags it `v<VERSION>` — the
`v` prefix matches every existing tag — and pushes the branch and the tag.
Prompts first, and shows you exactly what it will commit.

### `publish`

Uploads the built artifacts, after re-verifying them and confirming the tag
exists. Rehearse against TestPyPI first if you want:

```sh
python scripts/release.py publish 0.7.0 --test-pypi
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple rhealpixdggs
python scripts/release.py publish 0.7.0
```

A released version number cannot be reused or replaced, which is why this
stage asks twice over: once for the tag, once at the prompt.

### `announce`

Creates a **draft** GitHub release for `v<VERSION>` via `gh` — titled
`v<VERSION>`, matching the tag and every earlier release. Nothing is
publicly visible until you run `finalize`.

The body is the version's `CHANGES.rst` section, followed by GitHub's
auto-generated notes (the "What's Changed" PR list, new contributors, and
the full-changelog link). The changelog part is converted from
reStructuredText to Markdown on the way — ``literals`` become backticks and
`--` becomes an em dash; **bold**, bullet lists and `#123` issue references
already mean the same thing in both, and GitHub links the issue references
itself. Only that subset is handled, so if you reach for anything more
elaborate in the changelog, check the rendering.

The notes are printed for you to read before anything is created, and the
stage refuses to clobber a release that already exists. A version like
`0.7.0rc1` is marked as a prerelease automatically.

```sh
python scripts/release.py announce 0.7.0 --attach   # also upload the wheel and sdist
```

### `finalize`

Publishes the draft release that `announce` created, after showing you it
exists and prompting. Between the two stages, review the draft in the
browser — the auto-generated notes in particular — and edit it there if
anything needs touching up; `finalize` publishes whatever the draft says by
the time you run it.

## After publishing

The **conda package** is maintained at
[conda-forge/rhealpixdggs-feedstock](https://github.com/conda-forge/rhealpixdggs-feedstock);
its bot normally opens a version-bump PR automatically after the PyPI
release appears.

## Doing it by hand

If `scripts/release.py` is in the way, the underlying steps are:

```sh
# 1. edit pyproject.toml (version), CITATION.CFF (version, date-released)
./scripts/run_unittests.sh
./scripts/run_doctests.sh

# 2. tag
git commit -am "Release <VERSION>"
git tag v<VERSION>
git push
git push --tags

# 3. build
rm -rf dist/
poetry build

# 4. publish
poetry publish

# 5. draft GitHub release, with the CHANGES.rst entry as the body plus
#    GitHub's generated notes; then publish the draft once reviewed
gh release create v<VERSION> --title v<VERSION> --notes-file <notes.md> --generate-notes --draft
gh release edit v<VERSION> --draft=false
```
