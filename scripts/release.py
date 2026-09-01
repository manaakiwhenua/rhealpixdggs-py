#!/usr/bin/env python3
"""
Drive the rHEALPixDGGS release process described in RELEASING.md.

The release is split into five stages, each its own invocation. Nothing
irreversible happens without you asking for it by name, and every stage that
can be seen from outside the machine -- ``tag``, which pushes, ``publish``,
which uploads, and ``announce``, which posts the GitHub release -- also asks
for confirmation::

    python scripts/release.py check    0.7.0   # preflight only, changes nothing
    python scripts/release.py prepare  0.7.0   # bump versions, test, build, verify
    python scripts/release.py tag      0.7.0   # commit, tag, push
    python scripts/release.py publish  0.7.0   # upload to PyPI (--test-pypi first)
    python scripts/release.py announce 0.7.0   # GitHub release from the changelog

Run them in that order. Each stage re-runs the checks it depends on, so
stopping to fix something and starting again is safe.

Stdlib only, so it runs without installing anything. Requires Python 3.11+
(tomllib) and, for the later stages, Poetry 2.0+ -- this project's
pyproject.toml uses the PEP 621 [project] table, which older Poetry cannot
read.

Add --yes to skip the confirmation prompts, --dry-run to print the commands
a stage would run without running them.
"""

from __future__ import annotations

import argparse
import datetime
import json
import pathlib
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import tomllib
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
PYPROJECT = ROOT / "pyproject.toml"
CITATION = ROOT / "CITATION.CFF"
CHANGELOG = ROOT / "CHANGES.rst"
DIST = ROOT / "dist"

RELEASE_BRANCH = "master"
LICENSE_FILES = ("LICENSE", "LICENSE-LGPL", "LICENSE-MIT")
LICENSE_EXPRESSION = "LGPL-3.0-or-later OR MIT"
# Files whose copyright year range is maintained by hand (the docs compute
# theirs at build time).
COPYRIGHT_FILES = ("LICENSE", "LICENSE-MIT")
# The sdist ships tests/ and docs/ via an include pattern that overrides
# gitignore, so build leftovers under docs/ can stow away. See the note in
# pyproject.toml's [tool.poetry] section.
SDIST_MUST_INCLUDE = ("tests/", "docs/source/")
SDIST_MUST_NOT_INCLUDE = (".doctrees", "_build/", ".cache/", ".aux", ".idx", ".toc")

VERSION_RE = re.compile(r"^\d+\.\d+(\.\d+)?((a|b|rc)\d+)?$")


class Failure(Exception):
    """A check failed; the message is shown to the user and the stage stops."""


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

_dry_run = False
_assume_yes = False


def say(message: str) -> None:
    print(message)


def ok(message: str) -> None:
    print(f"  ok    {message}")


def note(message: str) -> None:
    print(f"  note  {message}")


def run(
    command: list[str], *, capture: bool = True, check: bool = True, cwd: pathlib.Path = ROOT
) -> str:
    """Run a command, echoing it when it changes something."""
    if not capture:
        say(f"  $ {' '.join(command)}")
        if _dry_run:
            return ""
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        capture_output=capture,
        check=False,
    )
    if check and result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()
        raise Failure(
            f"command failed ({result.returncode}): {' '.join(command)}"
            + (f"\n{detail}" if detail else "")
        )
    return (result.stdout or "").strip()


def git(*args: str, check: bool = True) -> str:
    return run(["git", *args], check=check)


def confirm(question: str) -> None:
    if _assume_yes or _dry_run:
        return
    reply = input(f"{question} [y/N] ").strip().lower()
    if reply not in ("y", "yes"):
        raise Failure("cancelled")


def version_key(version: str) -> tuple:
    """Sort key that orders prereleases before their final release."""
    match = re.match(r"^(\d+)\.(\d+)(?:\.(\d+))?(?:(a|b|rc)(\d+))?$", version)
    if not match:
        raise Failure(f"cannot parse version {version!r}")
    major, minor, patch, stage, serial = match.groups()
    stage_order = {"a": 0, "b": 1, "rc": 2, None: 3}[stage]
    return (
        int(major),
        int(minor),
        int(patch or 0),
        stage_order,
        int(serial or 0),
    )


# --------------------------------------------------------------------------
# reading and writing the version in the files that carry it
# --------------------------------------------------------------------------


def pyproject_version() -> str:
    with PYPROJECT.open("rb") as handle:
        return tomllib.load(handle)["project"]["version"]


def citation_version() -> str:
    match = re.search(r"^version:\s*(\S+)\s*$", CITATION.read_text(), re.M)
    if not match:
        raise Failure(f"no 'version:' line in {CITATION.name}")
    return match.group(1).strip("\"'")


def set_pyproject_version(version: str) -> bool:
    text = PYPROJECT.read_text()
    new_text, count = re.subn(
        r'(?m)^(version\s*=\s*")[^"]*(")', rf"\g<1>{version}\g<2>", text, count=1
    )
    if count != 1:
        raise Failure("could not find the version line in pyproject.toml")
    return write_if_changed(PYPROJECT, text, new_text)


def set_citation(version: str, released: datetime.date) -> bool:
    text = CITATION.read_text()
    new_text, count = re.subn(r"(?m)^(version:\s*).*$", rf"\g<1>{version}", text, count=1)
    if count != 1:
        raise Failure("could not find the version line in CITATION.CFF")
    new_text, count = re.subn(
        r"(?m)^(date-released:\s*).*$", rf"\g<1>{released.isoformat()}", new_text, count=1
    )
    if count != 1:
        raise Failure("could not find the date-released line in CITATION.CFF")
    return write_if_changed(CITATION, text, new_text)


def write_if_changed(path: pathlib.Path, old: str, new: str) -> bool:
    if old == new:
        return False
    if _dry_run:
        say(f"  would rewrite {path.name}")
        return True
    path.write_text(new)
    return True


# --------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------


def check_tools(*, need_poetry: bool) -> None:
    if sys.version_info < (3, 11):
        raise Failure("this script needs Python 3.11 or newer (tomllib)")
    if not need_poetry:
        return
    if shutil.which("poetry") is None:
        raise Failure("poetry is not on PATH; this project needs Poetry 2.0 or later")
    reported = run(["poetry", "--version"], check=False)
    match = re.search(r"(\d+)\.(\d+)", reported)
    if not match:
        raise Failure(
            "could not determine the Poetry version -- `poetry --version` said:\n"
            f"{reported or '(no output)'}"
        )
    if int(match.group(1)) < 2:
        raise Failure(
            f"Poetry 2.0 or later is required, found {match.group(0)}. Older versions "
            "cannot read this project's PEP 621 [project] table and fail with a "
            "cryptic \"'name'\" error."
        )
    ok(f"poetry {match.group(0)}")


def check_version_argument(version: str) -> None:
    if not VERSION_RE.match(version):
        raise Failure(
            f"{version!r} does not look like a release version "
            "(expected e.g. 0.7.0, or 0.7.0rc1)"
        )


def check_git_state(version: str, *, expect_bumped: bool) -> None:
    if git("rev-parse", "--is-inside-work-tree", check=False) != "true":
        raise Failure("not inside a git working tree")

    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    if branch != RELEASE_BRANCH:
        raise Failure(
            f"on branch {branch!r}, but releases are cut from {RELEASE_BRANCH!r}"
        )
    ok(f"on {RELEASE_BRANCH}")

    dirty_names = sorted(
        line[3:] for line in git("status", "--porcelain", "--untracked-files=no").splitlines()
    )
    if expect_bumped:
        allowed = {"pyproject.toml", "CITATION.CFF", CHANGELOG.name}
        unexpected = [name for name in dirty_names if name not in allowed]
        if unexpected:
            raise Failure(
                "unexpected uncommitted changes:\n  "
                + "\n  ".join(unexpected)
                + "\ncommit or stash them, then run this stage again"
            )
        ok("working tree carries only the release edits")
    elif dirty_names:
        raise Failure(
            "working tree is not clean:\n  "
            + "\n  ".join(dirty_names)
            + "\ncommit or stash before releasing"
        )
    else:
        ok("working tree clean")

    git("fetch", "--quiet", "--tags", "origin", check=False)
    upstream = git("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}", check=False)
    if upstream:
        behind, ahead = (
            git("rev-list", "--left-right", "--count", f"{upstream}...HEAD").split()
        )
        if int(behind):
            raise Failure(f"{behind} commit(s) behind {upstream}; pull first")
        if int(ahead) and not expect_bumped:
            raise Failure(f"{ahead} unpushed commit(s); push before releasing")
        ok(f"in step with {upstream}")

    tag = f"v{version}"
    if git("tag", "--list", tag):
        raise Failure(f"tag {tag} already exists locally")
    if git("ls-remote", "--tags", "origin", tag, check=False):
        raise Failure(f"tag {tag} already exists on origin")
    ok(f"tag {tag} is free")


def check_version_bump(version: str, *, expect_bumped: bool) -> None:
    current_pyproject = pyproject_version()
    current_citation = citation_version()

    if expect_bumped:
        if current_pyproject != version:
            raise Failure(
                f"pyproject.toml says {current_pyproject}, expected {version}; "
                "run the prepare stage first"
            )
        if current_citation != version:
            raise Failure(
                f"CITATION.CFF says {current_citation}, expected {version}; "
                "run the prepare stage first"
            )
        ok(f"pyproject.toml and CITATION.CFF both say {version}")
        return

    if current_pyproject != current_citation:
        note(
            f"pyproject.toml ({current_pyproject}) and CITATION.CFF "
            f"({current_citation}) disagree; both will be set to {version}"
        )
    if version_key(version) <= version_key(current_pyproject):
        raise Failure(
            f"{version} is not newer than the current version {current_pyproject}"
        )
    ok(f"{current_pyproject} -> {version}")


def changelog_section(version: str) -> str:
    """Return the body of ``version``'s section of CHANGES.rst."""
    lines = CHANGELOG.read_text().splitlines()

    def is_header(index: int) -> bool:
        return (
            index + 1 < len(lines)
            and lines[index].strip() != ""
            and set(lines[index + 1].strip()) == {"^"}
        )

    for index in range(len(lines) - 1):
        if lines[index].strip() == version and is_header(index):
            end = next(
                (j for j in range(index + 2, len(lines) - 1) if is_header(j)),
                len(lines),
            )
            return "\n".join(lines[index + 2 : end]).strip("\n")
    raise Failure(
        f"{CHANGELOG.name} has no section for {version}. Add one before releasing:\n"
        f"    {version}\n    {'^' * len(version)}\n"
        "with the release notes, including explicit **Breaking change:** callouts."
    )


def check_changelog(version: str) -> None:
    if not changelog_section(version).strip():
        raise Failure(f"the {version} section of {CHANGELOG.name} is empty")
    ok(f"{CHANGELOG.name} has release notes for {version}")


def rst_to_markdown(text: str) -> str:
    """
    Convert the small slice of reStructuredText this changelog actually uses.

    Handles ``literals`` and the -- em dashes; **bold**, *emphasis*, bullet
    lists and #123 issue references already mean the same thing in both, and
    GitHub turns the issue references into links by itself. Anything more
    elaborate would need a real converter, so keep the changelog to this
    subset or check the rendering before publishing.
    """
    spans: list[str] = []

    def stash(match: re.Match) -> str:
        spans.append(match.group(1))
        return f"\x00{len(spans) - 1}\x00"

    text = re.sub(r"``(.+?)``", stash, text, flags=re.S)
    text = re.sub(r"(?<=\s)--(?=\s)", "—", text)
    return re.sub(r"\x00(\d+)\x00", lambda m: f"`{spans[int(m.group(1))]}`", text)


def check_copyright_year() -> None:
    this_year = datetime.date.today().year
    for name in COPYRIGHT_FILES:
        path = ROOT / name
        match = re.search(r"Copyright \(c\) (\d{4})-(\d{4})", path.read_text())
        if not match:
            note(f"{name}: no 'Copyright (c) YYYY-YYYY' line found; check it by hand")
            continue
        if int(match.group(2)) < this_year:
            raise Failure(
                f"{name} copyright ends at {match.group(2)}, but it is {this_year}. "
                "Update the year range by hand (the docs compute theirs at build time)."
            )
    ok(f"copyright year ranges cover {this_year}")


def check_ci() -> None:
    """Advisory: ask GitHub whether the commit being released is green."""
    if shutil.which("gh") is None:
        note("gh not installed; skipping the CI check")
        return
    head = git("rev-parse", "HEAD")
    output = run(
        [
            "gh",
            "run",
            "list",
            "--branch",
            RELEASE_BRANCH,
            "--limit",
            "20",
            "--json",
            "headSha,conclusion,status,workflowName",
        ],
        check=False,
    )
    if not output.startswith("["):
        note("could not reach GitHub for the CI check; verify it by hand")
        return
    runs = [run_ for run_ in json.loads(output) if run_.get("headSha") == head]
    if not runs:
        note(f"no CI runs found for {head[:8]}; verify it by hand")
        return
    failed = [
        r for r in runs if r.get("status") == "completed" and r.get("conclusion") != "success"
    ]
    pending = [r for r in runs if r.get("status") != "completed"]
    if failed:
        raise Failure(
            "CI is not green on the commit being released:\n  "
            + "\n  ".join(f"{r['workflowName']}: {r['conclusion']}" for r in failed)
        )
    if pending:
        raise Failure(
            "CI is still running on the commit being released:\n  "
            + "\n  ".join(f"{r['workflowName']}: {r['status']}" for r in pending)
        )
    ok(f"CI green on {head[:8]} ({len(runs)} run(s))")


# --------------------------------------------------------------------------
# build and artifact verification
# --------------------------------------------------------------------------


def run_tests() -> None:
    say("\nRunning the tests and doctests")
    run(["bash", str(SCRIPTS / "run_unittests.sh")], capture=False)
    run(["bash", str(SCRIPTS / "run_doctests.sh")], capture=False)
    ok("tests and doctests pass")


def build(version: str) -> None:
    say("\nBuilding")
    if DIST.exists():
        say(f"  $ rm -rf {DIST.relative_to(ROOT)}")
        if not _dry_run:
            shutil.rmtree(DIST)
    run(["poetry", "build"], capture=False)
    if _dry_run:
        return
    ok(f"built into {DIST.relative_to(ROOT)}/")


def artifact_paths(version: str) -> tuple[pathlib.Path, pathlib.Path]:
    wheels = sorted(DIST.glob("*.whl"))
    sdists = sorted(DIST.glob("*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        raise Failure(
            f"expected exactly one wheel and one sdist in {DIST.name}/, found "
            f"{len(wheels)} wheel(s) and {len(sdists)} sdist(s). Clear the "
            "directory and build again."
        )
    wheel, sdist = wheels[0], sdists[0]
    normalised = version.replace("-", "_")
    for path in (wheel, sdist):
        if normalised not in path.name:
            raise Failure(f"{path.name} does not carry version {version}")
    return wheel, sdist


def verify_artifacts(version: str) -> None:
    say("\nVerifying the built artifacts")
    wheel, sdist = artifact_paths(version)
    ok(f"{wheel.name}")
    ok(f"{sdist.name}")

    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        metadata_name = next(
            (n for n in names if n.endswith(".dist-info/METADATA")), None
        )
        if metadata_name is None:
            raise Failure("the wheel has no METADATA file")
        metadata = archive.read(metadata_name).decode("utf-8")

    def metadata_field(field: str) -> list[str]:
        return re.findall(rf"(?mi)^{re.escape(field)}:\s*(.+?)\s*$", metadata)

    declared = metadata_field("Version")
    if declared != [version]:
        raise Failure(f"wheel METADATA Version is {declared}, expected ['{version}']")
    ok(f"METADATA Version: {version}")

    expression = metadata_field("License-Expression")
    if expression != [LICENSE_EXPRESSION]:
        raise Failure(
            f"wheel METADATA License-Expression is {expression}, "
            f"expected ['{LICENSE_EXPRESSION}']"
        )
    ok(f"METADATA License-Expression: {LICENSE_EXPRESSION}")

    shipped_licences = {pathlib.PurePosixPath(p).name for p in metadata_field("License-File")}
    missing = [name for name in LICENSE_FILES if name not in shipped_licences]
    if missing:
        raise Failure(
            "the wheel is missing license file(s): " + ", ".join(missing)
        )
    ok("all three license files are in the wheel")

    content_type = metadata_field("Description-Content-Type")
    if content_type != ["text/markdown"]:
        raise Failure(
            f"Description-Content-Type is {content_type}, expected ['text/markdown'] "
            "-- PyPI will not render the README otherwise"
        )
    ok("README ships as text/markdown")

    with tarfile.open(sdist) as archive:
        members = archive.getnames()
    stripped = ["/".join(name.split("/")[1:]) for name in members]
    for required in SDIST_MUST_INCLUDE:
        if not any(name.startswith(required) for name in stripped):
            raise Failure(f"the sdist is missing {required}")
    ok("sdist ships tests/ and docs/source/")

    strays = sorted(
        {name for name in stripped for bad in SDIST_MUST_NOT_INCLUDE if bad in name}
    )
    if strays:
        raise Failure(
            "the sdist has swept up build leftovers (Poetry's include patterns "
            "override gitignore):\n  " + "\n  ".join(strays[:20])
        )
    ok("sdist carries no build leftovers")


# --------------------------------------------------------------------------
# stages
# --------------------------------------------------------------------------


def stage_check(version: str) -> None:
    say(f"Preflight checks for {version}\n")
    check_version_argument(version)
    check_tools(need_poetry=True)
    check_git_state(version, expect_bumped=False)
    check_version_bump(version, expect_bumped=False)
    check_changelog(version)
    check_copyright_year()
    check_ci()
    say(f"\nReady to prepare {version}:  python scripts/release.py prepare {version}")


def stage_prepare(version: str) -> None:
    say(f"Preparing {version}\n")
    check_version_argument(version)
    check_tools(need_poetry=True)

    already_bumped = pyproject_version() == version
    check_git_state(version, expect_bumped=already_bumped)
    check_changelog(version)
    check_copyright_year()

    if already_bumped:
        note(f"already at {version}; leaving the version files alone")
    else:
        check_version_bump(version, expect_bumped=False)
        today = datetime.date.today()
        say("\nSetting the version")
        if set_pyproject_version(version):
            ok(f"pyproject.toml -> {version}")
        if set_citation(version, today):
            ok(f"CITATION.CFF -> {version}, released {today.isoformat()}")

    run_tests()
    build(version)
    if not _dry_run:
        verify_artifacts(version)

    say(
        f"\nPrepared {version}. Review the changes:\n"
        "    git diff\n"
        f"then:  python scripts/release.py tag {version}"
    )


def stage_tag(version: str) -> None:
    say(f"Tagging {version}\n")
    check_version_argument(version)
    check_git_state(version, expect_bumped=True)
    check_version_bump(version, expect_bumped=True)
    check_changelog(version)

    tag = f"v{version}"
    pending = git("status", "--porcelain", "--untracked-files=no")
    say("")
    if pending:
        say("Will commit:")
        for line in pending.splitlines():
            say(f"  {line}")
    else:
        note("nothing to commit; the version bump is already committed")
    say(f"Will tag {tag} and push it to origin/{RELEASE_BRANCH}.")
    confirm("\nThis pushes to GitHub. Continue?")

    if pending:
        run(["git", "commit", "-am", f"Release {version}"], capture=False)
    run(["git", "tag", tag], capture=False)
    run(["git", "push", "origin", RELEASE_BRANCH], capture=False)
    run(["git", "push", "origin", tag], capture=False)

    say(
        f"\nTagged and pushed {tag}. Once CI is green:\n"
        f"    python scripts/release.py publish {version} --test-pypi   # rehearsal\n"
        f"    python scripts/release.py publish {version}"
    )


def stage_publish(version: str, *, test_pypi: bool) -> None:
    target = "TestPyPI" if test_pypi else "PyPI"
    say(f"Publishing {version} to {target}\n")
    check_version_argument(version)
    check_tools(need_poetry=True)
    if pyproject_version() != version:
        raise Failure(
            f"pyproject.toml says {pyproject_version()}, expected {version}"
        )
    if not DIST.exists():
        raise Failure(f"no {DIST.name}/ directory; run the prepare stage first")
    verify_artifacts(version)

    if not test_pypi:
        tag = f"v{version}"
        if not git("tag", "--list", tag):
            raise Failure(f"tag {tag} does not exist; run the tag stage first")
        ok(f"tag {tag} exists")

    say("")
    confirm(
        f"Upload {version} to {target}? A released version number cannot be "
        "reused or replaced."
    )
    command = ["poetry", "publish"]
    if test_pypi:
        command += ["--repository", "testpypi"]
    run(command, capture=False)

    if test_pypi:
        say(
            "\nRehearsal uploaded. Check it installs:\n"
            "    pip install --index-url https://test.pypi.org/simple/ "
            "--extra-index-url https://pypi.org/simple rhealpixdggs\n"
            f"then:  python scripts/release.py publish {version}"
        )
    else:
        say(
            f"\nPublished {version}. Now:\n"
            f"    python scripts/release.py announce {version}\n"
            "conda-forge/rhealpixdggs-feedstock's bot normally opens a "
            "version-bump PR by itself once PyPI updates."
        )


def stage_announce(version: str, *, draft: bool, attach: bool) -> None:
    say(f"Creating the GitHub release for {version}\n")
    check_version_argument(version)
    if shutil.which("gh") is None:
        raise Failure(
            "gh is not on PATH; create the release by hand from the "
            f"{CHANGELOG.name} entry, or install the GitHub CLI"
        )

    tag = f"v{version}"
    if not git("ls-remote", "--tags", "origin", tag):
        raise Failure(
            f"tag {tag} is not on origin; run the tag stage first"
        )
    ok(f"tag {tag} is on origin")

    existing = run(["gh", "release", "view", tag, "--json", "url"], check=False)
    if existing.startswith("{"):
        raise Failure(
            f"a GitHub release for {tag} already exists: {json.loads(existing)['url']}\n"
            "Edit it in the browser, or delete it with "
            f"`gh release delete {tag}` and run this stage again."
        )

    notes = rst_to_markdown(changelog_section(version))
    if not notes.strip():
        raise Failure(f"the {version} section of {CHANGELOG.name} is empty")
    ok(f"release notes: {len(notes.splitlines())} lines from {CHANGELOG.name}")

    prerelease = bool(re.search(r"(a|b|rc)\d+$", version))
    if prerelease:
        note(f"{version} is a prerelease; marking it as one")

    assets: list[str] = []
    if attach:
        wheel, sdist = artifact_paths(version)
        assets = [str(wheel), str(sdist)]
        ok(f"attaching {wheel.name} and {sdist.name}")

    say("\n--- release notes ---")
    say(notes)
    say("--- end ---\n")

    if draft:
        note("creating this as a draft; publish it from the browser when ready")
    confirm(
        f"Create the {'draft ' if draft else ''}GitHub release {tag}"
        + (" and upload the artifacts?" if assets else "?")
    )

    # gh reads the body from a file, so the notes survive quoting intact.
    with tempfile.TemporaryDirectory() as workspace:
        notes_file = pathlib.Path(workspace) / f"release-notes-{version}.md"
        notes_file.write_text(notes)
        command = [
            "gh",
            "release",
            "create",
            tag,
            "--title",
            version,
            "--notes-file",
            str(notes_file),
        ]
        if draft:
            command.append("--draft")
        if prerelease:
            command.append("--prerelease")
        command += assets
        run(command, capture=False)

    say(
        f"\nGitHub release {tag} created"
        + (" as a draft; publish it when you are happy with it." if draft else ".")
    )


# --------------------------------------------------------------------------


def main(argv: list[str]) -> int:
    global _dry_run, _assume_yes

    # Accepted either side of the stage name, so both of these work:
    #   scripts/release.py --dry-run prepare 0.7.0
    #   scripts/release.py prepare 0.7.0 --dry-run
    # SUPPRESS keeps the subparser's default from clobbering a flag that was
    # given before the stage name; the values are read back with getattr.
    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument(
        "--yes",
        action="store_true",
        default=argparse.SUPPRESS,
        help="skip confirmation prompts",
    )
    shared.add_argument(
        "--dry-run",
        action="store_true",
        default=argparse.SUPPRESS,
        help="print the commands a stage would run, without running them",
    )

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[shared],
    )
    subparsers = parser.add_subparsers(dest="stage", required=True)

    for name, help_text in (
        ("check", "preflight checks only; changes nothing"),
        ("prepare", "bump the version files, run the tests, build, verify"),
        ("tag", "commit the bump, tag, and push"),
        ("publish", "upload the built artifacts"),
        ("announce", "create the GitHub release from the tag and the changelog"),
    ):
        sub = subparsers.add_parser(name, help=help_text, parents=[shared])
        sub.add_argument("version", help="the version being released, e.g. 0.7.0")
        if name == "publish":
            sub.add_argument(
                "--test-pypi",
                action="store_true",
                help="rehearse against TestPyPI instead of PyPI",
            )
        if name == "announce":
            sub.add_argument(
                "--draft",
                action="store_true",
                help="create it as a draft, to review before it goes live",
            )
            sub.add_argument(
                "--attach",
                action="store_true",
                help="upload the built wheel and sdist as release assets",
            )

    args = parser.parse_args(argv)
    # Keep stdout in step with stderr, so a failure message lands after the
    # lines it follows even when the output is piped to a file.
    sys.stdout.reconfigure(line_buffering=True)
    _dry_run = getattr(args, "dry_run", False)
    _assume_yes = getattr(args, "yes", False)
    if _dry_run:
        say("(dry run: nothing will be changed, pushed, or uploaded)\n")

    try:
        if args.stage == "check":
            stage_check(args.version)
        elif args.stage == "prepare":
            stage_prepare(args.version)
        elif args.stage == "tag":
            stage_tag(args.version)
        elif args.stage == "publish":
            stage_publish(args.version, test_pypi=args.test_pypi)
        elif args.stage == "announce":
            stage_announce(args.version, draft=args.draft, attach=args.attach)
    except Failure as failure:
        print(f"\nstopped: {failure}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\ncancelled", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
