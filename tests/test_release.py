"""
Tests for the release script's changelog-to-release-notes conversion.
"""

import importlib.util
import pathlib
import unittest

_SPEC = importlib.util.spec_from_file_location(
    "release", pathlib.Path(__file__).resolve().parents[1] / "scripts" / "release.py"
)
assert _SPEC is not None and _SPEC.loader is not None
release = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(release)


class RstToMarkdownTestCase(unittest.TestCase):
    def test_unwraps_paragraphs_and_keeps_paragraph_breaks(self):
        # GitHub renders every newline in release notes as a line break, so
        # hard-wrapped changelog paragraphs must come out on one line each,
        # with the blank line between entries preserved.
        rst = (
            "``Cell.centroid`` on dart cells integrates with a fixed rule\n"
            "instead of adaptive quadrature -- much faster (issue #120).\n"
            "\n"
            "``Cell.nw_vertex`` on dart cells picks the vertex in the\n"
            "plane (issue #122).\n"
        )
        self.assertEqual(
            release.rst_to_markdown(rst),
            "`Cell.centroid` on dart cells integrates with a fixed rule "
            "instead of adaptive quadrature — much faster (issue #120).\n"
            "\n"
            "`Cell.nw_vertex` on dart cells picks the vertex in the "
            "plane (issue #122).",
        )

    def test_lists_and_indented_blocks_keep_their_lines(self):
        rst = (
            "Two fixes:\n"
            "\n"
            "- the first one, which is described over\n"
            "  two lines\n"
            "- the second one\n"
            "\n"
            "Example::\n"
            "\n"
            "    literal line one\n"
            "    literal line two\n"
        )
        self.assertEqual(
            release.rst_to_markdown(rst),
            "Two fixes:\n"
            "\n"
            "- the first one, which is described over\n"
            "  two lines\n"
            "- the second one\n"
            "\n"
            "Example::\n"
            "\n"
            "    literal line one\n"
            "    literal line two",
        )

    def test_literal_spanning_a_wrap_survives(self):
        rst = "Call ``rhp_to_geo(index,\nplane=False)`` for the centroid.\n"
        self.assertEqual(
            release.rst_to_markdown(rst),
            "Call `rhp_to_geo(index, plane=False)` for the centroid.",
        )

    def test_released_changelog_sections_keep_one_paragraph_per_entry(self):
        # CHANGES.rst separates entries with blank lines so that announce
        # renders one paragraph per entry. Released sections are history and
        # should keep their structure: 0.8.0 has five entries, 0.7.1 two.
        for version, entries in (("0.8.0", 5), ("0.7.1", 2)):
            notes = release.rst_to_markdown(release.changelog_section(version))
            paragraphs = [p for p in notes.split("\n\n") if p.strip()]
            self.assertEqual(len(paragraphs), entries, version)
            for paragraph in paragraphs:
                self.assertNotIn("\n", paragraph, f"{version}: {paragraph[:80]}")
