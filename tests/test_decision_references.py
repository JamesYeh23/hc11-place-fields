"""Every decision number cited in the repository has to resolve to an entry.

Step 10b's commit message announced D10.6, D10.7 and D10.8, `docs/worklog.md`
cited D10.8, and none of the three had been written to `docs/decisions.md`. The
dangling cross-reference sat in a pushed public document for a day. Nothing in
the suite could have noticed, because a decision log is prose and prose does not
fail a test.

This is the smallest check that would have caught it: collect every ``Dx.y``
written anywhere in the tracked text of the repository, collect every ``Dx.y``
that `decisions.md` actually defines as a heading, and insist the first set sits
inside the second. It cannot tell whether an entry says what the citation claims
-- only that something is there to read. Commit messages are out of reach.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DECISIONS = ROOT / "docs" / "decisions.md"

#: A decision number as it is written in prose and in comments.
CITATION = re.compile(r"\bD(\d+)\.(\d+)\b")
#: A decision number as `decisions.md` defines one: a level-3 heading.
DEFINITION = re.compile(r"^### D(\d+)\.(\d+)\b", re.MULTILINE)

SEARCHED = ("docs/*.md", "*.md", "src/hc11/*.py", "scripts/*.py", "tests/*.py",
            "config/*.yaml")


def defined() -> set[str]:
    text = DECISIONS.read_text(encoding="utf-8")
    return {f"D{major}.{minor}" for major, minor in DEFINITION.findall(text)}


def cited() -> dict[str, set[str]]:
    """Every decision number in the tracked text, by the file that cites it."""
    out: dict[str, set[str]] = {}
    for pattern in SEARCHED:
        for path in sorted(ROOT.glob(pattern)):
            if path == Path(__file__):
                continue
            text = path.read_text(encoding="utf-8")
            if path == DECISIONS:
                text = DEFINITION.sub("", text)
            for major, minor in CITATION.findall(text):
                out.setdefault(f"D{major}.{minor}", set()).add(
                    str(path.relative_to(ROOT)))
    return out


def test_decisions_md_defines_at_least_one_entry():
    """A guard on the guard: a parser that finds nothing would pass everything."""
    assert len(defined()) > 50


def test_every_cited_decision_number_resolves_to_an_entry():
    known = defined()
    missing = {number: sorted(files) for number, files in cited().items()
               if number not in known}
    assert not missing, (
        "decision numbers cited with no entry in docs/decisions.md: "
        + "; ".join(f"{number} (cited in {', '.join(files)})"
                    for number, files in sorted(missing.items()))
        + ". Write the entry or correct the citation -- a reader who follows the "
          "reference has to find something.")


def test_no_decision_number_is_defined_twice():
    text = DECISIONS.read_text(encoding="utf-8")
    found = [f"D{major}.{minor}" for major, minor in DEFINITION.findall(text)]
    duplicates = sorted({number for number in found if found.count(number) > 1})
    assert not duplicates, f"defined more than once in decisions.md: {duplicates}"


@pytest.mark.parametrize("number", ["D10.6", "D10.7", "D10.8"])
def test_the_three_entries_that_dangled_are_present(number):
    """Named explicitly, because these are the ones that got away."""
    assert number in defined()
