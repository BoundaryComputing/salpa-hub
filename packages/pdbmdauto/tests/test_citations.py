"""Every citation in the package names a work the README's References list.

The References section of README.md is the package's one list of the works
it cites. A node's [citation] table in meta.toml, and the DOIs in the node's
README, must be works from that list, each with its own title.

The pdb2pqr node broke this in both places: it gave the title of the 2004
paper, "PDB2PQR: an automated pipeline for the setup of Poisson-Boltzmann
electrostatics calculations", with the DOI of the 2007 one,
10.1093/nar/gkm276 ("PDB2PQR: expanding and upgrading automated preparation
of biomolecular structures for molecular simulations"). The References cite
the 2004 paper as 10.1093/nar/gkh381.
"""

import re
from pathlib import Path

import pytest

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    tomllib = pytest.importorskip("tomli")

PACKAGE = Path(__file__).resolve().parent.parent
DOI = re.compile(r"10\.\d{4,9}/[^\s)\]>`'\"]*[^\s)\]>`'\".,;]")
# *italic*, not **bold**: the References set each title in italics.
ITALIC = re.compile(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", re.S)


def words(title):
    """A title as words, so case, dashes and punctuation do not count."""
    title = title.replace("–", "-").replace("—", "-").lower()
    return " ".join(re.sub(r"[^a-z0-9-]+", " ", title).split())


def references():
    """The References section as {DOI: the title set in italics before it}."""
    readme = (PACKAGE / "README.md").read_text(encoding="utf-8")
    section = readme.split("\n## References", 1)[1].split("\n## ", 1)[0]
    found = {}
    for entry in section.split("\n- "):
        for doi in DOI.finditer(entry):
            before = [t for t in ITALIC.finditer(entry) if t.end() <= doi.start()]
            found[doi.group()] = before[-1].group(1) if before else ""
    return found


def cited_nodes():
    for meta in sorted(PACKAGE.glob("*/meta.toml")):
        citation = tomllib.loads(meta.read_text(encoding="utf-8")).get("citation")
        if citation and citation.get("doi"):
            yield pytest.param(meta.parent, citation, id=meta.parent.name)


NODES = list(cited_nodes())


def test_the_check_finds_what_it_checks():
    """Without this, a parse that found nothing would pass everything below."""
    assert len(NODES) >= 3, NODES
    assert len(references()) >= 5, references()


@pytest.mark.parametrize("node_dir,citation", NODES)
def test_a_node_cites_a_listed_work_by_its_own_title(node_dir, citation):
    listed = references()
    doi = citation["doi"]
    assert doi in listed, f"{node_dir.name} cites {doi}, which the References do not"
    assert words(citation["reference"]) == words(listed[doi]), (
        f"{node_dir.name} gives {doi} the title {citation['reference']!r}; "
        f"the References give it {listed[doi]!r}"
    )


@pytest.mark.parametrize("node_dir,citation", NODES)
def test_a_node_readme_cites_what_its_meta_toml_cites(node_dir, citation):
    readme = node_dir / "README.md"
    if not readme.exists():
        pytest.skip(f"{node_dir.name} has no README")
    dois = set(DOI.findall(readme.read_text(encoding="utf-8")))
    unlisted = dois - set(references())
    assert not unlisted, f"{node_dir.name}/README.md cites {unlisted}, unlisted"
    assert citation["doi"] in dois, f"{node_dir.name}/README.md omits {citation['doi']}"
