"""Checks on the helper functions. Run with: python -m pytest tests -q"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd
import pytest
from rdkit import Chem

import chemspace as cs


def test_standardise_strips_salt():
    """Sodium benzoate should reduce to benzoic acid, neutralised."""
    mol = cs.standardise("[Na+].[O-]C(=O)c1ccccc1")
    assert Chem.MolToSmiles(mol) == "O=C(O)c1ccccc1"


@pytest.mark.parametrize("bad", ["", None, "not_a_molecule", "C(C"])
def test_standardise_rejects_garbage(bad):
    assert cs.standardise(bad) is None


def test_prepare_counts_losses():
    df = pd.DataFrame({"smiles": ["CCO", "OCC", "c1ccccc1", "xxx"]})
    out, audit = cs.prepare(df, label="t")
    assert audit["input records"] == 4
    assert audit["unparseable"] == 1
    assert audit["duplicate structures"] == 1  # CCO and OCC are the same molecule
    assert audit["retained"] == 2 == len(out)


def test_fsp3_endpoints():
    """Benzene is fully sp2, hexane fully sp3."""
    assert cs.describe(cs.standardise("c1ccccc1"))["Fsp3"] == 0.0
    assert cs.describe(cs.standardise("CCCCCC"))["Fsp3"] == 1.0


def test_lipinski_counts_violations():
    assert cs.lipinski_violations({"MW": 300, "cLogP": 2, "HBD": 1, "HBA": 4}) == 0
    assert cs.lipinski_violations({"MW": 700, "cLogP": 7, "HBD": 1, "HBA": 4}) == 2


def test_identical_molecules_have_similarity_one():
    fps = cs.fingerprints([cs.standardise("CCOc1ccccc1"), cs.standardise("CCOc1ccccc1")])
    assert cs.nearest_neighbour_similarity([fps[0]], [fps[1]])[0] == pytest.approx(1.0)


def test_unrelated_molecules_are_dissimilar():
    query = cs.fingerprints([cs.standardise("CCCCCCCCCC")])
    ref = cs.fingerprints([cs.standardise("c1ccc2[nH]ccc2c1")])
    assert cs.nearest_neighbour_similarity(query, ref)[0] < 0.2


def test_scaffold_of_acyclic_molecule_is_empty():
    assert cs.murcko_scaffold(cs.standardise("CCCCCC")) == ""


def test_scaffold_ignores_side_chains():
    """Toluene and ethylbenzene share the benzene scaffold."""
    a = cs.murcko_scaffold(cs.standardise("Cc1ccccc1"))
    b = cs.murcko_scaffold(cs.standardise("CCc1ccccc1"))
    assert a == b == "c1ccccc1"


def test_clustering_groups_near_neighbours_and_separates_chemotypes():
    """Two homologous alkanes belong together; an indole does not belong with them."""
    smiles = ["CCCCCCCCCC", "CCCCCCCCCCC", "c1ccc2[nH]ccc2c1"]
    fps = cs.fingerprints([cs.standardise(s) for s in smiles])
    assignment, clusters = cs.butina_clusters(fps, cutoff=0.4)
    assert assignment[0] == assignment[1]
    assert assignment[2] != assignment[0]
    assert len(clusters) == 2


def test_every_molecule_lands_in_exactly_one_cluster():
    smiles = ["CCCCCCCCCC", "CCCCCCCCCCC", "c1ccc2[nH]ccc2c1", "c1ccc2[nH]cc(C)c2c1"]
    fps = cs.fingerprints([cs.standardise(s) for s in smiles])
    assignment, clusters = cs.butina_clusters(fps, cutoff=0.4)
    assert len(assignment) == len(smiles)
    assert sum(len(c) for c in clusters) == len(smiles)
    assert set(assignment) == set(range(len(clusters)))
