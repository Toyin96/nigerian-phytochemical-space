"""Helper functions for the chemical space analysis of Nigerian medicinal plant compounds.

Every step the notebook performs lives here so that it can be tested and
re-run outside Jupyter. Nothing in this module makes a biological or
therapeutic claim: it computes structural descriptors and compares
distributions between sets of molecules.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina

RDLogger.DisableLog("rdApp.*")

SEED = 42

_fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_chooser = rdMolStandardize.LargestFragmentChooser()
_uncharger = rdMolStandardize.Uncharger()


# ---------------------------------------------------------------- loading

def standardise(smiles):
    """Parse SMILES, keep the largest fragment (removes salts and solvents), neutralise charges.

    Returns None when the string cannot be parsed, so callers can count and
    report how many records were lost rather than failing silently.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = _uncharger.uncharge(_chooser.choose(mol))
    except Exception:
        return None
    return mol if mol.GetNumHeavyAtoms() > 0 else None


def load_smi(path, sep="\t"):
    """Read a two-column .smi file (SMILES, identifier)."""
    rows = []
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split(sep)
            if not parts or not parts[0].strip():
                continue
            rows.append((parts[0].strip(), parts[1].strip() if len(parts) > 1 else ""))
    return pd.DataFrame(rows, columns=["smiles", "identifier"])


def prepare(df, smiles_col="smiles", label=""):
    """Standardise a frame of SMILES, drop failures and exact duplicates.

    Returns the cleaned frame and a dictionary recording what was removed,
    so the notebook can report attrition instead of hiding it.
    """
    n_input = len(df)
    out = df.copy()
    out["mol"] = out[smiles_col].map(standardise)
    n_unparsed = int(out["mol"].isna().sum())
    out = out.dropna(subset=["mol"]).copy()
    out["canonical_smiles"] = out["mol"].map(Chem.MolToSmiles)
    n_before_dedup = len(out)
    out = out.drop_duplicates(subset="canonical_smiles").reset_index(drop=True)
    audit = {
        "set": label,
        "input records": n_input,
        "unparseable": n_unparsed,
        "duplicate structures": n_before_dedup - len(out),
        "retained": len(out),
    }
    if label:
        out["set"] = label
    return out, audit


# ----------------------------------------------------------- descriptors

DESCRIPTORS = {
    "MW": Descriptors.MolWt,
    "cLogP": Descriptors.MolLogP,
    "HBD": rdMolDescriptors.CalcNumHBD,
    "HBA": rdMolDescriptors.CalcNumHBA,
    "TPSA": rdMolDescriptors.CalcTPSA,
    "RotB": rdMolDescriptors.CalcNumRotatableBonds,
    "Rings": rdMolDescriptors.CalcNumRings,
    "AromaticRings": rdMolDescriptors.CalcNumAromaticRings,
    "Fsp3": rdMolDescriptors.CalcFractionCSP3,
    "HeavyAtoms": lambda m: m.GetNumHeavyAtoms(),
}


def _count_atom(mol, symbol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == symbol)


def describe(mol):
    """Standard physicochemical and topological descriptors for one molecule."""
    row = {name: fn(mol) for name, fn in DESCRIPTORS.items()}
    row["Stereocentres"] = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                                         useLegacyImplementation=False))
    row["N count"] = _count_atom(mol, "N")
    row["O count"] = _count_atom(mol, "O")
    row["NO ratio"] = row["N count"] / row["O count"] if row["O count"] else np.nan
    return row


def descriptor_table(df):
    """Append descriptor columns to a prepared frame."""
    return pd.concat([df.reset_index(drop=True),
                      pd.DataFrame([describe(m) for m in df["mol"]])], axis=1)


def lipinski_violations(row):
    """Count Lipinski rule-of-five violations. Descriptive only: natural products
    routinely break these rules and remain orally active, so a high count is a
    note about the chemistry, not a verdict on a compound."""
    return int(row["MW"] > 500) + int(row["cLogP"] > 5) + int(row["HBD"] > 5) + int(row["HBA"] > 10)


# ---------------------------------------------------------- fingerprints

def fingerprints(mols):
    """Morgan fingerprints (radius 2, 2048 bits) as RDKit bit vectors."""
    return [_fpgen.GetFingerprint(m) for m in mols]


def fp_matrix(fps):
    """Bit vectors as a dense numpy array, for PCA and other numeric work."""
    arr = np.zeros((len(fps), fps[0].GetNumBits()), dtype=np.uint8)
    for i, fp in enumerate(fps):
        DataStructs.ConvertToNumpyArray(fp, arr[i])
    return arr


def nearest_neighbour_similarity(query_fps, reference_fps, batch=256):
    """For each query molecule, the Tanimoto similarity to its most similar
    molecule in the reference set. Used to measure how much of one library is
    already structurally represented in another."""
    best = np.zeros(len(query_fps))
    for i, fp in enumerate(query_fps):
        sims = DataStructs.BulkTanimotoSimilarity(fp, reference_fps)
        best[i] = max(sims) if sims else 0.0
    return best


# ------------------------------------------------------------- scaffolds

def murcko_scaffold(mol):
    """Bemis-Murcko scaffold as canonical SMILES. Empty string for acyclic molecules."""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return ""


def scaffold_summary(scaffolds):
    """Counts that describe how structurally varied a library's ring systems are."""
    series = pd.Series([s for s in scaffolds if s])
    counts = series.value_counts()
    n_acyclic = sum(1 for s in scaffolds if not s)
    return {
        "molecules": len(scaffolds),
        "acyclic (no scaffold)": n_acyclic,
        "unique scaffolds": int(counts.size),
        "scaffolds per 100 molecules": round(100 * counts.size / max(len(scaffolds), 1), 1),
        "singleton scaffolds": int((counts == 1).sum()),
        "share in singletons": round(float((counts == 1).sum() / max(counts.size, 1)), 3),
    }, counts


# ------------------------------------------------------------ clustering

def butina_clusters(fps, cutoff=0.4):
    """Butina clustering on Tanimoto distance.

    cutoff is a distance, so 0.4 groups molecules that are at least 0.6 similar.
    Returns the cluster assignment per molecule and the clusters themselves.
    """
    n = len(fps)
    distances = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        distances.extend(1.0 - s for s in sims)
    clusters = Butina.ClusterData(distances, n, cutoff, isDistData=True)
    assignment = np.zeros(n, dtype=int)
    for idx, cluster in enumerate(clusters):
        for member in cluster:
            assignment[member] = idx
    return assignment, clusters
