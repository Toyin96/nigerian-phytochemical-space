# Chemical space of Nigerian medicinal plant compounds

A reproducible analysis of the structural diversity of 1,813 compounds catalogued from
Nigerian medicinal plants, compared against two reference libraries from mainstream drug
discovery.

In 2018 I spent five months in four Hausa-Fulani communities in Lagos State documenting
medicinal plant use for my undergraduate thesis in Applied Botany at the University of Lagos.
I recorded 184 species, the parts used and how each remedy was prepared. The study ended with
a question fieldwork could not answer: what is actually inside these plants at the molecular
level? This repository is a first computational pass at that question.

**It is a descriptive structural analysis. It makes no claim that any compound is active,
safe or useful against any disease.**

## The question

> How structurally diverse are the compounds already catalogued from Nigerian medicinal
> plants, and how much of that chemical space is represented in the synthetic libraries that
> modern drug discovery screens against?

## Data

| Set | Compounds | Source |
|---|---|---|
| Nigerian medicinal plants | 1,813 | [NMPP-DB v1.0](https://nmppdb.com.ng) |
| Synthetic drug-like | 5,000 | random sample of ZINC 250k |
| Bioactive medicinal chemistry | 5,531 | ChEMBL EGFR set, via [TeachOpenCADD](https://github.com/volkamerlab/teachopencadd) |

The reference sets are comparators, not targets. The ChEMBL set is used only because it is a
well-curated sample of compounds medicinal chemists have actually worked on. Nothing is
screened against EGFR or any other protein.

## What the analysis found

After standardising structures, removing duplicates and computing descriptors, fingerprints
and scaffolds:

| | Nigerian plants | ZINC synthetic | ChEMBL bioactives |
|---|---|---|---|
| Median Fsp3 (three-dimensionality) | 0.73 | 0.39 | 0.21 |
| Median stereocentres | 1 | 1 | 0 |
| Unique Murcko scaffolds per 100 molecules | 28.7 | 82.3 | 33.7 |
| Scaffolds appearing in only one compound | 72.0% | 92.6% | 66.5% |

Fingerprint clustering (Butina, Tanimoto 0.6) splits the Nigerian set into 1,089 clusters
from 1,813 compounds, 80% of them singletons.

**Chemical space overlap.** For each Nigerian compound, the Tanimoto similarity to its most
similar molecule in each reference library:

| | vs ZINC | vs ChEMBL |
|---|---|---|
| Median nearest-neighbour similarity | 0.25 | 0.22 |
| No neighbour above 0.40 | 95.6% | 90.3% |
| Near-duplicate present (>= 0.85) | 0.0% | 0.5% |

## Sensitivity analysis: restricting by molecular weight

The Nigerian set skews small and simple: median molecular weight 222, and 610 of 1,813
compounds have no ring at all. Fatty acids and simple terpenes are abundant. That is a
confound, because an unbranched fatty acid scores a high Fsp3 for an uninteresting reason and
may lack a close neighbour in a drug-like library simply because nobody stocks fatty acids.

Repeating the analysis on the subset of every library within molecular weight 250 to 600 Da
with at least one ring. This is a sensitivity analysis rather than a matched design: the
subsets share a molecular weight window, but their distributions within it are not
statistically matched, so the test asks whether the findings survive the restriction rather
than removing size as a variable.

| | All compounds | MW-restricted |
|---|---|---|
| Median Fsp3 | 0.73 | 0.50 (ZINC 0.38, ChEMBL 0.21) |
| Median stereocentres | 1 | 3 (ZINC 1, ChEMBL 0) |
| No neighbour above 0.40 vs ZINC | 95.6% | 97.7% |
| No neighbour above 0.40 vs ChEMBL | 90.3% | 80.0% |

The restriction weakens the Fsp3 finding, strengthens the stereochemistry finding and leaves
the novelty finding essentially intact. The drop against ChEMBL is expected: medicinal chemistry
has drawn on natural product scaffolds for decades.

## Limitations

- **Structural novelty is not biological activity.** The overwhelming majority of novel
  structures do nothing of interest. Establishing otherwise needs assay data, which this
  analysis does not have and does not substitute for.
- **NMPP-DB reflects what has been isolated and published**, which is a biased sample of what
  these plants contain. Well-studied species are over-represented, and absence from the
  database is not evidence of absence from the plant.
- **One reference sample.** The ZINC comparator is 5,000 molecules from a single curated
  subset; a different library would shift the overlap numbers.
- **PCA on binary fingerprints** captures a small share of total variance, so the projection
  in the notebook illustrates the measured overlap rather than demonstrating it.
- **Stereocentre counts include unassigned centres**, which inflates the count for structures
  deposited without full stereochemistry.

## Repository layout

```
.
├── notebooks/01_chemical_space.ipynb   executed analysis with outputs
├── src/chemspace.py                    standardisation, descriptors, fingerprints, clustering
├── tests/test_chemspace.py             unit tests for the helpers
├── data/raw/                           input files (see below)
├── results/                            generated CSV tables
└── figures/                            generated plots
```

## Reproducing

```bash
pip install -r requirements.txt
python -m pytest tests -q
jupyter nbconvert --to notebook --execute --inplace notebooks/01_chemical_space.ipynb
```

`data/raw/` needs three files. Two download directly:

```bash
mkdir -p data/raw
curl -L -o data/raw/EGFR_compounds.csv \
  https://raw.githubusercontent.com/volkamerlab/teachopencadd/master/teachopencadd/talktorials/T001_query_chembl/data/EGFR_compounds.csv
curl -L -o data/raw/zinc250k.csv \
  https://raw.githubusercontent.com/aspuru-guzik-group/chemical_vae/master/models/zinc/250k_rndm_zinc_drugs_clean_3.csv
```

The third, `data/raw/nmppdb.smi`, is downloaded manually from
[nmppdb.com.ng/download](https://nmppdb.com.ng/download), choosing the SMILES format.

All randomness is seeded (`SEED = 42`), so the figures and tables reproduce exactly.

## Next steps

The honest next step is not a larger model but better data: linking these structures to
measured bioactivity where it exists, and treating absence of data as absence of data rather
than as a negative result.

## Author

Toyin Michael Onagoruwa. B.Sc. Applied Botany, University of Lagos. Software engineer,
currently working on backend and data systems in Lagos.

## Acknowledgements and licences

Compound data from the Nigerian Medicinal Plants and Phytochemicals Database (NMPP-DB),
ZINC and ChEMBL, each under its own terms; please cite the original sources rather than this
repository. The code here is released under the MIT licence.
