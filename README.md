# The Dynamics of Speed and Accuracy: Motor Coupling and the Slopes of the Classical Laws

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24184851.svg)](https://doi.org/10.5281/zenodo.24184851)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.13](https://img.shields.io/badge/Python-3.13-blue.svg)](requirements.txt)

Code and data for the paper. No data were collected for this study. Every measured value comes from a published table
or figure or from a released open dataset, as [data/REGISTRY.md](data/REGISTRY.md) records.

**Chulwook Park** · Department of Physical Education, Seoul National University · [ORCID 0000-0001-8714-5760](https://orcid.org/0000-0001-8714-5760)

---

## The law in one paragraph

The paper proposes a dynamic law of speed and accuracy, T = a + βH(t) log₂X. H(t) is the averaged information of
the coupling among the elements of the motor system after t steps of practice, in bits. Held at one step, the law
gives Hick's law, the Hick–Hyman law and Fitts' law, each with the slope βH(t). Across practice, coupling
concentrates, H(t) falls and every slope falls with it.

## Data sources

Every measured value in the paper comes from a published record or a released dataset. The table lists each source
and the file that holds what the paper takes from it.
[data/REGISTRY.md](data/REGISTRY.md) gives the page of every table, the second reading of every value, the licence and
the checksums.

| Source | What the paper takes | Held in |
| --- | --- | --- |
| Woodworth (1899), Tables I to IV | Error of aimed movements at ten rates, eyes open and shut | `data/transcribed/woodworth1899_tables_I-IV_right_hand.csv` |
| Hyman (1953), Figure 1, Tables 1 and 2, p. 194 | Reaction times of participant G. C. in 24 conditions | `data/read_from_figures/hyman1953_fig1_GC_*.csv`, `data/transcribed/hyman1953_*.csv` |
| Fitts (1954), Tables 1 to 3 | Movement times of reciprocal tapping in Table 1, with Tables 2 and 3 transcribed beside it | `data/transcribed/fitts1954_table*.csv` |
| Schmidt et al. (1979), Figures 8 and 9 | Nine mean effective widths | `data/read_from_figures/schmidt1979_fig8_fig9_distance_We.csv` |
| Mowbray and Rhoades (1959), Tables I to III | Choice times across practice, with those of Leonard (1954) as they table them | `data/transcribed/mowbray_rhoades1959_*.csv` |
| How We Type (Feit et al., 2016), CC BY-NC 4.0 | One row for each of 30 typists, with the slope of same-finger moves between letter keys | `data/open_extracted/p2_per_typist.csv`, `p4_per_typist.csv` |
| Online game (Stafford & Dewar, 2014), no data licence stated | Mean score over the first 30 plays, 30 values | `data/open_extracted/p5_mean_curve.csv` |
| Throwing practice (Ozkaya et al., 2018) and Kyokushin karate (Szczęsna et al., 2020), CC BY 4.0 | Values of the analyses in Supplemental Section S7 | `data/open_extracted/throwing_*.csv`, `karate_*.csv` |

The raw files of the four open datasets are not copied here. Each is obtained from its host and checked against
`data/checksums`.

## Key results and where to check them

Each row gives a result as the paper reports it, the stored file that holds it and the script that produces it.

| Result in the paper | Stored in | Produced by |
| --- | --- | --- |
| Minimal model, H after 2,000 steps: 4.319 bits under random, 4.08 under preferential and 3.72 under fitness-weighted attachment (Figure 4A) | `sim_out/sim_practice_result.txt` | `code/sim_practice.py` |
| Full model, H from 4.32 bits to 3.01 at step 1,000 and 1.65 at step 3,000 (Figure 5A) | `sim_out/sim_full_result.txt` | `code/sim_full.py` |
| Closed forms against the simulated change of H over 60 windows, r = .994 and slope 0.98 (Figure 5C) | `sim_out/sim_full_result.txt` | `code/sim_full.py` |
| Directed practice, slope from 94.7 ms per bit to 85.9, 80.3, 65.9 and 42.6 at step 1,000 | `sim_out/sim_lambda_result.txt` | `code/sim_lambda.py` |
| Selection as a race, 0.685 and 0.703 time units per bit against the predicted 0.693 (Figure S5) | `sim_out/sim_race_result.txt` | `code/sim_race.py` |
| Hick's law, RT = 227 + 148 log₂n, r = .993 (Figure 6A) | `notebooks/fig6_hick.ipynb` | `figures/fig6_hick.py` |
| Hick–Hyman law, RT = 217 + 148 H_S over 24 conditions, r = .985 (Figure 7A) | `notebooks/fig7_hyman.ipynb` | `figures/fig7_hyman.py` |
| Fitts' law, MT = 12.8 + 94.7 log₂(2D/W), r = .983 (Figure 8A) | `data/check_transcriptions_result.txt` | `code/check_transcriptions.py` |
| Impulse-variability equation, W_e = 2.25 + 0.030 D/MT, r = .976 (Figure 9A) | `sim_out/force_spread_result.txt` | `code/force_spread.py` |
| One simulated performer, slopes of 148, 148 and 94.7 ms per bit at step 0 and 103, 103 and 66.0 at step 1,000 (Table 2) | `sim_out/sim_static_cases_result.txt` | `code/sim_static_cases.py` |
| Practice tables, slope of 62, 51 and 45 ms per bit in Leonard's data and of 19 and 2 in Mowbray and Rhoades (Table 3) | `notebooks/figS3_practice_tables.ipynb` | `figures/figS3_practice_tables.py` |
| How We Type, 62 ms per bit in the slower half of the typists and 35 in the faster half (Table 3) | `analysis/p2_result.txt` | `analysis/p2_test.py` |
| Online game, 62.9 percent of 5,293 players and α = 0.275 (Figure S4B) | `analysis/p5_axon_result.txt` | `analysis/p5_axon.py` |

The plans and outputs of every analysis of an open dataset, as reported in Supplemental Section S7, are in
[analysis/](analysis/README.md).

## Quick start

```bash
pip install -r requirements.txt
python -I code/check_transcriptions.py
python -I code/sim_full.py sim_out
python figures/fig5_full_model.py sim_out
```

The first command checks every transcribed table and every value read from a figure, the second runs the full model
of practice (about two minutes), and the third draws Figure 5. Every command runs from the repository root. The
notebooks need Jupyter as well. The full list of commands is at the end of this page.

## Repository layout

| Path | What it holds |
| --- | --- |
| `code/` | The simulations of practice and of the classical laws as static cases, the selection of a command as a race among the elements, the network models, and the checks of the published values |
| `sim_out/` | The outputs of the simulations as run |
| `data/` | Values transcribed from published tables and figures, values extracted from the open datasets, checksums, and the registry of sources |
| `analysis/` | The plan, script and output of each analysis of an open dataset. Its [README](analysis/README.md) relates the labels of the plans to the paper. |
| `figures/` | One script and one image for each figure, and every panel in `figures/panels` |
| `notebooks/` | One notebook for each figure, with its code, the values it prints and the figure it shows |
| `tools/` | The script that builds the notebooks |
| [`MANIFEST.md`](MANIFEST.md) | Every file, and the open datasets that are cited but not copied here |

The analyses in `analysis/` need the raw files of each open dataset, obtained from its host and checked against
`data/checksums`.

## Licence

The code is released under the MIT licence. The values, figures and outputs made for the paper are released under
CC BY 4.0, and the values derived from the open datasets keep the terms of their source, as
[MANIFEST.md](MANIFEST.md#licence-of-the-data) sets out.

## Citation

Please cite the paper. The DOI [10.5281/zenodo.24184851](https://doi.org/10.5281/zenodo.24184851) cites every version of
this repository, and **Cite this repository**, on the right of this page, gives the citation of the code and data.

<details>
<summary>All commands, in order</summary>

```bash
python -I code/sim_practice.py > sim_out/sim_practice_result.txt   # under a minute, also writes sim_results.npz in the current folder
python -I code/sim_full.py sim_out                                 # about two minutes
python -I code/sim_lambda.py
python -I code/urn_check.py > sim_out/urn_check_result.txt         # under a minute
python -I code/sim_static_cases.py sim_out                         # about ten seconds
python -I code/sim_composition.py sim_out                          # after sim_full.py
python -I code/sim_flows.py sim_out                                # after sim_composition.py
python -I code/sim_networks.py sim_out
python -I code/network_models.py sim_out                           # about three and a half minutes
python -I code/force_spread.py sim_out                             # after sim_static_cases.py
python -I code/sim_race.py sim_out                                 # about three minutes
python -I code/check_transcriptions.py > data/check_transcriptions_result.txt

# CLONE is a clone of github.com/pcw8531/control-entropy-coordinative-structures at commit 0660a649f5e0
python figures/fig1_classical_laws.py
python -I figures/fig2_mechanism.py CLONE
python figures/fig3_practice_step.py
python figures/fig4_minimal_model.py
python figures/fig5_full_model.py sim_out
python figures/fig6_hick.py
python figures/fig7_hyman.py
python figures/fig8_fitts.py
python figures/fig9_impulse_variability.py sim_out
python figures/fig10_slope_reading.py sim_out
python figures/fig11_system_model.py sim_out
python -I figures/figS1_network_models.py CLONE
python figures/figS2_static_cases.py sim_out
python figures/figS3_practice_tables.py
python figures/figS4_open_records.py
python figures/figS5_selection_race.py sim_out                     # after sim_race.py
python figures/fig_p2_p4_howwetype.py
python figures/fig_throwing.py
python figures/fig_karate.py
python figures/fig_p5_axon.py
python tools/make_notebooks.py . CLONE                             # writes and runs the 16 notebooks
```

The commands of the analyses of the open datasets are in [analysis/README.md](analysis/README.md).

</details>
