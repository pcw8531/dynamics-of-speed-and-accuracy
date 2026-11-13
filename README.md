# The Dynamics of Speed and Accuracy: Motor Coupling and the Slopes of the Classical Laws

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.13](https://img.shields.io/badge/Python-3.13-blue.svg)](requirements.txt)

Code and data for the paper. No data were collected for this study. Every measured value comes from a published table
or figure or from a released open dataset, as [data/REGISTRY.md](data/REGISTRY.md) records.

**Chulwook Park** · Department of Physical Education, Seoul National University · [ORCID 0000-0001-8714-5760](https://orcid.org/0000-0001-8714-5760)

---

## Quick start

```bash
pip install -r requirements.txt
python -I code/check_transcriptions.py
python -I code/sim_full.py sim_out
python figures/fig5_full_model.py sim_out
```

The first command checks every transcribed table and every value read from a figure, the second runs the full model
of practice (about two minutes), and the third draws Figure 5. Every command runs from the repository root. The
notebooks need Jupyter as well.

## Repository layout

| Path | What it holds |
| --- | --- |
| `code/` | The simulations of practice and of the classical laws as static cases, the network models, and the checks of the published values |
| `sim_out/` | The outputs of the simulations as run |
| `data/` | Values transcribed from published tables and figures, values extracted from the open datasets, checksums, and the registry of sources |
| `analysis/` | The plan, script and output of each analysis of an open dataset. Its [README](analysis/README.md) relates the labels of the plans to the paper. |
| `figures/` | One script and one image for each figure, and every panel in `figures/panels` |
| `notebooks/` | One notebook for each figure, with its code and the values it prints |
| `tools/` | The script that builds the notebooks |
| [`MANIFEST.md`](MANIFEST.md) | Every file, and the open datasets that are cited but not copied here |

The analyses in `analysis/` need the raw files of each open dataset, obtained from its host and checked against
`data/checksums`.

## Licence

The code is released under the MIT licence. The values, figures and outputs made for the paper are released under
CC BY 4.0, and the values derived from the open datasets keep the terms of their source, as
[MANIFEST.md](MANIFEST.md#licence-of-the-data) sets out.

## Citation

Please cite the paper. **Cite this repository**, on the right of this page, gives the citation of the code and data.

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
python figures/fig_p2_p4_howwetype.py
python figures/fig_throwing.py
python figures/fig_karate.py
python figures/fig_p5_axon.py
python tools/make_notebooks.py . CLONE                             # writes notebooks/ and the same figures again
```

The commands of the analyses of the open datasets are in [analysis/README.md](analysis/README.md).

</details>
