# Release manifest

Every file behind the article falls into one of three groups. Only group 1 is uploaded with this repository. Group 2 is
cited by its DOI or address, with the SHA-256 of the exact files used, so anyone can obtain the same files and check
them. Group 3 is cited by repository and commit.

## 1. Uploaded with this repository

These are files made for this article and small enough to hold here. Values derived from the open datasets keep
the terms of their source, as LICENSE states.

- code/: sim_practice.py (minimal model, Figure 4), sim_full.py (full model, Figures 5 and 11, Section S5.3), sim_lambda.py
  (directed practice, Section S5.3), urn_check.py (Section S3.2), sim_static_cases.py (the classical laws in one
  simulated performer, Section S4.4), force_spread.py (the impulse-variability equation at two steps of practice,
  Section S4.5, Figure 9), sim_composition.py, sim_flows.py and sim_networks.py (where the coupling goes, Figures 4B
  and 5B), network_models.py (the five network models of Table S3, Section S2.1), check_transcriptions.py,
  read_figures.py, README_sim_full.md.
- sim_out/: the outputs of the simulations as run.
- figures/: one script and one image per figure, with the shared palette and type in style.py. Manuscript Figures 1
  to 11 (fig1_classical_laws, fig2_mechanism, fig3_practice_step, fig4_minimal_model, fig5_full_model, fig6_hick,
  fig7_hyman, fig8_fitts, fig9_impulse_variability, fig10_slope_reading, fig11_system_model), Supplemental Figures S1
  to S4 (figS1_network_models, figS2_static_cases, figS3_practice_tables, figS4_open_records), and the figures of the
  analyses of Section S7 as run (fig_p2_p4_howwetype, fig_throwing, fig_karate, fig_p5_axon). figures/panels/ holds
  every lettered panel of Figures 1, 2, 4 to 11 and S1 to S4 as its own image, and fig2_centre the unlettered strip of
  Figure 2. Figure 3 is one chart without panels and is placed whole.
- notebooks/: one notebook for each figure of the manuscript and the Supplemental Material, 15 in all, with the code
  of its script and the values that the code prints, and no stored plot.
- tools/: make_notebooks.py, which writes and runs the notebooks.
- data/transcribed/: values transcribed from published tables, each with its page, read twice. Woodworth (1899)
  Tables I to IV, Fitts (1954) Tables 1 to 3, Hyman (1953) Tables 1 and 2 and the component times of p. 194, Mowbray
  and Rhoades (1959) Tables I to III, Szczęsna et al. (2021) Table 1.
- data/read_from_figures/: values given only in figures, each read twice. Hyman (1953) Figure 1, participant G. C.,
  Experiments I to III (Figures 1B, 6 and 7 and Supplemental Figure S2), and Schmidt et al. (1979) Figures 8 and 9
  (Figures 1D and 9A).
- data/open_extracted/: every value the paper extracts from the open datasets of group 2, one table for each
  analysis, with a note on every column.
- data/checksums/: SHA-256 of the files used from group 2 (howwetype_motion_capture, karate_selected, throwing_files)
  and of the extracted tables (open_extracted).
- data/REGISTRY.md: every source, with its use, licence and checksums.
- data/check_transcriptions_result.txt: the output of code/check_transcriptions.py, which checks every transcribed
  table and every value read from a figure.
- analysis/: the plan of each analysis of an open dataset, written before it was run, its script, its output as run and its
  derived values per participant or per recording. Plans and run logs keep the labels of the draft of their date, and
  analysis/README.md relates them to the claims and equations of the article.
- README.md, MANIFEST.md, LICENSE, CITATION.cff, requirements.txt and .gitignore.

## 2. Cited, not uploaded

Open datasets, used under the terms of their release. Their raw files are not copied here, because of their size,
and because two of them limit redistribution, but every value the paper takes from them is held in
data/open_extracted.

| Source | Where | Licence | Files used |
|---|---|---|---|
| Online game (Stafford & Dewar, 2014b) | github.com/tomstafford/axongame, commit 58c163a895cf | MIT code, no data licence stated | data_by_cookie.json, SHA-256 in data/REGISTRY.md |
| How We Type (Feit et al., 2016b) | Zenodo 10.5281/zenodo.4047775 | CC BY-NC 4.0 | typing data and motion capture, SHA-256 in data/REGISTRY.md and data/checksums |
| Throwing practice (Ozkaya et al., 2018) | figshare 10.6084/m9.figshare.c.4017808 | CC BY 4.0 | per-subject CSV files and Trialinfo.xlsx, data/checksums/throwing_files.sha256 |
| Kyokushin karate (Szczęsna et al., 2020) | figshare 10.6084/m9.figshare.c.4981073 | CC BY 4.0 | 146 C3D files, data/checksums/karate_selected.sha256 |
| Typing performance (Dhakal et al., 2018) | userinterfaces.aalto.fi/136Mkeystrokes | research use with attribution | values reported in the paper only, no reanalysis |

The classical papers (Woodworth, 1899; Hick, 1952; Hyman, 1953; Fitts, 1954; Mowbray & Rhoades, 1959; Schmidt et al.,
1979) are cited, and only the values transcribed from them are held here, in group 1.

## 3. Cited by repository and commit

- Park (2026), the released code and data of the earlier paper: github.com/pcw8531/control-entropy-coordinative-structures,
  commit 0660a649f5e0 (Zenodo 10.5281/zenodo.21789488). Figure 2 and Supplemental Figure S1 import its generators, and
  the musculoskeletal value is taken from it.
- py-c3d, the C3D reader used in the karate test: github.com/EmbodiedCognition/py-c3d, commit b3060462bf35, MIT.
