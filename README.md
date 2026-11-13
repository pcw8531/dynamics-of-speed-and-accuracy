# The Dynamics of Speed and Accuracy: code and data

Code and data for the manuscript "The Dynamics of Speed and Accuracy: Motor Coupling and the Slopes of the Classical
Laws" by Chulwook Park. The repository holds the simulations of practice, the values transcribed from the published
records of the classical laws, the plan, script and output of each analysis of an open dataset, one script and one
notebook for each figure, and every output as run. No data were collected for this study. Every measured value
comes from a published table or figure or from a released dataset, as data/REGISTRY.md records.

## Folders

- code: the simulations. sim_practice.py runs the minimal model (Figure 4), sim_full.py the full model (Figure 5,
  Supplemental Section S5.3) and sim_lambda.py its runs at four strengths of directed practice. urn_check.py checks the
  Pólya-urn limit of Supplemental Section S3.2. sim_static_cases.py runs the classical laws of Hick, Hyman and Fitts
  in one simulated performer at steps 0 and 1,000 of the full model (Supplemental Section S4.4), and force_spread.py
  gives the impulse-variability equation at the same two steps (Section S4.5, Figure 9). sim_composition.py
  repeats the runs of Figures 4 and 5 with the same seeds and records where the coupling goes, sim_flows.py runs
  fitness-weighted formation and the fixed condition with selective elimination from six other starting compositions
  (Figures 4B and 5B), and sim_networks.py keeps every link of one run of each, drawn beside those panels.
  network_models.py builds the ring lattice, random, small-world, scale-free and fitness models on 1,000 elements and
  writes Table S3 of the Supplemental Material (Section S2.1). check_transcriptions.py checks every transcribed table
  and every value read from a figure against its own arithmetic and against published fits. read_figures.py repeats
  the automated reading of the two published figures, given the PDFs.
- data: REGISTRY.md lists every source, with its use, licence and checksums. transcribed holds values from published
  tables with their pages, read_from_figures holds values that the sources gave only in figures, each read twice,
  open_extracted holds every value the paper extracts from the open datasets of Supplemental Section S7, and checksums
  holds the SHA-256 of the files used from open datasets and of the extracted tables. check_transcriptions_result.txt
  is the output of code/check_transcriptions.py.
- sim_out: the outputs of the simulations as run.
- figures: one script and one image per figure, in the palette and type of style.py. fig1_classical_laws.py draws Figure
  1 from the transcribed tables and the values read from figures. fig2_mechanism.py draws Figure 2 from the generators
  of Park (2026), imported from that deposit at commit 0660a649f5e0, and from code/sim_full.py. fig3_practice_step.py
  draws Figure 3, one practice step as a flow chart in black, each box marked with the lines of Table 1 of the
  manuscript that carry it out. fig4_minimal_model.py and fig5_full_model.py draw Figures 4 and 5 from the simulations,
  their middle panels from sim_out/sim_flows.npz and sim_out/sim_networks.npz. fig6_hick.py, fig7_hyman.py,
  fig8_fitts.py and fig9_impulse_variability.py draw Figures 6 to 9, the classical laws and the impulse-variability
  equation as static cases. fig10_slope_reading.py draws Figure 10, Leonard's practice read on the path of averaged
  information, from the transcribed table, code/sim_full.py and sim_out/sim_static_cases.npz. fig11_system_model.py
  draws Figure 11, the law as a systems model of practice from its initial state, and the law of the full model at
  successive steps and over practice, with the even ring of step 0 and the network of the simulated performer at step
  1,000, rebuilt with code/sim_full.py from the run stored in sim_out/sim_static_cases.npz. figS1_network_models.py,
  figS2_static_cases.py, figS3_practice_tables.py and figS4_open_records.py draw Supplemental Figures S1 to S4, S1 the
  network models of code/network_models.py with the gap of averaged information and S4 the outputs of Section S7.
  fig_p2_p4_howwetype.py, fig_throwing.py, fig_karate.py and fig_p5_axon.py draw the analyses of Supplemental Section S7
  as run. Each manuscript and supplemental figure script with lettered panels also writes every panel as its own image
  to figures/panels, named by figure and letter (fig4_C.png), cut from the finished figure with the panels of a row cut
  to one height, and the manuscript places these panels one by one. Figure 2 is placed whole, and so is Figure 3, one
  chart without panels. Every figure is set in the regular serif of the text, with bold panel letters only, and narrow
  trend lines (style.bolden). FIG_WEIGHT=bold redraws them in bold type with thicker lines.
- analysis: the plan of each analysis of an open dataset, written before it was run, its script, its output as run and its
  derived values (Supplemental Section S7). analysis/README.md relates the labels of the plans to the claims and
  equations of the article.
- notebooks: one notebook for each figure of the manuscript and the Supplemental Material (fig1_classical_laws.ipynb
  to figS4_open_records.ipynb), with the code of its script in cells and the values that the code prints. No plot is
  stored in a notebook. Each writes the same image files as its script, and its last cell lists them.
  tools/make_notebooks.py writes and runs them.

## Rules

- No synthetic values in data. Every number in the paper traces to a file here, a transcribed table with its page, or
  a script.
- Each plan of an analysis of an open dataset was fixed before it was run, and every result is reported.

## Run, from the repository root

- python -I code/sim_practice.py > sim_out/sim_practice_result.txt   (under a minute, and it also writes
  sim_results.npz in the folder where it is run)
- python -I code/sim_full.py sim_out   (about two minutes)
- python -I code/sim_lambda.py
- python -I code/urn_check.py > sim_out/urn_check_result.txt   (under a minute)
- python -I code/sim_static_cases.py sim_out   (about ten seconds)
- python -I code/sim_composition.py sim_out   (about one and a half minutes, after sim_full.py)
- python -I code/sim_flows.py sim_out   (about two and a half minutes, after sim_composition.py)
- python -I code/sim_networks.py sim_out   (about half a minute)
- python -I code/network_models.py sim_out   (about three and a half minutes)
- python -I code/force_spread.py sim_out   (after sim_static_cases.py)
- python -I code/check_transcriptions.py > data/check_transcriptions_result.txt
- python figures/fig1_classical_laws.py
- python -I figures/fig2_mechanism.py <clone of github.com/pcw8531/control-entropy-coordinative-structures at 0660a649f5e0>
- python figures/fig3_practice_step.py
- python figures/fig4_minimal_model.py   (about half a minute)
- python figures/fig5_full_model.py sim_out
- python figures/fig6_hick.py, fig7_hyman.py and fig8_fitts.py   (about ten seconds each)
- python figures/fig9_impulse_variability.py sim_out
- python figures/fig10_slope_reading.py sim_out   (about ten seconds)
- python figures/fig11_system_model.py sim_out
- python -I figures/figS1_network_models.py <the same clone>
- python figures/figS2_static_cases.py sim_out
- python figures/figS3_practice_tables.py
- python figures/figS4_open_records.py
- python figures/fig_p2_p4_howwetype.py, fig_throwing.py, fig_karate.py and fig_p5_axon.py   (from the tables in
  analysis/)
- python tools/make_notebooks.py . <the same clone>   (about one and a half minutes, writes notebooks/ and the same
  figures again)

The analyses of the open datasets in analysis/ need the raw files of each dataset, obtained from its host and checked
against data/REGISTRY.md and data/checksums. analysis/README.md says how each is run.

## Requirements

Python 3.13 with the packages in requirements.txt, at the versions with which every figure, panel and notebook here
was produced. openpyxl and py-c3d are needed only for the analyses of the throwing and karate datasets.

## Licence

The code is released under the MIT licence. The values, figures and outputs made for this article are released under
CC BY 4.0, and the values derived from the open datasets keep the terms of their source. LICENSE gives the details.

## Citation

Please cite the paper. The Zenodo record of this repository gives a DOI for each release and one DOI for all
versions, and CITATION.cff holds the citation of the code and data.
