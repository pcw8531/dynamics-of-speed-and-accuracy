# Analyses of the open datasets

One file per test, written and fixed before the test is run. Each plan states the dataset (registry ID), the quantity computed, the prediction it tests, the exact statistic, and what result would count against the prediction.

Plans and their status:
- plan_P2.md, Prediction 2 in How We Type (fixed 4 November 2026, run once).
- plan_P4_howwetype.md, Prediction 4 in How We Type (fixed 4 November 2026, run once).
- plan_P5_axon.md, Prediction 5 in the Axon data (fixed 4 November 2026, run once).
- plan_throwing.md, coupling across throwing practice (fixed 4 November 2026, run once).
- plan_karate.md, Prediction 4 at the level of joints in Kyokushin karate (fixed 4 November 2026 at 14:15 UTC, before any data file was obtained, run once on 5 November 2026).

The results of each test are in its results and output files, and Supplemental Section S7 of the paper reports them.

## Running the analyses

Each script runs from the repository root on the raw files of its dataset, obtained from the host that data/REGISTRY.md names and checked against data/checksums. It prints its statistics, which the run kept as its output file, and writes its tables of values into this folder.

- python -I analysis/p2_test.py <typing folder of *_matched.txt> <keyboard_flat_coordinates.csv> <motion capture folder>, output p2_result.txt, table p2_per_typist.csv.
- python -I analysis/p4_howwetype.py <typing folder of *_matched.txt>, after p2_test.py, output p4_howwetype_result.txt, table p4_per_typist.csv.
- python -I analysis/p5_axon.py <data_by_cookie.json>, output p5_axon_result.txt, table p5_mean_curve.csv.
- python -I analysis/throwing_test.py <folder with S1.zip to S7.zip and Trialinfo.xlsx>, output throwing_result.txt, tables throwing_per_throw.csv and throwing_per_session.csv.
- python -I analysis/karate_test.py <folder with the athlete folders of Karate.zip> <clone of py-c3d at commit b3060462bf35>, output karate_result.txt, tables karate_per_recording.csv and karate_per_athlete.csv.

feasibility_counts_P2.py gives the counts that plan_P2.md cites (feasibility_P2_howwetype.md), and throwing_diagnostic_links.py the diagnostic reported in results_throwing.md. The usage of each is at its head. The figures of these runs are drawn from the tables by figures/fig_p2_p4_howwetype.py, fig_throwing.py, fig_karate.py and fig_p5_axon.py.

## Labels in the plans and in the article

The plans, scripts and outputs keep the labels of the draft of their date. The article has since renumbered its equations and no longer numbers its predictions. The labels correspond as follows.

- Prediction 2 in plan_P2.md, feasibility_P2_howwetype.md, feasibility_counts_P2.py, p2_test.py, results_howwetype.md and figures/fig_p2_p4_howwetype.py is the claim that the slope of a law is proportional to the averaged information of the coupling, b = βH. The ratio of slopes that plan_P2.md calls Equation 11 is Equation 10 of the article.
- Prediction 4 in plan_P4_howwetype.md, p4_howwetype.py, results_howwetype.md, figures/fig_p2_p4_howwetype.py, plan_throwing.md and plan_karate.md is the claim that lower averaged information goes with higher skill. karate_test.py, karate_result.txt, results_karate.md and figures/fig_karate.py, from the run of 5 November 2026, call the same claim Prediction 2, its label in the draft of that date.
- Prediction 5 in plan_P5_axon.md, p5_axon.py and results_axon.md is the claim that performance changes quickly and then slowly across practice. The curve that these files, p5_axon_result.txt and figures/fig_p5_axon.py call Equation 18 is Equation 16 of the article, and the column fit_eq18 of p5_mean_curve.csv holds its fit.
- Equations 6 and 7 in the plans, the shares of coupling and their averaged information in bits, are Equations 5 and 6 of the article.
- Equations 13 to 17 in plan_throwing.md, the dynamic claim of the law, are Equations 11 to 17 of the article.
- Section S2.10 in plan_karate.md, the change of H with each unit of coupling, is Supplemental Section S3.2 of the article, and the main-text section When Averaged Information Falls and When It Rises states that H is computed over one fixed set of elements.
