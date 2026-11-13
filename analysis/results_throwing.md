# Results in the throwing practice data (4 November 2026)

Run once under plan_throwing.md with the change recorded there at 13:42 UTC (run 13:43 to 13:47 UTC). Raw output: throwing_result.txt. Per-throw and per-session values: throwing_per_throw.csv, throwing_per_session.csv. Figure: figures/fig_throwing.png.

Data kept: 2,117 clean throws, 57 body markers, 171 links per network, ceiling log2 57 = 5.833 bits. Session mean H ranged from about 5.23 to 5.62 bits, 90 to 96 percent of the ceiling, so this sparse measure varies, unlike the dense measure used in How We Type.

## Primary test, H falls with practice: undecided
- Slopes of session mean H per session, averaged over the two arms, for participants 1 to 7: 0.0041, −0.0019, −0.0009, −0.0016, 0.0113, −0.0004, 0.0002 bits per session.
- Four of seven are negative and all are small. One-sided Wilcoxon p = .53, two-sided p = 1.00.

## Secondary tests, reported as run
- H and ball speed within participant: mean Spearman ρ per participant 0.096, −0.073, 0.194, −0.046, 0.452, 0.036, −0.025. The predicted negative relation is not seen (one-sided p = .85).
- Dominant against non-dominant arm: mean H was higher on the dominant arm in six of seven participants (for example 5.548 against 5.499 for participant 1). The predicted lower H on the dominant arm is not seen (one-sided p = .98).
- Practice check: ball speed rose across sessions on the dominant arm for six of seven participants (slopes −0.51 to 1.69 per session) and on the non-dominant arm for six of seven.

## Notes recorded with the results
- One ball speed in Trialinfo.xlsx reads 580 (participant 6, dominant arm, session 3w3d, throw 12) against 46 to 50 for the other throws of that session, which points to a typing error in the source. It is left as recorded. It barely moves the rank-based correlation and raises the practice-check slope for that arm.
- Diagnostic, not a test (throwing_diagnostic_links.py, 200 throws sampled with seed 20261104): about 21 percent of the 171 strongest links join two markers on the same body segment (thigh, shank, upper arm, forearm or foot clusters), against 2.3 percent of all possible pairs. These links are set by the rigid marker clusters and cannot change with practice, so part of the network describes the marker layout rather than coordination.
