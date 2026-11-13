# Analysis plan for the throwing practice data, fixed before any data were obtained (4 November 2026)

Source: Ozkaya et al. (2018), Three-dimensional motion capture data during repetitive overarm throwing practice, Scientific Data 5, 180272 (registry C1). Seven women with no throwing experience, 15 sessions over 5 weeks, 15 dominant-hand and 15 non-dominant-hand throws per session, marker trajectories at 200 Hz, ball speed of every throw from a speed gun, no feedback on performance. Written from the paper alone. No data file has been opened.

Claim tested: practice concentrates coupling among the moving parts of the body, so the averaged information H of the coupling network falls across sessions (the dynamic claim of the law, Equations 13 to 17), lower H goes with faster throws (Prediction 4), and the dominant side, practised over a lifetime, has lower H than the non-dominant side.

## Trials
- Included: trials marked as having no problem (letter O) in the Trial Information file, excluding the first three throws of each side in each session, which the authors classify as warm-up.
- Window: from the start event to the end event of each throw.

## The averaged information H of one throw
- Markers: the tracking markers present, with no missing value inside the window, in every included trial, so that all throws share one set of N markers.
- Signal: the speed of each marker, the magnitude of its frame-to-frame velocity after a fourth-order Butterworth low-pass filter at 15 Hz applied forward and backward to the positions.
- Coupling: the absolute correlation |r| between the speed series of each pair of markers over the window.
- Network: the 3N pairs with the largest |r| become links, all others are dropped. This keeps a mean degree of 6, the density of the networks in the earlier deposit (300 links on 100 elements, Park, 2026), and lets degrees differ, so that H can fall below its ceiling.
- Entropy: degree shares p_i = k_i / Σ k_j and H = −Σ p_i log2 p_i (Equations 6 and 7).

## Tests
- Primary, H falls with practice: for each participant and side, the least-squares slope of the session mean of H against session number (1 to 15). The two sides are averaged within each participant, giving seven slopes. One-sided Wilcoxon signed-rank test that the slope is below zero.
- Secondary, H and performance: for each participant and side, the Spearman correlation between H and ball speed over the included throws, averaged within participant. One-sided Wilcoxon test that it is below zero.
- Secondary, dominant against non-dominant: each participant's mean H over all sessions on each side. One-sided paired Wilcoxon test that the dominant side is lower.
- Practice check, reported as description: the slope of session mean ball speed against session number, for each side.

## Decision rule
- Support: the primary test gives p < .05 in the predicted direction.
- Against: the primary test gives p < .05 in the opposite direction (two-sided).
- Undecided: any other case. The secondary tests are reported with their own p values and do not change the verdict.

## Reported whatever the outcome
N, the number of throws kept for each participant and side, the seven slopes, the secondary statistics, and a figure of session mean H and ball speed across sessions for each participant and side. Any change to this plan after the data are read is reported as a change, with its reason.

## Change made at the structure check, before any coupling was computed (4 November 2026, 13:42 UTC)
Reason: the data files show that the marker rule as written drops a marker for every throw if it has a gap in any single throw. Applied to the 2,294 readable clean throws, it would keep only 29 of the 58 body markers and remove most of the throwing hand, which defeats the purpose of measuring coupling across the body.

New rule, replacing "Markers" above:
- Body markers only. Ball, FP1 and FP2 are not on the body and are left out.
- Markers present in every included throw are kept, except a marker whose gaps inside the window would remove more than 2 percent of the throws, which is dropped (R5MC only).
- A throw with a gap in any kept marker inside its window is excluded.
- Throws whose end event comes before the start event (three throws) are excluded, since no window can be formed.
Result of the rule (counts only): N = 57 markers, 2,117 throws, 14 or 15 sessions for every participant and side.
Everything else in the plan is unchanged. Ball speed is taken from Trialinfo.xlsx.
