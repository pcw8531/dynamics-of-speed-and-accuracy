# Analysis plan for Prediction 4 in How We Type, fixed before the Prediction 2 output was read

Written on 4 November 2026 while the Prediction 2 test was running and before its output was opened. The Prediction 2 output shows the mean H of the faster and slower halves, which bears on this test, so the plan had to be fixed first.

Claim tested: lower averaged information of the coupling that carries the task goes with faster and more accurate performance (Prediction 4). Under b = βH, a typist with lower H has a flatter slope and so shorter intervals for the same moves, and the law ties accuracy gained with practice to the same concentration of coupling.

Data: How We Type, Version 2.0, sentences condition, all 30 typists (registry D1).

Measures, one per typist:
- H: the same value as in the Prediction 2 test (strength shares of the |r| network of the 50 hand-marker speed series, Equations 6 and 7), taken from analysis/p2_per_typist.csv.
- Speed: the median inter-key interval, as in the Prediction 2 test. Shorter is faster.
- Accuracy: error corrections per character, the number of BackSpace keystrokes divided by the number of characters in the final transcribed sentences, summed over the sentences. This follows the measure of corrections per character used for the 136 million keystrokes sample (Dhakal et al., 2018). Fewer is more accurate.

Statistic: Spearman rank correlation across the 30 typists, ρ(H, speed) and ρ(H, corrections), each with a one-sided p value for the predicted positive sign and a two-sided p value for the opposite sign.

Decision rule:
- Support: both correlations positive with one-sided p < .05.
- Against: either correlation negative with two-sided p < .05.
- Undecided: any other case. Each correlation is also reported on its own.

Reported whatever the outcome: both correlations with their p values and 95 percent bootstrap intervals (10,000 resamples of typists, seed 20261104), and a figure of speed and corrections against H. Any change to this plan after the data are read is reported as a change, with its reason.

## Result (run once, 4 November 2026)
Undecided. ρ(H, speed) = 0.229, one-sided p = .112; ρ(H, corrections) = −0.214, two-sided p = .257. See results_howwetype.md.
