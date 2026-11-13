# Results in the Kyokushin karate data, Prediction 2 (5 November 2026)

Run once under plan_karate.md with the two changes recorded there (05:31 UTC, rate read from each file; 11:26 UTC, edge frames trimmed). Run at 11:28 UTC. Raw output: karate_result.txt. Per-recording and per-athlete values: karate_per_recording.csv, karate_per_athlete.csv. Figure: figures/fig_karate.png. Checksums of the 146 files used: data/checksums/karate_selected.sha256.

Data kept: all 146 recordings of S01 and S02, condition E01, trials T01 and T02, from all 37 athletes. Of the 74 athlete-technique cells, 72 hold two trials, and B0367 S01 and B0402 S02 hold one. No recording failed the validity rule. Edge frames were trimmed in three recordings, and no gap needed interpolation. One recording had a joint outside the 42 strongest links. Per-athlete H ranged from 3.608 to 3.734 bits (mean 3.689, SD 0.029), 94.8 to 98.1 percent of the ceiling log2 14 = 3.807 bits. Mean H was 3.668 bits for the reverse punch (S01) and 3.711 for the front kick (S02).

## Primary test, H falls as grade rises with age controlled: against
- Partial Spearman correlation of grade rank and H, age controlled: r = 0.331.
- One-sided p (predicted negative) = .976. Two-sided p = .048 from 10,000 permutations, seed 20261105.
- By the decision rule (against: r > 0 and two-sided p < .05) the test counts against the prediction. The p value sits close to the limit.

## Secondary results, reported as run and not used for the verdict
- (a) Spearman of grade rank and H without age: ρ = 0.115, two-sided p = .50.
- (b) Reverse punch alone: r = 0.354, two-sided p = .032 (against by the same rule). Front kick alone: r = 0.153, two-sided p = .37 (undecided).
- (c) H from the strength shares of the full absolute correlation network: mean 3.775 bits, r = 0.248, two-sided p = .14 (undecided).
- Grade and age: Spearman 0.866.

## Notes recorded with the results
- Post hoc check of the measure, made after the run and not part of the verdict (karate_posthoc_rate.txt). The four athletes recorded at 200 Hz include three dan grades. Resampling all 131 recordings at 250 Hz to 200 Hz changed H by +0.0001 bits on average (SD 0.0012, largest 0.013), against a spread of 0.029 bits between athletes, so the difference in rate does not account for the result.
- Descriptive: the Spearman correlation of H and age is −0.058.
- As in the How We Type test, H sits close to its ceiling. Here it spans 95 to 98 percent of the ceiling, wider than in How We Type (99.4 to 100) but still narrow.
