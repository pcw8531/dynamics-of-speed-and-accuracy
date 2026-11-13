# Results in How We Type, Predictions 2 and 4 (4 November 2026)

Both tests were run once, under plans fixed beforehand (plan_P2.md; plan_P4_howwetype.md, written at 07:06:11 UTC, before the Prediction 2 output first existed at 07:08:19 UTC). Raw output: p2_result.txt and p4_howwetype_result.txt. Per-typist values: p4_per_typist.csv. Figure: figures/fig_p2_p4_howwetype.png. A first launch at 07:05:47 UTC stopped before reading any data because a timing utility was missing, and produced no output.

## Prediction 2, verdict: against
- 2,641 moves, 30 typists, 50 markers.
- Faster half: b = 35.25 ms per bit (95 percent interval 23.24 to 47.90), H = 5.6324 bits, b/H = 6.26.
- Slower half: b = 61.78 ms per bit (41.34 to 70.80), H = 5.6344 bits, b/H = 10.97.
- ρ = 0.571, interval 0.369 to 0.943, which excludes 1.
- What the numbers show: the slope of the faster half is lower by 43 percent at the point estimates, with overlapping intervals, while H is the same to within 0.04 percent. Every typist's H lies between 5.6097 and 5.6433 bits, within 0.6 percent of the ceiling log2 50 = 5.6439, with a standard deviation of 0.0085 bits. Under b = βH with a common β, a slope ratio far from 1 needs an H ratio far from 1, and this H has almost no room to move.

## Prediction 4, verdict: undecided
- ρ(H, median inter-key interval) = 0.229, one-sided p = .112, interval −0.164 to 0.569. Predicted sign, not reliable.
- ρ(H, corrections per character) = −0.214, two-sided p = .257, interval −0.563 to 0.194. Opposite sign, not reliable.

## Limits recorded with the results
- The averaged information measured here, from the absolute correlations of all pairs of marker speeds, sits at the ceiling for every typist. A dense network of correlations among markers that share segments and hands spreads strength almost evenly, whereas the model concentrates links on a few elements. The measure was chosen in advance and is reported as chosen.
- Inter-key intervals come in clock steps of about 15.6 ms, so many median intervals tie at 125, 156 and 187 ms, and typists at the boundary of the two halves (155 and 156 ms) are not separable in speed.
- Per-typist slopes rest on 18 to 346 moves over about 1.5 bits of difficulty and range from −62 to 96 ms per bit.
