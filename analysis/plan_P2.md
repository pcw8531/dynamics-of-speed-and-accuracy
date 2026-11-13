# Analysis plan for Prediction 2, fixed before any analysis (4 November 2026)

Claim tested: across people who share the same rate of command processing, the slope of a speed-accuracy law is proportional to the averaged information of the coupling that carries the movement, b = βH. For two groups this gives b_faster / b_slower = H_faster / H_slower (Equation 11).

Data: How We Type, Version 2.0 (Feit et al., 2016), sentences condition, all 30 typists. Typing logs give the keystrokes, keys and fingers. Motion capture gives the 3-D positions of 52 hand markers at 240 frames per second. Nothing has been computed from either beyond the counts in feasibility_P2_howwetype.md. The code that runs this plan is analysis/p2_test.py, and it runs once, on the full data.

## The slope b
- Moves: two consecutive keystrokes in one sentence, made by the same finger, on two different lowercase letter keys (the 29 letter keys, each with its own position).
- Excluded: the first keystroke of a sentence, intervals of 0 ms or less, intervals over 1,000 ms, and any move with a BackSpace immediately before or after it.
- Difficulty: ID = log2(2D/W), with D the 3-D distance between key centres and W = 19 mm. A single W shifts the intercept only.
- Count under these rules (counts only, 4 November 2026): 2,641 moves, 18 to 346 per typist, median 67.
- Model: inter-key interval = a_i + b_g ID, by least squares, with one intercept for each typist i and one slope for each group g.

## The averaged information H
- Frames: from the first to the last keystroke of each sentence.
- Gaps of up to 24 frames (0.1 s) are filled by linear interpolation. Longer gaps stay missing.
- Positions are low-pass filtered at 10 Hz (fourth-order Butterworth, forward and backward) within each unbroken stretch of at least 240 frames. Each marker's speed is the magnitude of its frame-to-frame velocity.
- Markers: those with valid speed in at least 80 percent of frames in every typist, so all typists share one set of N markers.
- Coupling: w_jk = |r|, the absolute correlation between the speed series of markers j and k over frames where both are valid (at least 240 such frames), for every pair, with no threshold. In the model a link added again to the same pair raises both degrees, so a degree is a sum of link weights; strength is its empirical counterpart.
- Shares and entropy: s_j = Σ_k w_jk, p_j = s_j / Σ s, H = −Σ p_j log2 p_j (Equations 6 and 7). Group H is the mean over its typists.

## Groups
Typing speed is each typist's median inter-key interval over all keystrokes, leaving out the first keystroke of each sentence and intervals of 0 ms or less or over 1,000 ms. The 15 typists with the shortest median interval form the faster group, the other 15 the slower group.

## Decision rule
- Statistic: ρ = (b_faster / b_slower) / (H_faster / H_slower). The law predicts ρ = 1.
- Interval: 95 percent percentile bootstrap, 10,000 resamples of typists within each group, seed 20261104. Slopes, H and ρ are recomputed in every resample.
- Support: the interval includes 1 and lies within 0.80 to 1.25, the range conventionally used for equivalence of ratios.
- Against: the interval excludes 1.
- Undecided: any other case, or when the 95 percent interval of either group's slope includes 0.

## Reported whatever the outcome
ρ with its interval; b and H for each group; β = b/H for each group; the number of moves and markers kept; one figure of each typist's slope against H. Any change to this plan after the data are read is reported as a change, with its reason.

## Structure check of the motion capture, before the run (4 November 2026, counts only)
- 30 files, one per typist, all matching the typing logs, 240 frames per second, 51,620 to 108,438 frames each. SHA-256 checksums in data/checksums/howwetype_motion_capture.sha256.
- Layout: tab-separated, with one line of capture settings before the header. The reader in p2_test.py was adapted to this. No rule changed.
- No missing values in any file, so the gap rule has no effect. Exact zeros are below 0.002 percent of values.
- One typist (419655) has 50 markers, without the two thumb base markers (L_T1, R_T1). Under the marker rule these two are left out for everyone, so N = 50.
- One typist (120778) has motion capture for 34 of 50 sentences (882 of 1,301 keypresses). H for this typist uses those frames.

## Result (run once, 4 November 2026)
Against. ρ = 0.571, 95 percent interval 0.369 to 0.943. See results_howwetype.md.
