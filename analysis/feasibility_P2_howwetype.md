# Feasibility of the Prediction 2 test in How We Type (4 November 2026)

Purpose: check that the data can carry the test before the plan is written. Only counts and ranges were computed (analysis/feasibility_counts_P2.py). No slope, entropy or correlation was computed.

Data: Typing data Version 2.0, 30 typists, sentences condition (registry D1). The release also holds motion capture at 240 frames per second of 52 markers on both hands, with each keypress marked (Motion Capture.zip, larger than the 400 MB copy limit, not yet read).

Key coordinates: given in metres, with Escape at the origin. The key rows differ in z, so distances are three-dimensional. The median nearest-neighbour distance between letter keys is 19.0 mm, the standard key pitch.

Coordinate problems: the name "0" appears twice with two positions, and five positions are shared by several names ("exclam" and "1"; "quotedbl" and "2"; "equal", "9", "7" and "0"; "plus" and "question"; "underscore" and "minus"). The 29 letter keys (a to z with aring, adiaeresis, odiaeresis) each have one position of their own. The test therefore uses letter keys only.

Counts per typist, letter keys only, consecutive keystrokes within one sentence made by the same finger on two different keys (minimum, median, maximum):
- moves: 23, 76.5, 366 (3 typists have fewer than 30)
- distinct distances: 9, 18, 46
- distance between key centres: shortest 18.0 to 19.0 mm, longest 45 to 130 mm
- span of log2(2D/W): 1.25, 1.47, 2.80 bits

For comparison, counting every key with a unique position gives 34 to 554 moves and distances up to 186 mm, but those keys include the space bar, shift and others whose width differs from a letter key.

Key width: with one width W for all letter keys, log2(2D/W) = log2 D + log2(2/W), so W shifts the intercept only and does not change the slope.

Timer: 82 percent of inter-key intervals lie within 1.5 ms of a multiple of 15.625 ms (most common values 125, 187, 155, 156, 62, 93 ms). The clock therefore ticks in steps of about 15.6 ms. This adds noise of a few milliseconds to each interval but no bias to an average over many moves.

Conclusion: a Fitts-type slope over same-finger letter moves can be estimated, but the range of difficulty is narrow (about 1.5 bits for most typists) and some typists have few moves. Pooling moves within a group of typists, with an intercept for each typist, is more stable than one slope per typist. The plan must fix in advance the measure of averaged information, the grouping, the exclusions and the rule for judging the result.
