# Analysis plan: Prediction 4 at the level of joints, Kyokushin karate (registry C3)

Written on 4 November 2026 at 23:15 KST (14:15 UTC), before any file of this dataset has been downloaded or opened. Only the published
data descriptor has been read (Szczęsna, Błaszczyszyn & Pawlyta, 2021, Scientific Data 8, 13,
doi 10.1038/s41597-021-00801-5). Its Table 1 was transcribed and read twice, from the text layer and from the page
image, and the two readings agree (data/transcribed/szczesna2021_table1_participants.csv).

## Why this dataset and why joints

The law treats the elements of the motor system as its degrees of freedom. The two marker tests so far measured
coupling among skin markers, and many markers sit on one rigid body part, so part of their coupling is fixed by the
skeleton and cannot change with skill. Here the elements are joints. The release holds, for every recording, the joint
angles computed by Plug-In Gait from 39 markers at 250 Hz, for 37 athletes from 9th kyu to 4th dan. The theory of
Section S2.10 also says that H must be compared over one fixed set of elements. Every athlete has the same joints, so
the set is fixed by design.

## Prediction tested

Prediction 4. Lower averaged information of the coordination that carries a task goes with higher skill.
Here skill is the karate grade, so the law predicts that H falls as grade rises, once age is taken into account.

## Data used

- Techniques done by every athlete in the air condition (descriptor, Table 2): S01 Gyaku-Zuki (reverse punch) and
  S02 Mae-Geri (front kick). Kicks S03 to S05 are left out, because athletes who could not do them were not recorded,
  so using them would select by skill.
- Condition E01 (in the air) only, trials T01 and T02, every file named YYYY-MM-DD-CODE-S01-E01-T0J.c3d or
  YYYY-MM-DD-CODE-S02-E01-T0J.c3d.
- C3D files are read with py-c3d (github.com/EmbodiedCognition/py-c3d, commit b3060462bf35, MIT licence). Only the
  reader may be adapted to the layout of the released files. The rules below stay as written.

## Fixed rules

1. Elements, 14 joints: left and right hip, knee, ankle, shoulder, elbow and wrist, from the Plug-In Gait outputs
   LHipAngles, RHipAngles, LKneeAngles, RKneeAngles, LAnkleAngles, RAnkleAngles, LShoulderAngles, RShoulderAngles,
   LElbowAngles, RElbowAngles, LWristAngles, RWristAngles, and the spine and neck from LSpineAngles and LNeckAngles.
   The left and right spine and neck outputs describe the same joint with mirrored signs, so one of each is used.
   Segment angles relative to the room (pelvis, thorax, head, foot progression) are not joints and are not used.
   If one of the 14 outputs is absent from the release, that joint is dropped for every athlete, and this is decided
   from the parameter list of the first file read, before any value is computed.
2. Valid frames. A frame is invalid for a joint if any of its three angles is missing or exactly zero in all three
   components. A recording is used only if at least 80 percent of its frames are valid for all joints. Invalid frames
   inside a recording are filled by linear interpolation when the gap is at most 25 frames (0.1 s), and a recording
   with a longer gap is not used.
3. Signal. Each angle series is unwrapped (degrees), filtered with a fourth-order zero-lag Butterworth low-pass at
   15 Hz, and the angular speed of each joint is the norm of the first difference of its three angles times 250,
   in degrees per second. The whole recording is used, which holds three or four repetitions of the technique.
4. Coupling and H. For each recording, the absolute Pearson correlations between the 14 angular speed series form the
   coupling network. The 3N strongest pairs are kept (42 links for 14 joints, a mean degree of 6), the rule used in
   the throwing test and in the deposit of Park (2026). Ties at the cut are broken by a stable sort of pair index.
   H is the entropy in bits of the degree shares (Equations 6 and 7). Its ceiling is log2 14 = 3.807 bits.
5. Per athlete. H for a technique is the mean over that athlete's usable recordings of it. The athlete's H is the
   mean of the two technique values. An athlete with no usable recording of either technique is left out, and the
   number left out is reported.
6. Skill. Grade rank from Table 1 of the descriptor, 9th kyu = 1 up to 1st kyu = 9, then 1st dan = 10 up to
   4th dan = 13. Training years are not used, because two rows of Table 1 give training years equal to age
   (B0404, 29 and 29; B0405, 28 and 28), which is unlikely to be correct.

## Test and decision rule

Primary. Partial Spearman correlation between grade rank and H, controlling for age. Grade, H and age are ranked,
the ranks of grade and of H are each regressed on the rank of age, and r is the Pearson correlation of the two
residuals. Age is controlled because grade and age are strongly related in this sample (Spearman 0.87), since the
athletes include children of 10 to 14 years. The one-sided p value comes from 10,000 permutations of the H residuals,
seed 20261105.

- Support: r < 0 and one-sided p < .05.
- Against: r > 0 and two-sided p < .05.
- Undecided: otherwise.

Secondary, reported with the primary and not used for the verdict:
(a) the Spearman correlation of grade rank and H without control for age;
(b) the primary test for each technique separately;
(c) the primary test with H computed from strength shares of the full absolute correlation network, the measure of
    the How We Type test.

## Outputs

analysis/karate_result.txt (the run log), analysis/karate_per_recording.csv, analysis/karate_per_athlete.csv,
figures/fig_karate.py and its figure, and analysis/results_karate.md. The test is run once on the full data.
Any change to these rules after this plan is saved is written below, with its date and reason, before the run.

## Changes after saving

5 November 2026, 05:31 UTC, before any value was computed. Only the header and the channel labels of one file
(2017-01-31-B0367-S01-E01-T01.c3d) were read.
- The file header gives 200 frames per second, not the 250 Hz stated in the descriptor. The sampling rate is
  therefore read from the header of each file. The rules in time units stay as written: the low-pass cut at 15 Hz,
  gaps filled up to 0.1 s (20 frames at 200 Hz), and angular speed in degrees per second.
- All 14 joint outputs listed in rule 1 are present, so no joint is dropped.
- The figshare notice of 2 March 2023 corrects the names of three files, all in techniques S03 and S05, so none of
  them is used here.

5 November 2026, 11:26 UTC, before any value of H was computed. All 146 recordings that the plan selects (S01 and
S02, E01, T01 and T02, all 37 athletes) were read for a structure check that counted frames and gaps only. Two
athlete-technique cells hold one trial (B0367 S01 and B0402 S02), which rule 5 allows.
- 131 files are at 250 frames per second and 15 at 200 (all files of B0367 to B0370, recorded on 31 January 2017).
  The rate is read from each file, as recorded above, so the gap limit of 0.1 s is 25 or 20 frames.
- All 14 joint outputs are present in every file. Every recording has at least 99.2 percent of its frames valid for
  all 14 joints, and none has an invalid gap inside it.
- Three recordings begin or end with invalid frames, which rule 2 does not cover because they are not inside the
  recording: B0378 S01 T01 (1 leading frame), B0388 S01 T02 (18 leading) and B0401 S02 T02 (9 leading, 1 trailing).
  Such edge frames cannot be interpolated, so they are trimmed and the recordings are kept. No other rule changes.

