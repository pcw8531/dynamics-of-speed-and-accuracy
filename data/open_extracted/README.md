# Values extracted from the open datasets

This folder holds every value the paper takes from the four open datasets of Supplemental Section S7, one table for
each analysis, as written by the analysis scripts under the plans in the analysis folder. The scripts write the same
files to analysis/, and the copies here are identical (data/checksums/open_extracted.sha256). The raw files of the
datasets are not copied. Each is fetched from its host and checked against the SHA-256 listed in data/REGISTRY.md and
data/checksums, so that anyone can extract the same values again.

## How We Type (Feit, Weir and Oulasvirta, 2016a, 2016b), registry D1, CC BY-NC 4.0

p2_per_typist.csv. One row for each of the 30 typists: uid (the typist in the dataset), median_iki_ms (median
inter-key interval of the typist), H_bits (averaged information of the correlation network among the 50 hand markers),
slope_ms_per_bit (slope of inter-key interval on the index of difficulty for same-finger moves between letter keys)
and group (faster or slower half by median interval). Script analysis/p2_test.py, plan analysis/plan_P2.md. The slopes
of the two halves in Table 3 of the main text are fitted with one intercept for each typist, as p2_result.txt reports.

p4_per_typist.csv. The same 30 typists with corrections_per_char, the corrections per character typed. Script
analysis/p4_howwetype.py, plan analysis/plan_P4_howwetype.md.

## Throwing practice (Ozkaya et al., 2018), registry D2, CC BY 4.0

throwing_per_throw.csv. One row for each throw kept by the plan: subject, session, session_no, side (dominant DP or
non-dominant NDP arm), trial, H_bits (averaged information of the network that links the pairs of markers whose
speeds correlate most strongly within the window of the throw, as the plan defines it) and ball_speed. throwing_per_session.csv gives the mean H and ball speed of each subject, side
and session, with n throws. Script analysis/throwing_test.py, plan analysis/plan_throwing.md.

## Kyokushin karate (Szczęsna, Błaszczyszyn and Pawlyta, 2020, 2021), registry D3, CC BY 4.0

karate_per_recording.csv. One row for each of the 146 C3D recordings selected by the plan: file, athlete, technique,
trial, rate (the frame rate of the file), frames, valid (share of valid frames), trimmed, max_gap, H_degree and H_strength
(averaged information of the joint network by degree and by strength), isolated_joints and used. karate_per_athlete.csv
gives each of the 37 athletes with H for each technique, the mean over recordings, and age, grade, rank of grade and
years of training from Table 1 of Szczęsna et al. (2021) (data/transcribed/szczesna2021_table1_participants.csv).
Script analysis/karate_test.py, plan analysis/plan_karate.md.

## Online game (Stafford and Dewar, 2014a, 2014b), registry D4

p5_mean_curve.csv. The mean score of the 5,293 players with 30 to 100 plays at each of plays 1 to 30 (mean_score), with
the least-squares fit of Equation 16 written for a score that rises (fit_eq18, named by the equation number of the
draft of its date). Script analysis/p5_axon.py, plan analysis/plan_P5_axon.md. The repository of the game states no
data licence, so only these 30 means are given here, and the file of plays is fetched by its hash.
