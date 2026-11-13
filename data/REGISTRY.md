# Data registry

Every source behind a number in the paper, what was taken from it, where it is held and under what terms. No synthetic
or simulated values are stored in data/. Simulation outputs are in sim_out/, beside the code that writes them.

## Published tables and figures, values held in data/

W1. Woodworth (1899), Tables I to IV (pp. 34 and 35). Average error in mm for the right hand at ten metronome rates,
eyes open and shut. File: transcribed/woodworth1899_tables_I-IV_right_hand.csv, with notes in the .md file of the same
name. Transcribed and read a second time independently; all 80 cells agree. Used in Figure 1A and in the account of
Woodworth's two components.

W2. Hyman (1953). Reaction times are given only in his Figure 1 (p. 192). Participant G. C., all 24 conditions of
Experiments I to III, read twice, from the centroids of the plotted markers and by eye on gridded enlargements; the
readings differ by 1.0 ms at most in Experiment I and 3.3 ms in Experiments II and III. Files:
read_from_figures/hyman1953_fig1_GC_experiment1.csv and read_from_figures/hyman1953_fig1_GC_experiments2_3.csv. The
lights and probabilities of each condition (Tables 1 and 2, pp. 190 and 191) and the times of the two kinds of light
in Experiment II, Condition 3 (p. 194) are in transcribed/hyman1953_tables1_2.csv and
transcribed/hyman1953_p194_components.csv. Used in Figure 1B, Figures 6 and 7 and Supplemental Figure S2.

W3. Schmidt, Zelaznik, Hawkins, Frank and Quinn (1979). The nine mean effective widths of the short-movement-time
experiment (10, 20 and 30 cm by 140, 170 and 200 ms, four participants), shown both in Figure 8 (p. 426) and in
Figure 9 (p. 427), so each figure gives one reading; the readings differ by 0.16 mm at most. File:
read_from_figures/schmidt1979_fig8_fig9_distance_We.csv. Used in Figures 1D and 9A.

W4. Fitts (1954), Tables 1 to 3 (pp. 385 to 388). Reciprocal tapping, disc transfer and pin transfer, with index of
difficulty, mean movement time, errors and index of performance. Files: transcribed/fitts1954_table1.csv,
fitts1954_table2.csv and fitts1954_table3.csv, with notes in the .md files. Checked against the printed index of
difficulty and index of performance. Table 1 is used in Figures 1C and 8 and Supplemental Figure S2.

W5. Mowbray and Rhoades (1959), Tables I to III (pp. 19 to 21). Table I gives the three-choice and six-choice times of
Leonard (1954) as tabled by them, Table II the two-choice and four-choice means and standard deviations of one highly
practised participant by trial, and Table III the error percentages. Files: transcribed/mowbray_rhoades1959_tables1_3.csv
and transcribed/mowbray_rhoades1959_table2.csv, with notes in the .md files, each read twice from the page images. Used
in Figure 10 and Supplemental Figure S3.

W6. Szczęsna, Błaszczyszyn and Pawlyta (2021), Table 1, the participants of the karate dataset (age, gender, weight,
height, training years, grade). File: transcribed/szczesna2021_table1_participants.csv, read twice. Used in the karate
test of Supplemental Section S7.

Cited, with no values taken: Hick (1952), Newell, Carlton and Hancock (1984), Crossman and Goodeve (1983).

## Open datasets, cited and not copied

Each is obtained from its host. The values the paper extracts from each are in data/open_extracted, with a
note on every column. The SHA-256 of every file used is given here or in data/checksums, so the same files
can be checked.

D1. How We Type (Feit, Weir and Oulasvirta, 2016a, 2016b). Zenodo, doi 10.5281/zenodo.4047775 (Typing data Version 2.0),
CC BY-NC 4.0. Files used: Typing data v2.0.zip (386,519 bytes,
c35a17c2c6ffaf0c8a044d60e99b25a7e3cedfe16bd24eecbb70304c3b47c85d), Readme.txt (6,827 bytes,
4d62ebde494ead2a5ddbef0fafe4a3a49ff579e10c8216eaf5340f471221252c), Background.xlsx (16,224 bytes,
5686f7f5ea787761187de65323f4297184f384354bbbe49e13a526d540d370bd), keyboard_flat_coordinates.csv (1,594 bytes,
b2c4648d31b47af49a4a08873f9a7e9d380f732807b926f7336e903a27bf1a30), and the 30 motion capture files listed in
checksums/howwetype_motion_capture.sha256. Used in the How We Type tests of Supplemental Section S7.

D2. Throwing practice (Ozkaya et al., 2018). figshare, doi 10.6084/m9.figshare.c.4017808, CC BY 4.0. Files used:
checksums/throwing_files.sha256. Used in the throwing test of Supplemental Section S7.

D3. Kyokushin karate (Szczęsna, Błaszczyszyn and Pawlyta, 2020, 2021). figshare, doi 10.6084/m9.figshare.c.4981073,
CC BY 4.0. The 146 C3D files selected by the plan are listed in checksums/karate_selected.sha256. They are read with
py-c3d, github.com/EmbodiedCognition/py-c3d, commit b3060462bf35, MIT licence. Used in the karate test of
Supplemental Section S7.

D4. Online game (Stafford and Dewar, 2014a, 2014b). github.com/tomstafford/axongame, commit 58c163a895cf. File used:
data_by_cookie.json (185,523,121 bytes, 0901ce85f8967be7ff4622faf25ff415b4541259ad2dc3b30ef9f29b324046bd). The
repository states no data licence, so the file is fetched by hash and cited, not redistributed. Used in the practice
curve test of Supplemental Section S7.

D5. Typing performance of 168,594 people (Dhakal, Feit, Kristensson and Oulasvirta, 2018), released for research use
with attribution at userinterfaces.aalto.fi/136Mkeystrokes. Only values reported in the paper are cited, in the
introduction, and no reanalysis is made.

## Released code and data of earlier work

P1. Park (2026), Zenodo, doi 10.5281/zenodo.21789488, github.com/pcw8531/control-entropy-coordinative-structures at
commit 0660a649f5e0. Its generators make the four structures of Figure 2 and those of Supplemental Figure S1, and it
holds the averaged information of the human musculoskeletal network (Murphy et al., 2018, CC BY-SA 4.0), 7.7369 bits,
compared with 7.9876 bits for its random null, cited in the text.

## Reported values cited without reanalysis

Each of these is cited for what it reports. Xu et al. (2009) was checked against the full text.

- Bailey and Chen (1988)
- Ericsson, Krampe and Tesch-Römer (1993)
- Galna and Sparrow (2006)
- Horak and Nashner (1986)
- Salmoni, Schmidt and Walter (1984)
- Sparrow and Newell (1998)
- Wulf (2013)
- Xu et al. (2009)
- Yang, Pan and Gan (2009)
