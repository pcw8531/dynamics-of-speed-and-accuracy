# Values read from published figures

These values were not printed as numbers in the source papers. They were read from the figures and are labelled "read from figure" wherever they are used. Hyman's values appear in Figure 1B, in the analyses of Hick's and Hyman's laws (Figures 6 and 7) and in Supplemental Figure S2. Schmidt et al.'s values appear in Figure 1D only.

## Hyman (1953), Figure 1, participant G. C., Experiment I

Source: Hyman, R. (1953). Stimulus information as a determinant of reaction time. Journal of Experimental Psychology, 45(3), 188-196. Figure 1, p. 192, first panel. The PDF is held by the author and is not redistributed.

What was read: the eight open circles of Experiment I, one to eight equally likely alternatives. The information of each condition is fixed by the text (0.00, 1.00, 1.58, 2.00, 2.32, 2.58, 2.81 and 3.00 bits), so only reaction time was read.

How: the axes were calibrated from the tick marks of the 300 dpi page image, piecewise between neighbouring ticks. Reading A located each circle by the centroid of its enclosed white hole (code/read_figures.py). Reading B was made by eye on fourfold enlargements with a 10 ms grid. The stored value is the mean of the two readings.

Checks: the two readings differ by 1.0 ms at most. The information read from the figure differs from the text by 0.011 bits at most, which checks the horizontal calibration. A straight line through the stored values gives r = .993, against the r = .991 that Hyman reports for this participant in Experiment I. Hyman's line for this participant over all three experiments, RT = 212 + 153 H with r = .985, is printed in the figure. Figure 1B draws the least-squares line through Experiment I.

## Hyman (1953), Figure 1, participant G. C., Experiments II and III

Source: as above, the same panel and page image.

What was read: the eight open squares of Experiment II (unequal probabilities) and the eight open triangles of Experiment III (sequential dependence). The information of each condition is fixed by Hyman's Tables 1 and 2 (data/transcribed/hyman1953_tables1_2.csv), so only reaction time was read; the column of each marker was also read, to check the horizontal calibration.

How: the same calibration as Experiment I. Reading A located each marker by the centroid of its enclosed white hole (code/read_figures.py, function hyman_exp23). Square holes are taller than wide and hold 55 to 80 pixels; triangle holes hold 15 to 35 pixels. For each condition the hole nearest its column is taken. The triangle at 2.81 bits overlaps the circle of Experiment I at the same information, which splits its hole in two, and the two parts are merged by area. Reading B was made by eye on 12-fold enlargements with a grid. The stored value is the mean of the two readings.

Checks: the two readings differ by 3.3 ms at most. The information read from the figure lies within 0.053 bits of Tables 1 and 2. Straight lines through the stored values give r = .987 for Experiment II and r = .973 for Experiment III, against the .989 and .972 that Hyman reports for this participant (p. 191). A line through all 24 conditions gives RT = 217.3 + 147.8 H with r = .985, against the line RT = 212 + 153 H with r = .985 printed in the figure.

## Hyman (1953), p. 194, the two kinds of light in Condition 3 of Experiment II

Transcribed rather than read, in data/transcribed/hyman1953_p194_components.csv: the observed and predicted mean times of the frequent light (13/16, 0.30 bits), of the rare lights (3/16 in all, 4.00 bits) and of the condition (0.99 bits), as printed. Weighting the printed component means by 13/16 and 3/16 gives 358.3 and 364.1 ms, within 3 ms of the printed condition means of 361 and 363 ms.

## Schmidt et al. (1979), Figures 8 and 9, distance errors

Source: Schmidt, R. A., Zelaznik, H., Hawkins, B., Frank, J. S., & Quinn, J. T., Jr. (1979). Motor-output variability: A theory for the accuracy of rapid motor acts. Psychological Review, 86(5), 415-451. Figure 8, p. 426, and Figure 9, p. 427. The PDF is held by the author and is not redistributed.

What was read: the nine mean effective target widths in the direction of movement, for distances of 10, 20 and 30 cm and movement times of 140, 170 and 200 ms, four participants. The same nine means appear in both figures, against distance in Figure 8 and against average velocity in Figure 9.

How: the axes were calibrated from the tick marks of the 300 dpi page images. Filled markers were located by their centroids after removing thin lines (code/read_figures.py). Figure 9 gives the first reading, with both coordinates as plotted. Figure 8 gives the second reading of the width, with distance fixed by the design.

Checks: the two readings of the width differ by 0.10 mm on average and by 0.16 mm at most, Figure 9 lying slightly lower throughout. Distances read from Figure 8 differ from the design values by 0.06 cm at most. Eight of the nine velocities plotted in Figure 9 lie within 6 cm/s of distance divided by movement time. The ninth, 30 cm at 170 ms, is plotted at 200.8 cm/s against a nominal 176.5 cm/s. The caption calls the horizontal axis the average movement velocity, so the values are kept as plotted. A straight line through the Figure 9 values gives We = 2.25 + 0.030 A/MT with r = .976, against the 2.12 + 0.033 A/MT and r = .97 reported in the text. The line drawn in the original figure has a slope of 0.033. Figure 1D shows the Figure 9 values with the line reported in the text.
