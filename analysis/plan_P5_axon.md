# Analysis plan for Prediction 5 in the Axon game data, fixed before the analysis (4 November 2026)

Claim tested: performance changes quickly and then slowly across practice (Equation 18). For a score that rises with skill, Equation 18 reads S(t) = S∞ − R(1 + t)^(−α), with α > 0.

Data: data_by_cookie.json from Stafford and Dewar (2014), retrieved at commit 58c163a895cf (registry D3). Each machine identity holds numbered plays, each with a score, a date and an hour. Counts seen before this plan (no scores examined): 854,064 identities, 3,099,705 plays, 5,706 identities with at least 30 plays, 422 with at least 100. About 2 percent of identities have gaps in their play numbers.

Sample: identities with at least 30 and at most 100 plays. The upper bound follows the original analysis code, which names 101 attempts as a sensible limit because the regularity of the curve breaks down for identities with very many plays (sup_machinechurn.py). Plays are ordered by their number, and the first 30 are used.

Measures per identity: block means of score over plays 1 to 10 (m1), 11 to 20 (m2) and 21 to 30 (m3); early gain g1 = m2 − m1; later gain g2 = m3 − m2; deceleration d = g1 − g2.

Primary test: the share of identities with d > 0, with a one-sided binomial test against one half and a 95 percent Wilson interval.
- Support: the share exceeds one half with one-sided p < .05 and the median early gain g1 is above zero.
- Against: the share is below one half with two-sided p < .05, or the median g1 is zero or below.
- Undecided: any other case.

Reported whatever the outcome: the share and its interval, the medians of g1, g2 and d, and, as a description of the average curve, the least-squares fit of Equation 18 to the mean score at plays 1 to 30, with a 95 percent interval for α from 1,000 bootstrap resamples of identities (seed 20261104). Averaging can make a curve look like a power function (Heathcote et al., 2000), which is why the primary test is made on each identity's own curve.
