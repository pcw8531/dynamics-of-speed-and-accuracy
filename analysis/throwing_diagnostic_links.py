"""Diagnostic, not a test: which marker pairs form the 3N strongest links, in a fixed random sample of 200 analysed throws.
Usage: python -I analysis/throwing_diagnostic_links.py <folder with S1.zip ... S7.zip and Trialinfo.xlsx>"""
import sys, zipfile, re, collections, numpy as np, pandas as pd
from pathlib import Path
HERE = Path(__file__).resolve().parent
src = (HERE / "throwing_test.py").read_text()
src = src.split("# pass 1:")[0]          # helper functions only
exec(src)
t = pd.read_csv(HERE / "throwing_per_throw.csv")
samp = t.sample(200, random_state=20261104)
seg = lambda mk: re.sub(r"\d+$", "", mk)      # e.g. RTH1 -> RTH, LSK3 -> LSK
same = []; topseg = collections.Counter()
keep = None
for _, r in samp.iterrows():
    z = zipfile.ZipFile(D / f"S{r.subject}.zip")
    name = [n for n in z.namelist() if re.search(rf"/{r.session}/S{r.subject}_{r.side}_{int(r.trial):02d}_O\.csv$", n)][0]
    P = read_trial(z, name)
    if keep is None:
        keep = sorted(mk for mk in P if mk != "R5MC")
    P = {mk: P[mk] for mk in keep}
    names = sorted(P); S = []
    for mk in names:
        Q = filtfilt(b, a, P[mk], axis=0); S.append(np.linalg.norm(np.diff(Q, axis=0), axis=1) * FS)
    W = np.abs(np.corrcoef(np.array(S))); N = len(names); iu = np.triu_indices(N, 1)
    order = np.argsort(W[iu])[::-1][:3 * N]
    pairs = [(names[iu[0][i]], names[iu[1][i]]) for i in order]
    s = [seg(x) == seg(y) for x, y in pairs]; same.append(np.mean(s))
    for x, y in pairs:
        if seg(x) == seg(y): topseg[seg(x)] += 1
print("throws sampled", len(same), "share of the 3N links that join two markers on the same segment group: median %.3f, min %.3f, max %.3f" % (np.median(same), np.min(same), np.max(same)))
allpairs = [(x, y) for i, x in enumerate(keep) for y in keep[i+1:]]
print("share of all possible pairs that are same-group:", round(np.mean([seg(x) == seg(y) for x, y in allpairs]), 3))
print("most frequent same-group links:", topseg.most_common(10))
