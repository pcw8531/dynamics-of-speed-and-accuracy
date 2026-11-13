"""Feasibility counts for the Prediction 2 test in How We Type (Feit et al., 2016).

Counts and ranges only. No slope, entropy or correlation is computed here.
Usage: python -I analysis/feasibility_counts_P2.py <folder of *_matched.txt> <keyboard_flat_coordinates.csv>
Key coordinates are in metres in the release and are converted to millimetres. Distances are three-dimensional,
because the key rows differ in z.
"""
import sys, glob, numpy as np, pandas as pd
D, C = sys.argv[1], sys.argv[2]
c = pd.read_csv(C, sep="\t"); print("coord columns:", c.columns.tolist())
name = c.columns[0]
cnt = c[name].value_counts()
pos = c.groupby(["x", "y", "z"])[name].apply(list)
shared = set(n for l in pos if len(l) > 1 for n in l) | set(cnt[cnt > 1].index)
letters = [k for k in c[name] if (isinstance(k, str) and len(k) == 1 and k.isalpha()) or k in ("aring", "adiaeresis", "odiaeresis")]
print("letter keys in file:", len(letters), sorted(letters)); print("letters with shared position:", [k for k in letters if k in shared])
L = set(letters) - shared
xy = 1000 * c.drop_duplicates(name).set_index(name)[["x", "y", "z"]]
rows, q = [], []
for f in sorted(glob.glob(D + "/*_matched.txt")):
    t = pd.read_csv(f, sep="\t")
    k, fi, iki, st = t.key_symbol.values, t.finger.values, t.iki.values, t.stimulus.values
    q.append(iki[1:])
    m = (k[1:] != k[:-1]) & (fi[1:] == fi[:-1]) & (st[1:] == st[:-1]) & np.isin(k[1:], list(L)) & np.isin(k[:-1], list(L))
    i = np.where(m)[0] + 1
    d = np.linalg.norm(xy.loc[k[i]].values - xy.loc[k[i - 1]].values, axis=1) if len(i) else np.array([])
    rows.append((t.user_id.iloc[0], len(i), len(np.unique(np.round(d, 1))), d.min() if len(i) else np.nan, d.max() if len(i) else np.nan))
r = pd.DataFrame(rows, columns=["uid", "n_moves", "n_dist", "dmin", "dmax"])
print(r.describe().loc[["min", "50%", "max"]].round(1).to_string())
print("typists with fewer than 30 moves:", int((r.n_moves < 30).sum()))
a = np.concatenate(q).astype(float); a = a[a > 0]
res = np.abs(a / 15.625 - np.round(a / 15.625)) * 15.625
print("IKIs:", len(a), "share within 1.5 ms of a multiple of 15.625 ms:", round(float((res <= 1.5).mean()), 3))
print("most common IKI values:", pd.Series(a).value_counts().head(8).to_dict())

import itertools
P = xy.loc[sorted(L)].values
nn = [np.sort(np.linalg.norm(P - p, axis=1))[1] for p in P]
print("letter keys: nearest-neighbour distance median", round(float(np.median(nn)), 1), "mm")
r["ID_range"] = np.log2(2 * r.dmax / 19) - np.log2(2 * r.dmin / 19)
print("span of log2(2D/W) per typist (bits): min, median, max", r.ID_range.min().round(2), r.ID_range.median().round(2), r.ID_range.max().round(2))
