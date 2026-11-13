"""Reading A of the published figures of Hyman (1953) and Schmidt et al. (1979), used in Figures 1, 6 and 7.

Needs the published PDFs, which are not redistributed, and poppler's pdfimages.
Usage: python -I code/read_figures.py <Hyman 1953 PDF> <Schmidt et al. 1979 PDF>
Prints the values stored in data/read_from_figures. Reading B of Hyman's figure was made by eye on
gridded enlargements and is recorded in the CSVs, not computed here. Experiment I enters Figure 1B and Figures 6 and 7,
Experiments II and III enter Figure 7.
"""
import subprocess, sys, tempfile
from pathlib import Path
import numpy as np
from PIL import Image
from scipy import ndimage as ndi


def page_image(pdf, page):
    d = Path(tempfile.mkdtemp())
    subprocess.run(["pdfimages", "-f", str(page), "-l", str(page), "-png", pdf, str(d / "p")], check=True)
    return np.array(Image.open(sorted(d.glob("p-*.png"))[0]).convert("L")) < 128


def piecewise(p, P, V):
    """Linear between neighbouring ticks, extended linearly past the outer ticks."""
    P, V = np.asarray(P, float), np.asarray(V, float)
    if P[0] > P[-1]:
        P, V = P[::-1], V[::-1]
    i = int(np.clip(np.searchsorted(P, p) - 1, 0, len(P) - 2))
    return float(V[i] + (V[i + 1] - V[i]) * (p - P[i]) / (P[i + 1] - P[i]))


def filled_markers(B, radius=4, min_area=120, max_area=400):
    y, x = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    lab, n = ndi.label(ndi.binary_opening(B, structure=x * x + y * y <= radius * radius))
    out = []
    for k in range(1, n + 1):
        ys, xs = np.where(lab == k)
        if min_area <= len(xs) <= max_area:
            out.append((xs.mean(), ys.mean()))
    return out


def hyman(pdf):
    """Figure 1 (journal p. 192), participant G. C., Experiment I. Open circles located by the centroid of
    their enclosed white hole. Axis positions were measured on the 300 dpi page image."""
    B = page_image(pdf, 5)[500:1235, 150:930]
    rows, vals = [214.0, 391.0, 563.0, 734.5], [600, 400, 200, 0]
    cols, bits_ax = [148.0, 321.5, 495.5, 669.5], [0, 1, 2, 3]
    lab, n = ndi.label(~B)
    edge = set(np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]])))
    holes = []
    for k in set(range(1, n + 1)) - edge:
        ys, xs = np.where(lab == k)
        if 50 <= len(xs) <= 100 and abs(np.ptp(xs) - np.ptp(ys)) <= 2:
            holes.append((xs.mean(), ys.mean()))
    for b in [0.00, 1.00, 1.58, 2.00, 2.32, 2.58, 2.81, 3.00]:
        c = cols[0] + b * (cols[-1] - cols[0]) / 3
        cx, cy = min(holes, key=lambda h: abs(h[0] - c) + (0 if 140 < h[1] < 560 else 1e9))
        print(f"Hyman G.C. {b:.2f} bits: read {piecewise(cx, cols, bits_ax):.3f} bits, {piecewise(cy, rows, vals):.1f} ms")


HYMAN_II = [0.47, 0.72, 0.99, 1.39, 1.75, 2.16, 2.38, 2.75]       # Table 1, mean bits of each condition
HYMAN_III = [0.72, 0.92, 1.36, 1.79, 2.21, 1.00, 2.00, 2.81]      # Table 2, conditions 1 to 8


def hyman_exp23(B):
    """Figure 1, participant G. C., Experiments II (open squares) and III (open triangles), from the same page image
    and calibration as Experiment I. Each marker is located by the centroid of its enclosed white hole. Square holes
    are taller than wide, 55 to 80 pixels, filling more than 0.72 of their box. Triangle holes hold 15 to 35 pixels.
    Holes outside the plotting area (the key above row 160, the axis labels) are left out. For each condition the hole
    nearest its column is taken. Where a triangle overlaps another marker its hole splits in two, so fragments within
    14 pixels of the nearest one are merged by area."""
    rows, vals = [214.0, 391.0, 563.0, 734.5], [600, 400, 200, 0]
    cols = [148.0, 321.5, 495.5, 669.5]
    lab, n = ndi.label(~B)
    edge = set(np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]])))
    sq, tri = [], []
    for k in sorted(set(range(1, n + 1)) - edge):
        ys, xs = np.where(lab == k)
        w, h, area = np.ptp(xs) + 1, np.ptp(ys) + 1, len(xs)
        cx, cy = xs.mean(), ys.mean()
        if not (140 < cx and 160 < cy < 640):
            continue
        if 55 <= area <= 80 and area / (w * h) > 0.72 and h > w:
            sq.append((cx, cy, area))
        elif 15 <= area <= 35:
            tri.append((cx, cy, area))
    out = []
    for exp, bits, cand in (("II", HYMAN_II, sq), ("III", HYMAN_III, tri)):
        for b in bits:
            c = np.interp(b, [0, 1, 2, 3], cols)
            h0 = min(cand, key=lambda h: abs(h[0] - c))
            part = [h for h in cand if np.hypot(h[0] - h0[0], h[1] - h0[1]) < 14]
            wt = np.array([h[2] for h in part], float)
            cx = float(np.dot(wt, [h[0] for h in part]) / wt.sum())
            cy = float(np.dot(wt, [h[1] for h in part]) / wt.sum())
            rt = piecewise(cy, rows, vals)
            out.append((exp, b, cx, cy, rt))
            print(f"Hyman G.C. Exp {exp} {b:.2f} bits: hole at column {cx:.1f}, row {cy:.1f}"
                  f"{' (two fragments merged)' if len(part) > 1 else ''}, {rt:.1f} ms")
    return out


def schmidt(pdf):
    """Figures 8 (p. 426) and 9 (p. 427), distance errors, short movement times. The same nine means
    appear in both figures, so each figure gives one reading of the effective width."""
    B = page_image(pdf, 12)[300:1150, 420:1680]
    py = np.polyfit([29.5, 109.5, 184.5, 261.0, 336.5, 416.5, 496.0, 574.0, 650.5], range(9, 0, -1), 1)
    px = np.polyfit([208.0, 618.0, 1031.0], [10, 20, 30], 1)
    for x, y in sorted(filled_markers(B[:725, 57:]), key=lambda m: (round(m[0] / 50), m[1])):
        print(f"Schmidt Fig 8: A = {np.polyval(px, x + 57):.2f} cm, We = {np.polyval(py, y):.3f} mm")
    B = page_image(pdf, 13)[1700:2620, 380:1600]
    py = np.polyfit([12.5, 91.5, 167.0, 246.0, 323.0, 403.0, 480.5, 556.5, 636.5, 714.5], range(10, 0, -1), 1)
    px = np.polyfit([342.0, 603.0, 861.5, 1119.5], [50, 100, 150, 200], 1)
    for x, y in sorted(filled_markers(B[:788, 87:], min_area=150)):
        x += 87
        if 900 < x < 980 and 420 < y < 530:
            continue  # legend symbols
        print(f"Schmidt Fig 9: A/MT = {np.polyval(px, x):.1f} cm/s, We = {np.polyval(py, y):.3f} mm")


if __name__ == "__main__":
    hyman(sys.argv[1]); hyman_exp23(page_image(sys.argv[1], 5)[500:1235, 150:930]); schmidt(sys.argv[2])
