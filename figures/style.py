"""Type, palette and colour scale shared by every figure script, the same as in fig2_mechanism.py.

The palette and the diverging scale are those of Park (2026, Figures 1 to 6), unchanged. The scale reads
log2(p_i N): blue below an even share, pale yellow at exactly 1/N, yellow to red above it, on +-2.716 bits.
"""
import logging
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, rgb_to_hsv, hsv_to_rgb

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
mpl.rcParams.update({"font.family": ["Times New Roman", "Liberation Serif", "DejaVu Serif"], "font.size": 7.5,
                     "mathtext.fontset": "stix", "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.major.size": 1.8, "ytick.major.size": 1.8,
                     "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "axes.labelsize": 7.2, "savefig.dpi": 600,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.unicode_minus": True,
                     "legend.frameon": True, "legend.facecolor": "white", "legend.edgecolor": "none",
                     "legend.framealpha": 1.0, "legend.fancybox": False, "legend.borderpad": 0.2})
# Type and line weight. By default the text, numbers and mathematics of every plot are set in the regular serif of the
# text (Liberation Serif has the metrics of Times New Roman) and the trend lines are drawn narrow, at NARROW times the
# width given in each script; only the panel letters are bold. FIG_WEIGHT=bold sets bold type and thicker trend lines
# instead, kept as an option. Symbols the serif lacks fall back to STIX.
BOLD = os.environ.get("FIG_WEIGHT", "plain") == "bold"
WEIGHT = "bold" if BOLD else "normal"                                # weight of every text but the panel letters
_W = ":bold" if BOLD else ""
mpl.rcParams.update({"font.weight": WEIGHT, "axes.labelweight": WEIGHT, "axes.titleweight": WEIGHT,
                     "mathtext.fontset": "custom", "mathtext.rm": "Liberation Serif" + _W,
                     "mathtext.it": "Liberation Serif:italic" + _W, "mathtext.bf": "Liberation Serif:bold",
                     "mathtext.sf": "Liberation Serif" + _W, "mathtext.fallback": "stix"})
NARROW = 0.75                                                        # trend lines at 0.75 of their given width
THICK = 1.4 if BOLD else NARROW                                     # factor applied to the trend lines


def bolden(fig, factor=THICK, lo=0.95):
    """Set the width of the trend lines of every plot, once per figure, narrower by NARROW (or thicker by 1.4 with
    FIG_WEIGHT=bold): each line of an axes with ticks whose width is at
    least lo points (grid lines, reference rules and task drawings are thinner or have no ticks), and the matching
    lines of its legend."""
    if getattr(fig, "_bold_done", False):
        return
    fig._bold_done = True
    for ax in fig.axes:
        if not (len(ax.get_xticks()) or len(ax.get_yticks())):
            continue
        for ln in ax.get_lines():
            if ln.get_linestyle() not in ("None", "", " ") and ln.get_linewidth() >= lo:
                ln.set_linewidth(ln.get_linewidth() * factor)
        leg = ax.get_legend()
        for h in (getattr(leg, "legend_handles", None) or []) if leg is not None else []:
            if hasattr(h, "get_linestyle") and hasattr(h, "set_linewidth") and h.get_linestyle() not in ("None", "", " ") \
                    and h.get_linewidth() >= lo:
                h.set_linewidth(h.get_linewidth() * factor)

INK, GREY, GUIDE, RULE, TRACK = "#1A1A1A", "#8A8A8A", "#455A64", "#8C989E", "#E3E7E9"
LINK, EDGE = "#C3CBD0", "#37474F"
TONE = {"regular": "#B9C4C9", "small-world": "#93A8B1", "random": "#6E8C99", "scale-free": "#48707F"}
_HOT = plt.get_cmap("hot_r")(np.linspace(0.06, 0.94, 256))
_COOL = plt.get_cmap("Blues")(np.linspace(0.92, 0.06, 256))
CMAP = LinearSegmentedColormap.from_list("hot_div", np.vstack([_COOL, _HOT]))
VLIM = 2.7157337518864986
CNORM = Normalize(-VLIM, VLIM)


def share_colour(bits):
    """Colour of an element that sits `bits` = log2(p_i N) above (positive) or below (negative) an even share."""
    return CMAP(CNORM(bits))


def outline(rgba):
    """The colour of a block drawn in outline: the hue of its colour on the scale of Figure 2, made deep enough to
    show as a line on white (the pale yellow of an even share becomes a deep yellow). Figures 7 and 9."""
    h, s_, v = rgb_to_hsv(np.asarray(rgba[:3], float))
    if 0.12 < h < 0.19:                                              # yellow is drawn as gold, which shows on white
        h = 0.135
    s_, v = max(s_, 0.85), min(v, 0.92)
    rgb = hsv_to_rgb([h, s_, v])
    while 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] > 0.62:
        v *= 0.95; rgb = hsv_to_rgb([h, s_, v])
    return rgb


def demand_colour(x, lo, hi):
    """Colour of a task condition on the blue side of the same scale: the easiest condition (x = lo) takes the pale
    yellow of an even share, and the colour deepens toward blue as the demand of the task grows, the way the
    information of an element grows in Figure 2. Used for the data markers of Figures 1 and 6 to 8."""
    f = np.clip((np.asarray(x, float) - lo) / (hi - lo), 0.0, 1.0)
    return CMAP(CNORM(-VLIM * 0.95 * f))


class Sheet:
    """A figure laid out in inches, with panel letters set right before their titles and the pair centred."""

    def __init__(self, fw, fh):
        self.fw, self.fh = fw, fh
        self.fig = plt.figure(figsize=(fw, fh))
        save = self.fig.savefig

        def savefig(*a, **kw):                                       # trend line widths set before the first save
            bolden(self.fig)
            return save(*a, **kw)
        self.fig.savefig = savefig

    def ax(self, x, y, w, h):
        return self.fig.add_axes([x / self.fw, y / self.fh, w / self.fw, h / self.fh])

    def text(self, x, y, s, **kw):
        return self.fig.text(x / self.fw, y / self.fh, s, **kw)

    def head(self, letter, title, x_centre, y=None, gap=0.06):
        rend = self.fig.canvas.get_renderer(); y = self.fh - 0.14 if y is None else y
        widths = []
        for txt, kw in ((letter, dict(fontsize=9, fontweight="bold")), (title, dict(fontsize=8))):
            t_ = self.fig.text(0, 0, txt, **kw); widths.append(t_.get_window_extent(rend).width / self.fig.dpi); t_.remove()
        x0 = x_centre - (widths[0] + gap + widths[1]) / 2
        a = self.text(x0, y, letter, fontsize=9, fontweight="bold", ha="left", va="baseline")
        b = self.text(x0 + widths[0] + gap, y, title, fontsize=8, ha="left", va="baseline")
        self.heads = getattr(self, "heads", {}); self.heads[letter] = [a, b]

    def head_over(self, ax, letter, title, gap=0.06):
        bb = ax.get_position()
        self.head(letter, title, (bb.x0 + bb.x1) / 2 * self.fw, gap=gap)


def inset(ax, bounds, frame=False):
    """A drawing of the task inside a panel, without a frame unless frame=True."""
    a = ax.inset_axes(bounds)
    a.set_xlim(0, 10); a.set_ylim(0, 6); a.set_xticks([]); a.set_yticks([])
    for sp in a.spines.values():
        sp.set_visible(frame); sp.set_linewidth(0.5); sp.set_color(RULE)
    a.patch.set_facecolor("white"); a.patch.set_alpha(1.0); a.set_zorder(2.5)
    return a


# ---------------------------------------------------------------- ternary panels (Figures 4B and 5B)
R3 = np.sqrt(3.0) / 2.0


def tern_xy(c):
    """Shares (top, lower left, lower right), each row summing to 1, to the plane. Figures 4B and 5B pass the shares
    of the chance hubs, the peripheral elements and the fit elements in that order."""
    c = np.atleast_2d(np.asarray(c, float))
    return c[:, 2] + 0.5 * c[:, 0], R3 * c[:, 0]


def tern_frame(ax, corners=("chance hubs", "peripheral elements", "fit elements"), grid=(0.25, 0.5, 0.75),
               all_ticks=False, dotted=False, tick_fs=5.6, tick_color=None, corner_fs=6.4, inside=True, ticks=None,
               edge_titles=None, off=0.17, top=0.12, outline="solid", slide=0.0):
    """The triangle of shares (top, lower left, lower right) with guide lines for each share. With ticks, each share is
    read on the edge that runs to its own corner, the top share up the left edge, the lower left share leftward along
    the base and the lower right share down the right edge, each value with a short mark. With dotted, the guide lines
    are thin dotted lines as in Figure 2B, and with inside=False there are only short marks on the edges. edge_titles
    names the three shares at the middle of the left edge, the base and the right edge, each with an arrow toward the
    larger share; without it the names in corners are set at the corners. outline="dotted" draws the triangle itself
    as a dotted line (Figures 4B and 5B, with grid=() for no guide lines), and slide moves the titles of the left and
    right edges down along their edges, away from what is drawn above the triangle."""
    tick_color = GREY if tick_color is None else tick_color
    ax.set_xlim(-0.08, 1.08); ax.set_ylim(-0.20 if all_ticks else -0.16, R3 + top); ax.set_aspect("equal")
    ax.axis("off")
    if outline == "dotted":
        for (x0, y0), (x1, y1) in (((0, 0), (1, 0)), ((1, 0), (0.5, R3)), ((0.5, R3), (0, 0))):   # each side on its own
            ax.plot([x0, x1], [y0, y1], color=EDGE, lw=0.8, zorder=2, ls=(0, (0.1, 1.9)), dash_capstyle="round")
    else:
        ax.plot([0, 1, 0.5, 0], [0, 0, R3, 0], color=EDGE, lw=0.7, zorder=2, solid_joinstyle="round")
    gkw = dict(color=GRIDC, lw=0.45, ls=(0, (1.0, 1.6))) if dotted else dict(color=TRACK, lw=0.5)
    L = 0.025
    for g in grid:
        if inside:
            for a, b in (((g, 0, 1 - g), (g, 1 - g, 0)), ((0, g, 1 - g), (1 - g, g, 0)), ((0, 1 - g, g), (1 - g, 0, g))):
                x, y = tern_xy([a, b])
                ax.plot(x, y, zorder=1, **gkw)
        else:                                                        # short ticks on the edges in place of guide lines
            for c_, (dx, dy) in (([g, 1 - g, 0], (-1.0, 0.0)), ([0, g, 1 - g], (0.5, -R3)), ([1 - g, 0, g], (0.5, R3))):
                if not all_ticks and dx != -1.0:
                    continue
                x, y = tern_xy(c_)
                ax.plot([x[0], x[0] + L * dx], [y[0], y[0] + L * dy], color=EDGE, lw=0.6, zorder=2)
        if ticks is None:
            x, y = tern_xy([g, 1 - g, 0])
            ax.text(x[0] - 0.035, y[0], f"{g:g}", fontsize=tick_fs, ha="right", va="center", color=tick_color)
            if all_ticks:
                x, y = tern_xy([0, g, 1 - g])                        # the lower left share, along the base
                ax.text(x[0] + 0.018, y[0] - 0.022, f"{g:g}", fontsize=tick_fs, ha="center", va="top", color=tick_color)
                x, y = tern_xy([1 - g, 0, g])                        # the lower right share, along the right edge
                ax.text(x[0] + 0.03, y[0] + 0.03, f"{g:g}", fontsize=tick_fs, ha="left", va="center", color=tick_color)
    for v in (ticks or ()):                                          # ticks with short marks, the ends included
        L2 = 0.022
        x, y = tern_xy([v, 1 - v, 0]); x, y = x[0], y[0]             # the top share, up the left edge
        ax.plot([x, x - L2], [y, y], color=EDGE, lw=0.6, zorder=2)
        ax.text(x - L2 - 0.012, y, f"{v:.1f}", fontsize=tick_fs, ha="right", va="center", color=tick_color)
        x, y = tern_xy([0, v, 1 - v]); x, y = x[0], y[0]             # the lower left share, leftward along the base
        ax.plot([x, x + 0.5 * L2], [y, y - R3 * L2], color=EDGE, lw=0.6, zorder=2)
        ax.text(x + 0.5 * L2 + 0.006, y - R3 * L2 - 0.008, f"{v:.1f}", fontsize=tick_fs, ha="center", va="top",
                color=tick_color)
        x, y = tern_xy([1 - v, 0, v]); x, y = x[0], y[0]             # the lower right share, down the right edge
        ax.plot([x, x + 0.5 * L2], [y, y + R3 * L2], color=EDGE, lw=0.6, zorder=2)
        ax.text(x + 0.5 * L2 + 0.012, y + R3 * L2, f"{v:.1f}", fontsize=tick_fs, ha="left", va="center",
                color=tick_color)
    if edge_titles:                                                  # each title at the middle of its own edge, outside the
        t_left, t_base, t_right = edge_titles                        # 0.5 tick, with an arrow toward the larger share
        ax.text(0.25 - off * R3 - 0.5 * slide, R3 / 2 + off * 0.5 - R3 * slide, t_left + r" $\rightarrow$",
                fontsize=corner_fs, ha="center",
                va="center", rotation=60, rotation_mode="anchor")
        ax.text(0.5, -off + 0.03, r"$\leftarrow$ " + t_base, fontsize=corner_fs, ha="center", va="center")
        ax.text(0.75 + off * R3 + 0.5 * slide, R3 / 2 + off * 0.5 - R3 * slide, t_right + r" $\rightarrow$",
                fontsize=corner_fs, ha="center",
                va="center", rotation=-60, rotation_mode="anchor")
        return
    yc = -0.10 if all_ticks else -0.05
    ax.text(0.5, R3 + 0.04, corners[0], fontsize=corner_fs, ha="center", va="bottom")
    ax.text(-0.02, yc, corners[1], fontsize=corner_fs, ha="left", va="top")
    ax.text(1.02, yc, corners[2], fontsize=corner_fs, ha="right", va="top")


GRIDC = "#B0BAC0"
CLEAR = dict(bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"))  # text set over a grid line


def tb_axes(ax, xticks=None, yticks=None, xlim=None, ylim=None, grid="both", twin=None, twin_ticks=None):
    """Lines above and below the plot only, ticks at the minimum, middle and maximum of each axis, and thin dotted grid
    lines at the ticks, as in Figure 2B. twin is a second y axis on the right, given its own three ticks."""
    for a in [ax] + ([twin] if twin is not None else []):
        a.spines["left"].set_visible(False); a.spines["right"].set_visible(False)
        a.spines["top"].set_visible(True); a.spines["bottom"].set_visible(True)
        for sp in ("top", "bottom"):
            a.spines[sp].set_linewidth(0.6); a.spines[sp].set_color(INK); a.spines[sp].set_zorder(6)  # above a legend
    if xlim is not None: ax.set_xlim(*xlim)
    if ylim is not None: ax.set_ylim(*ylim)
    if xticks is not None: ax.set_xticks(xticks)
    if yticks is not None: ax.set_yticks(yticks)
    if twin is not None:
        twin.set_ylim(*twin.get_ylim())
        if twin_ticks is not None: twin.set_yticks(twin_ticks)
    from matplotlib.ticker import FuncFormatter
    def fmt_for(ticks):                                              # the same number of decimals on one axis
        dec = max((len(f"{v:g}".split(".")[1]) if "." in f"{v:g}" else 0) for v in ticks)
        return FuncFormatter(lambda v, _: f"{v:.{dec}f}".replace("-", "\u2212"))
    if xticks is not None: ax.xaxis.set_major_formatter(fmt_for(xticks))
    if yticks is not None: ax.yaxis.set_major_formatter(fmt_for(yticks))
    if twin is not None and twin_ticks is not None: twin.yaxis.set_major_formatter(fmt_for(twin_ticks))
    if grid:
        ax.grid(True, axis="both" if grid == "both" else grid, color=GRIDC, lw=0.4, ls=(0, (1.0, 1.6)), zorder=0)
        ax.set_axisbelow(True)
    return ax


def save_panels(sheet, panels, stem, pad=0.05, rows=None):
    """Save each panel as its own image, cut from the finished figure so that it matches the figure exactly.
    panels maps a letter to the axes and texts of that panel (its heading is added), or to a box in inches. Everything
    else in the figure is hidden while a panel given by its artists is saved, so a panel that reaches into the space
    of its neighbour is cut cleanly. The crop of each panel spans what remains across, and from the top to the bottom
    of all the panels of its row, so that the panels of a row have one height and line up when they are placed side by
    side. rows lists the letters of each row (by default all panels form one row). Files are written as
    <stem>_<letter>.png."""
    from matplotlib.transforms import Bbox
    fig = sheet.fig; fig.canvas.draw(); rend = fig.canvas.get_renderer()
    inv = fig.dpi_scale_trans.inverted(); out = []
    every = [a for a in list(fig.axes) + list(fig.texts) + list(fig.artists) + list(fig.lines) + list(fig.patches)
             if a.get_visible()]
    boxes = {}
    for letter, items in panels.items():                             # the extent of each panel on its own
        if isinstance(items, tuple):
            continue
        keep = list(items) + list(getattr(sheet, "heads", {}).get(letter, []))
        boxes[letter] = (keep, Bbox.union([a.get_tightbbox(rend) for a in keep]).transformed(inv).padded(pad))
    for row in (rows or [list(boxes)]):                              # one height for the panels of a row
        row = [r for r in row if r in boxes]
        if not row:
            continue
        y0, y1 = min(boxes[r][1].y0 for r in row), max(boxes[r][1].y1 for r in row)
        for r in row:
            bb = boxes[r][1]; boxes[r] = (boxes[r][0], Bbox.from_extents(bb.x0, y0, bb.x1, y1))
    for letter, items in panels.items():
        f = f"{stem}_{letter}.png"
        if isinstance(items, tuple):                                 # a box in inches, (x0, y0, x1, y1), nothing hidden
            for a in every:
                a.set_visible(True)
            fig.savefig(f, bbox_inches=Bbox.from_extents(*items)); out.append(f)
            continue
        keep, bb = boxes[letter]
        for a in every:
            a.set_visible(any(a is k for k in keep))
        fig.savefig(f, bbox_inches=bb)
        out.append(f)
    for a in every:
        a.set_visible(True)
    return out


def by_region(sheet, regions):
    """The axes and figure texts whose centres fall in each region, (x0, y0, x1, y1) in inches, for save_panels."""
    fig = sheet.fig; fig.canvas.draw(); rend = fig.canvas.get_renderer(); inv = fig.dpi_scale_trans.inverted()
    heads = [t for v in getattr(sheet, "heads", {}).values() for t in v]
    out = {k: [] for k in regions}
    for a in list(fig.axes) + list(fig.texts) + list(fig.artists):
        if not a.get_visible() or any(a is h for h in heads):
            continue
        bb = a.get_tightbbox(rend).transformed(inv); xc, yc = (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2
        for k, (x0, y0, x1, y1) in regions.items():
            if x0 <= xc < x1 and y0 <= yc < y1:
                out[k].append(a); break
    return out


# ---------------------------------------------------------------- contours, points and networks of Figures 4B and 5B
WARM = LinearSegmentedColormap.from_list("warm", [share_colour(v) for v in np.linspace(0.0, VLIM, 256)])


def h_colour(H, ceil):
    """Colour of an averaged information H on the warm side of the scale of Figure 2: pale yellow at the ceiling,
    deepening toward red as H falls (the colour log2 N - H bits has there)."""
    return WARM(np.clip((ceil - np.asarray(H, float)) / VLIM, 0, 1))


def tern_hfield(ax, sizes, levels, ceil, label_at=None, lw=0.8, fs=6.4, fill=True):
    """The averaged information that a composition would have if each group shared its coupling evenly among its own
    elements, H = -sum x_g log2 x_g + sum x_g log2 n_g, an upper bound on H at that composition, drawn as filled
    bands between the levels, each band in the colour that log2 N - H has at its middle on the warm side of the scale of
    Figure 2, with the levels as labelled lines. With fill=False there are no bands, and each line takes the colour of
    its own level, with a thin dark rim, on white."""
    import matplotlib.tri as mtri
    n = 160; pts = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            pts.append((i / n, j / n, 1 - i / n - j / n))
    c = np.clip(np.array(pts), 0, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -np.nansum(np.where(c > 0, c * np.log2(c), 0.0), 1)
    Hb = ent + c @ np.log2(np.asarray(sizes, float))
    x, y = tern_xy(c)
    tri = mtri.Triangulation(x, y)
    if fill:
        bounds = [Hb.min() - 1e-6] + list(levels) + [Hb.max() + 1e-6]
        cols = [h_colour(0.5 * (a + b), ceil) for a, b in zip(bounds[:-1], bounds[1:])]
        ax.tricontourf(tri, Hb, levels=bounds, colors=cols, zorder=0.5)
        cs = ax.tricontour(tri, Hb, levels=levels, colors=[EDGE], linewidths=lw, zorder=1.6)
    else:                                                            # each line in the colour of its level, thin dark rim
        import matplotlib.patheffects as pe
        cs = ax.tricontour(tri, Hb, levels=levels, colors=[h_colour(v, ceil) for v in levels], linewidths=lw, zorder=1.6)
        cs.set(path_effects=[pe.Stroke(linewidth=lw + 0.9, foreground=EDGE), pe.Normal()])
    ax.clabel(cs, fmt=lambda v: f"{v:g}", fontsize=fs, inline=True, inline_spacing=2.0, manual=label_at, colors=INK)
    return cs


def tern_dots(ax, c, every=100, ms=2.6, zorder=4):
    """A mean path as plain points every `every` steps, dark with a thin white rim so that points that overlap stay apart."""
    c = np.asarray(c)
    k = np.r_[np.arange(0, len(c) - 1, every), len(c) - 1]
    x, y = tern_xy(c[k])
    ax.plot(x[1:-1], y[1:-1], "o", ls="", ms=ms, color=INK, mec="white", mew=0.3, zorder=zorder)
    ax.plot(x[0], y[0], "o", ms=ms + 1.0, color="white", mec=INK, mew=0.8, zorder=zorder + 1)      # start, step 0
    ax.plot(x[-1], y[-1], "D", ms=ms + 1.2, color=INK, mec="white", mew=0.3, zorder=zorder + 2)   # finish, last step
    return x, y


def loop_layout(n, seed=5):
    """A smooth closed contour through a ring, so that it reads as a chain and not as a circle (Figure 2A)."""
    rng = np.random.default_rng(seed); th = 2 * np.pi * np.arange(n) / n
    r = 1 + sum(a * np.cos(m * th + ph) for m, a, ph in zip((2, 3, 4), (0.26, 0.18, 0.08), rng.uniform(0, 2 * np.pi, 3)))
    xy = np.c_[r * np.cos(th), r * np.sin(th)]; xy -= xy.mean(0); xy /= np.abs(xy).max()
    return xy * 0.95


def tree_layout(T, w, gap=225.0):
    """Rings by distance from the element with the largest share, each subtree given an angular span in proportion to
    its size, turned so that a gap between subtrees faces `gap` degrees (the layout of Figure 2D)."""
    root = max(T.nodes(), key=lambda z: (w[z], T.degree(z), -z))
    parent, order, depth, seen, queue = {root: None}, [root], {root: 0}, {root}, [root]
    while queue:
        u = queue.pop(0)
        for v in sorted(T.neighbors(u), key=lambda z: (-w[z], -T.degree(z), z)):
            if v not in seen:
                seen.add(v); parent[v] = u; depth[v] = depth[u] + 1; order.append(v); queue.append(v)
    kids = {n: [] for n in order}
    for n in order[1:]:
        kids[parent[n]].append(n)
    size = {}
    for n in reversed(order):
        size[n] = max(1, sum(size[c] for c in kids[n])) if kids[n] else 1
    span = {root: (0.0, 2 * np.pi)}
    for n in order:
        a0, a1 = span[n]; tot = sum(size[c] for c in kids[n]) or 1; acc = a0
        for c in kids[n]:
            wd = (a1 - a0) * size[c] / tot; span[c] = (acc, acc + wd); acc += wd
    dmax = max(depth.values()) or 1
    ang = {n: 0.5 * sum(span[n]) for n in order}; rad = {n: np.sqrt(depth[n] / dmax) for n in order}
    first = sorted(ang[c] for c in kids[root]) or [0.0]
    gaps = [((first[(i + 1) % len(first)] - first[i]) % (2 * np.pi), first[i]) for i in range(len(first))]
    wid, start = max(gaps); turn = np.deg2rad(gap) - (start + wid / 2)
    return {n: np.array([rad[n] * np.cos(ang[n] + turn), rad[n] * np.sin(ang[n] + turn)]) * 0.80 for n in order}


def hub_layout(A, k, p, size, seed=3):
    """Positions in the style of the hub network of Figure 2D. The hubs, elements holding at least 1.8 times an even
    share (at most four, at least the strongest one), sit near the middle, and every other element hangs on the hub
    with which it shares most links, fanned outward from that hub on a spoke that is longer the less coupling it
    holds. Elements without links sit loosely below. Returns the positions and the links drawn, each element's tie to
    its hub and the ties among the hubs."""
    n = len(k); rng = np.random.default_rng(seed)
    order = np.argsort(-k, kind="stable")
    hubs = [int(i) for i in order[:4] if p[i] >= 1.8 / n] or [int(order[0])]
    h = len(hubs)
    rad = lambda q: np.sqrt(size * (3 + 560 * q) / np.pi) / 13.25   # node radius in units of the inset
    rh = 0.0 if h == 1 else max(0.26, (2 * rad(p[hubs].max()) + 0.08) / (2 * np.sin(np.pi / h)))
    ang_h = {hb: np.deg2rad(90 + 360 * j / h + (30 if h == 3 else 0)) for j, hb in enumerate(hubs)}
    pos = {hb: rh * np.array([np.cos(ang_h[hb]), np.sin(ang_h[hb])]) for hb in hubs}
    leaves = [int(i) for i in order if k[i] > 0 and int(i) not in hubs]
    home = {i: max(hubs, key=lambda hb: (A[i, hb], k[hb])) for i in leaves}
    Iinfo = -np.log2(np.where(p > 0, p, 1.0))
    lo, hi = (Iinfo[leaves].min(), Iinfo[leaves].max()) if leaves else (0.0, 1.0)
    for hb in hubs:
        mine = [i for i in leaves if home[i] == hb]
        if not mine:
            continue
        span = np.deg2rad(min(330 if h == 1 else 360 / h + 90, 40 + 24 * len(mine)))
        centre = ang_h[hb] if h > 1 else np.deg2rad(90)
        for j, i in enumerate(mine):
            th = centre - span / 2 + span * (j + 0.5) / len(mine) + rng.uniform(-0.05, 0.05)
            f = 0.0 if hi == lo else (Iinfo[i] - lo) / (hi - lo)
            L = rad(p[hb]) + 0.14 + 0.30 * f + 0.16 * (j % 2) + rng.uniform(0.0, 0.04)
            pos[i] = pos[hb] + L * np.array([np.cos(th), np.sin(th)])
    lone = [i for i in range(n) if k[i] == 0]
    if lone:                                                         # elements without links, loosely outside
        R = max(np.hypot(*pos[i]) + rad(p[i]) for i in pos) + 0.14
        for j, i in enumerate(lone):
            th = np.deg2rad(200 + 140 * (j + 0.5) / len(lone)) + rng.uniform(-0.03, 0.03)
            pos[i] = (R + 0.10 * (j % 2)) * np.array([np.cos(th), np.sin(th)])
    P = np.array(list(pos.values())); P -= (P.max(0) + P.min(0)) / 2; sc = max(1.0, np.abs(P).max() / 1.12)
    pos = {i: P[j] for j, i in enumerate(pos)}
    pos = {i: v / sc for i, v in pos.items()}
    edges = [tuple(sorted((i, home[i]))) for i in leaves]
    edges += [tuple(sorted((hubs[a_], hubs[b_]))) for a_ in range(h) for b_ in range(a_ + 1, h) if A[hubs[a_], hubs[b_]] > 0]
    return pos, edges


def net_inset(ax, A, size=0.5, strong=0.25, alpha=0.92, seed=3):
    """A network in the reading and layout of Figure 2: node area by share p_i, node colour by log2(p_i N) on the scale
    of Figure 2, link width by the number of links between the pair. A ring is drawn as a closed chain (Figure 2A).
    Otherwise the hubs sit in the middle and every other element hangs on its hub, as in Figure 2D (hub_layout), and
    the links drawn are each element's tie to its hub, the ties among the hubs and the links of the strongest pairs,
    at least `strong` of the largest count. An element without links is a small open circle outside the network."""
    import networkx as nx
    A = np.asarray(A, float); n = len(A); k = A.sum(1); p = k / k.sum()
    ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.25, 1.25); ax.set_aspect("equal"); ax.axis("off")
    ax.patch.set_alpha(0.0)
    ring = np.all(k == k[0]) and np.all((A > 0).sum(1) == 2)
    if ring:
        G = nx.from_numpy_array(A); order = [0]
        while len(order) < n:
            order.append(next(v for v in G.neighbors(order[-1]) if v not in order))
        xy = loop_layout(n); pos = {v: xy[i] for i, v in enumerate(order)}
        edges = list(G.edges())
    else:
        mmax = A.max()
        pos, edges = hub_layout(A, k, p, size, seed)
        edges = sorted(set(edges) | {(int(u), int(v)) for u, v in np.argwhere(np.triu(A, 1) >= strong * mmax)})
    m = np.array([A[u, v] for u, v in edges]); mmax = m.max()
    for (u, v), mm in sorted(zip(edges, m), key=lambda z: z[1]):
        ax.plot(*zip(pos[u], pos[v]), color=LINK, alpha=alpha, lw=(0.35 + 1.7 * np.sqrt(mm / mmax)) * np.sqrt(size),
                solid_capstyle="round", zorder=2)
    for i in range(n):
        if k[i] == 0:
            ax.scatter(*pos[i], s=size * 6, color="white", edgecolors=GREY, linewidths=0.4, zorder=5, alpha=alpha)
        else:
            ax.scatter(*pos[i], s=size * (3 + 560 * p[i]), color=[share_colour(np.log2(p[i] * n))], edgecolors=EDGE,
                       linewidths=0.4, zorder=5, alpha=alpha)
    return pos


def inset_at(sheet, ax, xd, yd, w):
    """A square axes of side w inches centred on the data point (xd, yd) of ax."""
    fig = sheet.fig; fig.canvas.draw()
    X, Y = fig.dpi_scale_trans.inverted().transform(ax.transData.transform((xd, yd)))
    return sheet.ax(X - w / 2, Y - w / 2, w, w)


def net_above(sheet, ax, xd, yd, w, A, lines, mk, target, rad=0.0, turn=0.0):
    """A network inset of side w inches centred on the data point (xd, yd) of ax, with its label of two lines above
    it, the first line led by the marker of the same point in the panel, and a narrow dotted arrow from the network to
    the point `target` of ax that it comes from. The network is drawn slightly transparent where it meets the
    panel."""
    a_ = inset_at(sheet, ax, xd, yd, w); net_inset(a_, A, size=0.34)
    l1, l2 = lines.split("\n")
    t1 = a_.text(0.5, 1.30, l1, transform=a_.transAxes, fontsize=6.2, ha="center", va="center")
    a_.text(0.5, 1.10, l2, transform=a_.transAxes, fontsize=6.2, ha="center", va="center")
    sheet.fig.canvas.draw()
    bb = t1.get_window_extent(sheet.fig.canvas.get_renderer()).transformed(a_.transAxes.inverted())
    a_.plot(bb.x0 - 0.10, 1.30, mk, ms=3.4, color="white" if mk == "o" else INK, mec=INK,
            mew=0.8 if mk == "o" else 0.3, transform=a_.transAxes, clip_on=False)
    # the arrow starts at the edge of the network nearest the target and ends just short of the point
    X0, Y0 = ax.transData.transform((xd, yd)); X1, Y1 = ax.transData.transform(target)
    r = 0.40 * w * sheet.fig.dpi; ang = np.arctan2(Y1 - Y0, X1 - X0) + np.deg2rad(turn)
    start = ax.transData.inverted().transform((X0 + r * np.cos(ang), Y0 + r * np.sin(ang)))
    ax.annotate("", xy=target, xytext=start, zorder=6.5, annotation_clip=False,
                arrowprops=dict(arrowstyle="-|>,head_length=0.32,head_width=0.16", color=EDGE, lw=0.6, ls=(0, (1.2, 1.4)),
                                shrinkA=0, shrinkB=3.2, connectionstyle=f"arc3,rad={rad}"))
    return a_
