"""Write one Jupyter notebook for each figure of the paper, from the script of that figure, and run it.

Each notebook (notebooks/<script>.ipynb) holds the code of its script with three changes: the folders of the repository
are found from the notebook, nothing is written to disk, and the finished figure is shown in the notebook. Comments
are left out. A last cell sets the printed values beside the values of the paper.

This builder runs every code cell in order in one namespace, as a kernel would, and stores what each cell prints and
the figure it shows.
Usage: python tools/make_notebooks.py <repository root> <deposit of Park (2026a)>
"""
import ast, base64, contextlib, io, json, os, re, sys, tokenize
from pathlib import Path

REPO, DEPOSIT = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
MARK_RE = re.compile(r"# -{20,} (.*)")
LETTER_RE = re.compile(r"# ([A-Z]\. .*)")
NEEDS_DEPOSIT = ("fig2_mechanism", "figS1_network_models")
DEPOSIT_FILE = "deposit"
MODS = {"style": "figures/style.py", "sim_full": "code/sim_full.py", "sim_practice": "code/sim_practice.py",
        "force_spread": "code/force_spread.py", "network_models": "code/network_models.py"}

FIGS = [  # script, title, what it is, how it works
    ("fig1_classical_laws", "Figure 1. Classical laws and equations of speed and accuracy",
     "The four classical records, each with its line, and the square-root rule of force.",
     "Reads the transcribed tables and the values read from the published figures, fits one least-squares line to each record and draws the five panels."),
    ("fig2_mechanism", "Figure 2. How coupling sets the slope of speed and accuracy",
     "The dynamic law as a scheme, from the novice with even coupling (A) to the expert with hub coupling (D).",
     "Builds the four structures of 20 elements, runs the full model of practice 300 times for H(t), computes the law at A, at D and after practice, and draws the race under A and D."),
    ("fig3_practice_step", "Figure 3. One practice step of the simulations",
     "The practice step as a flow chart.",
     "Draws the boxes and arrows. Nothing is computed and no data are read."),
    ("fig4_minimal_model", "Figure 4. Minimal simulation of practice on a fixed set of 20 elements",
     "Averaged information under three rules of attachment, where the added coupling goes, and a late element that becomes a hub.",
     "Runs the minimal model 300 times for each rule and 200 times for the late element, and reads the stored shares and networks."),
    ("fig5_full_model", "Figure 5. Full model of practice with selective elimination and recruitment",
     "Averaged information with elimination and recruitment, where the coupling goes, and the closed forms against the simulated change.",
     "Reads the stored output of the full model and of the runs that follow the shares, and draws the three panels."),
    ("fig6_hick", "Figure 6. Hick's law as a static case of the dynamic law",
     "Participant G. C. of Hyman (1953), Experiment I, with the law at steps 0 and 1,000, the halvings of a choice, and the pairs of β and H that give one slope.",
     "Reads the eight reaction times, fits the line, runs the full model 300 times for H at step 1,000 and reads the simulated performer."),
    ("fig7_hyman", "Figure 7. The Hick–Hyman law as a static case of the dynamic law",
     "The 24 conditions of G. C. on one line, the two averaged informations of the law, and each light against the mean of its condition.",
     "Reads the 24 reaction times and the probabilities of the lights, fits one line and one line for each experiment, and reads the simulated performer."),
    ("fig8_fitts", "Figure 8. Fitts' law and Woodworth's two components as a static case of the dynamic law",
     "Fitts (1954) Table 1 with the law at steps 0 and 1,000, one aimed movement as initial adjustment and corrections, and two ways to a lower slope.",
     "Reads the table, fits the line of the one-ounce stylus, runs the full model 300 times for H at step 1,000 and reads the simulated performer."),
    ("fig9_impulse_variability", "Figure 9. The impulse-variability equation as a static case of the dynamic law",
     "The nine effective widths of Schmidt et al. (1979) with their line, the spread of force over the shares, and the correlation that keeps the slope.",
     "Reads the nine values and the shares of the simulated performer, fits the line and evaluates Equation 22."),
    ("fig10_slope_reading", "Figure 10. Leonard's practice read on the path of averaged information",
     "Leonard's three months placed on the path of H(t) of the full model by the ratio of their slopes.",
     "Reads the practice table, runs the full model 300 times, finds the steps at which H equals each month's value and draws the coupling at those steps."),
    ("fig11_system_model", "Figure 11. The dynamic law as a systems model of practice",
     "The loops of practice as a diagram, and the law at successive steps of practice.",
     "Draws the diagram, runs the full model for the path of H(t) and reads the stored network of one run."),
    ("figS1_network_models", "Figure S1. Network models and how far averaged information lies below its ceiling",
     "Three network rules, their degree distributions, and the gap between H and its ceiling.",
     "Builds the networks of Table S3 from their seeds, pools the degrees of 20 runs and evaluates Equations S2 and S3."),
    ("figS2_static_cases", "Figure S2. The classical laws in one simulated performer at two steps of practice",
     "Hick's, Hyman's and Fitts' tasks in one simulated performer at steps 0 and 1,000, and the ratio of the slopes.",
     "Reads the stored output of the simulated performer and fits a line to each task at each step."),
    ("figS3_practice_tables", "Figure S3. The two classical practice tables of the slope test",
     "The slope at each stage of practice in Leonard's table and in that of Mowbray and Rhoades, and the time saved against the bits of the choice.",
     "Reads the two transcribed tables and takes the differences between the set sizes, which differ by one bit."),
    ("figS4_open_records", "Figure S4. Skill and practice in two open records",
     "The slope of the slower and the faster typists of How We Type, and the mean score over the first 30 plays of the online game.",
     "Reads the outputs of the two analyses as run and draws them."),
    ("figS5_selection_race", "Figure S5. The selection of a command as a race among the elements",
     "The time to select an element against its information, the mean time of a command against H, and the wrong selections.",
     "Reads the stored output of the race simulation and fits the lines."),
]


PAPER = {
    "fig1_classical_laws": ["`MT = 12.8 + 94.7 ID`: the line of Fitts' Table 1 in the main text and in Figure 8A."],
    "fig2_mechanism": ["`scale-free H 3.9251`: 3.93 bits, in the note of Figure 2.",
                       "`practice: mean H 4.322 -> ... -> 3.005 (1000)`: from 4.32 to 3.01 bits.",
                       "`slopes A 94.70, D 86.00, t = 1000 65.85`: 94.7 and 65.9 ms per bit in the main text, 86.0 in Supplemental Section S2.2.",
                       "`T at log2 X = 4: 391.6`: the 391.6 ms of the 4-bit aiming task at step 0."],
    "fig3_practice_step": ["No value is printed. Compare the chart with Figure 3."],
    "fig4_minimal_model": ["`H(2000)` of 4.319, 4.079 and 3.717: 4.319, 4.08 and 3.72 bits.",
                           "`T(2000)` of 370.3 and 338.6 from `T(0)` 391.6: the times of the 4-bit aiming task.",
                           "`start (2, 2, 2)`, first share 0.15 to 0.466: the fit share from 0.15 to 0.47. Over the seven starts it ends between 0.277 and 0.685, the 0.28 to 0.69 of the text.",
                           "`87 percent of 200 runs`: the element recruited at step 50."],
    "fig5_full_model": ["`fixed H: 4.32 -> 3.01 (1000) -> 1.65 (3000)`, `no elimination 3.81, 3.66`, `recruitment 3.32 -> 4.05 (250) -> 3.70 -> 2.49`: the same values in the main text.",
                        "`windows: r = 0.994, slope 0.98, same sign 60 of 60`: Figure 5C.",
                        "`start (2, 2, 2)`, first share 0.688 at step 1,000 and 0.989 at step 3,000: 0.69 and 0.99 of all link ends. Over the seven starts it ends between 0.952 and 0.996, the 0.95 to 1.00 of the text."],
    "fig6_hick": ["`RT = 226.6 + 147.8 log2 n, r = 0.9926`: RT = 227 + 148 log₂n, r = .993.",
                  "`beta = 73.9 (N = 4), 34.20 (N = 20)`: 74 and 34 ms per bit.",
                  "`slope after practice 102.8`: 103 ms per bit."],
    "fig7_hyman": ["`all 24: RT = 217.3 + 147.8 H, r = 0.9851`: RT = 217 + 148 H_S, r = .985.",
                   "`slopes I 147.8, II 143.7, III 146.8`: 148, 144 and 147 ms per bit.",
                   "`two lights: (1.0, 339.6), (0.47, 306.5)`: 340 ms at 1 bit and 307 ms at 0.47 bits."],
    "fig8_fitts": ["`MT = 12.8 + 94.7 ID`: the line of the one-ounce stylus.",
                   "`residual by ID`, 72.5 at 1 bit: the easiest condition lies 73 ms above the line.",
                   "`slope after practice 65.85`: 65.9 ms per bit."],
    "fig9_impulse_variability": ["`W_e = 2.249 + 0.0303 D/MT`: W_e = 2.25 + 0.030 D/MT.",
                                 "`rbar keeping the slope -0.1644`: r̄ = −0.16."],
    "fig10_slope_reading": ["`Leonard slopes [62.0, 51.0, 45.0]`, `ratios [1.0, 0.8226, 0.7258]`: 62, 51 and 45 ms per bit, 82 and 73 percent.",
                            "`H targets [4.322, 3.555, 3.137]`: 3.56 and 3.14 bits.",
                            "`matched steps of the mean path [0, 593, 906]`: steps 593 and 906.",
                            "`time saved, three and six choices: (12.0, 29.0)`: 12 and 29 ms."],
    "fig11_system_model": ["`H(0) 4.322, H(1000) 3.005, H(3000) 1.645`: 4.32, 3.01 and 1.65 bits, as in Figure 5A."],
    "figS1_network_models": ["`continuum gap`, 3.4489, 1.3004 and 0.4427: 3.45, 1.30 and 0.443 bits (Section S2.2).",
                             "`exact growth gap` 0.7364, 0.6046, 0.5544 and `random, mean degree 2: 0.411`: 0.736, 0.605, 0.554 and 0.411 bits.",
                             "`preferential growth`, 0.7161 and 0.7304, and `growth with fitness`, 0.8702: 0.716, 0.730 and 0.870 bits.",
                             "`C = 2.2293, tail exponent 2.486, gap 0.5849`: C = 2.229, 2.49 and 0.585 bits.",
                             "`Figure 2 ensembles`, 0.0387, 0.3367 and 0.3966: 0.039, 0.337 and 0.397 bits."],
    "figS2_static_cases": ["`Hick 0.6971, Hyman 0.6971, Fitts 0.6975`, `H ratio 0.6954`, `L ratio 0.6972`: the ratio 0.697 against 0.695 in the main text."],
    "figS3_practice_tables": ["`Leonard slopes [62.0, 51.0, 45.0]` and `block slopes` 18.67 and 2.0: 62, 51, 45, 19 and 2 ms per bit (Table 3).",
                              "`time saved`, Leonard (12.0, 29.0) and Mowbray (75.0, 81.0): 12 and 29 ms, 75 and 81 ms."],
    "figS4_open_records": ["`slopes`, faster (35.25, 23.24, 47.9) and slower (61.78, 41.34, 70.8), `ratio 0.5706`: 35 [23, 48], 62 [41, 71] and 0.57 (Table 3).",
                           "`game 5293 0.6291 ('0.275', ...)`: 5,293 players, 62.9 percent and α = 0.275."],
    "figS5_selection_race": ["`A`, slope 0.6983, r = 0.9838: 0.698 time units per bit, r = .984 over 106 elements.",
                             "`B: theta = 3`, slope 0.6849, 190 and 132 ms, and `theta = 4`, slope 0.7026: 0.685 and 0.703 against 0.693, r = .9999 and .9995.",
                             "`B: even coupling`, slope 0.6861, r = 0.9997: 0.686 time units per bit.",
                             "`C`: wrong selections between 1.4 and 2.3 percent, and between 0.6 and 0.9 percent."],
}


def strip_comments(seg):
    """Remove comments from a statement, leaving strings as they are."""
    rows = seg.split("\n")
    cuts = {}
    try:
        for tok in tokenize.generate_tokens(io.StringIO(seg + "\n").readline):
            if tok.type == tokenize.COMMENT:
                cuts[tok.start[0] - 1] = tok.start[1]
    except tokenize.TokenError:
        return seg
    out = []
    for i, row in enumerate(rows):
        if i in cuts:
            row = row[:cuts[i]].rstrip()
            if not row.strip():
                continue
        out.append(row)
    return "\n".join(out)


def is_call(node, text):
    return isinstance(node, ast.Expr) and isinstance(node.value, ast.Call) and text in ast.unparse(node.value.func)


def transform(stem):
    """The code of the script as notebook cells, one for each section of the script."""
    src = (REPO / "figures" / f"{stem}.py").read_text(encoding="utf-8")
    lines = src.split("\n")
    body = ast.parse(src).body[1:]
    marks = {i + 1: MARK_RE.match(l).group(1) for i, l in enumerate(lines) if MARK_RE.match(l)}
    if not marks:
        marks = {i + 1: LETTER_RE.match(l).group(1) for i, l in enumerate(lines) if LETTER_RE.match(l)}
    cut_from = cut_to = None
    if stem == "fig2_mechanism":                                     # the saving of the figure and of its panels
        cut_from = next(i + 1 for i, l in enumerate(lines) if l.startswith("OUT = "))
        cut_to = next(n for n, t in marks.items() if t.startswith("the numbers the figure shows"))
    sections, last_end = {}, None
    for node in body:
        seg = ast.get_source_segment(src, node)
        if cut_from and cut_from <= node.lineno < cut_to:
            continue
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name in ("ROOT", "DEP") or (name in ("out", "OUT") and '"figures"' in seg):
                continue
        if any(is_call(node, t) for t in ("sys.path.insert", "mpl.use", "savefig", "save_panels", "mkdir")):
            continue
        if is_call(node, "print") and node.value.args and isinstance(node.value.args[0], ast.Constant) \
                and node.value.args[0].value == "written":
            rest = [a for a in node.value.args[1:] if not (isinstance(a, ast.Name) and a.id in ("out", "OUT"))]
            if not rest:
                continue
            seg = "print(" + ", ".join(ast.unparse(a) for a in rest) + ")"
        if isinstance(node, ast.ImportFrom) and node.module == "style":
            node.names = [a for a in node.names if a.name != "save_panels"]
            seg = ast.unparse(node)
        seg = strip_comments(seg)
        for a in (r'Path\(sys\.argv\[1\]\) if len\(sys\.argv\) > 1 else ROOT / "sim_out"',
                  r'Path\(sys\.argv\[1\] if len\(sys\.argv\) > 1 else "sim_out"\)'):
            seg = re.sub(a, 'ROOT / "sim_out"', seg)
        assert not re.search(r"__file__|sys\.argv|\.savefig\(|save_panels", seg), (stem, seg[:200])
        key = max([n for n in marks if n <= node.lineno], default=0)
        gap = last_end is not None and node.lineno - last_end > 1
        sections.setdefault(key, []).append(("\n" if gap and sections.get(key) else "") + seg)
        last_end = node.end_lineno
    return [(marks.get(k, ""), "\n".join(v)) for k, v in sorted(sections.items())]


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.strip("\n").splitlines(True)}


def code(text):
    return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
            "source": text.strip("\n").splitlines(True)}


SETUP = ('%matplotlib inline\n'
         '%config InlineBackend.figure_format = "retina"\n'
         '%config InlineBackend.print_figure_kwargs = {"bbox_inches": None}\n'
         'import logging, os, sys\n'
         'from pathlib import Path\n'
         'import matplotlib.pyplot as plt\n'
         'from IPython.display import display\n'
         'logging.getLogger("matplotlib").setLevel(logging.ERROR)\n'
         'plt.ioff()\n'
         'ROOT = Path.cwd().resolve()\n'
         'ROOT = ROOT.parent if ROOT.name == "notebooks" else ROOT\n'
         'os.chdir(ROOT)\n'
         'sys.path[:0] = [str(ROOT / "figures"), str(ROOT / "code")]\n')
SETUP_DEPOSIT = ('DEP = Path(os.environ.get("PARK2026_DEPOSIT", ROOT.parent / "control-entropy-coordinative-structures")).resolve()\n'
                 'assert (DEP / "code").is_dir(), f"no deposit of Park (2026a) at {DEP}"\n'
                 'sys.path.insert(0, str(DEP / "code"))\n')


def cells_of(stem, title, what, how, files):
    figure = "fig" if stem == "fig2_mechanism" else "S.fig"
    intro = f"# {title}\n\n**What it is.** {what}\n\n**How it works.** {how}\n\n"
    names = [f"`{f}`" if f != DEPOSIT_FILE else "`code/generators.py` of the deposit of Park (2026a)" for f in files]
    intro += "**Files read.** " + (", ".join(names) + "." if names else "None.")
    intro += "\n\n**To run.** Open this notebook from the `notebooks` folder of the repository, or from its root, and run all cells."
    if stem in NEEDS_DEPOSIT:
        intro += (" The structures of Park (2026a) come from its deposit: clone "
                  "github.com/pcw8531/control-entropy-coordinative-structures beside this repository and check out "
                  "commit 0660a649f5e0, or set the environment variable `PARK2026_DEPOSIT` to its path.")
    cells = [md(intro), code(SETUP + (SETUP_DEPOSIT if stem in NEEDS_DEPOSIT else ""))]
    for name, text in transform(stem):
        if not text.strip():
            continue
        if name:
            name = name.strip()
            cells.append(md("## " + name[0].upper() + name[1:]))
        cells.append(code(text))
    cells.append(md("## The figure"))
    cells.append(code(("" if stem == "fig2_mechanism" else "from style import bolden\nbolden(S.fig)\n")
                      + f"display({figure})\nplt.close(\"all\")"))
    cells.append(md("## In the paper\n\nPrinted here, then as the manuscript or the Supplemental Material gives it.\n\n"
                    + "\n".join("- " + line for line in PAPER[stem])))
    for i, c in enumerate(cells):
        c["id"] = f"cell-{i + 1:02d}"
    return cells


def run(stem, cells):
    """Run the code cells in one namespace from the notebooks folder. Returns the files of the repository they read."""
    import matplotlib
    matplotlib.rcdefaults()
    for m in list(MODS) + ["generators"]:
        sys.modules.pop(m, None)
    os.chdir(REPO / "notebooks")
    saved_path, shown, opened = list(sys.path), [], set()

    def display(fig):                                                # stands in for the display of the kernel
        buf = io.BytesIO(); fig.canvas.print_figure(buf, format="png", dpi=fig.dpi * 2, bbox_inches=None)
        w, h = fig.get_size_inches() * fig.dpi
        shown.append((buf.getvalue(), int(round(w)), int(round(h))))

    def hook(event, args):
        if event == "open" and isinstance(args[0], (str, os.PathLike)) and (args[1] is None or "r" in str(args[1])):
            p = Path(os.fspath(args[0])); p = (p if p.is_absolute() else Path.cwd() / p).resolve()
            if REPO in p.parents and p.suffix in (".csv", ".txt", ".npz"):
                opened.add(p.relative_to(REPO).as_posix())
    WATCH["hook"] = hook
    ns, n = {"__name__": "__main__", "display": display}, 0
    try:
        for c in cells:
            if c["cell_type"] != "code":
                continue
            n += 1; c["execution_count"] = n
            text = "\n".join(l for l in "".join(c["source"]).split("\n")
                             if not l.startswith("%") and l != "from IPython.display import display")
            buf = io.StringIO(); shown.clear()
            with contextlib.redirect_stdout(buf):
                exec(compile(text, f"<{stem} cell {n}>", "exec"), ns)
            c["outputs"] = []
            if buf.getvalue():
                c["outputs"].append({"output_type": "stream", "name": "stdout", "text": buf.getvalue().splitlines(True)})
            for png, w, h in shown:
                c["outputs"].append({"output_type": "display_data",
                                     "data": {"image/png": base64.b64encode(png).decode("ascii"), "text/plain": ["<Figure>"]},
                                     "metadata": {"image/png": {"width": w, "height": h}}})
        files = sorted(opened) + [f for m, f in MODS.items() if m in sys.modules]
        if "generators" in sys.modules:
            files.append(DEPOSIT_FILE)
    finally:
        WATCH["hook"] = None
        sys.path[:] = saved_path
        import matplotlib.pyplot as plt
        plt.close("all")
    return files


WATCH = {"hook": None}
sys.addaudithook(lambda e, a: WATCH["hook"](e, a) if WATCH["hook"] else None)
sys.dont_write_bytecode = True
os.environ["PARK2026_DEPOSIT"] = str(DEPOSIT); os.environ["MPLBACKEND"] = "Agg"
(REPO / "notebooks").mkdir(exist_ok=True)
for stem, title, what, how in FIGS:
    files = run(stem, cells_of(stem, title, what, how, []))          # first run: the files the notebook reads
    cells = cells_of(stem, title, what, how, files)
    assert run(stem, cells) == files, stem
    nb = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                       "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}
    (REPO / "notebooks" / f"{stem}.ipynb").write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"notebooks/{stem}.ipynb: {sum(c['cell_type'] == 'code' for c in cells)} code cells, {len(files)} files read")
