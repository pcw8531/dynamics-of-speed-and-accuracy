"""Write one Jupyter notebook for each figure of the paper, from the script of that figure, and run it.

Each notebook (notebooks/<script>.ipynb) has the description of the figure, a setup cell that finds the repository
root and gives the script the argument it takes on the command line, the code of the script split into cells at its
section marks, and a last cell that lists the finished figure and every panel it writes to figures/panels. The code
is the code of the script, with only three changes so that it runs in a notebook: the folder of the script becomes
ROOT / "figures", its parent becomes ROOT, and sys.argv[1] becomes ARGS[0].

This builder runs every code cell in order in one namespace, as a kernel would, and keeps what each cell prints, so
that the notebooks open with the values they report. No plot is stored in a notebook. The images are in figures/ and
figures/panels, and running a notebook again in Jupyter writes the same files.
Usage: python tools/make_notebooks.py <repository root> <deposit of Park (2026)>
"""
import contextlib, io, json, os, re, sys
from pathlib import Path

REPO, DEPOSIT = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
FIGS = [  # script, title, panels written to figures/panels, whether it needs the deposit of Park (2026)
    ("fig1_classical_laws", "Figure 1. Classical laws and equations of speed and accuracy, each with its task", "ABCDE", False),
    ("fig2_mechanism", "Figure 2. How coupling sets the slope of speed and accuracy", ["A", "B", "C", "D", "centre"], True),
    ("fig3_practice_step", "Figure 3. One practice step of the simulations", "", False),
    ("fig4_minimal_model", "Figure 4. Minimal simulation of practice on a fixed set of 20 elements", "ABC", False),
    ("fig5_full_model", "Figure 5. Full model of practice with selective elimination and recruitment", "ABC", False),
    ("fig6_hick", "Figure 6. Hick’s law as a static case of the dynamic law", "ABC", False),
    ("fig7_hyman", "Figure 7. The Hick–Hyman law as a static case of the dynamic law", "ABC", False),
    ("fig8_fitts", "Figure 8. Fitts’ law and Woodworth’s two components as a static case", "ABC", False),
    ("fig9_impulse_variability", "Figure 9. The impulse-variability equation as a static case of the dynamic law", "ABC", False),
    ("fig10_slope_reading", "Figure 10. Leonard’s practice read on the path of averaged information", "ABC", False),
    ("fig11_system_model", "Figure 11. The dynamic law as a systems model of practice", "AB", False),
    ("figS1_network_models", "Figure S1. Network models and how far averaged information lies below its ceiling", "ABCD", True),
    ("figS2_static_cases", "Figure S2. The classical laws in one simulated performer at two steps of practice", "ABCD", False),
    ("figS3_practice_tables", "Figure S3. The two classical practice tables of the slope test", "ABC", False),
    ("figS4_open_records", "Figure S4. Skill and practice in two open records", "AB", False),
]
MARK = "# ---------------------------------------------------------------- "


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.strip("\n").splitlines(True)}


def code(text):
    return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
            "source": text.strip("\n").splitlines(True)}


def notebook_cells(stem, title, panels, deposit):
    src = (REPO / "figures" / f"{stem}.py").read_text()
    doc, body = src.split('"""', 2)[1:]
    body = (body.replace("Path(__file__).resolve().parent.parent", "ROOT")
                .replace("Path(__file__).resolve().parent", '(ROOT / "figures")')
                .replace("len(sys.argv) > 1", "len(ARGS) > 0").replace("sys.argv[1]", "ARGS[0]")
                .replace('str((ROOT / "figures"))', 'str(ROOT / "figures")').replace("ROOT = ROOT\n", ""))
    assert "__file__" not in body and "sys.argv" not in body, stem
    run_line = [l for l in doc.splitlines() if l.startswith("Run from the repository root")]
    words = " ".join(l.strip() for l in doc.splitlines() if not l.startswith("Run from"))
    intro = (f"# {title}\n\n{words}\n\n"
             "**To run.** Open this notebook from the `notebooks` folder of the repository, or from its root, and run "
             "all cells. It needs Python 3 with the packages of `requirements.txt`, and the outputs of the "
             "simulations in `sim_out`. The same code is in `figures/" + stem + ".py`, and the two write the same "
             "files.")
    if run_line:
        intro += f" From a terminal: `{run_line[0].split(': ', 1)[1].strip()}`."
    setup = ('import os, sys\nfrom pathlib import Path\nimport matplotlib.pyplot as plt\n'
             'plt.ioff()                                                           # the figure goes to its file, and no plot is shown\n'
             'ROOT = Path.cwd().resolve()\nROOT = ROOT.parent if ROOT.name == "notebooks" else ROOT\n'
             'os.chdir(ROOT)                                                       # the script runs from the repository root\n')
    if deposit:
        intro += (" The structures of Park (2026) come from its deposit: clone "
                  "github.com/pcw8531/control-entropy-coordinative-structures beside this repository and check out "
                  "commit 0660a649f5e0, or set the environment variable `PARK2026_DEPOSIT` to its path.")
        setup += ('DEP = Path(os.environ.get("PARK2026_DEPOSIT", ROOT.parent / "control-entropy-coordinative-structures")).resolve()\n'
                  'assert (DEP / "code").is_dir(), f"no deposit of Park (2026) at {DEP}"\n'
                  'ARGS = [str(DEP)]                                                    # the argument of the script\n')
    else:
        setup += 'ARGS = [str(ROOT / "sim_out")]                                       # the argument of the script\n'
    cells = [md(intro), code(setup)]
    if MARK not in body:                                             # panels marked "# A. ..." instead
        body = re.sub(r"\n# ([A-Z])\. ", lambda m: "\n" + MARK + m.group(1) + ". ", body)
    parts = body.split(MARK)
    head = parts[0].strip("\n")
    if head:
        cells.append(code(head))
    for part in parts[1:]:
        name, rest = (part.split("\n", 1) + [""])[:2]
        cells.append(md("## " + name.strip()[0].upper() + name.strip()[1:]))
        cells.append(code(("# " if re.match(r"[A-Z]\. ", name.strip()) else MARK) + name.strip() + "\n" + rest))
    files = [f"figures/{stem}.png"] + [f"figures/panels/{stem.split('_')[0]}_{p}.png" for p in panels]
    cells.append(md("## Files written\n\nThe finished figure as it goes into the article"
                    + (", and each panel as it is placed on its own" if panels else ", placed whole") + ". Open them in "
                    "`figures` and `figures/panels`."))
    cells.append(code("for f in " + json.dumps(files) + ":\n"
                      "    print(f, \"written\" if (ROOT / f).is_file() else \"missing\")"))
    for i, c in enumerate(cells):
        c["id"] = f"cell-{i + 1:02d}"                                 # nbformat 4.5 gives every cell an id
    return cells, files


os.environ["PARK2026_DEPOSIT"] = str(DEPOSIT); os.environ["MPLBACKEND"] = "Agg"
(REPO / "notebooks").mkdir(exist_ok=True)
for stem, title, panels, deposit in FIGS:
    cells, files = notebook_cells(stem, title, panels, deposit)
    os.chdir(REPO / "notebooks")
    import matplotlib
    matplotlib.rcdefaults()                                          # a fresh kernel for each notebook
    saved_path = list(sys.path)
    ns, n = {"__name__": "__main__"}, 0
    for c in cells:
        if c["cell_type"] != "code":
            continue
        n += 1; c["execution_count"] = n
        src = "".join(c["source"])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            exec(compile(src, f"<{stem} cell {n}>", "exec"), ns)
        text = buf.getvalue().replace(str(REPO) + os.sep, "")             # paths as seen from the repository root
        if text:
            c["outputs"].append({"output_type": "stream", "name": "stdout", "text": text.splitlines(True)})
    import matplotlib.pyplot as plt
    plt.close("all"); sys.path[:] = saved_path
    for m in [m for m in sys.modules if m in ("style", "generators", "sim_full", "sim_practice", "sim_composition")]:
        del sys.modules[m]                                            # each notebook starts with a fresh kernel
    nb = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                       "language_info": {"name": "python", "pygments_lexer": "ipython3"}},
          "nbformat": 4, "nbformat_minor": 5}
    (REPO / "notebooks" / f"{stem}.ipynb").write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
    assert "image/png" not in json.dumps(nb) and all("id" in c for c in cells), stem
    print(f"notebooks/{stem}.ipynb: {n} code cells, {len(files)} files written")
