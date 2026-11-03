**🇺🇸 English** · [🇧🇷 Português](README.pt-BR.md)

# mmpbsa_plot

MM/PBSA plots for **g_mmpbsa** users (the version for the GROMOS force field).
Developed at the **LMM (Laboratório de Modelagem Molecular), UEFS**.

The script scans the subfolders, detects the files of each system and opens a menu to
choose what to plot. At the end it writes a `.txt` report with your choices and the
numerical results.

## Requirements

Python 3.8+ with `numpy`, `pandas` and `matplotlib`.

```bash
pip install numpy pandas matplotlib
```

## Expected layout

One subfolder per system (replicate, ligand, etc.), containing the `g_mmpbsa` output files:

```
my_project/
├── system_A/
│   ├── binding_energy.xvg            # energy vs. time (time + vdW, elec., polar, nonpolar, total)
│   ├── energy_summary.csv            # component averages ± SD
│   └── residues_energy_summary.csv   # per-residue decomposition
├── system_B/
│   └── ...
└── mmpbsa_plot.py
```

Not every folder needs all three files: the menu only offers what exists.
Energies are read in kJ/mol and shown in kcal/mol; time is shown in ns.

## Usage

```bash
python mmpbsa_plot.py                 # interactive menu, current folder
python mmpbsa_plot.py /path/to/data   # another root folder
```

The first menu asks for the **language** (Português do Brasil or English US). The choice
applies to everything: menus and messages, plot texts, table headers and the `.txt`
report. In Portuguese, axis numbers use a decimal comma. The last language used is saved
in `plot_style.ini`. To skip the question: `--lang pt` or `--lang en`.
(The `--help` text is always in English.)

Then the menu asks, in this order:

1. **Systems** to include (`1,3`, `1-3` or `all`).
2. **Label** of each system. The default is the folder name; Enter keeps it.
3. **Plots**: energy vs. time (total, vdW, electrostatic, polar, nonpolar),
   averages ± SD, per-residue decomposition.
4. **Time range** in ns (no default: you type start and end).
5. **Residues**: Top 5, Top 10 or All. *All is not recommended*: with hundreds of
   residues the bar chart is unreadable, so only the sequence-position profile and the
   complete CSV are generated.
6. **Ligand**: the script detects the ligand (the last residue that is not a standard
   amino acid) and asks whether to keep it in the ranking. If kept, it usually dominates
   the top.
7. **Colors**: comma-separated hex (`#0072B2,#D55E00`) or Enter for the default palette.
8. **Organization** of the output folders (see below).
9. **Formats** (png, pdf, svg), CSV export of the data, and the interactive window.

### Output

By default everything goes to `MM_PBSA_plots/`, organized by plot type:

```
MM_PBSA_plots/
├── energy_vs_time/     # energy over time (one plot per component)
├── summary/            # component averages ± SD
├── per_residue/        # per-residue decomposition (one per system)
├── tables/             # CSVs: rankings, "mean ± SD" tables, plot data
└── mmpbsa_report_YYYYMMDD_HHMMSS.txt
```

In Portuguese the folders are named `energia_vs_tempo/`, `medias/`, `por_residuo/` and
`tabelas/`. Other layouts are available: **by file format** (`png/`, `pdf/`, `svg/`,
handy when submitting) or **everything in one folder**.

### What the averages plot shows (`energy_components_mean`)

Bars with the mean of each binding-energy term, per system, read from
`energy_summary.csv`: vdW, electrostatic, polar solvation, nonpolar solvation and total
(total = sum of the four). The average is over the **whole trajectory** and the error bar
is the **standard deviation across frames** (the energy fluctuation), not the uncertainty
of the mean. State this in the figure caption. The entropic term (−TΔS) is not included.

## Editing labels, text and font

On the first run a `plot_style.ini` is created in the root folder. Edit it and run again:

```ini
[names]
# accepts matplotlib mathtext (comments only on their own line)
system_A = $\it{Mtb}$PS - KS1

[text_en]
energy_ylabel = ENERGY (kcal/mol)
residue_title = Per-residue energy decomposition
comp_vdw = vDW

[font]
family = Arial
size = 12
```

The texts live in two blocks, `[text_en]` and `[text_pt]`, one per language.
A `plot_style.ini` from earlier versions is upgraded automatically on the first run
(the labels you defined are kept).
The default font is **Arial**. If it is not installed (common on Linux), the script uses
Liberation Sans, which has the same metrics, and records this in the report.

For fine-tuning outside the script:

- **SVG and PDF** keep the text editable. Open them in Inkscape or Illustrator
  (for SVG the font must be installed on the machine that opens it).
- **Interactive window** (menu option): the *Figure options* button in the matplotlib
  toolbar changes title, labels, limits and colors before saving.
- **CSV** of the data: to redo the plot in Origin, Excel or Grace.

## Running without the menu

Any option skips the corresponding question. With `-y` nothing is asked: it uses the
defaults (all systems, all formats, ligand dropped), so `--plots` is required, and so are
`--time` and `--top` when those plots are selected.

```bash
python mmpbsa_plot.py . --systems 1,3 --plots total,vdw,residues \
    --time 20 100 --top 10 --ligand drop --formats png,svg --csv \
    --lang en --layout type -y
```

`python mmpbsa_plot.py --help` lists all options.

## Notes

- The averages plot uses `energy_summary.csv`, i.e. the **whole trajectory**. The averages
  in the report for the chosen time range are computed from the `.xvg` and may differ.
- The `*_top_N_table.csv` file has the values in `mean ± SD` format.

## Authors

- Késsia Souza Santos: LMM (Laboratório de Modelagem Molecular), UEFS
