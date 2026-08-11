# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

MTEX example scripts — a standalone git repo that lives inside the MTEX source
tree at `mtex/examples`. Everything here is MATLAB `.m` code that runs inside
MATLAB with MTEX on the path. There is no build system, package manager, linter,
or test suite; correctness is checked by running a script and looking at the
figures it produces.

Every subdirectory (`ExGrains`, `ExSeismics`, `ExPlasticity`, `trueEbsd`, …)
holds small self-contained demos. There is no library here — no classes, and
almost no shared code. The two exceptions are local to their own directory:
`ExPlasticity/DisTypesFo.m` is a real function (`dSo = DisTypesFo(CS)`), and
`JAC-Creuziger/ColorMap4|8|20.m` are colormap scripts pulled in by the other
JAC scripts with `run('ColorMap4.m')`.

## Running things

Everything runs interactively in MATLAB, never from the shell:

```matlab
run('ExGrains/ExIceSphericity.m')     % a plain example script

addpath(genpath('<path to mtex-trueEbsd>'))   % external toolbox, see below
run('trueEbsd/example_WCCo.m')                % downloads its data via mtexdata
```

- Data comes either from `mtexdata <name>` (downloaded by MTEX, e.g.
  `mtexdata trueEbsdWCCo`) or from a local file next to the script.
- **`mtexdata` assigns into the *base* workspace**, so a script calling it
  cannot be driven from inside a wrapper function — `run` it from base, or the
  variable never appears.
- The published examples never depend on the current directory: they build
  absolute paths with `mtexExamplePath`, e.g.
  `fullfile(mtexExamplePath,'ExODFReconstruction','data','alt4_*.rw1')` or
  `path = [mtexExamplePath filesep 'ExGrains' filesep]`. Keep new examples to
  that convention. `JAC-Creuziger/` is the exception — see below.
- `trueEbsd/example_WCCo.m` needs the TrueEBSD toolbox, which is **not** part of
  MTEX or of this repo: <https://github.com/vtvivian/mtex-trueebsd> (branch
  `mtex7-compat` for MTEX 7), Apache-2.0. It also needs MATLAB R2024a+ and the
  Image Processing, Curve Fitting, and Statistics and Machine Learning
  toolboxes, and takes ~11 minutes to run.

## Documentation build

These scripts are published to HTML by the separate `mtex/makeDoc` toolbox (see
`../makeDoc/CLAUDE.md`), which runs MATLAB's `publish` over the `%%`-section
comments. That means the prose comments are user-facing documentation, not
incidental notes — `%%` starts a titled section, `%%%` a subsection, and
`|foo|` renders as inline code.

The table of contents is a two-level plain-text structure:

- `Examples.toc` lists directories (`<DirName> <Display Title>`).
- `<Dir>/<Dir>.toc` lists the scripts inside it (`<ScriptName> <Display Title>`).
- `<Dir>/<Dir>.m` is a stub holding the section title for that directory
  (`Examples.m` is the same stub for the top level).

Adding an example means adding the `.m` file *and* its `.toc` line. The
converse also holds: **a `.m` file in a listed directory is not published
unless its `.toc` line exists**, which is why `ExPlasticity/` contains
`SynEBSD.m` and `DisTypesFo.m` that no built page links to. Don't "fix" that
by adding them to the `.toc` — publishing them is a deliberate decision.

## Directories outside the documentation

Three directories have no `.toc` and are absent from `Examples.toc`, so nothing
in them reaches the built documentation. They are still part of the repo and
still expected to run.

- **`trueEbsd/`** — one demo needing an external toolbox; see below.
- **`phaseTransformation/`** — a single script,
  `Forsterite_to_WadsleyiteMainprice.m`.
- **`JAC-Creuziger/`** — a vendored NIST contribution (pole-figure inversion
  for TRIP steels), copied from
  <https://github.com/usnistgov/Texture-Sampling-PhaseMeasurement-BiasErrors>.
  It carries **its own `README.md` and `LICENSE.TXT`** — the NIST license, not
  the repo's root `LICENSE` — and its README is the authority on the data,
  the `.maa`/`.sum` file formats, and citation requirements. Treat it as
  third-party: fix things here only when asked, and preserve the attribution
  and disclaimer text.

  Its scripts are the one place where **the current directory matters**: they
  use relative paths (`run('ColorMap4.m')`, `pname = './ExperimentalData/…'`,
  `savepath = 'ExpFigures'`), so MATLAB must be `cd`-ed into `JAC-Creuziger/`.
  They write into `ExpFigures/`, `MtexData/`, `MtexDataHW/`, `ODFFigures/`,
  `ODFFiguresHW/`, all gitignored — generated output, never commit it. They
  also call `export_fig`, a third-party function that is not part of MTEX.

## trueEbsd

`trueEbsd/example_WCCo.m` is a single demo script, like everything else here.

It used to be different: this directory carried a **vendored copy of the whole
TrueEBSD toolbox** — `@trueEbsd` and `@distortedImg`, the `funcsV0` DIC code
with its committed Windows mex, and a `tools/` folder of coordinate helpers —
plus `example_copper.m` and a committed `html/`. That copy was a fork of an old
upstream snapshot: it had drifted to a **handle**-class `@trueEbsd` (methods
mutating in place, `job.calcShifts('fitErr')`) and cell-array indexing
(`job.resizedList{n}`), neither of which matches released TrueEBSD.

All of it was deleted. TrueEBSD is now consumed as the external toolbox it is:

- Upstream: <https://github.com/vtvivian/mtex-trueebsd>, Apache-2.0, v2.1.0.
  Use branch `mtex7-compat` with MTEX 7.
- `@trueEbsd` there is a **value** class — reassign, `job = calcShifts(job,…)`.
- The lists are **object arrays**, `job.resizedList(n)`, not cell arrays. The
  exception is `job.shifts`, a cell of `@pairShifts` arrays indexed
  `job.shifts{n}(m)` for hop `n`, distortion-model stage `m`.
- `ij2EbsdSquare` and friends live in that toolbox's `tools/`, which is why the
  script needs it on the path, not just for the classes.

Do not re-vendor it here. If the script breaks, the fix belongs upstream or in
the script, not in a local copy of the classes.

## Conventions and current state

- Comment style is heavy and explanatory by design. The prose comments *are*
  the published documentation, so write them for a reader, not as notes.
- `.asv` files are MATLAB autosaves — ignore them, never edit.
- `ExODFReconstruction/Knie_10/` is a data directory with thousands of files;
  avoid recursive searches through it.
- Input data lives in the repo next to the script that reads it (`PIL185.ctf`,
  `Al2O3-Corundum.cif`, `data/alt4_*.rw1`, …) and is committed. Anything a run
  *produces* is gitignored — `trueEbsd/*_out.mat`, the JAC output folders — so
  a dirty working tree full of those is normal, not something to clean up or
  commit.
