# Clean-up for ADR-001 (12 Oct 2026, Windows PC `D:\CLG\SEM-7\FINAL PROJECT`)

**Nothing was trashed or deleted.**

## Why there was nothing to remove

This PC had no working copy of the repository. The old `Scratch/` folder with the F0–F5 run outputs isn't on this
machine. I cloned a fresh copy from GitHub (`main` @ `fa15fa2`) into this folder. The clone contains only
tracked files (about 3 MB), and all of them are kept:

| Item | Size | Made by | Needed by ADR-001 | Decision |
|---|---|---|---|---|
| `src/aqf/` | 21 KB | core package | loaders, climatology (baselines) | keep |
| `Data Analysis/` | 78 KB | phase 1 code and report | data-quality rules referenced | keep |
| `Forecasting/src`, `PLAN.md`, `REPORT.md` | ~200 KB | phase 2 | frozen baselines (damped persistence, climatology, Ridge) | keep |
| `Forecasting/results/F5_final/` | 1.2 MB | F5, scored once on 2024 | the 2024 test record must stay documented | keep |
| `tools/`, `tests/`, `configs/`, `docs/`, `uv.lock` | small | project | yes | keep |

Run outputs (`Data Analysis/results/`, `Forecasting/results/F0..F4`, `data/processed/`) are gitignored and
rebuildable, so none exist in a fresh clone.

## Ambiguous items: kept and listed (not part of this repository)

| Path | Size | What it is | Why it was kept |
|---|---|---|---|
| `D:\CLG\SEM-7\FINAL PROJECT\AirQ\Dataset\` | 1.8 GB | byte-identical copy of `..\Dataset\Raw` (only its README differs) | it sits inside another git repo (`antonioroger2/AirQ`, no commits yet, LFS configured), so it may be the copy meant for that repo. Safe to remove if that repo isn't needed. |
| `D:\CLG\SEM-7\FINAL PROJECT\Dataset\Raw\` | 1.8 GB | TNAQ export (Tamil Nadu, 2018–2024), the data ADR-001 Phase 1 uses | data: never touched |
| `D:\FinalYear\` | 2.0 GB | other repo (`Adxrsh-17/AirQ-CV`) | not this project |
| `C:\Users\Admin\.gemini\antigravity\scratch\AirQ-CV-Integration\` | 83 MB | another tool's scratch clone of `Adxrsh-17/AirQ-CV` | not this project |
| `D:\CLG\SEM-7\PROJECT\MIP\` | — | blood-cell MIP course project | unrelated |

`trash-put` isn't available on Windows. If anything above should go, send it to the Recycle Bin (recoverable).
Never use `rm`.
