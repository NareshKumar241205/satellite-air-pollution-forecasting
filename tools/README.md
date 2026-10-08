# tools/

## `view_tiff.ipynb`: view any GeoTIFF

1. Open the notebook and choose the kernel **Python (Scratch .venv)**.
2. Set `IMAGE_PATH` in the settings cell.
3. Run All.

Settings: `BANDS`, `RGB` (Sentinel-2 true colour), `SAME_SCALE` (one colour scale for comparing bands),
`ZEROS_AS_MISSING` (raw Sentinel-5P zeros are missing retrievals, shown grey), `VMIN`/`VMAX`, `SAVE_TO`.

### Kernel not found in VS Code (Code - OSS)?
There were two causes on this machine, both fixed on 2026-10-08:
- **No kernel sources.** Code - OSS doesn't grant the Jupyter/Python extensions their proposed APIs. The log says
  `Extension 'ms-toolsai.jupyter' CANNOT use API proposal ...`.
  - Fix: `~/.vscode-oss/argv.json` now has `"enable-proposed-api": ["ms-toolsai.jupyter", "ms-toolsai.jupyter-renderers", "ms-python.python", "ms-python.vscode-python-envs"]`.
  - The backup is `argv.json.bak-2026-10-08`.
  - **Quit Code - OSS completely and start it again.** A window reload isn't enough.
- **`.venv` not auto-detected.** The Python extension's `pet` tool is missing (`spawn .../python-env-tools/bin/pet ENOENT`).
  The project's environment is therefore registered as a named kernel instead:
  ```bash
  uv run python -m ipykernel install --user --name aqf-scratch --display-name "Python (Scratch .venv)"
  ```

### Fallback without VS Code: JupyterLab in the browser
```bash
uv run --with jupyterlab jupyter lab tools/view_tiff.ipynb
```
It opens in the browser, and the kernel "Python (Scratch .venv)" is selected automatically.

### Pre-rendered outputs
`tools/outputs/` holds an executed copy (`view_tiff_executed.ipynb` / `.html`) and PNGs for a Sentinel-5P file,
a Sentinel-2 true-colour image and the three 2025 forecasts. It is regenerated locally and gitignored.
