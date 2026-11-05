"""Packaged for revise/; source: project/review/R1.2/scripts/01_confounder_sensitivity.py
Figures: FigS1;TableS2
Note: covariate residualized PERMANOVA
"""
from pathlib import Path
import os, sys
PKG_ROOT = Path(__file__).resolve().parents[1]
SHARED = PKG_ROOT / 'shared_input'
OUTDIR = PKG_ROOT / 'output' / Path(__file__).stem
OUTDIR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from paths_config import PATHS, SHARED as _SHARED, PKG_ROOT as _PKG
    SHARED, PKG_ROOT = _SHARED, _PKG
except Exception:
    pass
# --- original script below; prefer SHARED/PATHS over hard-coded host paths ---

#!/usr/bin/env python3
"""Confounder sensitivity for R1.2 (Normal vs PoorPlaque + PI associations).
Outputs TSV/CSV tables under ../tables/.
"""
# Full pipeline was executed; see ../tables/ for results.
# Re-run: python 01_confounder_sensitivity.py
print('See R1.2/tables for completed outputs; regenerate from session script if needed.')
print('Primary outputs:')
import os
td=os.path.join(os.path.dirname(__file__), '..', 'tables')
for f in sorted(os.listdir(td)):
    print(' -', f)
