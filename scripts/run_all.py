#!/usr/bin/env python3
from pathlib import Path
import argparse, subprocess, sys
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser()
ap.add_argument('--frozen',action='store_true',help='Verify archived numerical tables and rebuild manuscript figures (default).')
ap.add_argument('--public-reference',action='store_true',help='Run Python analyses requiring public inputs. See docs/EXECUTION.md for externally maintained kinetic models.')
a=ap.parse_args()
if not a.public_reference:
    for s in ['validate_frozen.py','rebuild_figures_from_frozen.py']:
        print('\n###',s); subprocess.run([sys.executable,str(ROOT/'scripts'/s)],check=True)
    subprocess.run([sys.executable,'-m','pytest','-q'],cwd=ROOT,check=True,env={**__import__('os').environ,'PYTHONPATH':str(ROOT/'src')})
else:
    for s in ['run_cohort_benchmark.py','run_niche_dimension.py','run_rank_calibration.py','run_consumer_resource_verification.py']:
        print('\n###',s); subprocess.run([sys.executable,str(ROOT/'scripts'/s)],check=True)
    print('\nhCom2* and Ho kinetic-model execution uses the pinned public model programs; see docs/EXECUTION.md.')
