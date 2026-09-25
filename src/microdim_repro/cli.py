from __future__ import annotations
import argparse
from pathlib import Path
import shutil

from .io import download_public_data
from .predictive import run_all_cohorts
from .niche import run_all_niche, run_sensitivity as run_niche_sensitivity
from .rank_calibration import run as run_calibration
from .consumer_resource import run as run_consumer
from .micom_proxy import run as run_micom
from .figures import figure1,figure2,figure3,figure4,figure_s1,figure_s2,figure_s3,figure_s4


def _refdir(): return Path(__file__).resolve().parent/'reference_results'

def reference_figures(outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True);r=_refdir()
    figure1(out/'Figure_1.png');figure2(r/'table_S4_cohort_benchmark.csv',out/'Figure_2.png')
    figure3(r/'sample_matched_primary.csv',out/'Figure_3.png')
    figure4(r/'fast_resource_reconstructed_curve.csv',out/'Figure_4.png')
    figure_s1(r/'sample_matched_primary.csv',out/'Figure_S1.png');figure_s2(r/'spectrum_sensitivity.csv',out/'Figure_S2.png')
    figure_s3(r/'table_S9_micom_dimensions.csv',out/'Figure_S3.png');figure_s4(r/'micom_cutoff_sensitivity.csv',out/'Figure_S4.png')

def main(argv=None):
    ap=argparse.ArgumentParser(prog='microdim-repro')
    sub=ap.add_subparsers(dest='cmd',required=True)
    p=sub.add_parser('download');p.add_argument('--data-dir',default='data');p.add_argument('--overwrite',action='store_true')
    p=sub.add_parser('cohorts');p.add_argument('--data-dir',default='data');p.add_argument('--outdir',default='results/cohorts')
    p=sub.add_parser('niche');p.add_argument('--data-dir',default='data');p.add_argument('--outdir',default='results/niche');p.add_argument('--subsamples',type=int,default=100)
    p=sub.add_parser('niche-sensitivity');p.add_argument('--data-dir',default='data');p.add_argument('--outdir',default='results/niche_sensitivity');p.add_argument('--subsamples',type=int,default=100)
    p=sub.add_parser('calibration');p.add_argument('--outdir',default='results/calibration')
    p=sub.add_parser('consumer-resource');p.add_argument('--outdir',default='results/consumer_resource')
    p=sub.add_parser('micom');p.add_argument('--data-dir',default='data');p.add_argument('--outdir',default='results/micom')
    p=sub.add_parser('figures');p.add_argument('--outdir',default='results/figures')
    p=sub.add_parser('all');p.add_argument('--data-dir',default='data');p.add_argument('--results-dir',default='results');p.add_argument('--skip-download',action='store_true')
    args=ap.parse_args(argv)
    if args.cmd=='download': download_public_data(args.data_dir,overwrite=args.overwrite)
    elif args.cmd=='cohorts': run_all_cohorts(args.data_dir,args.outdir)
    elif args.cmd=='niche': run_all_niche(args.data_dir,args.outdir,n_subsamples=args.subsamples)
    elif args.cmd=='niche-sensitivity': run_niche_sensitivity(args.data_dir,args.outdir,n_subsamples=args.subsamples)
    elif args.cmd=='calibration': run_calibration(args.outdir)
    elif args.cmd=='consumer-resource': run_consumer(args.outdir)
    elif args.cmd=='micom': run_micom(Path(args.data_dir)/'micom/exchanges_d0_high_fiber.csv',Path(args.data_dir)/'micom/vendor_diversity.csv',args.outdir)
    elif args.cmd=='figures': reference_figures(args.outdir)
    elif args.cmd=='all':
        d=Path(args.data_dir);r=Path(args.results_dir)
        if not args.skip_download: download_public_data(d)
        run_all_cohorts(d,r/'cohorts');run_all_niche(d,r/'niche');run_calibration(r/'calibration');run_consumer(r/'consumer_resource')
        run_micom(d/'micom/exchanges_d0_high_fiber.csv',d/'micom/vendor_diversity.csv',r/'micom');reference_figures(r/'figures')
