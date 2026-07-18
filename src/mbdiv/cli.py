"""mbdiv CLI entry point."""

import argparse
import os
import sys
import time

from .config import PipelineConfig
from .step1_merge import run_step1
from .step2_normalize import run_step2
from .step3_alpha import run_step3
from .step4_beta import run_step4
from .step5_top10 import run_step5


def parse_args():
    parser = argparse.ArgumentParser(
        prog="mbdiv",
        description="Microbiome Diversity Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("raw_data", nargs="?", default=None,
                        help="Species abundance Excel file")
    parser.add_argument("meta_data", nargs="?", default=None,
                        help="Sample metadata Excel file")

    parser.add_argument("-r", "--raw", dest="raw_flag", default=None,
                        help="Species abundance Excel (alias for positional)")
    parser.add_argument("-m", "--meta", dest="meta_flag", default=None,
                        help="Sample metadata Excel (alias for positional)")

    parser.add_argument("-o", "--output", default="result",
                        help="Output directory (default: result)")
    parser.add_argument("-t", "--top", type=int, default=10,
                        help="Top-N species count (default: 10)")
    parser.add_argument("-c", "--config", default="",
                        help="YAML config file (advanced options)")

    args = parser.parse_args()

    args.raw_data = args.raw_flag or args.raw_data
    args.meta_data = args.meta_flag or args.meta_data

    if not args.raw_data or not args.meta_data:
        parser.error(
            "Both raw data and metadata files are required.\n"
            "  Usage: mbdiv species.xlsx metadata.xlsx\n"
            "     or: mbdiv -r species.xlsx -m metadata.xlsx"
        )

    return args


def build_config(args) -> PipelineConfig:
    cfg = PipelineConfig()

    if args.config and os.path.exists(args.config):
        try:
            import yaml
            with open(args.config, "r") as f:
                yml = yaml.safe_load(f) or {}
            for key, val in yml.items():
                if hasattr(cfg, key):
                    if isinstance(getattr(cfg, key), list) and isinstance(val, str):
                        setattr(cfg, key, [v.strip() for v in val.split(",")])
                    else:
                        setattr(cfg, key, val)
            print(f"  Config loaded: {args.config}")
        except ImportError:
            print("  Warning: PyYAML not installed, ignoring -c/--config")

    cfg.raw_data = args.raw_data
    cfg.meta_data = args.meta_data
    cfg.output_dir = args.output
    cfg.top_n = args.top

    return cfg


def validate_inputs(cfg: PipelineConfig):
    errors = []
    if not os.path.exists(cfg.raw_data):
        errors.append(f"Raw data file not found: {cfg.raw_data}")
    if not os.path.exists(cfg.meta_data):
        errors.append(f"Metadata file not found: {cfg.meta_data}")
    if errors:
        for e in errors:
            print(f"  ERROR: {e}")
        sys.exit(1)

    missing = []
    for pkg in ["pandas", "numpy", "scipy", "matplotlib", "seaborn", "openpyxl"]:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    if missing:
        print("Missing required packages:")
        for p in missing:
            print(f"  pip install {p}")
        sys.exit(1)

    import importlib.util
    if importlib.util.find_spec("skbio") is None:
        print("  Note: scikit-bio not installed -> PCoA uses the manual engine")
        print("        (to use scikit-bio set pcoa_engine: 'skbio' in config and install scikit-bio)")


def run_pipeline(cfg: PipelineConfig):
    start = time.time()

    print("=" * 60)
    print("  Microbiome Diversity Analysis Pipeline (mbdiv v1.0)")
    print("=" * 60)
    print(f"  Input:    {cfg.raw_data}")
    print(f"  Metadata: {cfg.meta_data}")
    print(f"  Output:   {cfg.output_dir}")
    print(f"  Top-N:    {cfg.top_n}")
    print("=" * 60)

    os.makedirs(cfg.output_dir, exist_ok=True)

    merged_path = run_step1(cfg)
    rel_path = run_step2(cfg, merged_path)
    alpha_path = run_step3(cfg, merged_path)
    dist_path, coords_path, permanova_p = run_step4(cfg, rel_path)
    bacteria_rel_path, top10_path = run_step5(cfg, merged_path)

    elapsed = time.time() - start

    report_path = os.path.join(cfg.output_dir, "pipeline_report.txt")
    with open(report_path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("Pipeline Report\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Input:\n  Raw data:  {cfg.raw_data}\n  Metadata:  {cfg.meta_data}\n\n")
        f.write(f"Output:\n  Base dir:  {cfg.output_dir}\n\n")
        f.write("Results:\n")
        f.write(f"  Step 1 - Merged species:  {merged_path}\n")
        f.write(f"  Step 2 - Relative abund: {rel_path}\n")
        f.write(f"  Step 3 - Alpha diversity: {alpha_path}\n")
        f.write(f"  Step 4 - Beta distance:   {dist_path}\n")
        f.write(f"  Step 4 - PCoA coords:     {coords_path}\n")
        if permanova_p is not None:
            f.write(f"  Step 4 - PERMANOVA p:     {permanova_p}\n")
        f.write(f"  Step 5 - Bacteria abund:  {bacteria_rel_path}\n")
        f.write(f"  Step 5 - Top species:     {top10_path}\n\n")
        f.write(f"Elapsed: {elapsed:.1f}s\n")

    print("\n" + "=" * 60)
    print(f"  Pipeline complete! ({elapsed:.1f}s)")
    print(f"  Report: {report_path}")
    print("=" * 60)

    return report_path


def main():
    args = parse_args()
    cfg = build_config(args)
    validate_inputs(cfg)
    run_pipeline(cfg)


if __name__ == "__main__":
    main()
