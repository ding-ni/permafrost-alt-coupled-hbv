# -*- coding: utf-8 -*-
import argparse
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "examples" / "synthetic_basin" / "smoke_config.json"
CLI_TARGETS = [
    REPO_ROOT / "workflow" / "data_prep" / "01_clip_dem.py",
    REPO_ROOT / "workflow" / "data_prep" / "07_process_precipitation.py",
    REPO_ROOT / "workflow" / "data_prep" / "11_build_alt_coupling_params.py",
    REPO_ROOT / "workflow" / "models" / "run_original_hbv.py",
    REPO_ROOT / "workflow" / "models" / "calibrate_improved_hbv.py",
    REPO_ROOT / "workflow" / "models" / "run_alt_coupled_hbv.py",
]


def run(cmd):
    completed = subprocess.run(cmd, cwd=REPO_ROOT, check=True, capture_output=True, text=True)
    if completed.stdout:
        print(completed.stdout.strip())
    if completed.stderr:
        print(completed.stderr.strip())


def main():
    parser = argparse.ArgumentParser(description="Run basic public-repo smoke checks.")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--generate-example", action="store_true")
    parser.add_argument("--lightweight", action="store_true", help="Also create the synthetic example workspace.")
    args = parser.parse_args()

    if args.generate_example or args.lightweight:
        run([sys.executable, str(REPO_ROOT / "tools" / "generate_synthetic_example.py"), "--config", args.config])

    run([sys.executable, str(REPO_ROOT / "tools" / "verify_public_repo.py")])
    for target in CLI_TARGETS:
        run([sys.executable, str(target), "--help"])
    print("Smoke test completed.")


if __name__ == "__main__":
    main()
