from pathlib import Path
from datetime import datetime
import shutil
import subprocess

project = Path(__file__).resolve().parent.parent
raw = project / "data/raw/kang18_8vs8.rds"

if not raw.is_file():
    raise SystemExit("Başlangıç veri dosyası bulunamadı.")

names = [
    "check_environment.R",
    "01_inspect_kang_metadata.R",
    "02_audit_cell_coverage.R",
    "03_build_monocyte_pseudobulk.R",
    "04_filter_and_check_design.R",
    "05_fit_primary_deseq2.R",
    "06_pseudobulk_sample_qc.R",
    "07_plot_paired_genes.R",
    "08_leave_one_donor_out.R",
    "09_summarize_lodo_genes.R",
    "10_review_sensitive_genes.R",
]

for name in names:
    if not (project / "scripts" / name).is_file():
        raise SystemExit(f"Script bulunamadı: {name}")

if shutil.which("Rscript") is None:
    raise SystemExit("Rscript bulunamadı.")

stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
run_dir = project / "repro_checks" / stamp

for folder in [
    "scripts", "docs", "logs", "data/raw", "data/processed",
    "results/tables", "results/figures"
]:
    (run_dir / folder).mkdir(parents=True, exist_ok=True)

shutil.copy2(raw, run_dir / "data/raw" / raw.name)

for name in names:
    shutil.copy2(project / "scripts" / name, run_dir / "scripts" / name)

(project / "repro_checks/latest_run.txt").write_text(
    str(run_dir) + "\n"
)

print(f"Run directory: {run_dir}", flush=True)

for name in names:
    log = run_dir / "logs" / f"{Path(name).stem}.log"
    print(f"START: {name}", flush=True)

    with log.open("w") as handle:
        result = subprocess.run(
            ["Rscript", f"scripts/{name}"],
            cwd=run_dir,
            stdout=handle,
            stderr=subprocess.STDOUT,
            env=None,
        )

    if result.returncode != 0:
        print(f"FAILED: {name}\nLog: {log}", flush=True)
        lines = log.read_text(errors="replace").splitlines()
        print("\n".join(lines[-40:]), flush=True)
        raise SystemExit(result.returncode)

    print(f"COMPLETED: {name}", flush=True)

(run_dir / "execution_completed.txt").write_text(
    "All analysis scripts completed successfully.\n"
)

print("\nALL ANALYSIS SCRIPTS COMPLETED.", flush=True)
print("Numerical comparison is still required.", flush=True)
