import os
import subprocess
import sys

def run_command(cmd, shell=True):
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
    return result.returncode

def check_file(filepath):
    return os.path.exists(filepath) and os.path.getsize(filepath) > 0

rna_runs = [f"SRR833099{i}" for i in range(0, 9)]
wgbs_runs = [f"SRR83496{i}" for i in range(15, 24)]

os.makedirs("results/fastq/raw", exist_ok=True)
os.makedirs("temp/sra", exist_ok=True)
os.makedirs("temp/fastq", exist_ok=True)

failures = []

rna_completed = 0
for run in rna_runs:
    out_file = f"results/fastq/raw/{run}_1.fastq.gz"
    if check_file(out_file):
        print(f"Skipping {run}, output exists.")
        rna_completed += 1
        continue
    print(f"Processing {run} (RNA-seq SINGLE)...")
    temp_sra = f"temp/sra/{run}"
    temp_fastq = f"temp/fastq/{run}"
    os.makedirs(temp_sra, exist_ok=True)
    os.makedirs(temp_fastq, exist_ok=True)
    if run_command(f"prefetch {run} -O {temp_sra}") != 0:
        failures.append(run)
        continue
    sra_file = f"{temp_sra}/{run}/{run}.sra"
    if not os.path.exists(sra_file):
        sra_file = f"{temp_sra}/{run}.sra"
    if run_command(f"fasterq-dump --split-files {sra_file} -O {temp_fastq}") != 0:
        failures.append(run)
        continue
    fq_options = [f"{temp_fastq}/{run}_1.fastq", f"{temp_fastq}/{run}.fastq"]
    fq_path = next((opt for opt in fq_options if os.path.exists(opt)), None)
    if not fq_path:
        failures.append(run)
        continue
    if run_command(f"pigz -c {fq_path} > {out_file}") != 0:
        failures.append(run)
        continue
    rna_completed += 1

wgbs_completed = 0
for run in wgbs_runs:
    out1, out2 = f"results/fastq/raw/{run}_1.fastq.gz", f"results/fastq/raw/{run}_2.fastq.gz"
    if check_file(out1) and check_file(out2):
        print(f"Skipping {run}, outputs exist.")
        wgbs_completed += 1
        continue
    print(f"Processing {run} (WGBS PAIRED)...")
    temp_sra, temp_fastq = f"temp/sra/{run}", f"temp/fastq/{run}"
    os.makedirs(temp_sra, exist_ok=True); os.makedirs(temp_fastq, exist_ok=True)
    if run_command(f"prefetch {run} -O {temp_sra}") != 0:
        failures.append(run); continue
    sra_file = f"{temp_sra}/{run}/{run}.sra"
    if not os.path.exists(sra_file): sra_file = f"{temp_sra}/{run}.sra"
    if run_command(f"fasterq-dump --split-files {sra_file} -O {temp_fastq}") != 0:
        failures.append(run); continue
    fq1, fq2 = f"{temp_fastq}/{run}_1.fastq", f"{temp_fastq}/{run}_2.fastq"
    if not (os.path.exists(fq1) and os.path.exists(fq2)):
        failures.append(run); continue
    if run_command(f"pigz -c {fq1} > {out1}") != 0 or run_command(f"pigz -c {fq2} > {out2}") != 0:
        failures.append(run); continue
    wgbs_completed += 1

print(f"\n--- Summary ---\nRNA-seq single-end files: {rna_completed}\nWGBS pairs: {wgbs_completed}")
if failures: print(f"Failures: {', '.join(failures)}")
else: print("No failures.")
print("\nFiles in results/fastq/raw:")
for f in sorted(os.listdir("results/fastq/raw")): print(f"results/fastq/raw/{f}")
