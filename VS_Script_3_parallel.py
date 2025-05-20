import os
import subprocess
import glob
import csv
import shutil
import concurrent.futures
from datetime import datetime

# Configuration constants
RECEPTOR_FILE = "/home/tuhin/drugbank_3D_PDBQT/final_IPK1.pdbqt"
LIGAND_DIR = "/home/tuhin/drugbank_3D_PDBQT/Mols_2"
OUTPUT_DIR = "/home/tuhin/drugbank_3D_PDBQT/drugbank_output_IPK1_1"
CONF_FILE = "/home/tuhin/drugbank_3D_PDBQT/VS_Config_IPK1.txt"
VINA_PATH = "/home/tuhin/drugbank_3D_PDBQT/autodock_vina_1_1_2_linux_x86/bin/vina"
FAULTY_DIR = "/home/tuhin/drugbank_3D_PDBQT/faulty_structs_4"
ERROR_LOG_FILE = os.path.join(OUTPUT_DIR, "vina_error_log.csv")

# Create output and faulty directories if needed
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FAULTY_DIR, exist_ok=True)

# Initialize error log file
if not os.path.exists(ERROR_LOG_FILE):
    with open(ERROR_LOG_FILE, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Timestamp", "Ligand", "Error"])

def run_vina_docking(receptor, ligand, config, out_pdbqt, out_log):
    """Run Vina docking process."""
    try:
        cmd = f"{VINA_PATH} --receptor {receptor} --ligand {ligand} --config {config} --out {out_pdbqt} --cpu 1"
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            timeout=500,
            stdout=open(out_log, 'w'),
            stderr=subprocess.STDOUT
        )
        return True
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Subprocess error: {e}")
    except subprocess.TimeoutExpired:
        raise TimeoutError("Timeout (500s) exceeded")

def process_ligand(idx, total, ligand_file):
    """Process a single ligand file."""
    ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
    out_pdbqt = os.path.join(OUTPUT_DIR, f"{ligand_name}_out.pdbqt")
    out_log = os.path.join(OUTPUT_DIR, f"{ligand_name}_log.txt")

    if os.path.exists(out_pdbqt) and os.path.exists(out_log):
        print(f"Skipping {ligand_name}: {idx+1}/{total} (already processed)")
        return True

    print(f"Processing {ligand_name}: {idx+1}/{total}")
    try:
        success = run_vina_docking(RECEPTOR_FILE, ligand_file, CONF_FILE, out_pdbqt, out_log)
        if not success:
            raise RuntimeError("Unknown docking failure")
        return True
    except Exception as e:
        # Log the error
        with open(ERROR_LOG_FILE, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([datetime.now().isoformat(), ligand_name, str(e)])

        # Copy faulty ligand to fault directory
        dest_path = os.path.join(FAULTY_DIR, os.path.basename(ligand_file))
        shutil.copy2(ligand_file, dest_path)
        print(f"Copied faulty ligand to: {dest_path}")
        return False

def parse_vina_results(log_file):
    """Extract lowest binding affinity from Vina log file."""
    lowest_affinity = None
    with open(log_file, 'r') as f:
        for line in f:
            if line.strip().startswith('1 '):
                first_result = line.split()
                if len(first_result) >= 2:
                    lowest_affinity = float(first_result[1])
                break
    return lowest_affinity

if __name__ == "__main__":
    ligand_files = sorted(glob.glob(os.path.join(LIGAND_DIR, "*.pdbqt")))
    total = len(ligand_files)

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=min(os.cpu_count(), 64)
    ) as executor:
        futures = {
            executor.submit(process_ligand, idx, total, lig): lig
            for idx, lig in enumerate(ligand_files)
        }

        for future in concurrent.futures.as_completed(futures):
            ligand_file = futures[future]
            try:
                future.result()
            except Exception as e:
                print(f"Exception processing {ligand_file}: {e}")

    # Post-processing results
    results = []
    for log_file in glob.glob(os.path.join(OUTPUT_DIR, "*_log.txt")):
        ligand_name = os.path.splitext(os.path.basename(log_file))[0].replace("_log", "")
        binding_affinity = parse_vina_results(log_file)
        if binding_affinity is not None:
            try:
                drugbank_id = ligand_name.split('_')[1]
            except IndexError:
                drugbank_id = ligand_name  # fallback
            results.append((drugbank_id, binding_affinity))

    # Sort and save results
    results.sort(key=lambda x: x[1])
    csv_file = os.path.join(OUTPUT_DIR, "docking_results_drugbank_.csv")
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Rank", "DrugBank ID", "Binding Affinity (kcal/mol)"])
        for idx, (drugbank_id, binding_affinity) in enumerate(results, start=1):
            writer.writerow([idx, drugbank_id, binding_affinity])

    print(f"Results have been saved to {csv_file}")
    print(f"Any errors have been logged to {ERROR_LOG_FILE}")
