import subprocess
import sys
import os

INPUT_FASTA = sys.argv[1]

MAFFT_PATH = "mafft"


def run_mafft(input_fasta, mafft_path):
    if not os.path.exists(input_fasta):
        print(f"Error: file not found: {input_fasta}")
        sys.exit(1)

    output_fasta = f"aligned_{input_fasta}"

    cmd = [mafft_path, "--auto", input_fasta]
    print(f"Running MAFFT on {input_fasta}...")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print("MAFFT ERROR:\n", result.stderr)
        sys.exit(1)

    with open(output_fasta, "w") as f:
        f.write(result.stdout)

    print(f"\n✅ Alignment complete: {output_fasta}")
    return output_fasta


if __name__ == "__main__":
    run_mafft(INPUT_FASTA, MAFFT_PATH)
