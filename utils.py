import numpy as np
import yaml
from collections import Counter
import subprocess

from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline


def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def calculate_length_percentiles(records, lower_percentile, upper_percentile):
    """Calculate length thresholds based on percentiles."""
    lengths = [len(record.seq) for record in records]
    lower_threshold = np.percentile(lengths, lower_percentile)
    upper_threshold = np.percentile(lengths, upper_percentile)
    return lower_threshold, upper_threshold


def filter_reads_by_length(records, lower_threshold, upper_threshold):
    """Filter reads based on length thresholds."""
    filtered_records = []
    for record in records:
        length = len(record.seq)
        if lower_threshold <= length <= upper_threshold:
            filtered_records.append(record)
    return filtered_records


def downsample_reads(records, target_count=500):
    """Randomly downsample reads to target count."""
    if len(records) <= target_count:
        print(
            f"Warning: Number of reads ({len(records)}) is less than or equal to target count ({target_count}). No downsampling performed."
        )
        return records

    # return random.sample(records, target_count)
    return records[:target_count]


def run_mafft(input_file, output_file):
    """Run MAFFT alignment on the input sequences."""
    mafft_cline = MafftCommandline(
        input=input_file, thread=4, retree=2, adjustdirection=True
    )
    stdout, stderr = mafft_cline()

    with open(output_file, "w") as f:
        f.write(stdout)

    return output_file


def run_flye_polish(consensus_fasta, reads_fastq, output_dir, threads=4):
    """Run Flye polish on the consensus sequence."""
    cmd = [
        "flye",
        "--polish-target",
        str(consensus_fasta),
        "--nano-raw",
        str(reads_fastq),
        "--iterations",
        "1",
        "--out-dir",
        str(output_dir),
        "-t",
        str(threads),
    ]

    print("Running Flye polish...")
    subprocess.run(cmd, check=True)

    # Return path to polished consensus
    polished_consensus = output_dir / "polished_1.fasta"
    return polished_consensus


def calculate_consensus(alignment_file):
    """Calculate consensus sequence from the alignment."""
    fasta = AlignIO.read(alignment_file, "fasta")
    consensus_seq = ""
    msa_len = len(fasta[0])
    for i in range(msa_len):
        consensus_base = Counter([read.seq[i] for read in fasta]).most_common(1)[0][0]
        consensus_seq += consensus_base
    seq = consensus_seq.replace("-", "").upper()
    return seq
