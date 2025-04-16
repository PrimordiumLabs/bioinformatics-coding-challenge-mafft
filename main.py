import argparse
from pathlib import Path
import numpy as np
import yaml
import random
from collections import Counter

from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

TMP_DIR = Path("./tmp")
TMP_DIR.mkdir(exist_ok=True)


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


def main(config):
    """Run the filtering pipeline with the given configuration."""
    # Get paths from config
    input_fastq = Path(config["input"]["fastq"])
    output_dir = Path(config["input"]["output_dir"])
    intermediate_dir = output_dir / "intermediate"

    # Create output directory if it doesn't exist
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    # Read all records first to calculate percentiles
    records = list(SeqIO.parse(input_fastq, "fastq"))

    # Calculate length thresholds
    lower_threshold, upper_threshold = calculate_length_percentiles(
        records,
        config["filtering"]["lower_percentile"],
        config["filtering"]["upper_percentile"],
    )

    print(
        f"Length thresholds: lower={lower_threshold:.2f}, upper={upper_threshold:.2f}"
    )

    # Filter records
    filtered_records = filter_reads_by_length(records, lower_threshold, upper_threshold)

    # Write filtered records to output file
    filtered_fastq = intermediate_dir / f"filtered_{input_fastq.name}"
    SeqIO.write(filtered_records, filtered_fastq, "fastq")

    print(f"Original number of reads: {len(records)}")
    print(f"Filtered number of reads: {len(filtered_records)}")
    print(f"Filtered reads written to: {filtered_fastq}")

    # Downsample reads
    downsampled_records = downsample_reads(filtered_records)
    downsampled_fastq = intermediate_dir / f"downsampled_{input_fastq.name}"
    SeqIO.write(downsampled_records, downsampled_fastq, "fastq")
    print(f"Downsampled to {len(downsampled_records)} reads")
    print(f"Downsampled reads written to: {downsampled_fastq}")

    # Convert FASTQ to FASTA for MAFFT
    downsampled_fasta = intermediate_dir / f"downsampled_{input_fastq.stem}.fasta"
    SeqIO.convert(downsampled_fastq, "fastq", downsampled_fasta, "fasta")

    # Run MAFFT alignment
    aligned_fasta = intermediate_dir / f"aligned_{input_fastq.stem}.fasta"
    print("Running MAFFT alignment...")
    run_mafft(downsampled_fasta, aligned_fasta)
    print(f"Alignment written to: {aligned_fasta}")

    # Calculate consensus and write as FASTA using BioPython
    consensus = calculate_consensus(aligned_fasta)
    consensus_record = SeqRecord(
        seq=Seq(consensus),
        id=f"consensus_{input_fastq.stem}",
        description="Consensus sequence from MAFFT alignment",
    )
    consensus_file = output_dir / f"consensus_{input_fastq.stem}.fasta"
    SeqIO.write(consensus_record, consensus_file, "fasta")
    print(f"Consensus sequence written to: {consensus_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter FASTQ reads by length percentiles"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default="config.yml",
        help="Path to configuration file (default: config.yml)",
    )
    args = parser.parse_args()

    # Load configuration and run pipeline
    config = load_config(args.config)
    main(config)
