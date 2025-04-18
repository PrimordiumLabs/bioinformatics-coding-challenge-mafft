import argparse
from pathlib import Path
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from utils import (
    calculate_length_percentiles,
    filter_reads_by_length,
    downsample_reads,
    run_mafft,
    calculate_consensus,
    run_flye_polish,
    load_config,
)

# TODO: move things to config
# TODO: for the ends, be more permissive when making consensus


def run_pipeline(
    input_fasta: Path,
    output_dir: Path,
    min_read_percentile: float = 50,
    max_read_percentile: float = 95,
    read_cap: int = 500,
) -> Path:
    # Get paths from config
    intermediate_dir = output_dir / "intermediate"
    step = 0

    # Create output directory if it doesn't exist
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    # Read all records first to calculate percentiles
    records = list(SeqIO.parse(input_fasta, "fasta"))

    # Calculate length thresholds
    lower_threshold, upper_threshold = calculate_length_percentiles(
        records,
        min_read_percentile,
        max_read_percentile,
    )

    print(
        f"Length thresholds: lower={lower_threshold:.2f}, upper={upper_threshold:.2f}"
    )

    # Filter records
    filtered_records = filter_reads_by_length(records, lower_threshold, upper_threshold)
    filtered_fasta = intermediate_dir / f"{step}_filtered_{input_fasta.name}"
    SeqIO.write(filtered_records, filtered_fasta, "fasta")
    step += 1

    print(f"Original number of reads: {len(records)}")
    print(f"Filtered number of reads: {len(filtered_records)}")
    print(f"Filtered reads written to: {filtered_fasta}")

    # Downsample reads
    downsampled_records = downsample_reads(filtered_records, read_cap)
    downsampled_fasta = intermediate_dir / f"{step}_downsampled_{input_fasta.name}"
    SeqIO.write(downsampled_records, downsampled_fasta, "fasta")
    step += 1
    print(f"Downsampled to {len(downsampled_records)} reads")
    print(f"Downsampled reads written to: {downsampled_fasta}")

    # Run MAFFT alignment
    aligned_fasta = intermediate_dir / f"{step}_aligned_{input_fasta.stem}.fasta"
    print("Running MAFFT alignment...")
    run_mafft(downsampled_fasta, aligned_fasta)
    print(f"Alignment written to: {aligned_fasta}")
    step += 1

    # Calculate consensus and write as FASTA using BioPython
    consensus = calculate_consensus(aligned_fasta)
    consensus_record = SeqRecord(
        seq=Seq(consensus),
        id=f"consensus_{input_fasta.stem}",
        description="Consensus sequence from MAFFT alignment",
    )
    consensus_file = intermediate_dir / f"{step}_consensus_{input_fasta.stem}.fasta"
    SeqIO.write(consensus_record, consensus_file, "fasta")
    print(f"Consensus sequence written to: {consensus_file}")
    step += 1

    # Run Flye polish
    flye_output_dir = intermediate_dir / f"{step}_flye_polish"
    flye_output_dir.mkdir(exist_ok=True)
    polished_consensus = run_flye_polish(
        consensus_file, downsampled_fasta, flye_output_dir, threads=4
    )
    print(f"Polished consensus written to: {polished_consensus}")
    step += 1

    # Copy polished consensus to output directory
    polished_consensus_file = output_dir / "assembly.fasta"
    shutil.copy(polished_consensus, polished_consensus_file)
    print(f"Polished consensus written to: {polished_consensus_file}")

    return polished_consensus_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter FASTA sequences by length percentiles"
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to input FASTA file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default="output",
        help="Path to output directory (default: output)",
    )

    # Configs

    # Length bounds
    parser.add_argument(
        "--min-read-percentile",
        type=float,
        default=50,
        help="Minimum read percentile (default: 50)",
    )
    parser.add_argument(
        "--max-read-percentile",
        type=float,
        default=95,
        help="Maximum read percentile (default: 95)",
    )

    # Read cap
    parser.add_argument(
        "--read-cap",
        type=int,
        default=2000,
        help="Read cap (default: 2000)",
    )

    args = parser.parse_args()

    # Load configuration and run pipeline
    run_pipeline(
        args.input,
        args.output,
        args.min_read_percentile,
        args.max_read_percentile,
        args.read_cap,
    )
