"""
NOTE:
- should we explain the library prep in the rapid barcode protocol so they expect to see tagmented reads from a pcr product?
"""

import argparse
from pathlib import Path
import shutil
from typing import List
from collections import Counter
from tqdm import tqdm

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


def filter_background_reads(records: List[SeqRecord]) -> List[SeqRecord]:
    KMER_SIZE = 15
    COMMON_KMER_CUTOFF = 500
    MIN_OVERLAPPING_KMER_THRESHOLD = 50
    print(
        "Filtering background reads based on common kmers. Kmer size:",
        KMER_SIZE,
        "Common kmer cutoff:",
        COMMON_KMER_CUTOFF,
        "Minimum overlapping kmer threshold:",
        MIN_OVERLAPPING_KMER_THRESHOLD,
    )

    # Get kmer counts per read
    read_kmer_counts = {}
    total_kmers = []
    for record in tqdm(records):
        read_kmers = []
        for i in range(len(record.seq) - KMER_SIZE + 1):
            kmer = record.seq[i : i + KMER_SIZE]
            read_kmers.append(kmer)
            total_kmers.append(kmer)
        read_kmer_counts[record.id] = Counter(read_kmers)

    # Get total unique kmers
    total_kmer_counts = Counter(total_kmers)
    print("Total unique kmers:", len(total_kmer_counts))

    # Identify the most common kmers
    most_common_kmers_counts = total_kmer_counts.most_common(COMMON_KMER_CUTOFF)
    most_common_kmers = set([kmer for kmer, _ in most_common_kmers_counts])
    print("Most common kmers:", len(most_common_kmers))

    # For each read, how many of its kmers are "globally common"
    read_common_kmer_overlap = {}
    for read_id, kmer_counts in tqdm(read_kmer_counts.items()):
        read_common_kmer_overlap[read_id] = len(
            set(kmer_counts.keys()) & most_common_kmers
        )

    # Filter reads with too few common kmers
    filtered_reads = [
        record
        for record in records
        if read_common_kmer_overlap[record.id] >= MIN_OVERLAPPING_KMER_THRESHOLD
    ]

    return filtered_reads


def run_pipeline(input_fastq: Path, output_dir: Path, config: dict) -> Path:
    # Get paths from config
    intermediate_dir = output_dir / "intermediate"
    step = 0

    # Create output directory if it doesn't exist
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    # Read all records first to calculate percentiles
    records = list(SeqIO.parse(input_fastq, "fastq"))
    print(f"Original number of reads: {len(records)}")

    # Filter background reads
    filtered_records = filter_background_reads(records)
    filtered_fastq = intermediate_dir / f"{step}_filtered_background_{input_fastq.name}"
    SeqIO.write(filtered_records, filtered_fastq, "fastq")
    step += 1
    print(f"Filtered background reads: {len(filtered_records)}")
    print(f"Filtered background reads written to: {filtered_fastq}")

    # Calculate length thresholds
    lower_threshold, upper_threshold = calculate_length_percentiles(
        filtered_records,
        config["filtering"]["lower_percentile"],
        config["filtering"]["upper_percentile"],
    )
    print(
        f"Length thresholds: lower={lower_threshold:.2f}, upper={upper_threshold:.2f}"
    )

    # Filter records
    filtered_records = filter_reads_by_length(
        filtered_records, lower_threshold, upper_threshold
    )
    filtered_fastq = intermediate_dir / f"{step}_filtered_{input_fastq.name}"
    SeqIO.write(filtered_records, filtered_fastq, "fastq")
    step += 1
    print(f"Filtered number of reads: {len(filtered_records)}")
    print(f"Filtered reads written to: {filtered_fastq}")

    # Downsample reads
    downsampled_records = downsample_reads(filtered_records)
    downsampled_fastq = intermediate_dir / f"{step}_downsampled_{input_fastq.name}"
    SeqIO.write(downsampled_records, downsampled_fastq, "fastq")
    step += 1
    print(f"Downsampled to {len(downsampled_records)} reads")
    print(f"Downsampled reads written to: {downsampled_fastq}")

    # Convert FASTQ to FASTA for MAFFT
    downsampled_fasta = (
        intermediate_dir / f"{step}_downsampled_{input_fastq.stem}.fasta"
    )
    SeqIO.convert(downsampled_fastq, "fastq", downsampled_fasta, "fasta")
    step += 1

    # Run MAFFT alignment
    aligned_fasta = intermediate_dir / f"{step}_aligned_{input_fastq.stem}.fasta"
    print("Running MAFFT alignment...")
    run_mafft(downsampled_fasta, aligned_fasta)
    print(f"Alignment written to: {aligned_fasta}")
    step += 1

    # Calculate consensus and write as FASTA using BioPython
    consensus = calculate_consensus(aligned_fasta)
    consensus_record = SeqRecord(
        seq=Seq(consensus),
        id=f"consensus_{input_fastq.stem}",
        description="Consensus sequence from MAFFT alignment",
    )
    consensus_file = intermediate_dir / f"{step}_consensus_{input_fastq.stem}.fasta"
    SeqIO.write(consensus_record, consensus_file, "fasta")
    print(f"Consensus sequence written to: {consensus_file}")
    step += 1

    # Run Flye polish
    flye_output_dir = intermediate_dir / f"{step}_flye_polish"
    flye_output_dir.mkdir(exist_ok=True)
    polished_consensus = run_flye_polish(
        consensus_file, downsampled_fastq, flye_output_dir, threads=4
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
        description="Filter FASTQ reads by length percentiles"
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to input FASTQ file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default="output",
        help="Path to output directory (default: output)",
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
    run_pipeline(args.input, args.output, config)
