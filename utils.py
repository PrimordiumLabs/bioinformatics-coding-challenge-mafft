import numpy as np
import yaml
from collections import Counter
import subprocess
import random

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align.Applications import MafftCommandline


def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def calculate_length_percentiles(records, lower_percentile=50, upper_percentile=95):
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


def downsample_reads(records, target_count=2000):
    """Randomly downsample reads to target count."""
    if len(records) <= target_count:
        print(
            f"Warning: Number of reads ({len(records)}) is less than or equal to target count ({target_count}). No downsampling performed."
        )
        return records

    return random.sample(records, target_count)
    # return records[:target_count]


def run_mafft(input_file, output_file):
    """Run MAFFT alignment on the input sequences."""
    mafft_cline = MafftCommandline(
        input=input_file, thread=4, retree=2, adjustdirection=True
    )
    stdout, stderr = mafft_cline()

    with open(output_file, "w") as f:
        f.write(stdout)

    return output_file


def run_flye_polish(consensus_fasta, reads_fasta, output_dir, threads=4):
    """Run Flye polish on the consensus sequence."""
    cmd = [
        "flye",
        "--polish-target",
        str(consensus_fasta),
        "--nano-raw",
        str(reads_fasta),
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


# def is_near_ends(read_len, read_idx, end_threshold=0.15):
#     start_threshold = end_threshold * read_len
#     end_threshold = read_len - start_threshold

#     return read_idx < start_threshold or read_idx > end_threshold


def calculate_consensus(alignment_file):
    """Calculate consensus sequence from the alignment."""
    fasta = AlignIO.read(alignment_file, "fasta")
    consensus_seq = ""
    msa_len = len(fasta[0])
    for i in range(msa_len):
        base_counts = Counter([read.seq[i] for read in fasta]).most_common()
        consensus_base = base_counts[0][0]

        # In tagmented reads, most of the reads will be missing one of the ends
        # if consensus_base == "-" and is_near_ends(msa_len, i):
        #     consensus_base = base_counts[1][0]
        consensus_seq += consensus_base
    seq = consensus_seq.replace("-", "").upper()
    return seq


def get_alignment(
    seq1: str | Seq,
    seq2: str | Seq,
    context_size=10,
    combine_threshold: int | None = 10,
    max_edit_dist_print_len: int = 50,
) -> tuple[int, str]:
    """
    Align two sequences and return str of BLAST-like diffs.

    `combine_threshold` is space between mismatches/indels required
    to combine into a single diff -- if None, no combining.
    `print_output` determines whether to print the output or return it as a string.
    """
    try:
        import edlib
    except ImportError:
        raise ImportError("`edlib` must be installed to use this function.")

    if isinstance(seq1, Seq) and isinstance(seq2, Seq):
        seq1 = str(seq1)
        seq2 = str(seq2)
    elif isinstance(seq1, str) and isinstance(seq2, str):
        pass
    else:
        raise ValueError("seq1 and seq2 must be the same type (str or Seq)")

    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # Perform the alignment
    result = edlib.align(seq1, seq2, mode="NW", task="path")

    # Extract the edit distance
    edit_distance = result["editDistance"]

    if edit_distance == 0 or edit_distance > max_edit_dist_print_len:
        return (edit_distance, "")

    alignments = edlib.getNiceAlignment(result, seq1, seq2)
    query_aligned = alignments["query_aligned"]
    target_aligned = alignments["target_aligned"]

    # Prepare visualization line
    visualization_line = "".join(
        "|" if q == t else " " if q == "-" or t == "-" else "·"
        for q, t in zip(query_aligned, target_aligned)
    )

    # Detect and combine closely spaced differences
    diff_regions = []
    current_diff_start = None
    current_diff_end = -1  # Initialize current_diff_end before the loop

    for i, (q, t) in enumerate(zip(query_aligned, target_aligned)):
        if q != t:  # Found a difference
            if current_diff_start is None:
                current_diff_start = i  # Start of a new diff region
            current_diff_end = i
        elif current_diff_start is not None:
            if combine_threshold is not None and (
                i - current_diff_end <= combine_threshold
            ):
                continue  # Extend the current diff region
            else:
                # Close off the current diff region and reset
                diff_regions.append((current_diff_start, current_diff_end))
                current_diff_start = None

    # Ensure the last diff region is added if the sequence ends with a diff
    if current_diff_start is not None:
        diff_regions.append((current_diff_start, current_diff_end))

    output = []
    # Print combined differences with context and improved visualization
    for start, end in diff_regions:
        context_start = max(0, start - context_size)
        context_end = min(len(query_aligned), end + context_size + 1)
        if start != end:
            output.append(f"Pos. {start} → {end}:")
        else:
            output.append(f"Pos. {start}:")
        output.append(f"Ref: {query_aligned[context_start:context_end]}")
        output.append(f"     {visualization_line[context_start:context_end]}")
        output.append(f"Seq: {target_aligned[context_start:context_end]}\n")

    output_str = "\n".join(output)

    return (edit_distance, output_str)
