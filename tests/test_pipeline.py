import sys

sys.path.append("..")

from pathlib import Path
from Bio import SeqIO
import edlib

from main import run_pipeline

# from solution import run_pipeline
from utils import load_config, get_alignment


def is_same_sequence(seq1: str, seq2: str) -> bool:
    return edlib.align(seq1, seq2, mode="NW", task="path")["editDistance"] < 2


def test_pipeline_1():
    fasta = Path("data/test1/test1_reads.fasta")
    ref = Path("data/test1/test1_ref.fasta")
    output = Path("data/test1/test1_output")
    output.mkdir(parents=True, exist_ok=True)

    output = run_pipeline(fasta, output)
    assert output.exists()

    output_seq = SeqIO.read(output, "fasta")
    ref_seq = SeqIO.read(ref, "fasta")

    if not is_same_sequence(output_seq.seq, ref_seq.seq):
        dist, alignment = get_alignment(
            output_seq.seq, ref_seq.seq, max_edit_dist_print_len=1000
        )
        print(
            f"Alignment (edit distance: {dist})\n",
            alignment,
        )
        assert False, "Sequences are not the same"


if __name__ == "__main__":
    test_pipeline_1()
    print("Tests passed!")
