import sys

sys.path.append("..")

from pathlib import Path
from Bio import SeqIO
import edlib

from main import run_pipeline
from utils import load_config


def is_same_sequence(seq1: str, seq2: str) -> bool:
    return edlib.align(seq1, seq2, mode="NW", task="path")["editDistance"] < 2


def test_pipeline():
    fastq = Path("data/test1/test1_reads.fastq")
    ref = Path("data/test1/test1_ref.fasta")
    output = Path("data/test1/test1_output")
    output.mkdir(parents=True, exist_ok=True)

    config = load_config(Path("../config.yml"))
    output = run_pipeline(fastq, output, config)
    assert output.exists()

    output_seq = SeqIO.read(output, "fasta")
    ref_seq = SeqIO.read(ref, "fasta")

    assert is_same_sequence(output_seq.seq, ref_seq.seq)


if __name__ == "__main__":
    test_pipeline()
    print("Tests passed!")
