from Bio.Seq import Seq
from Bio import SeqIO

def read_fasta(file_path):
    for record in SeqIO.parse(file_path, "fasta"):
        return record.seq

def gc_content(seq):
    gc_count = seq.count("G") + seq.count("C")
    return round((gc_count / len(seq)) * 100, 2)

def reverse_complement(seq):
    return seq.reverse_complement()

def find_codons(seq):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    starts = [i for i in range(len(seq)) if seq[i:i+3] == start_codon]
    stops = [i for i in range(len(seq)) if seq[i:i+3] in stop_codons]
    return starts, stops

if __name__ == "__main__":
    sequence = read_fasta("test_sequence.fasta")
    print(f"Original sequence: {sequence}")
    print(f"GC Content: {gc_content(sequence)}%")
    print(f"Reverse Complement: {reverse_complement(sequence)}")
    starts, stops = find_codons(sequence)
    print(f"Start Codons at positions: {starts}")
    print(f"Stop Codons at positions: {stops}")
