from collections import Counter

def read_fasta(file_path):
    with open(file_path, "r") as f:
        header = f.readline().strip()
        sequence = "".join(line.strip() for line in f if not line.startswith(">"))
    return header, sequence

def process_fasta(file_path):
    header, sequence = read_fasta(file_path)
    length = len(sequence)
    alphabet = set(sequence)
    freq = Counter(sequence)

    print(f"\nHeader: {header}")
    print(f"Sequence length: {length}")
    print(f"Alphabet: {', '.join(sorted(alphabet))}")
    print("Relative frequencies:")
    for base in sorted(alphabet):
        percentage = (freq[base] / length) * 100
        print(f"  {base}: {percentage:.2f}%")

if __name__ == "__main__":
    fasta_file = input("Enter path to FASTA file: ").strip()
    process_fasta(fasta_file)
