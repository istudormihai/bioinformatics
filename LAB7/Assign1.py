import matplotlib.pyplot as plt

def find_repetitions(dna_sequence, min_length=6, max_length=10):
    all_repetitions = {}
    sequence_upper = dna_sequence.upper()
    
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence_upper) - length + 1):
            pattern = sequence_upper[i:i + length]
            if all(nucleotide in 'ATGC' for nucleotide in pattern):
                if pattern not in all_repetitions:
                    all_repetitions[pattern] = []
                all_repetitions[pattern].append(i)
    
    repetitions = {pattern: positions for pattern, positions in all_repetitions.items() 
                   if len(positions) > 1}
    
    return repetitions


def main():
    dna_sequence = """
    AGTTGAAGTTGTGGAATATCTTGTTTATGCAACATGATTGTATCATGGTTTTGGCTCAAACATGCAATTA
    TTTAGTTAACAGTACATGATATTCTTTTGTTAATGTGATAACCCCTGAGGTTCCTTTTTTACAGTGCAGG
    AAGAGGCGAAGGACGGGGGAAGATATTTAAAGAGCTTCAGACACGGATTTTGCACTTCCCTCCACAATGT
    CTGCACAGACTCTGCTGCACCTCTTGGCTCTGGTGGCCTTTTTCTTTGCAGGATCCGATGCACAACTTCA
    TGAATGTGGCTTAGCTTCTCCCAACTTCAAAATAGTTGGAGGTCAGGACGGCTCACCTGGAAGTTGGCCC
    TGGCAGGTGAAGCTTTATGGCCCGTTCGGGTGCGCAGGCTCCCTGATCAACAAAAATTGGGTTCTGACTG
    CAGCTCACTGCGTCTCTAGAACAACTCCGTCCATGTGGACGGTGGTTCTGGGGCAGCAGAATCTGAGCGA
    CACACAGATGAGAACCGGAGTGAGAAGGAACGTTGGAAGAATCATCGTACACCCCAAGTTCAATCCCTCC
    CCCATCGACAACGACATTGCGTTGCTCAGAAAAAACAACAAAATAGTAATTAATATTTAAGGAAAAACAT
    TTTAAATGTTTGATTTTAATATTTACTTGTTTGCTCTTCTGTTATT
    """

    dna_sequence = "".join(char for char in dna_sequence if char.upper() in 'ATGC')
    
    print(f"\nSequence length: {len(dna_sequence)} nucleotides")
    print("Searching for repetitions...")

    repetitions = find_repetitions(dna_sequence)

    print(f"Found {len(repetitions)} repeated patterns\n")

    sorted_patterns = sorted(repetitions.items(), 
                             key=lambda x: len(x[1]), 
                             reverse=True)
    
    for pattern, positions in sorted_patterns[:20]:
        print(f"{pattern} ({len(pattern)}bp): {len(positions)} occurrences")
    
    if sorted_patterns:
        top_n = 20
        top_patterns = [p for p, _ in sorted_patterns[:top_n]]
        top_counts = [len(pos) for _, pos in sorted_patterns[:top_n]]

        plt.figure(figsize=(12, 6))
        plt.bar(range(len(top_patterns)), top_counts, color='teal')
        plt.xticks(range(len(top_patterns)), top_patterns, rotation=90)
        plt.xlabel("DNA Pattern (6–10 bp)")
        plt.ylabel("Occurrences")
        plt.title(f"Top {top_n} Most Frequent DNA Repetitions")
        plt.tight_layout()
        plt.show()
    else:
        print("No repetitions found.")

if __name__ == "__main__":
    main()
