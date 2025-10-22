from collections import Counter
import matplotlib.pyplot as plt

genetic_code = {
    'UUU':'F','UUC':'F','UUA':'L','UUG':'L',
    'CUU':'L','CUC':'L','CUA':'L','CUG':'L',
    'AUU':'I','AUC':'I','AUA':'I','AUG':'M',
    'GUU':'V','GUC':'V','GUA':'V','GUG':'V',
    'UCU':'S','UCC':'S','UCA':'S','UCG':'S',
    'CCU':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACU':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCU':'A','GCC':'A','GCA':'A','GCG':'A',
    'UAU':'Y','UAC':'Y','UAA':'Stop','UAG':'Stop',
    'CAU':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAU':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAU':'D','GAC':'D','GAA':'E','GAG':'E',
    'UGU':'C','UGC':'C','UGA':'Stop','UGG':'W',
    'CGU':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGU':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGU':'G','GGC':'G','GGA':'G','GGG':'G'
}

def read_fasta(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return seq

def count_codons(dna_seq):
    rna_seq = dna_seq.upper().replace('T', 'U')
    codons = [rna_seq[i:i+3] for i in range(0, len(rna_seq)-2, 3)]
    codons = [c for c in codons if len(c) == 3]
    return Counter(codons)

def plot_top_codons(counter, title):
    top10 = counter.most_common(10)
    codons, counts = zip(*top10)
    plt.bar(codons, counts, color='skyblue', edgecolor='black')
    plt.title(title)
    plt.xlabel("Codon")
    plt.ylabel("Frequency")
    plt.show()
    return top10

def count_amino_acids(counter):
    aa_count = Counter()
    for codon, freq in counter.items():
        aa = genetic_code.get(codon, '')
        if aa not in ('Stop', ''):
            aa_count[aa] += freq
    return aa_count


covid_fasta = "covid.fna"     
influenza_fasta = "influenza1.fna" 

# 1. Read both sequences
covid_seq = read_fasta(covid_fasta)
influenza_seq = read_fasta(influenza_fasta)

# 2. Count codons
covid_codons = count_codons(covid_seq)
influenza_codons = count_codons(influenza_seq)

# 3. Plot top 10 codons for each
top10_covid = plot_top_codons(covid_codons, "Top 10 Codons - COVID-19")
top10_influenza = plot_top_codons(influenza_codons, "Top 10 Codons - Influenza")

# 4. Compare common codons
covid_top_set = set(c[0] for c in top10_covid)
influenza_top_set = set(c[0] for c in top10_influenza)
common_codons = covid_top_set & influenza_top_set

print("\nCommon most frequent codons between COVID-19 and Influenza:")
print(common_codons if common_codons else "None found")

# 5. Top 3 amino acids per genome
covid_aa = count_amino_acids(covid_codons)
influenza_aa = count_amino_acids(influenza_codons)

print("\nTop 3 amino acids in COVID-19 genome:")
for aa, freq in covid_aa.most_common(3):
    print(f"{aa}: {freq}")

print("\nTop 3 amino acids in Influenza genome:")
for aa, freq in influenza_aa.most_common(3):
    print(f"{aa}: {freq}")
