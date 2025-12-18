import pandas as pd
import numpy as np

motifs = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT"
]

alphabet = ["A", "C", "G", "T"]
motif_length = len(motifs[0])

count_matrix = pd.DataFrame(
    0, index=alphabet, columns=range(motif_length)
)

for motif in motifs:
    for i, base in enumerate(motif):
        count_matrix.loc[base, i] += 1

print("\nCount Matrix")
print(count_matrix)

pseudocount = 1
weight_matrix = count_matrix + pseudocount

print("\nWeight Matrix (with pseudocounts)")
print(weight_matrix)

column_sums = weight_matrix.sum(axis=0)
freq_matrix = weight_matrix / column_sums

print("\nRelative Frequency Matrix")
print(freq_matrix.round(3))

null_prob = 0.25
log_likelihood_matrix = np.log(freq_matrix / null_prob)

print("\nLog-Likelihood Matrix")
print(log_likelihood_matrix.round(3))

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

def score_window(window, ll_matrix):
    score = 0.0
    for i, base in enumerate(window):
        score += ll_matrix.loc[base, i]
    return score

results = []

for i in range(len(S) - motif_length + 1):
    window = S[i:i + motif_length]
    score = score_window(window, log_likelihood_matrix)
    results.append((i, window, score))

print("\nSliding Window Scores")
for pos, window, score in results:
    print(f"Position {pos:2d} | {window} | Score = {score:.3f}")

scores = [r[2] for r in results]
mean_score = np.mean(scores)
std_score = np.std(scores)

threshold = mean_score + std_score

print("\nSignal Detection")
print(f"Mean score     : {mean_score:.3f}")
print(f"Std deviation  : {std_score:.3f}")
print(f"Threshold      : {threshold:.3f}")

print("\nPutative exon–intron borders:")
found = False
for pos, window, score in results:
    if score > threshold:
        print(f"  -> Position {pos:2d} | {window} | Score = {score:.3f}")
        found = True

if not found:
    print("  No strong exon–intron border detected.")

print("\n=== Conclusion ===")
if found:
    print(
        "High log-likelihood scores indicate that the sequence S\n"
        "contains regions consistent with an exon–intron boundary.\n"
        "These regions match the learned splice-site motif profile."
    )
else:
    print(
        "No significant log-likelihood peaks were detected.\n"
        "The sequence S likely does not contain an exon–intron border."
    )
