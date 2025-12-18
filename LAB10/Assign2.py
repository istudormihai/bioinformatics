# Download 10 influenza genomes. Adapt your application from the previous assignment in order to scan each genome for possible motifs.
# For each genome make a chart that shows the signal with the most likely locations of real functional motifs.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

motifs = [
    "GTCATTACTA",
    "ACACAATAGA",
    "GCGAGGGGTG",
    "GGGGGGGGGG",
    "TTTTTTTTTT",
    "AATCCAAAGA",
    "AAGAACATAA",
    "AGGGTTCAGG",
    "CTATTGTCTT"
]

nucleotides = ['A', 'C', 'G', 'T']
ntIndex = {nt: i for i, nt in enumerate(nucleotides)}

motifLength = len(motifs[0])

countMatrix = np.zeros((4, motifLength))
for motif in motifs:
    for position, nucleotide in enumerate(motif):
        countMatrix[ntIndex[nucleotide], position] += 1

pseudoCount = 1
countMatrixPc = countMatrix + pseudoCount
frequencyMatrix = countMatrixPc / countMatrixPc.sum(axis=0)
nullModel = 0.25
logLikelihoodMatrix = np.log(frequencyMatrix / nullModel)

def readFasta(filePath):
    seq = ""
    with open(filePath, 'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq += line.strip().upper()
    return seq

def scanSequence(sequence, pwmMatrix, motifLength):
    scores = []
    for i in range(len(sequence) - motifLength + 1):
        window = sequence[i:i + motifLength]
        score = 0
        for pos, nt in enumerate(window):
            if nt in ntIndex:
                score += pwmMatrix[ntIndex[nt], pos]
        scores.append(score)
    return scores

def getTopPeaks(scores, topN=5):
    scoresArr = np.array(scores)

    topIndices = scoresArr.argsort()[-topN:][::-1]
    return topIndices, scoresArr[topIndices]

fig, axes = plt.subplots(5, 2, figsize=(16, 22))
axes = axes.flatten()

for i in range(1, 11):
    fastaFile = f"Influenza{i}.fasta"
    if not os.path.exists(fastaFile):
        print(f"File {fastaFile} not found, skipping...")
        continue

    genomeSeq = readFasta(fastaFile)
    scores = np.array(scanSequence(genomeSeq, logLikelihoodMatrix, motifLength))
    topIndices, topScores = getTopPeaks(scores, topN=5)

    ax = axes[i - 1]

    ax.plot(scores, color="#1f77b4", linewidth=1.2, label="PWM score")
    ax.fill_between(
        range(len(scores)),
        scores,
        color="#1f77b4",
        alpha=0.25
    )

    ax.scatter(
        topIndices,
        topScores,
        color="#d62728",
        marker="^",
        s=80,
        zorder=3,
        label="Top motif hits"
    )

    threshold = scores.mean() + scores.std()
    ax.axhline(
        threshold,
        color="black",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
        label="Mean + 1σ"
    )

    for idx, score in zip(topIndices, topScores):
        ax.annotate(
            str(idx),
            (idx, score),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center",
            fontsize=8
        )

    ax.set_title(f"{fastaFile}", fontsize=11, fontweight="bold")
    ax.set_xlabel("Genome position")
    ax.set_ylabel("Log-likelihood score")

    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(fontsize=8, loc="upper right")

plt.tight_layout()
plt.show()

