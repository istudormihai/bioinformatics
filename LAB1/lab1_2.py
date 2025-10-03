sequence = "ATTTCGCCGATA"

s = set(sequence)

frequencies = {}

for ch in sequence: 
    frequencies[ch] = frequencies.get(ch, 0) + 1

for el in sorted(s):
    percentage = percentage = (frequencies[el] / len(sequence)) * 100
    print(f"{el} - {percentage:.2f}%")