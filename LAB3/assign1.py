import math
def formula1(dna):
    dna = dna.upper()
    g = dna.count('G')
    c = dna.count('C')
    a = dna.count('A')
    t = dna.count('T')
    tm = 4 * (g + c) + 2 * (a + t)
    return tm

def formula2(dna, na_conc=0.05):
    dna = dna.upper()
    length = len(dna)
    gc_percent = (dna.count('G') + dna.count('C')) / length * 100
    tm = 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - (600 / length)
    return tm

if __name__ == "__main__":
    dna = input("Enter DNA sequence: ").strip()
    tm1 = formula1(dna)
    tm2 = formula2(dna)
    print(f"\nFirst formula: {tm1:.2f} °C")
    print(f"Second formula: {tm2:.2f} °C")