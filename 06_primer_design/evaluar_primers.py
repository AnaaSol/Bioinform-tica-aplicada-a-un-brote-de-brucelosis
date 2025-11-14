#!/usr/bin/env python3
import csv
import re

def gc_content(seq):
    seq = seq.upper()
    g = seq.count("G")
    c = seq.count("C")
    return 100 * (g + c) / len(seq)

def has_bad_repeats(seq):
    """Detecta repeticiones largas que afectan PCR (AAAAA, CCCCC, etc)"""
    return bool(re.search(r"(A{5,}|T{5,}|G{5,}|C{5,})", seq))

def basic_dimer_check(seq1, seq2):
    """Chequeo muy simple de dimerización (complementos de 4 nt consecutivos)"""
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    complement = str.maketrans("ATGC", "TACG")
    seq2_rc = seq2.translate(complement)[::-1]

    # Ventana deslizante de 4 nt
    for i in range(len(seq1) - 3):
        kmer = seq1[i:i+4]
        if kmer in seq2_rc:
            return True
    return False

def score_primers(row):
    score = 100  # arrancamos desde 100, restamos penalizaciones

    fw = row["Primer_Fw"]
    rv = row["Primer_Rv"]
    tm_fw = float(row["Tm_Fw"])
    tm_rv = float(row["Tm_Rv"])
    size = int(row["Tamaño_Producto"])

    # === ΔTm ===
    delta_tm = abs(tm_fw - tm_rv)
    if delta_tm > 5:
        score -= 25
    elif delta_tm > 3:
        score -= 10

    # === GC% ideal entre 40–60% ===
    gc_fw = gc_content(fw)
    gc_rv = gc_content(rv)

    for gc in [gc_fw, gc_rv]:
        if gc < 35 or gc > 70:
            score -= 15
        elif gc < 40 or gc > 60:
            score -= 5

    # === Runs de bases ===
    if has_bad_repeats(fw):
        score -= 10
    if has_bad_repeats(rv):
        score -= 10

    # === Dímeros básicos ===
    if basic_dimer_check(fw, rv):
        score -= 15

    # === Tamaño del producto ===
    if size < 120 or size > 300:
        score -= 20  # penalización fuerte (ideal PCR diagnóstica 120–300 bp)

    if size > 1000:
        score -= 40  # amplicón demasiado grande

    return score

# ============================
# MAIN SCRIPT
# ============================

with open("res_primers.tab", "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    results = []

    for row in reader:
        s = score_primers(row)
        row["Score"] = s
        results.append(row)

# Ordenar de mejor a peor
results_sorted = sorted(results, key=lambda x: x["Score"], reverse=True)

# Mostrar ranking
print("\n=== RANKING DE PRIMERS ===")
for r in results_sorted:
    print(f"{r['Gen']}: Score={r['Score']}, Amplicón={r['Tamaño_Producto']} bp, ΔTm={abs(float(r['Tm_Fw']) - float(r['Tm_Rv'])):.1f}")

best = results_sorted[0]

print("\n=== MEJOR GEN PARA KIT DIAGNÓSTICO ===")
print(f"Gen: {best['Gen']}")
print(f"Score: {best['Score']}")
print(f"Amplicón: {best['Tamaño_Producto']} bp")
print(f"Tm Fw: {best['Tm_Fw']}°C")
print(f"Tm Rv: {best['Tm_Rv']}°C")
print(f"Secuencia Fw: {best['Primer_Fw']}")
print(f"Secuencia Rv: {best['Primer_Rv']}")

