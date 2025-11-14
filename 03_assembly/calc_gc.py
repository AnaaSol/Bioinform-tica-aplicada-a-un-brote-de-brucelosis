from Bio import SeqIO
from Bio.SeqUtils import GC

# Leer FASTA
record = list(SeqIO.parse("scaffold_final.fasta", "fasta"))[0]

# Estadísticas
size = len(record.seq)
gc_content = GC(record.seq)

print("=== ESTADÍSTICAS DEL GENOMA ===")
print(f"Tamaño: {size:,} pb")
print(f"%G+C: {gc_content:.2f}%")
print()
print("=== COMPARACIÓN ===")
print("B. suis (genoma completo): 3.3 Mb, 57.2% GC")
print(f"Este fragmento: {size/1000:.1f} kb, {gc_content:.2f}% GC")

