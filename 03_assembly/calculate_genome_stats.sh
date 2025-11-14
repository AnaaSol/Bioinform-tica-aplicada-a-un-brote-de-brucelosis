#!/bin/bash

  FASTA="scaffold_final.fasta"

  echo "=== ESTADÍSTICAS DEL GENOMA ENSAMBLADO ==="
  echo ""

  # Tamaño
  SIZE=$(grep -v "^>" $FASTA | tr -d '\n' | wc -c)
  echo "Tamaño del genoma: $SIZE pb"

  # Obtener secuencia sin saltos de línea
  SEQ=$(grep -v "^>" $FASTA | tr -d '\n')

  # Contar bases
  A=$(echo $SEQ | grep -o "A" | wc -l)
  T=$(echo $SEQ | grep -o "T" | wc -l)
  G=$(echo $SEQ | grep -o "G" | wc -l)
  C=$(echo $SEQ | grep -o "C" | wc -l)

  echo ""
  echo "Composición de bases:"
  echo "  A: $A ($(awk "BEGIN {printf \"%.2f\", ($A/$SIZE)*100}")%)"
  echo "  T: $T ($(awk "BEGIN {printf \"%.2f\", ($T/$SIZE)*100}")%)"
  echo "  G: $G ($(awk "BEGIN {printf \"%.2f\", ($G/$SIZE)*100}")%)"
  echo "  C: $C ($(awk "BEGIN {printf \"%.2f\", ($C/$SIZE)*100}")%)"

  # Calcular %GC
  GC=$((G + C))
  GC_PERCENT=$(awk "BEGIN {printf \"%.2f\", ($GC/$SIZE)*100}")

  echo ""
  echo "CONTENIDO G+C: $GC_PERCENT%"
  echo ""

  # Comparación con B. suis reportado
  echo "=== COMPARACIÓN CON Brucella suis ==="
  echo "Genoma completo B. suis (referencia NCBI):"
  echo "  Tamaño: ~3.3 Mb (3,315,175 pb)"
  echo "  %G+C: ~57.2%"
  echo ""
  echo "Fragmento ensamblado en este TP:"
  echo "  Tamaño: $SIZE pb (~$(awk "BEGIN {printf \"%.2f\", ($SIZE/3315175)*100}")% del genoma)"
  echo "  %G+C: $GC_PERCENT%"

  # Análisis
  DIFF=$(awk "BEGIN {printf \"%.2f\", $GC_PERCENT - 57.2}")
  echo ""
  echo "Diferencia en %G+C: $DIFF%"

  if (( $(echo "$DIFF < 2 && $DIFF > -2" | bc -l) )); then
      echo "✓ El %G+C es consistente con B. suis"
  else
      echo "⚠ El %G+C difiere del reportado para B. suis"
  fi

