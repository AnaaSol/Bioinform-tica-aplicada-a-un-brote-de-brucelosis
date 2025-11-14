  #!/bin/bash

  echo "=== ENSAMBLADO HÍBRIDO - ESTADÍSTICAS ==="
  echo ""

  echo "Scaffolds:"
  grep -c "^>" scaffolds.fasta
  echo ""

  echo "Contig más largo:"
  grep "^>" scaffolds.fasta | sed 's/.*length_//' | sed 's/_cov.*//' | sort -rn | head -1
  echo ""

  echo "Cobertura promedio:"
  grep "^>" scaffolds.fasta | sed 's/.*cov_//' | awk '{sum+=$1; count++} END {print sum/count}'
  echo ""

  echo "Longitud total:"
  grep "^>" scaffolds.fasta | sed 's/.*length_//' | sed 's/_cov.*//' | awk '{sum+=$1} END {print sum " pb"}'
  echo ""

  echo "Detalle de scaffolds:"
  grep "^>" scaffolds.fasta | sed 's/>//' | column -t
