# Bioinformatica aplicada a un brote de brucelosis

Análisis de un brote de Brucella suis en Argentina (Diciembre 2024). El proyecto involucra ensamblado genómico, anotación y análisis comparativo usando datos de secuenciación Sanger e Illumina para caracterizar el patógeno y desarrollar herramientas diagnósticas.

### ESCENARIO

  - Patógeno: Brucella suis (bacteria Gram-negativa causante de brucelosis/fiebre de Malta)
  - Brote: Afectando ganado y humanos en Argentina
  - Datos disponibles:
    - Secuencias Sanger (150-300 lecturas) del aislado bacteriano
    - Lecturas paired-end Illumina (55,000-70,000 lecturas)
  - Objetivo: Análisis genómico completo para comprensión epidemiológica y desarrollo de diagnóstico PCR

### ESTRUCTURA DEL PROYECTO

Investigar y documentar:
  - Características del genoma (tamaño, %GC, estructura cromosoma/plásmidos)
  - Mecanismos de patogénesis y factores de virulencia
  - Impacto en salud pública y consecuencias económicas

Secuencias Sanger:
  - Procesar electroferogramas (formato ab1)
  - Eliminar regiones de baja calidad y contaminación de vector
  - Herramientas: Phred, Lucy2

Lecturas Illumina:
  - Solo evaluación de calidad (NO procesar)
  - Herramientas: FastQC

Ensamblado Sanger:
  - Generar contig(s) y reportar longitud/cobertura
  - Herramienta: CAP3

Ensamblado híbrido:
  - Combinar contigs Sanger + lecturas Illumina
  - Herramienta: SPAdes (modo --trusted-contigs)
  - Calcular tamaño del genoma y %GC
  - Comparar con B. suis de referencia reportado

Predicción de genes y anotación:
  - Predicción de genes ab initio
  - Pipeline de anotación integrado (ej: Prokka, RAST)
  - Visualización en navegador de genomas (Artemis/IGV)
  - Refinamiento de anotación funcional vía BLAST + dominios conservados
  - Generar tabla gen-función

Análisis comparativo: Seleccionar una proteína y realizar:
  - PSI-BLAST: 3 iteraciones, top 10 hits de otras especies
  - Alineamiento Múltiple: Clustal X (tu proteína + 10 hits)
  - Análisis de Dominios: InterPro para dominios conservados y términos GO
  - Logo de Secuencia: WebLogo para regiones conservadas
  - Árbol Filogenético: Construir y analizar relaciones evolutivas de B. suis

Final: Proponer gen óptimo para kit diagnóstico PCR basado en unicidad y características de amplificación

### HERRAMIENTAS REQUERIDAS

  - Calidad: Phred, Lucy2, FastQC
  - Ensamblado: CAP3, SPAdes
  - Anotación: Prokka/RAST, BLAST, InterPro
  - Visualización: Artemis, IGV
  - Filogenética: PSI-BLAST, Clustal X, software de filogenia
  - Análisis de secuencias: WebLogo
  - Programación: Perl con módulos BioPerl

