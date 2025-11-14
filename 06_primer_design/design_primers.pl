#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::GFF;

# Script para diseño de primers PCR para genes predichos
# Autor: TP Integrador Brucella suis
# Fecha: 2025-11-13

# Verificar argumentos
if (@ARGV != 2) {
    die "Uso: $0 <archivo_fasta> <archivo_gff>\n";
}

my ($fasta_file, $gff_file) = @ARGV;

# Verificar que los archivos existan
die "Error: No se encuentra $fasta_file\n" unless -e $fasta_file;
die "Error: No se encuentra $gff_file\n" unless -e $gff_file;

print STDERR "Leyendo genoma desde $fasta_file...\n";
print STDERR "Leyendo genes desde $gff_file...\n";

# Leer secuencia genómica
my $seqio = Bio::SeqIO->new(-file => $fasta_file, -format => 'fasta');
my %genome_seqs;
while (my $seq = $seqio->next_seq) {
    $genome_seqs{$seq->id} = $seq;
}

# Abrir archivos de salida
open(my $primers_fh, '>', 'res_primers.tab') or die "No se puede crear res_primers.tab: $!\n";
open(my $products_fh, '>', 'res_prod.fa') or die "No se puede crear res_prod.fa: $!\n";

# Escribir encabezado del archivo de primers
print $primers_fh "Gen\tPrimer_Fw\tTm_Fw\tTa_Fw\tPrimer_Rv\tTm_Rv\tTa_Rv\tTamaño_Producto\n";

# Leer archivo GFF
my $gffio = Bio::Tools::GFF->new(-file => $gff_file, -gff_version => 3);

my $gene_count = 0;

while (my $feature = $gffio->next_feature()) {
    # Procesar solo features de tipo CDS
    next unless $feature->primary_tag eq 'CDS';

    # Obtener información del gen
    my $seq_id = $feature->seq_id;
    my $start = $feature->start;
    my $end = $feature->end;
    my $strand = $feature->strand;  # +1 o -1

    # Obtener nombre del gen desde los atributos
    my $gene_name;
    if ($feature->has_tag('locus_tag')) {
        ($gene_name) = $feature->get_tag_values('locus_tag');
    } elsif ($feature->has_tag('ID')) {
        ($gene_name) = $feature->get_tag_values('ID');
    } else {
        $gene_name = "gene_${gene_count}";
    }

    # Verificar que tengamos la secuencia del contig
    unless (exists $genome_seqs{$seq_id}) {
        print STDERR "Advertencia: No se encuentra secuencia para $seq_id (gen $gene_name)\n";
        next;
    }

    my $genome_seq = $genome_seqs{$seq_id};

    # Diseñar primers según la hebra
    my ($primer_fw, $primer_rv, $product_seq);

    if ($strand == 1) {
        # Gen en hebra positiva (+)
        # Primer Fw: 5 nt antes del inicio del gen, 20 nt total
        my $fw_start = $start - 5;
        $fw_start = 1 if $fw_start < 1;  # No salir del contig
        my $fw_end = $fw_start + 19;

        # Primer Rv: 5 nt después del fin del gen, 20 nt total, complemento reverso
        my $rv_end = $end + 5;
        $rv_end = $genome_seq->length if $rv_end > $genome_seq->length;
        my $rv_start = $rv_end - 19;

        # Extraer secuencias
        $primer_fw = $genome_seq->subseq($fw_start, $fw_end);
        my $rv_seq_forward = $genome_seq->subseq($rv_start, $rv_end);
        $primer_rv = reverse_complement($rv_seq_forward);

        # Secuencia del producto PCR
        $product_seq = $genome_seq->subseq($fw_start, $rv_end);

    } else {
        # Gen en hebra negativa (-)
        # Primer Fw: 5 nt después del fin del gen (en coord genómicas), complemento reverso
        my $fw_end = $end + 5;
        $fw_end = $genome_seq->length if $fw_end > $genome_seq->length;
        my $fw_start = $fw_end - 19;

        # Primer Rv: 5 nt antes del inicio del gen, 20 nt total
        my $rv_start = $start - 5;
        $rv_start = 1 if $rv_start < 1;
        my $rv_end = $rv_start + 19;

        # Extraer secuencias
        my $fw_seq_forward = $genome_seq->subseq($fw_start, $fw_end);
        $primer_fw = reverse_complement($fw_seq_forward);
        $primer_rv = $genome_seq->subseq($rv_start, $rv_end);

        # Secuencia del producto PCR (hebra negativa, complemento reverso)
        my $prod_forward = $genome_seq->subseq($rv_start, $fw_end);
        $product_seq = reverse_complement($prod_forward);
    }

    # Calcular Tm (Temperatura de Melting)
    my $tm_fw = calculate_tm($primer_fw);
    my $tm_rv = calculate_tm($primer_rv);

    # Calcular Ta (Temperatura de Annealing)
    my $ta_fw = $tm_fw - 5;
    my $ta_rv = $tm_rv - 5;

    # Calcular tamaño del producto
    my $product_size = length($product_seq);

    # Escribir resultados en archivo de primers
    printf $primers_fh "%s\t%s\t%.1f\t%.1f\t%s\t%.1f\t%.1f\t%d\n",
        $gene_name, $primer_fw, $tm_fw, $ta_fw, $primer_rv, $tm_rv, $ta_rv, $product_size;

    # Escribir producto PCR en archivo multifasta
    print $products_fh ">$gene_name size=$product_size\n";
    print $products_fh format_fasta_sequence($product_seq), "\n";

    $gene_count++;
}

close($primers_fh);
close($products_fh);

print STDERR "\nProcesamiento completado:\n";
print STDERR "  Genes procesados: $gene_count\n";
print STDERR "  Archivo de primers: res_primers.tab\n";
print STDERR "  Archivo de productos: res_prod.fa\n";

#############################
# Subrutinas
#############################

# Calcular Tm usando fórmula: Tm = 4(G+C) + 2(A+T)
sub calculate_tm {
    my ($seq) = @_;
    $seq = uc($seq);

    my $a = ($seq =~ tr/A//);
    my $t = ($seq =~ tr/T//);
    my $g = ($seq =~ tr/G//);
    my $c = ($seq =~ tr/C//);

    my $tm = 4 * ($g + $c) + 2 * ($a + $t);

    return $tm;
}

# Calcular complemento reverso
sub reverse_complement {
    my ($seq) = @_;
    $seq = uc($seq);
    $seq = reverse($seq);
    $seq =~ tr/ATGC/TACG/;
    return $seq;
}

# Formatear secuencia FASTA a 60 caracteres por línea
sub format_fasta_sequence {
    my ($seq) = @_;
    my $formatted = '';
    my $length = length($seq);

    for (my $i = 0; $i < $length; $i += 60) {
        my $chunk = substr($seq, $i, 60);
        $formatted .= $chunk . "\n";
    }

    return $formatted;
}

__END__

=head1 NOMBRE

design_primers.pl - Diseño de primers PCR para genes predichos de Brucella suis

=head1 SINOPSIS

  perl design_primers.pl <archivo_fasta> <archivo_gff>

=head1 DESCRIPCIÓN

Script para diseño automático de primers PCR para todos los genes predichos
en un ensamblado genómico. Genera primers forward y reverse de 20 nt que
amplifican cada gen completo, comenzando 5 nt antes y terminando 5 nt después
de las coordenadas del gen.

=head1 ARGUMENTOS

  archivo_fasta    Archivo FASTA con la secuencia genómica ensamblada
  archivo_gff      Archivo GFF3 con los genes predichos

=head1 SALIDA

  res_primers.tab  Tabla con información de primers (TSV)
  res_prod.fa      Secuencias de productos PCR (multiFASTA)

=head1 AUTOR

TP Integrador - Brucella suis

=cut
