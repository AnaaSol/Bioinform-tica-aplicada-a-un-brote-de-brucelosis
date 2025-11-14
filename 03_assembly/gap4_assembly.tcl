# Gap4 TCL script for assembling Sanger sequences
# Based on Practica7_2024.pdf methodology

# Open/create database
set db_name "brucella_sanger"
set db_version "0"

# Create new database
if {[catch {gap_open $db_name -create -version $db_version} io]} {
    puts "Error creating database: $io"
    exit 1
}

# Assembly parameters (from practice defaults)
set assemble_params {
    -max_pmismatch 30.0
    -min_match 20
    -min_overlap 20
    -word_length 4
    -max_dashes 10
    -min_conf 8.0
    -use_avg_insert 1
    -max_read_length 20000
    -align_mode 1
}

# Read file list from pregap.passed
set passed_file "/home/ana/Desktop/Brucelosis_TP_Integrador/02_pregap4_processing/pregap.passed"
if {[catch {open $passed_file r} fh]} {
    puts "Error opening $passed_file: $fh"
    exit 1
}

set file_list [split [read $fh] "\n"]
close $fh

# Filter empty lines
set file_list [lsearch -all -inline -not -exact $file_list ""]

puts "Found [llength $file_list] files to assemble"

# Normal shotgun assembly
puts "Starting assembly..."
if {[catch {
    eval directed_assembly \
        -io $io \
        -files $file_list \
        {*}$assemble_params \
        -align_which 2 \
        -enter_all 1 \
        -permit_joins 1
} result]} {
    puts "Assembly error: $result"
    gap_close -io $io
    exit 1
}

puts "Assembly completed"
puts "Result: $result"

# Get assembly statistics
set num_contigs [db_info num_contigs -io $io]
puts "Number of contigs: $num_contigs"

# Close database
gap_close -io $io

puts "Assembly finished successfully"
exit 0
