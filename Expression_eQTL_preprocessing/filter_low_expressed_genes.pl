#!/usr/bin/env perl
use strict;
use warnings;

sub usage {
    return <<"USAGE";
Usage:
  perl filter_low_expressed_genes.pl <TPM_matrix> <COUNT_matrix> [min_TPM] [min_TPM_fraction] [min_COUNT] [min_COUNT_fraction] [output_mode]

Filter genes that pass both TPM and raw-count expression thresholds.

Input matrix format:
  Tab-delimited text with a header line. Column 1 is gene ID/name; remaining
  columns are samples. TPM and COUNT matrices must use the same gene IDs.

Arguments:
  min_TPM             Default: 0.5
  min_TPM_fraction    Default: 0.2
  min_COUNT           Default: 6
  min_COUNT_fraction  Default: 0.2
  output_mode         1 = output filtered TPM matrix, 2 = output filtered COUNT matrix. Default: 1

Example:
  perl filter_low_expressed_genes.pl genes.TPM.txt genes.fragments.txt 0.1 0.2 6 0.2 2 > genes.filtered.count.txt
USAGE
}

my (
    $matrix_exp_TPM,
    $matrix_exp_COUNT,
    $min_value_TPM,
    $min_percent_TPM,
    $min_value_COUNT,
    $min_percent_COUNT,
    $output_mode
) = @ARGV;

die usage() if !defined $matrix_exp_TPM || !defined $matrix_exp_COUNT;
$min_value_TPM    = 0.5 if !defined $min_value_TPM;
$min_percent_TPM  = 0.2 if !defined $min_percent_TPM;
$min_value_COUNT  = 6   if !defined $min_value_COUNT;
$min_percent_COUNT = 0.2 if !defined $min_percent_COUNT;
$output_mode      = 1   if !defined $output_mode;

die "ERROR: TPM matrix not found: $matrix_exp_TPM\n"       if !-e $matrix_exp_TPM;
die "ERROR: COUNT matrix not found: $matrix_exp_COUNT\n"   if !-e $matrix_exp_COUNT;
die "ERROR: output_mode must be 1 (TPM) or 2 (COUNT)\n"    if $output_mode !~ /\A[12]\z/;

my @passed_tpm   = load_passing_genes($matrix_exp_TPM,   $min_value_TPM,   $min_percent_TPM);
my @passed_count = load_passing_genes($matrix_exp_COUNT, $min_value_COUNT, $min_percent_COUNT);
my %keep = map { $_ => 1 } intersect_lists(\@passed_tpm, \@passed_count);

my $output_matrix = $output_mode == 1 ? $matrix_exp_TPM : $matrix_exp_COUNT;
print_filtered_matrix($output_matrix, \%keep);

sub load_passing_genes {
    my ($filename, $min_value, $min_fraction) = @_;
    my @genes;

    open my $fh, '<', $filename or die "ERROR: Cannot open $filename: $!\n";
    my $header = <$fh>;
    die "ERROR: Empty matrix: $filename\n" if !defined $header;
    chomp $header;
    my @header_fields = split /\t/, $header;
    my $num_samples = @header_fields - 1;
    die "ERROR: Matrix must contain at least one sample column: $filename\n" if $num_samples < 1;

    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /\A\s*\z/;
        my ($gene_name, @values) = split /\t/, $line;
        if (passes_expression_threshold(\@values, $min_value, $num_samples, $min_fraction)) {
            push @genes, $gene_name;
        }
    }
    close $fh;

    return @genes;
}

sub passes_expression_threshold {
    my ($values, $min_value, $num_samples, $min_fraction) = @_;
    my $pass_count = 0;

    for my $value (@$values) {
        $pass_count++ if defined $value && $value ne '' && $value >= $min_value;
        return 1 if $pass_count / $num_samples >= $min_fraction;
    }

    return 0;
}

sub intersect_lists {
    my ($list_a, $list_b) = @_;
    my %seen_a = map { $_ => 1 } @$list_a;
    return grep { exists $seen_a{$_} } @$list_b;
}

sub print_filtered_matrix {
    my ($filename, $keep) = @_;

    open my $fh, '<', $filename or die "ERROR: Cannot open $filename: $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        if ($. == 1) {
            print "$line\n";
            next;
        }
        my ($gene_name) = split /\t/, $line, 2;
        print "$line\n" if exists $keep->{$gene_name};
    }
    close $fh;
}
