#!/usr/bin/env perl
use strict;
use warnings;

sub usage {
    return <<"USAGE";
Usage:
  perl expression_matrix_to_phenotype_bed.pl <gene_bed> <expression_matrix> > phenotype.bed

Create a tensorQTL-compatible phenotype BED from a gene BED file and an
expression matrix.

Inputs:
  gene_bed           BED6-like file: chrom, start, end, gene_id, score, strand
  expression_matrix  Tab-delimited matrix with gene ID in column 1 and sample
                     IDs in the header.

Output:
  Header plus BED-like rows:
    #chr  start  end  phenotype_id  sample1  sample2 ...

Notes:
  - Only autosomes 1-22 are retained.
  - chr prefixes should be removed before running if genotypes use numeric
    chromosome names.
  - The phenotype position is the midpoint of the input BED interval, matching
    the historical expression eQTL pipeline behavior.
USAGE
}

my ($bed_file, $expression_file) = @ARGV;
die usage() if !defined $bed_file || !defined $expression_file;
die "ERROR: BED file not found: $bed_file\n"              if !-e $bed_file;
die "ERROR: Expression matrix not found: $expression_file\n" if !-e $expression_file;

my %autosomes = map { $_ => 1 } 1..22;
my %gene_to_position = load_bed_to_hash($bed_file, \%autosomes);

open my $expr_fh, '<', $expression_file or die "ERROR: Cannot open $expression_file: $!\n";
my $header = <$expr_fh>;
die "ERROR: Empty expression matrix: $expression_file\n" if !defined $header;
chomp $header;
my ($id_header, @sample_ids) = split /\t/, $header;
die "ERROR: Expression matrix must contain at least one sample column\n" if @sample_ids < 1;
print join("\t", "#chr", "start", "end", "phenotype_id", @sample_ids), "\n";

while (my $line = <$expr_fh>) {
    chomp $line;
    next if $line =~ /\A\s*\z/;

    my ($gene_id, @values) = split /\t/, $line;
    next if !defined $gene_id || !exists $gene_to_position{$gene_id};
    print join("\t", $gene_to_position{$gene_id}, $gene_id, @values), "\n";
}
close $expr_fh;

sub load_bed_to_hash {
    my ($filename, $autosomes) = @_;
    my %gene_to_position;

    open my $bed_fh, '<', $filename or die "ERROR: Cannot open $filename: $!\n";
    while (my $line = <$bed_fh>) {
        chomp $line;
        next if $line =~ /\A\s*\z/;
        my ($chrom, $start, $end, $gene_id) = split /\t/, $line;
        next if !defined $gene_id;

        $chrom =~ s/\Achr//;
        next if !exists $autosomes->{$chrom};

        my $midpoint_1based = int(($start + $end) / 2 + 1);
        $gene_to_position{$gene_id} = join "\t", $chrom, $midpoint_1based - 1, $midpoint_1based;
    }
    close $bed_fh;

    return %gene_to_position;
}
