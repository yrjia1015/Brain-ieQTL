#!/usr/bin/env perl
use strict;
use warnings;

sub usage {
    return <<"USAGE";
Usage:
  perl gtf_to_gene_bed.pl <annotation.gtf> [feature_label] > genes.bed

Convert gene features from a GTF/GFF-like annotation file to BED6.

Output columns:
  chrom  start0  end  gene_id  feature_label  strand

The gene identifier is selected from the first available attribute in:
  gene_id, ID, Name, gene_name

Example:
  perl gtf_to_gene_bed.pl gencode.v37lift37.genes.gtf gene > genes.hg19.bed
  sed -i -e 's/^chr//' genes.hg19.bed
USAGE
}

my ($gtf_file, $feature_label) = @ARGV;
die usage() if !defined $gtf_file;
$feature_label = "gene" if !defined $feature_label || $feature_label eq "";
die "ERROR: Annotation file not found: $gtf_file\n" if !-e $gtf_file;

gtf_to_bed($gtf_file, $feature_label);

sub gtf_to_bed {
    my ($filename, $feature_label) = @_;

    open my $gtf_fh, '<', $filename or die "ERROR: Cannot open $filename: $!\n";
    while (my $line = <$gtf_fh>) {
        chomp $line;
        next if $line =~ /\A#/ || $line =~ /\A\s*\z/;

        my @fields = split /\t/, $line;
        next if @fields < 9;
        my ($chrom, $type, $start, $end, $strand, $attributes) = @fields[0, 2, 3, 4, 6, 8];
        next if $type ne "gene";

        my %attribute = parse_attributes($attributes);
        my $gene_id = first_existing_attribute(\%attribute, qw(gene_id ID Name gene_name));
        next if !defined $gene_id || $gene_id eq "";

        print join("\t", $chrom, $start - 1, $end, $gene_id, $feature_label, $strand), "\n";
    }
    close $gtf_fh;
}

sub first_existing_attribute {
    my ($attribute, @keys) = @_;
    for my $key (@keys) {
        return $attribute->{$key} if exists $attribute->{$key};
    }
    return;
}

sub parse_attributes {
    my ($attribute_string) = @_;
    my %attribute;

    for my $field (split /;\s*/, $attribute_string) {
        next if $field eq "";
        my ($key, $value);
        if ($field =~ /\A([^=\s]+)=("?)(.*?)\2\z/) {
            ($key, $value) = ($1, $3);
        } elsif ($field =~ /\A(\S+)\s+"?([^"]+)"?\z/) {
            ($key, $value) = ($1, $2);
        } else {
            next;
        }
        $attribute{$key} = $value;
    }

    return %attribute;
}
