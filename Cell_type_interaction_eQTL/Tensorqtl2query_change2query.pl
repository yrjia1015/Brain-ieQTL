use strict;
use warnings;

# 1st file: gene reference
# 2nd file: SNP reference
# 3rd file: TensorQTL results
# 4th argument: output file path
my ($fileA, $fileB, $fileC, $fileD) = @ARGV;

# Build hash table AA
my %AA;
open(my $fhA, '<', $fileA) or die "Cannot open file A: $!";
while (my $line = <$fhA>) {
    chomp $line;
    my @elements = split(/\t/, $line);
    my $key = $elements[3];
    my $value = join("\t", @elements);
    $AA{$key} = $value;
}
close $fhA;

# Build hash table BB
my %BB;
open(my $fhB, '<', $fileB) or die "Cannot open file B: $!";
while (my $line = <$fhB>) {
    chomp $line;
    my @elements = split(/\t/, $line);
    my $key = $elements[1];
    my $value = join("\t", @elements);
    $BB{$key} = $value;
}
close $fhB;

open(my $fhC, '<', $fileC) or die "Cannot open file C: $!";
open(my $fhD, '>', $fileD) or die "Cannot open file D: $!";
while (my $line = <$fhC>) {
    chomp $line;
    if ($. == 1) {
        my $new_line = "$line\tProbe_Chr\tProbe_bp\tProbe_bp_no\tprobe\tGene\tOrientation\tChr\trsid\tno\tBP\tA1\tA2";
        print $fhD "$new_line\n";
    } else {
        my @elements = split(/\t/, $line);
        if (exists $AA{$elements[0]}) {
            $line .= "\t$AA{$elements[0]}";
        }
        if (exists $BB{$elements[1]}) {
            $line .= "\t$BB{$elements[1]}";
        }
        my @new_elements = split /\t/, $line;
        my $new_element_27 = $new_elements[26];
        my $new_element_28 = $new_elements[27];
        if ($new_element_27 !~ /^[A-Z]+$/ || $new_element_28 !~ /^[A-Z]+$/) {
            next;
        }
        print $fhD "$line\n";
    }
}
close $fhC;
close $fhD;