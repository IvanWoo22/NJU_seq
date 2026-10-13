#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use File::Basename;
use File::Temp qw(tempfile);
@ARGV == 2 or die "Usage: signature_count.pl sites.tsv output.pdf\n";

sub BASE_CHANGE {
    my $BASE = shift;
    if ($BASE eq "T") {
        return ("U");
    }
    else {
        return ($BASE);
    }
}

sub METHYL {
    my $BASE = shift;
    return ($BASE . "m");
}

my $filedir = dirname(__FILE__);
my $envdir = $ENV{'PWD'};

open(my $IN, "<", $ARGV[0]);
my ($OUT, $temp_path) = tempfile(UNLINK => 1);

my %sig;
while (<$IN>) {
    chomp;
    next if /^#/ || /^(?:Chromosome|Position)\t/;
    my (undef, undef, undef, $lef, $nm, $rgt, undef) = split;
    ($lef, $nm, $rgt) = map { my $b = uc($_ // ""); $b =~ tr/U/T/; $b } ($lef, $nm, $rgt);
    next unless "$lef$nm$rgt" =~ /^[ACGT]{3}$/;
    if (exists($sig{ $lef . $nm . $rgt })) {
        $sig{ $lef . $nm . $rgt }++;
    }
    else {
        $sig{ $lef . $nm . $rgt } = 1;
    }
}

my $data_table;
foreach my $nm (qw(A G C T)) {
    foreach my $lef (qw(A G C T)) {
        foreach my $rgt (qw(A G C T)) {
            if (exists($sig{ $lef . $nm . $rgt })) {
                $data_table .=
                    BASE_CHANGE($lef)
                        . METHYL(BASE_CHANGE($nm))
                        . BASE_CHANGE($rgt)
                        . "\t$sig{$lef . $nm . $rgt}\n";
            }
            else {
                $data_table .=
                    BASE_CHANGE($lef)
                        . METHYL(BASE_CHANGE($nm))
                        . BASE_CHANGE($rgt) . "\t0\n";
            }
        }
    }
}

print $OUT ("$data_table");

close($OUT);
close($IN);
system('Rscript', "$filedir/signature_count.R", $temp_path, $ARGV[1]) == 0
    or die "Signature plot failed (R exit status $?)\n";
