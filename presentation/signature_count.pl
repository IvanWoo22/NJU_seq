#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use File::Basename;

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
open(my $OUT, ">", $envdir . "/" . $ARGV[0] . ".tmp");

my %sig;
while (<$IN>) {
    chomp;
    my (undef, undef, undef, $lef, $nm, $rgt, undef) = split;
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

system(
    "Rscript $filedir/signature_count.R $envdir/$ARGV[0].tmp $envdir/$ARGV[1]");

system("rm $envdir/$ARGV[0].tmp");

__END__
