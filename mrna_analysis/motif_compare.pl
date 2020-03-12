#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open(my $MIRNA_FH, "<", $ARGV[0]);
open(my $NM_FH, "<", $ARGV[1]);

my %mirna_motif;
while(<$MIRNA_FH>){
    chomp;
    my($name, $motif) = split /\t/;
    if(exists($mirna_motif{$motif})){
        $mirna_motif{$motif} .= ",$name";
    }else{
        $mirna_motif{$motif} = $name;
    }
}
close($MIRNA_FH);
while(<$NM_FH>){
    chomp;
    my (undef,undef,undef,undef,undef,undef,undef,undef,$motif) = split/\t/;
    if(exists($mirna_motif{$motif})){
        my $mirna_num = split(/,/,$mirna_motif{$motif});
        print("$_\t$mirna_num\t$mirna_motif{$motif}\n");
    }else{
        print("$_\t0\tNULL\n");
    }
}
close($NM_FH);
__END__