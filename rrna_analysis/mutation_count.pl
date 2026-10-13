#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

@ARGV == 3 or die "Usage: mutation_count.pl reference.fa five_columns.tsv reference_id\n";
open my $REF,    "<", $ARGV[0];
open my $IN_SAM, "<", $ARGV[1];

sub CIGAR_CONVERT {
    my $CIGAR = shift;
    $CIGAR =~ /^(?:[1-9][0-9]*[=XMID])+$/ or die "Unsupported CIGAR\n";
    my $SEQ   = shift;
    my @SEQ   = split( "", $SEQ );
    my @ALIGNMENT;
    while ( $CIGAR =~ /^([0-9]+)/ ) {
        my $LENGTH = $1;
        $CIGAR =~ s/^([0-9]+)//;
        if ( $CIGAR =~ /^[=XM]/ ) {
            $CIGAR =~ s/^[=XM]//;
            foreach ( 1 .. $LENGTH ) {
                push @ALIGNMENT, shift(@SEQ);
            }
        }
        elsif ( $CIGAR =~ /^I/ ) {
            $CIGAR =~ s/^I//;
            foreach ( 1 .. $LENGTH ) {
                shift(@SEQ);
            }
        }
        elsif ( $CIGAR =~ /^D/ ) {
            $CIGAR =~ s/^D//;
            foreach ( 1 .. $LENGTH ) {
                push @ALIGNMENT, "";
            }
        }
    }
    @SEQ == 0 && !grep { !defined $_ } @ALIGNMENT or die "CIGAR/sequence length mismatch\n";
    return \@ALIGNMENT;
}

my @site;
my @site_base;

readline($REF);
my $ref;
while (<$REF>) {
    s/\r?\n//;
    /^>/ and die "Expected one FASTA reference\n";
    $ref .= uc $_;
}
close($REF);

@site_base = split( "", $ref );

foreach ( 0 .. $#site_base ) {
    ${ $site[$_] }{"A"} = 0;
    ${ $site[$_] }{"G"} = 0;
    ${ $site[$_] }{"C"} = 0;
    ${ $site[$_] }{"T"} = 0;
    ${ $site[$_] }{"N"} = 0;
    ${ $site[$_] }{""} = 0;
}

while (<$IN_SAM>) {
    chomp;
    my @row = split /\t/;
    @row == 5 or die "Expected five alignment columns\n";
    if ( $row[1] eq $ARGV[2] ) {
        my ( undef, undef, $position, $cigar, $seq ) = split( "\t", $_ );
        my $alignment = CIGAR_CONVERT( $cigar, uc $seq );
        my @alignment = @{$alignment};
        $position =~ /^[1-9][0-9]*$/ && $position + @alignment - 1 <= @site_base or die "Alignment outside reference\n";
        foreach ( 0 .. $#alignment ) {
            ${ $site[ $position - 1 + $_ ] }{ $alignment[$_] } += 1;
        }
    }
}
close($IN_SAM);

# Stable ten columns: position, ref, A,G,C,T,N, consensus, status, deletions.
foreach my $i (0 .. $#site_base) {
    my $counts = $site[$i];
    my @bases = sort { $counts->{$b} <=> $counts->{$a} || $a cmp $b } qw(A G C T N);
    my $max_count = $counts->{$bases[0]};
    my @ties = grep { $counts->{$_} == $max_count } @bases;
    my ($consensus, $status) = ('N', 'uncovered');
    if ($max_count > 0) {
        $consensus = @ties == 1 ? $ties[0] : 'N';
        $status = $consensus eq 'N' ? 'ambiguous' : $consensus eq $site_base[$i] ? 'match' : 'mut';
    }
    print join("\t", $i+1, $site_base[$i], @{$counts}{qw(A G C T N)}, $consensus, $status, $counts->{''}), "\n";
}
