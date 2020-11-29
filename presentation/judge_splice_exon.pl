#!/usr/bin/env perl
use warnings;
use strict;
use autodie;
use AlignDB::IntSpan;

open( my $IN1, "<", $ARGV[0] );

my ( %ase_s, %ase_e, %ase_id, %ase_gene );
while (<$IN1>) {
    chomp;
    my @tmp = split "\t";
    $tmp[0] =~ s/chr//;
    push( @{ $ase_s{ $tmp[0] . $tmp[5] } },  $tmp[1] );
    push( @{ $ase_e{ $tmp[0] . $tmp[5] } },  $tmp[2] );
    push( @{ $ase_id{ $tmp[0] . $tmp[5] } }, $tmp[4] );
    if ( exists( $ase_gene{ $tmp[0] . $tmp[5] . $tmp[4] } ) ) {
        $ase_gene{ $tmp[0] . $tmp[5] . $tmp[4] } += $tmp[2] - $tmp[1] + 1;
    }
    else {
        $ase_gene{ $tmp[0] . $tmp[5] . $tmp[4] } = $tmp[2] - $tmp[1] + 1;
    }
}
close($IN1);

my ( %cse_s, %cse_e, %cse_id, %cse_gene );
open( my $IN2, "<", $ARGV[1] );
while (<$IN2>) {
    chomp;
    my @tmp = split "\t";
    $tmp[0] =~ s/chr//;
    push( @{ $cse_s{ $tmp[0] . $tmp[5] } },  $tmp[1] );
    push( @{ $cse_e{ $tmp[0] . $tmp[5] } },  $tmp[2] );
    push( @{ $cse_id{ $tmp[0] . $tmp[5] } }, $tmp[4] );
    if ( exists( $cse_gene{ $tmp[0] . $tmp[5] . $tmp[4] } ) ) {
        $cse_gene{ $tmp[0] . $tmp[5] . $tmp[4] } += $tmp[2] - $tmp[1] + 1;
    }
    else {
        $cse_gene{ $tmp[0] . $tmp[5] . $tmp[4] } = $tmp[2] - $tmp[1] + 1;
    }
}
close($IN2);

my ( $ase_times, $cse_times, $ase_ranges, $cse_ranges ) = ( 0, 0, 0, 0 );
my %gene;
open( my $IN3, "<", $ARGV[2] );
while (<$IN3>) {
    chomp;
    my @tmp = split "\t";
    my ( $ase, $cse ) = ( 0, 0 );
    foreach my $no ( 0 .. $#{ $ase_s{ $tmp[0] . $tmp[2] } } ) {
        if (    ( $tmp[1] >= ${ $ase_s{ $tmp[0] . $tmp[2] } }[$no] )
            and ( $tmp[1] <= ${ $ase_e{ $tmp[0] . $tmp[2] } }[$no] ) )
        {
            $ase++;
            unless (
                exists(
                    $gene{
                            $tmp[0]
                          . $tmp[2]
                          . ${ $ase_id{ $tmp[0] . $tmp[2] } }[$no]
                    }
                )
              )
            {
                $ase_ranges +=
                  $ase_gene{ $tmp[0]
                      . $tmp[2]
                      . ${ $ase_id{ $tmp[0] . $tmp[2] } }[$no] }
                  if (
                    exists(
                        $ase_gene{
                                $tmp[0]
                              . $tmp[2]
                              . ${ $ase_id{ $tmp[0] . $tmp[2] } }[$no]
                        }
                    )
                  );
                $cse_ranges +=
                  $cse_gene{ $tmp[0]
                      . $tmp[2]
                      . ${ $ase_id{ $tmp[0] . $tmp[2] } }[$no] }
                  if (
                    exists(
                        $cse_gene{
                                $tmp[0]
                              . $tmp[2]
                              . ${ $ase_id{ $tmp[0] . $tmp[2] } }[$no]
                        }
                    )
                  );
                $gene{  $tmp[0]
                      . $tmp[2]
                      . ${ $ase_id{ $tmp[0] . $tmp[2] } }[$no] } = 1;
            }
        }
    }
    foreach my $no ( 0 .. $#{ $cse_s{ $tmp[0] . $tmp[2] } } ) {
        if (    ( $tmp[1] >= ${ $cse_s{ $tmp[0] . $tmp[2] } }[$no] )
            and ( $tmp[1] <= ${ $cse_e{ $tmp[0] . $tmp[2] } }[$no] ) )
        {
            $cse++;
            unless (
                exists(
                    $gene{
                            $tmp[0]
                          . $tmp[2]
                          . ${ $cse_id{ $tmp[0] . $tmp[2] } }[$no]
                    }
                )
              )
            {
                $ase_ranges +=
                  $ase_gene{ $tmp[0]
                      . $tmp[2]
                      . ${ $cse_id{ $tmp[0] . $tmp[2] } }[$no] }
                  if (
                    exists(
                        $ase_gene{
                                $tmp[0]
                              . $tmp[2]
                              . ${ $cse_id{ $tmp[0] . $tmp[2] } }[$no]
                        }
                    )
                  );
                $cse_ranges +=
                  $cse_gene{ $tmp[0]
                      . $tmp[2]
                      . ${ $cse_id{ $tmp[0] . $tmp[2] } }[$no] }
                  if (
                    exists(
                        $cse_gene{
                                $tmp[0]
                              . $tmp[2]
                              . ${ $cse_id{ $tmp[0] . $tmp[2] } }[$no]
                        }
                    )
                  );
                $gene{  $tmp[0]
                      . $tmp[2]
                      . ${ $cse_id{ $tmp[0] . $tmp[2] } }[$no] } = 1;
            }
        }
    }
    if ( $ase + $cse > 0 ) {
        $ase_times += $ase / ( $ase + $cse );
        $cse_times += $cse / ( $ase + $cse );
    }
}

print "$ase_times\n$cse_times\n$ase_ranges\n$cse_ranges\n";

__END__