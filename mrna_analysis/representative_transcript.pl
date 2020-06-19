#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use Getopt::Long;

Getopt::Long::GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'geneid=s'  => \my $geneid,
    'transid=s' => \my $transid,
) or Getopt::Long::HelpMessage(1);

sub TAGS {
    my $INFO = shift;
    my $TSL;
    if ( $INFO =~ m/transcript_support_level=([1-5]);/ ) {
        $TSL = $1;
    }
    else {
        $TSL = 6;
    }
    my @TMP  = ( $INFO =~ m/tag=([A-Z,a-z,0-9,_,\,]*);/ );
    my @TAG  = split( /,/, $TMP[0] );
    my $MANE = 0;
    $MANE = grep /^MANE_Select$/i, @TAG;
    my $AP;
    if ( $INFO =~ m/appris_principal_([1-5])/ ) {
        $AP = $1;
    }
    else {
        $AP = 6;
    }
    my $BASIC = 0;
    $BASIC = grep /^basic$/i, @TAG;
    return ( $MANE, $AP, $BASIC, $TSL );
}

sub GET_ID {
    my ( $FEATURE, $INFO ) = @_;
    my $ID;
    if ( $INFO =~ m/$FEATURE(\w+)\.?[0-9]*;/ ) {
        $ID = $1;
    }
    else {
        warn("There is a problem in $INFO;\n");
    }
    return $ID;
}

my %gene_info;
while (<STDIN>) {
    chomp;
    my ( undef, undef, undef, $start, $end, undef, undef, undef, $info ) =
      split /\t/;
    my $trans_name = GET_ID( $transid, $info );
    my $gene_name  = GET_ID( $geneid,  $info );
    @{ $gene_info{$gene_name}{$trans_name} } = TAGS($info);
    my $length = $end - $start;
    push( @{ $gene_info{$gene_name}{$trans_name} }, $length );
}

foreach my $gene ( sort( keys(%gene_info) ) ) {
    my ( %ap, %basic, %tsl, %length );
    my $mane = 0;
    foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) ) {
        $mane++ if ${ $gene_info{$gene}{$trans} }[0] == 1;
        $ap{$trans}     = ${ $gene_info{$gene}{$trans} }[1];
        $basic{$trans}  = ${ $gene_info{$gene}{$trans} }[2];
        $tsl{$trans}    = ${ $gene_info{$gene}{$trans} }[3];
        $length{$trans} = ${ $gene_info{$gene}{$trans} }[4];
    }
    if ( $mane == 1 ) {
        foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) ) {
            print("$trans\n") if ${ $gene_info{$gene}{$trans} }[0] == 1;
        }
    }
    else {
        if ( $mane > 2 ) {
            foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) ) {
                if ( ${ $gene_info{$gene}{$trans} }[0] == 0 ) {
                    delete( $gene_info{$gene}{$trans} );
                    delete( $ap{$trans} );
                    delete( $basic{$trans} );
                    delete( $tsl{$trans} );
                    delete( $length{$trans} );
                }
            }
        }
        my @ap_sort = sort { $ap{$a} <=> $ap{$b} } keys %ap;
        if ( ( $#ap_sort == 0 ) or ( $ap{ $ap_sort[0] } < $ap{ $ap_sort[1] } ) )
        {
            print("$ap_sort[0]\n");
        }
        else {
            foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) ) {
                if ( $ap{$trans} > $ap{ $ap_sort[0] } ) {
                    delete( $gene_info{$gene}{$trans} );
                    delete( $ap{$trans} );
                    delete( $basic{$trans} );
                    delete( $tsl{$trans} );
                    delete( $length{$trans} );
                }
            }
            my $basic = 0;
            foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) ) {
                $basic++ if $basic{$trans} == 1;
            }
            if ( $basic == 1 ) {
                foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) ) {
                    print("$trans\n") if $basic{$trans} == 1;
                }
            }
            else {
                if ( $basic > 2 ) {
                    foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) )
                    {
                        if ( $basic{$trans} == 0 ) {
                            delete( $gene_info{$gene}{$trans} );
                            delete( $ap{$trans} );
                            delete( $basic{$trans} );
                            delete( $tsl{$trans} );
                            delete( $length{$trans} );
                        }
                    }
                }
                my @tsl_sort = sort { $tsl{$a} <=> $tsl{$b} } keys %tsl;
                if (   ( $#tsl_sort == 0 )
                    or ( $tsl{ $tsl_sort[0] } < $tsl{ $tsl_sort[1] } ) )
                {
                    print("$tsl_sort[0]\n");
                }
                else {
                    foreach my $trans ( sort( keys( %{ $gene_info{$gene} } ) ) )
                    {
                        if ( $tsl{$trans} > $tsl{ $tsl_sort[0] } ) {
                            delete( $gene_info{$gene}{$trans} );
                            delete( $ap{$trans} );
                            delete( $basic{$trans} );
                            delete( $tsl{$trans} );
                            delete( $length{$trans} );
                        }
                    }
                    my @length_sort =
                      sort { $length{$a} <=> $length{$b} } keys %length;
                    print("$length_sort[0]\n");
                }
            }
        }
    }
}

