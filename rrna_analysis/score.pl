#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use POSIX;
use File::ReadBackwards;

sub SCORE {
    my ( $TR_END_COUNT, $NC_END_COUNT, $TR_TOTAL, $NC_TOTAL ) = @_;
    my @SCORE;
    my @END_COR;

    $SCORE[0]   = 0;
    $SCORE[1]   = 0;
    $END_COR[0] = 0;
    $END_COR[1] = 0;

    for my $CURRENT ( 2 .. $#$TR_END_COUNT - 2 ) {
        my ( $N_END, $N_END_P1, $T_END, $T_END_P1, $END_COR );

        $T_END    = $TR_END_COUNT->[$CURRENT] || 1;
        $T_END_P1 = $TR_END_COUNT->[ $CURRENT + 1 ];
        $END_COR  = POSIX::ceil(
            $TR_END_COUNT->[ $CURRENT + 1 ] * $NC_TOTAL / $TR_TOTAL );

        $N_END    = $NC_END_COUNT->[$CURRENT] || 1;
        $N_END_P1 = $NC_END_COUNT->[ $CURRENT + 1 ];

        my $SCORE = $T_END_P1 / $T_END - $N_END_P1 / $N_END;
        push( @SCORE,   $SCORE );
        push( @END_COR, $END_COR );
    }

    push( @SCORE,   0 );
    push( @SCORE,   0 );
    push( @END_COR, 0 );
    push( @END_COR, 0 );

    for my $TER_START ( 0 .. 19 ) {
        $SCORE[$TER_START]               = 0;
        $END_COR[$TER_START]             = 0;
        $SCORE[ $#SCORE - $TER_START ]   = 0;
        $END_COR[ $#SCORE - $TER_START ] = 0;
    }

    return ( \@SCORE, \@END_COR );
}

my ( @site, @base, @end_count, @score, @end_count_cor, @total );
$total[0] = 0;
my $nc_content = do {
    open( my $fh, "<", $ARGV[0] );
    local $/;
    <$fh>;
};
my @nc_lines = split( "\n", $nc_content );

foreach my $line (@nc_lines) {
    chomp($line);
    my ( $site, $base, undef, $end_count ) = split /\t/, $line;
    push( @site,              $site );
    push( @base,              $base );
    push( @{ $end_count[0] }, $end_count );
    $total[0] += $end_count;
}

foreach my $sample ( 1 .. $#ARGV ) {
    my $tr_content = do {
        open( my $fh, "<", $ARGV[$sample] );
        local $/;
        <$fh>;
    };
    my @tr_lines = split( "\n", $tr_content );

    $total[$sample] = 0;
    foreach my $line (@tr_lines) {
        chomp($line);
        my ( undef, undef, undef, $end_count ) = split /\t/, $line;
        push( @{ $end_count[$sample] }, $end_count );
        $total[$sample] += $end_count;
    }
    ( $score[$sample], $end_count_cor[$sample] ) =
      SCORE( $end_count[$sample], $end_count[0], $total[$sample], $total[0] );
}

foreach my $site ( 0 .. $#site - 1 ) {
    print("$site[$site]\t$base[$site]\t$end_count[0][$site]\t");
    my ( $SoaS, $SoaC ) = ( 0, 0 );
    my $SoaS_tv = 30 * $#ARGV;
    my $SoaC_tv = 3 * $#ARGV;
    foreach my $sample ( 1 .. $#ARGV ) {
        $SoaC += $end_count_cor[$sample][$site];
        $SoaS += $score[$sample][$site];
        my $score_format = sprintf( "%.2f", $score[$sample][$site] );
        print("$end_count[$sample][$site]\t$score_format\t");
    }
    my ( $FC_LINE, $SC_LINE ) = ( "meet", "meet" );
    if ( $SoaS < $SoaS_tv ) {
        $SC_LINE = "unmeet";
    }
    if ( ( $SoaC < $SoaC_tv * $end_count[0][ $site + 1 ] ) or ( $SoaC == 0 ) ) {
        $FC_LINE = "unmeet";
    }
    my $SoaC_format        = $end_count[0][ $site + 1 ] || 1;
    my $fold_change_format = sprintf( "%.2f", $SoaC / $SoaC_format );
    my $score_format       = sprintf( "%.2f", $SoaS );
    print("$fold_change_format\t$FC_LINE\t$score_format\t$SC_LINE\n");
}

__END__
