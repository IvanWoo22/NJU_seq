#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use POSIX;
use File::Basename;
use File::Spec;

sub SCORE {
    my ( $TR_END_COUNT, $NC_END_COUNT, $TR_TOTAL, $NC_TOTAL ) = @_;
    my @SCORE = (0) x scalar(@$TR_END_COUNT);
    my @END_COR = @SCORE;
    return (\@SCORE, \@END_COR) unless $TR_TOTAL > 0 && $NC_TOTAL > 0;
    for my $current (2 .. $#$TR_END_COUNT - 2) {
        my $t_end = $TR_END_COUNT->[$current] || 1;
        my $n_end = $NC_END_COUNT->[$current] || 1;
        $SCORE[$current] = $TR_END_COUNT->[$current+1]/$t_end - $NC_END_COUNT->[$current+1]/$n_end;
        $END_COR[$current] = POSIX::ceil($TR_END_COUNT->[$current+1] * $NC_TOTAL / $TR_TOTAL);
    }
    for my $i (0 .. $#SCORE) {
        if ($i < 20 || $i > $#SCORE-20) { $SCORE[$i] = $END_COR[$i] = 0; }
    }

    return ( \@SCORE, \@END_COR );
}

@ARGV >= 2 or die "Usage: score.pl NC.count TR.count [TR2.count ...]\n";
my @sample_names;

for my $i ( 0 .. $#ARGV ) {
    my $file = $ARGV[$i];
    my $dir = File::Basename::dirname($file);
    my $prefix = File::Basename::basename($dir);

    push @sample_names, $prefix;
}

my ( @site, @base, @end_count, @score, @end_count_cor, @total );
my $normalization_ok = 1;
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
    defined($end_count) && $site == @site+1 && $end_count =~ /^\d+$/
        or die "NC must contain contiguous 1-based positions and nonnegative counts\n";
    push( @site,              $site );
    push( @base,              $base );
    push( @{ $end_count[0] }, $end_count );
    $total[0] += $end_count;
}

@site > 0 or die "Empty reference count table\n";
foreach my $sample ( 1 .. $#ARGV ) {
    my $tr_content = do {
        open( my $fh, "<", $ARGV[$sample] );
        local $/;
        <$fh>;
    };
    my @tr_lines = split( "\n", $tr_content );

    @tr_lines == @nc_lines or die "Reference lengths differ between NC and TR\n";
    my $i = 0;
    $total[$sample] = 0;
    foreach my $line (@tr_lines) {
        chomp($line);
        my ($pos, $base_tr, undef, $end_count) = split /\t/, $line;
        defined($end_count) && $pos == $site[$i] && $base_tr eq $base[$i] && $end_count =~ /^\d+$/
            or die "TR reference/coordinate/count mismatch at row " . ($i+1) . "\n";
        $i++;
        push( @{ $end_count[$sample] }, $end_count );
        $total[$sample] += $end_count;
    }
    ( $score[$sample], $end_count_cor[$sample] ) =
      SCORE( $end_count[$sample], $end_count[0], $total[$sample], $total[0] );
    if ($total[$sample] == 0 || $total[0] == 0) {
        warn "Normalization unavailable for $sample_names[$sample]: NC_TOTAL=$total[0], TR_TOTAL=$total[$sample]; joint candidates suppressed\n";
        $normalization_ok = 0;
    }
}

# Legacy output uses zero/unmeet for an unavailable joint analysis. The Rust
# unified schema emits explicit status and NA instead of these compatibility values.
unless ($normalization_ok) {
    for my $sample (1 .. $#ARGV) {
        $score[$sample] = [(0) x scalar(@site)];
        $end_count_cor[$sample] = [(0) x scalar(@site)];
    }
}

print "#Position\tBase\t${sample_names[0]}_count";
foreach my $s ( 1 .. $#ARGV ) {
    print "\t${sample_names[$s]}_count\t${sample_names[$s]}_score";
}
print "\tFoldChange\tFC_criteria\tTotalScore\tScore_criteria\n";

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
