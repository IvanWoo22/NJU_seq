#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use POSIX;

sub SCORE {
    my ( $TR_END_COUNT, $NC_END_COUNT, $TR_TOTAL, $NC_TOTAL ) = @_;

    my %SCORE;
    my %END_COR;

    for my $CURRENT ( keys %$TR_END_COUNT ) {
        my ( $CHR, $DIR, $POS ) = split( /\t/, $CURRENT );
        my ( $FORMAL_POS, $END_COR );

        if ( $DIR eq "+" ) {
            $FORMAL_POS = $POS - 1;
        }
        else {
            $FORMAL_POS = $POS + 1;
        }

        my ( $N_END, $N_END_P1, $T_END, $T_END_P1 );

        $T_END = $TR_END_COUNT->{ $CHR . "\t" . $DIR . "\t" . $FORMAL_POS }
          // 1;
        $N_END = $NC_END_COUNT->{ $CHR . "\t" . $DIR . "\t" . $FORMAL_POS }
          // 1;

        $T_END_P1 = $TR_END_COUNT->{$CURRENT};
        $END_COR =
          POSIX::ceil( $TR_END_COUNT->{$CURRENT} * $NC_TOTAL / $TR_TOTAL );

        $N_END_P1 = $NC_END_COUNT->{$CURRENT} // 1;

        my $SCORE = $T_END_P1 / $T_END - $N_END_P1 / $N_END;
        $SCORE{$CURRENT}   = $SCORE;
        $END_COR{$CURRENT} = $END_COR;
    }

    return ( \%SCORE, \%END_COR );
}

my ( @end_count, @score, @end_count_cor, @total, %info, %all_site_id,
    %score_all );

open( my $IN_NC, "<", $ARGV[0] );

while ( my $line = <$IN_NC> ) {
    chomp $line;
    my @tmp = split /\t/, $line;
    my $id  = $tmp[0] . "\t" . $tmp[2] . "\t" . $tmp[1];
    $end_count[0]{$id} = $tmp[8];
}

close($IN_NC);

for my $sample ( 1 .. $#ARGV ) {
    $total[0] = 0;

    open( my $IN_TR, "<", $ARGV[$sample] );

    while ( my $line = <$IN_TR> ) {
        chomp $line;
        my @tmp = split /\t/, $line;
        my $id  = $tmp[0] . "\t" . $tmp[2] . "\t" . $tmp[1];

        if ( exists $all_site_id{$id} ) {
            $all_site_id{$id}++;
        }
        else {
            $info{$id}        = join( "\t", @tmp[ 3 .. 6 ] );
            $all_site_id{$id} = 1;
        }

        $end_count[$sample]{$id} = $tmp[8];
        $total[$sample] += $tmp[8];

        if ( exists $end_count[0]{$id} ) {
            $total[0] += $end_count[0]{$id};
        }
    }

    close($IN_TR);

    my ( $score_ref, $end_count_cor_ref ) = SCORE(
        \%{ $end_count[$sample] },
        \%{ $end_count[0] },
        $total[$sample], $total[0]
    );
    $score[$sample]         = $score_ref;
    $end_count_cor[$sample] = $end_count_cor_ref;
}

for my $id ( keys %all_site_id ) {
    if ( $all_site_id{$id} == $#ARGV ) {
        my ( $chr, $dir, $pos ) = split( /\t/, $id );
        my $NC_END_COUNT = exists $end_count[0]{$id} ? $end_count[0]{$id} : 0;
        my ( $SoaS, $SoaC ) = ( 0, 0 );
        my $SoaStv = 30 * $#ARGV;
        my $SoaCtv = 3 * $#ARGV;

        for my $sample ( 1 .. $#ARGV ) {
            $SoaC += $end_count_cor[$sample]{$id};
            $SoaS += $score[$sample]{$id};
        }

        if ( $SoaS >= $SoaStv && $SoaC >= $SoaCtv * $NC_END_COUNT ) {
            my $key = "$chr\t$pos\t$dir\t$info{$id}\t$NC_END_COUNT";
            $NC_END_COUNT = 1 if $NC_END_COUNT == 0;
            my $fc  = $SoaC / $NC_END_COUNT;
            my $fca = $fc / $#ARGV;
            for my $sample ( 1 .. $#ARGV ) {
                $key .=
"\t$end_count[$sample]{$id}\t$end_count_cor[$sample]{$id}\t$score[$sample]{$id}";
            }
            $key .= "\t$fc\t$fca";
            $score_all{$key} = $SoaS;
        }
    }
}

foreach my $key ( sort { $score_all{$b} <=> $score_all{$a} } keys %score_all ) {
    my $ScoreAve = $score_all{$key} / $#ARGV;
    print("$key\t$score_all{$key}\t$ScoreAve\n");
}

__END__
