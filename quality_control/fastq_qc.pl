#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Long;
use PerlIO::gzip;

my $dirname = dirname(__FILE__);

sub TAIL_BASE_COUNT {
    my %SEQ = @_;
    my ( $A_COUNT, $G_COUNT, $C_COUNT, $T_COUNT, $TOTAL_COUNT ) =
      ( 0, 0, 0, 0, 0 );
    foreach ( keys(%SEQ) ) {
        chomp( my $SEQ_STRING = $_ );
        $TOTAL_COUNT += $SEQ{$_};
        my @SEQ = split( //, $SEQ_STRING );
        if ( $SEQ[-1] eq "A" ) {
            $A_COUNT += $SEQ{$_};
        }
        elsif ( $SEQ[-1] eq "G" ) {
            $G_COUNT += $SEQ{$_};
        }
        elsif ( $SEQ[-1] eq "C" ) {
            $C_COUNT += $SEQ{$_};
        }
        elsif ( $SEQ[-1] eq "T" ) {
            $T_COUNT += $SEQ{$_};
        }
    }
    return ( $A_COUNT, $G_COUNT, $C_COUNT, $T_COUNT, $TOTAL_COUNT );
}

sub HEAD_BASE_COUNT {
    my %SEQ = @_;
    my ( $A_COUNT, $G_COUNT, $C_COUNT, $T_COUNT, $TOTAL_COUNT ) =
      ( 0, 0, 0, 0, 0 );
    foreach ( keys(%SEQ) ) {
        chomp( my $SEQ_STRING = $_ );
        $TOTAL_COUNT += $SEQ{$_};
        my @SEQ = split( //, $SEQ_STRING );
        if ( $SEQ[0] eq "A" ) {
            $A_COUNT += $SEQ{$_};
        }
        elsif ( $SEQ[0] eq "G" ) {
            $G_COUNT += $SEQ{$_};
        }
        elsif ( $SEQ[0] eq "C" ) {
            $C_COUNT += $SEQ{$_};
        }
        elsif ( $SEQ[0] eq "T" ) {
            $T_COUNT += $SEQ{$_};
        }
    }
    return ( $A_COUNT, $G_COUNT, $C_COUNT, $T_COUNT, $TOTAL_COUNT );
}

sub BODY_BASE_COUNT {
    my %SEQ = @_;
    my ( $A_COUNT, $G_COUNT, $C_COUNT, $T_COUNT, $TOTAL_COUNT ) =
      ( 0, 0, 0, 0, 0 );
    foreach ( keys(%SEQ) ) {
        chomp( my $SEQ_STRING = $_ );
        my @SEQ = split( //, $SEQ_STRING );
        foreach my $BASE (@SEQ) {
            if ( $BASE eq "A" ) {
                $A_COUNT     += $SEQ{$_};
                $TOTAL_COUNT += $SEQ{$_};
            }
            elsif ( $BASE eq "G" ) {
                $G_COUNT     += $SEQ{$_};
                $TOTAL_COUNT += $SEQ{$_};
            }
            elsif ( $BASE eq "C" ) {
                $C_COUNT     += $SEQ{$_};
                $TOTAL_COUNT += $SEQ{$_};
            }
            elsif ( $BASE eq "T" ) {
                $T_COUNT     += $SEQ{$_};
                $TOTAL_COUNT += $SEQ{$_};
            }
        }
    }
    return ( $A_COUNT, $G_COUNT, $C_COUNT, $T_COUNT, $TOTAL_COUNT );
}

my (
    @head_a,      @head_g,              @head_c,   @head_t,
    @head_sum,    @tail_a,              @tail_g,   @tail_c,
    @tail_t,      @tail_sum,            @body_a,   @body_g,
    @body_c,      @body_t,              @body_sum, @length_range,
    @reads_count, %length_distribution, @short,    @long
);

foreach my $sample ( 0 .. $#ARGV - 1 ) {
    open( my $in_fh, "<:gzip", $ARGV[$sample] );
    my %seq;
    $reads_count[$sample] = 0;
    while (<$in_fh>) {
        $reads_count[$sample]++;
        chomp( my $seq_string = <$in_fh> );
        readline($in_fh);
        readline($in_fh);
        if ( exists( $seq{$seq_string} ) ) {
            $seq{$seq_string}++;
        }
        else {
            $seq{$seq_string} = 1;
        }
    }

    (
        $body_a[$sample], $body_g[$sample], $body_c[$sample],
        $body_t[$sample], $body_sum[$sample]
    ) = BODY_BASE_COUNT(%seq);
    (
        $head_a[$sample], $head_g[$sample], $head_c[$sample],
        $head_t[$sample], $head_sum[$sample]
    ) = HEAD_BASE_COUNT(%seq);
    (
        $tail_a[$sample], $tail_g[$sample], $tail_c[$sample],
        $tail_t[$sample], $tail_sum[$sample]
    ) = TAIL_BASE_COUNT(%seq);

    my $short = 200;
    my $long  = 1;
    foreach my $seq_string ( keys(%seq) ) {
        my $length = length($seq_string);
        if ( exists( $length_distribution{ $ARGV[$sample] . "\t" . $length } ) )
        {
            $length_distribution{ $ARGV[$sample] . "\t" . $length } +=
              $seq{$seq_string};
        }
        else {
            $length_distribution{ $ARGV[$sample] . "\t" . $length } =
              $seq{$seq_string};
        }
        $short = $length if ( $short > $length );
        $long  = $length if ( $long < $length );
    }
    $length_range[$sample] = $short . " - " . $long;
    $short[$sample]        = $short;
    $long[$sample]         = $long;
}

my $name = join( "\t", @ARGV[ 0 .. $#ARGV - 2 ] );

my $body_a   = join( "\t", @body_a );
my $body_g   = join( "\t", @body_g );
my $body_c   = join( "\t", @body_c );
my $body_t   = join( "\t", @body_t );
my $body_sum = join( "\t", @body_sum );

my $head_a   = join( "\t", @head_a );
my $head_g   = join( "\t", @head_g );
my $head_c   = join( "\t", @head_c );
my $head_t   = join( "\t", @head_t );
my $head_sum = join( "\t", @head_sum );

my $tail_a   = join( "\t", @tail_a );
my $tail_g   = join( "\t", @tail_g );
my $tail_c   = join( "\t", @tail_c );
my $tail_t   = join( "\t", @tail_t );
my $tail_sum = join( "\t", @tail_sum );

my $reads_count  = join( "\t", @reads_count );
my $length_range = join( "\t", @length_range );

open( my $body,        ">", $ARGV[-2] . "/" . $ARGV[-1] . "_body.tsv" );
open( my $head,        ">", $ARGV[-2] . "/" . $ARGV[-1] . "_head.tsv" );
open( my $tail,        ">", $ARGV[-2] . "/" . $ARGV[-1] . "_tail.tsv" );
open( my $summary,     ">", $ARGV[-2] . "/" . $ARGV[-1] . "_summary.tsv" );
open( my $length_dist, ">", $ARGV[-2] . "/" . $ARGV[-1] . "_length.tsv" );

print $body ("$name\n$body_a\n$body_g\n$body_c\n$body_t\n$body_sum\n");

print $head ("$name\n$head_a\n$head_g\n$head_c\n$head_t\n$head_sum\n");

print $tail ("$name\n$tail_a\n$tail_g\n$tail_c\n$tail_t\n$tail_sum\n");

print $summary ("$name\n$reads_count\n$length_range\n");

foreach ( sort { $a cmp $b } keys %length_distribution ) {
    print $length_dist ("$_\t$length_distribution{$_}\n");
}

close($body);
close($head);
close($tail);
close($summary);
close($length_dist);

system(
    "Rscript $dirname/picture_draw.R \
    $ARGV[-2]/$ARGV[-1]_body.tsv \
    $ARGV[-2]/$ARGV[-1]_head.tsv \
    $ARGV[-2]/$ARGV[-1]_tail.tsv \
    $ARGV[-2]/$ARGV[-1]_length.tsv \
    $ARGV[-2]/$ARGV[-1]_summary.tsv \
    $ARGV[-1].pdf"
);

system("mv $ARGV[-1].pdf $ARGV[-2]/$ARGV[-1].pdf");

system("tar zcvf $ARGV[-1]_fastqc.tar.gz $ARGV[-1]_*.tsv");

__END__
