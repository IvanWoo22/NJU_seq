#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

my ( @line, @gene, @gene_range, @gene_position, @mark );
my $i = 0;
while (<>) {
    chomp;
    push( @line, $_ );
    my ( $chr, undef, $type, $start, $end, undef, $dir, undef, undef ) =
      split /\t/;
    if ( $type eq "gene" ) {
        push( @gene, $i );
        my $position = $chr . "\t" . $dir;
        push( @gene_position, $position );
        my $range = AlignDB::IntSpan->new;
        $range->add_range( $start, $end );
        push( @gene_range, $range );
    }
    $i++;
}

# Compare all same-strand intervals, including nested and nonadjacent genes.
my %groups;
push @{$groups{$gene_position[$_]}}, $_ for 0 .. $#gene;
@mark = (0) x scalar(@gene);
for my $group (values %groups) {
    my ($max_end, $max_id);
    for my $id (sort { $gene_range[$a]->min <=> $gene_range[$b]->min } @$group) {
        my ($start, $end) = ($gene_range[$id]->min, $gene_range[$id]->max);
        if (defined($max_end) && $start <= $max_end) {
            $mark[$id] = $mark[$max_id] = 1;
        }
        if (!defined($max_end) || $end > $max_end) { ($max_end, $max_id) = ($end, $id); }
    }
}

foreach ( 0 .. $#gene ) {
    if ( ( $mark[$_] == 0 ) and ( $_ < $#gene ) ) {
        foreach my $t ( $gene[$_] .. $gene[ $_ + 1 ] - 1 ) {
            print("$line[$t]\n");
        }
    }
    elsif ( ( $mark[$_] == 0 ) and ( $_ == $#gene ) ) {
        foreach my $t ( $gene[$_] .. $#line ) {
            print("$line[$t]\n");
        }
    }
}

__END__
