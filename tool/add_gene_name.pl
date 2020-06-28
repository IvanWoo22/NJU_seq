#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;

Getopt::Long::GetOptions(
    'help|h' => sub { Getopt::Long::HelpMessage(0) },
    'id=s'   => \my $id_feature,
    'name=s' => \my $name_feature,
    'file=s' => \my $file,
    'col=s'  => \my $col,
) or Getopt::Long::HelpMessage(1);

sub GET_INFO {
    my ( $ID_FEATURE, $NAME_FEATURE, $INFO ) = @_;
    my ( $ID, $NAME );
    if ( $INFO =~ m/$ID_FEATURE(\w+)\.?[0-9]*/ ) {
        $ID = $1;
    }
    else {
        warn("There is a problem in $INFO;\n");
    }
    if ( $INFO =~ m/$NAME_FEATURE([A-Z,a-z,0-9,\-,\.,\(,\), ,\[,\],%,:,\/,\+,',`]+);/ ) {
        $NAME = $1;
    }
    else {
        warn("There is a problem in $INFO;\n");
    }
    return ( $ID, $NAME );
}

my %gene_info;

while (<STDIN>) {
    chomp;
    my ( undef, undef, undef, undef, undef, undef, undef, undef, $info ) =
      split /\t/;
    my ( $id, $name ) = GET_INFO( $id_feature, $name_feature, $info );
    $gene_info{$id} = $name;
}

open( my $REF_FH, "<", $file );
while (<$REF_FH>) {
    chomp;
    my @line = split "\t";
    my @id   = split "/", $line[6];
    my $name;
    if ( $#id > 0 ) {
        my @name = map { $gene_info{$_} } @id;
        $name = join( "/", @name );
    }
    else {
        $name = $gene_info{ $id[0] };
    }
    my $out =
        join( "\t", @line[ 0 .. $col - 2 ] ) . "\t"
      . $name . "\t"
      . join( "\t", @line[ $col - 1 .. $#line ] );
    print "$out\n";
}
close($REF_FH);
__END__
