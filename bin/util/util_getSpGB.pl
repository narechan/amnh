#! /usr/bin/perl -w

=head1 NAME

util_getSpGB.pl

=head1 SYNOPSIS

  util_getSpGB.pl

Options:

 --help     Show brief help and exit
 --list     List of accessions to query

=head1 DESCRIPTION

Given a list of genbank accessions, finds and writes the sequence

=head1 SEE ALSO

perl.

=head1 AUTHOR

Apurva Narechania

=head1 COPYRIGHT

Copyright (c) 2007 Cold Spring Harbor Laboratory

This library is free software;  you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

####SETUP####

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::GenBank;
use Bio::SeqIO;

my ($help, $list, $outdir);
GetOptions(
	   'h|help'       => \$help,
	   'l|list=s'     => \$list,
	   );

pod2usage(2) if $help;

`mkdir -p $outdir`;

####MAIN####

my $gb = Bio::DB::GenBank->new(-format => 'Fasta');

open (LIST, "$list");
my $sequence = {};
while (my $line = <LIST>){
    chomp $line;
    
    print STDERR "Working on $line\n";

    # get the record from remote genbank
    my $seq_obj = $gb->get_Seq_by_acc($line);
    my $id = $seq_obj->id();
    my $desc = $seq_obj->desc();
    
    print "$id\t$desc\n";
    
}
