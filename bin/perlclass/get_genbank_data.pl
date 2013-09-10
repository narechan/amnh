#! /usr/bin/perl -w

=head1 NAME

get_genbank_data.pl

=head1 SYNOPSIS

  get_genbank_data.pl

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

# load libraries
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::GenBank;

# get commandline options
my ($help, $list);
GetOptions(
	   'h|help'       => \$help,
	   'l|list=s'     => \$list,
	   );

pod2usage(2) if $help;

# instantiate genbank object
my $gb = new Bio::DB::GenBank(-format => 'Fasta');

# foreach item on the list get the 
# sequence record by accession
open (LIST, "$list");
while (my $line = <LIST>){
    chomp $line;
    
    print STDERR "Working on $line\n";
    
    # get the sequence object and associated data
    my $seq_obj = $gb->get_Seq_by_acc($line);
    my $id      = $seq_obj->id();
    my $desc    = $seq_obj->desc();
    my $seq     = $seq_obj->seq();
    
    # print in fasta format
    print ">$id\t$desc\n$seq\n";
}
close (LIST);
