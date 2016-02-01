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

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::GenBank;

my ($help, $list, $outdir);
GetOptions(
	   'h|help'       => \$help,
	   'l|list=s'     => \$list,
	   'o|outdir=s'   => \$outdir,
	   );

pod2usage(2) if $help;

`mkdir -p $outdir`;

my $gb = new Bio::DB::GenBank(-format => 'Fasta');

open (LIST, "$list");
while (my $line = <LIST>){
    chomp $line;
    my ($name, $acc) = split (/\s/, $line);
    
    print STDERR "Working on $line\n";

    my $seq_obj = $gb->get_Seq_by_acc($line);
#    my $seq_obj = $gb->get_Seq_by_acc($acc); 
    my $id = $seq_obj->id();
    my $desc = $seq_obj->desc();
    my $seq = $seq_obj->seq();
    
    #print ">$id\t$desc\n$seq\n";
    open (F, ">$outdir/$name");
    print F ">$name\n$seq\n";
    close (F);
}
close (LIST);
