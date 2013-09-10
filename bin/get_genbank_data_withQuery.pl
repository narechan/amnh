#! /usr/bin/perl -w

=head1 NAME

get_genbank_data.pl

=head1 SYNOPSIS

  get_genbank_data.pl

Options:

 --help     Show brief help and exit
 --list     List of accessions to query
 --outdir   Is your outdir

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

my $queries = {};
open (L, "$list");
while (my $line = <L>){
    chomp $line;
    $line =~s/\[accn\]//g;
    my ($acc, $id, $query) = split (/\t/, $line);
    $queries->{$query} = $acc;
}
close (L);
    
my $gb = new Bio::DB::GenBank;

foreach my $query (keys %$queries){
    print STDERR "Working on $query\n";
    
    my ($start, $end) = split (/\-/, $query);
    
    # specific to WGS records names
    my ($stchars, $ststring) = $start =~ /^(.{5})(.*)$/s;
    my ($endchars, $endstring) = $end =~ /^(.{5})(.*)$/s;

    open (FO, ">$outdir/$queries->{$query}");
    my $counter = $ststring;
    until ($counter > $endstring){
	my $acc = $stchars . $counter;
	
	my $seq_obj = $gb->get_Seq_by_acc($acc);
	my $id = $seq_obj->id();
	my $seq = $seq_obj->seq();

	print FO ">$id\n$seq\n";
	$counter++;
    }
    close (FO);
}


=head
foreach my $query (keys %$queries){
    next if ($query eq "AACK01000001:AACK01000140[accn]");
    my $q = Bio::DB::Query::GenBank->new
	(-query   => $query,
	 -db      => 'nucleotide');
    
    my $seqio = $gb->get_Stream_by_query($q);
 
    open (FO, ">$outdir/$queries->{$query}.fa");
    while(my $seq = $seqio->next_seq){
	my $id = $seq->id();
	my $sequence = $seq->seq();
	
	print FO ">$id\n$sequence\n";
    }
    close (FO);
}
=cut
