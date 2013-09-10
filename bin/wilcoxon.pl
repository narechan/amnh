#!/usr/bin/perl -w

=head1 NAME

wilcoxon.pl

=head1 SYNOPSIS

  wilcoxon.pl

Options:

 --help        Show brief help and exit
 --file1           Is your first distribution of AUCs
 --file2           Is your second distribution of AUCs

=head1 DESCRIPTION

Given two distributions -- calculates the nonpara probability
that two samples are from the same distribution

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2008 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####TODO:
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Test::WilcoxonRankSum;

my ($help, $file1, $file2);
GetOptions(
    'h|help'          => \$help,
    'file1|a=s'      => \$file1,
    'file2|b=s'        => \$file2,
    ) or pod2usage;

pod2usage if $help;

for my $option ($file1, $file2){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

#####MAIN#####

# parse the AUC files and store data
# they have to have the same number of nodes
# and each node index must correspond
my $data1 = parse ($file1);
my $data2 = parse ($file2);

foreach my $node (sort {$a <=> $b} keys %$data1){

    my $wilcoxtest = Statistics::Test::WilcoxonRankSum->new();
    $wilcoxtest->load_data ($data1->{$node}, $data2->{$node});
    my $prob = $wilcoxtest->probability();

    print "$node\t$prob\n";
}

#####SUBS#####

sub parse {
    my $file1 = shift;

    my $data = {};
    open (A, "$file1");
    while (my $line = <A>){
	chomp $line;
	
	my @data = split (/\t/, $line);
	my $sample = shift @data;
	
	my $counter = -1;
	foreach my $element (@data){
	    $counter++;
	    
	    push @{$data->{$counter}}, $element;
	}
    }
    close (A);
    
    return ($data);
}
