#!/usr/bin/perl
=head1 NAME                                                                       
                                                                                  
boids_pp.pl
                                                                                  
=head1 SYNOPSIS                                                                   
                                                                                  
  boids_pp.pl --                                                           
    cut out boids given x and y ranges
                                                                                  
Options:                                                                          
                                                                                  
 --help        Show brief help and exit                                           
 --infile      Is the set of boids coordinates you want to analyze
 --xrange      Is the xrange coordinates to cut out; comma separated, no spaces
 --yrange      Is the yrange coordinates to cut out; comma separated, no spaces   

=head1 DESCRIPTION                                                                
                                                                                  
cut boids out of an image
                                                                                  
=head1 AUTHOR                                                                     
                                                                                  
Apurva Narechania                                                                 
anarechania *a|t* amnh.org                                                        
                                                                                  
=head1 COPYRIGHT                                                                  
                                                                                  
Copyright (c) 2012 American Museum of Natural History                             
                                                                                  
This library is free software;                                                    
you can redistribute it and/or modify                                             
it under the same terms as Perl itself.                                           
                                                                                  
=cut                                                                              

# ----------------------------------------------------                            

#####                                                                             
#####                                                                             

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($help, $infile, $xrange, $yrange);
GetOptions(
    'h|help'          => \$help,
    'i|infile=s'      => \$infile,
    'x|xrange=s'      => \$xrange,
    'y|yrange=s'      => \$yrange,
	   ) or pod2usage;

pod2usage if $help;

for my $option ($infile, $xrange, $yrange){
    (warn ("Missing a required option\n") and pod2usage)
        unless ($option);
}

####MAIN####

# parse the boids coordinates file and find relavent groups
my $boidskept = {};
my ($x1, $x2) = split (/\,/, $xrange);
my ($y1, $y2) = split (/\,/, $yrange);
open (B, "$infile");
while (my $line = <B>){
    chomp $line;
    my ($gene, $x, $y) = split (/\t/, $line);
    if (($x >= $x1) and ($x <= $x2) and ($y >= $y1) and ($y <= $y2)){
	$boidskept->{$gene} = 1;
	print "$gene\n";
    }
    else {
	next;
    }
}
close (B);
