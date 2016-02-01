#!/usr/bin/perl
=head1 NAME                                                                          
                                                                                     
pp_pa.pl                                                                      
                                                                                     
=head1 SYNOPSIS                                                                      
                                                                                     
  pp_pa.pl --                                                                 
    Collates flocks across multiple replicates of the clusterflock procedure to assess their 
    robustness. Produces a flock matrix where each column comprises a replicate and each row 
    an entity. Entities are assigned characters corresponding to their flock membership for that replicate.
                                                                                     
Options:                                                                             
                                                                                     
 --help        Show brief help and exit                                              
 --flockdir    Is a directory containing clusterflock output directories, one per replicate calculated 
 --indir       Is the directory of flocked entities (ortholog fastas)
 --burnstart   Is the first frame to mine
 --burnend     Is the last frame to mine
                                                                                     
=head1 DESCRIPTION                                                                   

Create a flock matrix to assess robustness across multiple clusterflock replicates
Note that the most common settings for burnstart and burnend are burnstart=brunend=index_of_final_frame.
So if your simulations last 1000 frames, both burnstart and burnend would be set to 1000 so only the final 
frame in each replicate is analyzed.

=head1 AUTHOR                                                                        
                                                                                     
Apurva Narechania                                                                    
anarechania *a|t* amnh.org                                                           
                                                                                     
=head1 COPYRIGHT                                                                     
                                                                                     
Copyright (c) 2015 American Museum of Natural History                                
                                                                                     
This library is free software;                                                       
you can redistribute it and/or modify                                                
it under the same terms as Perl itself.                                              
                                                                                     
=cut                                                                                 

# ----------------------------------------------------   

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($help, $flockdir, $indir, $burnstart, $burnend); #where indir is the directory from flocking
GetOptions(
    'h|help'          => \$help,
    'f|flockdir=s'    => \$flockdir,
    'i|indir=s'       => \$indir,
    's|startburn=s'   => \$burnstart,
    'e|endburn=s'     => \$burnend,
    ) or pod2usage;

pod2usage if $help;

for my $option ($indir, $flockdir, $burnstart, $burnend){
    (warn ("Missing a required option\n") and pod2usage)                            
        unless ($option);                                                           
}                                                                                   

####MAIN####                                                                         

# create a symbol library (each letter of the alphabet)
my $library = {};
my $c = -1;
my @letters;
for my $letter ("A".."Z"){
#my @letters = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');
#foreach my $letter (@letters){
    $c++;
    $library->{$c} = $letter;
    push (@letters, $letter);
    last if ($c == 26);
}
my $symbolstring = join " ", @letters;

# parse out the boids from fasta file names                                          
opendir (D, "$indir");
my @charsets = sort (readdir (D));
shift @charsets;
shift @charsets;

# parse the flocking info
opendir (F, "$flockdir");
my @repdirs = sort (readdir (F));
shift @repdirs;
shift @repdirs;
closedir (F);
my $repdirs = @repdirs;

my $pa = {};
foreach my $repdir (@repdirs){
    for (my $i = $burnstart; $i <= $burnend; $i++){
#	print STDERR "$repdir\t$i\n";
	
	# if the optics file exists, parse it,
	# otherwise insert missing data for a replicate
	# that failed at this point
	if (-e "$flockdir/$repdir/flocks/$i.optics"){
	    open (I, "$flockdir/$repdir/flocks/$i.optics");
	    while (my $line = <I>){
		chomp $line;
		my ($flock, $gene) = split (/\t/, $line);
		if ($flock > 26){
#		if ($flock > 20){
		    print STDERR "Too many flocks: not enough characters!\n";
		    die;
		}
		else{
		    push (@{$pa->{$gene}}, $library->{$flock});
#		    push (@{$pa->{$gene}}, $flock);
		}
	    }
	    close (I);
	}
	else {
	    foreach my $charset (@charsets){
		push (@{$pa->{$charset}}, "?");
	    }
	}
    }
}

my $ntax = 0;
foreach my $tax (sort keys %$pa){
    $ntax++;
}

print "#NEXUS\n";
print "BEGIN DATA;\n";
print "DIMENSIONS NTAX=$ntax NCHAR=$repdirs;\n";
print "FORMAT SYMBOLS=\"$symbolstring\";\n";
print "MATRIX\n";
my $freqs = {};
foreach my $tax (sort keys %$pa){
    my @paarray = @{$pa->{$tax}};

    my @paarray2;
    my $a = 0;
    my $b = 0;
    foreach my $pa (@paarray){
	if ($pa eq "A"){
	    $a++;
	    push (@paarray2, "A");
	}
	else {
	    $b++;
	    push (@paarray2, "B");
	}
    }
    my $freqa = $a / ($a + $b);
    $freqs->{$freqa}++;

    if ($freqa <= 0.8){
	print STDERR "$tax\t$a\t$b\t$freqa\n";
    } 
   
    my $pastring = join "", @paarray;
    my $pastring2 = join "", @paarray2;
    print "$tax\t$pastring\n";
#    print "$tax\t$pastring2\n";
}
print ";\n";
print "END;\n";

for (my $i = 0; $i <= 100; $i++){
    my $ith = $i / 100;
    if (exists ($freqs->{$ith})){
#	print STDERR "$ith\t$freqs->{$ith}\n";
    }
    else {
#	print STDERR "$ith\t0\n";
    }
}
