#!/usr/bin/perl -w

=head1 NAME

R_wrapper.pl

=head1 SYNOPSIS

  R_wrapper.pl -- this program wraps around an R script 
    that requires only one data file, if R is invoked like this:

    Invoke % R --slave --args [infile] < mm.R

    Multiple infiles are run through the R script using this
    wrapper.

Options:

 --help        Show brief help and exit
 --indir       Contains file to be analyzed
 --outdir      Is your output dir
 --script      Is your R script to run

R must be in your path

=head1 DESCRIPTION

Run R script on a series of input files

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

#####
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;

my ($help, $indir, $outdir, $script);
GetOptions(
    'h|help'          => \$help,
    'i|indir=s'       => \$indir,
    'o|outdir=s'      => \$outdir,
    's|script=s'      => \$script,
    ) or pod2usage;

pod2usage if $help;

for my $option ($indir, $outdir, $script){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/results`;
`mkdir -p $outdir/stderr`;
`mkdir -p $outdir/subfiles`;

#####MAIN#####

# get all the files
opendir (I, "$indir");
my @files = sort (readdir (I));
shift @files;
shift @files;
closedir (I);

# cycle through the files
foreach my $file (@files){
    open (F, "$indir/$file");
    
    # store all the data in the file
    my %data;
    while (my $line = <F>){
	chomp $line;
	my @line = split (/\t/, $line);
	$data{$line[0]} = [@line];
    }

    # create subfiles with the last column chopped off
    my $counter = 0;
    while (%data){
	$counter++;
	my $covcache;

	open (W, ">$outdir/subfiles/$file.$counter");
	foreach my $cov (sort {$a <=> $b} keys %data){
	    $covcache = $cov;
	    my $string = join "\t", @{$data{$cov}};
	    print W "$string\n";
	}
	close (W);

	delete ($data{$covcache});
    }
    close (F);
}
    
# get all the subfiles                                                                                 
opendir (S, "$outdir/subfiles");
my @subfiles = sort (readdir (S));
shift @subfiles;
shift @subfiles;
closedir (S);


foreach my $subfile (@subfiles){
    warn "R $subfile\n";
    R_run ($subfile, $outdir, $script);
}

# get all the results
opendir (R, "$outdir/results");
my @rfiles = sort (readdir (R));
shift @rfiles;
shift @rfiles;
closedir (R);

# tabulate the Vmax and Km
foreach my $rfile (@rfiles){
    my @rfile = split (/\./, $rfile);
    open (RF, "$outdir/results/$rfile");
    while (my $line = <RF>){
	chomp $line;
	if ($line =~m/^Vm/){
	    my @data = split (/\s+/, $line);
	    print "$rfile[1]\t$rfile[2]\t$data[1]\n";
	}
    }
    close (RF);
}
	
#####SUBS#####

sub R_run{
    my $file  = shift;
    my $out   = shift;
    my $script = shift;

    # R
    my $r = "R";
    $r .= " --slave";
    $r .= " --args $out/subfiles/$file";
    $r .= " < $script";
    $r .= " 1>$out/results/$file.R_result";
    `$r 2>$out/stderr/$file.R_stderr`;
}
