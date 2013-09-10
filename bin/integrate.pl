#!/usr/bin/perl

# Purpose: Accept XY coordinates and integrate using trapezoidal rule of 
#          straight points and those interpolated using a loess curve.
#          There must be the same number of X and Y values.

# Options:
#          -i input file

# use the following modules 
use Getopt::Std;
use Math::Spline;

# HARD CODED location of the loess R script
my $loess = "/home/apurva/bin/loess.R";

# get inputs                                                                    
my %opts = ();
getopts ('i:', \%opts);
my $infile = $opts{'i'};

# integrate the straight coordinates
my ($x, $y) = parse ($infile);
my $area = integrate_trap ($x, $y);
print "$area\n";

# integrate the loess coordinates
`R --slave --args 0.75 $infile $infile.Rout < $loess`;
my ($xl, $yl) = parse ("$infile.Rout");
my $areal = integrate_trap ($xl, $yl);
print "$areal\n";

###SUBS###

# parse the file
sub parse{
    my $infile = shift;
    
    my @x;
    my @y;
    open (FILE, "$infile");
    while (my $line = <FILE>){
	chomp $line;
	my ($x, $y) = split (/\t/, $line);
	push (@x, $x);
	push (@y, $y);
    }
    close (FILE);

    return (\@x, \@y);
}

# integrate with the trapezoidal rule
sub integrate_trap{
    my $x = shift;
    my $y = shift;

    my $points = @$x;
    my $integration;
    for (my $i = 0; $i < $points - 1; $i++){
	my $height = @$x[$i+1] - @$x[$i];
	my $area   = 0.5 * $height * (@$y[$i+1] + @$y[$i]);
	$integration += $area;
    }
    
    return ($integration);
}

