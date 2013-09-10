#!/usr/bin/perl
    
use strict;
use warnings;
use Chart::Gnuplot;
use Getopt::Long;
use Pod::Usage;

my ($help, $itlog, $dimmensions, $highfile, $outgif);
GetOptions(
    'h|help'          => \$help,
    'i|itlog=s'      => \$itlog,
    'd|dimmensions=s'      => \$dimmensions,
    'h|highlights=s'  => \$highfile,
    'o|outgif=s'     => \$outgif,
           ) or pod2usage;

pod2usage if $help;


# create the gnuplot object for drawing the frame                            
my $chart = Chart::Gnuplot->new(output => $outgif,
#				xtics => undef,
#				ytics=>undef,
				xtics => {
#				    labels    => [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250],
				    labels   => [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500],
				},
				ytics => {
#				    labels    => [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250],
                                    labels   => [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500],
				},
				xrange=>[0, $dimmensions],
				yrange=>[0, $dimmensions],
#				grid=>'on',
#				grid   => {
#				    width => 100
#				    },
				size=>'square 2');

my $highlights = {};
if ($highfile){
    open (H, "$highfile");
    while (my $line = <H>){
        chomp $line;
        $highlights->{$line} = 1;
    }
    close (H);
}

my @x;
my @y;
my @xhigh;
my @yhigh;
open (F, "$itlog");
while (my $line = <F>){
    chomp $line;
    my ($boid, $x, $y) = split (/\t/, $line);
    
    if (exists($highlights->{$boid})){
	push (@xhigh, $x);
	push (@yhigh, $y);
    }
    else {
	push (@x, $x);
	push (@y, $y);
    }
}
close (F);

if ($highfile){
    my $data = Chart::Gnuplot::DataSet->new(xdata=>\@x,
					    ydata=>\@y,
					    style=>"points",
					    pointtype=>"fill-circle",
					    pointsize=>0.5);
    my $datahigh = Chart::Gnuplot::DataSet->new(xdata=>\@xhigh,
						ydata=>\@yhigh,
						style=>"points",
						pointtype=>"fill-circle",
						pointsize=>0.5);
    
    $chart->plot2d($data, $datahigh);
}
else {
    my $data = Chart::Gnuplot::DataSet->new(xdata=>\@x,
                                            ydata=>\@y,
                                            style=>"points",
                                            pointtype=>"fill-circle",
                                            pointsize=>0.5);
    $chart->plot2d($data);
}
