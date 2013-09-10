#!/usr/bin/perl
    
use strict;
use warnings;
use GD::Graph::Cartesian;
#use GD::Graph::Polar;
use Graphics::ColorNames;

# create and pop GD object for drawing the frame                            
my $gdobj=GD::Graph::Cartesian->new (height=>800,
				     width=>800,
				     borderx=>5,
				     bordery=>5,
				     minx=>0,
				     miny=>0,
				     maxx=>$ARGV[1],
				     maxy=>$ARGV[2],
				     );
#my $highlights = {};
#if ($ARGV[3]){
#    open (H, "$ARGV[3]");
#    while (my $line = <H>){
#	chomp $line;
#	$highlights->{$line} = 1;
#    }
#    close (H);
#}

open (F, "$ARGV[0]");
while (my $line = <F>){
    chomp $line;
    my ($boid, $x, $y) = split (/\t/, $line);
#    if (exists($highlights->{$boid})){
#	$gdobj->color("red");
#    }
#    else {
#	$gdobj->color("black");
#    }
#    ($gdobj->color("red")) if (exists($highlights->{$boid}));
#    $gdobj->color("blue");
#    $gdobj->color("red");
    $gdobj->addPoint  ($x => $y);
#    $gdobj->color("black");
    $gdobj->addString ($x => $y, $boid);# if ($boid == 1);
#    $gdobj->color("black");
}
close (F);

print $gdobj->draw;
