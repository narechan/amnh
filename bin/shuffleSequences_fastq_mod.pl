#!/usr/bin/perl

$filenameA = $ARGV[0];
$filenameB = $ARGV[1];
$filenameOut = $ARGV[2];

open $FILEA, "< $filenameA";
open $FILEB, "< $filenameB";

open $OUTFILE, "> $filenameOut";

while(<$FILEA>) {
    chomp $_;
    
    my ($acc, $junk) = split (/\s/, $_);
    #    $_ =~s/\s/\-/g;
#    print $OUTFILE "$_#0/1\n";
    print $OUTFILE "$acc#0/1\n";
    $_ = <$FILEA>;
    print $OUTFILE $_; 
    $_ = <$FILEA>;
#    chomp $_;
    print $OUTFILE $_; #"$_#0/1d\n"; 
    $_ = <$FILEA>;
    print $OUTFILE $_; 
    
    $_ = <$FILEB>;
    chomp $_;
#    $_ =~s/\s/\-/g;
#    print $OUTFILE "$_#0/2\n"; 
    print $OUTFILE "$acc#0/2\n";
    $_ = <$FILEB>;
    print $OUTFILE $_;
    $_ = <$FILEB>;
#    chomp $_;
    print $OUTFILE $_; #"$_#0/2\n";
    $_ = <$FILEB>;
    print $OUTFILE $_;
}
