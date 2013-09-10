#!/usr/bin/perl

$filenameA = $ARGV[0];
$filenameOut1 = $ARGV[1];
$filenameOut2 = $ARGV[2];

open $FILEA, "< $filenameA";

open $OUTFILE1, "> $filenameOut1";
open $OUTFILE2, "> $filenameOut2";

while(<$FILEA>) {
    chomp $_;

    if ($_ =~m/1$/){
	print $OUTFILE1 "$_\n";
	$_ = <$FILEA>;
	print $OUTFILE1 $_; 
	$_ = <$FILEA>;
        print $OUTFILE1 $_;
	$_ = <$FILEA>;
	print $OUTFILE1 $_;
    }
    else {
	print $OUTFILE2 "$_\n";;
        $_ = <$FILEA>;
        print $OUTFILE2 $_;
        $_ = <$FILEA>;
        print $OUTFILE2 $_;
        $_ = <$FILEA>;
        print $OUTFILE2 $_;
    }
}
close ($FILEA);
close ($OUTFILE1);
close ($OUTFILE2);

