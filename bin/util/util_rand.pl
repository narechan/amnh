#!/usr/bin/perl

my @pairs = ("aba", "ab", "a", "b", "c");
while (1){
    my $rand = int(rand(@pairs));
    print "$rand $pairs[$rand]\n";
}
