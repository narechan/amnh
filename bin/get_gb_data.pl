#! /usr/bin/perl -w

use strict;
use Bio::DB::GenBank;

my $gb = new Bio::DB::GenBank(-format => 'Fasta');

my $seq_obj = $gb->get_Seq_by_acc($ARGV[0]);
my $id = $seq_obj->id();
my $desc = $seq_obj->desc();
my $seq = $seq_obj->seq();
print ">$id\t$desc\n$seq\n";
