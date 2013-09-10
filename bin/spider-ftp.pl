#!/usr/bin/perl -w

#####SETUP#####

use strict;
use Net::FTP;
use Getopt::Long;

my ($outdir, $configfile);
GetOptions(
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
           );

#####MAIN#####

my $config = parse_config($configfile);

my $f = Net::FTP->new($config->{'START'});
$f->login($config->{'USER'}, $config->{'PASS'});

$f->cwd('metagenomes/0/');
$f->get('4441026.3.processed.tar.gz');
#$, = "\n";
#my @files = $f->dir;
print yellow;

#####SUBS######

sub parse_config{
    my $file = shift;

    my %config;
    open (F, "$file");
    while (my $line = <F>){
        chomp $line;
        my ($key, $value) = split (/\=/, $line, 2);
        $config{$key} = $value;

    }
    close (F);

    return (\%config);
}
