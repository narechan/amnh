#!/usr/bin/perl -w

#####SETUP#####

use strict;
use WWW::Mechanize;
use Getopt::Long;

my ($outdir, $configfile);
GetOptions(
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$configfile,
           );

#####MAIN#####

my $config = parse_config($configfile);

my $mech = WWW::Mechanize->new(autocheck => 1);
$mech->get($config->{'START'});

my @links = $mech->find_all_links( url_regex =>qr/$config->{'REGEX'}/);

for my $link ( @links ) {
    my $url = $link->url_abs;
    my $filename = $url;
    $filename =~ s[^.+/][];

    print "Fetching $url";
    $mech->get( $url, ':content_file' => $filename );

    print "   ", -s $filename, " bytes\n";
}

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
