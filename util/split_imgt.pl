#!/usr/bin/perl

use lib "./";
use strict;
use fasta;

my $id;
my $info;
my $sequence;

if(open(my $in, '<', $ARGV[0]))
{
    while((($id, $info, $sequence) = fasta::ReadFasta($in)) && ($id ne ""))
    {
        my @fields = split(/\|/, $info);
        if(!($fields[1] =~ /^TRAC/))
        {
            my $hl     = ($fields[1] =~ /^IGH/)?'heavy':'light';
        
            my $species = $fields[2];
            $species =~ s/\s/-/g;
            $species =~ s/_.*//;
            my $domain  = "\L$fields[4]";
            $domain =~ s/-region//;
            $domain =~ s/ch.*/c/;
            my $filename = "${species}_${hl}_${domain}_dna.faa";
            print "$filename\n";
        }
    }
}

