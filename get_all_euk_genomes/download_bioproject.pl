#!/usr/bin/perl


use strict;

use LWP::Simple;
# my $query = 'chimpanzee[orgn]+AND+biomol+mrna[prop]';
my $query = '72423[BioProject]';
my $db = 'protein';

# http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&sort=&query_key=2&qty=18384&filter=all


#assemble the esearch URL
my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";
#my $url = $base . "esearch.fcgi?db=nuccore&term=$query&usehistory=y";



#post the esearch URL
my $output = get($url);
print $output;

#parse WebEnv, QueryKey and Count (# records retrieved)
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);


#open output file for writing
open(OUT, ">chimp.fna") || die "Can't open file!\n";


#retrieve data in batches of 500
my $retmax = 500;
for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
        my $efetch_url = $base ."efetch.fcgi?db=$db&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
        my $efetch_out = get($efetch_url);
        print OUT "$efetch_out";
}
close OUT;

