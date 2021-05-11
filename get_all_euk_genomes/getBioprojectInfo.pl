#!/usr/bin/perl

use strict;

use Bio::DB::Taxonomy;

my $db = Bio::DB::Taxonomy->new(-source => 'entrez');
use LWP::Simple;

#my $genomefile = "/data/projects/kinetochore/data2/genomes/genomelist.txt";
my $genomefile = $ARGV[0];

open(GENOMEFILE, $genomefile);

my @taxon_prjids = ();
my %taxon_info;
my @labels = ();
my %label2row;
my @dbs = ("protein", "nucleotide");

my $debug = 1;


while(my $line = <GENOMEFILE>)
{
	chomp($line);
	if( $line =~ /^\#/)
	{
		$line =~ s/\#//;
		@labels = split("\t", $line);
		my $lid = 0;
		foreach my $label (@labels)
		{
			$label2row{$label} = $lid;
			$lid++;
		}
	} else
	{
		my @thisarr = split("\t", $line);
		$taxon_info{$thisarr[2]} = \@thisarr;
		push(@taxon_prjids, $thisarr[2]);
	}
}

my %taxoncache = {};

foreach my $prjid (@taxon_prjids)
{
	my $lineprint = '';

	my $taxonid = $taxon_info{$prjid}->[1];

	print join("\t", @{$taxon_info{$prjid}}) . "\t";
	my $protein_count = 0;
	my $nt_count = 0;
	my $nt_size = 0;

	foreach my $db (@dbs)
	{
		if($db eq 'nucleotide' && $protein_count != 0)
		{
			$nt_count = getProjectCount($db, $prjid);
		#	$nt_size = getProjectSum($db, $prjid);

		} elsif($db eq 'protein')
		{
			$protein_count = getProjectCount($db, $prjid);
		} else
		{
		}

	}
	print join("\t", $protein_count, $nt_count, $nt_size);
	print "\n";

}

my $query = '72423[BioProject]';
my $db = 'protein';

# http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&sort=&query_key=2&qty=18384&filter=all


sub getProjectCount
{
	my $db = shift;
	my $query = shift;

	$query = $query . '[BioProject]';

	#assemble the esearch URL
	my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
	my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y&DbFrom=bioproject";

	#post the esearch URL
	my $output = get($url);

	if($debug)
	{
		print "=======================================\n$output\n==========================================\n";
	}
	#parse WebEnv, QueryKey and Count (# records retrieved)
	my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
	my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
	my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

	return $count;
}


sub getProjectSum
{
	my $db = shift;
	my $query = shift;

	$query = $query . '[BioProject]';

	#assemble the esearch URL
	my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
	my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

	#post the esearch URL
	my $output = get($url);

	# print "=======================================\n$output\n==========================================\n";
	#parse WebEnv, QueryKey and Count (# records retrieved)
	my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
	my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
	my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

	my $retmax = 500;
	my $sum = 0;

	if($count > 0)
	{
		for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
			my $efetch_url = $base ."efetch.fcgi?db=$db&WebEnv=$web";
			$efetch_url .= "&query_key=$key&retstart=$retstart";
			$efetch_url .= "&retmax=$retmax&rettype=summary&retmode=xml";
			my $efetch_out = get($efetch_url);
			my @lengths = $efetch_out =~ /<Seq-inst_length>(\d+)<\/Seq-inst_length>/g;
			$sum += $_ for @lengths;
		#	print join(":", @lengths) . "\n";
			#print  "$efetch_out";
		}
	}


	return $sum;
}

