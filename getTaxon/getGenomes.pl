#!/usr/bin/perl

use Bio::DB::Taxonomy;
use LWP::Simple;
use strict;

my $db = Bio::DB::Taxonomy->new(-source => 'entrez');


my $taxon_cache_file = ".taxoninfo.txt";
my $bioproject_cache_file = ".bioprojectinfo.txt";
my $default_save_dir = "genomes";
mkdir($default_save_dir);

my $use_cache = 1;
my %taxids;
my %hasRefseq;
my %taxoncache;
my %type;

# Get all bioprojects

my $ftp_summary = 'ftp://ftp.ncbi.nlm.nih.gov/bioproject/summary.txt';
my $summary_text = '';
my %downloaded;

if(-e $bioproject_cache_file && ! -z $bioproject_cache_file)
{
	my $filesize = -s $bioproject_cache_file;
	print "Using bioproject cache file: " . -$filesize . " bytes\n";
	open(SUMTXT, $bioproject_cache_file);
	while(<SUMTXT>)
	{
		$summary_text .= $_;
	}
} else
{
	print "Downloading bioproject summary file.\n";
	my $summary_text = get($ftp_summary);
	print $summary_text;
	open(SUMTXT, ">", $bioproject_cache_file);
	print SUMTXT $summary_text;
	close(SUMTXT);
	my $filesize = -s $bioproject_cache_file;
	print "Done downloading file size: " . $filesize . " bytes\n";
}

if(-e $taxon_cache_file)
{
	my $num_records = 0;
	open(TAXONCACHE, $taxon_cache_file);
	while(<TAXONCACHE>)
	{
		my $line = $_;
		chomp($line);
		my ($taxid, @adesc) = split("\t", $line);
		my $desc = join("\t", @adesc);
		$taxoncache{$taxid} = $desc;
		$num_records++;
	}
	close(TAXONCACHE);
	print "$num_records taxon records loaded from cache.\n";
}

open(TAXONCACHE, ">", $taxon_cache_file);
select((select(TAXONCACHE), $|=1)[0]);
while(my ($taxid, $desc) = each(%taxoncache))
{
	print TAXONCACHE join("\t", $taxid, $desc) . "\n";
}


my @all_bioprojects = split("\n", $summary_text);

# Only look at Genomes
# Field 6, "Genome sequencing" or "RefSeq Genome"


foreach my $bioproject (@all_bioprojects)
{
	my ($organism_name, $taxid, $prjacc, $prjid, $prjtype, $prjdata, $date) = split("\t", $bioproject);
	if($prjdata eq "Genome sequencing" || $prjdata eq "RefSeq Genome")
	{
		my $shortName =  getShortName($organism_name);
		print "Checking $organism_name [$taxid] ($shortName):";
		# Am I a Euk?
		if(isEuk($taxid))
		{
			print " is a Eukaryote";
			$taxids{$taxid} = $organism_name;
			if($prjdata eq 'RefSeq Genome')
			{
				print " has refseq genome";
				$hasRefseq{$taxid} = 1;
			} else
			{
				print " is not a refseq genome";
			}
		} else
		{
			print " is not a Eukaryote";
		}
		print "\n";
	} else
	{
		# Do nothing
	}

}


# Get all taxon ids and names:

# If there is a refseq genome, download that
# "protein bioproject"[Filter] AND $organism[orgn] AND refseq

foreach my $taxid (keys %hasRefseq)
{
	# Filename, $taxid, $orgname, $type
	my $organism_name = $taxids{$taxid};
	my $filename = $default_save_dir . "/" . getShortName($organism_name) . ".fasta";
	my $fastaSeq = getFasta($organism_name, "refseq");

	if(length($fastaSeq)  > 0)
	{
		if(! -e $filename || -z $filename)
		{
			open(FILE, ">", $filename);
			print FILE $fastaSeq;
			close(FILE);
		}
		$downloaded{$taxid};
	}

}

# If no refseq genome, download other non-refseq
# "protein bioproject"[Filter] AND $organism[orgn] AND txid$taxon_id

# Determine naming scheme First letter of phyla, species name, remove special chars and spaces

sub getTaxonName
{
	my $taxonid = shift;
	my $lineprint = '';

	if(exists($taxoncache{$taxonid}))
	{
		return $taxoncache{$taxonid};
	}

	my $taxon = $db->get_taxon(-taxonid => $taxonid);
	if(!defined($taxon))
	{
		return undef;
	}
	my @common_names = $taxon->common_names();

	$lineprint .= $taxon->division() . "\t";
	$lineprint .= join(",", @common_names) . "\t";
	
	my $t = $taxon;
	my @all_taxa = ();
	while(defined($t))
	{
		unshift(@all_taxa,  @{$t->name('scientific')});
		my $a = $db->ancestor($t);
		$t = $a;
	}
	$lineprint .= join(",", @all_taxa);
	$taxoncache{$taxonid} = $lineprint;
	print TAXONCACHE join("\t", $taxonid, $lineprint) . "\n";
	return $lineprint;

}

sub isEuk
{
	my $taxonid = shift;
	if( getTaxonName($taxonid) =~ /,Eukaryota/)
	{
		return 1;
	} else
	{
		return 0;
	}

}

sub getShortName
{
	my $name = shift;

	$name =~ s/\.//g;
	$name =~ s/\-//g;
	$name =~ s/\://g;
	$name =~ s/\+//g;
	$name =~ s/\(//g;
	$name =~ s/\)//g;
	$name =~ s/\\//g;
	$name =~ s/\///g;
	$name =~ s/\,//g;
	$name =~ s/\=//g;
	$name =~ s/\#//g;
	$name =~ s/\'//g;
	$name =~ s/\`//g;
	$name =~ s/\"	//g;

	my @parts = split(" ", $name);
	my $shortname = '';
	$shortname .= uc(substr($parts[0], 0, 1));
	shift(@parts);
	$shortname .= join("_", @parts);
	return $shortname;
}

sub getFasta
{
	my $taxid = shift;
	my $organism_name = shift;
	my $extra = shift;

	my $db = "protein";
	my $seq_text = '';
	
	my $query = '"protein bioproject"[Filter] AND ' . $organism_name . '[orgn]';
	if(defined($taxid))
	{
	 	$query .= ' AND txid' . $taxid;
	}

	if(defined($extra))
	{
		$query .= " AND $extra";
	}

	#assemble the esearch URL
	my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
	my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

	my $output = get($url);


	#parse WebEnv, QueryKey and Count (# records retrieved)
	my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
	my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
	my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);


	#retrieve data in batches of 500
	my $retmax = 500;
	for (my $retstart = 0; $retstart < $count; $retstart += $retmax) 
	{
	        my $efetch_url = $base ."efetch.fcgi?db=$db&WebEnv=$web";
	        $efetch_url .= "&query_key=$key&retstart=$retstart";
	        $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
	        my $efetch_out = get($efetch_url);
	        $seq_text .= "$efetch_out";
	}


	return $seq_text;
}














