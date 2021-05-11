#!/usr/bin/perl


use strict;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Data::Dumper;
use Bio::Tools::GFF;
use WWW::Curl::Simple;
use IO::String;

my $q = $ARGV[0];
my $gff_file = $ARGV[1];
my $fasta_file = $ARGV[2];
my $curl = WWW::Curl::Simple->new();

my ($chr, $coords) = split(':', $q);
my ($start, $end) = split('\.\.', $coords);
$start =~ s/\,//g;
$end =~ s/\,//g;
# print "Getting " . join("\t", $chr, $coords, $start, $end) . "\n";

my $seq;
my $gffio;
my @feat;

# http://toxodb.org/cgi-bin/gbrowse/toxodb/?plugin=FastaDumper;plugin_action=Go;view_start=1340976;view_stop=1345975


# Get FASTA file
#
if(!defined($fasta_file))
{
	my $fastalink = join("", 
		"http://toxodb.org/cgi-bin/gbrowse/toxodb/?plugin=FastaDumper;plugin_action=Go;FastaDumper.format=text;",
		"name=$chr:$start..$end;",
		"view_start=$start;view_stop=$end;"
		);
	my $res = $curl->get($fastalink);
    	if ($res->is_success) 
	{
		my $fastatxt = $res->decoded_content;
		# print $fastatxt;
		my $iofasta = IO::String->new($fastatxt);
		my $fasta = Bio::SeqIO->new(-format=>'fasta', -fh=>$iofasta);
		$seq = $fasta->next_seq();
    	} else
	{
		die("Can not connect and get fasta file!");
	}
#	open(TMP, 'tmp.fasta');
	
} else
{
	my $fasta = Bio::SeqIO->new(-format=>'fasta', -file=>$fasta_file);
	$seq = $fasta->next_seq;
}
if(!defined($gff_file))
{
	my $gfflink = join("", 
			"http://toxodb.org/cgi-bin/gbrowse/toxodb/?plugin=TrackDumper;plugin_action=Go;plugin_button=Cancel&plugin_button=Configure&plugin_button=Go&TrackDumper.mode=selected&TrackDumper.version=3&TrackDumper.region=All&TrackDumper.disposition=view&TrackDumper.print_config=off&plugin_button=Cancel&plugin_button=Configure&plugin_button=Go;",
			"name=$chr:$start..$end;",
			"view_start=$start;view_stop=$end;"
			);
	my $res = $curl->get($gfflink);
    	if ($res->is_success) 
	{
		my $gfftxt = $res->decoded_content;
		my $gffedited_text = "";
		my $past_header = 0;
		for (split /^/, $gfftxt) 
		{
			if($_ =~ /##gff-version 3/)
			{
				$past_header = 1;
			}
			if($past_header)
			{
				$gffedited_text .= $_;
			}
		}
		# Get rid of the text up until the gff-version 3 line
		# print $gffedited_text;
		my $iogff = IO::String->new($gffedited_text);
		$gffio = Bio::Tools::GFF->new(-fh=>$iogff, -gff_version => 3);
    	} else
	{
		die("Can not connect and get fasta file!");
	}

} else
{
	$gffio = Bio::Tools::GFF->new(-file => $gff_file, -gff_version => 3);
}


while(my $feat = $gffio->next_feature())
{
	my @labels;
	my $new_start = $feat->start-$start+1;
	my $new_end = $feat->end-$start+1;
	if($new_start < 1 && $new_end < 1)
	{
		next;
	}
	if($new_start > $seq->length && $new_end > $seq->length)
	{
		next;
	}
	if($new_start < 1)
	{
		$new_start = 1;
		if($feat->strand == 1)
		{
			push(@labels, "start truncated");
		} else
		{
			push(@labels, "end truncated");
		}
	}
	if($new_end > $seq->length)
	{
		$new_end = $seq->length;
		if($feat->strand == 1)
		{
			push(@labels, "end truncated");
		} else
		{
			push(@labels, "start truncated");
		}
	}
	#print join(":", $feat->start, $feat->end, $new_start, $new_end) . "\n";
	$feat->set_attributes(-start=>$new_start, -end=>$new_end);
	if($feat->has_tag("Name"))
	{
		push(@labels, join(" ", $feat->get_tag_values("Name")));
	}
	if($feat->has_tag("product"))
	{
		push(@labels, join(" ", $feat->get_tag_values("product")));
	}

	if(scalar @labels > 0)
	{
		$feat->add_tag_value("label", join(" - ", @labels));
		#print Dumper($feat);
	}
	push(@feat, $feat);
}

$seq->add_SeqFeature(@feat);

my $outstring;
my $ioout = IO::String->new($outstring);
my $seqout = Bio::SeqIO->new(-format => 'genbank', -fh=>$ioout);
#my $seqout = Bio::SeqIO->new(-format=>'genbank', -file=>'>outfile.gb');
$seqout->write_seq($seq);

print $outstring;

