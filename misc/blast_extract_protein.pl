#!/usr/bin/perl

use Getopt::Long;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SeqIO;
use Bio::AlignIO;
use IO::String;

use strict;

my $dbs;
my $query_file;

my $options = GetOptions(       "database=s", \$dbs,
                                "query_file=s", \$query_file
			);

my @params = (     'database'    => "$dbs",
		   'program'     => "tblastn"
		   );

my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
my $seqio = Bio::SeqIO->new(-file=>$query_file, -format=>'fasta');
while(my $seq = $seqio->next_seq)
{

	my $blast_report = $factory->blastall($seq);

	my $this_result = $blast_report->next_result;

	while( my $hit = $this_result->next_hit()) 
	{
		while(my $hsp = $hit->next_hsp())
		{

			# Get the HSP
			print join("\t", $hit->name(), $hsp->evalue(),  $hsp->hit->start, $hsp->hit->end, $hsp->hit->strand, $hsp->hit->length) . "\n";
			print join("\n", $hsp->query_string, $hsp->homology_string, $hsp->hit_string) . "\n";

			# Now I want to take this area and expand out and find the gene that covers it
			my $contig_seq = get_sequence($dbs, $hit->name());
			if(defined($contig_seq))
			{
				print $contig_seq->display_id . "\n";

				# Pull out the sequences from start to end
				my $sub_seq; 
				my $trans_seq;
				if($hsp->hit->strand == -1)
				{
					$trans_seq = $contig_seq->trunc($hsp->hit->start, $hsp->hit->end)->revcom->translate;
					$sub_seq = $contig_seq->trunc($hsp->hit->start, $hsp->hit->end)->revcom;
				} else
				{
					$trans_seq = $contig_seq->trunc($hsp->hit->start, $hsp->hit->end)->translate;
					$sub_seq = $contig_seq->trunc($hsp->hit->start, $hsp->hit->end);
				}
				print $trans_seq->seq() . "\n\n";



			} else
			{
				print "No SEQ\n";
			}
			

			
		}
	}
}



sub get_sequence
{
	my $db = shift;
	my $seqname = shift;

	my $seq = `fastacmd -d $db -p F -s $seqname`;
	my $iofh = IO::String->new($seq);
	my $sio = Bio::SeqIO->new(-format=>'fasta', -fh=>$iofh);
	my $seq = $sio->next_seq();
	if(defined($seq))
	{
		return $seq;
	} else
	{
		return undef;
	}
}
