#!/usr/bin/perl


use strict;
use Bio::SeqIO;
use File::Temp qw/tempfile tempdir/;

# This file will take a file of orfs in fasta format and print out their orfid and where they localize using IPSORT's algorithm
# As well as inserting the results into the database
# Usage: ./parse_targetp.pl giardia giardia_orfs.fasta


my $signalp_bin = 'targetp -N -s 0.0 -t 0.65 -o 0.52';

my $in_file = $ARGV[0];

my $orfs = Bio::SeqIO->new('-file'=> $in_file, '-format'=>'Fasta');
my $mito_out = Bio::SeqIO->new(-file=>">" . $in_file . '.mito', -format=>'fasta');
my ($seq_fh, $seq_filename) = tempfile();
my ($output_fh, $output_filename) = tempfile();

while(my $orf = $orfs->next_seq() )
{
	my $seqout = Bio::SeqIO->new(-file=>">" . $seq_filename, -format=>'fasta');
	$seqout->write_seq($orf);
	my ($name) = $orf->display_id;

	my $output;
	$output = system("$signalp_bin $seq_filename > $output_filename");

	open(ANALYSIS, $output_filename);

	my $full_report = '';
	my $score = 0;
	my $header = '';
	my $length = '';
	my $mtp = '';
	my $sp = '';
	my $other = '';
	my $loc = '';
	my $rc = '';

	while(<ANALYSIS>)
	{
		my $line = $_;
		$full_report .= $line;
		
		if($line =~ /^--/)
		{
			$line = <ANALYSIS>;
			$full_report .= $line;
			($header, $length, $mtp, $sp, $other, $loc, $rc) = $line =~ /^(.+)\s+(\d+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([A-Za-z\_\?]+)\s+(\d+)\s*$/;
			chomp($header, $length, $mtp, $sp, $other, $loc, $rc);

			if($loc =~  /\_/)
			{
				$loc = 'other';
			} elsif($loc =~ /\*/)
			{
				$loc = 'undetermined';
			} elsif($loc =~ /\?/)
			{
				$loc = 'undetermined';
			}

			my $score = 0;
			if($loc eq 'other')
			{
				$score = $other;
			}elsif($loc eq 'M')
			{
				$score = $mtp;
			} elsif($loc eq 'S')
			{
				$score = $sp;
			} elsif($loc eq 'undetermined')
			{
				$score = 0;
			} else
			{
				$loc = 'undetermined';
				$score = 0;
			}
			
			print join("\t", $name, 'targetp',$loc, $score) . "\n";
			if($loc eq 'M')
			{
				$mito_out->write_seq($orf);
			}

			while(<ANALYSIS>)
       			{ 		
				$full_report .= $_;
			}
		}

	}
}



