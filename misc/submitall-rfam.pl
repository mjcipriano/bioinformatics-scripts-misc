#!/usr/bin/perl


use strict;

use Bio::SeqIO;
use Bio::Seq;

my $infile = 'giardia12';
my $infile = 'unused_reads.fas';

my $seqs = Bio::SeqIO->new(-file=>$infile, -format=>"fasta");

while(my $seq = $seqs->next_seq())
{
	my $name = $seq->display_id;
	open(FASTA, ">", $seq->display_id . ".fas");
	open(SUBMIT, ">", $seq->display_id . ".sh");
	print FASTA ">" . $seq->display_id . "\n" . $seq->seq() . "\n";
	close(FASTA);
	print SUBMIT '#!/bin/csh
#$ -j y
#$ -o /xraid/habitat/mcipriano/rfam/5S/rfam.out
#$ -N rfams_' . $seq->display_id . '
#$ -cwd

cd /xraid/habitat/mcipriano/rfam/5S;
/opt/BioBrew/infernal/bin/cmsearch 5S.cm ' . $seq->display_id . '.fas > ' . $seq->display_id . '.rfam;
echo "DONE";';

	close(SUBMIT);
	open(FASTA, ">>", $seq->display_id . ".fas");

	for(1..499)
	{
		if($seq = $seqs->next_seq())
		{
			print FASTA ">" . $seq->display_id . "\n" . $seq->seq() . "\n";
		}
	}
	close(FASTA);
	system('chmod 755 /xraid/habitat/mcipriano/rfam/5S/' . $name . '.sh');
	system('cd /xraid/habitat/mcipriano/rfam/5S/;qsub ' . $name . '.sh');
}
