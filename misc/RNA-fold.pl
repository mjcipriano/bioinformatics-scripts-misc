#!/usr/bin/perl


use strict;

use Bio::SeqIO;

my $input_fasta = $ARGV[0];

my $out_dir = `pwd`;
chomp($out_dir);

my ($root_name) = $input_fasta =~ /^(.+)\./;

print $root_name . "\n\n";

my $seqs = Bio::SeqIO->new( -file=>$input_fasta, -format=>'fasta');

while (my $seq = $seqs->next_seq)
{
	my $seq_name = $seq->display_id;
	my $range = $seq->desc;
	my ($range_name) = $range =~ /(\d+\-\d+)/;
	my $new_name = $seq_name . '_' . $range_name;
	my $sequence =  $seq->seq;
	my $string = ">" . $new_name . "\n" . $sequence . "\n";

	my $dir = "$out_dir/$root_name";
	mkdir($dir);
	system("cd $dir;echo '$string' | /xraid/bioware/linux/ViennaRNA/bin/RNAfold > $new_name.fold");
	my $short_name = substr($new_name, 0, 12);
	system("cd $dir;mv $short_name" . "_ss.ps" . " $new_name.ps");
}
