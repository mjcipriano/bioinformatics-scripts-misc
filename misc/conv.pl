#!/usr/bin/perl

use Bio::SeqIO;


$in  = Bio::SeqIO->new(-file => "chrI.embl",
                         -format => 'embl');
$out = Bio::SeqIO->new(-file => ">chrI.gb",
                         -format => 'genbank');

while ( my $seq = $in->next_seq() ) {$out->write_seq($seq); }

#while (my $seq = $in->next_seq) {
#		for my $feat ($seq->top_SeqFeatures) {
#			print $feast->gff_string,"\n";
#		}
#	}

