#!/usr/bin/perl

use Bio::Matrix::PSM::IO;


my $in_filename = 'GABP_transfac.transfac';
my $out_filename = 'GABP_mast.mast';

my $psmin = Bio::Matrix::PSM::IO->new(-format=>'transfac', -file=>$in_filename) or die ("Huh?");
my $psmout = Bio::Matrix::PSM::IO->new(-format=>'masta', -file=>$out_filename, -mtype=>'SEQ');

my $psm = $psmin->next_psm();
print $psm->width();

# my $psm_header = $psm->header;
# print $psm_header{IC};
# print $psm_header{sites};

print $psm;

# $psm->id("Test");
 $psmout->write_psm($psm);

