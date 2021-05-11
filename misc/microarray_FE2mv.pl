#!/usr/bin/perl

use strict;

# Takes as input a filename and outputs filename.mev and filename.ann


my $file = $ARGV[0];

my $probe_desc_file = 'agilent_probe_annotation.txt';

my $probename_desc;
my $probename_symbol;
my $probename_genename;
my $probename_acc;
my $probename_sys;

open(PROBEDES, $probe_desc_file);
while(<PROBEDES>)
{
	my $line = $_;
	chomp($line);
	my @line = split("\t", $line);

	my $probe = $line[1];
	$probename_desc->{$probe} = $line[5];
	$probename_symbol->{$probe} = $line[2];
	$probename_acc->{$probe} = $line[4];
	$probename_sys->{$probe} = $line[0];
	$probename_genename->{$probe} = $line[3];
}
# Now get the TIGR TC annotations

my $tc_tog;
my $tc_go_term;
my $tc_tentative_annotation;
my $tc_et_name;
my $tc_oligo;

my $tc_file = 'zf_tc_annotation.txt';
open(TCFILE, $tc_file);
my $header = <TCFILE>;
while(<TCFILE>)
{
	my $line = $_;
	chomp($line);
	my ($tc, $tog, $go, $tc_et, $ten_annotation, $oligo) = split("\t", $line);
	$tc_tog->{$tc} = $tog;
	$tc_go_term->{$tc} = $go;
	$tc_tentative_annotation->{$tc} = $ten_annotation;
	$tc_et_name->{$tc} = $tc_et;
	$tc_oligo->{$tc} = $oligo;
}

open(INFILE, $file) or die ("No file to open!\n");
# Open the out files
open(MEV, ">", $file . ".mev");
print MEV join("\t", "UID", "IA", "IB", "R", "C", "MR", "MC", "BkgA", "BkgB", "SDA", "SDB", "SDBkgA", "SDBkgB", "MedA", "MedB", "MNA", "MNB") . "\n";

open(ANN, ">", $file . ".ann");
print ANN join("\t", "UID", "R", "C", "accessions", "location", "control_type", "probe_name", "gene_name", "systematic_name", "description", "description2", "gene_symbol", "gene_name2", "additional_acc", "systematic_name", "tc", "tog", "tc_go", "tc_et_name", "tc_annotation", "tc_oligo", "organism_uid") . "\n";	



# Search until the feature line 
my $go = 1;

while($go)
{
	my $row = <INFILE>;
	if($row =~ /^FEATURES/)
	{
		$go = 0;
	}
}

my $uid_hash;
while(<INFILE>)
{
	my $row = $_;
	# skip rows that do not start with DATA
	if( !($row =~ /^DATA/) )
	{
		next;
	}
	chomp($row);
	my @line = split("\t", $row);

	my $uid = $line[8];
	my $row = $line[2];
	my $col = $line[3];
	my $accessions = $line[4];
	my $location = $line[5];
	my $control_type = $line[9];
	my $probe_name = $line[10];
	my $gene_name = $line[11];
	my $systematic_name = $line[12];
	my $description = $line[13];
	my $pos_x = $line[14];
	my $pos_y = $line[15];
	my $IA = $line[23]; # Green Signal
	my $IB = $line[24]; # Red Signal
	my $SA = $line[31]; # Num Pixels
	my $MNA = $line[33];
	my $MNB = $line[34];
	my $MedA = $line[35];
	my $MedB = $line[36];
	my $SDA = $line[37];
	my $SDB = $line[38];
	my $BkgA = $line[41];
	my $BkgB = $line[42];
	my $SDBkgA = $line[45];
	my $SDBkgB = $line[46];
	my $NetA = $line[93];
	my $NetB = $line[94];

	if($control_type != 0)
	{
		next;
	}
	# Create MEV file

	print MEV join("\t", "$uid", "$IA", "$IB", "$row", "$col", "1", "1", "$BkgA", "$BkgB", "$SDA", "$SDB", "$SDBkgA", "$SDBkgB", "$MedA", "$MedB", "$MNA", "$MNB") . "\n";
	
#	$uid_hash->{$uid} = join("\t", "$uid", "$row", "$col", "$accessions", "$location", "$control_type", "$probe_name", "$gene_name", "$systematic_name", "$description")
	# "tc", "tog", "tc_go", "tc_annotation"
	# Get the tc number
	my $tc_num;
	if($probename_acc->{$probe_name} =~ /tc\|\w+/)
	{
		($tc_num) = $probename_acc->{$probe_name} =~ /tc\|(\w+)/;
	}
	print ANN join("\t", "$uid", "$row", "$col", "$accessions", "$location", "$control_type", "$probe_name", "$gene_name", "$systematic_name", "$description", $probename_desc->{$probe_name}, $probename_symbol->{$probe_name}, $probename_genename->{$probe_name}, $probename_acc->{$probe_name}, $probename_sys->{$probe_name}, $tc_num, $tc_tog->{$tc_num}, $tc_go_term->{$tc_num}, $tc_et_name->{$tc_num}, $tc_tentative_annotation->{$tc_num}, $tc_oligo->{$tc_num}, $uid ) . "\n";	
}

