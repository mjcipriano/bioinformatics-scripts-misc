#!/usr/bin/perl
 
use strict;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Root::IO;
 
my $organism = $ARGV[0];

my $rfam_file = $ARGV[1];
my $fasta_input = $ARGV[2];

 
my $cmsearch_bin = '/xraid/bioware/linux/infernal-0.55/src/cmsearch'; 

# Make a temporary input file

#system("mkdir temp");

#system("$cmsearch_bin $rfam_file $fasta_input > $outfile");



# Open the Rfam directory

my $rfam_dir = "/blastdb/rfam";
my $seed_file = "Rfam.seed";
my $cwd = `pwd`;
chomp($cwd);

my $infile = $ARGV[0];
if(!$infile)
{
	die("No Input file specified\n");
}

open(RFAM, $rfam_dir . "/" . $seed_file);

my $rfam_hash;
my $rfam_id_hash;
my $rfam_desc_hash;
my $rfam_wsize_hash;
my $rfam_comment_hash;

my $temp_ac;
my $temp_id;
my $temp_desc;
my $temp_comment;
my $temp_wsize;
my $last_id = '0';

while(<RFAM>)
{
	if($_ =~ /^\#=GF\ AC/)
	{
		($temp_ac) = $_ =~ /^\#=GF\ AC\s+(.+)$/;
		$rfam_hash->{$temp_ac} = 1;
	} elsif($_ =~ /^\#=GF\ ID/)
	{
		($temp_id) = $_ =~ /^\#=GF\ ID\s+(.+)$/;
		$rfam_id_hash->{$temp_ac} = $temp_id;
	} elsif($_ =~ /^\#=GF\ DE/)
	{
		($temp_desc) = $_ =~ /^\#=GF\ DE\s+(.+)$/;
		$rfam_desc_hash->{$temp_ac} = $temp_desc;
	} elsif($_ =~ /^\#=GF\ BM\s+cmsearch/)
	{
		($temp_wsize) = $_ =~ /^\#=GF\ BM\s+cmsearch.+\-W\ (\d+)/;
		$rfam_wsize_hash->{$temp_ac} = $temp_wsize;
	}
	
}
while(my ($key, $val) = each(%$rfam_hash))
{
#	print join("\t", $key, $rfam_id_hash->{$key}, $rfam_desc_hash->{$key}, $rfam_wsize_hash->{$key}) . "\n";
 
	my $shell_name = $cwd . '/' . $rfam_id_hash->{$key} . '.sh';
 	open(SUBMIT, ">", $shell_name);
	print SUBMIT '#!/bin/csh
#$ -j y
#$ -o ' . $cwd . '/rfam.out
#$ -N ' . $rfam_id_hash->{$key} . '
#$ -cwd
 
cd ' . $cwd . ';
' . $cmsearch_bin . ' ' . $rfam_dir . "/" . $key . '.cm ' . $infile . ' > ' . $rfam_id_hash->{$key}  . '.rfam;
rfam_parse.pl ' . $rfam_id_hash->{$key} . ' ' . $rfam_id_hash->{$key}  . '.rfam;
echo "DONE";';

	close(SUBMIT);
 
	system("chmod 755 $shell_name");
#	system("cd $cwd;qsub $shell_name");
}

