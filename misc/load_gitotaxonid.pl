#!/usr/bin/perl


use DBD::mysql;
use RegTransBase;
use Getopt::Long;
use File::Temp qw /tempfile tempdir/;

use strict;



my $regobj = RegTransBase->new();
my $conn = $regobj->connection();

my $seqdb_dbh = $conn->get_seqdb_dbh();




my $gi_nt_file = '/bioware/taxon_tools/gi_taxid_nucl.dmp';
my $gi_aa_file = '/bioware/taxon_tools/gi_taxid_prot.dmp';
my $notlocal = 0;

my $options = GetOptions (      "nuc=s"=>\$gi_nt_file,
				"prot=s"=>\$gi_aa_file,
				"notlocal"=>\$notlocal
                );

my $tempdir = tempdir(CLEANUP=>1);

if(!defined($gi_nt_file))
{
	system("cd $tempdir;wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz;gunzip gi_taxid_nucl.dmp.gz");
	$gi_nt_file = "$tempdir/gi_taxid_nucl.dmp";

}

if(!defined($gi_aa_file))
{
	system("cd $tempdir;wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz;gunzip gi_taxid_prot.dmp.gz");
	$gi_nt_file = "$tempdir/gi_taxid_prot.dmp";
}

my $insert_sth = $seqdb_dbh->prepare("insert into gitotaxid (gi, ncbi_taxon_id) VALUES (?, ?)");
my $del_sth = $seqdb_dbh->prepare("delete from gitotaxid where gi = ?");
my $load_infile_sth = $seqdb_dbh->prepare("load data local infile ? into table gitotaxid");
my $del_all_sth = $seqdb_dbh->prepare("truncate table gitotaxid");




# Delete the table
$del_all_sth->execute();

# Now process both files

my $num_nt = process_file($gi_nt_file);

print "Loaded $num_nt nucleotide records!\n";

my $num_aa = process_file($gi_aa_file);
print "Loaded $num_aa protein records!\n";



sub process_file
{
	my $filename = shift;

	
	if($notlocal)
	{
		open(FILE, $filename) or die("FILE " . $filename . " NOT FOUND!\n");

		my $count = 0;
		while(<FILE>)
		{
			my $line = $_;
			chomp $line;
			my ($gi, $taxid) = split("\t", $line);
		
			# First delete this gi entry
			$del_sth->execute($gi);
			# Now insert this entry
			$insert_sth->execute($gi, $taxid);
			$count++;
		}
		# Done
		return $count;


	} else
	{
		$load_infile_sth->execute($filename);
		return 1;
	}

	return 0;
	
}



