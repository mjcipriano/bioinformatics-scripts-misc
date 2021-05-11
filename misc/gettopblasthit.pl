#!/usr/bin/perl

use strict;
    
# Connect to the database
use DBI;
    
    
my $driver = "mysql";
my $hostname = "mib";
my $port = "3306";
my $user = "gid";
my $password = "NOPE";
my $database = "antonospora01";
    
my $dsn = "DBI:$driver:database=$database;host=$hostname;port=$port";
my $dbh = DBI->connect($dsn, $user, $password);

my $orfsh = $dbh->prepare("select orfid from orfs where delete_fg = 'N'");
 my $top_blasth = $dbh->prepare("select hit_name, description, evalue, accession_number, gi from blast_results where sequence_type_id = 2 AND db IN (2, 3) AND algorithm = 3 AND idname = ? AND evalue < 1e-2 order by evalue limit 1");

$orfsh->execute();

while(my $orf_row = $orfsh->fetchrow_hashref)
{
	$top_blasth->execute($orf_row->{orfid});
	if($top_blasth->rows > 0)
	{
		my $bh = $top_blasth->fetchrow_hashref;
		print join("\t", $orf_row->{orfid}, $bh->{evalue}, $bh->{accession_number}, $bh->{gi}, $bh->{hit_name}, $bh->{description}) . "\n"; 
	} else
	{
		print $orf_row->{orfid} . "\t" . "No Blast Hit\n";
	}
}
