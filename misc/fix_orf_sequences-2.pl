#!/usr/bin/perl
 
use Mbl;
use Bio::Seq;
  
my $mbl = Mbl::new(undef, 't_cruzi');
my $dbh = $mbl->dbh;

my $update = 0;


my $get_orfid_2 = $dbh->prepare("select distinct orfid from annotation where notes like '%codon_start:2%'");
my $get_orfid_3 = $dbh->prepare("select distinct orfid from annotation where notes like '%codon_start:3%'");
my $get_orf_info = $dbh->prepare("select orfid, sequence, start, stop, direction from orfs where orfid = ?");
my $update_orf = $dbh->prepare("update orfs set sequence = ?, start = ?, stop = ? where orfid = ?");


$get_orfid_2->execute;




while(my $row = $get_orfid_2->fetchrow_hashref)
{
	print $row->{orfid} . "\n";

	$get_orf_info->execute($row->{orfid});
	my $orfinfo = $get_orf_info->fetchrow_hashref or die("No orf found for " . $row->{orfid});

	# These are off by one on their start codon, so delete the first leter of the sequence in their sequence column, and move their actual start by one
	my $new_seq = substr($orfinfo->{sequence}, 1);
	my $new_start = $orfinfo->{start};
	my $new_stop = $orfinfo->{stop};
#	print $orfinfo->{sequence} . "\n\n" . $new_seq . "\n";

	# Now if the direction is +, add one to the start location
	if($orfinfo->{direction} eq "+")
	{
		$new_start++;
	} elsif($orfinfo->{direction} eq '-')
	{
		$new_stop--;
	} else
	{
		die("No direction!");
	}


	# ok, now we can update the records
	if($update)
	{
		$update_orf->execute($new_seq, $new_start, $new_stop, $row->{orfid});
	}

}

$get_orfid_3->execute;

while(my $row = $get_orfid_3->fetchrow_hashref)
{
	print $row->{orfid} . "\n";

	$get_orf_info->execute($row->{orfid});
	my $orfinfo = $get_orf_info->fetchrow_hashref or die("No orf found for " . $row->{orfid});

	# These are off by one on their start codon, so delete the first 2 leters of the sequence in their sequence column, and move their actual start by 2
	my $new_seq = substr($orfinfo->{sequence}, 2);
	my $new_start = $orfinfo->{start};
	my $new_stop = $orfinfo->{stop};
#	print $orfinfo->{sequence} . "\n\n" . $new_seq . "\n";

	# Now if the direction is +, add one to the start location
	if($orfinfo->{direction} eq "+")
	{
		$new_start++;
		$new_start++;
	} elsif($orfinfo->{direction} eq '-')
	{
		$new_stop--;
		$new_stop--;
	} else
	{
		die("No direction!");
	}


	# ok, now we can update the records
	if($update)
	{
		$update_orf->execute($new_seq, $new_start, $new_stop, $row->{orfid});
	}

}


