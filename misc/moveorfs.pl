#!/usr/bin/perl
 

# This script will check all of the old orfs and try and
# Map them to a new assembly and insert the results into the orfs table
 
use Bio::SeqIO;
use Bio::Seq;
 
# Connect to the database
use DBI;
 
use strict;
 
my $driver = "mysql";
my $hostname = "mib";
my $port = "3306";
my $user = "gid";
my $password = "NOPE";
my $old_database = "giardia";
my $new_database = "giardia090403";
 
my $odsn = "DBI:$driver:database=$old_database;host=$hostname;port=$port";
my $odbh = DBI->connect($odsn, $user, $password);
 
my $odrh = DBI->install_driver("mysql");

my $ndsn = "DBI:$driver:database=$new_database;host=$hostname;port=$port";
my $ndbh = DBI->connect($ndsn, $user, $password);
my $ndrh = DBI->install_driver("mysql");

my $base_dir = "/var/bio/ARACHNE/";
my $working_dir = "gl080403";
my $run_dir = "run1";
my $root_dir = $base_dir . '/' . $working_dir . '/';

# New fasta bases file
my $fasta_bases_file = $root_dir . '/' . $run_dir . '/' . 'assembly.bases';

my $debug = 1;

# determine the max insert id for new orfs and set it

my $maxidh = $odbh->prepare("select max(orfid)+1 as next_orf from orfs");
$maxidh->execute();
my $maxrow_hash = $maxidh->fetchrow_hashref;
my $next_id = $maxrow_hash->{next_orf};
print "Next new orfid is $next_id \n";
my $insidh = $ndbh->prepare("set insert_id=" . $next_id);
$insidh->execute();

my $get_orf_query = "select ORFid, sequence, annotation, annotation_type, source, delete_fg, delete_reason, start, stop, direction, attributes, old_orf, orf_name from orfs order by orfid";
my $sth = $odbh->prepare($get_orf_query);
$sth->execute;

# iterate through all old orfs
while(my $this_row = $sth->fetchrow_hashref)
{
	if($debug)
	{
		print "Checking:\t" . $this_row->{ORFid} . "\t";
	}

	# now iterate through all sequences and check if this sequence is in one of the contigs

	my $sequences = Bio::SeqIO->new('-file'         => "<$fasta_bases_file",
        	                        '-format'       => "fasta");

	# Create a sequence object for the query string
	my $sequence_string_db = $this_row->{sequence};

	# Clean up the query string
	
        $sequence_string_db =~ s/[^ATGC]//ig;

	my $this_orf_sequence = Bio::Seq->new ( -display_id 	=> $this_row->{ORFid},
						-seq		=> $sequence_string_db);
	my $sequence_to_check = lc($this_orf_sequence->seq);
	my $sequence_to_check_rc = lc($this_orf_sequence->revcom()->seq);
	my @hit_array ;
	while (my $seq = $sequences->next_seq)
	{
		my $last_index = 0;
                my $this_index = -1;
		my $num_hits = 0;
		my $start = 0;
		my $stop = 0;
	
		if($debug)
		{
		#	print "\t" . $seq->display_id . "\n";
		}
		my $checked_first = 0;
		my $dir = "+";
		# check both this sequence and its reverse complement;
		while($checked_first < 2)
		{
			$last_index = 0;
			$this_index = -1;
			if($checked_first == 0)
			{
				$this_index = index(lc($seq->seq), $sequence_to_check);
				$dir = "+";
			} else
			{
				$this_index = index(lc($seq->seq), $sequence_to_check_rc);
				$dir = "-";
			}

			# while we have a hit
			while($this_index > -1)
			{
				$num_hits++;
				
				$start = $this_index +1; # index is zero based, so add one
				$stop = $start + length($this_orf_sequence->seq);

				if(length($this_orf_sequence->seq)%3 == 2)
				{
					if($dir eq '+')
					{
						$stop = $start + length($this_orf_sequence->seq) + 1;
						
					} else
					{
						$start = $start -1;					
					}
				}
				if($debug)
				{
	                                print  "ORF \t"
	                                . $this_row->{ORFid}
	                                . " matches \t" . $seq->display_id
	                                . " location:\t" . $start . "\t"
	                                .  $stop . "\t" . $dir . "\n";
				}
				

				push @hit_array, [$this_row->{ORFid}, $seq->display_id, $start, $stop, $dir];
				$last_index = $this_index;
	 			if($checked_first == 0)
	                        {
					$this_index = index(lc($seq->seq), $sequence_to_check, $last_index+1);
				} else
				{
					$this_index = index(lc($seq->seq), $sequence_to_check_rc, $last_index+1);

				}
			} # End while this db row sequence has checked a particular contig
			$checked_first++;
		} # End while this db row sequence reveresed has checked a particular contig
	} # End checking db row sequence to all contigs

	# Now process the @hit_array
	if((scalar @hit_array) == 0)
	{
		# Expire this orf with reason code, "disconnected"

		if($debug)
		{
		  print "ORF\t" . $this_row->{ORFid} . "\t Not found!\n";
		}
		my $insert_query = "insert into orfs 
			(ORFid, sequence, annotation, annotation_type,source,	
			 delete_fg, delete_reason, old_orf, attributes, orf_name) VALUES (" .
			"'" . $this_row->{ORFid} . "'," .
                        "'" . $this_orf_sequence->seq . "'," .
                        $ndbh->quote($this_row->{annotation}) . "," .
                        "'" . $this_row->{annotation_type} . "'," .
                        "'" . $this_row->{source} . "'," .
                        "'Y'," .
                        "'disconnected', 'Y'," .
                        "'" . $this_row->{attributes} . "'," .
                        "'" . $this_row->{orf_name} . "'" .
			")"; 
	 		my $insert_del_ref = $ndbh->prepare($insert_query);
                        $insert_del_ref->execute;

			
	} elsif((scalar @hit_array) > 1)
	{
		if($debug)
		{
			print "ORF\t" . $this_row->{ORFid} . "\t Multiple HITS!\n";
		}

		# we must then expire the ORFid and create new ORFid's
		foreach my $this_hit(@hit_array)
		{
			#print "old ORF\t" . $this_hit->[0] . "\t now in " . $this_hit->[1] . "\n";
	                my $insert_query = "insert into orfs
                        (ORFid, sequence, annotation, annotation_type,source, old_orf,
                         delete_fg, contig, start, stop, direction, attributes, orf_name) VALUES (" .
                        "'" . $next_id . "' ," .
                        "'" . $this_orf_sequence->seq . "'," .
                        $ndbh->quote($this_row->{annotation}) . "," .
                        "'" . $this_row->{annotation_type} . "'," .
                        "'" . $this_row->{source} . "'," .
                        "'" . 'Y' . "'," .
                        "'" . 'N' . "'," .
                        "'" . $this_hit->[1] . "'," .
                        "'" . $this_hit->[2] . "'," .
                        "'" . $this_hit->[3] . "'," .
                        "'" . $this_hit->[4] . "'," . 
                        "'" . $this_row->{attributes} . "'," .
                        "'" . $this_row->{orf_name} . "'" . 
			")";
			my $insert_ref = $ndbh->prepare($insert_query);
			$insert_ref->execute;
			$next_id++;
			my $last_insert_id = $ndbh->{'mysql_insertid'};
			my $insert_reassign_query = 
			"insert into orf_reassign (old_orf, new_orf) VALUES (" .
			"'" . $this_row->{ORFid} . "'," .
			"'" . $last_insert_id . "')";
			my $insert_reassign = $ndbh->prepare($insert_reassign_query);
			$insert_reassign->execute;
				
		}
	} else
	{
		# it is the only hit, so we can just use this orf id
		my $delete_flag = "";
		my $delete_reason = $this_row->{delete_fg};
		if($this_row->{delete_fg} eq "no")
		{
			$delete_flag = "N";
			$delete_reason = "";
		} else
		{
			$delete_flag = "Y";
		}
		my $insert_query = "insert into orfs
                        (ORFid, sequence, annotation, annotation_type,source,
                         delete_fg, delete_reason, contig, start, stop, direction, attributes, orf_name) VALUES (" .
                        "'" . $this_row->{ORFid} . "'," .
                        "'" . $this_orf_sequence->seq . "'," .
                        $ndbh->quote($this_row->{annotation}) . "," .
                        "'" . $this_row->{annotation_type} . "'," .
                        "'" . $this_row->{source} . "'," .
                        "'" . $delete_flag . "'," .
                        "'" . $delete_reason . "'," .
                        "'" . $hit_array[0][1] . "'," .
                        "'" . $hit_array[0][2] . "'," .
                        "'" . $hit_array[0][3] . "'," .
                        "'" . $hit_array[0][4] . "'," . 
                        $ndbh->quote($this_row->{attributes})  . "," .
                        "'" . $this_row->{orf_name} . "'" .
			")"; 
			my $insert_ref = $ndbh->prepare($insert_query);
			$insert_ref->execute;
	}

	
} # end checking all db rows




