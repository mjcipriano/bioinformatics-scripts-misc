#!/usr/bin/perl


use Bio::SeqIO;
use Bio::Seq;
#use Bio::Factory::EMBOSS;

# Connect to the database
use DBI;

use strict;

#my $organism = 'giardiabac03';
#my $organism = 'giardianobac03';
#my $organism = 'giardia';
my $organism  = 'giardia090403';
#my $organism = 'nosema';
#my $organism = 't_brucei';

my $driver = "mysql";
my $hostname = "mib";
my $port = "3306";
my $user = "gid";
my $password = "NOPE";
my $database = $organism;

my $dsn = "DBI:$driver:database=$database;host=$hostname;port=$port";

# File/script locations


my $db_database_name = $organism . "db";
my $db_supercontig_database_name = $organism . "sc";
my $db_supercontig_read_database_name = $organism . "screads";

#my $base_dir = "/var/bio/ARACHNE";
#my $working_dir = "gl063003";
#my $run_dir = "run1";

#my $base_dir = "/xraid/ARACHNE/Arachne_data";
#my $working_dir = "gl080403";
#my $run_dir = "run1";

#my $base_dir = "/xraid/ARACHNE/Arachne_data";
#my $working_dir = "gl080403";
#my $run_dir = "run2";

#my $base_dir = "/xraid/ARACHNE/Arachne_data";
#my $working_dir = "nl080703";
#my $run_dir = "run1";

my $base_dir = "/xraid/licor2/giardia";
my $working_dir = "gl090403";
my $run_dir = "run1";


my $root_dir = $base_dir . '/' . $working_dir . '/';
my $xml_file = $root_dir . "/traceinfo/reads.xml";
my $assembly_unplaced_file = $root_dir . '/' . $run_dir . '/' . 'assembly.unplaced';
my $assembly_links_file = $root_dir . '/' . $run_dir . '/' . 'assembly.links';
my $assembly_reads_file = $root_dir . '/' . $run_dir . '/' . 'assembly.reads';
my $fasta_bases_file = $root_dir . '/' . $run_dir . '/' . 'assembly.bases';
my $fasta_reads_bases_file =  $root_dir . '/' . 'reads.fasta';
my $fasta_reads_quality_file =  $root_dir . '/' . 'reads.qual';

my $gff_output_file = "gff/" . $organism . ".gff";
my $gff_supercontig_output_file = "gff/" . $organism . "_supercontig.gff";
my $gff_supercontig_read_output_file = "gff/" . $organism . "_supercontig_read.gff";
my $gff_supercontig_orf_output_file = "gff/" . $organism . "_orf_supercontig.gff";
my $gff_orf_output_file = "gff/" . $organism . "_orf.gff";
my $gff_supercontig_sage_output_file =  "gff/" . $organism . "_supercontig_sage.gff";
my $gff_sage_output_file =  "gff/" . $organism . "_sage.gff";
my $gff_repeat_file = "gff/" . $organism . '_repeat.gff';
my $gff_supercontig_repeat_file = "gff/" . $organism . 'supercontig__repeat.gff';



my $temp_file_name = "gff/" . $organism . "_temp.gff";

my $orf_input_file = "gff/" . $organism . "_orfs.gff";
my $database_bulk_loader = "../bulk_load_gff_nowarn.pl";
my $database_incremental_loader = "../load_gff.pl";
my $ace_dir = $root_dir . '/' . $run_dir . '/' . 'acefiles';
my $sage_dir = '/xraid/licor2/gsage/';
my $sage_file = 'sageanalysis/results/091203.alllib.tags';
my $supercontig_fasta_file = 'fasta/' . $organism . "_supercontig.fas";


my $blast_bin_dir = '/var/bio/blast/bin/';
my $blast_db_dir = '/var/bio/blast/data/';
my $repeatFinder_bin = '/usr/local/bin/repeatFinder';

# ORGANISM.fasta file should exist in the tables directory which contains known genes to be used to train glimmer and Codon Preference with
my $glimmer_bin_dir = 'glimmer/';




my $debug = 0;
# Turn off/on parts of this script

## Putting this in here so I do not make a mistake
my $drop_schema				= 0;


my $create_schema			= 1;
#my $drop_schema				= 0;
my $parse_xml      			= 1;
my $parse_unplaced 			= 1;
my $parse_links    			= 1;
my $parse_reads    			= 1;
my $parse_reads_bases 			= 1;
my $parse_reads_qual 			= 1;
my $create_modified_fasta		= 1; # will create a supercontig fasta file to use against the supercontig view
my $move_orfs_from_old			= 1; # calls moveorfs.pl script (edit that script first)
my $find_orfs_glimmer			= 1;
my $delete_invalid_orfs			= 1;
my $parse_ace	   			= 0;
my $create_supercontig_orf_from_input 	= 0; # will read in a tab delimited file and create a gff file, not used any more
my $run_orf_tests 			= 1; # will run tests on all orfs to give P/F or scores to determine if they are real orfs
my $create_sage_from_file 		= 1; # to import sage data from a tab delimited file
my $map_sage_to_db			= 1;
my $map_sage_to_orf 			= 1; # first test run, do not use anymore
my $map_sage_to_orf_secondary 		= 1; # this will dump data into a temp table to be used for next part
my $map_sage_to_orf_tert		= 1; # this will insert into orftosage for determined tag->orf mappings
my $find_repeat				= 1;
my $create_supercontig_orf_from_db      = 1;
my $create_sage_from_db 		= 1; # this will create sage.gff files
my $output_gff     			= 1;
my $create_supercontig_map 		= 1;
my $create_supercontig_read_map 	= 1;

my $create_blast_db			= 1;
my $load_db        			= 1; # this will load the database with this organisms gff files

my $minimum_gap_fg 			= 1;
my $minimum_gap_length			= 100;

my $transcript_tail = 15;

# include package
use XML::DOM;

if($drop_schema)
{
	print "Start: Dropping Schema\n";
        # I Really don't want to drop this one ##SAFETY##
        if($organism eq "giardia")
        {
        	exit;
        }

        system("echo \"drop database $organism\" | mysql");
        system("echo \"drop database $db_database_name\" | mysql");
	system("echo \"drop database $db_supercontig_database_name\" | mysql");
	system("echo \"drop database $db_supercontig_read_database_name\" | mysql");
	print "End: Dropping Schema";
}

if($create_schema)
{
	print "Start: Creating Schema\n";
	system("echo \"create database $organism\" | mysql");
	system("echo \"create database $db_database_name\" | mysql");
	system("echo \"create database $db_supercontig_database_name\" | mysql");
	system("echo \"create database $db_supercontig_read_database_name\" | mysql");


	system("mysql -D $organism < sql/create_script.sql");

	system("echo \"grant select on $db_database_name.* to nobody\@localhost\" | mysql");
	system("echo \"grant select on $db_supercontig_database_name.* to nobody\@localhost\" | mysql");
	system("echo \"grant select on $db_supercontig_read_database_name.* to nobody\@localhost\" | mysql");

	system("echo \"grant all privileges on $organism.* to $user\@'%' identified by '$password'\" | mysql");
	system("echo \"grant all privileges on $db_database_name.* to $user\@'%' identified by '$password'\" | mysql");
	system("echo \"grant all privileges on $db_supercontig_database_name.* to $user\@'%' identified by '$password'\" | mysql");
	system("echo \"grant all privileges on $db_supercontig_read_database_name.* to $user\@'%' identified by '$password'\" | mysql");
	system("echo \"flush privileges\" | mysql");
	print "End: Creating Schema\n";
}




# NOW We can connect to the database we just created

my $dbh = DBI->connect($dsn, $user, $password);
                                                                                                                                                                                    
my $drh = DBI->install_driver("mysql");





if($parse_xml)
{
	
	print "Start: Parse reads.xml file\n";

	# instantiate parser
	my $xp = new XML::DOM::Parser();
	
	# parse and create tree
	my $doc = $xp->parsefile($xml_file) or die ("Can not open reads.xml file!");
	
	# get root node (trace_volume) 
	my $root = $doc->getDocumentElement();
	
	
	# get children (trace)
	my @my_reads = $root->getChildNodes();
	
	# iterate through trace_volume list
	foreach my $node (@my_reads)
	{
		# if element node
		if ($node->getNodeType() == 1)
		{
	                my $trace_name = '';
	                my $center_name = '';
	                my $plate_id = '';
	                my $well_id = '';
	                my $template_id = '';
	                my $library_id = '';
	                my $trace_end = '';
	                my $trace_direction = '';
	
			my @children = $node->getChildNodes();
	
			# iterate through child nodes
			foreach my $item (@children)
			{
				# check element name
				if (lc($item->getNodeName) eq "trace_name")
				{
					$trace_name = $item->getFirstChild()->getData;
				}elsif (lc($item->getNodeName) eq "center_name")
				{
					$center_name =  $item->getFirstChild()->getData;
				}elsif (lc($item->getNodeName) eq "plate_id")
				{
					$plate_id = $item->getFirstChild()->getData;
				}elsif (lc($item->getNodeName) eq "well_id")
				{
	                                $well_id = $item->getFirstChild()->getData;
				}elsif (lc($item->getNodeName) eq "template_id")
	                        {
	                                $template_id = $item->getFirstChild()->getData;
	                        }elsif (lc($item->getNodeName) eq "library_id")
	                        {
	                                $library_id = $item->getFirstChild()->getData;
	                        }elsif (lc($item->getNodeName) eq "trace_end")
	                        {
	                                $trace_end = $item->getFirstChild()->getData;
	                        }elsif (lc($item->getNodeName) eq "trace_direction")
	                        {
	                                $trace_direction = $item->getFirstChild()->getData;
	                        }
	
			}
	        my $query = 'insert into reads (
	                    read_name,
	                    center_name,
	                    plate_id,
	                    well_id,
	                    template_id,
	                    library_id,
	                    trace_end,
	                    trace_direction
	                  ) VALUES (' .
	               "'" . $trace_name  . "'" . ', ' . 
	               "'" . $center_name . "'" . ', ' . 
	               "'" . $plate_id    . "'"  . ', ' . 
	               "'" . $well_id     . "'" . ', ' .
	               "'" . $template_id . "'" . ', ' .
	               "'" . $library_id  . "'" . ', ' .
	               "'" . $trace_end   . "'" . ', ' .
	               "'" . $trace_direction . "'" . 
	               ')';
	         my $sth = $dbh->prepare($query);
	         $sth->execute;
	         $sth->finish;
		print ".";
		}
	}
        print "End: Parse reads.xml file\n";

}	
#######DONE WITH READS.XML FILE###############

#######START ASSEMBLY.UNPLACED FILE#############

if($parse_unplaced)
{
	print "Start: Parse assembly.unplaced file\n";

	open(ASSM, "$assembly_unplaced_file") or die ("Can not open assembly.unplaced");
	while (<ASSM>)
	{
	    my $line = $_;
	    my @this_line = split(" ", $line);
	    my $update_status_query = "UPDATE reads 
	                            SET status = '$this_line[1]' 
	                            WHERE read_name = '$this_line[0]'";
	    my $sth = $dbh->prepare($update_status_query);
	    $sth->execute;
	    $sth->finish;
	    print ".";
	
	}
	close(ASSM);
        print "End: Parse assembly.unplaced file\n";
}
#######DONE ASSEMBLY.UNPLACED FILE#############

#######START ASSEMBLY.LINKS FILE#############

if($parse_links)
{
	print "Start: Parse assembly.links file\n";
	open(LINKS, "$assembly_links_file") or die ("Can not open assembly.links");
	
	#super_id
	#num_bases_in_super
	#num_contigs_in_super
	#ordinal_num_of_contig_in_supercontig_id
	#contig_number
	#length_of_contig
	#estimated_gap_before_contig
	#estimated_gap_after_contig
	
	while (<LINKS>)
	{
	    my $line = $_;
	    my @this_line = split(" ", $line);
	
	    my $update_links_query = "insert into links
	                            (super_id, bases_in_super, contigs_in_super, ordinal_number,
	                            contig_number, contig_length, gap_before_contig, gap_after_contig)
	                            VALUES ('$this_line[0]', 
	                                    '$this_line[1]',
	                                    '$this_line[2]',
	                                    '$this_line[3]',
	                                    '$this_line[4]',
	                                    '$this_line[5]',
	                                    '$this_line[6]',
	                                    '$this_line[7]')";
	
	    if($this_line[0] ne '#super_id')
	    {	
		    my $sth = $dbh->prepare($update_links_query);
		    $sth->execute;
		    $sth->finish;
		    print ".";
	    }
	
	}
	close(LINKS);
	print "Completed parsing assembly.links file\n";
                                                                                                                                                                                   
	# Now determine where contigs are within the supercontig with or without minimum gap length 
                                                                                                                                                                                    
        my $query = '
                select distinct
                super_id
                FROM
                links';
        my $links_result = $dbh->prepare($query);
        $links_result->execute;
                                                                                                                                                                                    
                                                                                                                                                                                    
        while(my $links_array = $links_result->fetchrow_hashref)
        {
                                                                                                                                                                                    
                # We need to find out the total modified supercontig length And update where the contig starts in the supercontig
                my $start_super_base = 1;
                my $running_end = 0;
                my $contq = "select contig_number, contig_length, ordinal_number, gap_before_contig from links where super_id = '" . $links_array->{super_id} . "' ORDER BY ordinal_number";
                                                                                                                                                                                    
                my $conth = $dbh->prepare($contq);
                $conth->execute();
                while(my $this_contig = $conth->fetchrow_hashref)
                {
                        my $this_start = 0;
                        if($this_contig->{gap_before_contig} < 0)
                        {
                                $this_start = $running_end + $minimum_gap_length;
                        } else
                        {
                                $this_start = $running_end + $this_contig->{gap_before_contig};
                        }
                        my $updq = "update links set
                                modified_contig_start_base = '" . $this_start . "'
                                where contig_number = '" . $this_contig->{contig_number} . "'";
                        my $updh = $dbh->prepare($updq);
                        $updh->execute();
                        $running_end = $this_start + $this_contig->{contig_length};
                }
                my $updq = "update links set
                                modified_bases_in_super = '" . $running_end . "'
                                where super_id = '" . $links_array->{super_id} . "'";
                        my $updh = $dbh->prepare($updq);
                        $updh->execute();
        }






	# Now without minimum gap length
	my $query = '
	        select contig_number,
	        super_id,
	        bases_in_super,
	        contigs_in_super,
	        ordinal_number,
	        contig_length,
	        gap_before_contig,
	        gap_after_contig,
	        contig_start_super_base,
	        modified_contig_start_base,
	        modified_bases_in_super
	        FROM
	        links ORDER BY super_id, ordinal_number';
	my $links_result = $dbh->prepare($query);
	$links_result->execute;
                                                                                                                                                                                    
	my $last_super_id = '';
	my $super_running_total = 0;



	while(my $links_array = $links_result->fetchrow_hashref)
	{
                                                                                                                                                                                    
		if($links_array->{super_id} != $last_super_id)
		{
			$super_running_total = 0;
		}
                                                                                                                                                                                    
		my $start_val = $super_running_total + $links_array->{gap_before_contig};
		my $end_val = $start_val + $links_array->{contig_length};
		$super_running_total = $end_val;
		$last_super_id = $links_array->{super_id};

		my $update_base_query = "UPDATE links set contig_start_super_base = '" . $start_val . "'
	                    WHERE contig_number = '" . $links_array->{contig_number} ."'";
		my $update_base_result = $dbh->prepare($update_base_query);
		$update_base_result->execute;
	}
        print "End: Parse assembly.links file\n";

}

#######END ASSEMBLY.LINKS FILE#############

#######START ASSEMBLY.READS FILE#############
if($parse_reads)
{
	print "Start: Parse assembly.reads file\n";
	open(READS, "$assembly_reads_file") or die ("Can not open assembly.reads");

	#read_name
	#read_status
	#read_len_untrim
	#first_base_of_trim
	#read_len_trim
	#contig_number
	#contig_length
	#trim_read_in_contig_start
	#trim_read_in_contig_stop
	#orientation
	#read_pair_name
	#read_pair_status
	#read_pair_contig_number
	#observed_insert_size
	#given_insert_size
	#given_insert_std_dev
	#observed_inserted_deviation
	
	
	while (<READS>)
	{
	    my $line = $_;
	    my @this_line = split("\t", $line);
            foreach my $this_var (@this_line) 
	    {
		if($this_var eq "")
		{
		  $this_var = "NULL";
		} else
		{
		  $this_var = "'" . $this_var . "'";
		}
	    }
	    my $update_reads_query = "insert into reads_assembly 
					(read_name,
					read_status,
					read_len_untrim,
					first_base_of_trim,
					read_len_trim,
					contig_number,
					contig_length,
					trim_read_in_contig_start,
					trim_read_in_contig_stop,
					orientation,
					read_pair_name,
					read_pair_status,
					read_pair_contig_number,
					observed_insert_size,
					given_insert_size,
					given_insert_std_dev,
					observed_inserted_deviation)
	                            VALUES ($this_line[0],
	                                    $this_line[1],
	                                    $this_line[2],
	                                    $this_line[3],
	                                    $this_line[4],
	                                    $this_line[5],
	                                    $this_line[6],
	                                    $this_line[7],
	                                    $this_line[8],
	                                    $this_line[9],
	                                    $this_line[10],
	                                    $this_line[11],
	                                    $this_line[12],
	                                    $this_line[13],
	                                    $this_line[14],
	                                    $this_line[15],
	                                    $this_line[16])";
	                                                                                                 
	    my $sth = $dbh->prepare($update_reads_query);
	    $sth->execute;
	    $sth->finish;
	    print ".";

	}
	close(READS);
        print "End: Parse assembly.reads file\n";

}

#######END ASSEMBLY.READS FILE#############


#######START READS.BASES##############
if($parse_reads_bases)
{

        print "Start: Parse assembly.bases file\n";

	my $in  = Bio::SeqIO->new('-file' => "$fasta_reads_bases_file",
                         '-format' => 'Fasta');

	while ( my $seq = $in->next_seq() )
	{
		my $insert_query = "insert into reads_bases (read_name, bases) VALUES ( '" 
			. $seq->id . "', '" .  $seq->seq . "')";
 
		my $sth = $dbh->prepare($insert_query);
		$sth->execute;
	}
        print "End: Parse assembly.bases file\n";

}


#######END READS.BASES FILE#############


#######START CREATE SUPERCONTIG FASTA##########

if($create_modified_fasta)
{

        print "Start: Create modified fasta file\n";


	open(SUPERFASTA, ">$supercontig_fasta_file");

	my $superquery = 'select distinct super_id from links order by super_id';
	my $superh = $dbh->prepare($superquery);
	$superh->execute();

	while(my $supercon = $superh->fetchrow_hashref)
	{
	
		print SUPERFASTA '>supercontig_' . $supercon->{super_id} . "\n";

		my $query = "
	                select contig_number,
	                super_id,
	                bases_in_super,
	                contigs_in_super,
	                ordinal_number,
	                contig_length,
	                gap_before_contig,
	                gap_after_contig,
	                contig_start_super_base,
	                modified_contig_start_base,
	                modified_bases_in_super
	                FROM
	                links WHERE
			super_id = '" . $supercon->{super_id} . "'
			ORDER BY ordinal_number";
		my $linksh = $dbh->prepare($query);
		$linksh->execute();

		while(my $this_contig = $linksh->fetchrow_hashref)
		{
			my $sequences = Bio::SeqIO->new('-file'         => "<$fasta_bases_file",
                                                '-format'       => "fasta");
			my $contig_found = 0;
			my $contig_bases = '';
			while( (my $fasta_contig = $sequences->next_seq) && (!$contig_found) )
			{
				if($fasta_contig->display_id eq ('contig_' . $this_contig->{contig_number}) )
				{
					$contig_found = 1;
					$contig_bases = $fasta_contig->seq;
				}
			
			}
			my $num_placeholder = 0;
			if($this_contig->{gap_before_contig} < 0)
			{
				$num_placeholder = $minimum_gap_length;
			}else
			{
				$num_placeholder = $this_contig->{gap_before_contig};
			}
			for(my $i = 0;$i < $num_placeholder;$i++)
			{
				print SUPERFASTA 'N';
			}
			print SUPERFASTA $contig_bases;
			print SUPERFASTA "\n";
			
			
		}
	}

        print "End: Create modified fasta file\n";


}

#######END CREATE SUPERCONTIG FASTA##########

if($move_orfs_from_old)
{
	system('moveorfs.pl');
}

#######START FIND ORFS VIA GLIMMER############

if($find_orfs_glimmer)
{

        print "Start: Discover ORFs via GLIMMER2\n";

	system($glimmer_bin_dir . 'build-icm -f < tables/' . $organism . '.fasta > tables/' . $organism . '_glimmer.bin');
        my $sequences = Bio::SeqIO->new('-file'         => "<$fasta_bases_file",
                                        '-format'       => "fasta");

	my $tmp_file = 'temp/' . $organism . '_tmp.fas'; 
        while (my $seq = $sequences->next_seq)
        {
		open(TMPFILE, ">", $tmp_file);
	
		print TMPFILE ">" . $seq->display_id . "\n" . $seq->seq();
		close(TMPFILE);
		system($glimmer_bin_dir . '/glimmer2 ' . $tmp_file . ' tables/' . $organism . '_glimmer.bin > temp/' . $organism . '_glimmer.out');

		open (ORFILE, "<", 'temp/' . $organism . '_glimmer.out');
		while (defined (my $line = <ORFILE>)) 
		{
			if ($line =~/\[/) 
			{
 				$line=~/\s+\d+\s+(\d+)\s+(\d+)\s+\[([-+])\d\s.+/;
				my $start = $1;
				my $stop = $2;
				my $dir = $3;
				my $orf_seq = '';
				if($start > $stop)
				{
					my $tmpend = $stop;
					$stop = $start;
					$start = $tmpend;
				}
				if($dir eq "+" && $seq->length() >= $stop+3)
				{
					$stop = $stop+3;
				}

				if($dir eq "-" && $start > 3)
				{
					$start = $start - 3;
				}
				# print $seq->display_id . "\tGLIMMER\tORF\t$1\t$2\t.\t$3\t0\n";
				if($dir eq "-")
				{
                                        $orf_seq = $seq->trunc($start, $stop)->revcom()->seq() . "\n";
				} else
				{
                                        $orf_seq = $seq->subseq($start, $stop) . "\n";
				}
				# First see if we already have this orf in the database 
		                my $chk_query = "select orfid from orfs 
                                        WHERE
                                        contig = '"     . $seq->display_id   . "' AND
                                        start <= '"     . $start    . "' AND
                                        stop >= '"      . $stop     . "' AND
                                        direction = '"  . $dir . "' AND
                                        (start % 3) = '"  . ($start%3). "'";
				my $chk_h = $dbh->prepare($chk_query);
				$chk_h->execute();
				if($chk_h->rows() > 0)
				{
					# do nothing
				} elsif( ($stop - $start) < 200)
				{
					# do nothing
				} elsif( ($stop - $start) > 20000)
				{
					# do nothing
				} else
				{

					my $ins_query = "insert into orfs (sequence, source, delete_fg, contig, start, stop, direction, old_orf) VALUES ( " . 
						"'" . $orf_seq . "'," .
	                                        "'" . 'GLIMMER2' . "'," .
	                                        "'" . 'N' . "'," .
	                                        "'" . $seq->display_id . "'," .
	                                        "'" . $start . "'," .
	                                        "'" . $stop . "'," .
	                                        "'" . $dir . "'," .
	                                        "'" . 'N' . "')" ;
					my $insh = $dbh->prepare($ins_query);
					$insh->execute();
				}
			}
		}
	
	}
        print "End: Discover ORFs via GLIMMER2\n";


}
#######END FIND ORFS VIA GLIMMER############


#######START FIND INVALID ORFS############
if($delete_invalid_orfs)
{


	# First delete all called orfs that are too large (> 20000)

	my $delh = $dbh->prepare("update orfs set delete_fg = 'Y', delete_reason='invalid size' where stop-start > 20000 AND delete_fg = 'N'");
	$delh->execute();

	# Now check for and set delete_fg for orfs that are duplicates of other orfs
	# or are entierly within another orf. If one of them is an old orf
	# keep that one
	 
	# Now we will mark all duplicate orfs. We will always keep the lower numbered orf
	# as the correct orf, and we will insert into the orf_reassign table to keep track of
	# these orfs
	 
	my $last = 1;
	my $last_orf_id = 0;
	 
	while($last)
	{
	        # pick one orfs higher then the last orf
	        my $orf_query = "SELECT orfid,
	                        sequence,
	                        annotation,
	                        annotation_type,
	                        source,
	                        delete_fg ,
	                        delete_reason,
	                        contig,
	                        start,
	                        stop,
	                        direction
	                        FROM orfs
	                        where delete_fg = 'N' AND
	                        orfid > " . $last_orf_id . " order by orfid LIMIT 1";
	        my $sth = $dbh->prepare($orf_query);
	        $sth->execute;
	 
	        my $this_row = $sth->fetchrow_hashref;
	 
	        if($this_row)
	        {
	                $last_orf_id = $this_row->{orfid};
	 
	                my $find_query = "select orfid from orfs WHERE
	                                contig = '" . $this_row->{contig}       . "' AND
	                                start  = '" . $this_row->{start}        . "' AND
	                                stop  = '" . $this_row->{stop}          . "' AND
	                                direction  = '" . $this_row->{direction}. "' AND
	                                orfid != '" . $this_row->{orfid} . "' AND delete_fg = 'N'";
	                my $findh = $dbh->prepare($find_query);
	                $findh->execute;
	                print "ORF\t" . $this_row->{orfid} . "\t duplicates: " . $findh->rows . "\n";
	 
	                while(my $find_row = $findh->fetchrow_hashref)
	                {
	                        my $insert_query = "INSERT into orf_reassign (old_orf, new_orf) VALUES
	                                                ('" . $find_row->{orfid} . "', '" .
	                                                $this_row->{orfid} . "')";
	                        my $update_query = "UPDATE orfs set delete_fg = 'Y',
	                                                delete_reason = 'duplicate'
	                                                where orfid = '" . $find_row->{orfid} . "' AND delete_fg = 'N'";
	                        my $insh = $dbh->prepare($insert_query);
	                        $insh->execute;
	                        my $updh = $dbh->prepare($update_query);
	                        $updh->execute;
	                }
	 
	        } else
	        {
	                $last = 0;
	        }
	 
	}
	 
	 
	# we will select all orfs which are not deleted that are greater then the last orf
	# we looked at and then grab one record and update all orfs that follow the dulicate
	# criteria
	 
	my $last = 1;
	my $last_orf_id = 0 ;
	 
	while($last)
	{
	        my $orf_query = "SELECT orfid,
	                        sequence,
	                        annotation,
	                        annotation_type,
	                        source,
	                        delete_fg ,
	                        delete_reason,
	                        contig,
	                        start,
	                        stop,
	                        direction
	                        FROM orfs
	                        where delete_fg = 'N' AND
	                        orfid > '" . $last_orf_id . "' order by orfid LIMIT 1";
	        #print $orf_query . "\n\n";
	        my $sth = $dbh->prepare($orf_query);
	        $sth->execute;
	 
	        my $this_row = $sth->fetchrow_hashref;
	 
	        if($this_row)
	        {
	                $last_orf_id = $this_row->{orfid};
	                my $update_query = "UPDATE orfs set delete_fg = 'Y',
	                                        delete_reason = 'within other'
	                                        WHERE
        	                                contig = '"     . $this_row->{contig}   . "' AND
	                                        start >= '"     . $this_row->{start}    . "' AND
        	                                stop <= '"      . $this_row->{stop}     . "' AND
	                                        direction = '"  . $this_row->{direction}. "' AND
	                                        (start % 3) = '"  . ($this_row->{start}%3). "' AND
        	                                delete_fg = 'N' AND
	                                        orfid != '"     . $this_row->{orfid}    . "'";
	                my $updh = $dbh->prepare($update_query);
	                $updh->execute;
	                #print $update_query;
	                #print "\n\n\n";
	                print "ORF\t" . $this_row->{orfid} . "\t within: " . $updh->rows;
	                print "\n";
	        }else
	        {
	                $last = 0;
	        }
	 
	}# end if ($last)
	
	
}
                           
######END REMOVE INVALID ORFS##############


if($parse_ace)
{

 ## THIS IS UNFINISHED!!!!!!!!####################

  my $ace_file  = $ace_dir . '/' . 'supercontig.ace.1';

  open(ACE, "$ace_file") or die ("Can not open $ace_file");
  my @whole_acefile = <ACE>;

  my $this_file = join("\n", @whole_acefile);
  my @contigs = split("\nCO ", $this_file);

  #$contigs[0] is garbage 
  #$contigs[1-?] is each contig in this ace file
  
  my $num_contigs = scalar @contigs;
  print "Parsing $num_contigs Contigs in file $ace_file\n";
  for(my $i = 1;$i < $num_contigs;$i++)
  {
    
  }


}




if($create_supercontig_orf_from_input)
{
        open(ORF, "$orf_input_file") or die ("Can not open $orf_input_file");
	open(ORFSUPER, '>', "$gff_supercontig_orf_output_file") or die ("Can not open $gff_supercontig_orf_output_file");
        while (<ORF>)
        {
            my $line = $_;
            my @this_line = split("\t", $line);
	    my $contig_number = $this_line[0];
	    my $source = $this_line[1];
            my $feature = $this_line[2];
            my $start = $this_line[3];
            my $end = $this_line[4];
	    my $score = $this_line[5];
	    my $strand = $this_line[6];
	    my $frame = $this_line[7];
	    my $attr = $this_line[8];

	    my @contig_number_arr = split("_", $contig_number);
	    $contig_number = $contig_number_arr[1];
            my $get_contig_pos_query = "select contig_number, super_id, contig_start_super_base, contig_length, modified_contig_start_base from links where contig_number = '" . $contig_number . "'";
            my $sth = $dbh->prepare($get_contig_pos_query);
            $sth->execute;
	    my $links_array = $sth->fetchrow_hashref;
	    my $new_start = 0;
	    my $new_end = 0;

	    if($minimum_gap_fg)
	    {
	      $new_start = $start + $links_array->{modified_contig_start_base};
              $new_end = $end +  $links_array->{modified_contig_start_base};
  	    } else
	    {
	      $new_start = $start + $links_array->{contig_start_super_base};
	      $new_end = $end +  $links_array->{contig_start_super_base};
	    }
	    print ORFSUPER "supercontig_" . $links_array->{super_id} . "\t" .
			   $source 	. "\t" .
		   	   $feature 	. "\t" .
			   $new_start	. "\t" .
			   $new_end	. "\t" .
			   $score	. "\t" .
			   $strand	. "\t" .
			   $frame	. "\t" .
			   $attr	. "\n";
                                                                                                 
        }




}

if($create_sage_from_file)
{
        print "Start: Parsing Sage Input File\n";

  open(SAGEFILE, "$sage_dir" . "$sage_file") or die("Can not open sage file");

  my $first_line = 1;
  my @library_names;

  # Delete from sage_tags and sage_results;
  $dbh->prepare('delete from sage_tags')->execute();
  $dbh->prepare('delete from sage_results')->execute();

  while(<SAGEFILE>)
  {
    chop;
    my $line = $_;
    my @sage_array = split("\t", $line);

    # Get library names from first line of array
    if($first_line)
    {
      my $count = 0;
      foreach my $lib_name (@sage_array)
      {
        if($count > 1)
	{
	  push(@library_names, $lib_name)
	}

        $count++;
      }  
      $first_line = 0;
    } else {

      # insert the tag values, if they are already in the database they should not insert
      # since they are flaged as unique in the database
      my $insert_tag_query = "insert into sage_tags (tagID, sequence) 
				VALUES (" . 
				"'" . $sage_array[0] . "', " .
				"'" . $sage_array[1] . "')";
      my $sth = $dbh->prepare($insert_tag_query);
      $sth->execute;
  
      # Start at -2 so we can skip the first two values (tagID, sequence)
      my $count = -2;
   
      # iterate through the remaining libraries
      foreach my $lib (@sage_array)
      {
        if($count >= 0)
        {
	  my $insert_result_query = "insert into sage_results (tagID, library, result)
	  			VALUES ('" . $sage_array[0] . "', " .
  	  			"'" . $library_names[$count] . "', " .
				"'" . $lib . "')";
	  $sth = $dbh->prepare($insert_result_query);
          $sth->execute;
	}
        $count++;			
      }
    }
			
  }
        print "End: Parsing Sage Input File\n";

}

if($map_sage_to_db)
{

        print "Start: Mapping Sage Tags to Genome\n";

# Iterate through all old sage tags

	# Delete old tagmap
	$dbh->prepare('delete from tagmap')->execute();
 
	my $get_tag_query = "select tagid, sequence from sage_tags";
	my $sth = $dbh->prepare($get_tag_query);
	$sth->execute;
 
	while(my $this_row = $sth->fetchrow_hashref)
	{
	        if($debug)
	        {
	                print "Checking:\t" . $this_row->{tagid};
	        }
	 
	        # now iterate through all sequences and check if this sequence is in one of the contigs
	 
	        my $sequences = Bio::SeqIO->new('-file'         => "<$fasta_bases_file",
	                                        '-format'       => "fasta");
	 
	        # Create a sequence object for the query string
	        my $sequence_string_db = $this_row->{sequence};
	 
	        $sequence_string_db =~ s/[^ATGC]//ig;
	                                                                                                                                                                                 
	        my $this_tag_sequence = Bio::Seq->new ( -display_id     => $this_row->{tagid},
	                                                -seq            => $sequence_string_db);
	        my $sequence_to_check = lc($this_tag_sequence->seq);
	        my $sequence_to_check_rc = lc($this_tag_sequence->revcom()->seq);
	        my @hit_array ;
	        while (my $seq = $sequences->next_seq)
	        {
	                my $last_index = 0;
	                my $this_index = -1;
	                my $num_hits = 0;
	                                                                                                                                                                                 
	                if($debug)
	                {
	                        print "\t" . $seq->display_id . "\n";
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
	 
	                                print   "TAG \t"
	                                . $this_row->{tagid}
	                                . " matches \t" . $seq->display_id
	                                . " location:\t" . $this_index . "\t"
	                                .  (length($this_tag_sequence->seq) + $this_index) . "\t" . $dir . "\n";
	                                push @hit_array, [$this_row->{tagid}, $seq->display_id, $this_index+1, (length($this_tag_sequence->seq) + $this_index), $dir];
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
	                # Tag not found, so give a message
 
	                if($debug)
	                {
	                  print "TAG\t" . $this_row->{tagid} . "\t Not found!\n";
	                }
	                my $insert_query = "insert into tagmap
	                        (tagID, contig, start, stop, direction
	                         ) VALUES (" .
	                        "'" . $this_row->{tagid} . "'," .
	                        "NULL ," .
	                        "NULL ," .
	                        "NULL ," .
	                        "NULL " .
	                         ")";
                        my $insert_ref = $dbh->prepare($insert_query);
                        $insert_ref->execute;
                        #print $insert_query . "\n";
 
 
	        } else
	        {
	               if( (scalar @hit_array) && ($debug) )
	               {
	                        print "TAG\t" . $this_row->{tagid} . "\t Multiple HITS!\n";
	               }
	 
	                # We must now map it to each location
	 
	               foreach my $this_hit(@hit_array)
	               {
	                       my $insert_query = "insert into tagmap
	                       (tagID, contig, start, stop, direction
	                        ) VALUES (" .
	                       "'" . $this_row->{tagid} . "'," .
	                       "'" . $this_hit->[1] . "'," .
	                       "'" . $this_hit->[2] . "'," .
	                       "'" . $this_hit->[3] . "'," .
	                       "'" . $this_hit->[4] . "')";
	                       my $insert_ref = $dbh->prepare($insert_query);
	                       $insert_ref->execute;
	                       #print $insert_query . "\n";
	 
	                } # end foreach tag hit
	       	} # end else hit array == 0
	
	 
	} # end checking all db rows
        print "End: Mapping Sage Tags to Genome\n";

}



if($run_orf_tests)
{

	# This will run testcode, and genescan on each orf in the database and insert those values inside the database

	# Get all orfs in the database and their sequence
	my $query = "
        	select orfid,
		sequence,
	        contig,
	        start,
	        stop,
	        direction
	        FROM
	        orfs where sequence != '' AND sequence is not null AND stop-start >=200";
	my $sth = $dbh->prepare($query);
	$sth->execute;
	
	while(my $this_row = $sth->fetchrow_hashref)
	{
		# Create a temp file to store the sequence in
		my $tempoutfile = "temp/outfile.tmp";

		open(TMPFILE, ">temp/tmp.seq") or die ("Can not open temp file for testcode prediction");
		print TMPFILE ">temp\n";
		print TMPFILE $this_row->{sequence};
		
#		system('/usr/local/bin/tcode -sequence temp/tmp.seq -outfile temp/outfile.tmp -window 200');

		# now parse this file
#		open (RESULT, "$tempoutfile");
#		my $get_line = 0;
#		my $total_frames = 0;
#		my $total_hits = 0;
#		my $running_frame_score = 0;
#		my $test_result = '';
#		while(<RESULT>)
#		{
#			my $this_line = $_;
#			my $is_coding = '';
#
#			if($get_line)
#			{
#				my $this_score = 0;
#				(undef, undef, $this_score, $is_coding) = split(" ", $this_line);
#				if($is_coding ne '')
#				{
#					$total_frames++;
#					$running_frame_score = $running_frame_score + $this_score;
#				}
#				if($is_coding eq 'Coding')
#				{
#					$total_hits++;
#				}
#			}
#			if(substr($this_line, 0, 7) eq "  Start")
#			{
#				$get_line = 1;
#			}
#		}		
#		if($running_frame_score/$total_frames >= 1)
#		{
#			$test_result = 'P';
#			#print $this_row->{orfid} . "\t" . $test_result . "\n";
#		} else
#		{
#			$test_result = 'F';
 #                       #print $this_row->{orfid} . "\t" . $test_result . "\n";
#		}
#		
#		close(RESULT);	


		# Now check Testcode using the testcode_unix package
		
		system('echo "temp/tmp.seq 100 1" | ../testcode_unix');

		my $testcode_result = 0;
		my $testcode_values = 0;
		my $testcode_average = 0;
		my $test_result = 'F';
		open(RESULTTCODE, "temp/tmp_tcode_frd.dat");
		while(<RESULTTCODE>)
		{
			$testcode_result = $testcode_result + $_;
			$testcode_values++;
		}

		$testcode_average = $testcode_result / $testcode_values;
		if($testcode_average >= 9.7)
		{
			$test_result = 'P';
		}else
		{
			$test_result = 'F';
		}

		# Now check the GeneScan output from the testcode_unix package

                my $genescan_result = 0;
                my $genescan_values = 0;
                my $genescan_average = 0;
		my $gene_result = 'F';

		open(RESULTGSCAN, "temp/tmp_gscan.dat");
                while(<RESULTGSCAN>)
                {
                        $genescan_result = $genescan_result + $_;
                        $genescan_values++;
                }

                $genescan_average = $genescan_result / $genescan_values;
                if($genescan_average >= 4.0)
                {
                        $gene_result = 'P';
                }else
                {
                        $gene_result = 'F';
                }




		my $cusage_program = "/usr/local/bin/chips -seqall temp/tmp.seq -outfile temp/outfile.tmp";
                system($cusage_program);

                # now parse this file
                open (RESULT, "$tempoutfile");
		my $codon_usage_score = 0;

                while(<RESULT>)
                {
                        my $this_line = $_;
                        (undef, undef, $codon_usage_score) = split(" ", $this_line);
                }
                close(RESULT);

		
		# Now check codon usage /adaptation index

		my $cai_program = "/usr/local/bin/cai -seqall temp/tmp.seq -outfile temp/outfile.tmp -cfile tables/" . $organism . ".cod -sbegin1 1";
		system($cai_program);
	
                # now parse this file
                open (RESULT, "$tempoutfile");
                my $cai_result_1 = 0;
                while(<RESULT>)
                {
                        my $this_line = $_;
			(undef, undef, undef, $cai_result_1) = split(" ", $this_line);
		}
		close(RESULT);

		# Second Frame
                my $cai_program = "/usr/local/bin/cai -seqall temp/tmp.seq -outfile temp/outfile.tmp -cfile tables/" . $organism . ".cod -sbegin1 2";
                system($cai_program);

                # now parse this file
                open (RESULT, "$tempoutfile");
                my $cai_result_2 = 0;
                while(<RESULT>)
                {
                        my $this_line = $_;
                        (undef, undef, undef, $cai_result_2) = split(" ", $this_line);
                }
                close(RESULT);

		# Third Frame
                my $cai_program = "/usr/local/bin/cai -seqall temp/tmp.seq -outfile temp/outfile.tmp -cfile tables/" . $organism . ".cod -sbegin1 3";
                system($cai_program);

                # now parse this file
                open (RESULT, "$tempoutfile");
                my $cai_result_3 = 0;
                while(<RESULT>)
                {
                        my $this_line = $_;
                        (undef, undef, undef, $cai_result_3) = split(" ", $this_line);
                }
                close(RESULT);



		my $cai_test = 'F';

		if( ($cai_result_1 > $cai_result_2) && 
		    ($cai_result_1 > $cai_result_3) )
		{
			$cai_test = 'P';
		}else
		{
			$cai_test = 'F';
		}
		#print $cai_result_1 . "\n";
                #print $cai_result_2 . "\n";
                #print $cai_result_3 . "\n" . "\n\n";
		


                my $udquery = 	"UPDATE orfs set TestCode = '" 	. $test_result . "',
				 TestScore = '" 		. $testcode_average . "',
                                 GeneScan = '" 			. $gene_result . "',
                                 GeneScanScore = '" 		. $genescan_average . "',
                                 CodonUsage = '"             . $codon_usage_score . "',
                                 CodonPreferenceScore = '"             . $cai_result_1 . "',
                                 CodonPreference = '" 		. $cai_test . "' 
				WHERE orfid = " . $this_row->{orfid};

                my $udh = $dbh->prepare($udquery);
                $udh->execute;


	}


}

if($map_sage_to_orf)
{
	# first find all tags that have at least a 10 value in at least 1 library
	my $sage_query = "select distinct tagid from sage_results where result >= 10";
	my $sth = $dbh->prepare($sage_query);
	$sth->execute;
	my $count = 0;

	while(my $this_row = $sth->fetchrow_hashref)
	{
		# Now find out if this tag maps to only 1 part of the genome
		my $check_map_query = "select count(*) as number_map from tagmap where tagid = '" . $this_row->{tagid} . "'";
		my $ckh = $dbh->prepare($check_map_query);
		$ckh->execute;
		my $number_maps = $ckh->fetchrow_hashref->{number_map};
		if($number_maps == 1)
		{
			# Get all the information about this tag

			my $tag_query = "select tagid, contig, start, stop, direction from tagmap where tagid = '" . $this_row->{tagid} . "'";
			my $tagh = $dbh->prepare($tag_query);
			$tagh->execute;
			
			my $tag = $tagh->fetchrow_hashref;
			# We map to only one part of the genome, so now try and find a good orf close by
			#  A valid orf is a orf that ends close to our start, is going in the same direction as us
			#  And is by a nla3 site towards it's end. If we only map at one place and there is only one orf
			#  Then we can just assign outselves to that orf
			my $get_orf_query = "select orfid, start, stop, contig, TestCode, CodonScore, GeneScan from orfs 
						WHERE
						    delete_fg = 'N'
						AND contig = '" . $tag->{contig} . "'
						AND start-15 < '" . $tag->{start} . "'
						AND stop+15 > '" . $tag->{stop} . "'" .
						"AND direction = '" . $tag->{direction} . "'";

			                        my $get_orf_query = "select orfid, start, stop, contig, TestCode, CodonScore, GeneScan from orfs
                                                WHERE
                                                    delete_fg = 'N'
                                                AND contig = '" . $tag->{contig} . "'
                                                AND ( (direction = '+' AND '+' = '" . $tag->{direction} . "' 
						AND start < '" . $tag->{start} . "'
                                                AND stop+15 > '" . $tag->{stop} . "') OR (" .
						" direction = '-' AND '-' = '" . $tag->{direction} . "'
						AND start-15 < '" . $tag->{start} . "' AND stop > '" . $tag->{stop} . "') )";

			my $findh = $dbh->prepare($get_orf_query);
			$findh->execute;
			if($findh->rows == 1)
			{
				my $found_orf = $findh->fetchrow_hashref;
				print "Sage Tag\t" . $this_row->{tagid} . "\t Maps ORF\t" . $found_orf->{orfid} . "\n";
				$count++;
			}
		}
	}
	print "\n\n Total $count\n";
}

if($map_sage_to_orf_secondary)
{

	# Delete the sage_temp table
	my $query = 'delete from sage_temp';
	$dbh->prepare($query)->execute();

	# Grab only those that map to the genome at leat once
	my $query = 'select distinct tagid from tagmap where contig is not null order by tagid';
	my $sth = $dbh->prepare($query);
	$sth->execute;
	
	
	while(my $tag_outer = $sth->fetchrow_hashref)
	{
		# Now check if this tag is unique or if it maps to multiple places
		my $unique_genome = 0;
		my $unique_transcript = 0;
		my $primary = 0;

		my $num_query = "select count(*) as num_map from tagmap where tagid = '" . $tag_outer->{tagid} . "'";
		my $numh = $dbh->prepare($num_query);
		$numh->execute();
		my $num_map = $numh->fetchrow_hashref->{num_map};
	
		
		# If it maps to only one location, then we set it as unique in transcript and unique in genome
		if($num_map == 1)
		{
			$unique_genome = 1;
			$unique_transcript = 1;
		}
		
		# iterate through all of the tags that this maps to
		my $tag_query = "select distinct tagid, contig, start, stop, direction, assignment, id from tagmap where contig is not null AND tagid = '" . $tag_outer->{tagid} . "'";
		my $taginh = $dbh->prepare($tag_query);
		$taginh->execute();

		while(my $tag = $taginh->fetchrow_hashref)
		{
			# Now get all orf transripts that are overlapping this sage tag (either direction)
	
			my $orf_query = "select distinct orfid, start, stop, contig, direction from orfs
                                                WHERE
                                                    delete_fg = 'N'
                                                AND contig = '" . $tag->{contig} . "'
                                                AND start-$transcript_tail <= '" . $tag->{start} . "'
                                                AND stop+$transcript_tail >= '" . $tag->{stop} . "'";
			my $orfh = $dbh->prepare($orf_query);
			$orfh->execute;
			while(my $orf_hit = $orfh->fetchrow_hashref)
			{
				# Now check if it is at the primary NlaIII site
				my $sequences = Bio::SeqIO->new('-file'         => "<$fasta_bases_file",
	                                        '-format'       => "fasta");
				my $transcript = '';
				my $prim_nlasite = 0;
				my $prim_rev_nlasite = 0;
				
				# Create an array of arrays for each orf hit
				#  orfid,start,stop,TranscriptSequence, PrimaryNlaIII site, PrimaryReverseNlaIII site
				my @orf_array;
				my $found_contig = 0;
				# Iterate through the sequences file to find the contig we are interested in
				while ( (my $seq = $sequences->next_seq) && (!$found_contig) )
				{
					if($seq->display_id eq $tag->{contig})
					{
					        $found_contig = 1;
						
						# Grab start -> stop+$transcript_tail
						if($orf_hit->{direction} eq '+')
						{
						        my $stop = $orf_hit->{stop};
							my $start = $orf_hit->{start};
							my $offset = 0;
						        if($orf_hit->{stop}+$transcript_tail > $seq->length())
							{
							  $stop = $seq->length();
							} else
							{
							  $stop = $orf_hit->{stop} + $transcript_tail;
							}

							if($tag->{direction} ne $orf_hit->{direction})
							{
							    $offset = 3;
							}
							$transcript = $seq->trunc($start,$stop);
							$prim_nlasite = rindex($transcript->seq(), 'CATG');
							if($prim_nlasite == -1)
							{
								print "ORF\t" . $orf_hit->{orfid} . "\t No Primary Sense NlaIII site\n";
							} else
							{
								$prim_nlasite = $prim_nlasite + $orf_hit->{start} + $offset;
							}

							# Now find the primary reverse nla3 site, but first get the hypot reverse transcript 
							# which is the start of the orf -$transcript_tail -> stop
							if( ($orf_hit->{start} - $transcript_tail) < 1)
							{
							    $start = 1;
							}else
							{
							    $start = $orf_hit->{start} - $transcript_tail;
							}
							if($orf_hit->{stop} > $seq->length)
							{
							    $stop = $seq->length();
							} else
							{
								$stop = $orf_hit->{stop};
							}

							$transcript = $seq->trunc($start, $stop);
                                                        $prim_rev_nlasite = index($transcript->seq(), 'CATG');
							if($prim_rev_nlasite == -1)
                                                        {
                                                                print "ORF\t" . $orf_hit->{orfid} . "\t No Primary AntiSense NlaIII site\n";
                                                        } else
                                                        {
                                                                $prim_rev_nlasite = $prim_rev_nlasite + $start + $offset;
                                                        }


							
						} else # Grab start-$transcript_tail -> stop
						{
					
						        my $start = 0;
							my $stop = 0;
							my $offset = 0;
						        if($orf_hit->{start}- $transcript_tail <= 1)
							{
							  $start = 1;
							} else
							{
							  $start = $orf_hit->{start} - $transcript_tail;
							}
							if($orf_hit->{stop} > $seq->length())
							{
								$stop = $seq->length();
							} else
							{
								$stop = $orf_hit->{stop};
							}
							$transcript = $seq->trunc($start, $stop);
							$prim_nlasite = index($transcript->seq(), 'CATG') ;

							if($tag->{direction} eq $orf_hit->{direction})
							{
							    $offset = 3;
							}
                                                        if($prim_nlasite == -1)
                                                        {
                                                                print "ORF\t" . $orf_hit->{orfid} . "\t No Primary Sense NlaIII site\n";
                                                        } else
                                                        {
                                                                $prim_nlasite = $prim_nlasite + $start + $offset;
                                                        }
							$start = $orf_hit->{start};
							if( ($orf_hit->{stop} + $transcript_tail) > $seq->length() )
							{
							    $stop = $seq->length();
							}else
							{
							    $stop = $orf_hit->{stop} + $transcript_tail;
							}
							$transcript = $seq->trunc($start, $stop);

                                                        $prim_rev_nlasite = rindex($transcript->seq(), 'CATG');
                                                        if($prim_rev_nlasite == -1)
                                                        {
                                                                print "ORF\t" . $orf_hit->{orfid} . "\t No Primary Antisense NlaIII siter\n";
                                                        } else
                                                        {
                                                                $prim_rev_nlasite = $prim_rev_nlasite + $start + $offset;
                                                        }

						} # End else
						
						my $tagstart = 0;
						if($tag->{direction} eq "-")
						{
						  $tagstart = $tag->{stop};
						} else
						{
						  $tagstart = $tag->{start};
						}
						my $code = '';
						if( ($orf_hit->{direction} eq $tag->{direction}) && ($tagstart == $prim_nlasite) )
						{
						  $code = 'Primary Sense Tag';
						} elsif( ($orf_hit->{direction} ne $tag->{direction}) && ($tagstart == $prim_rev_nlasite) )
                                                {				 
						    $code = 'Primary Antisense Tag';

                                                } elsif( $orf_hit->{direction} eq $tag->{direction} )
						{
						    $code = 'Alternate Sense Tag';
						} else
						{
						    $code = 'Alternate Antisense Tag';
						}
						

						print "TAG\t" . $tag->{tagid} . "\t" . $tagstart . "\t" . $tag->{direction} . "\tORF\t" . $orf_hit->{orfid} . "\t" . $orf_hit->{direction} . "\t" . $prim_nlasite . "\t" . $prim_rev_nlasite . "\t" . $code . "\t" ;
						my $tag_ins_query = "insert into sage_temp (tagid, start, direction, orfid, orf_direction, tagmapid, tagtype) VALUES (".
						"'" .  $tag->{tagid} . "'," .
						"'" .  $tagstart . "'," .
						"'" .  $tag->{direction} . "'," .
						"'" .  $orf_hit->{orfid} . "'," .
						"'" .  $orf_hit->{direction} . "'," .
                                                "'" .  $tag->{id} . "'," .
						"'" .  $code . "')" ;
						my $taginsh = $dbh->prepare($tag_ins_query);
						$taginsh->execute();
						
						if($unique_genome)
						{
						    print "UG\t";
						} else
						{
						    print "NG\t";
						}

						if($unique_transcript)
						{
						    print "UT\t";
						}else
						{
						    print "NT\t";
						}
						print "\n";
					
					} # end if this is the correct contig
				} # end while I am iterating through all sequences
			} # end while I am done with all orfs that overlap this sagetag

		} # end if this tag maps to only one part of the genome

	} # end iterate through all tags
}


if($map_sage_to_orf_tert)
{

	# Delete the orftosage table;

	$dbh->prepare('delete from orftosage')->execute();

	# Now iterate through all tags again
	# Check how many orf hits a tag has, 
		#if it has only one, mark it, otherwise we need to make a decision
		# if it has more then one, and only one of them is a primary NlaIII site, mark that orf as the correct orf
		# If none lie on a primary NlaIII site and only one is on an alternate forward site, mark that --NOT DONE
		# If all are on alternate reverse sites, and only one is on a primary alternate, mark that     -- NOT DONE

	my $query = "select distinct tagid from tagmap";
	my $tagh = $dbh->prepare($query);
	$tagh->execute();
	my $continue = 1;
	my $ins_query = "insert into orftosage (orfid, tagid, tagtype, tagmapid, unique_genome_fg, unique_trans_fg) VALUES (?, ?, ?, ?, ?, ?)";
        my $insh = $dbh->prepare($ins_query);
	while(my $tags = $tagh->fetchrow_hashref)
	{
		$continue = 1;
		my $tagquery = "select tagid, start, direction, orfid, orf_direction, tagtype, tagmapid from sage_temp where tagid = '" . $tags->{tagid} . "'";
		my $tagmaph = $dbh->prepare($tagquery);
		$tagmaph->execute();

		# check if this tag is unique to the genome
		my $unique_genome = 0;

		my $ugh = $dbh->prepare("select count(*) as num_tag from tagmap where tagid = '" . $tags->{tagid} . "'");
		$ugh->execute;
		my $ugnum = $ugh->fetchrow_hashref;
		if($ugnum->{num_tag} == 1)
		{
			$unique_genome = 1;
		}
		

		# If we only map to one place and one orf, this will only return one record
		if($tagmaph->rows == 0)
		{
			# Do nothing
		} elsif($tagmaph->rows == 1)
		{
			# Now get the record and update 
			my $this_row = $tagmaph->fetchrow_hashref;
			$insh->execute($this_row->{orfid}, $this_row->{tagid}, $this_row->{tagtype}, $this_row->{tagmapid}, $unique_genome, 1);
			
		} else
		{
			# Check if only one of these is a primary sense tag
			my $checkh = $dbh->prepare("select orfid, tagmapid from sage_temp where tagid = '" . $tags->{tagid} . "' AND tagtype = 'Primary Sense Tag'");
			$checkh->execute();
			if($checkh->rows == 0)
			{
				# This is not a primary sense tag for any orf
				$continue = 1;
			} elsif($checkh->rows > 1)
			{
				# We can not make any observations on this due to it maping to multiple primaries
				$continue = 1;
			} elsif($checkh->rows == 1)
			{
				my $this_row = $checkh->fetchrow_hashref;
				$insh->execute($this_row->{orfid}, $tags->{tagid}, 'Primary Sense Tag', $this_row->{tagmapid}, 0, 0);
				$continue = 0;
			}

		} # END ELSE THIS IS IN ZERO OR ONE TIMES ONLY

	} # END WHILE ITERATING THROUGH TAGS
 

}



if($find_repeat)
{

	# Run the repeat finder program
	system($repeatFinder_bin . ' -in ' . $fasta_bases_file . ' -out ' . 'temp/' . $organism . '.repeat ' . $organism);
	
	# Parse the file and place into database
	open(REPEAT, 'temp/' . $organism . '.repeat') or die ("Can not open $organism.repeat");;
	open(REPEATSUPEROUT, ">", $gff_supercontig_repeat_file );
        open(REPEATOUT, ">", $gff_repeat_file );

	while(<REPEAT>)
	{
		my $line = $_;
		if(substr($line, 0, 1) eq '#')
		{
			# Do nothing
		} else
		{
			chomp($line);
			my ($contig, $rep_name, $start, $stop, $class) = split("\t", $line);
		        my $get_contig_pos_query = "select contig_number, super_id,
		                        contig_start_super_base, contig_length, modified_contig_start_base
		                        from links where concat('contig_', contig_number) = '" . $contig . "'";
		        my $posh = $dbh->prepare($get_contig_pos_query);
		        $posh->execute;
		        my $super_result = $posh->fetchrow_hashref;

		        my $new_start = 0;
		        my $new_stop = 0;
		        if($minimum_gap_fg)
		        {
		                $new_start = $start + $super_result->{modified_contig_start_base};
		                $new_stop  = $stop +  $super_result->{modified_contig_start_base};
		        } else
		        {
		                $new_start = $start + $super_result->{contig_start_super_base};
		                $new_stop  = $stop +  $super_result->{contig_start_super_base};
		        }

		        print REPEATSUPEROUT       "supercontig_" . $super_result->{super_id}      . "\t" .
		                        'ClosureRepeatFinder-3.7'     . "\t" .
		                        'repeat_region'                   . "\t" .
		                        $new_start              . "\t" .
		                        $new_stop               . "\t" .
		                        '.'                     . "\t" .
		                        '.'  			. "\t" .
		                        '.'                     . "\t" .
		                        "Sequence $rep_name" 	. "\n";
		        print REPEATOUT $contig     		. "\t" .
                                        'ClosureRepeatFinder-3.7'     . "\t" .
                                        'repeat_region'         . "\t" .
		                        $start       		. "\t" .
		                        $stop       		. "\t" .
		                        '.'                     . "\t" .
		                        '.'  			. "\t" .
		                        '.'                     . "\t" .
		                        "Sequence $rep_name"	. "\n";

		        $posh->finish;
		}
	}
	
}


if($create_blast_db)
{

	# We need to create a fasta file first, then call the formatdb command on it

	# Contigs first, this is the same as our assembly.bases
	system('cp ' . $fasta_bases_file . ' ' . $blast_db_dir . '/' . $organism );
        system($blast_bin_dir . 'formatdb -t ' . $organism . ' -i ' . $blast_db_dir . '/' . $organism . ' -p F -o T');

	

	# ORFs as nt sequences
                                                                                                                                                                                  
        open(ORFFILE, '>', $blast_db_dir . '/' . $organism . '_orfs_nt');
                                                                                                                                                                                  
        my $query = "select orfid, sequence from orfs where delete_fg = 'N'";
                                                                                                                                                                                  
        my $sth = $dbh->prepare($query);
        $sth->execute;
                                                                                                                                                                                  
        while(my $orf = $sth->fetchrow_hashref)
        {
                print ORFFILE ">" . $orf->{orfid} . "\n";
                # load this sequence into a sequence object
                my $this_sequence_obj = Bio::Seq->new( -id  => 'my_sequence',
                                -seq =>$orf->{sequence});
                                                                                                                                                                                  
                                                                                                                                                                                  
                print ORFFILE $this_sequence_obj->seq();
                print ORFFILE "\n";
        }
        close ORFFILE;

	system($blast_bin_dir . 'formatdb -t ' . $organism . '_orfs_nt -i ' . $blast_db_dir . '/' . $organism . '_orfs_nt -p F -o T');



	# Translated ORFs

	open(ORFFILE, '>', $blast_db_dir . '/' . $organism . '_orfs_aa');

	my $query = "select orfid, sequence from orfs where delete_fg = 'N'";
                                                                                                                                                                                  
	my $sth = $dbh->prepare($query);
	$sth->execute;

	while(my $orf = $sth->fetchrow_hashref)
	{
	        print ORFFILE ">" . $orf->{orfid} . "\n";
	        # load this sequence into a sequence object
	        my $this_sequence_obj = Bio::Seq->new( -id  => 'my_sequence',
	                        -seq =>$orf->{sequence});
	                                                                                                                                                                                  
	        print ORFFILE $this_sequence_obj->translate()->seq();
	        print ORFFILE "\n";
	}
	close ORFFILE;

        system($blast_bin_dir . 'formatdb -t ' . $organism . '_orfs_aa -i ' . $blast_db_dir . '/' . $organism . '_orfs_aa -p T -o T');


	# Unused Reads

        open(READFILE, '>', $blast_db_dir . '/' . $organism . '_unused_reads_nt');
                                                                                                                                                                                  
        my $query = "select reads.read_name as readname, reads_assembly.read_name, reads_bases.bases 
			FROM 
			reads_bases, 
			reads left join reads_assembly ON reads.read_name=reads_assembly.read_name 
			WHERE reads.read_name = reads_bases.read_name AND reads_assembly.read_name is null";
                                                                                                                                                                                  
        my $sth = $dbh->prepare($query);
        $sth->execute;
                                                                                                                                                                                  
        while(my $read = $sth->fetchrow_hashref)
        {
                print READFILE ">" . $read->{readname} . "\n";
                # load this sequence into a sequence object
                my $this_sequence_obj = Bio::Seq->new( -id  => 'my_sequence',
                                -seq =>$read->{bases});

                print READFILE $this_sequence_obj->seq();
                print READFILE "\n";
        }
        close READFILE;
                                                                                                                                                                                  
        system($blast_bin_dir . 'formatdb -t ' . $organism . '_unused_reads_nt -i ' . $blast_db_dir . '/' . $organism . '_unused_reads_nt -p F -o T');
                                                                                                                                                                                  


}
#######START OUTPUT GFF GILE#################
if($output_gff)
{
	print "Started creating $gff_output_file file\n";
	open(GFF,'>', "$gff_output_file") or die ("Can not open $gff_output_file file");

	# First Annotate each contig
	my $gff_query = "select distinct contig_number, contig_length from reads_assembly";
	my $sth = $dbh->prepare($gff_query);
        $sth->execute;
        while(my $this_row = $sth->fetchrow_hashref)
        {
          print GFF "contig_" . $this_row->{contig_number} . "\t" .
	            "ARACHNE" . "\t" .
                    "contig" . "\t" .
                    "1" . "\t" .
                    $this_row->{contig_length} . "\t" .
		    "." . "\t" .
                    "." . "\t" .
                    "." . "\t" .
                    "Sequence \"contig_" . $this_row->{contig_number} . "\"" . "\n"; 
	  print ".";
        }
	$sth->finish;
	
        # Now Annotate each supercontig
        my $super_query = "select distinct super_id, contigs_in_super, 
		      contig_number, ordinal_number, contig_length from links";
        $sth = $dbh->prepare($super_query);
        $sth->execute;
        while(my $this_row = $sth->fetchrow_hashref)
        {
          print GFF "contig_" . $this_row->{contig_number} . "\t" .
                    "ARACHNE" . "\t" .
                    "supercontig" . "\t" .
                    "1" . "\t" .
                    $this_row->{contig_length} . "\t" .
                    "+" . "\t" .
                    "." . "\t" .
                    "." . "\t" .
                    "Sequence \"supercontig_" . $this_row->{super_id} 
		        ." (" . $this_row->{ordinal_number} . "/" . $this_row->{contigs_in_super} . ")\"" . "\n";
          print ".";
        }
        $sth->finish;

        # Now Annotate each supercontig with link
        my $super_query = "select distinct super_id, contigs_in_super,
                      contig_number, ordinal_number, contig_length from links";
        $sth = $dbh->prepare($super_query);
        $sth->execute;
        while(my $this_row = $sth->fetchrow_hashref)
        {
          print GFF "contig_" . $this_row->{contig_number} . "\t" .
                    "ARACHNE" . "\t" .
                    "supercontig" . "\t" .
                    "1" . "\t" .
                    $this_row->{contig_length} . "\t" .
                    "+" . "\t" .
                    "." . "\t" .
                    "." . "\t" .
                    "Sequence \"supercontig_" . $this_row->{super_id}
                        . '"' . "\n";
          print ".";
        }
        $sth->finish;


	#Now Annotate each read
        my $gff_query = "SELECT distinct
			reads.read_name,
			r_assem.contig_number,
			reads.center_name,
			r_assem.trim_read_in_contig_start,
			r_assem.trim_read_in_contig_stop,
			r_assem.orientation,
			reads.plate_id,
			reads.well_id,
			reads.template_id,
			reads.library_id,
			r_assem.read_pair_name,
			r_assem.read_pair_contig_number,
			r_assem.observed_insert_size,
			r_assem.given_insert_size,
			r_assem.given_insert_std_dev,
			r_assem.observed_inserted_deviation,
			links.super_id as super_contig_number,
			links.contig_start_super_base
		      FROM reads,
			reads_assembly r_assem,
			links
		      WHERE reads.read_name = r_assem.read_name
			AND links.contig_number = r_assem.contig_number";
        $sth = $dbh->prepare($gff_query);
        $sth->execute;
        while(my $this_row = $sth->fetchrow_hashref)
        {
          print GFF "contig_" . $this_row->{contig_number} . "\t" .
                    "read" . "\t";

       # Check if this read is special :
       # 1- it's partner read is unplaced and not missing(it has a partner)
       # 2- it's partner read is on another contig
       # 3- it's partner read is on another supercontig
       # 4- it's partner is unplaced and missing
       # 5- it's partner is in the next contig and the gap is positive
       # 6- it's partner is in the next contig and the gap is negative
       # 7- it's partner is not in the next ordinal contig


          if($this_row->{read_pair_contig_number} == $this_row->{contig_number}) #they are in the same contig
	  {
	     print GFF "read" . "\t";

	  }
	  elsif($this_row->{read_pair_contig_number} == "") # there is no partner read
          {
            my $check_for_partner_query = "select read_name from reads where template_id = '".
                $this_row->{template_id} . "' AND read_name != '" . $this_row->{read_name} . "'"; 
	    my $check_for_partner_result = $dbh->prepare($check_for_partner_query);
            $check_for_partner_result->execute;
            my $num_rows_ret = $check_for_partner_result->rows;
            if($num_rows_ret == 0)
	    {
              print GFF "missing-partner" . "\t";
	    } else
   	    {
	      print GFF "unplaced-partner" . "\t";
            }
          } else #they are in different contigs, now check if they are in same supercontig
	  {
            my $check_partner_super_query = "
                        SELECT distinct
                        reads.read_name,
                        r_assem.contig_number,
                        r_assem.read_pair_name,
                        r_assem.read_pair_contig_number,
                        links.super_id as super_contig_number,
			r_assem.orientation
                      FROM reads,
                        reads_assembly r_assem,
                        links
                      WHERE reads.read_name = r_assem.read_name
                        AND links.contig_number = r_assem.contig_number
                        AND reads.read_name = '" . $this_row->{read_pair_name} . "'" ;
		my $inner_check = $dbh->prepare($check_partner_super_query);
		$inner_check->execute;
	        my $inner_row = $inner_check->fetchrow_hashref;

		if($inner_row->{super_contig_number} == $this_row->{super_contig_number})
		{
		  
		  # Now check if the two contigs overlap
################################################################
		  my $query = '
	        select contig_number,
	        super_id,
	        bases_in_super,
	        contigs_in_super,
	        ordinal_number,
	        contig_length,
	        gap_before_contig,
	        gap_after_contig
	        FROM
	        links  
	        WHERE
		super_id = ' .  $this_row->{super_contig_number} .
	        ' ORDER BY super_id, ordinal_number';
		  my $links_result = $dbh->prepare($query);
		  $links_result->execute;
		  my $contig_one = $this_row->{contig_number};
		  my $contig_two = $this_row->{read_pair_contig_number};
		  my $read_one_dir = $this_row->{orientation};
		  my $read_two_dir = $inner_row->{orientation};
		
		  my $contig_one_start = 0;
		  my $contig_two_start = 0;
		  my $contig_one_end = 0;
		  my $contig_two_end = 0;
	
		  my $last_super_id = '-1';
		  my $super_running_total = 0;
		  while(my $links_array = $links_result->fetchrow_hashref)
		  {
                                                                                              
		    if($links_array->{super_id} != $last_super_id)
		    {
		      $super_running_total = 0;
		    }
		                                                                                              
		    my $start_val = $super_running_total + $links_array->{gap_before_contig};
		    my $end_val = $start_val + $links_array->{contig_length};
		    $super_running_total = $end_val;
		    $last_super_id = $links_array->{super_id};
		
		    if($links_array->{contig_number} == $contig_one)
		    {
		      $contig_one_start = $start_val;
		      $contig_one_end   = $end_val;
		    }elsif($links_array->{contig_number} == $contig_two)
		    {
		      $contig_two_start = $start_val;
		      $contig_two_end   = $end_val;
		    }
	                                                                                              
		  }
		  $links_result->finish;

###############################################################
		  my $gap_type = "";
		  if( ($read_one_dir eq "+") && ($read_two_dir eq "-") )
		  {
		    if($contig_one_end < $contig_two_start)  
		    {
		      $gap_type = "partner-different-contig-positive-gap";
		    } else
		    {
                      $gap_type = "partner-different-contig-negative-gap";
		    }
		  } elsif( ($read_one_dir eq "-") && ($read_two_dir eq "+") )
		  {
                    if($contig_one_start > $contig_two_end)
                    {
                      $gap_type = "partner-different-contig-positive-gap";
                    } else
                    {
                      $gap_type = "partner-different-contig-negative-gap";
                    }
		  } elsif( $read_one_dir eq $read_two_dir)
		  {
		    # There is something wierd about this read, pairs in same direction
		    $gap_type = "partner-different-contig-same-direction";
		  } else
		  {
		    $gap_type = "partner-different-contig-exception";
		  }

		  print GFF $gap_type . "\t";
		} else
          	{
           	   print GFF "partner-different-supercontig" . "\t";
          	} 
          }
          print GFF $this_row->{trim_read_in_contig_start} . "\t" .
                    $this_row->{trim_read_in_contig_stop} . "\t" .
                    "." . "\t" .
                    $this_row->{orientation} . "\t" .
                    "." . "\t" .
                    "Sequence \"" . $this_row->{read_name} . "\"" . 
		    ' ; ReadPair "' . $this_row->{read_pair_name} . '" ; ReadPairContig "' .
		    $this_row->{read_pair_contig_number} . '" ; ReadContig "' . $this_row->{contig_number} . 
                    '" ; TemplateID "' . $this_row->{template_id} .
		    '" ; GivenInsertSize "' . $this_row->{given_insert_size} . 
		    '" ; GivenInsertStdDev "' . $this_row->{given_insert_std_dev} . 
		    '" ; ObservedInsertSize "' . $this_row->{observed_insert_size} . 
		    '" ; ObservedInsertStdDev "' . $this_row->{observed_inserted_deviation} . '"' .
		    "\n";
	  print ".";

        }
        $sth->finish;
	print "Completed creating $gff_output_file file\n";
#######END OUTPUT GFF GILE##############

}
if($create_supercontig_map)
{

  print "Started creating $gff_supercontig_output_file file\n";
  open(SUPER,'>', "$gff_supercontig_output_file") or die ("Can not open $gff_supercontig_output_file file");


  # Place Supercontigs

  my $query = '
        select distinct
        super_id,
        bases_in_super
        FROM
        links';
  my $links_result = $dbh->prepare($query);
  $links_result->execute;
                                                                                
                                                                                
  while(my $links_array = $links_result->fetchrow_hashref)
  {


    print SUPER     "supercontig_" . $links_array->{super_id}  . "\t" .
                    "ARACHNE" . "\t" .
                    "supercontig" . "\t" .
                    "1" . "\t" .
                    $links_array->{bases_in_super} . "\t" .
                    "." . "\t" .
                    "." . "\t" .
                    "." . "\t" .
                    "Sequence \"supercontig_" . $links_array->{super_id} . "\"" . "\n";
  }
  $links_result->finish;




	# Place Contigs

	my $query = '
		select contig_number,
		super_id,
		bases_in_super,
		contigs_in_super,
		ordinal_number,
		contig_length,
		gap_before_contig,
		gap_after_contig,
		contig_start_super_base,
		modified_contig_start_base,
		modified_bases_in_super
		FROM
		links ORDER BY super_id, ordinal_number';
	my $links_result = $dbh->prepare($query);
	$links_result->execute;

	my $last_super_id = '';
	my $super_running_total = 0;
        while(my $links_array = $links_result->fetchrow_hashref)
        {
		my $start_val = 0;
		my $end_val = 0;
                if($minimum_gap_fg)
                {
                        $start_val = $links_array->{modified_contig_start_base};
                        $end_val = $start_val + $links_array->{contig_length};
                                                                                                                                                                                    
                } else
                {
                        $start_val = $links_array->{contig_start_super_base};
                        $end_val = $start_val + $links_array->{contig_length};
                }
                        print SUPER "supercontig_" . $links_array->{super_id} . "\t" .
                            "ARACHNE" . "\t" .
                            "contig" . "\t" .
                            $start_val . "\t" .
                            $end_val . "\t" .
                            "." . "\t" .
                            "." . "\t" .
                            "." . "\t" .
                            "Sequence \"contig_" . $links_array->{contig_number} . '"' . "\n";
                        my $update_base_query = "UPDATE links set contig_start_super_base = '" . $start_val . "'
                                WHERE contig_number = '" . $links_array->{contig_number} ."'";
                        my $update_base_result = $dbh->prepare($update_base_query);
                        $update_base_result->execute;
                                                                                                                                                                                    
	}

	$links_result->finish;


}

if($create_supercontig_read_map)
{

	print("Starting to create $gff_supercontig_read_output_file\n");
	open(SUPERGFF, '>', "$gff_supercontig_read_output_file") or die ("Can not open $gff_supercontig_read_output_file");

	my $query = '
	        select distinct
	        super_id,
	        bases_in_super,
		modified_bases_in_super
	        FROM
	        links';
	my $links_result = $dbh->prepare($query);
	$links_result->execute;
                                                                                                 
                                                                                                 
	while(my $links_array = $links_result->fetchrow_hashref)
	{
		my $stop = 0;

		if($minimum_gap_fg)
		{
			$stop = $links_array->{modified_bases_in_super}
		} else
		{
			$stop = $links_array->{bases_in_super}
		}
    print SUPERGFF     "supercontig_" . $links_array->{super_id}  . "\t" .
                    "ARACHNE" . "\t" .
                    "supercontig" . "\t" .
                    "1" . "\t" .
                    $stop . "\t" .
                    "." . "\t" .
                    "." . "\t" .
                    "." . "\t" .
                    "Sequence \"supercontig_" . $links_array->{super_id} . "\"" . "\n";
                                                                                                 
	}

	$links_result->finish;


  # Place Contigs
                                                                                                 
        # Place Contigs
                                                                                                                                                                                    
        my $query = '
                select contig_number,
                super_id,
                bases_in_super,
                contigs_in_super,
                ordinal_number,
                contig_length,
                gap_before_contig,
                gap_after_contig,
                contig_start_super_base,
                modified_contig_start_base,
                modified_bases_in_super
                FROM
                links ORDER BY super_id, ordinal_number';
        my $links_result = $dbh->prepare($query);
        $links_result->execute;
                                                                                                                                                                                    
        my $last_super_id = '';
        my $super_running_total = 0;
        while(my $links_array = $links_result->fetchrow_hashref)
        {
		my $start_val = 0;
		my $end_val = 0;
                if($minimum_gap_fg)
                {
                        $start_val = $links_array->{modified_contig_start_base};
                        $end_val = $start_val + $links_array->{contig_length};
                                                                                                                                                                                    
                } else
                {
                        $start_val = $links_array->{contig_start_super_base};
                        $end_val = $start_val + $links_array->{contig_length};
                }
                        print SUPERGFF "supercontig_" . $links_array->{super_id} . "\t" .
                            "ARACHNE" . "\t" .
                            "contig" . "\t" .
                            $start_val . "\t" .
                            $end_val . "\t" .
                            "." . "\t" .
                            "." . "\t" .
                            "." . "\t" .
                            "Sequence \"contig_" . $links_array->{contig_number} . '"' . "\n";
                        my $update_base_query = "UPDATE links set contig_start_super_base = '" . $start_val . "'
                                WHERE contig_number = '" . $links_array->{contig_number} ."'";
                        my $update_base_result = $dbh->prepare($update_base_query);
                        $update_base_result->execute;

        }

        $links_result->finish;

	#Now Annotate each read
	my $gff_query = "SELECT distinct
			reads.read_name,
			r_assem.contig_number,
			reads.center_name,
			r_assem.trim_read_in_contig_start,
			r_assem.trim_read_in_contig_stop,
			r_assem.orientation,
			reads.plate_id,
			reads.well_id,
			reads.template_id,
			reads.library_id,
			r_assem.read_pair_name,
			r_assem.read_pair_contig_number,
                        r_assem.observed_insert_size,
                        r_assem.given_insert_size,
                        r_assem.given_insert_std_dev,
                        r_assem.observed_inserted_deviation,
			links.super_id as super_contig_number,
			links.contig_start_super_base,
			links.modified_contig_start_base
		      FROM reads,
			reads_assembly r_assem,
			links
		      WHERE reads.read_name = r_assem.read_name
			AND links.contig_number = r_assem.contig_number";
	my $sth = $dbh->prepare($gff_query);
	$sth->execute;
	while(my $this_row = $sth->fetchrow_hashref)
	{
		print SUPERGFF "supercontig_" . $this_row->{super_contig_number} . "\t" .
                	"read" . "\t";

	       # Check if this read is special :
	       # 1- it's partner read is unplaced and not missing(it has a partner)
	       # 2- it's partner read is on another contig
	       # 3- it's partner read is on another supercontig
	       # 4- it's partner is unplaced and missing
	       # 5- it's partner is in the next contig and the gap is positive
	       # 6- it's partner is in the next contig and the gap is negative
	       # 7- it's partner is not in the next ordinal contig


		if($this_row->{contig_number} == $this_row->{read_pair_contig_number}) #they are in the same contig
		{
	     		print SUPERGFF "read" . "\t";

	  	} elsif($this_row->{read_pair_contig_number} eq "") # there is no partner read
	        {
            		my $check_for_partner_query = "select read_name from reads where template_id = '".
			$this_row->{template_id} . "' AND read_name != '" . $this_row->{read_name} . "'"; 
			my $check_for_partner_result = $dbh->prepare($check_for_partner_query);
		        $check_for_partner_result->execute;
		        my $num_rows_ret = $check_for_partner_result->rows;
		        if($num_rows_ret == 0)
			{
		              print SUPERGFF "missing-partner" . "\t";
		   	} else
   	    		{
			      print SUPERGFF "unplaced-partner" . "\t";
            		}
		} else #they are in different contigs, now check if they are in same supercontig
	  	{
            		my $check_partner_super_query = "
                        	SELECT distinct
	                        reads.read_name,
	                        r_assem.contig_number,
	                        r_assem.read_pair_name,
	                        r_assem.read_pair_contig_number,
	                        links.super_id as super_contig_number,
				r_assem.orientation
	                      FROM reads,
	                        reads_assembly r_assem,
	                        links
	                      WHERE reads.read_name = r_assem.read_name
	                        AND links.contig_number = r_assem.contig_number
	                        AND reads.read_name = '" . $this_row->{read_pair_name} . "'" ;
			my $inner_check = $dbh->prepare($check_partner_super_query);
			$inner_check->execute;
		        my $inner_row = $inner_check->fetchrow_hashref;

			if($inner_row->{super_contig_number} == $this_row->{super_contig_number})
			{
		  
				# Now check if the two contigs overlap

				my $query = '
			        select contig_number,
			        super_id,
			        bases_in_super,
			        contigs_in_super,
			        ordinal_number,
			        contig_length,
			        gap_before_contig,
			        gap_after_contig
			        FROM
			        links  
			        WHERE
				super_id = ' .  $this_row->{super_contig_number} .
			        ' ORDER BY super_id, ordinal_number';
				  my $links_result = $dbh->prepare($query);

				$links_result->execute;
				my $contig_one = $this_row->{contig_number};
				my $contig_two = $this_row->{read_pair_contig_number};
				my $read_one_dir = $this_row->{orientation};
				my $read_two_dir = $inner_row->{orientation};
		
				my $contig_one_start = 0;
				my $contig_two_start = 0;
				my $contig_one_end = 0;
				my $contig_two_end = 0;
	
				my $last_super_id = '-1';
				my $super_running_total = 0;
				while(my $links_array = $links_result->fetchrow_hashref)
		  		{
                                                                                              
			    		if($links_array->{super_id} != $last_super_id)
			    		{
			      			$super_running_total = 0;
			    		}
		                                                                                              
			    		my $start_val = $super_running_total + $links_array->{gap_before_contig};
			    		my $end_val = $start_val + $links_array->{contig_length};
			    		$super_running_total = $end_val;
			    		$last_super_id = $links_array->{super_id};
		
			    		if($links_array->{contig_number} == $contig_one)
			    		{
			      			$contig_one_start = $start_val;
			      			$contig_one_end   = $end_val;
			    		}elsif($links_array->{contig_number} == $contig_two)
			    		{
			      			$contig_two_start = $start_val;
			      			$contig_two_end   = $end_val;
			    		}
		                                                                                              
			  	}
			  	
				$links_result->finish;

###############################################################
				my $gap_type = "";
			  	if( ($read_one_dir eq "+") && ($read_two_dir eq "-") )
			  	{
			    		if($contig_one_end < $contig_two_start)  
			    		{
			      			$gap_type = "partner-different-contig-positive-gap";
			    		} else
			    		{
	                      			$gap_type = "partner-different-contig-negative-gap";
			    		}
			  	} elsif( ($read_one_dir eq "-") && ($read_two_dir eq "+") )
			  	{
	                    		if($contig_one_start > $contig_two_end)
	                    		{
	                      			$gap_type = "partner-different-contig-positive-gap";
	                    		} else
	                    		{
	                      			$gap_type = "partner-different-contig-negative-gap";
	                    		}
			  	} elsif( $read_one_dir eq $read_two_dir)
			  	{
			    		# There is something wierd about this read, pairs in same direction
			    		$gap_type = "partner-different-contig-same-direction";
			  	} else
			  	{
			    		$gap_type = "partner-different-contig-exception";
			  	}
	
			  	print SUPERGFF $gap_type . "\t";
			} else
	          	{
	           		print SUPERGFF "partner-different-supercontig" . "\t";
	          	} 
		}
	
		my $read_start_val = 0;
		my $read_stop_val = 0;
		if($minimum_gap_fg)
		{
			$read_stop_val = $this_row->{trim_read_in_contig_stop} + $this_row->{modified_contig_start_base};
	            	$read_start_val = $this_row->{modified_contig_start_base} + $this_row->{trim_read_in_contig_start};
		} else
		{
			$read_stop_val = $this_row->{trim_read_in_contig_stop} + $this_row->{contig_start_super_base};
			$read_start_val = $this_row->{contig_start_super_base} + $this_row->{trim_read_in_contig_start};
		}
	
	        print SUPERGFF $read_start_val . "\t" .
	                    $read_stop_val . "\t" .
	                    "." . "\t" .
	                    $this_row->{orientation} . "\t" .
	                    "." . "\t" .
	                    "Sequence \"" . $this_row->{read_name} . "\"" .
                    	    ' ; ReadPair "' . $this_row->{read_pair_name} . '" ; ReadPairContig "' .
	                    $this_row->{read_pair_contig_number} . '" ; ReadContig "' . $this_row->{contig_number} .
                    '" ; TemplateID "' . $this_row->{template_id} .
                    '" ; GivenInsertSize "' . $this_row->{given_insert_size} .
                    '" ; GivenInsertStdDev "' . $this_row->{given_insert_std_dev} .
                    '" ; ObservedInsertSize "' . $this_row->{observed_insert_size} .
                    '" ; ObservedInsertStdDev "' . $this_row->{observed_inserted_deviation} . '"' .
		    "\n";
		print ".";
		#$innner_check->finish;
	}
	
        $sth->finish;
	print "Completed creating output.gff file\n";

	#######END OUTPUT SUPERGFF GILE##############
	                                                                                                 


}
if($create_supercontig_orf_from_db)
{

  my $query = "
        select orfid, 
	contig,
	start,
	stop,
	direction,
	attributes,
	TestCode,
	CodonPreference,
	CodonPreferenceScore,
	CodonUsage,
	TestScore,
	GeneScan,
	GeneScanScore,
	source
        FROM
        orfs where delete_fg = 'N' ";
  my $result = $dbh->prepare($query);
  $result->execute;


  open(ORFDBSCFILE, ">$gff_supercontig_orf_output_file");
  open(ORFDBFILE, ">$gff_orf_output_file");


  while(my $this_row = $result->fetchrow_hashref)
  {
	#my $attributes = $this_row->{attributes};


	#if($attributes eq '')
	#{
		# Do nothing
	#}else 
	#{
	#	$attributes = 'old_' . $attributes;
	#}
	
	# find best blast hit

	my $best_hit = '';
	my $blast_query = "select  hit_name, description, evalue from blast_results where sequence_type_id = 2 AND db = 3 AND algorithm = 3 AND idname = '" . $this_row->{orfid} . "' order by evalue limit 1";
	my $blasth = $dbh->prepare($blast_query);
	$blasth->execute;

	if($blasth->rows == 0)
	{
		$best_hit = 'No significant SwissProt Hit';
	} else
	{
		my $blast_row = $blasth->fetchrow_hashref;
		if($blast_row->{evalue} <= 1e-3)
		{
			$best_hit = $blast_row->{hit_name} . '|' . $blast_row->{description} . '|' . $blast_row->{evalue};
		} else
		{
			$best_hit = 'No significant SwissProt Hit';
		}	
	}

	# Find best NR hit

        my $best_nr_hit = '';
        my $blast_query = "select  hit_name, description, evalue from blast_results where sequence_type_id = 2 AND db = 2 AND algorithm = 3 AND idname = '" . $this_row->{orfid} . "'  AND (description not like '%ATCC 50803%' OR description like '%gb|%') order by evalue limit 1";
        my $blasth = $dbh->prepare($blast_query);
        $blasth->execute;

        if($blasth->rows == 0)
        {
                $best_nr_hit = 'No significant nr Hit';
        } else
        {
	  my $blast_row = $blasth->fetchrow_hashref;
	  if ($blast_row->{evalue} <= 1e-10)
	    {
	      $best_nr_hit = $blast_row->{hit_name} . '|' . $blast_row->{description} . '|' . $blast_row->{evalue};
	    } else
	      {
		$best_nr_hit = 'No significant nr Hit';
	      }
        }

        my $best_pfam_hits = '';
        my $pfam_query = "select  hit_name, description, evalue from blast_results where sequence_type_id = 2 AND db = 4 AND algorithm = 4 AND idname = '" . $this_row->{orfid} . "' AND evalue <= 1e-3 order by evalue";
        my $blasth = $dbh->prepare($pfam_query);
        $blasth->execute;

        if($blasth->rows == 0)
        {
                $best_pfam_hits = ' Pfam_ls "No significant Pfam Hit" ;';
        } else
	{
	  while(my $blast_row = $blasth->fetchrow_hashref)
	  {
	    $best_pfam_hits = $best_pfam_hits .  ' Pfam_ls "' . $blast_row->{hit_name} . '|' . $blast_row->{description} . '|' . $blast_row->{evalue} . '" ;';
          }
        }

	$best_nr_hit =~ s/\"/\\"/ig;
	$best_nr_hit =~ s/\'/\\'/ig;

	my $attributes = 'Sequence "' . $this_row->{orfid} . 
				'" ; TestCode "' 	. $this_row->{TestCode} . 
                                '" ; TestCodeScore "'        . $this_row->{TestScore} .
                                '" ; GeneScan "'        . $this_row->{GeneScan} .
                                '" ; GeneScanScore "'        . $this_row->{GeneScanScore} .
				'" ; CodonPreference "' . $this_row->{CodonPreference} .
                                '" ; CodonPreferenceScore "' . $this_row->{CodonPreferenceScore} .
                                '" ; CodonUsage "' . $this_row->{CodonUsage} .
                                '" ; nr "' . $best_nr_hit .
                                '" ; SwissProt "' . $best_hit . '" ;' .
				     $best_pfam_hits;

	my $get_contig_pos_query = "select contig_number, super_id, 
			contig_start_super_base, contig_length, modified_contig_start_base
			from links where concat('contig_', contig_number) = '" . $this_row->{contig} . "'";
        my $posh = $dbh->prepare($get_contig_pos_query);
        $posh->execute;
        my $super_result = $posh->fetchrow_hashref;

	my $new_start = 0;
	my $new_stop = 0;
	if($minimum_gap_fg)
	{
		$new_start = $this_row->{start} + $super_result->{modified_contig_start_base};
                $new_stop  = $this_row->{stop} +  $super_result->{modified_contig_start_base};
	} else
	{
        	$new_start = $this_row->{start} + $super_result->{contig_start_super_base};
	        $new_stop  = $this_row->{stop} +  $super_result->{contig_start_super_base};
	}

	print ORFDBSCFILE	"supercontig_" . $super_result->{super_id} 	. "\t" .
			$this_row->{source} 	. "\t" .
			'ORF'		    	. "\t" .
			$new_start		. "\t" .
			$new_stop		. "\t" .
                        '.'                     . "\t" .
			$this_row->{direction}	. "\t" .
                        '.'                     . "\t" .
			$attributes		. "\n";
        print ORFDBFILE $this_row->{contig}     . "\t" .
                        $this_row->{source}     . "\t" .
                        'ORF'                   . "\t" .
                        $this_row->{start}       . "\t" .
                        $this_row->{stop}       . "\t" .
                        '.'                     . "\t" .
                        $this_row->{direction}  . "\t" .
                        '.'                     . "\t" .
                        $attributes             . "\n";

	$posh->finish;

  }



}


if($create_sage_from_db)
{
 my $query = '
        select tagid,
	id,
        contig,
        start,
        stop,
        direction
        FROM
        tagmap where contig is not null';
  my $result = $dbh->prepare($query);
  $result->execute;


  open(SAGEDBSCFILE, ">$gff_supercontig_sage_output_file");
  open(SAGEDBFILE, ">$gff_sage_output_file");


  while(my $this_row = $result->fetchrow_hashref)
  {
	my $get_library_results_query = "select tagid, library, result from sage_results where tagid = '" . $this_row->{tagid} . "'";
	my $libh = $dbh->prepare($get_library_results_query);
	$libh->execute;
	my $library_results = '';

	while(my $lib_row = $libh->fetchrow_hashref)
	{
		$library_results = $library_results . " ; library_" . $lib_row->{library} . " " . $lib_row->{result};
	}

	# Now determine if we have mapped ourselves to an orf, and if so what kind of mapping it is
	my $sageorfq = "select orfid, tagid, tagmapid, tagtype, unique_genome_fg, unique_trans_fg from orftosage where tagmapid = '" . $this_row->{id} . "' AND tagid = '" . $this_row->{tagid} . "'";
	my $sage_orfh = $dbh->prepare($sageorfq);
	$sage_orfh->execute();

	my $sageorf_results = '';
	if($sage_orfh->rows > 0)
	{
		my $sageorf_row = $sage_orfh->fetchrow_hashref;
		$sageorf_results = ' Orf ' . $sageorf_row->{orfid} . ' ; TagType "' . $sageorf_row->{tagtype} . '" ; UniqueGenome ' . $sageorf_row->{unique_genome_fg} .
					' ; UniqueTranscript ' . $sageorf_row->{unique_trans_fg} ;
	}

        my $get_contig_pos_query = "select contig_number, super_id,
                        contig_start_super_base, contig_length, modified_contig_start_base
                        from links where concat('contig_', contig_number) = '" . $this_row->{contig} . "'";
        my $posh = $dbh->prepare($get_contig_pos_query);
        $posh->execute;
        my $super_result = $posh->fetchrow_hashref;
	my $new_start = 0;
	my $new_stop = 0;
	if($minimum_gap_fg)
	{
                $new_start = $this_row->{start} + $super_result->{modified_contig_start_base};
                $new_stop  = $this_row->{stop} +  $super_result->{modified_contig_start_base};
	} else
	{
        	$new_start = $this_row->{start} + $super_result->{contig_start_super_base};
	        $new_stop  = $this_row->{stop} +  $super_result->{contig_start_super_base};
	}

        print SAGEDBSCFILE       "supercontig_" . $super_result->{super_id}      . "\t" .
                        'sage'	     		. "\t" .
                        'sagetag'               . "\t" .
                        $new_start              . "\t" .
                        $new_stop               . "\t" .
                        '.'                     . "\t" .
                        $this_row->{direction}  . "\t" .
                        '.'                     . "\t" .
                        'sagetag ' . $this_row->{tagid} . $library_results . ' ; ' . $sageorf_results      . "\n";
        print SAGEDBFILE $this_row->{contig}     . "\t" .
                        'sage'		         . "\t" .
                        'sagetag'                   . "\t" .
                        $this_row->{start}       . "\t" .
                        $this_row->{stop}       . "\t" .
                        '.'                     . "\t" .
                        $this_row->{direction}  . "\t" .
                        '.'                     . "\t" .
                        'sagetag ' . $this_row->{tagid} . $library_results . ' ; ' . $sageorf_results       . "\n";

        $posh->finish;
  }

}

if($load_db)
{

  # Load the main gff files
  
  # To make things a bit faster, concat all the files into one large temp file
  # so we can do a bulk load rather then an incremental load

  # First process organismdb  
  system("cat $gff_output_file $gff_orf_output_file $gff_sage_output_file $gff_repeat_file > $temp_file_name");

  # Load the tempfile into organismdb
  system("$database_bulk_loader --database $db_database_name --fasta $fasta_bases_file $temp_file_name");

  # Next process organismsc
  if($minimum_gap_fg)
  {
    system("$database_bulk_loader --database $db_supercontig_database_name --fasta $supercontig_fasta_file $gff_supercontig_output_file");

  } else
  {
    system("$database_bulk_loader --database $db_supercontig_database_name $gff_supercontig_output_file");
  }


  # Next process organismscreads

  system("cat $gff_supercontig_read_output_file $gff_supercontig_orf_output_file $gff_supercontig_sage_output_file  $gff_supercontig_repeat_file > $temp_file_name");

  if($minimum_gap_fg)
  {
    system("$database_bulk_loader --database $db_supercontig_read_database_name --fasta $supercontig_fasta_file $temp_file_name");
  } else
  {
    system("$database_bulk_loader --database $db_supercontig_read_database_name $temp_file_name"); 
  }



}

$dbh->disconnect;
# end
