#!/usr/bin/perl


use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;
use DBI;
use XML::DOM;

use strict;

my %config;

my $organism = 'giardia03';

my $organism_name = 'giardia';
my $species_name = 'giardia lamblia';
my $release_date = '10/26/03';
my $release_version = '0.3';

# Use to move orf data between versions
my $old_orf_database = 'giardia02'; 

# Mysql database settings
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

# Change below to reflect your settings

my $base_dir = "/var/bio/ARACHNE";
my $working_dir = "gl101603";
my $run_dir = "run1";

my $root_web_dir = '/var/www/html/';
my $web_php_dir = 'bdb';

my $web_dir = $root_web_dir . $organism . '/';
my $download_dir = $web_dir . 'download/';
my $gbrowse_conf_dir = '/var/www/gbrowse/beta/';

my $root_dir = $base_dir . '/' . $working_dir . '/';
my $xml_file = $root_dir . "/traceinfo/reads.xml";
my $assembly_unplaced_file = $root_dir . '/' . $run_dir . '/' . 'assembly.unplaced';
my $assembly_links_file = $root_dir . '/' . $run_dir . '/' . 'assembly.links';
my $assembly_reads_file = $root_dir . '/' . $run_dir . '/' . 'assembly.reads';
my $fasta_bases_file = $root_dir . '/' . $run_dir . '/' . 'assembly.bases';
my $fasta_reads_bases_file =  $root_dir . '/' . 'reads.fasta';
my $fasta_reads_quality_file =  $root_dir . '/' . 'reads.qual';


# Try not to change these
my $gff_contig_dir = "load/$organism/contig/";
my $gff_supercontig_dir = "load/$organism/supercontig/";
my $gff_output_file = $gff_contig_dir . "contig.gff";
my $gff_supercontig_read_output_file = $gff_supercontig_dir . "supercontig_read.gff";
my $gff_supercontig_orf_output_file = $gff_supercontig_dir . "orf_supercontig.gff";
my $gff_orf_output_file = $gff_contig_dir . "orf.gff";
my $gff_supercontig_sage_output_file =  $gff_supercontig_dir . "supercontig_sage.gff";
my $gff_sage_output_file =  $gff_contig_dir . "sage.gff";
my $gff_repeat_file = $gff_contig_dir . "repeat.gff";
my $gff_supercontig_repeat_file = $gff_supercontig_dir . "supercontig_repeat.gff";
my $gff_trna_file = $gff_contig_dir . "trna.gff";
my $gff_supercontig_trna_file = $gff_supercontig_dir . "supercontig_trna.gff";
my $gff_coverage_file = $gff_contig_dir . "coverage.gff";
my $gff_supercontig_coverage_file = $gff_supercontig_dir . "supercontig_coverage.gff";

my $modified_reads_fasta_file = "load/" . $organism . "/trimed_reads.fasta";
my $supercontig_fasta_file = 'load/' . $organism . "/supercontig/supercontig.fasta";


# These below can be changed

my $orf_input_file = "gff/" . $organism . "_orfs.gff"; # Used when you want to convert a orf gff file from contig view to supercontig view

my $database_bulk_loader = "../bp_bulk_load_gff.pl";
my $database_incremental_loader = "../bp_load_gff.pl";
my $ace_dir = $root_dir . '/' . $run_dir . '/' . 'acefiles';
my $sage_dir = '/xraid/licor2/gsage/sageanalysis/results/';
my $sage_file = '10103.alllib.tags';

my $blast_bin_dir = '/var/bio/blast/bin/';
my $blast_db_dir = '/var/bio/blast/data/';
my $emboss_db_dir = '/biodb/emboss/';
my $repeatFinder_bin = '/usr/local/bin/repeatFinder';
my $trnascan_bin = '/usr/local/bin/tRNAscan-SE';

# ORGANISM.fasta file should exist in the tables directory which contains known genes to be used to train glimmer and Codon Preference with
my $fasta_train_file = "tables/$organism" . ".fasta";
my $glimmer_train_output = "tables/$organism/glimmer_out.bin";
my $cusp_codon_file = "tables/$organism/cusp_out.cod";
my $glimmer_bin_dir = '../glimmer/';
my $build_icm_bin = $glimmer_bin_dir . "build-icm";
my $glimmer2_bin = $glimmer_bin_dir . "glimmer2";
my $emboss_dir = "/usr/local/bin/";
my $cai_bin = $emboss_dir . "cai";
my $cusp_bin = $emboss_dir . "cusp";
my $chips_bin = $emboss_dir . "chips";
my $mummer_bin = "/var/bio/src/MUMmer3.08/mummer";

my $debug = 0;
# Turn off/on parts of this script

my $drop_schema				= 0;


my $create_schema			= 0;
my $create_directories			= 0;
my $parse_xml      			= 0;
my $parse_unplaced 			= 0;
my $parse_links    			= 0;
my $parse_reads    			= 0;
my $parse_reads_bases 			= 0;
my $parse_reads_qual 			= 0;
my $create_modified_fasta		= 0; # will create a supercontig fasta file to use against the supercontig view
my $create_modified_reads_fasta		= 0; # will create a fasta file of used reads with the unused potions of the reads cut out
my $move_orfs_from_old			= 0; # will move orfs over from an older version of the database and search in the new assembly
my $find_orfs_glimmer			= 0; # will use glimmer to search for new orfs in this assembly
my $delete_invalid_orfs			= 0; # this will check for duplications and invalid orfs and mark them as deleted
my $parse_ace	   			= 0; # Not finished - will parse an ace file and extract information on read locations
my $create_supercontig_orf_from_input 	= 0; # will read in a tab delimited file and create a gff file, not used any more
my $create_sage_from_file 		= 0; # to import sage data from a tab delimited file
my $run_orf_tests 			= 0; # will run tests on all orfs to give P/F or scores to determine if they are real orfs
my $map_sage_to_db			= 0;
my $map_sage_to_db_mummer		= 0;
my $map_sage_to_orf_secondary 		= 0; # this will dump data into a temp table to be used for next part
my $map_sage_to_orf_tert		= 0; # this will insert into orftosage for determined tag->orf mappings
my $find_repeat				= 0;
my $find_trna				= 0;
my $calculate_coverage			= 0; # Calculate the read overlap and create a gff file for a graph	
my $create_read_gff			= 0;
my $create_supercontig_orf_from_db      = 0;
my $create_sage_from_db 		= 1; # this will create sage.gff files

my $create_search_db			= 0; # Will create blast and emboss format databases
my $create_download_files		= 0;
my $create_web_files			= 0;
my $load_db        			= 0; # this will load the database with this organisms gff files

my $minimum_gap_fg 			= 1;
my $minimum_gap_length			= 100;

my $minimum_orf_length			= 150;
my $maximum_orf_length			= 18000;
my $remove_orf_false_start		= 1;
my $remove_orf_no_stop			= 1;

my $transcript_tail = 15;


# Process command line options


my $config_file = '';

#GetOptions('organism'=> \$organism);
#GetOptions('organism_name'=> \$organism_name);

GetOptions("configuration=s" => \$config_file);
print $config_file;
# Process Configuration File

if($config_file ne '')
{
	open(CONFIG, $config_file) or die("Can not open Configuration File: $config_file\n");

	while(<CONFIG>)
	{
		my $line = $_;
		chomp($line);
		# First check if this is a comment line
		if($line =~ m/^#/)
		{
			next;
		}
		# Now check to make sure this is a correct configuration line (it must have an equal sign and something before and after the equal sign

		if(!($line =~ /\w+\=\w+/))
		{
			next;
		}

		# Now process the line
		my ($variable, $value) = split("=", $line);
		$($variable) = $value;
	}
}
print $organism . "\n";

die;

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
	print "End: Dropping Schema\n";
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



if($create_directories)
{
        print "Start: Creating Directories\n";
	system("mkdir load");
	system("mkdir load/$organism");
        system("mkdir load/$organism/contig");
        system("mkdir load/$organism/supercontig");
	system("mkdir tables");
	system("mkdir tables/$organism");
        print "End: Creating Directories\n";

}

if($parse_xml)
{
	
	print "Start: Parse $xml_file file\n";

	# instantiate parser
	my $xp = new XML::DOM::Parser();
	
	# parse and create tree
	my $doc = $xp->parsefile($xml_file) or die ("Can not open $xml_file file!");
	
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
					if($item->getFirstChild())
					{
	                                	$template_id = $item->getFirstChild()->getData;
					}
	                        }elsif (lc($item->getNodeName) eq "library_id")
	                        {
					if($item->getFirstChild())
					{
	                                	$library_id = $item->getFirstChild()->getData;
					}
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
	                $dbh->quote($trace_name)  . ', ' . 
	                $dbh->quote($center_name) . ', ' . 
	                $dbh->quote($plate_id)    . ', ' . 
	                $dbh->quote($well_id)     . ', ' .
	                $dbh->quote($template_id) . ', ' .
	                $dbh->quote($library_id)  . ', ' .
	                $dbh->quote($trace_end)   . ', ' .
	                $dbh->quote($trace_direction) .  
	               ')';
	         my $sth = $dbh->prepare($query);
	         $sth->execute;
	         $sth->finish;
		}
	}
        print "End: Parse $xml_file file\n";

}	
#######DONE WITH READS.XML FILE###############

#######START ASSEMBLY.UNPLACED FILE#############

if($parse_unplaced)
{
	print "Start: Parse $assembly_unplaced_file file\n";

	open(ASSM, "$assembly_unplaced_file") or die ("Can not open $assembly_unplaced_file");
	my $update_status_query = "UPDATE reads SET status = ? WHERE read_name = ?";
        my $sth = $dbh->prepare($update_status_query);

	while (<ASSM>)
	{
	    my $line = $_;
	    my @this_line = split(" ", $line);

	    $sth->execute($this_line[1], $this_line[0]);
	    $sth->finish;
	
	}
	close(ASSM);
        print "End: Parse $assembly_unplaced_file file\n";
}
#######DONE ASSEMBLY.UNPLACED FILE#############

#######START ASSEMBLY.LINKS FILE#############

if($parse_links)
{
	print "Start: Parse $assembly_links_file file\n";
	open(LINKS, "$assembly_links_file") or die ("Can not open $assembly_links_file");
	
	#super_id
	#num_bases_in_super
	#num_contigs_in_super
	#ordinal_num_of_contig_in_supercontig_id
	#contig_number
	#length_of_contig
	#estimated_gap_before_contig
	#estimated_gap_after_contig
        my $update_links_query = "insert into links
                                (super_id, bases_in_super, contigs_in_super, ordinal_number,
                                contig_number, contig_length, gap_before_contig, gap_after_contig)
                                VALUES (?, ?, ?, ?, ?, ?, ?, ?)";
        my $sth = $dbh->prepare($update_links_query);
	
	while (<LINKS>)
	{
	    my $line = $_;
	    my @this_line = split(" ", $line);
	
	    if($this_line[0] ne '#super_id')
	    {	
		    $sth->execute($this_line[0], $this_line[1], $this_line[2], $this_line[3], $this_line[4], $this_line[5], $this_line[6], $this_line[7]);
	    }
	
	}

	$sth->finish;
	close(LINKS);
	print "Completed parsing $assembly_links_file file\n";
                                                                                                                                                                                   
	# Now determine where contigs are within the supercontig with or without minimum gap length 
                                                                                                                                                                                    
        my $query = '
                select distinct
                super_id
                FROM
                links';
        my $links_result = $dbh->prepare($query);
        $links_result->execute;

        my $contq = "select contig_number, contig_length, ordinal_number, gap_before_contig from links where super_id = ? ORDER BY ordinal_number";
        my $conth = $dbh->prepare($contq);

        my $updq = "update links set
                    modified_contig_start_base = ?
                    where contig_number = ?";
	my $updh = $dbh->prepare($updq);

        my $updqbases = "update links set
                    modified_bases_in_super = ?
                    where super_id = ?";
	my $updbases = $dbh->prepare($updqbases);
                                                                                                                                                                             
                                                                                                                                                                                    
        while(my $links_array = $links_result->fetchrow_hashref)
        {
                                                                                                                                                                                    
                # We need to find out the total modified supercontig length And update where the contig starts in the supercontig
                my $start_super_base = 1;
                my $running_end = 0;

                $conth->execute($links_array->{super_id});
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

                        $updh->execute($this_start, $this_contig->{contig_number});

                        $running_end = $this_start + $this_contig->{contig_length};
                }

                        $updbases->execute($running_end, $links_array->{super_id});
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
	print "Start: Parse $assembly_reads_file file\n";
	open(READS, "$assembly_reads_file") or die ("Can not open $assembly_reads_file");

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
                                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
	my $sth = $dbh->prepare($update_reads_query);
	while (<READS>)
	{
	    my $line = $_;
            chomp($line);
            $line =~ s/[\r\n]+$//;
                                                                                                                                                                                                                                                    
            my @this_line = split("\t", $line);
            while(scalar @this_line < 17)
            {
                push @this_line, '';
            }

            foreach my $this_var (@this_line) 
	    {
		if($this_var eq "")
		{
		  $this_var = undef;
		}

	    }

	    $sth->execute( $this_line[0],$this_line[1],$this_line[2],$this_line[3],$this_line[4],$this_line[5],$this_line[6],$this_line[7],$this_line[8],$this_line[9],$this_line[10],$this_line[11],
				$this_line[12],$this_line[13],$this_line[14],$this_line[15],$this_line[16]);

	}

	$sth->finish;
	close(READS);
        print "End: Parse $assembly_reads_file file\n";

}

#######END ASSEMBLY.READS FILE#############


#######START READS.BASES##############
if($parse_reads_bases)
{

        print "Start: Parse $fasta_reads_bases_file file\n";

	my $in  = Bio::SeqIO->new('-file' => "$fasta_reads_bases_file",
                         '-format' => 'Fasta');

	my $insert_query = "insert into reads_bases (read_name, bases) VALUES (?, ?)";
        my $sth = $dbh->prepare($insert_query);

	while ( my $seq = $in->next_seq() )
	{
		my $sth = $dbh->prepare($insert_query);
		$sth->execute($seq->id, $seq->seq);
	}
	$sth->finish;
        print "End: Parse $fasta_reads_bases_file file\n";

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
if($create_modified_reads_fasta)
{
        print "Start: Create modified reads fasta file\n";

	my $sequences = Bio::SeqIO->new('-file'         => "<$fasta_reads_bases_file",
					'-format'       => "fasta");
	open(READFASTA, ">", $modified_reads_fasta_file);
        while(my $this_read =  $sequences->next_seq)
	{
		my $disp_id = $this_read->display_id;
		$disp_id =~ s/ //g;
		my $query = 'select read_name, first_base_of_trim, read_len_trim from reads_assembly where read_name = ' . $dbh->quote($disp_id);
		my $rh = $dbh->prepare($query);
		$rh->execute();
		if(my $result = $rh->fetchrow_hashref)
		{
			my $start = $result->{first_base_of_trim} + 1;
			my $end  = $result->{read_len_trim} + 1 + $result->{first_base_of_trim};
			if($end > $this_read->length())
			{
				$end = $this_read->length();
			}
			print READFASTA ">" . $result->{read_name} . "\n";
			print READFASTA $this_read->subseq($start, $end) . "\n";
		}
		
	}

        print "End: Create modified reads fasta file\n";
}
if($move_orfs_from_old)
{
        print "Start: Moving orfs from database: $old_orf_database to database: $organism\n";

	my $old_database = $old_orf_database;
	my $new_database = $organism;
  
	my $odsn = "DBI:$driver:database=$old_database;host=$hostname;port=$port";
	my $odbh = DBI->connect($odsn, $user, $password);
  
	my $odrh = DBI->install_driver("mysql");
 
	my $ndsn = "DBI:$driver:database=$new_database;host=$hostname;port=$port";
	my $ndbh = DBI->connect($ndsn, $user, $password);
	my $ndrh = DBI->install_driver("mysql");
 
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
	 
	        my $this_orf_sequence = Bio::Seq->new ( -display_id     => $this_row->{ORFid},
	                                                -seq            => $sequence_string_db);
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
	                #       print "\t" . $seq->display_id . "\n";
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
	                        (ORFid, 
				sequence, 
				annotation, 
				annotation_type,
				source,
	                        delete_fg, 
				delete_reason, 
				old_orf, 
				attributes, 
				orf_name) VALUES (" .
	                        "'" . $this_row->{ORFid} . "'," .
	                        "'" . $this_orf_sequence->seq . "'," .
	                        $ndbh->quote($this_row->{annotation}) . "," .
	                        "'" . $this_row->{annotation_type} . "'," .
	                        "'" . $this_row->{source} . "'," .
	                        "'Y'," .
	                        "'disconnected', 
				'Y'," .
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
	                my $delete_reason = $this_row->{delete_reason};
	                if($this_row->{delete_fg} eq "no")
	                {
	                        $delete_flag = "N";
	                        $delete_reason = "";
	                } else
	                {
	                        $delete_flag = "Y";
	                }
	                my $insert_query = "insert into orfs
	                        (ORFid, 
				sequence, 
				annotation, 
				annotation_type,
				source,
	                        delete_fg, 
				delete_reason, 
				contig, 
				start, 
				stop, 
				direction, 
				attributes, 
				orf_name) VALUES (" .
	                        "'" . $this_row->{ORFid} . "'," .
	                        "'" . $this_orf_sequence->seq . "'," .
	                        $ndbh->quote($this_row->{annotation}) . "," .
	                        "'" . $this_row->{annotation_type} . "'," .
	                        "'" . $this_row->{source} . "'," .
	                        "'" . 'N' . "'," .
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

        print "End: Moving orfs from database: $old_database to database: $organism\n";

}

#######START FIND ORFS VIA GLIMMER############

if($find_orfs_glimmer)
{

        print "Start: Discover ORFs via GLIMMER2\n";

	system("$build_icm_bin -f < $fasta_train_file > $glimmer_train_output");

        my $sequences = Bio::SeqIO->new('-file'         => "<$fasta_bases_file",
                                        '-format'       => "fasta");

	my $tmp_file = 'temp/' . $organism . '_tmp.fas'; 
        while (my $seq = $sequences->next_seq)
        {
		open(TMPFILE, ">", $tmp_file);
	
		print TMPFILE ">" . $seq->display_id . "\n" . $seq->seq();
		close(TMPFILE);
		system($glimmer2_bin . ' ' . $tmp_file . ' ' . $glimmer_train_output . "-g $minimum_orf_length -X -o 0 -p 0 > temp/" . $organism . '_glimmer.out');

		open (ORFILE, "<", 'temp/' . $organism . '_glimmer.out');
		my $chk_query = "select orfid from orfs
                                  WHERE
                                  contig = ? AND
                                  start = ? AND
                                  stop = ? AND
                                  direction = ? AND
                                  (start % 3) = ?";
                my $chk_h = $dbh->prepare($chk_query);

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
				$chk_h->execute( $seq->display_id, $start, $stop, $dir, $start%3);

				if($chk_h->rows() > 0)
				{
					# do nothing
				} elsif( ($stop - $start) < $minimum_orf_length)
				{
					# do nothing
				} elsif( ($stop - $start) > $maximum_orf_length)
				{
					# do nothing
				} elsif(substr($orf_seq,0,3) ne 'ATG' &&  $remove_orf_false_start)
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

        print "Start: Delete invalid Orfs\n";

	# First delete all called orfs that are too large (> 20000)
        print "       Deleting Orfs of invalid size\n";

	my $delh = $dbh->prepare("update orfs set delete_fg = 'Y', delete_reason='invalid size' where ((stop-start > " . $maximum_orf_length  . ") OR (stop-start < " . $minimum_orf_length . ")) AND delete_fg = 'N'");
	$delh->execute();

	# Now check for false starts and stop codons within reading frame
        print "       Deleting Orfs with false stops and stop codons within reading frame\n";

	my $get_orf_query = "SELECT orfid,
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
	                        where sequence is not NULL and contig is not null ";
	my $sth = $dbh->prepare($get_orf_query);
	$sth->execute;
	while(my $this_row = $sth->fetchrow_hashref)
	{
	        # Bring the sequence into a bioperl object
	        my $this_orf_sequence = Bio::Seq->new ( -display_id     => $this_row->{ORFid},
	                                                -seq            => $this_row->{sequence}
	                                                );
 
	        if($this_orf_sequence->subseq(1,3) ne 'ATG' && $remove_orf_false_start)
	        {
	                #print "Orf\t" . $this_row->{orfid} . " \t is NOT a valid translating sequence.\tReason:\t" . "false start\n";
	                my $update_query = "UPDATE orfs set delete_fg='Y', delete_reason='false start'
	                                WHERE orfid = '" . $this_row->{orfid} . "'";
	                my $udq = $dbh->prepare($update_query);
	                $udq->execute;
	        } else
	        {
	                # Now check for a stop in middle
	                my $start = 1;
	                my $valid = 1;
	                my $not_valid_reason = '';
	                my $seq_length = $this_orf_sequence->length();
	                while( ($start <= ($seq_length -3)) && ($valid) )
	                {
	                        my $check_codon = $this_orf_sequence->subseq($start, $start+2);
				#print $check_codon . "\n";
	                        if( ($check_codon eq 'TAA') ||
	                            ($check_codon eq 'TAG') ||
	                            ($check_codon eq 'TGA') )
	                        {
	                                #print "STOP";
	                                # We have a stop codon in the middle of our sequence
	                                $valid = 0;
	                                $not_valid_reason = 'mid stop';
	                        }
	                        $start = $start+3;
	                }
	                                                                                                                                                                                                                                                      
	                # check that it ends in a stop codon
	                my $end_seq =  substr($this_row->{sequence}, -3);
	                                                                                                                                                                                                                                                      
	                if($valid && $remove_orf_no_stop)
	                {
	                        if( ($end_seq eq 'TAA') ||
	                            ($end_seq eq 'TAG') ||
	                            ($end_seq eq 'TGA') )
	                        {
	                                # We are good or we already have invalidated the sequence
	                        }else
	                        {
	                                # we are missing a stop codon, so not a valid orf
	                                $not_valid_reason = 'no stop';
	                                $valid = 0;
	                        }
	                }
	                                                                                                                                                                                                                                                      
	                # Get mark orfs that are < minimum_orf_length or > maximum_orf_length as invalid size
	                if($valid)
	                {
	                        my $size = length($this_row->{sequence});
	                        if( ($size < $minimum_orf_length) || ($size > $maximum_orf_length) )
	                        {
	                                $not_valid_reason = 'invalid size';
	                                $valid = 0;
	                        }
	                }
	                # if it is still valid, we are ok, otherwise print out the reason
	                if($valid)
	                {
	                        #print "Orf\t" . $this_row->{orfid} . " \t is a valid translating sequence.\n";
	                } else
	                {
	                         #print "Orf\t" . $this_row->{orfid} .
	                        #" \t is NOT a valid translating sequence.\tReason:\t" .
	                        #$not_valid_reason . "\n";
	                        my $update_query = "UPDATE orfs set delete_fg='Y',
	                                delete_reason='" . $not_valid_reason . "'
	                                WHERE orfid = '" . $this_row->{orfid} . "'";
	                        #print $update_query . "\n";
	                        my $udq = $dbh->prepare($update_query);
	                        $udq->execute;
	                } # end if $valid
	        } # end else ne 'ATG'
	} # end while
                                                                                                                                                                                                                                                      

	# Now check for and set delete_fg for orfs that are duplicates of other orfs
	# or are entierly within another orf. If one of them is an old orf
	# keep that one
	 
	# Now we will mark all duplicate orfs. We will always keep the lower numbered orf
	# as the correct orf, and we will insert into the orf_reassign table to keep track of
	# these orfs
        print "       Deleting Duplicated Orfs or Orfs internal to other Orfs\n";

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
	my $update_query = "UPDATE orfs set delete_fg = 'Y',
		delete_reason = 'within other'
		WHERE
		contig = ? AND
		start >= ? AND
		stop <= ? AND
		direction = ? AND
		(start % 3) = ? AND
		delete_fg = 'N' AND
		orfid != ?";
	my $updh = $dbh->prepare($update_query);
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
                                orfid > ? order by orfid LIMIT 1";
                #print $orf_query . "\n\n";
	my $sth = $dbh->prepare($orf_query);

 
	while($last)
	{
		if($debug >= 5)
		{
	        	print $orf_query . "\n\n";
		}
	        $sth->execute($last_orf_id);;
	 
	        my $this_row = $sth->fetchrow_hashref;
	 
	        if($this_row)
	        {
	                $last_orf_id = $this_row->{orfid};

	                $updh->execute($this_row->{contig}, $this_row->{start}, $this_row->{stop}, $this_row->{direction}, $this_row->{start}%3,  $this_row->{orfid});
			if($debug >= 5)
			{
		                print $update_query;
		                print "\n\n";
			}
	                print "ORF\t" . $this_row->{orfid} . "\t within: " . $updh->rows;
	                print "\n";
	        }else
	        {
	                $last = 0;
	        }
	 
	}# end if ($last)

        print "End: Delete invalid Orfs\n";

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
        print "Start: Convert Contig Orf gff file to Supercontig Orf gff file\n";

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

        print "Start: Convert Contig Orf gff file to Supercontig Orf gff file\n";
}

if($create_sage_from_file)
{
        print "Start: Parsing Sage Input File\n";

  open(SAGEFILE, "$sage_dir" . "$sage_file") or die("Can not open sage file $sage_file");

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

    my $insert_tag_query = "insert into sage_tags (tagID, sequence)
                                VALUES (?, ?)";
    my $insert_tagh = $dbh->prepare($insert_tag_query);
    my $insert_result_query = "insert into sage_results (tagID, library, result)
                                VALUES (?, ?, ?)";
    my $insert_resulth = $dbh->prepare($insert_result_query);


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

      $insert_tagh->execute($sage_array[0], $sage_array[1]);
  
      # Start at -2 so we can skip the first two values (tagID, sequence)
      my $count = -2;
   
      # iterate through the remaining libraries
      foreach my $lib (@sage_array)
      {
        if($count >= 0)
        {
          $insert_resulth->execute( $sage_array[0], $library_names[$count], $lib);
	}
        $count++;			
      }
    }
			
  }
        print "End: Parsing Sage Input File\n";

}
if($map_sage_to_db_mummer)
{

        # Delete old tagmap
        $dbh->prepare('delete from tagmap')->execute();

	# First create a temporary fasta file to be used with mummer
	my $temp_sage_fasta = "temp/$organism.sage.fasta";

	open(SAGEFASTA, ">", $temp_sage_fasta);

        my $get_tag_query = "select tagid, sequence from sage_tags order by tagid";
        my $sth = $dbh->prepare($get_tag_query);
        $sth->execute;
	
	while(my $this_row = $sth->fetchrow_hashref)
	{
		print SAGEFASTA ">" . $this_row->{tagid} . "\n" . $this_row->{sequence} . "\n";

	}
	close(SAGEFASTA);

	# Now run mummer
	system($mummer_bin . " -b -c $temp_sage_fasta $fasta_bases_file > temp/$organism.mummer.out");

	# Now parse the output

	my $dir = "+";
	my $contig = "";

	open(MUMMEROUT, "temp/$organism.mummer.out");

	my $insert_query = "insert into tagmap
				(tagID, contig, start, stop, direction
				) VALUES ( ?, ?, ?, ?, ?)";
                               my $insert_ref = $dbh->prepare($insert_query);
	
	while(<MUMMEROUT>)
	{
		my $line = $_;
		if($line =~ /^>/) # This is a header line for a new contig
		{
			($contig) = $line =~ /^> ([\w\.]+)/;
			if($line =~ /Reverse$/)
			{
				$dir = "-";
			} else
			{
				$dir = "+";
			}
		} else # This is a match result line
		{
			my $stop = '';
			my (undef, $tagid, undef, $start, $size) = split(/\s+/, $line);
			if($dir eq "-")
			{
				$stop = $start;
				$start = $stop - $size;
			} else
			{
				$stop = $start + $size;
			}
			$insert_ref->execute($tagid, $contig, $start, $stop, $dir);
		
		}
	}


}
if($map_sage_to_db)
{

        print "Start: Mapping Sage Tags to Genome\n";

	# Iterate through all old sage tags

	# Delete old tagmap
	$dbh->prepare('delete from tagmap')->execute();
 
	my $get_tag_query = "select tagid, sequence from sage_tags order by tagid";
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
	        my @hit_array;
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
	 
	                                print   join("\t", "TAG ",
	                                $this_row->{tagid},
	                                " matches " , $seq->display_id,
	                                " location:", $this_index,
	                                (length($this_tag_sequence->seq) + $this_index), $dir) . "\n";
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
	                } # End while this db row sequence reversed has checked a particular contig
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
	                        $dbh->quote($this_row->{tagid}) . ', ' .
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
	                       $dbh->quote($this_row->{tagid}) . "," .
	                       $dbh->quote($this_hit->[1]) . "," .
	                       $dbh->quote($this_hit->[2]) . "," .
	                       $dbh->quote($this_hit->[3]) . "," .
	                       $dbh->quote($this_hit->[4]) . ")";
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

	# create the codon usage table
	system("$cusp_bin -sequence $fasta_train_file -outfile $cusp_codon_file");

	# Get all orfs in the database and their sequence
	my $query = "
        	select orfid,
		sequence,
	        contig,
	        start,
	        stop,
	        direction
	        FROM
	        orfs where sequence != '' AND sequence is not null";
	my $sth = $dbh->prepare($query);
	$sth->execute;


        my $udquery =   "UPDATE orfs set TestCode = ? ,
                         TestScore = ? ,
                         GeneScan = ? ,
                         GeneScanScore = ? ,
                         CodonUsage = ? ,
                         CodonPreferenceScore = ? ,
                         CodonPreference = ?  
                         WHERE orfid = ?";
                                                                                                                                                                                                                                                       
	my $udh = $dbh->prepare($udquery);

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




		my $cusage_program = "$chips_bin -seqall temp/tmp.seq -outfile temp/outfile.tmp";
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

		my $cai_program = "$cai_bin -seqall temp/tmp.seq -outfile temp/outfile.tmp -cfile $cusp_codon_file -sbegin1 1";
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
                my $cai_program = "$cai_bin -seqall temp/tmp.seq -outfile temp/outfile.tmp -cfile $cusp_codon_file -sbegin1 2";
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
                my $cai_program = "$cai_bin -seqall temp/tmp.seq -outfile temp/outfile.tmp -cfile $cusp_codon_file -sbegin1 3";
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
		
                $udh->execute($test_result, $testcode_average, $gene_result, $genescan_average, $codon_usage_score, $cai_result_1, $cai_test, $this_row->{orfid});


	}


}

sub find_sequence_in_fasta {

	my $sequence = shift;
	my $name = shift;
	my @hit_array;
	my $sequences = Bio::SeqIO->new('-file'         => "<$fasta_bases_file",
        				'-format'       => "fasta");
	# change sequence to only ATGC
	$sequence =~ s/[^ATGC]//ig;
	my $sequence_obj = Bio::Seq->new ( -display_id     => $name,
                                           -seq            => $sequence);
	my $sequence_to_check = lc($sequence_obj->seq);
	my $sequence_to_check_rc = lc($sequence_obj->revcom()->seq);
	
	while(my $seq = $sequences->next_seq)
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
	                                                                                                                                                                                                                            
				print   "ID \t"
				. $name
				. " matches \t" . $seq->display_id
				. " location:\t" . $this_index . "\t"
				.  (length($sequence_obj->seq) + $this_index) . "\t" . $dir . "\n";
				push @hit_array, [$name, $seq->display_id, $this_index+1, (length($sequence_obj->seq) + $this_index), $dir];
				$last_index = $this_index;
				if($checked_first == 0)
				{
					$this_index = index(lc($seq->seq), $sequence_to_check, $last_index+1);
				} else
				{
					$this_index = index(lc($seq->seq), $sequence_to_check_rc, $last_index+1);
				}
			}
		} # END while checked first

		$checked_first++;
	} # End while this db row sequence has checked a particular contig

	return \@hit_array;
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

if($find_trna)
{
	system($trnascan_bin . ' -o temp/' . $organism . '.trna ' . $fasta_bases_file);
	open(TRNAFILE, "temp/" . $organism . '.trna');
	open(TRNAGFF, ">", $gff_trna_file);
        open(TRNASUPER, ">", $gff_supercontig_trna_file);
	my $count = 1;
	while(<TRNAFILE>)
	{
		# Skip the first 3 lines
		if($count >= 4)
		{
			my $line = $_;
			chomp($line);
			my ($contig, $ordinal, $start, $stop, $type, $anti_codon, $intron_start, $intron_end, $score) = split("\t", $line);
			my $dir = '';

			if($start > $stop)
			{
				$dir = '-';
			} else
			{
				$dir = '+';
			}

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
                                                                                                                                                                                                                                                      
                        print TRNASUPER      join("\t", "supercontig_" . $super_result->{super_id}, 'tRNAscan-SE-1.23', 'tRNA', $new_start, $new_stop ,$score,$dir ,'0',"tRNA $type ; anticodon $anti_codon" ) . "\n";
                        print TRNAGFF   join("\t", $contig, 'tRNAscan-SE-1.23', 'tRNA', $start, $stop ,$score,$dir ,'0',"tRNA $type ; anticodon $anti_codon" ) . "\n";
                                                                                                                                                                                                                                                      
                        $posh->finish;

		}
		$count++;
	}


	
}

if($calculate_coverage)
{

	open(COVER, ">", $gff_coverage_file);
	open(SUPERCOVER, ">", $gff_supercontig_coverage_file);

	# Iterate through each of the contigs, find out the overlap at each base
	my $sequences = Bio::SeqIO->new('-file'         => "<$fasta_bases_file",
					'-format'       => "fasta");

	my $check_query = 'select count(*) as num_overlap from reads_assembly 
			   where trim_read_in_contig_start <= ? 
			   AND trim_read_in_contig_stop >= ?
			   AND contig_number = ?';
	my $sth = $dbh->prepare($check_query);
	#my @result_array;

        my $get_contig_pos_query = "select contig_number, super_id,
                                    contig_start_super_base, contig_length, modified_contig_start_base
                                    from links where  contig_number = ?";

        my $posh = $dbh->prepare($get_contig_pos_query);

        while(my $seq = $sequences->next_seq)
        {
		my $length = $seq->length();
		my $pointer = 1;
		my ($contig_id) = $seq->display_id =~ /(\d+)/;
		my $last_coverage_val = 0;
		my $start_coverage_val = 1;
		my $stop_coverage_val = 1;
		
		# find out the starting coverage value
		$sth->execute($pointer, $pointer, $contig_id);
		my $result = $sth->fetchrow_hashref();
		my $last_coverage_val = $result->{num_overlap};
	
		while($pointer <= $length)
		{
			$sth->execute($pointer, $pointer, $contig_id);
			my $result = $sth->fetchrow_hashref();
			if( ($last_coverage_val == $result->{num_overlap}) && ($pointer != $length) )
			{
				$stop_coverage_val++;
			} else
			{
	                        $posh->execute($contig_id);
	                        my $super_result = $posh->fetchrow_hashref;


	                        my $new_start = 0;
	                        my $new_stop = 0;
	                        if($minimum_gap_fg)
	                        {
	                                $new_start = $start_coverage_val + $super_result->{modified_contig_start_base};
	                                $new_stop  = $stop_coverage_val +  $super_result->{modified_contig_start_base};
	                        } else
	                        {
	                                $new_start = $start_coverage_val + $super_result->{contig_start_super_base};
	                                $new_stop  = $stop_coverage_val +  $super_result->{contig_start_super_base};
	                        }

				
				#push(@result_array, ($contig_id, $start_coverage_val, $stop_coverage_val, $last_coverage_val));
				print COVER join("\t", $contig_id, 'count', 'coverage', $start_coverage_val, $stop_coverage_val, $last_coverage_val, '.', '.') . "\n";
                                print SUPERCOVER join("\t", $contig_id, 'count', 'coverage', $new_start, $new_stop, $last_coverage_val, '.', '.') . "\n";
				$last_coverage_val =  $result->{num_overlap};
				$start_coverage_val = $pointer;
				$stop_coverage_val = $pointer;

			}
			$pointer++;
		}

		
	}


}
if($create_search_db)
{

	# We need to create a fasta file first, then call the formatdb and dbifasta command on it

	system("mkdir $emboss_db_dir/$organism");

	# Contigs first, this is the same as our assembly.bases
	system('cp ' . $fasta_bases_file . ' ' . $blast_db_dir . '/' . $organism );
        system($blast_bin_dir . 'formatdb -t ' . $organism . ' -i ' . $blast_db_dir . '/' . $organism . ' -p F -o T');
	system("mkdir $emboss_db_dir/$organism/contigs");
        system("cp $fasta_bases_file $emboss_db_dir/$organism/contigs/$organism" );
	system("cd $emboss_db_dir/$organism/contigs/;dbifasta -idformat simple -directory . -filenames $organism -dbname $organism_name" . "_contigs -release $release_version -date $release_date");


	

       # Supercontigs for use with gbrowse
        system('cp ' . $supercontig_fasta_file . ' ' . $blast_db_dir . '/' . $organism . '_supercontig' );
        system($blast_bin_dir . 'formatdb -t ' . $organism . '_supercontig -i ' . $blast_db_dir . '/' . $organism . '_supercontig -p F -o T');

	# ORFs as nt sequences
	my $orf_nt_file = $blast_db_dir . '/' . $organism . '_orfs_nt';
        open(ORFFILE, '>', $orf_nt_file);
                                                                                                                                                                                  
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
        system("mkdir $emboss_db_dir/$organism/orfs_nt");
	system("cp $orf_nt_file $emboss_db_dir/$organism/orfs_nt/" );
        system("cd $emboss_db_dir/$organism/orfs_nt/;dbifasta -idformat simple -directory . -filenames $organism" . "_orfs_nt -dbname $organism_name" . "_orfs_nt -release $release_version -date $release_date");




	# Translated ORFs
	my $orf_aa_file = $blast_db_dir . '/' . $organism . '_orfs_aa';
	open(ORFFILE, '>', $orf_aa_file);

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
        system("mkdir $emboss_db_dir/$organism/orfs_aa");
        system("cp $orf_aa_file $emboss_db_dir/$organism/orfs_aa/" );
        system("cd $emboss_db_dir/$organism/orfs_aa/;dbifasta -idformat simple -directory . -filenames $organism" . "_orfs_aa -dbname $organism_name" . "_orfs_aa -release $release_version -date $release_date");


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

if($create_download_files)
{
	# Create directories
        system('mkdir ' . $web_dir);
        system('mkdir ' . $download_dir);

	# Create fasta file for orfs, translated orfs, used shotgun reads, unused reads, contigs
	
	#ORFs
	open(ORFFILE, ">", $download_dir . 'orfs.fas');

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
	
	system('gzip -f ' . $download_dir . 'orfs.fas');

        # Translated ORFs
                                                                                                                                                                                                                                                    
        open(ORFFILE, '>', $download_dir. 'orfs_aa.fas');
                                                                                                                                                                                                                                                    
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
        system('gzip -f ' . $download_dir . 'orfs_aa.fas');

        # Unused Reads
                                                                                                                                                                                                                                                    
        open(READFILE, '>', $download_dir . 'unused_reads.fas');
                                                                                                                                                                                                                                                    
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

        system('gzip -f ' . $download_dir . 'unused_reads.fas');

	# Used Reads
        open(READFILE, '>', $download_dir . 'used_reads.fas');
                                                                                                                                                                                                                                                    
        my $query = "select reads.read_name as readname, reads_assembly.read_name, reads_bases.bases
                        FROM
                        reads_bases,
                        reads left join reads_assembly ON reads.read_name=reads_assembly.read_name
                        WHERE reads.read_name = reads_bases.read_name AND reads_assembly.read_name is not null";
                                                                                                                                                                                                                                                    
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
                                                                                                                                                                                                                                                    
        system('gzip -f ' . $download_dir . 'used_reads.fas');

	# All Reads
        system('cp ' . $fasta_reads_bases_file . ' ' . $download_dir . 'all_reads.fas');
        system('gzip -f ' . $download_dir . 'all_reads.fas');

	# Contigs
	system('cp ' . $fasta_bases_file . ' ' . $download_dir . 'contigs.fas');
	system('gzip -f ' . $download_dir . 'contigs.fas');


}

if($create_web_files)
{
	# For each web file that must be created, open the file, perform a search of !ORGANISM! and replace with $organism
	system('cp -f template/html/blast.php ' . $web_dir);
        system('cp -f template/html/bottom.html ' . $web_dir);
        system('cp -f template/html/index.php ' . $web_dir);
        system('cp -f template/html/style.css ' . $web_dir);
        system('cp -f template/html/top_header.php ' . $web_dir);
        system('cp -f template/html/blast.html ' . $web_dir);
        system('cp -f template/html/top.html ' . $web_dir);
        system('cp -f template/contig.conf ' . $gbrowse_conf_dir . $organism . '.conf');
        system('cp -f template/supercontig.conf ' . $gbrowse_conf_dir . $organism . 'screads.conf');

	my @modify_array ;
	push @modify_array, $web_dir . 'blast.html';
        push @modify_array, $web_dir . 'top.html';
        push @modify_array, $gbrowse_conf_dir . $organism . 'screads.conf';
        push @modify_array, $gbrowse_conf_dir . $organism . '.conf';


	my $search = '!ORGANISM!';
	my $replace = $organism;

	foreach my $file (@modify_array)
	{
		print "editing $file\n";
		my @contents;
		open(FH, "+< $file")                 or die "Opening: $!";
		while (<FH>) {
			my $line = $_;
			$line =~ s/$search/$replace/g;
			push(@contents,$line);
		}

		# change ARRAY here
		seek(FH,0,0)                        or die "Seeking: $!";
		print FH @contents                    or die "Printing: $!";
		truncate(FH,tell(FH))               or die "Truncating: $!";
		close(FH)                           or die "Closing: $!";
	}

}
### BEGIN OUTPUT GFF

if($create_read_gff)
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
        }
        $sth->finish;


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
			links.contig_length,
			links.modified_contig_start_base
		      FROM reads,
			reads_assembly r_assem,
			links
		      WHERE reads.read_name = r_assem.read_name
			AND links.contig_number = r_assem.contig_number";
	my $sth = $dbh->prepare($gff_query);
	$sth->execute;

	# Query used in while loop to get information about a read's partner
                my $check_partner_super_query = "
		SELECT distinct
		reads.read_name,
		r_assem.contig_number,
		r_assem.read_pair_name,
		r_assem.read_pair_contig_number,
		links.super_id as super_contig_number,
		links.contig_start_super_base,
		links.contig_length,
		r_assem.orientation
		FROM reads,
		reads_assembly r_assem,
		links
		WHERE reads.read_name = r_assem.read_name
		AND links.contig_number = r_assem.contig_number
		AND reads.read_name = ?";
	my $inner_check = $dbh->prepare($check_partner_super_query);


	# Query used in while loop to check if a partner exists
	my $check_for_partner_query = 'select read_name from reads where template_id = ?  AND read_name != ?';
	my $check_for_partner_result = $dbh->prepare($check_for_partner_query);

	while(my $this_row = $sth->fetchrow_hashref)
	{
		my $read_type = '';

	       # Check if this read is special :
	       # 1- it's partner read is unplaced and not missing(it has a partner)
	       # 2- it's partner read is on another contig
	       # 3- it's partner read is on another supercontig
	       # 4- it's partner is unplaced and missing
	       # 5- it's partner is in the next contig and the gap is positive
	       # 6- it's partner is in the next contig and the gap is negative
	       # 7- it's partner is not in the next ordinal contig

                $inner_check->execute($this_row->{read_pair_name});
                my $inner_row = $inner_check->fetchrow_hashref;
		

		# Do I have a partner?
		if($this_row->{read_pair_contig_number} eq "") # there is no partner read
                {
                        $check_for_partner_result->execute($this_row->{template_id}, $this_row->{read_name});
                        my $num_rows_ret = $check_for_partner_result->rows;
                        if($num_rows_ret == 0)
                        {
                              $read_type =  'missing-partner';
                        } else
                        {
                              $read_type = 'unplaced-partner';
                        }
		} elsif($this_row->{contig_number} eq $this_row->{read_pair_contig_number}) #they are in the same contig
		{
			if($this_row->{orientation} eq $inner_row->{orientation})
			{
				$read_type = 'partner-different-contig-same-direction';
			} else
			{
	     			$read_type = 'read';
			}
		} elsif($inner_row->{super_contig_number} ne $this_row->{super_contig_number}) # They are in different supercontigs
		{
			$read_type = 'partner-different-supercontig';
		} else
	  	{
				# Now check if the two contigs overlap

				my $contig_one = $this_row->{contig_number};
				my $contig_two = $this_row->{read_pair_contig_number};
				my $read_one_dir = $this_row->{orientation};
				my $read_two_dir = $inner_row->{orientation};
		
				my $contig_one_start = $this_row->{contig_start_super_base};
				my $contig_two_start = $inner_row->{contig_start_super_base};
				my $contig_one_end = $this_row->{contig_start_super_base} + $this_row->{contig_length};
				my $contig_two_end =  $inner_row->{contig_start_super_base} + $inner_row->{contig_length};
	
			  	if( ($read_one_dir eq "+") && ($read_two_dir eq "-") )
			  	{
			    		if($contig_one_end <= $contig_two_start)  
			    		{
			      			$read_type = "partner-different-contig-positive-gap";
			    		} else
			    		{
	                      			$read_type = "partner-different-contig-negative-gap";
			    		}
			  	} elsif( ($read_one_dir eq "-") && ($read_two_dir eq "+") )
			  	{
	                    		if($contig_one_start >= $contig_two_end)
	                    		{
	                      			$read_type = "partner-different-contig-positive-gap";
	                    		} else
	                    		{
	                      			$read_type = "partner-different-contig-negative-gap";
	                    		}
			  	} elsif($this_row->{orientation} eq $inner_row->{orientation})
		                {
		                	$read_type = "partner-different-contig-same-direction";
		                } else
			  	{
			    		$read_type = "partner-different-contig-exception";
			  	}
	
			} 
	
		my $read_super_start_val = 0;
		my $read_super_stop_val = 0;
		if($minimum_gap_fg)
		{
			$read_super_stop_val = $this_row->{trim_read_in_contig_stop} + $this_row->{modified_contig_start_base};
	            	$read_super_start_val = $this_row->{modified_contig_start_base} + $this_row->{trim_read_in_contig_start};
		} else
		{
			$read_super_stop_val = $this_row->{trim_read_in_contig_stop} + $this_row->{contig_start_super_base};
			$read_super_start_val = $this_row->{contig_start_super_base} + $this_row->{trim_read_in_contig_start};
		}
	
	        print SUPERGFF  join("\t", "supercontig_" . $this_row->{super_contig_number}, 
                        	"read", 
				$read_type,
				$read_super_start_val, 
	                    	$read_super_stop_val,
	                    	'.',
	                    	$this_row->{orientation},
	                    	'.',
	                    	"Sequence \"" . $this_row->{read_name} . "\"" .
				' ; ReadPair "' . $this_row->{read_pair_name} . '" ; ReadPairContig "' .
	                    	$this_row->{read_pair_contig_number} . '" ; ReadContig "' . $this_row->{contig_number} .
                    		'" ; TemplateID "' . $this_row->{template_id} .
                    		'" ; GivenInsertSize "' . $this_row->{given_insert_size} .
                    		'" ; GivenInsertStdDev "' . $this_row->{given_insert_std_dev} .
                    		'" ; ObservedInsertSize "' . $this_row->{observed_insert_size} .
                    		'" ; ObservedInsertStdDev "' . $this_row->{observed_inserted_deviation} . '"') .
				"\n";

		print GFF join("\t", "contig_" . $this_row->{contig_number},
				"read",
				$read_type,
				$this_row->{trim_read_in_contig_start},
				$this_row->{trim_read_in_contig_stop}.
				'.',
				$this_row->{orientation},
				'.',
				"Sequence \"" . $this_row->{read_name} . "\"" .
				' ; ReadPair "' . $this_row->{read_pair_name} . '" ; ReadPairContig "' .
				$this_row->{read_pair_contig_number} . '" ; ReadContig "' . $this_row->{contig_number} .
				'" ; TemplateID "' . $this_row->{template_id} .
				'" ; GivenInsertSize "' . $this_row->{given_insert_size} .
				'" ; GivenInsertStdDev "' . $this_row->{given_insert_std_dev} .
				'" ; ObservedInsertSize "' . $this_row->{observed_insert_size} .
				'" ; ObservedInsertStdDev "' . $this_row->{observed_inserted_deviation} . '"') .
				"\n";
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
        my $blast_query = "select  hit_name, description, evalue from blast_results where sequence_type_id = 2 AND db = 2 AND algorithm = 3 AND idname = " . $dbh->quote($this_row->{orfid}) . "  AND (description not like '%ATCC 50803%' OR description like '%gb|%') order by evalue limit 1";
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
        my $pfam_query = "select  hit_name, description, evalue from blast_results where sequence_type_id = 2 AND db = 4 AND algorithm = 4 AND idname = " . $dbh->quote($this_row->{orfid}) . " AND evalue <= 1e-3 order by evalue";
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

	my $attributes = 'Orf "' . $this_row->{orfid} . 
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
			from links where concat('contig_', contig_number) = " . $dbh->quote($this_row->{contig}) ;
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

  # Get the percentages in each library also, first get the totals for each library

  	my $percent_query = 'select library, sum(result) as total_tags from sage_results group by library';
  	my $percent_h = $dbh->prepare($percent_query);
  	$percent_h->execute();

	# Now take each grand total value and put it into a hash
        my %total_hash;
        while(my $percent_row = $percent_h->fetchrow_hashref)
        {
                $total_hash{$percent_row->{library}} = $percent_row->{total_tags};
        }


  while(my $this_row = $result->fetchrow_hashref)
  {
	my $get_library_results_query = "select tagid, library, result from sage_results where tagid = '" . $this_row->{tagid} . "'";
	my $libh = $dbh->prepare($get_library_results_query);
	$libh->execute;
	my $library_results = '';


	while(my $lib_row = $libh->fetchrow_hashref)
	{
		$library_results = $library_results .   " ; library_" . $lib_row->{library} . " " . $lib_row->{result} . 
						     	" ; librarypercent_" . $lib_row->{library} . " " . $lib_row->{result}/$total_hash{$lib_row->{library}}*100 ;
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
  
  # First process organismdb  

  system("$database_bulk_loader --create --database $db_database_name $gff_output_file $gff_orf_output_file $gff_sage_output_file $gff_repeat_file $gff_trna_file $modified_reads_fasta_file");


  # Next process organismscreads

  if($minimum_gap_fg)
  {
    system("$database_bulk_loader --create --database $db_supercontig_read_database_name $gff_supercontig_read_output_file $gff_supercontig_orf_output_file $gff_supercontig_sage_output_file $gff_supercontig_repeat_file $gff_supercontig_trna_file  $modified_reads_fasta_file $supercontig_fasta_file");
  } else
  {
    system("$database_bulk_loader --create --database $db_supercontig_read_database_name $gff_supercontig_read_output_file $gff_supercontig_orf_output_file $gff_supercontig_sage_output_file $gff_supercontig_repeat_file $gff_supercontig_trna_file  $modified_reads_fasta_file ");
  }



}

$dbh->disconnect;
# end
