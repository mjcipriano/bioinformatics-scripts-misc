#!/usr/bin/perl 

use strict;
use Bio::Tools::Run::RemoteBlast;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Tools::BPlite::Sbjct;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Root::IO;
use DBI;

$ENV{PATH} = '/bin:/usr/bin:/usr/local/bin:/var/bio/blast/bin:/var/bio/bin:/usr/local/bio/bin:/usr/local/bio/blast';
BEGIN {$ENV{BLASTDIR} = '/usr/local/bio/blast';}
BEGIN {$ENV{BLASTDATADIR} = '/blastdb/'; }


my $driver = "mysql";
my $hostname = "alien.mbl.edu";
my $port = "3306";
my $user = "gid";
my $password = "NOPE";
my $database = "giardia";
my $dsn = "DBI:$driver:database=$database;host=$hostname;port=$port";
my $dbh = DBI->connect($dsn, $user, $password);
my $drh = DBI->install_driver("mysql");

my $hit_to_store = 10;

my $db = '/blastdb/nt';
my $e_val= '.01';
my $prog = 'blastn';


#		'readmethod' => 'SearchIO' 
my @params = (  'p' => $prog,
		'd' => $db,
                'e' => $e_val
);

my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

my $query = "
select rb.read_name, MID(rb.bases, ra.first_base_of_trim+1, read_len_trim)
from reads_bases rb,
blast ra
where rb.read_name = ra.read_name
LIMIT 1
";
# print "Starting\n";
my $keep_going = 1;

while($keep_going)
{
  my $sth = $dbh->prepare($query);
  $sth->execute;
  if($sth->rows == 0)
  {
    	$keep_going = 0;
    	# Do Nothing
  } else
  {
    	my $this_row = $sth->fetchrow_arrayref;
	my $delete_query = "delete from blast where read_name = '" . $this_row->[0] . "'";
         print $delete_query . "\n";
	 print $this_row->[1] . "\n";
        my $delhand = $dbh->prepare($delete_query);
        #$delhand->execute;

    	my $input = Bio::Seq->new(-id=>$this_row->[0] , '-seq' => $this_row->[1] );
	#Blast a sequence against a database:
    	my $r = $factory->blastall($input);

	while ( my @rids = $factory->each_rid ) 
	{
		foreach my $rid ( @rids ) 
		{
			my $rc = $factory->retrieve_blast($rid);
			if( !ref($rc) ) 
			{
				if( $rc < 0 ) 
				{
					$factory->remove_rid($rid);
				}
				sleep 5;
			} else 
			{ 
				my $current_hit = 0;
				while( (my $result = $rc->next_result()) 
					&& ($current_hit < $hit_to_store) )
				{
					$current_hit++;
					$factory->remove_rid($rid);
					while(my $hit = $result->next_hit)
					{ 
						my $hsp = $hit->next_hsp;
						my $insert_query = "insert into blast_results (
				        		read_name,
					        	score,
							read_start,
							read_end,
							hit_start,
							hit_end,
							hit_name,
							accession_number,
							description )
							VALUES (" .
	                                                "'" . $result->query_name() . "'," .
							"'" . $hsp->score . "'," .
	                                                "'" . $hsp->query->start . "'," .
	                                                "'" . $hsp->query->end . "'," .
	                                                "'" . $hsp->hit->start . "'," .
	                                                "'" . $hsp->hit->end . "'," .
	                                                "'" . $hit->name . "'," .
	                                                "'" . $hit->accession . "'," .
	                                                "'" . $hit->description . "')";
  					my $insert_handler = $dbh->prepare($insert_query);
  					$insert_handler->execute;

					 print $insert_query; 

					} # END HSP		
				} # END CURRENT RESULT WHILE
			}  # END RESULT IF STATEMENT
		}  # END FOREACH RESULT ENTRY
	}  # END WHILE THERE IS A RESULT ENTRY
  }  # END ELSE  THERE IS SOMETHING LEFT TO DO

}  # END WHILE KEEP GOING
