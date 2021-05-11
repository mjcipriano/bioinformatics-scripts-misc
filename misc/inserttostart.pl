#!/usr/bin/perl

use Bio::Seq;
use Bio::Tools::Run::RemoteBlast;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO::Writer::HTMLResultWriter;
use Bio::Tools::BPlite::Sbjct;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Root::IO;
use DBI;
use Bio::DB::GFF;



$ENV{PATH} = '/bin:/usr/bin:/usr/local/bin:/var/bio/blast/bin:/var/bio/bin:/usr/local/bio/bin:/usr/local/bio/blast';
BEGIN {$ENV{BLASTDIR} = '/usr/local/bio/blast';}
BEGIN {$ENV{BLASTDATADIR} = '/blastdb/'; }
   
# Connect to the database
   
use strict;
   
my $driver = "mysql";
my $hostname = "mib.mbl.edu";
my $port = "3306";
my $user = "gid";
my $password = "NOPE";
my $database = "giardia04";
   
my $dsn = "DBI:$driver:database=$database;host=$hostname;port=$port";
my $dbh = DBI->connect($dsn, $user, $password);
   
my $drh = DBI->install_driver("mysql");
 




my $hits_to_store = 20;
 
 
my $e_val= '1e-40';
my $prog = 'blastx';
 
 
my $writerhtml = new Bio::SearchIO::Writer::HTMLResultWriter();
 
 
my $query = " select s.id, s.idname, s.sequence_type_id, s.db_id, s.algorithm_id, s.sequence, s.translate, a.name as algorithm_name, db.name as database_name from sequence_search s, algorithms a, db
WHERE db.id = s.db_id
AND a.id = s.algorithm_id
LIMIT 1
";
 
my $sth = $dbh->prepare($query);
 
my $current_query = "insert into current_search (idname, sequence_type_id, db_id, algorithm_id) VALUES (?, ?, ?, ?)";
my $inshand = $dbh->prepare($current_query);
 
my $delete_query = "delete from sequence_search where id = ?";
my $delhand = $dbh->prepare($delete_query);
 
my $del_result_query = "delete from blast_results
                                where idname = ?
                                AND sequence_type_id = ?
                                AND db=?
                                AND algorithm = ?";
 
my $delrshand = $dbh->prepare($del_result_query);
 
my $del_report_query = "delete from blast_report_full
                                where idname = ?
                                AND sequence_type_id = ?
                                AND db_id=?
                                AND algorithm_id = ?";
 
my $delrphand = $dbh->prepare($del_report_query);
my $insert_query = "insert into blast_results (
                        idname,
                        sequence_type_id,
                        score,
                        evalue,
                        read_start,
                        read_end,
                        hit_start,
                        hit_end,
                        hit_name,
                        accession_number,
                        description,
                        algorithm,
                        db,
                        gaps,
                        frac_identical,
                        frac_conserved,
                        query_string,
                        hit_string,
                        homology_string,
                        hsp_rank,
                        hsp_strand,
                        hsp_frame
                       )
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
my $insert_handler = $dbh->prepare($insert_query);
                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                     
my $insert_full_string = "insert into blast_report_full
                          (idname,
                           sequence_type_ip->query->end),
                           db_id,
                           algorithm_id,
                           report)
                           VALUES (?, ?, ?, ?, ?)";
my $insert_full_handler = $dbh->prepare($insert_full_string);
                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                     
my $delhand_c_query = "delete from current_search where idname = ? AND sequence_type_id = ? AND db_id = ? AND algorithm_id = ?";
my $delhand_c = $dbh->prepare($delhand_c_query);
                                                                                                                                                                                                                                                     
my $keep_going = 1;




my $insert_search_query = 'insert into sequence_search (idname, sequence_type_id, db_id, algorithm_id, sequence, translate) VALUES (?, ?, ?, ?, ?, ?)';
my $insert_search_h = $dbh->prepare($insert_search_query);

my $bases_file = 'blastme.fas';

# Load all of the intergenic spaces into a sequence object;
my $sequences = Bio::SeqIO->new('-file'         => "$bases_file",
                                '-format'       => "fasta");

print "Starting check\n";
open(OUT, ">", "outfile.tab");

while( my $spacer = $sequences->next_seq)
{
	# Now grab the sequence 
	$insert_search_h->execute($spacer->id, 3, 2, 2, $spacer->seq(), 'N');

}




