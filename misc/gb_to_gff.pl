#!/usr/bin/perl -w
use strict;

use Bio::Tools::GFF;
use Bio::SeqIO;
use XML::DOM;
use Bio::Seq;
use DBI;
 
my $organism = 't_brucei';
  
my $driver = "mysql";
my $hostname = "mib";
my $port = "3306";
my $user = "gid";
my $password = "NOPE";
my $database = $organism;
  
my $dsn = "DBI:$driver:database=$database;host=$hostname;port=$port";
  
# File/script locations
  
my $gff_file_name = $organism . '.gff';
my $fasta_file_name = $organism . '.fas';
 
my $db_database_name = $organism . "db";
my $db_supercontig_database_name = $organism . "sc";
my $db_supercontig_read_database_name = $organism . "screads";
  
my $dbh = DBI->connect($dsn, $user, $password);
                                                                                                                                                                                 
     
my $drh = DBI->install_driver("mysql");
 
 

my ($seqfile) = @ARGV;
die("must define a valid seqfile to read") unless ( defined $seqfile && -r $seqfile);

my $seqio = new Bio::SeqIO(-format => 'genbank',
			   -file   => $seqfile);
my $count = 0;
my $file_name = '';
while( my $seq = $seqio->next_seq ) {
    $count++;
    # defined a default name
    my $fname = sprintf("%s.gff", $seq->display_id || "seq-$count");
    $file_name = $fname;
    my $gffout = new Bio::Tools::GFF(-file => ">$fname" ,
				     -gff_version => 2);
    
    foreach my $feature ( $seq->top_SeqFeatures() ) {
	$gffout->write_feature($feature);
    }
}
open(GFF, $file_name);
my $ins_query = 'insert into orfs (
                                  orfid,
                                  orf_name,
                                  annotation,
                                  annotation_type,
                                  source,
                                  delete_fg,
                                  delete_reason,
                                  contig,
                                  start,
                                  stop,
                                  direction) VALUES (?,?,?,?,?,?,?,?,?,?,?)';
my $ins_h = $dbh->prepare($ins_query);

while(<GFF>)
{
	my $this_line = $_;
	chomp($this_line);
	my ($contig, $source,$type, $start, $stop, undef, $dir, undef, $desc) = split(/\t/, $this_line);
	my $name = '';
	my $annot = '';
	my @desc_array = split(/;/, $desc);
	if(scalar @desc_array >= 2)
	{
	foreach my $feat(@desc_array)
	{
		$feat =~ s/^\s+//g;
		$feat =~ s/\s+$//g;
		if( substr($feat, 0, 4) eq 'gene')
		{
			my $sstart = index($feat, '"');
			my $send = index($feat,'"', $sstart+1);
			$name = substr($feat, $sstart+1, $send-$sstart-1);
		} elsif( substr($feat, 0, 7) eq 'product')
		{
                        my $sstart = index($feat, '"');
                        my $send = rindex($feat,'"');
                        $annot = substr($feat, $sstart, $send-$sstart+2);
		}
	}

	if($name ne '')
	{
			
		print $name . "\t" . $annot . "\n";
		$ins_h->execute(undef, $name, $annot, 'Sanger', $source, 'N', undef, $contig, $start, $stop, $dir);
	}
	}

}
