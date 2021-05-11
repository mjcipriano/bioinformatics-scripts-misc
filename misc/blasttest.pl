#!/usr/bin/perl

BEGIN {$ENV{BLASTDIR} = ' /usr/local/bio/blast/'; }
BEGIN {$ENV{PATH} = '/bin:/usr/bin:/usr/local/bin:/var/bio/blast/bin:/var/bio/bin:/usr/local/bio/blast'; }
use Getopt::Long;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Tools::BPlite::Sbjct;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Root::IO;
use strict;

my $dbs = '/blastdb/nt';
my $dir = "/tmp";
my @params = (     'database'    => "$dbs",
		   'program'     => "blastn"
		   );

my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
my $input = Bio::Seq->new( -id  => 'my_query',
                           -seq =>'TGGATTGGGCCTATCTATTCAAATAATTGCAATATTGAAAATAAT');


my $blast_report = $factory->blastall($input);

my $this_result = $blast_report->next_result;

while( my $hit = $this_result->next_hit()) 
{ 
   print "\thit name: ", $hit->name(), " significance: ", $hit->significance(), "\n";

  # Get the HSP
  my $this_hsp = $hit->hsp;
  print "Start: " . $this_hsp->hit->start . "\t" . "End:" . $this_hsp->hit->end . "\n"; 
}


