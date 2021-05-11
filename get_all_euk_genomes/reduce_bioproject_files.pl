#!/usr/bin/perl



use strict;

my @all_lines = ();
my %prjid2line;
my %taxon_projects;


while(<>)
{
        my $line = $_;
        chomp($line);
        my @vals = split("\t", $line);

	my $prjid = $vals[2];
	my $taxon_id = $vals[1];
	my $taxon_name = $vals[0];
	my $num_proteins = $vals[10];
	my $num_scaff = $vals[11];
	my $size_scaff = $vals[12];

	if($num_proteins < 100)
	{
		next;
	}
        push(@all_lines, $line);
        $prjid2line{$prjid} = $line;
	if(!exists $taxon_projects{$taxon_id})
	{
		$taxon_projects{$taxon_id} = ();
	}
        push(@{$taxon_projects{$taxon_id}}, $prjid);
        
        # Does it have any protein
}


foreach my $taxonid ( sort {$a <=> $b} keys %taxon_projects )
{
	my @projects = @{$taxon_projects{$taxonid}};
	
	if(scalar @projects == 1)
	{
		print $prjid2line{$projects[0]} . "\n";
	} else
	{
		# is one of them refseq?
		my $check_refseq = 0;
		foreach my $prjid (@projects)
		{
			my @fields = split("\t", $prjid2line{$prjid});
			my $rec_type = $fields[5];
			if($rec_type eq 'RefSeq Genome')
			{
				print $prjid2line{$prjid} . "\n";
				$check_refseq = 1;
				last; 
			}
				
		}
	}
#        print join(":", @{$taxon_projects{$taxonid}}) . "\n";
}

