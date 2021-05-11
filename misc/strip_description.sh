perl -e 'open(T, "orfs_aa.fas"); while(<T>) { if($_ =~ m/^>/) { ($orfid) = $_ =~ /^>(\d+)\ \|/; print ">" . $orfid . "\n";} else { print $_;} }' > orfs_aa_stripped.fas
