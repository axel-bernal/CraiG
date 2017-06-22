# ids are sneu_ID, so we get rid of that sneu at the beginning
perl -ne 'chomp(); if($_ =~ /^>sneu_(\S+)\s*/) { print "\n>$1\n";} else { print $_; }' | sed 1,1d | sed -e '$a\'