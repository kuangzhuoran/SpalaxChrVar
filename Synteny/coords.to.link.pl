#!/usr/bin/perl -w
use strict;
#explanation:this program is edited to
#edit by hewm;   Thu Jan  6 15:50:07 CST 2022
#Version 1.0    hewm@genomics.org.cn

die  "Version 1.0\t2022-01-06;\nUsage: $0 <InPut><Out>\n" unless (@ARGV ==2);

##############Befor  Start  , open the files ####################

open (IA,"$ARGV[0]") || die "input file can't open $!";

open (OA,">$ARGV[1]") || die "output file can't open $!" ;

################# Do what you want to do #######################
my $MinAlgLen=5000;
<IA>;
<IA>;
<IA>;
<IA>;
<IA>;
	while(<IA>)
	{
		chomp;
		my @inf=split;
	   next if  ($inf[6]<90);
		print OA "$inf[-2]\t$inf[0]\t$inf[1]\t$inf[-1]\t$inf[3]\t$inf[4]\n";
		
#               next if  ($inf[9]<$MinAlgLen);
	}