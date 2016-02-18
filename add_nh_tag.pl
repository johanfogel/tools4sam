#!/usr/bin/perl

use strict;

############
# Add NH-tag in SAM-file, based on the number of reported entries in the bam-file
#
# A bit memory intensive, but requires no specific sorting...
my @reads;
my $readCount;

while(<>){
	print and next if /^@/;
	chomp;
	my @r=split;
	

	$readCount->{$r[0]}++;

	push(@reads,$_);
}

foreach (@reads){
	my @r=split;
	#add the NH count
	my $nh="NH:i:".$readCount->{$r[0]};
	print "$_\t$nh\n";
}



