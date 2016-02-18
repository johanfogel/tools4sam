#!/usr/bin/env perl

use strict;
use Bio::SeqIO;

use Getopt::Std;

###################
# make calc 5' preference. Make it a ggplot compatible file

my %opts;
$opts{'m'}=18;
$opts{'x'}=30;



my $usage="
$0 [opts]
-f  file (read from stdin unless specified)
-m	min length (18)
-x	max length (30)
-d  print head
-s  sample name (otherwise equals filename)
-o  output file
-c  assume collapsed reads from fastx_collapser

-f and -s can be a comma (,) separated list
";

getopts("m:x:s:f:o:hdc",\%opts);


die $usage if $opts{'h'};


$opts{'f'}='stdin' unless $opts{'f'};

#take the filenames as samplenames unless specified
$opts{'s'}=$opts{'f'} unless $opts{'s'};

my @files=split(/,/,$opts{'f'});
my @samples=split(/,/,$opts{'s'});

my %tot_reads;
my @length=($opts{'m'},$opts{'x'});

my $skipped;

my $data;

my $si;


for(my $i=0;$i<@files;$i++){
    my $file=$files[$i];
    my $sample=$samples[$i];
    
    #check the input
    if($file eq 'stdin'){
        print STDERR "working with STDIN\n";
        $si=Bio::SeqIO->new(-fh=>\*STDIN ,-format=>'fasta');
    }elsif(-e $file){
        print STDERR "working with $file\n";
        if($file=~/gz$/){
            $si=Bio::SeqIO->new(-file=>"gunzip -c $file | ",-format=>'fasta');
        }else{
            $si=Bio::SeqIO->new(-file=>"< $file",-format=>'fasta');
        }
    }else{
        die"file $file was not found.\n";
    }


    while(my $seq=$si->next_seq){
        my $count=1;
        
        if($opts{'c'}){
            #assume input is collapsed using FASTX-Toolkit
            (my $id,$count)=split(/-/,$seq->id);
        }
        my $len=length($seq->seq);
        
        $skipped->{$sample}->{'long'}+=$count and next if $len>$length[-1];
        $skipped->{$sample}->{'short'}+=$count and next if $len<$length[0];
        
        my $five=substr($seq->seq,0,1);
        
        
        $data->{$len}->{$five}->{$sample}+=$count;	

        $tot_reads{$sample}+=$count;

    }
    $si->close;
}


#report the number of skipped entries
print STDERR "Number of skipped sequences\n";
foreach my $sample(sort(keys(%{$skipped}))){
	foreach my $cat(sort(keys(%{$skipped->{$sample}}))){
		print STDOUT "$sample\t$cat\t",$skipped->{$sample}->{$cat},"\n";
	}
}



my $openFile=open(my $fh ,"> $opts{'o'}");
$fh=*STDOUT unless $openFile;


print $fh "length\tbase\tsample\count\n" if $opts{'d'};
foreach my $l(sort(keys(%{$data}))){
	foreach my $m(sort(keys(%{$data->{$l}}))){
		foreach my $sample(sort(keys(%{$data->{$l}->{$m}}))){
			print $fh "$l\t$m\t$sample\t",$data->{$l}->{$m}->{$sample},"\n";
		}
	}
	
}
close $fh if $openFile;



