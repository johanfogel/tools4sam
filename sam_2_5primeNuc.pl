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
-f  file [required]
-m	min length (18)
-x	max length (30)
-d  print head
-s  sample name (otherwise equals filename)
-o  output file

-f and -s can be a comma (,) separated list
";

getopts("m:x:s:f:o:hd",\%opts);


die $usage if $opts{'h'};




#take the filenames as samplenames unless specified
$opts{'s'}=$opts{'f'} unless $opts{'s'};

my @files=split(/,/,$opts{'f'});
my @samples=split(/,/,$opts{'s'});

my %tot_reads;
my @length=($opts{'m'},$opts{'x'});

my $skipped;

my $data;




for(my $i=0;$i<@files;$i++){
    my $file=$files[$i];
    my $sample=$samples[$i];
    
    my $si;
    #check the input
    if($file eq 'stdin'){
        $si=*STDIN;
    }elsif(-e $file){
        print STDERR "working with $file\n";

        if($file=~/bam$/){
            #its a bamfile
            open(IN, "samtools view $file |");
            $si=*IN;
        }elsif($file=~/sam$/){
            #its a samfile
            open(IN, "samtools view -S $file |");
            $si=*IN;
        }else{
            die"$file is nor a sam or bam file\n";
        }
    }else{
        die"file $file was not found.\n";
    }


    while(<$si>){
        chomp;
        #assume input is collapsed using FASTX-Toolkit
        my @row=split;
        my $seq=$row[9];
        my $flag=$row[1];
        my $count=1;
        $count=1/$1 if /NH:i:(\d+)/;
        
        #reverse the sequence if 2nd strand
        if($flag & 16){
            $seq=reverse($seq);
            $seq=tr/ATCGatcg/TACGtacg/;
        }
        my $len=length($seq);
        

        
        $skipped->{$sample}->{'long'}+=$count and next if $len>$length[-1];
        $skipped->{$sample}->{'short'}+=$count and next if $len<$length[0];
        
        my $five=substr($seq,0,1);
        
         #DEBUGG
        #print"$flag\t$count\t$five\t$len\t$sample\n$seq\n";
        #print"$_\n";
        $data->{$len}->{$five}->{$sample}+=$count;	

        $tot_reads{$sample}+=$count;

    }
    $si->close;
}


#report the number of skipped entries
print STDERR "Number of skipped sequences\n";
foreach my $sample(sort(keys(%{$skipped}))){
	foreach my $cat(sort(keys(%{$skipped->{$sample}}))){
		print STDERR "$sample\t$cat\t",$skipped->{$sample}->{$cat},"\n";
	}
}



my $openFile=open(my $fh ,"> $opts{'o'}");
$fh=*STDOUT unless $openFile;


print $fh "length\tbase\tsample\tCount\n" if $opts{'d'};
foreach my $l(sort(keys(%{$data}))){
	foreach my $m(sort(keys(%{$data->{$l}}))){
		foreach my $sample(sort(keys(%{$data->{$l}->{$m}}))){
			print $fh "$l\t$m\t$sample\t",$data->{$l}->{$m}->{$sample},"\n";
		}
	}
	
}
close $fh if $openFile;

