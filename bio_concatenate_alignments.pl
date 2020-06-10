#!/usr/bin/perl
use Getopt::Std;
use File::Basename;

my $Version = "2020/6/9";
my $Program = basename $0;
my $Author = "Yan-Bo Sun (Yunnan University)";
my %opts;
getopts('hi:n:r:s:S:l:N:',\%opts);

die "
Program for concatenating FASTA alignments into super alignments
Contact: $Author
Version: $Version

Usage: perl $Program -i <FASTA Folder Name> [Other Options]

Options:
	-n INT   Number of FASTA files to concatenate [All]
	-r INT   Number of replicates [1]
	         (work with -n; Each replicate will generate a super
	         alignment containing -n input FASTA files)
	-s FLOAT Identity Cutoff (min.) for input FASTA alignments [0]
	         (FASTA with lower Identity will not be used)
	-S FLOAT Identity Cutoff (max.) for input FASTA alignments [1]
	         (FASTA with higher Identity will not be used)
	-l INT   Length Cutoff for input FASTA alignments [0]
	         (FASTA with shorter length will not be used)
	-N INT   Cutoff of Sequence Number in an input file [0]
	         (FASTA with less sequences not be used)
	-R       Sampling with Replacement [No]
	         (Only works with -n; Default = without Replacement)
	-h       Print this help information

Example:
	perl $Program -i in_alignments -s 0.85 -l 150
" if(!$opts{i} || $opts{h});

my $inFolder = $opts{i};#Filename should end with '.fas'
my @inFiles = glob "$inFolder/*.fa*"; 
my $numFile = defined $opts{n}?$opts{n}:0;
my $numRep = defined $opts{r}?$opts{r}:1;
my $cutSim = defined $opts{s}?$opts{s}:0;
my $cutSim2 = defined $opts{S}?$opts{S}:1;
my $cutLen = defined $opts{l}?$opts{l}:0;
my $cutSeq = defined $opts{N}?$opts{N}:0;
my $SamplingFlag = defined $opts{R}?1:0;
my $runLog = $inFolder.".concate.log";
open LOG,">$runLog";

#Print Parameters:
print STDERR "[Parameters]:\n";
print STDERR "    Number of files per replicate = $numFile\n";
print STDERR "    Number of replicate = $numRep\n";
print STDERR "    Cutoff of min.similarity = $cutSim\n";
print STDERR "    Cutoff of max.similarity = $cutSim2\n";
print STDERR "    Cutoff of alignment length = $cutLen\n";
print STDERR "    Cutoff of Sequence Number = $cutSeq\n";

#Perform concatenate replicates:
print STDERR "[Concatenate]:\n";
foreach my $tmpRep (1..$numRep){
	if ($numFile == 0){#All Files:
		my %hashIDSeq;
		foreach my $tmpFile (@inFiles){
			chomp $tmpFile;
			print STDERR "\r    [Rep-$tmpRep]: $tmpFile    ";
			my (@idList,@seqList,$tmpSeq,$tmpNum,$tmpLen,$meanIdentiry);
			open handle,"$tmpFile";
			while(<handle>){
				chomp;
				if(/\>(\S+)/){
					push (@idList,$1);
					push (@seqList,$tmpSeq) if ($tmpSeq ne "");
					$tmpSeq = "";
				}else{
					$tmpSeq .= $_;
				}
			}
			push (@seqList,$tmpSeq);
			close handle;
			
			if($cutSeq > 0){
				$tmpNum = scalar (@idList);
				if ($tmpNum < $cutSeq){
					print LOG "[Failed]\t$tmpFile contains $tmpNum seqs (< $cutSeq)\n" ; next;
				}
			}
			if($cutLen > 0){
				$tmpLen = length $seqList[0];
				if ($tmpLen < $cutLen){
					print LOG "[Failed]\t$tmpFile contains shorter seqs (< $cutLen)\n"; next;
				}
			}
			if($cutSim > 0 or $cutSim2 < 1){
				my ($identity,$replicate);
				foreach my $i (0..(@seqList-2)){
					foreach my $j (($i+1)..(@seqList-1)){
						$replicate++;
						$identity += &sub_pair_identity($seqList[$i],$seqList[$j]);
					}
				}
				$meanIdentiry = $replicate > 0?($identity/$replicate):0;
				if ($meanIdentiry < $cutSim){
					print LOG "[Failed]\t$tmpFile contains poor seqs (sim < $cutSim)\n"; next;
				}elsif($meanIdentiry > $cutSim2){
					print LOG "[Failed]\t$tmpFile contains high similar seqs (sim > $cutSim2)\n"; next;
				}
			}
			print LOG "[PASS]\t$tmpFile\t$meanIdentiry\t$tmpLen\t$tmpNum\n";
			foreach my $i (0..$#idList)
			{$hashIDSeq{$idList[$i]}.=$seqList[$i];}
		}
		#save to file:
		open OUT,">$inFolder\.concate\.all\.fas";
		while (my ($tmpID,$tmpSeq) = each %hashIDSeq){
			print OUT ">$tmpID\n$tmpSeq\n";
		}
		close OUT;
		print LOG "Saved to $inFolder\.concate\.all\.fas\n";
		print LOG "*"x50,"\n";
		last;
	}else{
		my (%hashIDList,%hashIDSeq);
		my $numFileReaded = 0;
		while($numFileReaded < $numFile){
			my $tmpIdx = rand($#inFiles);
			my $tmpFile = $inFiles[$tmpIdx];
			if($SamplingFlag == 0){
				if (exists $hashIDList{$tmpFile}){
					next;
				}else{
					$hashIDList{$tmpFile} = 0
				}
			}
			print STDERR "\r    [Rep-$tmpRep]: $numFileReaded ($tmpFile)    ";
			my (@idList,@seqList,$tmpSeq,$tmpNum,$tmpLen,$meanIdentiry);
			open handle,"$tmpFile";
			while(<handle>){
				chomp;
				if(/\>(\S+)/){
					push (@idList,$1);
					push (@seqList,$tmpSeq) if ($tmpSeq ne "");
					$tmpSeq = "";
				}else{
					$tmpSeq .= $_;
				}
			}
			push (@seqList,$tmpSeq) if ($tmpSeq ne "");
			close handle; #print STDERR "@seqList\n";
			
			if($cutSeq > 0){
				$tmpNum = scalar (@idList);
				if ($tmpNum < $cutSeq){
					print LOG "[Failed]\t$tmpFile contains $tmpNum seqs (< $cutSeq)\n" ; next;
				}
			}
			if($cutLen > 0){
				$tmpLen = length $seqList[0];
				if ($tmpLen < $cutLen){
					print LOG "[Failed]\t$tmpFile contains shorter seqs (< $cutLen)\n"; next;
				}
			}
			if($cutSim > 0 or $cutSim2 < 1){
				my ($identity,$replicate);
				foreach my $i (0..(@seqList-2)){
					foreach my $j (($i+1)..(@seqList-1)){
						$replicate++;
						my $tmp_value = &sub_pair_identity($seqList[$i],$seqList[$j]);
						$identity += $tmp_value;
					}
				}
				$meanIdentiry = $replicate > 0?($identity/$replicate):0;
				if ($meanIdentiry < $cutSim){
					print LOG "[Failed]\t$tmpFile contains bad seqs (sim < $cutSim)\n"; next;
				}elsif($meanIdentiry > $cutSim2){
					print LOG "[Failed]\t$tmpFile contains high similar seqs (sim > $cutSim2)\n"; next;
				}
			}
			$numFileReaded ++;
			print LOG "[PASS]\t$tmpFile\t$meanIdentiry\t$tmpLen\t$tmpNum\n";
			foreach my $i (0..$#idList)
			{$hashIDSeq{$idList[$i]}.=$seqList[$i];}
		}
		#save to file:
		open OUT,">$inFolder\.concate\.rep$tmpRep\.$numFile\.fas";
		while (my ($tmpID,$tmpSeq) = each %hashIDSeq){
			print OUT ">$tmpID\n$tmpSeq\n";
		}
		close OUT;
		print LOG "Saved to $inFolder\.concate\.rep$tmpRep\.$numFile\.fas\n";
		print LOG "*"x50,"\n";
		print STDERR "\n";
	}
}
print STDERR "[Done!]\n";
close LOG;

######################################################
sub sub_pair_identity{
	my ($seq1,$seq2) = @_;
	my $seqLength = length($seq1)-1;
	my ($identity,$blockLength);
	foreach my $i (0..$seqLength){
		my $base1 = uc substr($seq1,$i,1);
		my $base2 = uc substr($seq2,$i,1);
		if(($base1 ne '-') and ($base2 ne '-')){
			$blockLength++;
			$identity++ if ($base1 eq $base2);
		}
	}
	return $blockLength==0?0:scalar($identity/$blockLength);
}