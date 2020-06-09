#!/usr/bin/perl
use Getopt::Std;
use File::Basename;

my $Version = "2020/6/9";
my $Program = basename $0;
my $Author = "Yan-Bo Sun (Yunnan University)";
my %opts;
getopts('hi:f:r:s:l:n:',\%opts);

die "
Program for concatenating FASTA alignments into super alignments
Contact: $Author
Version: $Version

Usage: perl $Program -i <FASTA Folder Name> [Other Options]

Options:
	-f INT   Number of FASTA files to concatenate [All]
	-r INT   Number of replicates [1]
	         (work with -n; Each replicate will generate a super
	         alignment containing -n input FASTA files)
	-s FLOAT Identity Cutoff for input FASTA alignments [NULL]
	         (FASTA with lower Identity will not be used)
	-l INT   Length Cutoff for input FASTA alignments [NULL]
	         (FASTA with shorter length will not be used)
	-n INT   Cutoff of Sequence Number in an input file [NULL]
	         (FASTA with less sequences not be used)
	-h       Print this help information

Example:
	perl $Program -i in_alignments -s 0.85 -l 150
" if(!$opts{i} || $opts{h});

my $inFolder = $opts{i};#Filename should end with '.fas'
my @inFiles = glob "$inFolder/*.fa*"; 
my $numFile = defined $opts{f}?$opts{f}:0;
my $numRep = defined $opts{r}?$opts{r}:1;
my $cutSim = defined $opts{s}?$opts{s}:0;
my $cutLen = defined $opts{l}?$opts{l}:0;
my $cutSeq = defined $opts{n}?$opts{n}:0;
my $runLog = $inFolder.".concate.log";
open LOG,">$runLog";

#Print Parameters:
print STDERR "[Parameters]:\n";
print STDERR "    IDs to concatenate = $idFile\n";
print STDERR "    Number of files to concatenate = $numFile\n";
print STDERR "    Number of replicate = $numRep\n";
print STDERR "    Cutoff of alignment similarity = $cutSim\n";
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
			if($cutSim > 0){
				my ($identity,$replicate);
				foreach my $i (0..(@seqList-2)){
					foreach my $j (($i+1)..(@seqList-1)){
						$replicate++;
						$identity += &sub_pair_identity($seqList[$i],$seqList[$j]);
					}
				}
				$meanIdentiry = $replicate > 0?($identity/$replicate):0;
				if ($meanIdentiry < $cutSim){
					print LOG "[Failed]\t$tmpFile contains bad seqs (sim < $cutSim)\n"; next;
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
		my %hashIDSeq;
		my $numFileReaded = 0;
		while($numFileReaded < $numFile){
			my $tmpIdx = rand($#inFiles);
			my $tmpFile = $inFiles[$tmpIdx];
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
					print LOG "[Failed]\t$tmpFile contains $tmpNum seqs (< $cutSeq)\n"; next;
				}
			}
			if($cutLen > 0){
				$tmpLen = length $seqList[0];
				if ($tmpLen < $cutLen){
					print LOG "[Failed]\t$tmpFile contains shorter seqs (< $cutLen)\n"; next;
				}
			}
			if($cutSim > 0){
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