# cd /sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/11/tumor/virusfinder1/output/step3/b-align/crest/
#!/usr/bin/perl
use warnings;
use strict;

my $sample=$ARGV[0];#"11";
my $phe=$ARGV[1];#"tumor";


my $datadir=$ARGV[2];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/";

my $workdir=$datadir.$sample."/".$phe."/virusfinder1/output/step3/d-align/";

my $er=3;

#`module load samtools`;

my $pos_file=$workdir."crest/reads_dis_on_breakpoint_hsa";
print $pos_file,"\n";
open(my $FILE1,$pos_file) or die"I can't read the file";
my $npos=0;
my @pos;
for my $line(<$FILE1>){
	chomp $line;
	my @a=split(/\t/,$line);
	if ($a[2] ne $a[6]){
		$pos[2*$npos]=$a[0];
		$pos[2*$npos+1]=$a[1];
		$npos=$npos+1;		
	}
	
}
close($FILE1);
print "number of pos:$npos\n";


`rm temp201`;
`rm temp202`;

my $output=$workdir."crest/reads_dis_on_breakpoint_hsa";
for my $i(0..$npos-1){

	my @a=split(/\_/,$pos[2*$i]);
        my $u=$a[2]-$a[1]+1;
	my $v;
	my $w;
        if($pos[2*$i+1]-1000<1){
                $v=1;
        }
        else{
                $v=$pos[2*$i+1]-1000;
        }
        if($pos[2*$i+1]+1000>$u){
                $w=$u;
        }else{
                $w=$pos[2*$i+1]+1000;
        }
        my $cmd1="samtools view $workdir"."alignment.sorted.bam ".$pos[2*$i].":".$v."-".$w."> temp201";
        `$cmd1`;
        print $cmd1,"\n";


	my $cmd111="awk '\$2==".($pos[2*$i+1]+1)."' alignment.sorted.bam.sclip.txt | cut -f4 |uniq";	
	print $cmd111,"\n";
	my $cmd11=`$cmd111`;
	my @aa=split(/\n/,$cmd11);
#	print $cmd11;

	if($#aa>=0){

		my $na=$#aa+1;
#		print "na:$na\n";

		my $nr=int($na/$er);
		if($na%$er>0){
			$nr=$nr+1;
		}
#		print "nr:$nr\n";
	
		my @bb;	
		if($nr>1){
			for my $j(0..$nr-2){
				$bb[$j]=$aa[$j*$er];
				for my $k(1..$er-1){
					$bb[$j]=$bb[$j]."|".$aa[$j*$er+$k];
				}
			}		
	
		}
	
		$bb[$nr-1]=$aa[($nr-1)*$er];
		for my $k(1..($na-$er*($nr-1)-1)){
			$bb[$nr-1]=$bb[$nr-1]."|".$aa[($nr-1)*$er+$k];	
		}
	
		for my $j(0..$nr-1){

			my $cmd21="grep -E \"$bb[$j]\" temp201 | grep -v 'SA:Z:chrVirus' | awk '\$6 !~ /^[1-9][0-9]M\$|\\*/'| cut -f1 | sort | uniq";
			print $cmd21,"\n";
			my $cmd211=`$cmd21`;
			print $cmd211,"\n";
			my @ab=split(/\n/,$cmd211);
			if($#ab>=0){
				my $nb=$#ab+1;
#				print $ab[0],"\n";
		
		
				my $nbr=int($nb/$er);
				if($nb%$er>0){
					$nbr=$nbr+1;	
				}
#				print "nbr:$nbr\n";	
		
				my @bc;
				if($nbr>1){
					for my $k(0..$nbr-2){
						$bc[$k]=$ab[$k*$er];
						for my $l(1..$er-1){	
							$bc[$k]=$bc[$k]."|".$ab[$k*$er+$l];
						}
					}
				}
				
		
				$bc[$nbr-1]=$ab[($nbr-1)*$er];
#				print "ab-nbr:",$ab[$nbr-1],"\t",$bc[$nbr-1],"\n";
				for my $k(1..($nb-$er*($nbr-1)-1)){
#					print $bc[$nbr-1],"\t";
					$bc[$nbr-1]=$bc[$nbr-1]."|".$ab[($nbr-1)*$er+$k];		
#					print "$bc[$nbr-1]\n";
				}
		
				for my $k(0..$nbr-1){
#					print $bc[$k],"\n";
					my $cmd22="grep -E \"$bc[$k]\" alignment.sorted.bam.sclip.txt | cut -f4,5 >> temp202";
					print $cmd22,"\n";
					`$cmd22`;
				}
			}
		}
	}
}


my $cmd3="awk '{OFS=\"\\t\"; print \">\"\$1\"\\n\"\$2}' temp202  > temp.fa";
print $cmd3,"\n";
`$cmd3`;
#$cmd4="blat -minScore=25 -minIdentity=85 virus.2bit temp.fa temp.ps1";
my $cmd4="blat -tileSize=7 -stepSize=1 -out=psl -minScore=15 -noHead -maxIntron=1 virus.2bit temp.fa temp.ps1";
print $cmd4,"\n";
`$cmd4`;

my $numr=`wc -l temp.ps1`;
my @numrr=split(/\ +/,$numr);
print "numr:$numrr[0]\n";
if($numrr[0]>0){
#	$cmd5="grep -E \$(awk 'NR>5' temp.ps1 | cut -f10 | sort | uniq | paste -sd \"\|\" ) alignment.sorted.bam.sclip.txt | sort | uniq > extra_sclip_integraion";
	my $cmd5="grep -E \$(cut -f10 temp.ps1 | sort | uniq | paste -sd \"\|\" ) alignment.sorted.bam.sclip.txt | sort | uniq > extra_sclip_integraion";
	print $cmd5,"\n";
	`$cmd5`;
}

