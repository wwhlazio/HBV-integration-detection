#
$sample=$ARGV[0];#"11";
$phe=$ARGV[1];#"tumor";


$datadir=$ARGV[2];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/";

$workdir=$datadir.$sample."/".$phe."/virusfinder1/output/step3/d-align/";

`module load samtools`;

$er=3;

$pos_file=$workdir."crest/results-virus-loci.txt-fp";
print $pos_file,"\n";
open(LOGFILE,$pos_file) or die"I can't read the file";
$npos=0;
$u=0;
for $line(<LOGFILE>){
	if($u>0){
	chomp $line;
	@a=split(/\t/,$line);
	if ($a[0] ne "chrVirus"){
		$pos[2*$npos]=$a[0];
		$pos[2*$npos+1]=$a[1];
		$npos=$npos+1;		
	}
	else{
		$pos[2*$npos]=$a[3];
		$pos[2*$npos+1]=$a[4];
		$npos=$npos+1;	
	}
	}
	$u=$u+1;
	
}
close(LOGFILE);
print "number of pos:$npos\n";


$output=$workdir."crest/reads_dis_on_breakpoint_hsa";
#$cmd="cat $workdir"."crest/alignment.sorted.bam.chr*sclip.txt > $workdir"."crest/alignment.sorted.bam.chr.temp.sclip.txt";
#`$cmd`;
open(file,">",$output) or die"I can't write to the file";
for $i(0..$npos-1){
	@a=split(/\_/,$pos[2*$i]);
	$u=$a[2]-$a[1]+1;
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
	$nhbv=0;
	$nschbv=0;
	$nasc=0;	
	$nssc=0;
	$nnsc=0;

	$cmd1="samtools view $workdir"."alignment.sorted.bam ".$pos[2*$i].":".$v."-".$w."> temp201";
	`$cmd1`;
	print $cmd1,"\n";
	$cmd111="awk '\$2==".($pos[2*$i+1]+1)."' alignment.sorted.bam.sclip.txt | cut -f4 | uniq";
	print $cmd111,"\n";
	$cmd11=`$cmd111`;
	print $cmd11,"\n";
	@aa=split(/\n/,$cmd11);
	
	print $#aa,"\n";
	
	$na=$#aa+1;	

	if($na>0){	
	$nr=int($na/$er);
	if($#na%$er>0){
		$nr=$nr+1;
	}		

	print "nr:",$nr,"\n";

	if($nr>1){
	for $j(0..$nr-2){
		$bb[$j]=$aa[$j*$er];
		for $k(1..$er-1){
			$bb[$j]=$bb[$j]."|".$aa[$j*$er+$k];
		}
	}
	}	
	
	$bb[$nr-1]=$aa[($nr-1)*$er];
	for $k(1..($#aa-$er*($nr-1))){
		$bb[$nr-1]=$bb[$nr-1]."|".$aa[($nr-1)*$er+$k];
	}

	
	for $j(0..$nr-1){
		print $bb[$j],"\n";
		$cmd2="grep -E \"$bb[$j]\" temp201 | awk '\$6 !~ /^[1-9][0-9]M\$|\*/'  | awk '\$7 ~ \/chrVirus\/' | grep -v 'SA:Z:chrVirus' | wc -l";
		print $cmd2,"\n";
		$nhbv+=`$cmd2`;
		$cmd3="grep -E \"$bb[$j]\" temp201 |  awk '\$6 !~ /^[1-9][0-9]M\$|\*/' | grep 'SA:Z:chrVirus' | wc -l";
		print $cmd3,"\n";
		$nschbv+=`$cmd3`;
		$cmd4="grep -E \"$bb[$j]\" temp201 |  awk '\$6 !~ /^[1-9][0-9]M\$|\*/' | cut -f1 | sort | uniq | wc -l";
		print $cmd4,"\n";
		$nasc+=`$cmd4`;
		$cmd5="grep -E \"$bb[$j]\" temp201 |  awk '\$6 !~ /^[1-9][0-9]M\$|\*/' | awk '\$7 ~ \/\=\/' | grep -v 'SA:Z:chrVirus' | cut -f1 | sort | uniq | wc -l";
		print $cmd5,"\n";
		$nssc+=`$cmd5`;
		$cmd6="grep -E \"$bb[$j]\" temp201 |  awk '\$6 !~ /^[1-9][0-9]M\$|\*/' |awk '\$7 !~ \/\=\/' | grep -v 'SA:Z:chrVirus' | wc -l";
		print $cmd6,"\n";
		$nnsc+=`$cmd6`
	}



	#$cmd2="grep -E \$(awk '\$2==".($pos[2*$i+1]+1)."' alignment.sorted.bam.sclip.txt | cut -f4 | paste -sd \"\|\") temp201 | awk '\$6 !~ /^[1-9][0-9]M\$|\*/'  | awk '\$7 ~ \/chrVirus\/' | grep -v 'SA:Z:chrVirus' | wc -l"; #paired_end on hbv
	#print $cmd2,"\n";
	#$nhbv=`$cmd2`;
	#print $cmd2,"\n";
	#$cmd3="grep -E \$(awk '\$2==".($pos[2*$i+1]+1)."' alignment.sorted.bam.sclip.txt | cut -f4 | paste -sd \"\|\") temp201 |  awk '\$6 !~ /^[1-9][0-9]M\$|\*/' | grep 'SA:Z:chrVirus' | wc -l"; #sclipped read on hbv
	#print $cmd3,"\n";
	#$nschbv=`$cmd3`;
	#print $cmd3,"\n";
	#$cmd4="grep -E \$(awk '\$2==".($pos[2*$i+1]+1)."' alignment.sorted.bam.sclip.txt | cut -f4 | paste -sd \"\|\") temp201 |  awk '\$6 !~ /^[1-9][0-9]M\$|\*/' | cut -f1 | sort | uniq | wc -l";
	#print $cmd4,"\n";
	#$nasc=`$cmd4`;
	#print $cmd4,"\n";
	#$cmd5="grep -E \$(awk '\$2==".($pos[2*$i+1]+1)."' alignment.sorted.bam.sclip.txt | cut -f4 | paste -sd \"\|\") temp201 |  awk '\$6 !~ /^[1-9][0-9]M\$|\*/' | awk '\$7 ~ \/\=\/' | grep -v 'SA:Z:chrVirus' | cut -f1 | sort | uniq | wc -l"; #2*# paired-end on same chr
	#print $cmd5,"\n";
	#$nssc=`$cmd5`;
	#print $cmd5,"\n";
	#print "nssc:$nssc\n";
	#$cmd6="grep -E \$(awk '\$2==".($pos[2*$i+1]+1)."' alignment.sorted.bam.sclip.txt | cut -f4 | paste -sd \"\|\") temp201 |  awk '\$6 !~ /^[1-9][0-9]M\$|\*/' |awk '\$7 !~ \/\=\/' | grep -v 'SA:Z:chrVirus' | wc -l"; #other nosise
	#print $cmd6,"\n";
	#$nnsc=`$cmd6`;
	#print $cmd6, "\n";
	chomp $nschbv;
	chomp $nasc;
	chomp $nhbv;
	chomp $nssc;
	chomp $nnsc;
	print file $pos[2*$i]."\t",$pos[2*$i+1],"\t",$nschbv,"\t",$nhbv,"\t",$nssc,"\t",$nnsc-$nhbv,"\t",$nasc,"\n";
	}
}
close(file);
