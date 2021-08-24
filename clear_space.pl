#perl /sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/bin/clear_space.pl "200" "normal" "N";
$sample=$ARGV[0];#"200";
$status=$ARGV[1];#"normal";
$s1=$ARGV[2];#"N";

$filename="/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/".$sample."/".$status."/file_name";
$dir="/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/".$sample."/".$status."/virus_seq/split/";
$dir1="/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/".$sample."/".$status."/virusfinder/unmapped/";
$dir2="/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/".$sample."/".$status."/virusfinder/integ_step2/";
$dir3="/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/".$sample."/".$status."/virusfinder/integ_step3/";
$dir4="/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/".$sample."/".$status."/virusfinder1/align_hsa/";

$ns=0;
open(logfile,$filename) or die"I can't write to the file!";
for $line(<logfile>){
	@a=split(/\n/,$line);
	$sp[$ns]=$a[0];
	$ns=$ns+1;
}
close(logfile);

$clear_file="clear_space_".$sample."_".$status.".sh";
open(file,">",$clear_file) or die"I can't write to the file!";
print file "rm $dir3","alignment_$sample","$s1","_merge.tmp.sorted.bam\n";
for $i(0..$ns-1){
	$dir_u=$dir.$sp[$i];
	$dir1_u=$dir1.$sp[$i];
	$dir2_u=$dir2.$sp[$i];
	$dir3_u=$dir3.$sp[$i];
	$dir4_u=$dir4.$sp[$i];
	print file "rm $dir2_u/alignment_$sp[$i]_merge.tmp.sorted.bam\n";
	@temp=`ls $dir_u`;
	for $j(0..$#temp-1){
		@a=split(/\n/,$temp[$j]);
		if(length($a[0])==2){
			print file "rm $dir_u/$a[0]/file_not_passed*\n";
#			print file "rm $dir_u/$a[0]/file_passed_*.fq.gz\n";
			print file "rm $dir_u/$a[0]/$sp[$i]_*_split_$a[0].fastq\n";
			print file "rm $dir1_u/$a[0]/B.unmapped.bam\n";
			print file "rm $dir1_u/$a[0]/L.unmapped.bam\n";
			print file "rm $dir1_u/$a[0]/R.unmapped.bam\n";
			print file "rm $dir1_u/$a[0]/alignment.sam\n";
	#		print file "rm $dir1_u/$a[0]/alignment.bam\n";
			print file "rm $dir1_u/$a[0]/unmapped.bam\n";
	#		print file "rm $dir1_u/$a[0]/unmapped.*.fq\n";
	#		print file "rm $dir1_u/$a[0]/unmapped.*.fa\n";
			print file "rm $dir1_u/$a[0]/integ_step1/unmapped.*.sai\n";
	#		print file "rm $dir1_u/$a[0]/integ_step1/alignment.sorted.bam*\n";
		}
	}
	print file "rm $dir2_u/alignment_$sp[$i]_merge.tmp.sorted.bam\n";		
	print file "rm $dir2_u/alignment_$sp[$i]_merge.sorted.bam\n";
	print file "rm $dir2_u/alignment_$sp[$i]_merge.sorted.bam.bai\n";
	print file "rm $dir3_u/alignment_$sample","$s1","_merge.tmp.sorted.bam\n";
	#print file "rm $dir4_u/alignment_$sp[$i].bam\n";
			
}
close(file);
