#collect for virusfinder

#perl /sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/bin/collect_data.pl "/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/11/normal/virusfinder1/" "/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/11/normal/virusfinder/" "11N" "/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/11/normal/file_name";

$workdir=$ARGV[0];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/11/normal/virusfinder1/";
$data_dir=$ARGV[1];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/11/normal/virusfinder/";

$AS=$ARGV[2];#"11N";

$file_name=$ARGV[3];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/11/normal/file_name";



$data_dir1=$workdir."output/step1";
$data_dir2=$workdir."output/step2";
$daa_dir3=$workdir."output/step3/b-align";



open(logfile,$file_name) or die"I can't read the file!";
$ns=0;
for $line(<logfile>){
	@a=split(/\n/,$line);
	$sp[$ns]=$a[0];
	$ns=$ns+1;
}
close(logfile);
print "number of sample:$ns\n";

###################collct unmapped reads#############################

$bam_all_combine=$workdir."align_hsa/bam_combine_all.lsf";
$bam_combine_sub=$workdir."align_hsa/bam_combine_sub.sh";
open(file1,">",$bam_all_combine) or die"I can't write to the file!";
open(file2,">",$bam_combine_sub) or die"I can't write to the file!";
print file1 "#!/bin/bash\n";
print file1 "#BSUB -P acc_GTEX\n";
print file1 "#BSUB -q premium\n";
print file1 "#BSUB -J VF_combine_bam\n";
print file1 "#BSUB -R \"rusage[mem=10000]\"\n";
print file1 "#BSUB -W 8:00\n";
print file1 "#BSUB -m mothra\n";
print file1 "#BSUB -o $workdir","align_hsa/$AS","_combine.stdout\n";
print file1 "#BSUB -eo $workdir","align_hsa/$AS","_combine.stderr\n";
print file1 "#BSUB -L /bin/bash\n";

print file1 "module load samtools\n";
for $i(0..$ns-1){
	$spf=$data_dir."unmapped/$sp[$i]";
	`ls $spf/*_unmapped.stderr > temp`;
	open(logfile,"temp") or die"I can't read the file!";
	$nsf=0;
	for $line(<logfile>){
		@a=split(/\//,$line);
		@b=split(/\_/,$a[$#a]);
		if(length($b[0])==2){
			$ffs[$nsf]=$b[0];
			$nsf=$nsf+1;
		}
	}
	print "$nsf\n";	
	close(logfile);
	`mkdir $workdir\"align_hsa/\"$sp[$i]`;
	
	$bam_cmb_file=$workdir."align_hsa/$sp[$i]/$sp[$i]_combine_bam.lsf";
	print $bam_cmb_file,"\n";
	open(file,">",$bam_cmb_file) or die"I can't write to the file!";
        print file "#!/bin/bash\n";
        print file "#BSUB -P acc_GTEX\n";
        print file "#BSUB -q premium\n";
        print file "#BSUB -J VF_combine_bam_sub\n";
        print file "#BSUB -R \"rusage[mem=10000]\"\n";
        print file "#BSUB -W 2:00\n";
        print file "#BSUB -m mothra\n";
        print file "#BSUB -o $workdir","align_hsa/$sp[$i]/$sp[$i]","_unmapped.stdout\n";
        print file "#BSUB -eo $workdir","align_hsa/$sp[$i]/$sp[$i]","_unmapped.stderr\n";
        print file "#BSUB -L /bin/bash\n";
	
	for $j(0..$nsf-1){
		$unmapped1=$spf."/$ffs[$j]/unmapped.1.fq";
		$unmapped2=$spf."/$ffs[$j]/unmapped.2.fq";
		$unmapped11=$spf."/$ffs[$j]/unmapped.1.fa";
		$unmapped21=$spf."/$ffs[$j]/unmapped.2.fa";
		if(-e $unmapped1){
			`cat $unmapped1 >> $data_dir1/unmapped.1.fq`;	
			`cat $unmapped11 >> $data_dir1/unmapped.1.fa`
		} 	
		else{
			print "error! $unmapped1 not available!\n";
		}
		if(-e $unmapped2){
			`cat $unmapped2 >> $data_dir1/unmapped.2.fq;`;
			`cat $unmapped21 >> $data_dir1/unmapped.2.fa;`
		}
		else{
			print "error! $unmapped2 not available!\n";	
		}	
		$bamfile=$spf."/$ffs[$j]/alignment.bam";
		$bamfile1=$workdir."align_hsa/$sp[$i]/alignment_$sp[$i]_$ffs[$j].bam";
		`ln -s $bamfile $bamfile1`;
	}
	print file "cd $workdir","align_hsa/$sp[$i]\n";	
	print file "module load samtools\n";
	print file "samtools merge alignment_$sp[$i].bam *.bam\n";
	print file1 "ln -s $workdir","align_hsa/$sp[$i]/alignment_$sp[$i].bam $workdir","align_hsa/alignment_$sp[$i].bam\n";
	close(file);
	print file2 "bsub < $bam_cmb_file\n";	
}
print file1 "cd $workdir","align_hsa\n";
print file1 "samtools merge alignment_$AS.bam *.bam\n";
print file1 "ln -s alignment_$AS.bam $workdir","output/step1/alignment.bam\n";
close(file1);
close(file2);

$cmd = "ln -s $data_dir"."integ_step3/alignment_$AS","_merge.sorted.bam $data_dir3/alignment.sorted.bam";
`$cmd`;
$cmd = "ln -s $data_dir"."integ_step3/alignment_$AS","_merge.sorted.bam.bai $data_dir3/alignment.sorted.bam.bai";
`$cmd`;
$cmd = "echo \"fake\" > $data_dir1/B.unmapped.bam";
`$cmd`;
$cmd = "echo \"fake\" > $data_dir1/L.unmapped.bam";
`$cmd`;
$cmd = "echo \"fake\" > $data_dir1/R.unmapped.bam";
`$cmd`;

#`ln -s $data_dir/output/alignment_$AS","_merge.sorted.bam $data_dir3/alignment.sorted.bam`;
#`ln -s $data_dir/output/alignment_$AS","_merge..sorted.bam.bai $data_dir3/alignment.sorted.bam.bai`;


