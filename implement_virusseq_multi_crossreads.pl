#split the file and filter the reads including duplicates and low quality reads



#perl /sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/bin/implement_virusseq_multi_crossreads.pl  1000 4000 20000 0:30 8:00 22N alloc /sc/orga/projects/zhuj05a/Wenhui/HBV/VirusSeq/Mosaik_bin /sc/orga/projects/zhuj05a/Wenhui/HBV/VirusSeq/Mosaik_JumpDb /sc/orga/projects/zhuj05a/Wenhui/HBV/VirusSeq/VirusSeq_Script /sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/22/normal/virus_seq/ /sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/22/normal/ 5000000;

$nspl=$ARGV[12];#5000000;
$mem1=$ARGV[0];#1000;
$mem2=$ARGV[1];#4000;
$mem3=$ARGV[2];#20000;
$time1=$ARGV[3];#"0:30";
$time2=$ARGV[4];#"8:00";
$AS=$ARGV[5];#"22N";
$que=$ARGV[6];#"alloc";



$bin_dir=$ARGV[7];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/VirusSeq/Mosaik_bin";
$ref_dir=$ARGV[8];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/VirusSeq/Mosaik_JumpDb";
$script_dir=$ARGV[9];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/VirusSeq/VirusSeq_Script";

$workdir=$ARGV[10];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/22/normal/virus_seq/";
$datadir=$ARGV[11];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/22/normal/";
$pipe_file1=$workdir."pipe_split.sh";
$pipe_filter=$workdir."pipe_filter.sh";


open(file1,">",$pipe_file1) or die "I can't write to the file";
open(file2,">",$pipe_filter) or die "I can't write to the file";

 

$sample_file=$datadir."\/"."file_name";
print $sample_file,"\n";



$ns=0;
open(logfile,$sample_file) or die"I can't read the file!";
for $line(<logfile>){
	@a=split(/\n/,$line);
	$sample[$ns]=$a[0];
	$ns=$ns+1;
}
print "number of file:$ns\n";


system("mkdir $workdir\"split\"");
system("mkdir $workdir\"filter\"");


for $j(0..$ns-1){

#############################################################split the input fastq######################################################################
	$sample_split=$workdir."split/$sample[$j]_split.lsf";
	open(file,">",$sample_split) or die"I can't write to the file!"; 
	print file "#!/bin/bash\n";
	print file "#BSUB -P acc_GTEX\n";
	print file "#BSUB -q $que\n";
	print file "#BSUB -J VS_split\n";
	print file "#BSUB -R \"rusage[mem=$mem2]\"\n";
	print file "#BSUB -W $time1\n";
	print file "#BSUB -m mothra\n";
	print file "#BSUB -o $workdir","split/$sample[$j]_split.stdout\n";
	print file "#BSUB -eo $workdir","split/$sample[$j]_split.stderr\n";
	print file "#BSUB -L /bin/bash\n";
                                                                                       
	print file "cd $datadir\n";
	print file "mkdir $workdir","split/$sample[$j]\n";
	print file "cp -r $sample[$j]_*.fastq.gz $workdir","split/$sample[$j]\n";
	print file "cd $workdir","split/$sample[$j]\n";
	print file "gunzip $sample[$j]_1.fastq.gz\n";
	print file "gunzip $sample[$j]_2.fastq.gz\n";
	print file "split -l $nspl $sample[$j]_1.fastq $sample[$j]_1_split_\n";
	print file "split -l $nspl $sample[$j]_2.fastq $sample[$j]_2_split_\n";
	print file "rm $sample[$j]_1.fastq\n";
	print file "rm $sample[$j]_2.fastq\n";
	close(file);
	print file1 "bsub < $sample_split\n";
 

#####################################################filter###########################################################################
	$sample_filter=$workdir."filter/$sample[$j]_filter.lsf";
	$cmd="mkdir $workdir"."filter/$sample[$j]";
	system("$cmd");
	open(file,">",$sample_filter) or die"I can't write to the file!";
	print file "#!/bin/bash\n";
	print file "#BSUB -P acc_GTEX\n";
	print file "#BSUB -q $que\n";
	print file "#BSUB -J VS_filter\n";
	print file "#BSUB -R \"rusage[mem=$mem1]\"\n";
	print file "#BSUB -W $time1\n";
	print file "#BSUB -m mothra\n";
	print file "#BSUB -o $workdir","filter/$sample[$j]_filter.stdout\n";
	print file "#BSUB -eo $workdir","filter/$sample[$j]_filter.stderr\n";
	print file "#BSUB -L /bin/bash\n";

	$samp_dirr=$workdir."split/".$sample[$j]."/";
	system("pwd");
	$cmd="ls ".$samp_dirr." > file_tmp";
	print $cmd,"\n";
	system("$cmd");
	open(logfile,"file_tmp") or die"I can't read the file\n";
	$nsf=0;
	for $line(<logfile>){
		@a=split(/\n/,$line);
		@b=split(/\_/,$a[0]);
		if($#b>0){
			if(length($b[$#b])==2){
				$ffs[$nsf]=$b[$#b];
				$nsf=$nsf+1;
			}
		}
	}
	close(logfile);
	print "number of splits for $sample[$j]:$nsf\n";
	
	for $i(0..$nsf-1){
		$sample_filter_sub=$workdir."filter/$sample[$j]/$sample[$j]_split_$ffs[$i]_filter.lsf";
		open(file1x,">",$sample_filter_sub) or die"I can't write to the file!";
		print file1x "#!/bin/bash\n";
		print file1x "#BSUB -P acc_GTEX\n";
		print file1x "#BSUB -q $que\n";
		print file1x "#BSUB -J VS_filter_sub\n";
		print file1x "#BSUB -R \"rusage[mem=$mem1]\"\n";
		print file1x "#BSUB -W $time1\n";
		print file1x "#BSUB -m mothra\n";
		print file1x "#BSUB -o $workdir","filter/$sample[$j]/$ffs[$i]_filter.stdout\n";
		print file1x "#BSUB -eo $workdir","filter/$sample[$j]/$ffs[$i]_filter.stderr\n";
		print file1x "#BSUB -L /bin/bash\n";
		
		print file1x "mkdir $workdir"."split/$sample[$j]/$ffs[$i]\n";
		print file1x "mv $samp_dirr","$sample[$j]_1_split_$ffs[$i] $workdir"."split/$sample[$j]/$ffs[$i]/$sample[$j]_1_split_$ffs[$i].fastq\n";
		print file1x "mv $samp_dirr","$sample[$j]_2_split_$ffs[$i] $workdir"."split/$sample[$j]/$ffs[$i]/$sample[$j]_2_split_$ffs[$i].fastq\n";
		print file1x "cd $workdir"."split/$sample[$j]/$ffs[$i]\n";
 		print file1x "perl /hpc/users/wangm08/packages/PRINSEQ/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $sample[$j]_1_split_$ffs[$i].fastq -fastq2 $sample[$j]_2_split_$ffs[$i].fastq -out_format 3 -log logfile";
                print file1x "  -derep 14 -min_qual_mean 20 -ns_max_p 10 -out_good file_passed -out_bad file_not_passed\n";
        	print file1x "rm file_not_passed*\n";
		print file1x "mv file_passed_1.fastq file_passed_1.fq\n";
	        print file1x "mv file_passed_2.fastq file_passed_2.fq\n";
		close(file1x);
		print file "bsub < $sample_filter_sub\n";
	}
	close(file);
	print file2 "bsub < $sample_filter\n";
}

close(file1);
close(file2);
