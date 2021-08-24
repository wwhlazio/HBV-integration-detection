#check the results of split and filtering
#perl check_pipe_line_v1.pl "/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/180/tumor/virus_seq/" "/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/180/tumor/file_name"

$workdir=$ARGV[0];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/180/tumor/virus_seq/";
$file_name=$ARGV[1];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/180/tumor/file_name";

$ns=0;
open(logfile,$file_name) or die"I can't read the file!";
for $line(<logfile>){
	@a=split(/\n/,$line);
	$sp[$ns]=$a[0];
	$ns=$ns+1;	
}
print "number of sample:$ns\n";
close(logfile);

######################################check split###############################################

$redo_split=$workdir."redo_split.sh";
open(file,">",$redo_split) or die"I can't write to the file!";
$stderr_dir=$workdir."split/";
for $i(0..$ns-1){
	$stderr=$stderr_dir.$sp[$i]."_split.stderr";
	open(logfile,$stderr) or die"I can't read the file!";
	@ff=<logfile>;
	if($#ff>0){
		print "error of split with $sp[$i]\n";
		print file "bsub < $workdir","split/$sp[$i]_split.lsf\n"; 
	}
	close(logfile);
}
close(file);

#####################################3check filter#####################################
$redo_filter=$workdir."redo_filter.sh";
open(file,">",$redo_filter) or die"I can't write to the file!";
$stderr_dir=$workdir."filter/";
$stderr_dirx=$workdir."split/";
for $i(0..$ns-1){
	$stderr_dir1=$stderr_dir."$sp[$i]/";
	$stderr_dirx1=$stderr_dirx."$sp[$i]/";
	system("ls $stderr_dir1*_filter.stderr > temp");
	open(logfile,"temp") or die"I can't read the file!";
	for $line(<logfile>){
		@a=split(/\//,$line);
		@b=split(/\_/,$a[$#a]);
		if(length($b[0])==2){
			$fc=$stderr_dirx1.$b[0]."/";
			$fc1=$fc."file_passed_1.fq";
			$fc2=$fc."file_passed_2.fq";
			if((!-e $fc1 || !-s $fc1)||(!-e $fc2 || !-s $fc2)){
				print "error of filter with $sp[$i] and $b[0]\n";
				print file "bsub < $stderr_dir1","$sp[$i]_split_$b[0]","_filter.lsf\n";
			}
		}
	}
	close(logfile);
}
close(file);
