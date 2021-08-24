#perl /sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/bin/check_vf_pipeline_v1.pl "/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/22/normal/virusfinder/" "/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/22/normal/file_name" 

$workdir=$ARGV[0];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/22/normal/virusfinder/";
$file_name=$ARGV[1];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/22/normal/file_name";

open(logfile,$file_name) or die "I can't read the file";
$ns=0;
for $line(<logfile>){
	@a=split(/\n/,$line);
	$sp[$ns]=$a[0];
	$ns=$ns+1;
}
close(logfile);

####################################icheck unmapped#################

$stderr_dir=$workdir."unmapped/";
$redo_unmapped=$workdir."redo_unmapped.sh";
open(file,">",$redo_unmapped) or die"I can't write to the file!";
for $i(0..$ns-1){
	$stderr_dir1=$stderr_dir.$sp[$i];
	`ls $stderr_dir1/*_unmapped.stderr > temp`;
	open(logfile,"temp") or die"I can't write to the file!";
	for $line(<logfile>){
		@a=split(/\//,$line);
		@b=split(/\_/,$a[$#a]);
		if(length($b[0])==2){
			$fc=$stderr_dir1."/$b[0]/";
			$fc1=$fc."alignment.bam";	
			$fc2=$fc."unmapped.1.fq";
			$fc3=$fc."unmapped.2.fq";
			if( (!-e $fc1 || ! -s $fc1) || (!-e $fc2 || ! -s $fc2) || (!-e $fc2 || ! -s $fc2) ){
				print -e $fc1,"\t",-e $fc2,"\t",-e $fc3,"\n";		
				print -s $fc1,"\t",-s $fc2,"\t",-s $fc3,"\n";
				print file "bsub < $workdir","unmapped/$sp[$i]/$b[0]/$sp[$i]_$b[0]_unmapped.lsf\n";
				print "effor with unmapped for $sp[$i] and $b[0]\n";
			}		
		}
	}
	close(logfile);
}
close(file);

######################integ_step1#####################


$stderr_dir=$workdir."unmapped/";
$redo_unmapped=$workdir."redo_integ_step1.sh";
open(file,">",$redo_unmapped) or die"I can't write to the file!";
for $i(0..$ns-1){
        $stderr_dir1=$stderr_dir.$sp[$i];
        `ls $stderr_dir1/*.stderr > temp`;
        open(logfile,"temp") or die"I can't write to the file!";
        for $line(<logfile>){
		@a=split(/\//,$line);
                @b=split(/\_/,$a[$#a]);
		if(length($b[0])==2){
			$fc=$stderr_dir1."/$b[0]/integ_step1/";	
			$fc1=$fc."alignment.sorted.bam";
			$fc2=$fc."alignment.sorted.bam.bai";
#			print $fc1,"\n";
#			print $fc2,"\n";
			if((!-e $fc1 || !-s $fc1)||(!-e $fc2 || !-s $fc2)){
				$u1= -e $fc1;
				$u2= -e $fc2;
				$v1= -s $fc1;
                                $v2= -s $fc2;
				print $u1,"\t",$u2,"\n";
                                print $v1,"\t",$v2,"\n";
				print file "bsub < $workdir","integ_step1/$sp[$i]/$b[0]/$sp[$i]_$b[0]_integ_step1.lsf\n"; 
				print "effor with integ_step1 for $sp[$i] and $b[0]\n";
			}
		}
        } 
        close(logfile);
}
close(file);
################################################check integ_step2#####################

$stderr_dir=$workdir."integ_step2/";
$redo_unmapped=$workdir."redo_integ_step2.sh";
open(file,">",$redo_unmapped) or die"I can't write to the file!";
for $i(0..$ns-1){
	$fc1=$stderr_dir.$sp[$i]."/alignment.sorted.bam";
	$fc2=$stderr_dir.$sp[$i]."/alignment.sorted.bam.bai";	
	if((!-e $fc1 || !-s $fc1)||(!-e $fc2 || !-s $fc2)){
		print "effor with integ_step2 for $sp[$i]\n";
                print file "bsub < $workdir","integ_step2/$sp[$i]/integ_step2_$sp[$i].lsf\n";
	}      
}
close(file);

print "BE careful on the reported errors!\n";
