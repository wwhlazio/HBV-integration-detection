

$sample=$ARGV[0];#"11";
$phe=$ARGV[1];#"tumor";


$datadir=$ARGV[2];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/";

$workdir=$datadir.$sample."/".$phe."/virusfinder1/output/step3/d-align/crest/";

$reads_result=$workdir."reads_dis_on_breakpoint_hsa";
$sclipt_result=$workdir."extra_sclip_integraion";

print $reads_result,"\n";

open(logfile1,$reads_result) or die"I can't read the file!";
if(-e $sclipt_result){
	open(logfile2,$sclipt_result) or die"I can't read the file";
}
$output=$workdir."remap_sclipt_result";
#print $output,"\n";
open(file,">",$output) or die"I can't write to the file!";


$npos=0;
for $line(<logfile1>){
	@a=split(/\t/,$line);
	$pos[3*$npos]=$a[0];
	$pos[3*$npos+1]=$a[1];
	$pos[3*$npos+2]=$a[2];
	$npos++;
} 

#print "number of all input postion is :$npos\n";
close(logfile1);

if(-e $sclipt_result){
$nsp=0;
for $line(<logfile2>){
	@a=split(/\t/,$line);	
	$sp[2*$nsp]=$a[0];
	$sp[2*$nsp+1]=$a[1];
	$nsp++;
}
close(logfile2);


$rnsp=0;
$rsp[0]=$sp[0];
$rsp[1]=$sp[1]-1;
$rsp[2]=1;
for $i(1..$nsp-1){
	$l=0;
	for $j(0..$rnsp){
		if($rsp[3*$j+1]==($sp[2*$i+1]-1)){
			$rsp[3*$j+2]++;
			$l=1;
		}	
	}
	if($l==0){
		$rnsp++;
		$rsp[3*$rnsp]=$sp[2*$i];
		$rsp[3*$rnsp+1]=$sp[2*$i+1]-1;
		$rsp[3*$rnsp+2]=1;
	}
	
}
#print "rnsp:$rnsp\n";

#for $i(0..$rnsp){
#	print $rsp[3*$i],"\t",$rsp[3*$i+1],"\t",$rsp[3*$i+2],"\n";
#}

#print "\n\n";

#$output=$workdir."remap_sclipt_result";
#print $output,"\n";
#open(file,">",$output) or die"I can't write to the file!";
for $i(0..$npos-1){
	$u=0;
	for $j(0..$rnsp){
		if($rsp[3*$j] eq $pos[3*$i] && $rsp[3*$j+1] == $pos[3*$i+1]){
			$u=$rsp[3*$j+2];		
		}
	}	
#	print $pos[3*$i],"\t",$pos[3*$i+1],"\t\|",$pos[3*$i+2],"\|\t",$u,"\n";
	if($pos[3*$i+2]>0 || $u>0){
		@a=split(/\_/,$pos[3*$i]);
#		print file $pos[3*$i],"\t",$pos[3*$i+1],"\t",$pos[3*$i+2],"\t",$u,"\n";	
		print file $a[0],"\t",$a[1]+$pos[3*$i+1]-1,"\t",$pos[3*$i+2],"\t",$u,"\n";	
	}
}
}
else{
	for $i(0..$npos-1){
		if($pos[3*$i+2]>0){
			@a=split(/\_/,$pos[3*$i]);
#			print file $pos[3*$i],"\t",$pos[3*$i+1],"\t",$pos[3*$i+2],"\t0\n";
			print file $a[0],"\t",$a[1]+$pos[3*$i+1]-1,"\t",$pos[3*$i+2],"\t0\n";
		}
	}
}
close(file);
