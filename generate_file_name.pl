#perl "/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/bin/generate_file_name.pl" "/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/22/normal/"
$dir=$ARGV[0];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/22/normal/";

$filename=$dir."file_name";

`ls $dir > temp`;

open(logfile,"temp") or die"I can't read the file!";
$nf=0;
for $line(<logfile>){
	@a=split(/\_/,$line);
	if($#a>0){

		if($nf==0){
			$fi[0]=$a[0];
			$nf=$nf+1;

		}
		else{
			$l=0;
			for $i(0..$nf-1){
				if($a[0] eq $fi[$i]){
					$l=1;	
				}
			}
			if($l==0){
				$fi[$nf]=$a[0];
				print $a[0],"\n";
				$nf=$nf+1;
			}
		}
	}	
}
print "number of files:$nf\n";
close(logfile);

open(file,">",$filename) or die"I can't write to the file!";
for $i(0..$nf-1){
	if($fi[$i] ne "file"){
		print file $fi[$i],"\n";
	}	
}
close(file);
