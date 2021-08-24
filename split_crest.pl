#!/usr/bin/perl
use strict;
use warnings;


my $dir="/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/";
my $sample="22"; #sample
my $phe="tumor"; #phenotype
my $step="b-align"; # could be d-align or b-align
my $nsub=1; # number of posssible region submit to crest, could be 1

my $que="premium";
my $mem="10000";
my $time="1:00";
my $flank_region_size=500;
my $cap3_bin="/sc/orga/projects/zhuj05a/Wenhui/HBV/VirusFinder/VirusFinder2.0/bin/cap3";
my $blat_bin="/sc/orga/projects/zhuj05a/Wenhui/HBV/VirusFinder/VirusFinder2.0/bin/blat";
my $read_len="90";  # read length


my $dir_cur=$dir.$sample."/".$phe."/virusfinder1/output/step3/$step/";

my $svdetect_file=$dir_cur."SVDetect/results-virus-loci.txt";

my $workdir=$dir_cur."crest/split/";

my $bin_dir="/sc/orga/projects/zhuj05a/Wenhui/HBV/VirusFinder/VirusFinder2.0/bin/";
if(!-d $workdir){
	`mkdir $workdir`;
}

my $all_sh=$dir_cur."crest/split/split_lsf.sh";

#collect all the non repeat region from svdetect
print $svdetect_file,"\n";
open(my $FILEX,"<",$svdetect_file) or die"I can't read the file!";
my $nreg=0;
my @reg;
my @tmp;
while(<$FILEX>){
	chomp;
	my @a=split(/\t/,$_);
	my $l=0;
	for my $i(0..$nreg-1){
		if($a[3] eq $tmp[$i]){
			$l=1;
			last;	
		}	
	}	
	if($l==0){
		$tmp[$nreg]=$a[3];	
		$reg[$nreg]=$_;
		$nreg++;
	}
}
close($FILEX);
print "number of total reagions:$nreg\n";
###########################################################


##############################generate the input region file for each split############################
my $u=int($nreg/$nsub);
for my$i(0..$u-1){
	my$workdir_sub=$workdir."split$i/";
	if(!-d $workdir_sub){
		`mkdir $workdir_sub`;
	}
	my $outsubfile=$workdir_sub."results-virus-loci.txt";
	open(my $FILE,">",$outsubfile) or die "I can't write to the file!";
	for my $j(0..$nsub-1){
		print $FILE $reg[$i*$nsub+$j],"\n";	
	}
	close($FILE);
}
if($u*$nsub < $nreg){
	my $workdir_sub=$workdir."split".$u."/";
	if(!-d $workdir_sub){
		`mkdir $workdir_sub`;
	}
	my $outsubfile=$workdir_sub."results-virus-loci.txt";
	open(my $FILE,">",$outsubfile) or die"I can't write to the file!";
	for my $i($nsub*$u..$nreg-1){
		print $FILE $reg[$i],"\n";
	}
	$u=$u+1;
}
#############################################################################
#my @tarchr;
#my $v=0;
#for $i(0..$nreg-1){
#	my @a=split(/\t/,$reg[$i]);
#	$tarchr[$v]=$a[3];
#	$v++;	
#}

#for $i(0..$u-1){
#	my $outsclipt=$workdir_sub."alignment.sorted.bam.sclip.txt";
#	my $outcover=$workdir_sub."alignment.sorted.bam.cover";
	
#	$cmd1="cat $dir_cur/crest/alignment.sorted.bam.$tarchr[$i].*.slip.txt > $outsclipt";
#	$cmd2="cat $dir_cur/crest/alignment.sorted.bam.$tarchr[$i].*.cover > $outcover";
#	`$cmd1`;
#	`$cmd2`;
#	$cmd3="cat $dir_cur/crest/alignment.sorted.bam.$tarchr[$i].chrVirus.slip.txt > $outsclipt";	
#	$cmd4="cat $dir_cur/crest/alignment.sorted.bam.$tarchr[$i].chrVirus.cover > $outcover";
#	`$cmd3`;
#	`$cmd4`;
#}

############################################################################get all the sub split lsf file##############################
open(my $FILE, ">",$all_sh) or die"I can't write to the file!";
for my $i(0..$u-1){
	my $lsf_file=$workdir.$step."_crest_sub_$i.lsf";
	my $stdout= $workdir."split_".$i.".stdout";
	if(-e $stdout){
		`rm $stdout`;
	}
	my $stderr= $workdir."split_".$i.".stderr";
        if(-e $stderr){
                `rm $stderr`;
        }
	open(my $FILEX,">",$lsf_file) or die"I can't write to the file!";
	print $FILEX "#!/bin/bash\n";
        print $FILEX "#BSUB -P acc_GTEX\n";
        print $FILEX "#BSUB -q $que\n";
        print $FILEX "#BSUB -J $sample","_$phe","_$i","_crest\n";
        print $FILEX "#BSUB -R \"rusage[mem=$mem]\"\n";
        print $FILEX "#BSUB -W $time\n";
        print $FILEX "#BSUB -m mothra\n";
        print $FILEX "#BSUB -o $workdir","split_$i.stdout\n";
        print $FILEX "#BSUB -eo $workdir","split_$i.stderr\n";
        print $FILEX "#BSUB -L /bin/bash\n";
	print $FILEX "export PERL5LIB=/hpc/users/wangm08/packages/trinity/trinityrnaseq_r2012-06-08/PerlLib:/hpc/users/wangm08/packages/BIO_DB_SAM/lib/perl5:/hpc/packages/minerva-common/CPAN/5.10.1/lib/perl5\n";
        print $FILEX "export JAVA_HOME=/hpc/packages/minerva-common/java/1.6.0_30/jdk1.6.0_30\n";
        print $FILEX "export PATH=\"/hpc/packages/minerva-common/java/1.6.0_30/jdk1.6.0_30/bin:\$PATH\"\n";
        print $FILEX "export _JAVA_OPTIONS=-Xmx1G\n";
        print $FILEX "export SAMTOOLS=/hpc/users/wangm08/packages/samtools/samtools-0.1.18\n";
        print $FILEX "export PATH=\"/hpc/users/wangm08/packages/samtools/samtools-0.1.18:\$PATH\"\n";
        print $FILEX "cd /sc/orga/projects/zhuj05a/Wenhui/HBV/VirusFinder/VirusFinder2.0\n";
        print $FILEX "./sys_check.pl\n";	
	print $FILEX "cd $bin_dir\n";
	print $FILEX "perl CREST.pl -f $dir_cur","crest/alignment.sorted.bam.cover -d $dir_cur","alignment.sorted.bam --ref_genome $dir_cur","hg19+virus.fa --boundary $workdir","split$i/results-virus-loci.txt --flankregion $flank_region_size -t $dir_cur","crest/hg19+virus.2bit -s $dir_cur","crest/virus.2bit --cap3 $cap3_bin --blat $blat_bin -l $read_len -o $workdir","split$i/ -p results --rmdup  --sensitive --min_sclip_reads 1\n";

	close($FILEX);
	print $FILE "bsub < $lsf_file\n";
}

close($FILE);
###############################################################################################################################################
