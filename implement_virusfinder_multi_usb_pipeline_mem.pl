#modify to use the bwa mem
#"the idea of split the reference chromosome will not work"
#use sub function the firt time"
#formalize the pipeline for reusefulness


#perl implement_virusfinder_multi_usb_pipeline_mem.pl "/sc/orga/projects/zhuj05a/Wenhui/HBV/VirusFinder/VirusFinder2.0/" "2:00" "8:00" "48:00" "1000" "8000" "25000" "alloc" "/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/11/tumor" "/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/11/tumor/virus_seq/split/" "/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/11/tumor/virusfinder/" "11T" "/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/82/tumor/configure_file_82T" "/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/82/tumor/hbv_ref_ncbi/HBV.fa" "/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/82/tumor/hbv_ref_ncbi/bwa_hg19_hbv_index" 


sub sub_pipe_line;
sub sub_pipe_line_new;
sub sub_pipe_line_integ_combine;


$vf_bin=$ARGV[0];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/VirusFinder/VirusFinder2.0/";

$time1=$ARGV[1];#"2:00";
$time2=$ARGV[2];#"8:00";
$time3=$ARGV[3];#"48:00";

$mem1=$ARGV[4];#"1000";
$mem2=$ARGV[5];#"8000";
$mem3=$ARGV[6];#"25000";

$que=$ARGV[7];#"alloc";

$data_dir=$ARGV[8];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/wgs_88/11/tumor";				#different for each data
$data_file=$data_dir."/file_name";								
$data_split=$ARGV[9];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/11/tumor/virus_seq/split/";	#different for each data
$workdir=$ARGV[10];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/11/tumor/virusfinder/";		#different for each data
$AS=$ARGV[11];#"11T";

$temp_config=$ARGV[12];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/82/tumor/configure_file_82T";  #unchanged
$hbv_ref=$ARGV[13];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/82/tumor/hbv_ref_ncbi/HBV.fa";	#unchanged
$dir_hg19_hbv_bwa_index=$ARGV[14];#"/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/82/tumor/hbv_ref_ncbi/bwa_hg19_hbv_index";	#unchanged


$pipe_unmapped=$workdir."extract_unmapped_reads";
$pipe_integ_step1=$workdir."call_integ_step1";
$pipe_integ_step2=$workdir."combine_integ_step2";
$pipe_integ_step3=$workdir."combine_integ_all.lsf";

#`cd $datadir`;
#`ls -l *.faste.gz > file_name`;

`mkdir $workdir\"unmapped\"`;
`mkdir $workdir\"integ_step1\"`;
`mkdir $workdir\"integ_step2\"`;
`mkdir $workdir\"integ_step3\"`;


open(file1,">",$pipe_unmapped) or die"I can't write to the file!";
open(file2,">",$pipe_integ_step1) or die"I can't write to the file!";
open(file3,">",$pipe_integ_step2) or die"I can't write to the file!";
open(file4,">",$pipe_integ_step3) or die"I can't write to the file!";


open(logfile,$temp_config) or die"I can't read the file!";
@template=<logfile>;


print $data_file,"\n";
open(logfile,$data_file) or die"I can't read the file!";
$ns=0;
for $line(<logfile>){
	@a=split(/\n/,$line);	
	$fp[$ns]=$a[0];
	$ns=$ns+1;	
}
print "number of files $ns\n";


#step1: splite the inout fastq files.
#use the split files from virus seq. The ones after filtering duplication.

#step2: align the reads to human and get the unmapped reads

print file4 "#!/bin/bash\n";
print file4 "#BSUB -P acc_GTEX\n";
print file4 "#BSUB -q $que\n";
print file4 "#BSUB -J VF_combine3\n";
print file4 "#BSUB -R \"rusage[mem=$mem2]\"\n";
print file4 "#BSUB -W $time2\n";
print file4 "#BSUB -m mothra\n";
print file4 "#BSUB -o $workdir","integ_step3/","integ_step3.stdout\n";
print file4 "#BSUB -eo $workdir","integ_step3/","integ_step3.stderr\n";
print file4 "#BSUB -L /bin/bash\n";
print file4 "export PERL5LIB=/hpc/users/wangm08/packages/trinity/trinityrnaseq_r2012-06-08/PerlLib:/hpc/users/wangm08/packages/BIO_DB_SAM/lib/perl5:/hpc/packages/minerva-common/CPAN/5.10.1/lib/perl5\n";
print file4 "export JAVA_HOME=/hpc/packages/minerva-common/java/1.6.0_30/jdk1.6.0_30\n";
print file4 "export PATH=\"/hpc/packages/minerva-common/java/1.6.0_30/jdk1.6.0_30/bin:\$PATH\"\n";
print file4 "export _JAVA_OPTIONS=-Xmx1G\n";
print file4 "export SAMTOOLS=/hpc/users/wangm08/packages/samtools/samtools-0.1.18\n";
print file4 "export PATH=\"/hpc/users/wangm08/packages/samtools/samtools-0.1.18:\$PATH\"\n";
print file4 "cd /sc/orga/projects/zhuj05a/Wenhui/HBV/VirusFinder/VirusFinder2.0\n";
print file4 "./sys_check.pl\n";
print file4 "cd $workdir","integ_step3\n";
print file4 "rm -rf $workdir","integ_step3/*.bam\n";
print file4 "rm -rf $workdir","integ_step3/*.bam.bai\n";

for $i(0..$ns-1){
	$cmd="mkdir $workdir"."unmapped/$fp[$i]/";
	system("$cmd");
	$cmd="rm $data_split".$fp[$i]."/temp";
	system("$cmd");
	$cmd="ls $data_split"."$fp[$i]/ > $data_split"."$fp[$i]/temp";
	system("$cmd");
	$file_split=$data_split.$fp[$i]."/temp";
	open(logfile,$file_split) or die"I can't read the file!";
	$nsf=0;
	for $line(<logfile>){
		@a=split(/\n/,$line);	
		$fps[$nsf]=$a[0];
		$nsf=$nsf+1;
	}
	$nsf=$nsf-1;	
	
	$align_for_unmap=$workdir."unmapped/$fp[$i]/align_for_unmap_$fp[$i].lsf";
	sub_pipe_line($align_for_unmap,$workdir,$fp[$i],$mem1,$mem2,$time1,$time2,$data_split,$nsf,@fps);
	print file1 "bsub < $align_for_unmap\n"; 

	$cmd="mkdir $workdir"."integ_step1/$fp[$i]/";
	system("$cmd");	
	
	$integ_step1=$workdir."integ_step1/$fp[$i]/integ_step1_$fp[$i].lsf";
	sub_pipe_line_new($integ_step1,$workdir,$fp[$i],$mem1,$mem2,$time1,$time2,$data_split,$nsf,"integ_step1","integ_step1",@fps);
	print file2 "bsub < $integ_step1\n";	

	$cmd="mkdir $workdir"."integ_step2/$fp[$i]";
	system("$cmd");
	
	$integ_step2=$workdir."integ_step2/$fp[$i]/integ_step2_$fp[$i].lsf";
	sub_pipe_line_integ_combine($integ_step2,$workdir,$fp[$i],$mem1,$time1,$nsf,"integ_step2",@fps);
	print file3 "bsub < $integ_step2\n";
	
	print file4 "ln -s $workdir","integ_step2/$fp[$i]/alignment.sorted.bam ","alignment_$fp[$i]_merge.sorted.bam\n"; 
	print file4 "ln -s $workdir","integ_step2/$fp[$i]/alignment.sorted.bam.bai ","alignment_$fp[$i]_merge.sorted.bam.bai\n";
}	
print file4 "module load samtools\n";
print file4 "samtools merge alignment_$AS","_merge.tmp.sorted.bam *.bam\n";
print file4 "samtools sort alignment_$AS","_merge.tmp.sorted.bam alignment_$AS","_merge.sorted\n";
print file4 "samtools index alignment_$AS","_merge.sorted.bam\n";

 



close(file1);
close(file2);
close(file3);
close(file4);



#sub_pip_line($workdir,$fp[$i],$mem1,$mem2,$time1,$time2,$data_split,$nsf,@fps);
sub sub_pipe_line{
	my $align_for_unmap= @_[0];
	my $workdir = @_[1];
	my $fp = @_[2];
	my $mem1 = @_[3];
	my $mem2 = @_[4];
	my $time1 = @_[5];
	my $time2 = @_[6];
	my $data_split = @_[7];
	my $nsf = @_[8];
	my @fps;

	for $i(9..scalar(@_)-1){
		$fps[$i-9]=@_[$i];
	}

        my $align_for_unmap=$workdir."unmapped/$fp/align_for_unmap_$fp.lsf";
        open(file,">",$align_for_unmap) or die"I can't write to the file!";
        print file "#!/bin/bash\n";
        print file "#BSUB -P acc_GTEX\n";
        print file "#BSUB -q $que\n";
        print file "#BSUB -J VF_unmap\n";
        print file "#BSUB -R \"rusage[mem=$mem1]\"\n";
        print file "#BSUB -W $time1\n";
        print file "#BSUB -m mothra\n";
        print file "#BSUB -o $workdir","unmapped/$fp","_unmapped.stdout\n";
        print file "#BSUB -eo $workdir","unmapped/$fp","_unmapped.stderr\n";
        print file "#BSUB -L /bin/bash\n";
	
        for my $j(0..$nsf-1){
                my $cmd="mkdir $workdir"."unmapped/$fp[$i]/$fps[$j]/"; 
                system("$cmd");

                my $cff=$workdir."unmapped/$fp/$fps[$j]/config_$fp_$fps[$j]";
                open(filex,">",$cff) or die "I can't write to the file";
                @template_tmp=@template;
                $template_tmp[0]="fastq1 = $data_split"."$fp/$fps[$j]/file_passed_1.fq\n";
                $template_tmp[1]="fastq2 = $data_split"."$fp/$fps[$j]/file_passed_2.fq\n";
                print filex @template_tmp;
                close(filex);

                my $align_for_unmap_sub=$workdir."unmapped/$fp/$fps[$j]/$fp"."_$fps[$j]_unmapped.lsf";
                open(filex,">",$align_for_unmap_sub) or die"I can't write to the file!";
                print filex "#!/bin/bash\n";
                print filex "#BSUB -P acc_GTEX\n";
                print filex "#BSUB -q $que\n";
                print filex "#BSUB -J VF_unmap_sub\n";
                print filex "#BSUB -R \"rusage[mem=$mem2]\"\n";
                print filex "#BSUB -W $time2\n";
                print filex "#BSUB -m mothra\n";
                print filex "#BSUB -o $workdir","unmapped/$fp/$fps[$j]","_unmapped.stdout\n";
                print filex "#BSUB -eo $workdir","unmapped/$fp/$fps[$j]","_unmapped.stderr\n";
                print filex "#BSUB -L /bin/bash\n";
		print filex "export PERL5LIB=/hpc/users/wangm08/packages/trinity/trinityrnaseq_r2012-06-08/PerlLib:/hpc/users/wangm08/packages/BIO_DB_SAM/lib/perl5:/hpc/packages/minerva-common/CPAN/5.10.1/lib/perl5\n";
		print filex "export JAVA_HOME=/hpc/packages/minerva-common/java/1.6.0_30/jdk1.6.0_30\n";
		print filex "export PATH=\"/hpc/packages/minerva-common/java/1.6.0_30/jdk1.6.0_30/bin:\$PATH\"\n";
		print filex "export _JAVA_OPTIONS=-Xmx1G\n";
		print filex "export SAMTOOLS=/hpc/users/wangm08/packages/samtools/samtools-0.1.18\n";
		print filex "export PATH=\"/hpc/users/wangm08/packages/samtools/samtools-0.1.18:\$PATH\"\n";
		print filex "cd /sc/orga/projects/zhuj05a/Wenhui/HBV/VirusFinder/VirusFinder2.0\n";
		print filex "./sys_check.pl\n";
                print filex "cd $data_split","$fp/$fps[$j]/\n";
		print filex "rm $workdir","unmapped/$fp/$fps[$j]/unmapped*\n";
                print filex "rm $workdir","unmapped/$fp/$fps[$j]/alignment*\n";
                print filex "rm $workdir","unmapped/$fp/$fps[$j]","_unmapped.stdout\n";
                print filex "rm $workdir","unmapped/$fp/$fps[$j]","_unmapped.stderr\n";
                print filex "perl $vf_bin","preprocess_new.pl -c $cff -o $workdir","unmapped/$fp/$fps[$j]/\n";
                close(filex);
                print file "bsub < $align_for_unmap_sub\n";
        }	
	close(file);
}


sub sub_pipe_line_new{
        my $align_for_unmap= @_[0];
        my $workdir = @_[1];
        my $fp = @_[2];
        my $mem1 = @_[3];
        my $mem2 = @_[4];
        my $time1 = @_[5];
        my $time2 = @_[6];
        my $data_split = @_[7];
        my $nsf = @_[8];
	my $label_lsf = @_[9];
	my $label_dir = @_[10];
        my @fps;

        for $i(11..scalar(@_)-1){
                $fps[$i-11]=@_[$i];
        }

	open(file,">",$align_for_unmap) or die"I can't write to the file!";
        print file "#!/bin/bash\n";
        print file "#BSUB -P acc_GTEX\n";
        print file "#BSUB -q $que\n";
        print file "#BSUB -J VF_call_int\n";
        print file "#BSUB -R \"rusage[mem=$mem1]\"\n";
        print file "#BSUB -W $time1\n";
        print file "#BSUB -m mothra\n";
        print file "#BSUB -o $workdir","$label_dir/$fp","_$label_dir.stdout\n";
        print file "#BSUB -eo $workdir","$label_dir/$fp","_$label_dir.stderr\n";
        print file "#BSUB -L /bin/bash\n";

        for my $j(0..$nsf-1){
                my $cmd="mkdir $workdir"."$label_dir/$fp[$i]/$fps[$j]/";
                system("$cmd");

                my $cff=$workdir."$label_dir/$fp/$fps[$j]/config_$fp_$fps[$j]";
                open(filex,">",$cff);
                @template_tmp=@template;
                $template_tmp[0]="fastq1 = $data_split"."$fp/$fps[$j]/file_passed_1.fq\n";
                $template_tmp[1]="fastq2 = $data_split"."$fp/$fps[$j]/file_passed_2.fq\n";
		$template_tmp[8]="bwa_bin = /hpc/packages/minerva-common/bwa/0.7.8/bwa-0.7.8/bwa\n";
                print filex @template_tmp;
                close(filex);

                my $align_for_unmap_sub=$workdir."$label_dir/$fp/$fps[$j]/$fp"."_$fps[$j]_$label_dir.lsf";
		open(filex,">",$align_for_unmap_sub) or die"I can't write to the file!";
                print filex "#!/bin/bash\n";
                print filex "#BSUB -P acc_GTEX\n";
                print filex "#BSUB -q $que\n";
                print filex "#BSUB -J VF_call_int_sub\n";
                print filex "#BSUB -R \"rusage[mem=$mem2]\"\n";
                print filex "#BSUB -W $time1\n";
                print filex "#BSUB -m mothra\n";
                print filex "#BSUB -o $workdir","$label_dir/$fp/$fps[$j]","_$label_dir.stdout\n";
                print filex "#BSUB -eo $workdir","$label_dir/$fp/$fps[$j]","_$label_dir.stderr\n";
                print filex "#BSUB -L /bin/bash\n";
                print filex "export PERL5LIB=/hpc/users/wangm08/packages/trinity/trinityrnaseq_r2012-06-08/PerlLib:/hpc/users/wangm08/packages/BIO_DB_SAM/lib/perl5:/hpc/packages/minerva-common/CPAN/5.10.1/lib/perl5\n";
                print filex "export JAVA_HOME=/hpc/packages/minerva-common/java/1.6.0_30/jdk1.6.0_30\n";
                print filex "export PATH=\"/hpc/packages/minerva-common/java/1.6.0_30/jdk1.6.0_30/bin:\$PATH\"\n";
                print filex "export _JAVA_OPTIONS=-Xmx1G\n";
                print filex "export SAMTOOLS=/hpc/users/wangm08/packages/samtools/samtools-0.1.18\n";
                print filex "export PATH=\"/hpc/users/wangm08/packages/samtools/samtools-0.1.18:\$PATH\"\n";
                print filex "cd /sc/orga/projects/zhuj05a/Wenhui/HBV/VirusFinder/VirusFinder2.0\n";
                print filex "./sys_check.pl\n";
		print filex "cd $workdir","unmapped/$fp/$fps[$j]\n";
		print filex "mkdir $label_dir\n";
		print filex "cd $label_dir\n";
		print filex "rm $workdir","unmapped/$fp/$fps[$j]/$label_dir/*\n";
		print filex "rm $workdir","$label_dir/$fp/$fps[$j]","_$label_dir.stdout\n";
                print filex "rm $workdir","$label_dir/$fp/$fps[$j]","_$label_dir.stderr\n";
		print filex "ln -s $dir_hg19_hbv_bwa_index/hg19+virus.fa hg19+virus.fa\n";
		print filex "ln -s $dir_hg19_hbv_bwa_index/hg19+virus.fa.amb hg19+virus.fa.amb\n";
		print filex "ln -s $dir_hg19_hbv_bwa_index/hg19+virus.fa.ann hg19+virus.fa.ann\n";
		print filex "ln -s $dir_hg19_hbv_bwa_index/hg19+virus.fa.bwt hg19+virus.fa.bwt\n";
		print filex "ln -s $dir_hg19_hbv_bwa_index/hg19+virus.fa.pac hg19+virus.fa.pac\n";
		print filex "ln -s $dir_hg19_hbv_bwa_index/hg19+virus.fa.sa hg19+virus.fa.sa\n";
		print filex "$vf_bin","detect_integration1_2.pl -c $cff -v $hbv_ref --fq1 $workdir","unmapped/$fp/$fps[$j]/unmapped.1.fq --fq2 $workdir","unmapped/$fp/$fps[$j]/unmapped.2.fq -o $workdir","unmapped/$fp/$fps[$j]/$label_dir\n";
                print filex $spc_cmd;
		close(filex);
                print file "bsub < $align_for_unmap_sub\n";
        }
	close(file);
}


sub sub_pipe_line_integ_combine{
        my $align_for_unmap= @_[0];
        my $workdir = @_[1];
        my $fp = @_[2];
        my $mem1 = @_[3];
        my $time1 = @_[4];
        my $nsf = @_[5];
        my $label_dir = @_[6];
        my @fps;

        for $i(7..scalar(@_)-1){
                $fps[$i-7]=@_[$i];
        }

        open(file,">",$align_for_unmap) or die"I can't write to the file!";
        print file "#!/bin/bash\n";
        print file "#BSUB -P acc_GTEX\n";
        print file "#BSUB -q $que\n";
        print file "#BSUB -J VF_comb2\n";
        print file "#BSUB -R \"rusage[mem=$mem1]\"\n";
        print file "#BSUB -W $time1\n";
        print file "#BSUB -m mothra\n";
        print file "#BSUB -o $workdir","$label_dir/$fp","_$label_dir.stdout\n";
        print file "#BSUB -eo $workdir","$label_dir/$fp","_$label_dir.stderr\n";
        print file "#BSUB -L /bin/bash\n";
	print file "cd $workdir","$label_dir/$fp\n";
	print file "rm -rf $workdir","$label_dir/$fp/*.bam\n";
	print file "rm -rf $workdir","$label_dir/$fp/*.bam.bai\n";
	print file "module load samtools\n";

        for my $j(0..$nsf-1){
		print file "ln -s $workdir","unmapped/$fp/$fps[$j]/integ_step1/alignment.sorted.bam alignment_$fp_$fps[$j].sorted.bam\n";
		print file "ln -s $workdir","unmapped/$fp/$fps[$j]/integ_step1/alignment.sorted.bam.bai alignment_$fp_$fps[$j].sorted.bam.bai\n";
        }
	print file "samtools merge alignment_$fp","_merge.tmp.sorted.bam *.bam\n";
#	print file "samtools sort alignment_$fp","_merge.tmp.sorted.bam alignment_$fp"."_merge.sorted\n";
#	print file "samtools index alignment_$fp","_merge.sorted.bam\n"; 
	print file "samtools view -h alignment_$fp","_merge.tmp.sorted.bam | sed 's/gi|21326584|ref|NC_003977.1|/chrVirus/g' | samtools view -bS - > alignment.sorted.bam\n";
        print file "samtools index alignment.sorted.bam\n";
        print file "rm alignment_$fp","_merge.tmp.sorted.bam\n";
	close(file);	
}

