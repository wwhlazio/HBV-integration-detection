#!/usr/bin/perl
use strict;
use warnings;

my $dir="/sc/orga/projects/zhuj05a/Wenhui/HBV/script/wgs_88/";
my $sample="22";
my $phe="tumor";
my $step="d-align";

my $dir_cur=$dir.$sample."/".$phe."/virusfinder1/output/step3/$step/crest/";
my $cmd="ls $dir_cur"."split/*";
my $tmp=`$cmd`;
my @tmp1=split(/\n/,$tmp);
my $nsp=0;
my @sp;
####################################################remove the old file in crest/##########################
my $out1=$dir_cur."results.predSV.txt";
my $out2=$dir_cur."results.predSV.txt-filter";
my $out3=$dir_cur."results.predSV.txt2";

if(-s $out1){
	`rm $out1`;
}
if(-s $out2){
	`rm $out2`;
}
if(-s $out3){
	`rm $out3`;
}
###################################################################################################################



###################################################check the result of split  if some thing wrong will generate error meeesage#####################

my $all_sh=$dir_cur."split/split_lsf_redo.sh";
my $comb_sh=$dir_cur."split/combin_sub_crest.sh";
open(my $FILE,">",$all_sh) or die"I can't write to the file!";
open(my $FILE1,">",$comb_sh) or die"I can't write to the file!";
while(<@tmp1>){
	my $stdoutfile=$_;
	my @a=split(/\./,$_);
	if($a[$#a] eq "stdout"){
		my @b=split(/\//,$a[$#a-1]);
		my @c=split(/\_/,$b[$#b]);	
		$sp[$nsp]=$c[$#c];
		my $cmd="grep \"Successfully completed\" $stdoutfile";
		my $tmp=`$cmd`;
		if(!length($tmp)){
			print "ERROR:please check split $nsp!";
			print $FILE "bsub < $dir_cur","split/$step","_crest_sub_$nsp.lsf\n";	
		}
		else{	
			my $in1=$dir_cur."split/split$sp[$nsp]/results.predSV.txt";
			my $in2=$dir_cur."split/split$sp[$nsp]/results.predSV.txt-filter";
			my $in3=$dir_cur."split/split$sp[$nsp]/results.predSV.txt2";
			print $in1,"\n";
			print $in2,"\n";
			print $in3,"\n";
			if(-s $in1){
				print $FILE1 "cat  $in1 >> $out1\n";
			}
			if(-s $in2){
                                print $FILE1 "cat  $in2 >> $out2\n";
                        }
			if(-s $in3){
                                print $FILE1 "cat  $in3 >> $out3\n";
                        }
		}
		$nsp++;
	}	
}
close($FILE);
close($FILE1);
##################################################################################################################
