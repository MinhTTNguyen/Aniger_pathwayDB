# August 3rd 2017
# Create gmt file from GO assigned by AspGD

#! /usr/perl/bin -w
use strict;

my $filein_GO="AspGD_GO_Aniger_CBS_BP_no_IEA_ND_RCA_with_proteinID.txt";
my $fileout_gmt="AspGD_GO_Aniger_CBS_BP_no_IEA_ND_RCA_with_proteinID.gmt";

open(In,"<$filein_GO") || die "Cannot open file $filein_GO";
my %hash_geneset;
while (<In>)
{
	$_=~s/\s*$//;
	my @cols=split(/\t/,$_);
	my $geneid=$cols[2];
	my $goid=$cols[5];
	my $goterm=$cols[6];
	my $goterm_uc=uc($goterm);
	my $geneset_name=$goterm_uc."%GOBP%".$goid."\t".$goterm;
	if($hash_geneset{$geneset_name}){$hash_geneset{$geneset_name}=$hash_geneset{$geneset_name}."\t".$geneid;}
	else{$hash_geneset{$geneset_name}=$geneid;}
}
close(In);
open(Out,">$fileout_gmt") || die "Cannot open file $fileout_gmt";
while (my ($k, $v)=each(%hash_geneset))
{
	my @geneids=split(/\t/,$v);
	my %hash_temp;
	foreach my $id (@geneids){$hash_temp{$id}++;}
	my @geneids_nr=keys(%hash_temp);
	my $gene_ids=join("\t",@geneids_nr);
	print Out "$k\t$gene_ids\n";
}
close(Out);