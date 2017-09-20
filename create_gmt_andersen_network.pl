# August 2nd 2017
# create gmt file for Andersen's pathway network

#! /usr/perl/bin -w
use strict;

my $filein_andersen="Andersen_pathways_CBS_IDs.txt";
my $fileout_gmt="Andersen_pathways_CBS.gmt";

open(In,"<$filein_andersen") || die "Cannot open file $filein_andersen";
my %hash_pathway_cbsids;
while (<In>)
{
	$_=~s/^\s*//;
	$_=~s/\s*$//;
	my @cols=split(/\t/,$_);
	my $pathway_id=$cols[0];
	my $pathway=$cols[1];
	my $cbs_ids=$cols[2];
	$cbs_ids=~s/\s*//g;
	#$cbs_ids=~s/\;/\t/;
	my @arr_cbsids=split(/;/,$cbs_ids);
	$cbs_ids=join("\t",@arr_cbsids);
	$pathway=~s/^\s*//;
	$pathway=~s/\s*$//;
	
	my $pathway_uc=uc($pathway);
	my $geneset_name=$pathway_uc."%Andersen%".$pathway_id;
	my $pwy=$geneset_name."\t".$pathway;
	if ($hash_pathway_cbsids{$pwy}){$hash_pathway_cbsids{$pwy}=$hash_pathway_cbsids{$pwy}."\t".$cbs_ids;}
	else{$hash_pathway_cbsids{$pwy}=$cbs_ids;}
}
close(In);

open(Out,">$fileout_gmt") || die "Cannot open file $fileout_gmt";
while (my ($k, $v)=each(%hash_pathway_cbsids))
{
	my @all_ids=split(/\t/,$v);
	my %hash_temp;
	foreach my $id (@all_ids){$hash_temp{$id}++;}
	my @ids_nr=keys(%hash_temp);
	my $id_nr=join("\t",@ids_nr);
	#print "\n$k\n$id_nr\n";exit;
	print Out "$k\t$id_nr\n";
}
close(Out);