# August 2nd 2017
# Replace CBS IDs by NRRL3 IDs in gmt file
# The corresponding list was from Erin

#! /usr/perl/bin -w
use strict;

my $file_gmt="AspGD_GO_Aniger_CBS_BP_no_IEA_ND_RCA_with_proteinID.gmt";
my $file_corresponding_IDs="corresponding_NRRL3IDs_of_CBS_IDs.txt";
my $fileout="AspGD_GO_Aniger_NRRL3_BP_no_IEA_ND_RCA.gmt";

#=============================================================================================#
open(Corresponding,"<$file_corresponding_IDs") || die "Cannot open file $file_corresponding_IDs";
my %hash_cbs_nrrl3;
while (<Corresponding>)
{
	$_=~s/\s*$//;
	my @cols=split(/\t/,$_);
	my $cbs_id=$cols[0];
	my $nrrl3_id=$cols[1];
	if ($nrrl3_id)
	{
		if ($hash_cbs_nrrl3{$cbs_id}){$hash_cbs_nrrl3{$cbs_id}=$hash_cbs_nrrl3{$cbs_id}."\t".$nrrl3_id;}
		else{$hash_cbs_nrrl3{$cbs_id}=$nrrl3_id;}
	}
}
close(Corresponding);

### filter duplicated nrrl3 ids of each cbs ids
while (my ($k, $v)=each(%hash_cbs_nrrl3))
{
	my %temp_hash;
	my @nrrl3_ids=split(/\t/,$v);
	foreach my $nrrl3_id (@nrrl3_ids){$temp_hash{$nrrl3_id}++;}
	my @nrrl3_ids_nr=keys(%temp_hash);
	my $nrrl3_id_nr=join("\t",@nrrl3_ids_nr);
	$hash_cbs_nrrl3{$k}=$nrrl3_id_nr;
}
#=============================================================================================#



#=============================================================================================#
open(GMT,"<$file_gmt") || die "Cannot open file $file_gmt";
open(Out,">$fileout") || die "Cannot open file $fileout";
while (<GMT>)
{
	$_=~s/\s*$//;
	my @cols=split(/\t/,$_);
	my $pathway_name=shift(@cols);
	my $pathway_desc=shift(@cols);
	my %hash_nrrl3_ids;
	foreach my $cbs_id (@cols)
	{
		my $corr_nrrl3_ids=$hash_cbs_nrrl3{$cbs_id};
		my @arr_nrrl3_ids=split(/\t/,$corr_nrrl3_ids);
		foreach my $each_nrrl3 (@arr_nrrl3_ids){$hash_nrrl3_ids{$each_nrrl3}++;}
	}
	my @arr_nrrl3_ids_nr=keys(%hash_nrrl3_ids);
	my $nrrl3_ids_nr=join("\t",@arr_nrrl3_ids_nr);
	print Out "$pathway_name\t$pathway_desc\t$nrrl3_ids_nr\n";
}
close(GMT);
close(Out);
#=============================================================================================#

