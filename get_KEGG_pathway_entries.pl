# August 1st 2017
# Download all KEGG pathway entries related to Aspergillus niger

#! /usr/perl/bin -w
use strict;
use File::Fetch;

my $kegg_org_code="ang"; # organism code of A. niger in KEGG is "ang"
my $folder_pwy_entries="KEGG_Aniger_pwy_entries";
mkdir $folder_pwy_entries;
chdir $folder_pwy_entries;

my $pathway_list_url='http://rest.kegg.jp/list/pathway/'.$kegg_org_code;
my $file_pathway_list=File::Fetch->new(uri=>$pathway_list_url);
my $fileout_pathway_list=$file_pathway_list -> fetch() || die $file_pathway_list->error;

open(List,"<$fileout_pathway_list") || die "Cannot open file $file_pathway_list";
while (<List>)
{
	$_=~s/\s*$//;
	my @cols=split(/\t/,$_);
	my $pathway_id=$cols[0];
	$pathway_id=~s/path\://;
	$pathway_id=~s/\s*//g;
	my $entry_url='http://rest.kegg.jp/get/'.$pathway_id;
	my $file_entry=File::Fetch->new(uri=>$entry_url);
	my $fileout_entry=$file_entry->fetch() || die $file_entry->error;
}
close(List);

