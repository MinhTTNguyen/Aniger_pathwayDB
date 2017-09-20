# August 7th 2017
# Create pathologic file for A.niger NRRL3

#! /usr/perl/bin -w
use strict;

my $path="/home/mnguyen/Research/PathwayTools_v21";
my $filein_annotation="Annotation_table_Supplemental_Table_S6_NRRL3_functional_predictions_20_03_2017.txt";
my $filein_gff3="NRRL3_model_changes_Aug31_KM.gff3";
my $folderout="Aniger_NRRL3_7Aug2017_Pathologic";

#======================================================================================================================#
# read gff3 file
# get STARTBASE, ENDBASE, INTRON
# get chromosomal ID to determine the output file

my %hash_geneID_scaffoldID; # this hash contains all genes
my %hash_geneID_exons;
my %hash_geneID_start;
my %hash_geneID_end;
my $gene_ID="";
open(GFF,"<$path/$filein_gff3") || die "Cannot open file $filein_gff3";
while (<GFF>)
{
	$_=~s/\s*$//;
	if ($_!~/^\#/)
	{
		my @cols=split(/\t/,$_);
		my $chr=$cols[0];
		my $seq_feature=$cols[2];
		my $start=$cols[3];
		my $end=$cols[4];
		my $strand=$cols[6];
		my $id_info=$cols[8];
		
		# Get gene IDs. Some NRRL3 gene IDS have letters and other symbols: NRRL3_00001.1, NRRL3_00001a
		if ($seq_feature eq "gene")
		{
			$gene_ID=$id_info;
			$gene_ID=~s/^.+Name\=//;
			$gene_ID=~s/\;.*$//;

			$hash_geneID_scaffoldID{$gene_ID}=$chr;
			#unless ($chr) {print "\n$gene_ID";}
			if ($strand=~/\+/)
			{
				$hash_geneID_start{$gene_ID}=$start;
				$hash_geneID_end{$gene_ID}=$end;
			}else
			{
				$hash_geneID_start{$gene_ID}=$end;
				$hash_geneID_end{$gene_ID}=$start;
			}
		}elsif($seq_feature eq "exon")
		{
			my $exon_position="";
			if ($strand=~/\+/){$exon_position=$start."-".$end;}
			else{$exon_position=$end."-".$start;}
			
			if ($hash_geneID_exons{$gene_ID}){$hash_geneID_exons{$gene_ID}=$hash_geneID_exons{$gene_ID}.";".$exon_position;}
			else{$hash_geneID_exons{$gene_ID}=$exon_position;}
		}
	}
}
$gene_ID="";
close(GFF);
#======================================================================================================================#



#======================================================================================================================#
# create hash containing intron information
my %hash_geneID_introns;# this hash only contains gene harbouring introns
while(my ($k_geneID,$v_exons) = each (%hash_geneID_exons))
{
	#print "\n$k_geneID\t$v_exons\n";exit;
	my $gene_start=$hash_geneID_start{$k_geneID};
	my $gene_end=$hash_geneID_end{$k_geneID};
	if ($gene_start<$gene_end) # if gene is on plus strand
	{
		if ($v_exons=~/\;/) #multiple exons, i.e., gene contains at least one intron
		{
			my @exons=split(/;/,$v_exons);
			my $exon_number=scalar(@exons);
			
			my $first_exon=$exons[0];
			my $intron_start="";
			if ($first_exon=~/.+\-(.+)/)
			{
				my $first_exon_end=$1;
				$intron_start=$first_exon_end+1;
			}else{print "\nError (line ".__LINE__."): Position of the first exon is not as described!\nFirst exon position:$first_exon";exit;}
			
			for (my $i=1;$i<$exon_number;$i++)
			{
				my $exon_position=$exons[$i];
				if ($exon_position=~/(.+)\-(.+)/)
				{
					my $exon_start=$1;
					my $exon_end=$2;
					my $intron_end=$exon_start-1;
					my $intron_position=$intron_start."-".$intron_end;
					
					$intron_start=$exon_end+1;
					
					if ($hash_geneID_introns{$k_geneID}){$hash_geneID_introns{$k_geneID}=$hash_geneID_introns{$k_geneID}.";".$intron_position;}
					else{$hash_geneID_introns{$k_geneID}=$intron_position;}
				}else{print "\nError (line ".__LINE__."): Position of the exon is not as described!\nExon position:$exon_position";exit;}
			}
		}
	}else # if gene is on minus strand
	{
		my @exons=split(/;/,$v_exons);
		my $exon_number=scalar(@exons);
		my $first_exon=$exons[$exon_number-1];
		my $intron_start="";
		if ($first_exon=~/.+\-(.+)/)
		{
			my $first_exon_end=$1;
			$intron_start=$first_exon_end-1;
		}else{print "\nError (line ".__LINE__."): Position of the first exon is not as described!\nFirst exon position:$first_exon";exit;}
		
		for (my $i=$exon_number-2;$i>=0;$i--)
		{
			my $exon_position=$exons[$i];
			if ($exon_position=~/(.+)\-(.+)/)
			{
				my $exon_start=$1;
				my $exon_end=$2;
				my $intron_end=$exon_start+1;
				my $intron_position=$intron_start."-".$intron_end;
				
				$intron_start=$exon_end-1;
				
				if ($hash_geneID_introns{$k_geneID}){$hash_geneID_introns{$k_geneID}=$hash_geneID_introns{$k_geneID}.";".$intron_position;}
				else{$hash_geneID_introns{$k_geneID}=$intron_position;}
			}else{print "\nError (line ".__LINE__.")s: Position of the exon is not as described!\nExon position:$exon_position";exit;}
		}
	}
}
#======================================================================================================================#


#======================================================================================================================#
# read annotation table and print out of files
mkdir $folderout;
open(Annotation,"<$path/$filein_annotation") || die "Cannot open file $filein_annotation";
while (<Annotation>)
{
	if ($_!~/^\#/)
	{
		$_=~s/\s*$//;
		my @cols=split(/\t/,$_);
		my $gene_id=$cols[0];

		my $function=$cols[3];
		my $evidence=$cols[4];
		my $evidence_source=$cols[5];
		my $ecs=$cols[7];
		my $swissprot_id=$cols[11];
		my $aspGD_id=$cols[13];
		my $sgd_id=$cols[17];
		my $mycoCLAP_id=$cols[22];
		my $sp=$cols[32];

		my $chr=$hash_geneID_scaffoldID{$gene_id};
		my $fileout=$chr.".pf";
		open(Out,">>$path/$folderout/$fileout") || die "Cannot open file $fileout";

		my $gene_name=$gene_id; #use gene ID as the name of gene in PGDB would make it easier for searching genes
		print Out "ID\t$gene_id\n";
		print Out "NAME\t$gene_id\n";
		print Out "STARTBASE\t$hash_geneID_start{$gene_id}\n";
		print Out "ENDBASE\t$hash_geneID_end{$gene_id}\n";

		my $all_intron=$hash_geneID_introns{$gene_id};
		if ($all_intron)
		{
			if ($all_intron=~/\;/)
			{
				my @introns=split(/;/,$all_intron);
				foreach my $intron (@introns){print Out "INTRON\t$intron\n";}
			}else{print Out "INTRON\t$all_intron\n";}
		}

		print Out "PRODUCT-TYPE\tP\n";

		if (($function=~/hypothetical protein/) or ($function=~/uncharacterized protein/)){$function="ORF";}
		print Out "FUNCTION\t$function\n";
		if ($sp eq "Y"){print Out "FUNCTION-COMMENT\tSP:Y\n";}

		########################### DBLINKS ###########################
		if($evidence_source eq "CAZy"){print Out "FUNCTION-COMMENT\t$evidence\n";}
		elsif($evidence_source eq "Literature")
		{
			if ($evidence=~/\|/)
			{
				$evidence=~s/\s*//g;
				my @arr=split(/\|/,$evidence);
				foreach my $pubmed_id (@arr){print Out "FUNCTION-CITATION\t$pubmed_id\n";}
			}else{print Out "FUNCTION-CITATION\t$evidence\n";}
		}else
		{
			if ($evidence_source eq "SwissProt"){print Out "DBLINK\tSP:$swissprot_id\n";}
			elsif ($evidence_source eq "SGD"){print Out "DBLINK\tSGD:$sgd_id\n";}
			elsif($evidence_source eq "AspGD"){print Out "DBLINK\tAspGD:$aspGD_id\n";}
			elsif ($evidence_source eq "mycoCLAP"){print Out "DBLINK\tmycoCLAP:$mycoCLAP_id\n";}
			elsif($evidence_source eq "FungiDB"){print Out "DBLINK\t$evidence\n";}
			elsif($evidence_source eq "InterPro")
			{
				if ($evidence=~/\|/)
				{
					my @arr=split(/\|/,$evidence);
					foreach my $ipr (@arr)
					{
						$ipr=~s/^\s*//;$ipr=~s/\s*$//;
						$ipr=~s/\:.+$//;
						print Out "DBLINK\tINTERPRO:$ipr\n";
					}
				}else{$evidence=~s/\:.*$//;print Out "DBLINK\tINTERPRO:$evidence\n";}
			}elsif($evidence_source eq "PFAM")
			{
				if ($evidence=~/\|/)
				{
					my @arr=split(/\|/,$evidence);
					foreach my $pf (@arr)
					{
						$pf=~s/^\s*//;$pf=~s/\s*$//;
						$pf=~s/\:.+$//;
						print Out "DBLINK\tPFAM:$pf\n";
					}
				}else{$evidence=~s/\:.*$//;print Out "DBLINK\tPFAM:$evidence\n";}
			}else
			{
				$evidence=~s/\s*//g;
				if ($evidence)
				{
					if ($evidence=~/\|/)
					{
						my @arr=split(/\|/,$evidence);
						foreach my $cdd (@arr)
						{
							$cdd=~s/^\s*//;$cdd=~s/\s*$//;
							$cdd=~s/^CDD\://;
							print Out "DBLINK\tCDD:$cdd\n";
						}
					}else{$evidence=~s/^CDD\://;print Out "DBLINK\tCDD:$evidence\n";}
				}
			}
		}
		########################### DBLINKS ###########################
		
		
		########################### EC ###########################
		if ($ecs)
		{
			$ecs=~s/\s*//g;
			if ($ecs=~/\|/)
			{
				my @arr_ecs=split(/\|/,$ecs);
				foreach my $ec (@arr_ecs){print Out "EC\t$ec\n";}
			}else{print Out "EC\t$ecs\n";}
		}
		########################### EC ###########################
		
		print Out "//\n";
		close(Out);
	}
}
close(Annotation);
#======================================================================================================================#