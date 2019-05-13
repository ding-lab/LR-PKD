#!/usr/bin/perl -w

open(IN, "$ARGV[0].results");
open(DA, "$ARGV[0].readcount");
open(OUT, ">$ARGV[0].output.tsv");
print OUT "Hugo_Symbol\tNCBI_Build\tChr\tStart\tEnd\tEffect\tType\tRef\tAlt\tNovelty\tNucleotide\tAminoAcid\tGene\tFeature\tFeature_type\tSIFT\tPolyPhen\tCLIN_SIG\tLONGRANGER_FILTER\tGT\tPS\tBX\tTotal_count\tRef_count\tAlt_count\tVAF\tLRPKD_FILTER\n";

while(<IN>)
{
	chomp;
	@l = split(/\t/,);
	$match = <DA>;
	chomp $match;
	@m = split(/\t/, $match);
	if($l[6] eq "SNP")
	{
		$ref=$l[7];
		$alt=$l[8];
	}elsif($l[6] eq "INS")
	{
		$ref=$m[2];
		$alt="+$l[8]";
	}elsif($l[6] eq "DEL")
	{
		$ref=$m[2];
		$alt="-$l[7]";
	}
	$nref="NA";
	$nalt="NA";
	for($i=5;$i<=$#m;$i++)
	{
		@d = split(/\:/,$m[$i]);
		if($d[0] eq $ref)
		{
			$nref=$d[1];
		}elsif($d[0] eq $alt)
		{
			$nalt=$d[1];
		}
	}
	$vaf=($nref eq "NA" || $nalt eq "NA") ? "NA" : sprintf("%.3f", $nalt/($m[3]))*100;
	if (index($l[19], "/") != -1 || $l[20] eq -1 || $nref eq "NA" || $nalt eq "NA" || $vaf < 20 || length(split(';', $l[21])) < 2) {
		$LRPKD_FILTER = "NOT_PASS"
	} else {
		$LRPKD_FILTER = "PASS"
	}
	print OUT "$_\t$m[3]\t$nref\t$nalt\t$vaf\%\t$LRPKD_FILTER\n";
}
