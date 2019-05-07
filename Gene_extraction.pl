#!/usr/bin/perl -w
#use strict;
#use warnings;

# Variants in coding region
open(IN, "<Gene_interest.txt");
while(<IN>)
{
	chomp;
	$h1{$_} = "";
}

open(DA, "<input.vep.GTPSBX.maf");
open(OUT1, ">Genes_interest_coding.results");
<DA>;<DA>;
while(<DA>)
{
	chomp;
	@l = split(/\t/,);
	next if ($l[8] eq "Splice_Region" || $l[8] eq "3'Flank" || $l[8] eq "3'UTR" || $l[8] eq "5'Flank" || $l[8] eq "5'UTR" || $l[8] eq "IGR" || $l[8] eq "Intron" || $l[8] eq "RNA" || $l[8] eq "Silent");
	if(exists $h1{$l[0]})
	{
		print OUT1 "$l[0]\t$l[3]\t$l[4]\t$l[5]\t$l[6]\t$l[8]\t$l[9]\t$l[10]\t$l[12]\t$l[13]\t$l[34]\t$l[36]\t$l[47]\t$l[48]\t$l[49]\t$l[71]\t$l[72]\t$l[85]\n";
	}
}

close IN;
close OUT1;
close DA;

# All types of variants
open(IN, "<Gene_interest.txt");
while(<IN>)
{
        chomp;
        $h1{$_} = "";
}

open(DA, "<input.vep.GTPSBX.maf");
open(OUT1, ">Genes_interest_all.results");
<DA>;<DA>;
while(<DA>)
{
        chomp;
        @l = split(/\t/,);
        if(exists $h1{$l[0]})
        {
                print OUT1 "$l[0]\t$l[3]\t$l[4]\t$l[5]\t$l[6]\t$l[8]\t$l[9]\t$l[10]\t$l[12]\t$l[13]\t$l[34]\t$l[36]\t$l[47]\t$l[48]\t$l[49]\t$l[71]\t$l[72]\t$l[85]\n";
        }
}

close IN;
close OUT1;
close DA;
