#!/usr/bin/perl

### perl script to parse exons from UCSC refFlat.txt file


if (!defined @ARGV) {
die  "Correct Syntax is: refflat_parse.perl <fileneame> \n\nPlease supply a file to parse.  The output file will be the input file name with \".intronBoundary\".exonsgaps\n\n";
}

$infile = shift @ARGV;

$outfile = $infile . ".intronBoundary.exonsgaps";


open INFILE, "<$infile";
open OUTFILE, ">$outfile";


while (<INFILE>) {
	chomp;
	if ($_ =~ /^\#/) {
		next;
	}
	
	($genename, $name, $chrom, $strand, $txStart, $txEnd, $cdsStart, 
	$cdsEnd, $exonCount, $exonStarts, $exonEnds)=split(/\t/,$_,11);
	
	@exstarts=split(/,/,$exonStarts);
	@exends=split(/,/,$exonEnds);
	
	$count=0;
	while ($count <$exonCount) {
		    $ex_num=$count+1;
			
			if ($ex_num == 1) {
					$introns[$count]="NA";
					$intronsStart[$count]="NA";
					$intronsEnd[$count]="NA";
				}
			if ($ex_num >1) {
					$introns[$count]=$exstarts[$count] - $exends[$count-1];
					$intronsStart[$count]=$exends[$count-1] + 1;
					$intronsEnd[$count]=$exstarts[$count] - 1;
				}
		
			print OUTFILE "$genename\t$name\t$chrom\t$strand\t$txStart\t$txEnd\t$cdsStart\t$cdsEnd\t$ex_num\t$exstarts[$count]\t$exends[$count]\t$introns[$count]\t$intronsStart[$count]--$intronsEnd[$count]\n";
					
			$count++;
	}
}
