#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_j $opt_i $opt_f $opt_t $opt_c $opt_s $opt_p $opt_l $opt_r);
getopts("j:i:f:t:c:s:p:l:r");

my $iterations = 5;
my $frag = 150;
my $fdr = 0.001;
my $background2 = "K562_input_tags.txt";
my $overpeaks = "K562_inputs.txt_duplications.gff";
my $normalize = "K562_inputs.txt_corrected.txt";
my $lob = 1200;
my $join = 0;

die "
Enter tester file
options:
-j \t distance between peaks, allowing for the peaks to be joined. Default is $join
-r \t blur data when studying histone modifications.
-l \t blur length for histone modifications Default is $lob\. 
-i \t permutations. Default is $iterations
-f \t minimum chromatin fragment size(minimum peak size). Default is $frag. Lowest is 100bp. Highest is 250bp. Increase in increments of 30.
-t \t background file, including 30bp fragments that are actually sequenced, representing the sequenceable region of the genome
-c \t corrected file
-s \t fdr. Default is 0.001
-p \t gff file containing all peaks found in the background model

" unless ($ARGV[0]);

######################## MAIN ########################

my $name = $ARGV[0];

if ($opt_i){
    if ($opt_i > 10){
        die "must be 10 or less iterations at this time\n"
    }
    $iterations = $opt_i;
}

$join = $opt_j if $opt_j;
$background2 = $opt_t if $opt_t;
$frag = $opt_f if $opt_f;
$fdr = $opt_s if $opt_s;
$overpeaks = $opt_p if $opt_p;
$normalize = $opt_c if $opt_c;

#print out file that contains iterations of randomized, but biased datapoints, combine it with  our test file, sort, bin, and test bins for significance        
my $test_count = file_count ($name, $name); #test file count
my $tag_count = $test_count;
my $run = 2;
my $pcount = 0;
my $tagcount = 0;
my $sublob = 0;

binning ($name, "tempbin_$name\_test", "1");
duplic ("tempbin_$name\_test", $overpeaks);
system ("cat tempbin_$name\_test $normalize | sort -k 1,1 -k 2n > temp_freqs_$name");
system ("rm tempbin_$name\_test");
binsig ("temp_freqs_$name", $name); #outputs file $name\_temp_adjusted, which is to be used for the remainder of analysis
system ("rm temp_freqs_$name");

if ($opt_r){
    if ($opt_l){
	$lob = $opt_l;
    }
    $sublob = int($lob/30);
    blur ("$name\_temp_adjusted", $sublob);
}

#print the data after it has been manipulated
printit ("$name\_temp_adjusted", "$name\_adjusted.sgr");

my $bg_count = file_count ($background2, $background2);
my ($bg2) = random_chop ($tag_count, $bg_count, $background2, $name);
binning ($bg2, "tempbin_$name\_background"); #background2 is a file of only TAGS, not extended tags. It shows where in the genome can be sequenced.
system ("rm $bg2");
$bg_count = file_count ("tempbin_$name\_background", "tempbin_$name\_background"); #determine the number of bins into which your tags may fall in the background model.
my $co = 2;

#make ordered background file
open (BGF, "tempbin_$name\_background"); #take the bins and number them
open (BGFI, ">$background2\_count");
my $bgfc = 0;
while (<BGF>){
    my ($cc, $ss, $vv) = split;
    print BGFI "$bgfc\t$cc\t$ss\n";
    ++$bgfc;
}
close BGF;
close BGFI;
system ("sort $background2\_count >$background2\_count_temp");
system ("rm $background2\_count");

while ($run){
    if ($run eq 1){ # change to file with common peak regions removed
	system ("rid_bgPeaks.pl $name\_signifpeaks.gff $name >temp_num_$name");
	$tag_count = file_count ("temp_num_$name", "temp_num_$name");
	system ("rm temp_num_$name");
	$co -= 2;
	if ($co < 3){
	    $co = 3;
	}
    }

    my @signs = ("plus", "minus");
    for (my $k = 0; $k < $iterations; ++$k){
	open (OT, ">tempbinS_$name\_$k");
	for (my $j = 0; $j < $tag_count; ++$j){ #go through all of your tags and determine in which new bin they will reside
	    my $nu = int(rand($bg_count));
	    my $si = int(rand(@signs));
	    if ($signs[$si] eq "minus"){
		$nu -= 6;
	    }
	    print OT "$nu\n";
	}
	close OT;

	system ("sort tempbinS_$name\_$k > tempbinS2_$name\_$k");
	system ("mv tempbinS2_$name\_$k tempbinS_$name\_$k");
	system ("join tempbinS_$name\_$k $background2\_count_temp | awk \'{OFS = \"\t\"; print \$2, \$3, \$3+200}\' >tempbinS2_$name\_$k");
	system ("sort -k 1,1 -k 2n tempbinS2_$name\_$k > tempbinS_$name\_$k");
	system ("rm tempbinS2_$name\_$k");
	binning ("tempbinS_$name\_$k", "tempbin_$name\_$k");
	system ("rm tempbinS_$name\_$k");
	if ($opt_r){
	    blur ("tempbin_$name\_$k", $sublob);
	}
    }
    print "Background 1 files complete\n";
    
    #Now test to determine FDR
    my $confirm = 1;
    my $fdr_cut = 0;
    
    while ($confirm > $fdr_cut && (($fdr_cut - $confirm) < 1)){
	++$co;
	filter ("$name\_temp_adjusted", $co); #filters out bins that fall in acceptable category
	($pcount, $tagcount) = get_peaks ($name, "$name\_temp_adjusted_filtered.txt", $frag, "$name\_signifpeaks.gff", "final");   #This takes peaks from the hash data
	system ("rm $name\_temp_adjusted_filtered.txt");
	$fdr_cut = $pcount * $fdr;
	$confirm = $fdr_cut - 1;
	my $peaks = 0;
	my $counter = 0;
	my $coun = 0;
	
	while ($counter < $iterations){ #run through all ten rand files to check signif. If number of peaks in background excedes number of peaks alloted for any of the iterations, then go back and add one more tag to significance threshold.
	    filter ("tempbin_$name\_$counter", $co); #filters out bins that fall in acceptable category
	    ($confirm) = get_peaks ($name, "tempbin_$name\_$counter\_filtered.txt", $frag, "temp_$name", $fdr_cut);   #This takes peaks from the hash data. Returns # of peaks 
	    system ("rm tempbin_$name\_$counter\_filtered.txt");
	    ++$counter;
	    $coun += $confirm;
	    print "$co\t$pcount\t$fdr_cut\t$confirm\n";
	}
	$confirm = $coun/$iterations;
    }
    print "FDR determined cutoff: $co\n";
    --$run;
}

if ($opt_j){
    join_peaks ("$name\_signifpeaks.gff", $join);
}

system ("rm $name\_temp_adjusted temp*$name\*");
peak_stats ("$name\_signifpeaks.gff");
print "Summary file created and program complete\n";

############################################ SUBROUTINES ##########################################


sub join_peaks {
    my ($file, $join) = @_;
    my $chr = "NA";
    my $start = 0;
    my $end = 0;
    my @values = ();

    open (IN, $file);
    open (OUT, ">tempsmush_$file");

    while (<IN>){
	my (@line) = split;
	if ($line[0] ne $chr || ($line[3] > ($end + $join))){ #if the next peak is on a new chromosome or far away from the last peak, print the last peak in its entirety, if exists, and move on.
	    if ($values[0]){
		my @use = sort {$b <=> $a} @values;
		print OUT "$chr\t$line[1]\t$line[2]\t$start\t$end\t$use[0]\t+\t.\t$chr\_$start\n";
	    }
	    $chr = $line[0];
	    $start = $line[3];
	    $end = $line[4];
	    @values = ($line[5]);
	}
	if (($line[0] eq $chr) && ($line[3] <= ($end + $join))){ # if the next peak is very close to the last, join them.
	    $end = $line[4];
	    push (@values, $line[5]);
	}
    }
    close IN;
    close OUT;
    system ("mv tempsmush_$file $file");
}



sub printit {
    my ($input, $output) = @_;

    open (SM, $input);
    open (PR, ">temp_$input\_removed");
    my $chr = "NA";

    while (<SM>){
	my ($ch, $st, $va) = split;
	if ($ch ne $chr){
	    close PR;
	    open (PR, ">$ch\_$output");
	    $chr = $ch;
	}
	print PR "$ch\t$st\t$va\n";
    }
    close PR;
    system ("rm temp_$input\_removed");
}

sub random_chop {
    my ($test_count, $bg_count, $bg, $name) = @_;

    my $freq = $test_count/$bg_count;
    open (IN, $bg);
    my $bname = "tempsort_$name";
    open (O1, ">$bname");

    while (<IN>){
	my $line = $_;
	chomp $line;
	if (rand(1) < $freq){
	    print O1 "$line\n";
	}
    }
    close IN;
    close O1;
    return ($bname);
}


sub binsig {
    my ($infile, $name) = @_;

    open (SORT, $infile);
    open (OUTER, ">$name\_temp_adjusted");
    my $chromb = "chr";
    my $blocks = 0;
    my $blocke = 29;
    my $tester = 1;
    my $fval;
    
    while (<SORT>){
	my (@sorted) = split;
	my $chromo = $sorted[0];
	my $start = $sorted[1];
	my $value = $sorted[2];
	if ($sorted[3]){
	    $tester = $value;
	} else {
	    $fval = $value/$tester;
	    print OUTER "$chromo\t$start\t$fval\n";
	    $tester = 1;
	}
    }
    close SORT;
    close OUTER;
}
	

sub duplic {
    my ($bin, $fold) = @_;

    open (OUT, ">2_$bin");
    open (CHR, $fold);
    my @peaks = <CHR>;
    close CHR;
    my ($pc, $d, $d1, $ps, $pe, $pv) = split ("\t", $peaks[0]);

    
    open (BINN, $bin);

    while (<BINN>){
	my ($chr, $start, $value) = split;
	if ($chr eq $pc && $pe < $start && $peaks[0]){
	    shift @peaks;
	    if ($peaks[0]){
		($pc, $d, $d1, $ps, $pe, $pv) = split ("\t", $peaks[0]);
	    } else {
		$pc = "chr";
		$ps = 0;
		$pe = 0;
		$pv = 0;
	    }
	}

	if ($chr eq $pc && $ps <= $start && $pe >= $start){
	    $value /= $pv;
	}
	print OUT "$chr\t$start\t$value\n";
    }
    close CHR;
    close BINN;
    close OUT;
    system ("mv 2_$bin $bin");
}

sub peak_stats {
    my ($file) = @_;

    open (IN, $file);
    my $peak_count = file_count ($file, $file);
    my @heights = ();
    my $hit = 0;
    my $wide = 0;
    
    while (<IN>){
	my ($chr, $field1, $field2, $start, $end, $value, $sign, $dot, $field3) = split;
	$wide += ($end - $start);
	$hit += $value;
	push (@heights, $value);
    }
    close IN;
    
    my @sheight = sort {$a <=> $b} @heights;
    my $mid = int($peak_count/2);
    my $low = $sheight[0];
    my $high = $sheight[$peak_count-1];
    my $median = $sheight[$mid];
    my $peak_avg = $hit/$peak_count;
    my $width = $wide/$peak_count;
    my $speak_avg = sprintf "%.2f", $peak_avg;
    my $swidth = sprintf "%.2f", $width;
    
    $file =~ s/.gff//;
    open (OUT, ">Summary_$file\.txt");
    print OUT "Number of peaks:\t$peak_count\nAverage peak height:\t$speak_avg\nMedian peak height:\t$median\nHighest peak:\t$high\nLowest peak:\t$low\nAverage peak width:\t$swidth\n";
    close OUT;
}

sub filter {
    my ($file, $cut) = @_;
    open (FIL, $file);
    
    open (IN, ">$file\_filtered.txt");
    
    while (<FIL>){
	my $fili = $_;
	chomp $fili;
	my ($ch, $sta, $va) = split ("\t", $fili);
	
	if ($va >= $cut){
	    print IN "$ch\t$sta\t$va\n";
	} else {
	    next;
	}
    }
    close FIL;
    close IN;
}
    

sub binning {
    my ($infile, $filename, $type) = @_;

    open (SORT, $infile);
    my $chromb = "chr";
    my $blocks = 0;
    my $blocke = 29;
    my @sts = ();
    my @ends = ();
    my $tag = 0;
    my $count = 0;

    open (OF2, ">$filename");

    while (<SORT>){
	my (@sorted) = split;
	my $chromo = $sorted[0];
	my $start = $sorted[1];
	my $end = $sorted[2];

	        #if just starting, initiate bin
	unless ($tag) {
	    $blocks = ((int($start/30)) * 30);;
	    $blocke = $blocks + 29;
	    $chromb = $chromo;
	    $tag = 1;
	    if ($type){
		open (OF3, ">$chromb\_$infile\_rawbinning.sgr");
	    }
	}

	        # if the tag is already in the initiated bin, add it
	if (($chromo eq $chromb) and ((($blocks >= $start) and ($blocks <= $end)) or (($blocke >= $start) and ($blocke <= $end)))){
	    push (@sts, $start);
	    push (@ends, $end);

	                # elsif tag is outside initiated bin, evaluate eveything thus far and create new bin
	} elsif (($start > $blocke) or ($chromo ne $chromb)){
	    while (($start > $blocke || $chromo ne $chromb) && @sts){ #
		$count = scalar @sts;
		print OF2 "$chromb\t$blocks\t$count\n";
		if ($type){
		    my $av = int(($blocks + $blocke)/2);
		    print OF3 "$chromb\t$av\t$count\n";
		}
		 		    
		while (@ends && ($ends[0] <= $blocks + 30)){
		    shift @sts;
		    shift @ends;
		}
                #### we've finished analyzing this block. Time to move on to the next one
		$blocks += 30;
		$blocke = $blocks + 29;
	    }
	                # if after this analysis, the next tag is still larger, or in the next chromosome, and hasn't been included in the set, change the block coordinates
	    if (($start > $blocke) or ($chromo ne $chromb)){
		if ($chromo ne $chromb && $type){
		    close OF3;
		    open (OF3, ">$chromo\_$infile\_rawbinning.sgr");
		}
		
		$chromb = $chromo;
		$blocks = ((int($start/30)) * 30);
		$blocke = $blocks + 29;
	    }
	    push (@sts, $start);
	    push (@ends, $end);
	}
    }
    close SORT;
    close OF2;
    if ($type){
	close OF3;
    }
}

sub blur {
    my ($file, $blur) = @_;
    open (IN, $file);
    open (OUT, ">$file\_temp");
    my $ch2 = 0;
    my $st2 = 0;
    my $va2 = 0;
    my $avg_val = 0;
    my $num = -1;
    my $num2 = 0;
    my @trials = ();

    my $half_blur = int($blur/2);
    
    while (<IN>){
	my ($line) = $_;
	++$num;
	my $values = 0;
	my ($ch3, $st3, $va3) = split ("\t", $line);

	if ($num > $half_blur){
	    if ($num < $blur){
		push (@trials, $line);
		my ($chr, $start, $value) = split ("\t", $trials[$num2]);
		foreach my $fi (@trials){
		    my ($chr, $st, $va) = split ("\t", $fi);
		    $values += $va;
		}
		$avg_val = int($values/@trials);
		print OUT "$chr\t$start\t$avg_val\n";
		++$num2;
	    } else {
		my ($ch2, $st2, $va2) = split ("\t", $trials[$half_blur]);
		if ($ch3 eq $ch2 || $st3 <= ($st2 + ($half_blur * 30 + 1000))){
		    push (@trials, $line);
		    foreach my $fi (@trials){
			my ($chr, $st, $va) = split ("\t", $fi);
			$values += $va;
		    }
		    $avg_val = int($values/@trials);
		    print OUT "$ch2\t$st2\t$avg_val\n";
		    shift @trials;
		} elsif ($ch3 ne $ch2 || ($ch3 eq $ch2 && $st3 > ($st2 + ($half_blur * 30 + 1000)))){
		    while ($trials[$half_blur]){
			my ($chr1, $start1, $value1) = split ("\t", $trials[$half_blur]);
			foreach my $fi (@trials){
			    my ($chr, $st, $va) = split ("\t", $fi);
			    $values += $va;
			}
			$avg_val = int($values/@trials);
			print OUT "$chr1\t$start1\t$avg_val\n";
			shift @trials;
			$values = 0;
		    }
		    $num = 1;
		    $num2 = 0;
		    @trials = ($line);
		}
	    }
	} else {
	    push (@trials, $line);
	}
    }
    close IN;
    close OUT;
    system ("mv $file\_temp $file");
}


sub get_peaks {
    my ($name, $proc, $frag, $filename, $type) = @_;
    open (OT, ">$filename");
    open (OTS, ">Summet_$filename");

    my $redname;
    my $finish = 0;
    if ($type and $type eq "final") {
	$finish = 1;
	$type = 0;
	my $redname = "$name\_redbin.sgr";
	open (OT3, ">$redname");
    }

    my @endvalue = ();
    my @mids = ();
    my $middle;
    my $val1;
    my $chrom1;
    my $start1;
    my $end1;
    my $pcount = 0;
    my $tagcount = 0;
    my @pbins = ();
    my $a1;
    my $line;
    
    open (IN, $proc);
    my $marker = 1;
    
    while (<IN>) {
	my $totline = $_;
	chomp $totline;
	my @windows = split ("\t", $totline);
	if ($marker){
	    $chrom1 = $windows[0];
	    $start1 = $windows[1];
	    $end1 = $windows[1] + 29;
	    $middle = int(($start1 + $end1)/2);
	    $val1 = $windows[2];
	    $a1 = int(($start1 + $end1)/2);
	    push (@endvalue, $val1);
	    push (@mids, $middle);
	    push (@pbins, $totline);
	    $marker = 0;
	    next;
	} else {
	    my $chromer2 = $windows[0];
	    my $starter2 = $windows[1];
	    my $ender2 = $windows[1] + 29;
	    my $val2 = $windows[2];
	    my $a2 = int(($starter2 + $ender2)/2);
	    
	    if (($chromer2 ne $chrom1) or ($starter2 > $end1 + 1)){  #if the next chromosome is not the same as the previous or the next sequence is too far away from the previous, make the peak and/or continue
		if (($end1 - $start1) >= $frag) { #if peak is long enough, keep
		    ++$pcount;
		    unless ($finish){
			if ($type && $pcount > $type){ #in case of bg, return if more peaks than signif
			    return $pcount;
			}
		    }
		    my $top = 0;
		    my $topnum = -1;
		    foreach my $vs (@endvalue){
			if ($vs > $top){
			    $top = $vs;
			    ++$topnum;
			}
		    }
		    my $summit = $mids[$topnum];
		    print OT "$chrom1\tSole_search_$name\tTagPeak\t$start1\t$end1\t$top\t+\t.\t$chrom1\_$start1\_$end1\n";
		    print OTS "$chrom1\tSole_search_summit_$name\tSummet\t$summit\t$summit\t$top\t+\t.\t$chrom1\_$start1\_$end1\n";
		    if ($finish){
			foreach my $pot (@pbins){ #this is file with peak bins like peggy wanted
			    print OT3 "$pot\n";
			}
		    }
		    $tagcount += $top;
		}
		@endvalue = ($val2);
		$chrom1 = $chromer2;
		$start1 = $starter2;
		$end1 = $ender2;
		$middle = int(($start1 + $end1)/2);
		@mids = ($middle);
		$line = "$chromer2\t$a2\t$val2";
		@pbins = ($line);
	    } elsif (($chromer2 eq $chrom1) and ($starter2 <= $end1 + 1)){
		push (@endvalue, $val2);
		$middle = int(($starter2 + $ender2)/2);
		push (@mids, $middle);
		$line = "$chromer2\t$a2\t$val2";
		push (@pbins, $line);
		$end1 = $ender2;
	    }
	}
    }
    close OT;
    close IN;
    close OTS;
    if ($finish){
	close OT3;
    }
    return ($pcount, $tagcount);
}

sub file_count {
    my ($name, $out) = @_;
    my $test_count;

    my $rn = rand();
    system ("wc $name > temp$rn\_$out");
    my @file = open(INTEMP, "temp$rn\_$out");   
    while (<INTEMP>){
	my @counts = split;
	$test_count = $counts[0];
    }
    close INTEMP;
    system ("rm temp$rn\_$out");
    return $test_count;
}


