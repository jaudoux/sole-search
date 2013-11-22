#!/usr/bin/perl
use strict; 
use warnings;
use Getopt::Std;
use vars qw($opt_a);
getopts('a:');

my $length = 150;
my $frag1 = 2000;
my $frag2 = 10000;
my $spacing1 = 2000;
my $spacing2 = 20000;
my $alpha = 0.001;

die "
Enter an input file to get its common peaks
options:
-a alpha value. options: 0.1, 0.05, 0.01, 0.001, 0.0001
 
" unless ($ARGV[0]);

if ($opt_a){
    $alpha = $opt_a;
}

###Part 1: bin input into 200bp bins and print raw data for user visualization
my $background = $ARGV[0];
binning ($background, "tempbin_$background\_test", "200");

###Part 2: Determine Duplication and deletion events
my ($bave) = duplifind ("tempbin_$background\_test", "20"); #blurs binned data by 20 bins and produces tempbin_$background\_test_temp file

# print this to a file so that people can see the data smeared.
printit ("tempbin_$background\_test_temp", "$background\_smear.sgr");

# call deletions
reverse_filter ("tempbin_$background\_test_temp", $frag2, $background);

# call large duplications
my $peak_co = file_count ("tempbin_$background\_test_temp", "tempbin_$background\_test_temp");
my $average = $bave/$peak_co;
my $duplicated = 3*($average); #duplicated regions require 3-fold higher than background
my ($filtered) = filter ("tempbin_$background\_test_temp", $duplicated);
get_peaks ($filtered, $frag2, $spacing2, "tempbin_$background\_peaks2.gff", $average);

# call narrow duplication events 
system ("perl rid_bgPeaks.pl tempbin_$background\_peaks2.gff tempbin_$background\_test >temp2_$background\_smear.sgr");
$peak_co = file_count ("temp2_$background\_smear.sgr", "temp2_$background\_smear.sgr");
$average = $bave/$peak_co;
my $non_unique = 4*($average);
($filtered) = filter ("temp2_$background\_smear.sgr", $non_unique); 
get_peaks ($filtered, $frag1, $spacing1, "temp_$background\_peaks1.gff", $average);
system ("cat temp_$background\_peaks1.gff tempbin_$background\_peaks2.gff | sort -k 1,1 -k 4n >$background\_duplications.gff");
system ("rm *temp*$background*");

#reduce the background file by the duplication events, smear again, but not as wide, then estimate enrichment over background.
binning ($background, "tempbin_$background\_test", "30");
duplic ("tempbin_$background\_test", "$background\_duplications.gff");
$bave = duplifind ("tempbin_$background\_test", "34"); #blur 1020 bases
$peak_co = file_count ("tempbin_$background\_test", "tempbin_$background\_test");
$average = $bave/$peak_co;
my ($sd) = stadev ("tempbin_$background\_test", $average, $peak_co);
my ($z) = getalpha ($alpha);
my $cut = ($z * $sd) + $average;
my $corrected = correct ("tempbin_$background\_test", $average, $background, $cut);
limitat ($background, $corrected, $length, $alpha); 
system ("rm *temp*$background*");


####################################################
###################SUBROUTINES######################
####################################################

sub limitat {
    my ($name, $file, $frag, $alpha) = @_;

    open (OT3, ">$name\_corrected_$alpha\.sgr");

    my @endvalue = ();
    my $val1;
    my @endef;
    my $chrom1 = "NA";
    my $start1 = 1;
    my $end1 = 1;
    my @pbins = ();
    my $line;

    open (IN, $file);
    
    while (<IN>){
	my ($chromer2, $starter2, $val2, $namer2) = split;
	my $ender2 = $starter2 + 30;
	if (($chromer2 ne $chrom1) or ($starter2 > $end1 + 1)){  #if the next chromosome is not the same as the previous or the next sequence is too far away from the previous, make the peak and/or continue
	    if (($end1 - $start1) >= $frag) { #if peak is long enough, keep
		foreach my $pot (@pbins){ #this is file with peak bins like peggy wanted
		    print OT3 "$pot\n";
		}
	    }
	    @endvalue = ($val2);
	    $chrom1 = $chromer2;
	    $start1 = $starter2;
	    $end1 = $ender2;
	    @pbins = ();
	    $line = "$chromer2\t$starter2\t$val2\t$namer2";
	    push (@pbins, $line);
	} elsif (($chromer2 eq $chrom1) and ($starter2 <= $end1 + 1)){
	    push (@endvalue, $val2);
	    $line = "$chromer2\t$starter2\t$val2\t$namer2";
	    push (@pbins, $line);
	    $end1 = $ender2;
	    $chrom1 = $chromer2;
	}
    }
    close OT3;
}


sub getalpha {
    my ($alpha) = @_;
    my $z = 0;
    
    if ($alpha eq 0.1){
	$z = 1.29;
    } elsif ($alpha eq 0.05){
	$z = 1.65;
    } elsif ($alpha eq 0.01){
	$z = 2.33;
    } elsif ($alpha eq 0.001){
	$z = 3.00;
    } elsif ($alpha eq 0.0001){
	$z = 3.75;
    }
    return ($z);
}
    
sub stadev {
    my ($file, $avg, $num) = @_;

    print "$avg\t$num\t";
    open (IN, $file);
    my $added = 0;
    
    while (<IN>){
	my ($chr, $st, $val) = split;
	$added += (($val - $avg)*($val - $avg));
    }
    my $var = $added/($num - 1);
    my $sd = sqrt($var);
    print "$var\t$sd\n";
    return ($sd);
}

sub printit {
    my ($input, $output) = @_;

    open (SM, $input);
    open (PR, ">temp");
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
    close SM;
}


sub correct {
    my ($file, $avg, $background, $cut) = @_;

    open (COF, $file);
    my $output = "temp_$background\_correct";
    open (PRI, ">temp_$background\_correct");

    while (<COF>){
	my ($chr, $start, $value) = split;
	if ($value > $cut) { # if this value is larger than expected by chance in this background file, we will want to correct the foreground data
	    my $corrected = sprintf "%.4f", ($value/$avg);
	    print PRI "$chr\t$start\t$corrected\t$background\n";
	}
    }
    close COF;
    close PRI;
    return ($output);
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

sub duplifind {
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
    my $bave = 0;

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
		$bave += $avg_val;
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
		    $bave += $avg_val;
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
			$bave += $avg_val;
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
    return ($bave);
}

sub file_count {
    my ($name, $out) = @_;
    my $test_count;
    
    system ("wc $name > tempcount_$out");
    my @file = open(INTEMP, "tempcount_$out");   
    while (<INTEMP>){
        my @counts = split;
        $test_count = $counts[0];
    }
    close INTEMP;
    system ("rm tempcount_$out");
    return $test_count;
}

sub filter {
    my ($file, $cut) = @_;
    open (ING, $file);
    my @windtemp = ();


    while (<ING>){
        my (@line) = split;
        my $ch = $line[0];
        my $sta = $line[1];
        my $sto = $sta + 199;
        my $va = $line[2];

        if ($va >= $cut){
            push @windtemp, {
                "chrom" => $ch,
                "start" => $sta,
                "stop" => $sto,
                "value" => $va,
                "effect" => 1,
            };
        } else {
            next;
        }
    }
    close ING;
    
    push @windtemp, {
        "chrom" => "buffer",
        "start" => 0,
        "stop" => 0,
        "value" => 0,
        "effect" => 1,
    };
    return \@windtemp;
}

sub reverse_filter {
    my ($file, $frag, $name) = @_;
    open (RF, $file);
    open (RFO, ">$name\_deletions.gff");
    my $last_end = 0;

    while (<RF>){
        my (@line) = split;
        my $ch = $line[0];
        my $sta = $line[1];
        my $sto = $sta + 199;
        my $va = $line[2];
        
        if (($sta > ($last_end + $frag)) && ($sta < ($last_end + 2900000))){
            print RFO "$ch\tSoleSearch_deleted\tDeletions\t$last_end\t$sta\t1\t+\t.\tDeletions_$ch\:$last_end\n";
        }
        $last_end = $sto;
    }
}

sub binning {
    my ($infile, $filename, $bin) = @_;

    open (SORT, $infile);
    my $chromb = "chr";
    my $blocks = 0;
    my $blocke = $bin - 1;
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
            $blocks = ((int($start/$bin)) * $bin);;
            $blocke = $blocks + ($bin - 1);
            $chromb = $chromo;
            $tag = 1;
	    open (OF3, ">$chromb\_$infile\_background.sgr");
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
		my $av = int(($blocks + $blocke)/2);
		print OF3 "$chromb\t$av\t$count\n";

                while (@ends && ($ends[0] <= $blocks + $bin)){
                    shift @sts;
                    shift @ends;
                }
                                #### we've finished analyzing this block. Time to move on to the next one
                $blocks += $bin;
                $blocke = $blocks + ($bin - 1);
            }
                                    # if after this analysis, the next tag is still larger, or in the next chromosome, and hasn't been included in the set, change the block coordinates
            if (($start > $blocke) or ($chromo ne $chromb)){
		if ($chromo ne $chromb){
		    close OF3;
                    open (OF3, ">$chromo\_$infile\_rawbinning.sgr");
                }

                $chromb = $chromo;
                $blocks = ((int($start/$bin)) * $bin);
                $blocke = $blocks + ($bin - 1);
            }
            push (@sts, $start);
            push (@ends, $end);
        }
    }
    close SORT;
    close OF2;
    close OF3;
}

sub get_peaks {
    my ($windows, $frag, $space, $filename, $bin_mean) = @_;
   
    open (OT2, ">$filename");
    
    my @endvalue = ($$windows[0] -> {value});
    my $chrom1 = $$windows[0] -> {chrom};
    my $start1 = $$windows[0] -> {start};
    my $end1 = $$windows[0] -> {stop};
    my $num = @$windows;
    my $pcount = 0;
    my $tagcount = 0;
    
   
    for (my $p = 1; $p < $num; ++$p){   #going through the bin file to make peaks
        my $chromer2 = $$windows[$p] -> {chrom};
        my $starter2 = $$windows[$p] -> {start};
        my $ender2 = $$windows[$p] -> {stop};
        my $val2 = $$windows[$p] -> {value};
        
        if (($chromer2 ne $chrom1) or ($starter2 > $end1 + $space)){  #if the next chromosome is not the same as the previous or the next sequence is too far away from the previous, make the peak and/or continue
            my $cut = ((($end1 - $start1)*2)/(600)); #want to have at least 1/3 of the bins to be over cutoff value.
            if (($end1 - $start1) >= $frag && (@endvalue >= $cut)) { #if peak is long enough, keep
                ++$pcount;
                my @fin = sort {$b <=> $a} @endvalue;
        
                my $totfn = 0;
                foreach my $posib (@fin){
                    $totfn += $posib;
                }
                my $avbin = $totfn/@fin;
                my $fold = sprintf "%.2f", ($avbin/$bin_mean);
                print OT2 "$chrom1\tSoleSearch_duplicated\tNonUnique_Region_fold\t$start1\t$end1\t$fold\t+\t.\t$chrom1\_$start1\_$end1\n";
                $tagcount += $fin[0];
            }
            @endvalue = ($val2);
            $chrom1 = $chromer2;
            $start1 = $starter2;
            $end1 = $ender2;
        } elsif (($chromer2 eq $chrom1) and ($starter2 <= $end1 + $space)){
            push (@endvalue, $val2);
            $end1 = $ender2;
            $chrom1 = $chromer2;
        }
    }
    return ($pcount, $tagcount);
}

