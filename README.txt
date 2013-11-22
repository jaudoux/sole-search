Sole-Search v2

Two programs are required to determine statistically significant peaks using ChIP-seq data and Input data.

First, all raw ChiP data is combined into one file and raw Input data is combined into another file, and each are processed. Two files are produced: one represents the coordinates of the tags themselves (actual sequence that was read by the GA2) and the other represents the entire chromatin fragment whose end was sequenced (This is done by parsing large raw data files for uniquely mapped reads, obtaining the start coordinate, and adding 200bp in the orientation the tag was being sequenced. A shell program is included to show how this parsing occurs when starting with an eland_extended.txt file.
ex:
sh parse_export.txt <export_file> <output_name>

Files that come from this are:
Input_full-length.txt
Input_tags.txt
ChIP-File_full-length.txt
ChIP-File_tags.txt (****This file is not needed and can be eliminated from processing if desired****)

Program 1: normalize_sd.pl
This program processes the Input data (the file containing chr, start, start+200)
The available option is alpha value. In this case, the alpha value is used to determine regions of the input that are statistically enriched, in comparison to the rest of the background. If a region is enriched over background, it will be used later to normalize the ChIP-seq data. Therefore, a smaller alpha value will result in more peaks in the ChIP-seq (less of the test data has to be normalized) and a large alpha value will result in fewer peaks due to more regions being normalized. Alpha values that can be used are 0.1, 0.05, 0.01, 0.001, and 0.0001.
File that is used for input is:
Input_full-length.txt

ex. 
perl normalize_sd.pl -a 0.001 Input_full-length.txt

options:
-a alpha value. options: 0.1, 0.05, 0.01, 0.001, 0.0001

recommendations:
An alpha value of 0.01 or 0.001 is recommended for most transcription factors and an alpha value of 0.001 or 0.0001 is recommended for histone modifications or spreading factors.

Files that come from this program are:
1. Raw data visualization file:
chr*_Input_full-length.txt_rawbinning.sgr

2. Processed, smeared data visualization file
chr*_Input_full-length.txt_smear.sgr

3. gff file containing genomic duplication events
Input_full-length.txt_duplications.gff

4. gff file containing genomic deletion events
Input_full-length.txt_deletions.gff

5. corrected.txt file that includes bins with fold increase over average background.
Inputs_full-length.txt_corrected_0.0001.sgr

**The duplications file and corrected file will be used in the second program: the ChIP-seq data will be normalized by the fold increase of that region in the input genome.

Program 2: Sole-searchV2.pl
This program determines peaks from ChIP-seq data.
It requires four files corresponding to the background:
1. Input_tags.txt corresponding to the coordinates sequenced in the run
2. duplications file (Input_full-length.txt_duplications.gff)
3. corrected file (Inputs_full-length.txt_corrected_0.0001.sgr)
4. Raw extended tag data for ChIP (ChIP-File_full-length.txt)

ex.
perl Sole-searchV2_final.pl -t Input_tags.txt -c Input_full-length.txt_corrected_0.0001.sgr -s 0.01 -p Input_full-length.txt_duplications.gff ChIP-File_full-length.txt

options:
-j 	 distance between peaks, allowing for the peaks to be joined. Default is 0
-r 	 blur data when studying histone modifications.
-l 	 blur length for histone modifications Default is 1200. 
-i 	 permutations. Default is 5
-f 	 minimum chromatin fragment size(minimum peak size). Default is 150. Lowest is 100bp. Highest is 250bp. Increase in increments of 30.
-t 	 background file, including 30bp fragments that are actually sequenced, representing the sequenceable region of the genome
-c 	 corrected file
-s 	 fdr. Default is 0.001
-p 	 gff file containing all peaks found in the background model

recommendations:
When studying factors with >10000 binding sites, we recommend setting FDR to 0.001 or 0.0001. When studying histone modifications or if background is high, we recommend and FDR of 0.01 (or 0.1 in cases of high background or lots of spreading.)
When studying histone modifications, we suggest using the blur (-r) option and leave default blur length. We also suggest using the peak-joining option (-j) with a length of 300bp.
example for running the program for a histone modification:
perl Sole-searchV2_final.pl -t Input_tags.txt -c Input_full-length.txt_corrected_0.0001.sgr -s 0.01 -p Input_full-length.txt_duplications.gff -j 300 -r ChIP-File_full-length.txt

The program returns the following files:
1. Raw data for visualization (.sgr)
2. visualization file showing the ChIP-seq data after it has been normalized (adjusted.sgr)
3. peaks file with tag counts (signifpeaks.gff)
4. peak statistics summary file (Summary_)
