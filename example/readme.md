## Yeast example 
Example to test the scanPAV installation using yeast genome manually modified to include long indels plus snps for noise.

### Files:
  **yeast_data/S288C_reference_sequence_R64-1-1_20110203.fsa**   Yeast un-modified reference
  **yeast_data/S288c_with_indels_and_snps.fasta**   	     Yeast reference modified to include long indels and snps

  __yeast_data/indels_added.list__ 	  Long indels added to the modified assembly, simulated using the pipeline https://github.com/mlliou112/simulatesv 
  __yeast_data/snps_added.list__      SNPs added to the modified assembly, simulated using the pipeline https://github.com/mlliou112/simulatesv

#### Check if results are correct:
Compare your PAVs with the file in folder:

  __yeast_data/extracted_pavs__:   Folder that include the presence and absence PAVs from a properly working scanPAV run
  
	pavs_absent_in_S288c_with_indels_and_snps,   pavs_absent_in_S288c_with_indels_and_snps.fasta : vcf and fasta file for absence PAVs, corrisponding to deletions added by the simulation 

	pavs_present_in_S288c_with_indels_and_snps,  pavs_present_in_S288c_with_indels_and_snps.fasta : vcf and fasta file for presence PAVs, corrisponding to insertions added by the simulation


### Results:
  All the insertions and deletions are found by ScanPAV, but there are 2 false positives:
    False positive absence PAV:
	PAV_ref|NC_001140|_000381001_00001000 : occurring on the edge of a long insertion, part of the sequence is original, part is an insertion, bwa cannot map it completely
    False positive presence PAV:
	PAV_SV_ref|NC_001136|_000838001_00001000 : occurring on the edge of a long deletion, part of the sequence is from the original position, part is from after the deletion length,  bwa cannot map it as a unique sequence

    Both cases are two 1000 bp long PAVs, and they occurr on the edge of a real PAV. As mentioned on the main Readme.md, we recommend to filter out any PAV <= 1000 bp, as they are mostly noise.

