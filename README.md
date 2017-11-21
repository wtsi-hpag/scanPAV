# scanPAV
Pipeline to detect PAVs (presence/absence variations) in genome comparison using whole genome alignment.

Description

scanPAV is a pipeline to detect PAVs (presence/absence variations) using whole genome alignment. 

Say if you have two assemblies assembly_1 and assembly_2:

You want to find sequences presented in assembly_1 and missed in assembly_2:

This means assembly_1 = presence_assembly
           assembly_2 = absence_assembly

Run scanPAV

Using BWA (by default) 

/tmp/scanPAV/scanPAV -nodes 30 -align bwa -score 550 presence_assembly.fasta absence_assembly.fasta presence_absence.fasta > try.out

       nodes  (30)    - number of CPUs requested
       score  (550)   - smith-waterman alignment score 
       align  (bwa)   - using BWA as the sequence aligner 

Or using SMALT (https://sourceforge.net/projects/smalt/)
 
/tmp/scanPAV/scanPAV -nodes 30 -align smalt -score 550 presence_assembly.fasta absence_assembly.fasta presence_absence.fasta > try.out

This gives you the PAV file presence_absence.fasta in your working directory


Note: 1. you need to give the full path of the scanPAV input files or the files are in your working directory;
      2. presence_assembly.fasta and absence_assembly.fasta should be in your working directory;
      3. presence_assembly.fasta and absence_assembly.fasta are in full path.
      4. If you use bwa, you need to check if the binary verson provided in the package works:
         Type
         /tmp/scanPAV/scanPAV/scanPAV-bin/bwa
         If there are error messages, please copy a local working version to 
         /tmp/scanPAV/scanPAV/scanPAV-bin/
      5. The default aligner is bwa, but you also have the chance to use smalt, which is faster;
      6. The final results from smalt and bwa are relatively consistant, but there are some small differences

(3) Install

gunzip scanPAV-v1.1.tar.gz 
tar xvf scanPAV-v1.1.tar
make 

Please contact Zemin Ning ( zn1@sanger.ac.uk ) or Francesca Giordano ( fg6@sanger.ac.uk ) for any further information. 
 



a
