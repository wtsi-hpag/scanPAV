# scanPAV
Pipeline to detect presence/absence variations (PAVs) when comparing two genome assemblies, assembly_1 and assembly_2.

To find PAV sequences present in in assembly_1 (presence-assembly) but absent in assembly_2 (absence-assembly), the pipeline performs the following steps:

           1. the presence assembly assembly_1 contigs/scaffolds are shred 
                      in chunks of 1000 bases after removal of Ns
           2. the 1000 bases chunks are mapped against the absence assembly assembly_2
           3. groups of consequent 1000 bases chunks that are not find 
                      in the absence assembly assembly_2 are merged and printed out as 'absence PAVs'



### Download and Compile:
Requirements for compiling: gcc

	$ git clone https://github.com/SangerHpag/scanPAV.git
	$ cd scanPAV 
	$ make

(Tested with gcc-4.9.2, XXX) 

#### External packages
The genome aligner BWA (http://bio-bwa.sourceforge.net) and SMALT (http://www.sanger.ac.uk/science/tools/smalt-0) are included in the package.

#### Run:

           $ /full/path/to/scanPAV -nodes <nodes> -align <aligner> -score <sw-score> <presence.fasta> <absence.fasta>
           
           where:
             nodes:    number of CPUs requested  [ default = 30 ]
             sw-score: smith-waterman alignment score [ default = 550 ]
             aligner:  sequence aligner: bwa or smalt [ default = bwa ]
                              

#### Results
The PAV sequences will be in the file presence_absence.fasta in your working directory.

#### Some Notes:
           1. you need to give the full path of the scanPAV input files or the files are in your working directory;
           2. presence_assembly.fasta and absence_assembly.fasta should be in your working directory;
           3. presence_assembly.fasta and absence_assembly.fasta are in full path.
           4. If you use bwa, you need to check if the binary versin provided in the package works:
          
                    $ /full/path/to/scanPAV/scanPAV-bin/bwa
              This should print out the bwa help information. If this gives you error messages, 
              please link your own bwa installation to /full/path/to/scanPAV/scanPAV-bin/bwa :
                    $ ln -sf /full/path/to/my/own/bwa  /full/path/to/scanPAV/scanPAV-bin/bwa
              
           5. The default aligner is bwa, but you also have the chance to use smalt, which is faster;
           6. Results from smalt and bwa are relatively consistent, but some small differences 
              are to be expected due to different aligner sensitivities and internal parameters.
 
