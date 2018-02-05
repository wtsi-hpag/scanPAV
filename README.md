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
		$ ./install.sh
		
If everything compiled saccessfully you must see the final comment: 
		"Congrats: installation successful!"		

(Tested with gcc-4.9.2) 

#### External packages
The genome aligner BWA (http://bio-bwa.sourceforge.net) and SMALT (http://www.sanger.ac.uk/science/tools/smalt-0) are downloaded and compiled by scanPAV.

#### Run:

           $ /full/path/to/scanPAV/src/scanPAV -nodes <nodes> -align <aligner> -score <sw-score> \
	   	      </full/path/to/assembly_1.fasta> </full/path/to/assembly_2.fasta> \ 
		      <pavs_present_in_assembly_1.fasta>
           
           where:
	          /full/path/to/assembly_1.fasta: full path to the assembly file to be considered as "presence-assembly"
	     	  /full/path/to/assembly_2.fasta:  full path to the assembly file to be considered as "absence-assembly"
	     	  pavs_present_in_assembly_1.fasta:   output name for the pav sequences. 
	     		These are the PAVs present in assembly_1.fasta, but absent in assembly_2.fasta
	     
	       parameters:
             nodes:    number of CPUs requested  [ default = 30 ]
             sw-score: smith-waterman alignment score [ default = 550 ]
             aligner:  sequence aligner: bwa or smalt [ default = bwa ]
             
Please notice that scanPAV technically only finds sequences present in assembly_1.fasta that are 'absent' in assembly_2.fasta. To find the sequences that are absent in assembly_1.fasta but present in assembly_2.fasta, scanPAV
has to be re-run inverting assembly_1.fasta and assembly_2.fasta:

	   $ /full/path/to/scanPAV/src/scanPAV -nodes <nodes> -align <aligner> -score <sw-score> \
	   	      </full/path/to/assembly_2.fasta> </full/path/to/assembly_1.fasta> \ 
		      <pavs_absent_in_assembly_1.fasta> 
	   
	   	where: 	 
		   pavs_absent_in_assembly_1.fasta:  output name for the pav sequences. 
			These are the PAVs absent in assembly_1.fasta, but present in assembly_2.fasta

	
#### Results
The PAV sequences will be in the file pavs_present_in_assembly_1.fasta in your working directory. If you also run the scanPAV pipeline with assembly_1.fasta and assembly_2.fasta in reverse order, then you'll find also the file  pavs_absent_in_assembly_1.fasta in your working directory.

#### Some Notes on the aligners:
1. You can check if the bwa installation was succesfull by launching:
         
                    $ /full/path/to/scanPAV/src/scanPAV-bin/bwa
		    
   This should print out the bwa help information. If this gives you error messages, 
      you can try to re-install it (by running again the install.sh script) or you
      can link your own bwa installation into /full/path/to/scanPAV/src/scanPAV-bin/bwa :
	      
                    $ ln -sf /full/path/to/my/own/bwa  /full/path/to/scanPAV/src/scanPAV-bin/bwa
              
2. The default aligner is bwa, but you also have the chance to use smalt, which is faster;
3. Results from smalt and bwa are relatively consistent, but some small differences 
      are to be expected due to different aligner sensitivities and internal parameters.
 
