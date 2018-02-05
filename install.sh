#!/bin/bash


projdir=`pwd`

bindir=$projdir/src/scanPAV-bin/
mkdir -p $bindir
mkdir -p $projdir/src/log/

errs=0

##### Download and install BWA ######

echo "Downloading and installing BWA"
if [[ ! -s $bindir/bwa ]]; then

    if [[ ! -d $projdir/src/bwa ]]; then
	cd $projdir/src/
	git clone https://github.com/lh3/bwa.git &> $projdir/src/log/bwa_cloning.log
    fi

    if [[ ! -s $projdir/src/bwa/bwa ]]; then
	cd $projdir/src/bwa
	make &> $projdir/src/log/bwa_installation.log
    fi

    cp bwa $bindir
fi

if  [[ ! -s $bindir/bwa ]]; then
    echo " !! Error: bwa not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if bwa was downloaded properly:" $projdir/src/log/bwa_cloning.log 
    echo "   Check if the bwa was compiled properly:" $projdir/src/log/bwa_installation.log

    # Cleaning up
    cd $projdir/src
    rm -rf $projdir/src/bwa/bwa $bindir/bwa 
    
    errs=$(($errs+1))
else
    echo " BWA succesfully installed!"
    rm -rf $projdir/src/bwa/
fi


##### Download and install SMALT ######
echo; echo "Downloading and installing Smalt"
if [[ ! -s $bindir/smalt ]]; then
   
    if [[ ! -d $projdir/src/smalt-0.7.4 ]]; then
	cd $projdir/src/
	wget ftp://ftp.sanger.ac.uk/pub/resources/software/smalt/smalt-0.7.4.tgz &> $projdir/src/log/smalt_wget.log
	tar -xvzf smalt-0.7.4.tgz &> $projdir/src/log/smalt_untar.log
	rm -f smalt-0.7.4.tgz
    fi

    cp $projdir/src/smalt-0.7.4/smalt_x86_64 $bindir/smalt
fi
if  [[ ! -s $bindir/smalt ]]; then 
    echo " !! Error: smalt not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if smalt was downloaded properly:" $projdir/src/log/smalt_wget.log 
    echo "   Check if the folder was uncompressed properly:" $projdir/src/log/smalt_untar.log

    # Cleaning up	
    rm -rf $projdir/src/smalt-0.7.4/ $bindir/smalt

    errs=$(($errs+1))
else
    echo " Smalt succesfully installed!"
    rm -rf $projdir/src/smalt-0.7.4/
fi


###### Compile ScanPAV sources ######

echo; echo "Compiling scanPAV sources"

srcs=( scanPAV_seqout scanPAV_fasta scanPAV_shred scanPAV_confirm scanPAV_gapsize scanPAV scanPAV_bwasam scanPAV_bwatag scanPAV_bwagap scanPAV_clean scanPAV_cover scanPAV_process )

cd $projdir/src
make &> $projdir/src/log/sources_compilation.log

echo; echo "Checking installation:"
for src in "${srcs[@]}"; do
    if [[ ! -s $bindir/$src ]]; then 
        echo " !! Error: executable $src missing in $bindir"
	echo "    Please check for errors the log file:" $projdir/src/log/sources_*	
        errs=$(($errs+1))
    fi
done

if [  $errs -gt 0 ]; then echo; echo " ****  Errors occurred! **** "; echo; exit; 
else echo " Congrats: installation successful!"; fi




