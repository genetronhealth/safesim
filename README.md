The SafeSim package currently contains two tools: SafeMut and SafeMix.
The tools SafeMut and SafeMix are all aware of the unique molecular identifiers (UMIs) which are also known as molecular barcodes. 
The tool SafeMut spikes variants specified in a VCF file into an existing BAM file (which is not supposed to contain any of the specified variants) by directly modifying read sequences in a UMI-aware manner. 
The tool SafeMix mixes two BAM files in a predefined way to simulate variants at predefined fractions in a UMI-aware manner. 

# How to install

This package requires a compiler that supports the C++14 standard.
The Makefile in this directory compiles with g++, but the Makefile can be easily modified to use another compiler instead of g++ (for example, clang).
To install from scratch, please run: (./install-dependencies.sh && make clean && make all -j4 && make deploy). 
Please note that ./install-dependencies.sh requires bzip2 to decompress the downloaded files with the (.tar.bz2) extension.
This package depends on git 2.12+, htslib 1.6+ and bcftools 1.6+ (lower versions of htslib and bcftools may also work, but are not tested).
If these two dependencies were already installed, then install-dependencies.sh may not be necessary.
For trouble-shooting with the installation of htslib and bcftools, please check their official repositories at https://github.com/samtools/htslib and https://github.com/samtools/bcftools.
More specifically, if any error message containing "error while loading shared libraries" pops up, please use the command (./configure --disable-plugins --disable-bz2 --disable-lzma --disable-libcurl --disable-s3 --disable-largefile) to build the corresponding htslib first.
Although not required, it is highly recommmended that bcftools is installed at a system location (a location that can be found in the PATH environment variable).

In total, the installation should take about 4 minutes.

# How to use

Use the bin/safemut -h and bin/safemix -h commands to see the tool-specific usage help.

For UMI (unique molecular identifier, a.k.a. molecular barcode) to be detected, the read name (QNAME) in the input BAM file should be in the format of originalName#UMI.
For example, the UMI-labeled read name can be
 1. "H5G5ABBCC:4:1209:10114:63736#ACGTAACCA" (ACGTAACCA is the single-strand barcode) or 
 2. "H5G5ABBCC:1:3010:10412:33669#AGTA+TGGT" (AGTA+TGGT is the duplex barcode).
Please note that INFO/FA must be defined the header of the input VCF file in order to be effective, otherwise the default value of allele fraction is used by the simulation. 

# Other things

The word "safe" refers to the Safe-Sequencing System (Safe-SeqS) first described at https://doi.org/10.1073/pnas.1105422108 

For more information, please check the wiki.

