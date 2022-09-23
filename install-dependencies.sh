#!/usr/bin/env sh

# cmdline param 1: option for this script
# cmdline param 2: option for the configure script in bcftools

set -evx
currdir="${PWD}"
mkdir -p "${currdir}/bin/"

if [ $(echo "${1}" | grep skip-htslib | wc -l) -eq 0 ]; then
    mkdir -p "${currdir}/ext/"
    cd "${currdir}/ext/"
    if [ $(echo "${1}" | grep skip-downloading-htslib | wc -l) -eq 0 ]; then
        wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2
    fi
    tar -xvf htslib-1.11.tar.bz2
    mv "${currdir}/ext/htslib-1.11" "${currdir}/ext/htslib-1.11-lowdep"
    cd "${currdir}/ext/htslib-1.11-lowdep"
    ./configure -disable-plugins --disable-libcurl --disable-s3 --disable-largefile ${2} # --disable-bz2 and --disable-lzma are both for disabling CRAM files
    make -j 4
fi

if [ $(echo "${1}" | grep skip-fastq-tools | wc -l) -eq 0 ]; then
    mkdir -p "${currdir}/ext/"
    cd "${currdir}/ext/"
    if [ $(echo "${1}" | grep skip-downloading-fastq-tools | wc -l) -eq 0 ]; then
        wget https://github.com/dcjones/fastq-tools/archive/refs/tags/v0.8.3.tar.gz
    fi
    tar -xvf v0.8.3.tar.gz
    cd "${currdir}/ext/fastq-tools-0.8.3"
    ./autogen.sh
    ./configure
    make -j 4
    cp src/fastq-sort "${currdir}/bin/"
fi

