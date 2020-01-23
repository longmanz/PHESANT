#!/bin/bash

# Boot disk - change default disk size from 10GB to 200GB. (Debian GNU/Linux - though this doesn't matter too much).

sudo apt-get update
HOME=/home/dpalmer
INSTALL_DIR=/home/dpalmer/R/

# Install a bunch of packages to get R up and running on the cluster.
PKGS="bzip2 build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libpcre3-dev gfortran openjdk-8-jdk libxml2-dev"
for p in $PKGS; do
    sudo apt-get install -y $p
done

# Install R 3.4.
wget https://cran.r-project.org/src/base/R-3/R-3.4.1.tar.gz
tar -xvzf R-3.4.1.tar.gz
rm R-3.4.1.tar.gz

cd R-3.4.1
./configure --with-readline=no --with-x=no --prefix=$INSTALL_DIR
make
make install

# Add R to the path.
export PATH=$INSTALL_DIR/bin:$PATH

## Default repo
# local({r <- getOption("repos")
#        r["CRAN"] <- "http://cran.r-project.org" 
#        options(repos=r)
# })

# Install a bunch of packages:
Rscript -e 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)}); install.packages(c("data.table", "bit64", "optparse", "hyperSpec", "dplyr"))'

# Install git.
sudo apt-get install git
cd $HOME
# clone PHESANT library
git clone https://github.com/astheeggeggs/PHESANT.git

# I've assumed that the correctly parsed pheotype file to be passed to PHESANT is:
# /home/eatkinso/neale_lab_parsed_QC_Oct2019.tsv 

cd PHESANT/WAS

NUMPARTS=10
for i in `seq 1 $NUMPARTS`;
do
Rscript phenomeScan.r \
		--phenofile="../../../eatkinso/neale_lab_parsed_QC_Oct2019.tsv" \
		--variablelistfile="../variable-info/outcome_info_final_multi_ancestry_jan2020.tsv" \
		--datacodingfile="../variable-info/data-coding-ordinal-info-nov2019-update.txt" \
		--userId="userId" \
		--resDir="../../" \
		--out="multi_ancestry_jan_2020" \
		--partIdx="$i" \
		--numParts="${NUMPARTS}"
done

for i in `seq 1$NUMPARTS`;
do
	gzip -k "multi_ancestry_jan_2020.$i.tsv" 
done
