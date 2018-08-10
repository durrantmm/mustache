#!/usr/bin/env bash

echo "Deactivating any active conda environment..."
source deactivate;

echo "Removing mustache environment if already installed..."
conda env remove -n mustache --yes

echo "Installing mustache environment..."
conda env create -f envs/environment.yml
echo "Entering new mustache environment..."
source activate mustache

echo "Creating config.py..."
CURRENT_DIR="$(pwd -P)"
echo "BLASTDB = '${CURRENT_DIR}/data/mergem_blast_db/mergem'" > mustache/config.py
echo "HMMDB = '${CURRENT_DIR}/data/clusters.faa.hmm'" >> mustache/config.py
echo "FRAGGENESCAN = '${CURRENT_DIR}/bin/FragGeneScan1.30/run_FragGeneScan.pl'" >> mustache/config.py

echo "Unzipping mergem database":
rm -rf data/mergem_blast_db
tar -zxvf data/mergem_blast_db.tar.gz --directory data

echo "Unzipping the transposase cluster HMM file"
rm -rf data/clusters.faa.hmm
gunzip --stdout data/clusters.faa.hmm.gz > data/clusters.faa.hmm

echo "Unzipping and installing the FragGeneScan package..."
rm -rf bin/FragGeneScan1.30
tar -zxvf bin/FragGeneScan1.30.tar.gz --directory bin
cd bin/FragGeneScan1.30
make
cd ../../

echo "Running pip install"
pip install --editable .
source deactivate

echo "Installation Complete."

echo ""
echo "Before running mustache, activate the mustache environment with"
echo "> source activate mustache"
echo ""
echo "You can then run mustache by typing"
echo "> mustache [command]"
echo "Send any questions to mdurrant@stanford.edu"
