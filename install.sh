#!/bin/bash

conda create --name dash-cam --file spec-file.txt

# Install pbsim for simulator's model
git clone https://github.com/yukiteruono/pbsim2.git

# Install kraken2
git clone https://github.com/DerrickWood/kraken2.git
pushd ./kraken2
./install_kraken2.sh
popd