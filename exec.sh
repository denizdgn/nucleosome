#!/bin/bash
shopt -s expand_aliases
path_to_scripts=$(pwd)
chmod +x nucleus

echo '##### Required for nucleus #####' >> ~/.bashrc
echo 'conda activate nucleosome' >> ~/.bashrc
echo 'export LC_ALL=en_US.utf-8' >> ~/.bashrc
echo 'export LANG=en_US.utf-8' >> ~/.bashrc
echo "alias nucleus='python "''$path_to_scripts/nucleus''"'" >> ~/.bashrc
echo '##### Required for nucleus #####' >> ~/.bashrc


conda activate nucleosome