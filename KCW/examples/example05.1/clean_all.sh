#!/bin/bash

for dir in nspin1  nspin2  nspin4_noSOC_MAG  nspin4_noSOC_noMAG  nspin4_SOC_MAG  nspin4_SOC_noMAG; do 
cd $dir
 ./clean.sh
cd ../

done 