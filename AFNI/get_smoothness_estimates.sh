#!/bin/bash

#Study=Cambridge
Study=Beijing

Design=boxcar10_REML

start_directory=/home/andek/Research_projects/SingleSubject

#for Smooth in 1 2 3 4
for Smooth in 2
do

	if [ "$Smooth" -eq "1" ]; then
		Smoothing=4mm
	elif [ "$Smooth" -eq "2" ]; then
		Smoothing=6mm
	elif [ "$Smooth" -eq "3" ]; then
		Smoothing=8mm
	elif [ "$Smooth" -eq "4" ]; then
		Smoothing=10mm
	elif [ "$Smooth" -eq "5" ]; then
		Smoothing=12mm
	elif [ "$Smooth" -eq "6" ]; then
		Smoothing=14mm
	elif [ "$Smooth" -eq "7" ]; then
		Smoothing=16mm
	fi

	Subjects=0.0
	Significant=0.0
	FWE=0.0
	one=1.0

	touch Results/smoothness_estimates_${Study}_Smoothing${Smoothing}.txt

	# Loop over all subjects
	for i in /home/andek/Data/fcon1000/${Study}/*; do 

   		# Check if fMRI data exists for this directory
	    if [ -e ${i}/func/rest.nii.gz ]; then
	
			Subjects=$(echo "scale=3;$Subjects + $one" | bc)

   		    # Go to current directory
	        cd $i
	        # Get subject name
	        Subject=${PWD##*/}
			echo "-------------------------------"	
	   		echo "Processing" $Subject " " $Subjects
			echo "-------------------------------"	
	        # Go back to original directory
	        cd $start_directory

			# Check for significant clusters, calculate cluster size required
			cat Results/${Study}/${Smoothing}/${Design}/${Subject}.results/blur_est.${Subject}.1D
			fxyz=`cat Results/${Study}/${Smoothing}/${Design}/${Subject}.results/blur_est.${Subject}.1D`
			temp=()
			temp2=${fxyz[$((0))]}
			temp+=($temp2)
			xsmoothness=${temp[$((14))]}
			ysmoothness=${temp[$((15))]}
			zsmoothness=${temp[$((16))]}
			echo "X smoothness is $xsmoothness" 
			echo "Y smoothness is $ysmoothness" 
			echo "Z smoothness is $zsmoothness" 

			echo "$xsmoothness" >> Results/smoothness_estimates_${Study}_Smoothing${Smoothing}.txt

	    else
			echo "This directory does not contain any fMRI data!"
	    fi
	done
done



