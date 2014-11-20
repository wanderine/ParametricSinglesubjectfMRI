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

	# Loop over all subjects
	for i in /home/andek/Data/fcon1000/${Study}/*; do 

   		# Check if fMRI data exists for this directory
	    if [ -e ${i}/func/rest.nii.gz ]; then
				
   		    # Go to current directory
	        cd $i
	        # Get subject name
	        Subject=${PWD##*/}
			echo "-------------------------------"	
	   		echo "Processing" $Subject " " $Subjects
			echo "-------------------------------"	
	        # Go back to original directory
	        cd $start_directory

			rm $start_directory/Results/${Study}/${Smoothing}/${Design}/${Subject}.results/activity_map.nii.gz

			3dcalc -a "Results/${Study}/${Smoothing}/${Design}/${Subject}.results/stats.${Subject}_REML+tlrc[2]" -b "Results/${Study}/${Smoothing}/${Design}/${Subject}.results/full_mask.${Subject}+tlrc" -expr 'step(b)*a' -prefix $start_directory/Results/${Study}/${Smoothing}/${Design}/${Subject}.results/activity_map.nii.gz

	    else
			echo "This directory does not contain any fMRI data!"
	    fi
	done
done



