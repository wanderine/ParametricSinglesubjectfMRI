#!/bin/bash

# A help script to get the smoothness estimate of each dataset

clear

Study=Cambridge
DF=111

#Study=Beijing
#DF=217

# Design
Design=boxcar10


# Loop over different smoothing levels
for SmoothingLevel in 2
do
	date1=$(date +"%s")

	if [ "$SmoothingLevel" -eq "1" ] ; then
		Smoothing=4mm
	elif [ "$SmoothingLevel" -eq "2" ] ; then
		Smoothing=6mm
	elif [ "$SmoothingLevel" -eq "3" ] ; then
		Smoothing=8mm
	elif [ "$SmoothingLevel" -eq "4" ] ; then
		Smoothing=10mm
	elif [ "$SmoothingLevel" -eq "5" ] ; then
		Smoothing=12mm
	elif [ "$SmoothingLevel" -eq "6" ] ; then
		Smoothing=14mm
	elif [ "$SmoothingLevel" -eq "7" ] ; then
		Smoothing=16mm
	fi

	results_directory=/home/andek/Research_projects/RandomGroupAnalyses/Results/${Study}/${Smoothing}/${Design}

	Significant=0
	Subjects=0

	touch Results/smoothnessestimates_${Study}.txt

	# Loop over all subjects
	for i in /home/andek/Data/fcon1000/${Study}/* ; do

	    # Check if fMRI data exists for this directory
	    if [ -e ${i}/func/rest.nii.gz ]; then

			# Go to current directory
			cd $i
			# Get subject name
		   	Subject=${PWD##*/}
		    echo "Processing" $Subject
			# Go back to original directory
			cd $design_directory

			# Get smoothness 
			text=`smoothest -d ${DF} -m /home/andek/Research_projects/RandomGroupAnalyses/Results/${Study}/${Smoothing}/${Design}/${Subject}.feat/mask.nii.gz -r /home/andek/Research_projects/RandomGroupAnalyses/Results/${Study}/${Smoothing}/${Design}/${Subject}.feat/stats/res4d.nii.gz -V`

			#smoothest -d ${DF} -m /home/andek/Research_projects/RandomGroupAnalyses/Results/${Study}/${Smoothing}/${Design}/${Subject}.feat/mask.nii.gz -r /home/andek/Research_projects/RandomGroupAnalyses/Results/${Study}/${Smoothing}/${Design}/${Subject}.feat/stats/res4d.nii.gz -V

			temp=${text[$((0))]}
			values=()
			values+=($temp)
			xsmoothness=${values[$((106))]}

			#echo $temp
			echo "$xsmoothness" >> /home/andek/Research_projects/RandomGroupAnalyses/Results/smoothnessestimates_${Study}.txt 
			
	    fi

	done
done

