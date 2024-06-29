
#various metadata about the pulsars in order to run all these codes on them--------------
gptoolresults_path="/home/ygupta/work_area/garvit/gptool_results/"

rawdata_paths=("/data/test_data/Band4_test_april2021/B0329+54_bm1_pa_550_200_16_18apr2021.raw" "/data/test_data/Band4_test_april2021/B0740-28_bm1_pa_550_200_16_23apr2021.raw" "/data/test_data/32_089_18sep17/B1642-03_bm1_cdp_500_200_512_2_18sep17.raw.dat" "/data/37_016_data/37_016_band4_15mar20/J2145-0750_bm1_cdp_550_200_1024_2_15mar20.raw0")


pulsar_dir="/home/ygupta/work_area/garvit/"

pulsar_periods=(714.557690577 166.787764489 387.726184613 16.04999682533)

pulsar_taus=(0.16384 0.16384 0.01024 0.02048)

folding_offsets=(0.0 0.0 0.0 0.5)

pulses_phaseranges_min=(0.0 0.0 0.0 0.0)
pulses_phaseranges_max=(1.0 1.0 1.0 1.0)

stackplot_offsets=(1.0 0.1 0.4 0.05)
#---------------------------------------------------------------------

PS3="Enter Pulsar number: "

#using the SELECT construct for printing a menu of all pulsars and taking input from user
select pulsar in B0329+54 B0740-28 B1642-03 J2145-0750_msp

do 
	echo -e "\nSelected Pulsar: $pulsar\n"
	
	#running gptool on data to give the fullDM_filtered.gpt file
	echo -e "Running gptool for $pulsar data..\n"
	source /home/ygupta/Pulsar/scripts/presto_old_gptool4.3.5_meanval.bash
	cd "${gptoolresults_path}${pulsar}" 
	gptool -f ${rawdata_paths[$REPLY-1]} -tempo2
	
	#running folding.c code on the fullDM_filtered.gpt file
	echo -e "Folding $pulsar filtered dedispersed time series..\n"
	cd "${pulsar_dir}${pulsar}/folding"	
	gcc /home/ygupta/work_area/garvit/C_code/folding/folding.c -o /home/ygupta/work_area/garvit/C_code/folding/folding -lm -std=gnu99
	./../../C_code/folding/folding ${gptoolresults_path}${pulsar}"/fullDM_filtered.gpt" ${pulsar_periods[$REPLY-1]} ${pulsar_taus[$REPLY-1]} ${folding_offsets[$REPLY-1]}

	#using the gptool environment above changes the python paths to something else (using the PYTHONPATH variable and making it an environment variable). removing this variable using unset, restores the python paths that work correctly for the garvit_env python virtual environment
	unset PYTHONPATH

	#running the stacked_profiles.py code on the fullDM_filtrered.gpt file 
	echo -e "Running stacked_profiles.py code for $pulsar..\n"	
	source /home/ygupta/work_area/garvit/python_code/garvit_env/bin/activate
	cd "${pulsar_dir}${pulsar}/stacked_profiles"
	python3 /home/ygupta/work_area/garvit/python_code/stacked_profiles/stacked_profiles.py ${gptoolresults_path}${pulsar}"/fullDM_filtered.gpt" ${pulsar_periods[$REPLY-1]} ${pulsar_taus[$REPLY-1]} 50 ${pulses_phaseranges_min[$REPLY-1]} ${pulses_phaseranges_max[$REPLY-1]} ${stackplot_offsets[$REPLY-1]} ${folding_offsets[$REPLY-1]}  

done



