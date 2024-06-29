
rawdata_paths=("/data/test_data/Band4_test_april2021/B0329+54_bm1_pa_550_200_16_18apr2021.raw" "/data/test_data/Band4_test_april2021/B0740-28_bm1_pa_550_200_16_23apr2021.raw" "/data/test_data/32_089_18sep17/B1642-03_bm1_cdp_500_200_512_2_18sep17.raw.dat" "/data/37_016_data/37_016_band4_15mar20/J2145-0750_bm1_cdp_550_200_1024_2_15mar20.raw0")

#---------------------------------------------------------------------

#this script is to be run in the ~/work_area/garvit/bash_scripts/action_item5_11.09/ directory
PS3="Enter Pulsar number: "

#using the SELECT construct for printing a menu of all pulsars and taking input from user
select pulsar in B0329+54 B0740-28 B1642-03 J2145-0750_msp

do 
	echo -e "\nSelected Pulsar: $pulsar\n"
	
	sed "16 s/pulsar/$pulsar/" template_gptool.in > gptool.in
	
	echo -e "Running gptool for $pulsar data..\n"
	source /home/ygupta/Pulsar/scripts/presto_old_gptool4.3.5_meanval.bash
	gptool -f ${rawdata_paths[$REPLY-1]} -tempo2 

done


