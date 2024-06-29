#running gptool on data to give the fullDM_filtered.gpt file

source /home/ygupta/Pulsar/scripts/presto_old_gptool4.3.5_meanval.bash

cd /home/ygupta/work_area/garvit/gptool_results/B1642-03
gptool -f /data/test_data/32_089_18sep17/B1642-03_bm1_cdp_500_200_512_2_18sep17.raw.dat -tempo2

#folding the fullDM_filtered.gpt file

cd /home/ygupta/work_area/garvit/B1642-03/folding

gcc /home/ygupta/work_area/garvit/C_code/folding/folding.c -o /home/ygupta/work_area/garvit/C_code/folding/folding -lm -std=gnu99

./../../C_code/folding/folding /home/ygupta/work_area/garvit/gptool_results/B1642-03/fullDM_filtered.gpt 387.726184613 0.01024 

#using the gptool environment above changes the python paths to something else (using the PYTHONPATH variable and making it an environment variable). removing this variable using unset, restores the python paths that work correctly for the garvit_env python virtual environment
unset PYTHONPATH

#running the stacked_profiles.py code on the fullDM_filtrered.gpt file 

source /home/ygupta/work_area/garvit/python_code/garvit_env/bin/activate

cd /home/ygupta/work_area/garvit/B1642-03/stacked_profiles
python3 /home/ygupta/work_area/garvit/python_code/stacked_profiles/stacked_profiles.py /home/ygupta/work_area/garvit/gptool_results/B1642-03/fullDM_filtered.gpt 387.726184613 0.01024 50 0.0 1.0
