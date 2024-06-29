
from __future__ import print_function
import subprocess
import time
import shm  #to attach to system V shm used by gptool. This was specifically chosen to work in python 2.7
import struct
import signal
import re
import os
from datetime import datetime
import sys

#--------------------------------------------------------------
results_dirpath= "/common-h10/gpuuser/gbmon/"
get_gtaccode_command= "ssh gwbh6 /home/gpuuser/GWB/release/bin/get_gtac_code"
gptool_envfile= "/common-h10/gpuuser/pulsar/packages/scripts/presto_old_yg.bash"
gptool_location= "/common-h10/gpuuser/pulsar/packages/gptool_ver4.3.6_gbmon/"
read_shmhdr_command= "/home/gpuuser/GWB/release/bin/read_shm_hdr1"
get_beamno_command= "uname -a"
gptooldotin_master="/common-h10/gpuuser/pulsar/packages/gbmon_tests/gptool.in_master"
png_viewer_command= "/common-h10/gpuuser/pulsar/packages/gbmon_tests/gbmon_png_viewer"

#the key to the shm of das header. Taken from gmrt_newcorr.h in gptool_ver4.3.6
DAS_H_KEY = 1030

#in gptool when scan is on, the 'status' variable of the das header is set equal to DAS_START which simply holds the value 2 as can be seen in protocol.h in gptool_ver4.3.6. Same goes for DAS_STOP
DAS_START = 2
DAS_STOP = 3

#--------------------------------------------------------------

os.system("source "+gptool_envfile)

#If run from GUI, gtac code will be provided as an argument
gtac_code=""
if len(sys.argv)>1:
	gtac_code= sys.argv[1]
print(gtac_code)

#defining a custom print function, to simultaneously print on terminal as well as write to the log file of gbmon
class Logger:
 
    def __init__(self, filename):
        self.console = sys.stdout
        self.file = open(filename, 'a')
 
    def write(self, message):
        self.console.write(message)
        self.file.write(message)
 
    def flush(self):
        self.console.flush()
        self.file.flush()

#taylor the date and time into a particular format for gbmon
def formatDateTime(date_time):
	dateTime=str(date_time).split()
	date=dateTime[0].split('-')
	date.reverse()
	date=date[0]+date_time.strftime("%b")+date[2]
	date= date.lower()
	timestamp=dateTime[1].split('.')[0]

	return date, timestamp

#read the time from beam host shm header
def get_hdrtime():
	global gtac_code

	x=[read_shmhdr_command,"-c",gtac_code]
	a=subprocess.Popen(x, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	(out, err) = a.communicate()
	out=out.decode("utf-8").split("\n")

	timestamp=out[1].split()[-1]
	timestamp=timestamp[0:len(timestamp)-13]
	timestamp=timestamp.split(':')	
	timestamp="".join(timestamp)	

	return timestamp

#get the beam number using command 'uname -a'. it gives name of host machine like gwbh7 which we map to beam 1 and so on.
def get_beam_no():
	global beam_no

	os.system(get_beamno_command+" | awk '{print $2}' > tmp")
	host_no=open('tmp', 'r').read().split('.')[0][4:]  #.split('.')[0]
	os.remove('tmp')

	beam_no = int(host_no) - 6
	return beam_no

#gets observation parameters from shm header, decides which gptool.in to read and then rewrite gptool.in with changed parameters. Its called everytime gptool is run--------
def rewrite_gptooldotin():
	global scan_no, time_start, date, gtac_code, results_subdirpath, beam_no, results_dirpath
	
	print("Getting observation parameters and preparing gptool.in\n")

	#taking relevant parameters from the shm header-------------------
	x=[read_shmhdr_command,"-c",gtac_code]
	a=subprocess.Popen(x, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	(out, err) = a.communicate()
	out=out.decode("utf-8").split("\n")

	x=" ".join(x)
	gbmon_log.file.write(x+"\n")
	for line in out:
		gbmon_log.file.write(line+"\n")

	scan_no=int(out[0].split()[-1])
	time_start=get_hdrtime();
	date= formatDateTime(datetime.now())[0] #out[2].split()[-1]

	source= out[3].split()[-1] #"B0329+54"
	side_band=int(out[4].split()[-1])

	bandwidth=float(out[8].split()[-1])

	#gptool requires frequency given to be the lowest in the range in gptool.in-------
	rf=float(out[5].split()[-1])
	if(side_band==-1): rf=rf-bandwidth

	beam_type=out[6].split()[-1]
	gsb_stokes=int(out[7].split()[-1])
	
	channels=int(out[9].split()[-1])
	sampling_interval=float(out[10].split()[-1])/1000.0
	

	beam_no= get_beam_no()
	print("BEAM NO. ", beam_no)
	
	results_subdirpath= results_dirpath+"{}_{}_beam{}".format(date, gtac_code, beam_no)
	os.system("cd {} && mv gbmon.log_{}_{} gbmon.log_{}_{}_beam{}".format(results_dirpath, date ,gbmon_timestart ,date ,gbmon_timestart ,beam_no))	

	logfile="{}gbmon.log_{}_{}_beam{}".format(results_dirpath, date, gbmon_timestart, beam_no)
	gbmon_log=Logger(logfile)
	sys.stdout=gbmon_log

	#this if-else decides which gptool.in to use. If a gptool.in already exists in a time-stamped subdir in results_dirpath, (which means GUI is being used) use that. If not, then use gptool.in_beam<number> if it exists. If it doesnt then resort to gptool.in_master  
	if (os.path.isfile(results_subdirpath+"/gptool.in")): 
		string="mv {} {}{}_{}_{}_beam{}/".format(results_subdirpath, results_dirpath,  date,time_start, gtac_code,beam_no)
		os.system(string)

		results_subdirpath="{}{}_{}_{}_beam{}/".format(results_dirpath,  date,time_start, gtac_code,beam_no)

		f= open("{}gptool.in".format(results_subdirpath),'r')

	else:
		results_subdirpath="{}{}_{}_{}_beam{}/".format(results_dirpath,  date,time_start, gtac_code,beam_no)
		os.system("mkdir -p {}".format(results_subdirpath))

		if (os.path.isfile("{}gptool.in_beam{}".format(results_dirpath, beam_no))):
			f= open("{}gptool.in_beam{}".format(results_dirpath, beam_no),'r')

		else: 
			f= open(gptooldotin_master,'r')

	#editing existing gptool.in----------------------------
	params=f.read().split("\n")
	f.close()

	params[15]=source+"\t: Pulsar name"
	params[3]=beam_type+"\t\t: Beam mode"
	if(gsb_stokes==2): params[4]="0\t\t: Polarization mode (0-> intesity, 1-> stokes data)"
	elif(gsb_stokes==4): params[4]="1\t\t: Polarization mode (0-> intesity, 1-> stokes data)"
	params[9]=str(bandwidth)+"\t\t: Bandwidth(in Mhz)"
	params[11]=str(channels)+"\t\t: Number of channels"

	params[10]=str(side_band)+"\t\t: Sideband flag (-1-> decreasing +1-> increasing)"
	params[8]=str(rf)+"\t\t: Frequency band (lowest value in Mhz)"
	params[12]=str(sampling_interval)+"\t\t: Sampling Interval (in ms)"

	flag_channels=int(.03*channels)
	params[33]=str(flag_channels)+"\t\t: Number of channels to flag at band beginning"
	params[34]=str(flag_channels)+"\t\t: Number of channels to flag at band end"	

	smoothing_window=int(10*channels/bandwidth)	
	params[46]=str(smoothing_window)+"\t\t: Smoothing window size for bandshape normalization (in number of channels)"

	#if source name isnt a pulsar name, and if period and dm are '-1' in gptool.in, fix the values of DM and period----------
	period=params[16][:2]
	dm=params[17][:2]

	pulsar_regex= re.compile(r'^(B|J)(\d{4})(\+|\-)(\d{2}|\d{4})$')	
	mo=pulsar_regex.search(source)
	if mo==None:
		if period=="-1": period=100
		if dm=="-1": dm=0

		print("Source name is not a standard Pulsar name. \nUsing Period: {}ms and DM: {}".format(period,dm)) 
		
		params[16]=str(period)+"\t\t: Pulsar period (in milliseconds)"
		params[17]=str(dm)+"\t\t: DM (in pc/cc)"

	else:
		print("Source name is a standard Pulsar name. \n")
		if period=="-1": print("Using tempo-obtained period\n")
		else: print("Using Period: {} ms\n".format(period))

		if dm=="-1": print("Using psrcat-obtained DM\n")
		else: print("Using DM: {} pc/cc\n".format(dm))


	#prepare the required directories, if pgplot displau device is /png
	display_device= params[63][:4]	
	if (display_device== "/png" or display_device== "/PNG"): 
		os.system("mkdir -p {}realtime_beam{}/plots/".format(results_dirpath, beam_no))
		os.system("if test -z \"$(ls -A {0}realtime_beam{1}/plots/)\" ; then : ; else rm {0}realtime_beam{1}/plots/* ; fi ".format(results_dirpath, beam_no))

		params[63]= "{}realtime_beam{}/plots/png\t\t:pgplot display device type".format(results_dirpath, beam_no)
		os.system("mkdir -p {}/plots/".format(results_subdirpath))	

		#png_viewer=subprocess.Popen(png_viewer_command.split(" ")+[results_subdirpath+"/plots"])		
	
	params="\n".join(params)

	f=open("{}/gptool.in".format(results_subdirpath),'w')
	f.write(params)
	f.close()

#if no gtac code provided as command-line arg, get it by the get_gtaccode_command
def get_gtac_code():
	
	global gtac_code
	if gtac_code=="": 
		
		print("\nGTAC code not provided\nGetting GTAC code/Observation code\n")
		
		x=get_gtaccode_command.split(" ")
		
		a=subprocess.Popen(x, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		
		(out, err) = a.communicate()

		x=" ".join(x)
		print(x+"\n")
		
		for i in range(len(out)-2):
			if out[i:i+3]=="PRJ":
				gtac_code=out[i+6:-1]
				break

	
	print("GTAC code/Observation code being used: {}\n".format(gtac_code))	
	gbmon_log.file.flush()	

	return gtac_code


#handles Ctrl+C interruption and further options
def keyboardInterruptHandler(signal, frame):

	global flag1, flag2, date, time_start, time_stop, gbmon_timestart, gbmon_timestop;
	gbmon_log.file.seek(0,2)

	print("\n\nKeyboardInterrupt\n")	
	
	if gptool==None:
		i=input("\nEnter 0 to Continue\nEnter 1 to Exit GBMon\n")
		print("\nEnter 0 to Continue\nEnter 1 to Exit\nWaiting for input\n")
	
		gbmon_log.file.write(str(i)+"\n")
		if i==1: 
			print("\nExiting the gbmon\n")
			#gbmon_log.file.close()
			gbmon_timestop=formatDateTime(datetime.now())[1][:5]    #get_hdrtime();
			os.system("cd {} && mv gbmon.log_{}_{}_beam{} gbmon.log_{}_{}_{}_beam_no{}".format(results_dirpath,date,gbmon_timestart, beam_no, date, gbmon_timestart, gbmon_timestop, beam_no))
			
			exit(0)
		
		if i==0: 
			flag4=0
			flag2=0
			flag1=0

	else: 
		stop_gptool()
		i=input("\nEnter 0 to restart gptool\nEnter 1 to exit the GBMon\nWaiting for input\n")
		gbmon_log.file.write(str(i)+"\n")	    
	
		if i==1: 
			print("Exiting the gbmon\n")
			gbmon_timestop=formatDateTime(datetime.now())[1][:5] #get_hdrtime();
			os.system("cd {} && mv gbmon.log_{}_{} gbmon.log_{}_{}_{}".format(results_dirpath, date ,gbmon_timestart ,date ,gbmon_timestart ,gbmon_timestop))		
			exit(0)
		elif i==0: 
			start_gptool()
	
	gbmon_log.file.flush()


#move saved files to required directories
def save_files():
	global scan_no, time_start,time_stop, date, results_dirpath, results_subdirpath

	#os.chdir(results_subdirpath)
	#string="for i in bandshape.gpt profile_filtered.gpt profile_unfiltered.gpt gptool.in; do [ -f \"$i\" ] && mv \"$i\" \"$i\"_{}_scan{}_{}_{}; done".format(date,scan_no, time_start,time_stop) 
	#os.system(string)

	os.chdir(results_dirpath)
	string="mv {}_{}_{}_beam{} {}_{}_{}_{}_beam{}/".format(date,time_start,gtac_code, beam_no, date,time_start, time_stop,gtac_code,beam_no)
	os.system(string)

	results_subdirpath= results_dirpath+"{}_{}_{}_{}_beam{}".format(date,time_start,time_stop,gtac_code, beam_no)

	plots_dir = "{}/plots/".format(results_subdirpath)
	string="if [ -d {} ] ; then : ; cp {}realtime_beam{}/plots/* {}".format(plots_dir, results_dirpath, beam_no, results_subdirpath, plots_dir)
	os.system(string)


def start_gptool():
	global gptool, quiet_mode, scan_no, time_start, date, results_subdirpath

	rewrite_gptooldotin()

	os.chdir(results_subdirpath)	
	x=[gptool_location+"gptool" , "-tempo2", "-l", logfile,"-c", str(scan_no), "-r","-shmID", "1"]

	if quiet_mode==True: 
		FNULL = open(os.devnull, 'w')
		print("Starting gptool in quiet mode..\n")
		
		gptool= subprocess.Popen(x,stdout=FNULL, stderr=subprocess.STDOUT)
		x=" ".join(x)
		gbmon_log.file.write(x+"\n")
		gbmon_log.file.write("STDOUT= DEVNULL\n")

	else: 
		print("Starting gptool in verbose mode..\n")	
		gptool= subprocess.Popen(x)
		x=" ".join(x)
		gbmon_log.file.write(x+"\n")

	print("----------------------------------------\nScan No.: {}\nDate: {}\nStart Time: {}\n".format(scan_no, date, time_start))

	gbmon_log.file.flush()

def stop_gptool():
	global gptool, time_stop, png_viewer

	gbmon_log.file.seek(0,2)
	
	print("Stopping gptool and saving files..please wait\n")
	gptool.send_signal(signal.SIGINT)
	time.sleep(6) #wait for gptool to stop completely after Ctrl+C signal

	os.chdir(results_subdirpath)
	os.system('python2 {}plotgptoolsummary.py'.format(gptool_location))

	time_stop=get_hdrtime();
	print("Done\ngptool stop time: {}\n".format(time_stop))	
	save_files()
	gptool=None

	if png_viewer!=None:
		png_viewer.kill()
	png_viewer=None	

	gbmon_log.flush()

#print current time and waiting time while GBMon waits for shm to be created or scan to be started
def waiting_time(waitingstart_time):
	waiting_time= str(datetime.now()- waiting_starttime).split('.')[0]
	
	pos=gbmon_log.file.tell()
	print("Current Time: {}, Waiting Time: {}\r".format(formatDateTime(datetime.now())[1], waiting_time),end="")
	gbmon_log.file.seek(pos)
	
	gbmon_log.flush()
	
#for gbmon's logfile
date, gbmon_timestart= formatDateTime(datetime.now())  
gbmon_timestart=gbmon_timestart[:5]
logfile="{}gbmon.log_{}_{}".format(results_dirpath, date, gbmon_timestart)

#defining custom print function
gbmon_log=Logger(logfile)
sys.stdout=gbmon_log

gptool=None
png_viewer=None

quiet_mode=False
i=input("Enter 0 for quiet mode\nEnter 1 for verbose mode\n")
if i==0: quiet_mode=True
gbmon_log.file.write(str(i)+"\n")
gbmon_log.file.flush()


#these variables take care of the state of the observation
new_scan=False  #need not be a different source, can be the same source
prev_status=curr_status=DAS_STOP 
observation_stop=True
flag1=0 #to print "Observation has not been started.\nCould not find shared memory.\nWaiting.." once
flag2=0 # to print "Data aquisition is not in start mode. Waiting..\n" once
flag3=0 #to print Data aquisition is in start mode\n" once
flag4=0 # to print "\nLooking for ongoing observation\n" once

#looks out for Ctrl+C termination
signal.signal(signal.SIGINT, keyboardInterruptHandler)


while True: 
	time.sleep(1)	

	try: #checking for existence of shared memory continuously

		if (observation_stop==True and flag4==0): 
			print("\nLooking for ongoing observation\n")

		shmid=shm.getshmid(DAS_H_KEY)
		
		if (observation_stop==True): 
			print("Success. Found ongoing observation\n")

		gbmon_log.file.seek(0,2)
		if observation_stop==True:
			das_hdr = shm.memory(shmid)
			das_hdr.attach()
			observation_stop=False
			flag1=0

			gtac_code=get_gtac_code()

	except KeyError:
		if(flag1==0):
			waiting_starttime=datetime.now()
			print("Observation has not been started.\nCould not find shared memory.\nWaiting..")
			flag1=1 

		waiting_time(waiting_starttime)		
		observation_stop=True
		flag4=1
		continue


	#reading the shm header status continuously
	a=das_hdr.read(4,4)
	a=struct.unpack('i',a)
	
	if(a[0]!= DAS_START):
		prev_status=curr_status
		curr_status=DAS_STOP

		if (prev_status==DAS_START and curr_status==DAS_STOP): 
			print("\nData aquisition has been stopped\n")			
			stop_gptool()

		if flag2==0: 
			waiting_starttime=datetime.now()
			print("\nData aquisition is not in start mode. Waiting..\n")
			sys.stdout.flush()
			flag2=1
			flag3=0

		waiting_time(waiting_starttime)
		time.sleep(1)
		continue

	flag2=0
	gbmon_log.file.seek(0,2)
	if flag3==0: 
		print("Data aquisition is in start mode\n")
		flag3=1		

	prev_status=curr_status
	curr_status=DAS_START
	
	if (prev_status==DAS_STOP and curr_status==DAS_START): new_scan=True

	if new_scan==True: 
		gbmon_log.file.flush()
		start_gptool()	
		
	
	
	new_scan=False
		
	


	




