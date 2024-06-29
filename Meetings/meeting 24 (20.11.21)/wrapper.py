
import subprocess
import time
import shm
import struct

#rewrites gptool.in with changed parameters. Is called everytime scan is turned on
def rewrite_gptooldotin():

	#taking relevant parameters from the shm header
	a=subprocess.Popen(["/home/gpuuser/GWB/release/bin/read_shm_hdr","-c","TST2459"], stdout=subprocess.PIPE)

	(out, err) = a.communicate()
	out=out.decode("utf-8").split("\n")

	source=out[0].split()[-1]
	beam_type=out[1].split()[-1]
	gsb_stokes=int(out[2].split()[-1])
	bandwidth=float(out[3].split()[-1])
	channels=int(out[4].split()[-1])
	side_band=int(out[5].split()[-1])
	
	#gptool demands frequency given to be the lowest in the range
	rf=float(out[6].split()[-1])
	if(side_band==-1): rf=rf-bandwidth
	
	sampling_interval=float(out[7].split()[-1])/1000.0

	#editing existing gptool.in
	f= open("gptool.in",'r')
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

	params="\n".join(params)

	f=open("gptool.in",'w')
	f.write(params)
	f.close()

	return

#the key to the shm of das header. Taken from gmrt_newcorr.h
DAS_H_KEY = 1030

#in gptool when scan is on, the 'status' variable of the das header is set equal to DAS_START which simply holds the value 2 as can be seen in protocol.h. Same goes for DAS_STOP
DAS_START = 2
DAS_STOP = 3

#attaching to the das header shm
shmid=shm.getshmid(DAS_H_KEY)
das_hdr = shm.memory(shmid)
das_hdr.attach()

#these variables take care of the state of the observation
new_source=False  #need not be a different source, can be the same source
prev_status=curr_status=DAS_STOP 

while True:
	while True: 
		#reading the shm header status continuously
		a=das_hdr.read(4,4)
		a=struct.unpack('i',a)

		if(a[0]!= DAS_START):
			#new_source=False
			prev_status=curr_status
			curr_status=DAS_STOP
			print("DAS not in start mode\n")
			#p.terminate()

			if (prev_status==DAS_START and curr_status==DAS_STOP): p.terminate()
			time.sleep(1)
			break

		prev_status=curr_status
		curr_status=DAS_START
		
		#time.sleep(1)

		if (prev_status==DAS_STOP and curr_status==DAS_START): new_source=True

		if new_source==True: 
			rewrite_gptooldotin()
			p= subprocess.Popen(["/common-h10/gpuuser/pulsar/packages/gptool_ver4.5.1/gptool", "-r","-shmID", "1"])
		
		
		new_source=False
		
	


	




