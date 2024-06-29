#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<unistd.h>
#include<sys/shm.h>

//Generating simulated time series----------------------------------------------

int main(){

	float P = 500.12345; //in ms
	float pulse_width = 50.54321; //in ms
	float pulse_height= 1.0;	
	int N = 1000000;// #total samples in data
	float tau = 0.234;// #ms

	int num= (int)ceil((N* tau)/P);/// #no. of pulses/cycles. Last cycle might be incomplete
	int w= (int)round(pulse_width/tau);// #no. of samples in pulse width

	int index,p=1000; //#position of pulse begining (in terms of no. of samples) in every cycle
	
//create shared memory and fill time series-----------------------------------	
	// ftok to generate unique key
	key_t key = ftok("shm_folding",65);
	
	// shmget returns an identifier in shmid
	int shmid = shmget(key,sizeof(float)*N,0666|IPC_CREAT);
	
	// shmat to attach to shared memory
	float *time_series=(float*) shmat(shmid,NULL,0);
	
	float *array=time_series;
	for(int i=0;i<N;i++) *array++=0.2;
	
	for(int i=0;i<num;i++){
		index= (int)round(P*i/tau); //#index of every cycle
		for(int j=0;j<w;j++) time_series[index+p+j]=1.0;
	}
	
	shmdt(time_series);
	
	return 0;

}

