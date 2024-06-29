#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<unistd.h>
#include<sys/shm.h>
#include<stdbool.h>

//Generating simulated time series----------------------------------------------

int main(){

	float P = 500.12345; //in ms
	float pulse_width = 80.54321; //in ms
	float pulse_height= 1.0;	
	int N = 6000000;// #total samples in data
	float tau = 0.234;// #ms
	
	int samples=(int) round(P/tau); //samples in one period

	int num_pulses= (int)ceil((N* tau)/P);/// #no. of pulses/cycles. Last cycle might be incomplete
	int w= (int)round(pulse_width/tau);// #no. of samples in pulse width

	int index; //index of first sample of every cycle
	int p=1000; //#position of pulse begining (in terms of no. of samples) in every cycle
	
	int rewrite_cycles=0; //number of times shm has been written
	bool rewrite= false; //this is changed to True whenever the writer runs out of memory and need to rewrite
	
//create shared memory and fill time series-----------------------------------	
	// ftok to generate unique key
	key_t key = ftok("shm_folding",65);
	
	// shmget returns an identifier in shmid
	int mem_size=1000000; //in number of floats not in bytes
	int shmid = shmget(key,sizeof(float)*mem_size,0666|IPC_CREAT);

	// shmat to attach to shared memory
	float *time_series=(float*) shmat(shmid,NULL,0);
	
	for(int i=0;i<=mem_size;i++) time_series[i]=0.2;
	
	for(int i=0;i<num_pulses;i++){
	
		if(rewrite) {  //resets full array to 0.2
			for(int j=0;j<mem_size;j++) time_series[j]=0.2; //*time_series is set to 0.2 which stops folding until fully shm is written
			printf("SHM Rewritten %d times\n",rewrite_cycles);
		}
	
		index= ((int)round(P*i/tau))- rewrite_cycles*mem_size; //#index of every cycle
		
		int x=mem_size-index; 
		
		//below if-else code is for taking care of partial cycles being written at the end of shm and rest of the cycle to be written in the beginning of the rewriting
		if(x<samples){
			rewrite_cycles++;
			rewrite=true; 
			
			if(x>p && x<p+w){
				for(int j=0;j<x-p;j++) time_series[index+p+j]=1.0;
			}
			
			else if(x>p+w){
				for(int j=0;j<w;j++) time_series[index+p+j]=1.0;	
			}
			
			printf("memory fully written\n");
			
			*time_series=0.0; //this is the 'go' signal for folding
			sleep(5);
		}
		
		else if (index<samples){
		
			rewrite=false;
		
			if(index> samples-p-w && index<samples-p){
				for(int j=0;j<index-samples+p+w;j++) time_series[j]=1.0;
			}
			
			else if(index> samples-p){
				for(int j=0;j<w;j++) time_series[index-samples+p+j]=1.0;
			}
			
			for(int j=0;j<w;j++) time_series[index+p+j]=1.0;
		}
			
			
		else{
			rewrite=false;
			for(int j=0;j<w;j++) time_series[index+p+j]=1.0;
		}
	}
	
	shmdt(time_series);
	
	return 0;

}

