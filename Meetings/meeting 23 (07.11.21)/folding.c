#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<sys/shm.h>
#include<stdbool.h>
#include<unistd.h>

int main()
{
	int i=0;

	float *bins=(float*) malloc(sizeof(float)*1);
	double P,tau,bin_interval;
	int num_bins,bin_index, *bin_count=(int*) malloc(sizeof(int)*1);
	int buffer_binsize;

	float *folded=(float*) malloc(sizeof(float)*1),*temp=(float*) malloc(sizeof(float)*1),sum;

	FILE *gnuplot = popen("gnuplot -persistent", "w");
	
	FILE *foldedfile = fopen("folded_profile.txt", "w");
	int numcommands= 5;
	char title[100] = "set title 'Folded profile- "; 
	char *gnuplotcommands[] = {title,"set terminal png size 800,600", "set output 'folded_profile.png'","set xlabel 'phase'" ,"plot \"folded_profile.txt\" u 1:2 w lp"};


//retrieving the time series data from the shared memory ---------------------------
	
	//int N=10000000, 
	int mem_size=1000000;
	
	key_t key = ftok("shm_folding",65);

	// shmget returns an identifier in shmid
	int shmid = shmget(key,sizeof(float)*mem_size,0666);

	// shmat to attach to shared memory
	float *time_series=(float*) shmat(shmid,NULL,0);	
	
	
//Folding the time series---------------------------------------------------------

//info needed for folding

	P = 500.12345;
	tau=0.234;

	num_bins=(int)P/tau;
	bin_interval= (double)P/num_bins;
	bin_count=(int*)realloc(bin_count,num_bins*sizeof(int));
	for(i=0;i<num_bins;i++) bin_count[i]=0;
	
	//buffer_binsize keeps extra space for each bin in the 'bins' array
	buffer_binsize= 10000;//(int) (N*1.1)/num_bins;
	bins=(float*) realloc(bins, num_bins*buffer_binsize*sizeof(float)); 
	for(i=0;i<buffer_binsize*num_bins;i++) bins[i]=0.0;
	
	int start_over=0; //similar to 'rewrite' from the writer code.
	*time_series=0.2; // != 0 is a signal for the reader to wait

	//the actual folding loop
	for(i=0;i>-1;i++)
	{
	
	while(*time_series != 0.0) 
	{
		printf("Awaiting new data\n");
		sleep(1);
	}
	
	if(i-start_over*mem_size==1) printf("folding new shm data\n");
	
	if(i-start_over*mem_size>=mem_size) //if fully read shm then stop reading until *time_series is set to 0.0 again
	{
		start_over++;
		*time_series=0.2; //stop signal
	}
	
	bin_index =(int)round((fmod(i*tau,P))/bin_interval); //decides which bin the ith sample goes to
	if (bin_index == num_bins) bin_index=0;
	bin_count[bin_index]+=1;
	
	bins[bin_index*buffer_binsize+ bin_count[bin_index]-1]=time_series[i-start_over*mem_size];
	
	if(start_over==5) break;
	}
	
	//summing all values in each bin and taking average to get folded profile array	
	folded= (float*)realloc(folded,num_bins*sizeof(float));
	for(i=0;i<num_bins; i++) folded[i]=0.0;

	for(i=0;i<num_bins;i++)
	{
		sum=0.0;
		for(int j=0;j<bin_count[i];j++) sum+=bins[i*buffer_binsize+j];
		folded[i]=sum/bin_count[i];
	}


	//saving the folded profile in a textfile

	char header[500];
	sprintf(header, "#Period                           :%lf\n#Sampling Interval                :%lf\n#Total no. of lines is 3(including this line) and the numeric values at each line start from the 36th character\n",P,tau);
	fprintf(foldedfile, header);
	
	for(i=0;i<num_bins;i++)
	{
	fprintf(foldedfile, "%f %f\n", i/((float)num_bins), folded[i]);
	}
	
	// automatically saving the folded profile as png using gnuplot
	for(i=0;i<numcommands ;i++)
	{
	fprintf(gnuplot, "%s\n", gnuplotcommands[i]);
	}
	
	FILE *shotwell = popen("shotwell folded_profile.png","r");
	
	
	//detach from shared memory 
    	shmdt(time_series);
    
    	// destroy the shared memory
    	//shmctl(shmid,IPC_RMID,NULL);
		

	return 0;
}
