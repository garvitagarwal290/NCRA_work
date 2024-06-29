#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

int main(int argc, char* argv[])
{
	char *timeseries_path=(char*)malloc(1*strlen(argv[1]));
	FILE *f;
	int N=0,i=0;
	float val, *time_series= (float*) malloc(sizeof(float)*1);
	long int totalsize;

	float phase_offset=-1.0, *bins=(float*) malloc(sizeof(float)*1);
	double P,tau,bin_interval;
	int num_bins,bin_index, *bin_count=(int*) malloc(sizeof(int)*1);
	int buffer_binsize;

	float *folded=(float*) malloc(sizeof(float)*1),*temp=(float*) malloc(sizeof(float)*1),sum;
	int bin_offset;

	FILE *gnuplot = popen("gnuplot -persistent", "w");
	FILE *foldedfile = fopen("folded.txt", "w");
	int numcommands= 4;
	//char *gnuplotcommands[] = {"set title \"Folded\";set terminal png size 500,400;set output 'folded.png';plot \"folded.txt\" u 1:2 w lp"};
	char *gnuplotcommands[] = {"set title \"Folded\"","set terminal png size 800,600", "set output 'folded.png'" ,"plot \"folded.txt\" u 1:2 w lp"};


//loading the binary file data----------------------------------------------------

	strcpy(timeseries_path,argv[1]);
	f= fopen(timeseries_path,"rb");

	//calculating total number of samples N	
	fseek(f,0,SEEK_END);
	totalsize= ftell(f);
	N= totalsize/sizeof(float);
	
	//putting data in time_series
	fseek(f,0,SEEK_SET);
	time_series=(float*)realloc(time_series,sizeof(float)*N);
	for(i=0; i<N ;i++)
	{
		fread(&val, sizeof(float), 1, f);
		time_series[i]=val;
	}

	fclose(f);

//Folding the time series---------------------------------------------------------

//info needed for folding

	//commandline arguments	
	P = strtod(argv[2],NULL);
	tau=strtod(argv[3],NULL);
	if(argc==5) phase_offset=strtof(argv[4],NULL);

	
	num_bins=(int)P/tau;
	bin_interval= (double)P/num_bins;
	bin_count=(int*)realloc(bin_count,num_bins*sizeof(int));
	for(i=0;i<num_bins;i++) bin_count[i]=0;
	
	//buufer_binsize keeps extra space for each bin in the 'bins' array
	buffer_binsize= (int) (N*1.1)/num_bins;
	bins=(float*) realloc(bins, num_bins*buffer_binsize*sizeof(float)); 
	for(i=0;i<N;i++) bins[i]=0.0;

	//the actual folding loop
	for(i=0;i<N;i++)
	{
	bin_index =(int)round((fmod(i*tau,P))/bin_interval); //decides which bin the ith sample goes to
	if (bin_index == num_bins) bin_index=0;
	bin_count[bin_index]+=1;
	
	bins[bin_index*buffer_binsize+ bin_count[bin_index]-1]=time_series[i];
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

	//adjusting for the phase_offset in the folded profile
	if (phase_offset != -1)
	{
		bin_offset= (int)round(num_bins*phase_offset);
		temp= (float*) realloc(temp, num_bins*sizeof(float));

		for(i=0;i<num_bins - bin_offset; i++) temp[bin_offset + i] = folded[i];
		for(i=0;i<bin_offset; i++) temp[i] = folded[num_bins- bin_offset+i]; 
		folded=(float*) realloc(temp, num_bins*sizeof(float));	
	}	

	//saving the folded profile in a textfile

	char header[500];
	sprintf(header, "#Path                             :%s\n#Period                           :%lf\n#Sampling Interval                :%lf\n#phase_offset                     :%f\n#Total no. of lines is 5(including this line) and the numeric values at each line start from the 36th character\n",timeseries_path,P,tau,phase_offset);
	fprintf(foldedfile, header);
	
	for(i=0;i<num_bins;i++)
	{
	fprintf(foldedfile, "%f %f\n", i/((float)num_bins), folded[i]);
	}
	
	// automatically plotting the folded profile using gnuplot
	for(i=0;i<numcommands ;i++)
	{
	fprintf(gnuplot, "%s\n", gnuplotcommands[i]);
	}
		

	return 0;
}
