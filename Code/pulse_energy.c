#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include<stdbool.h>

//the format of the run command is as follows: ./[path-to-code's-executable-after-compiling] [-p num_onpulse_phaseranges] [-s num_offpulse stats phaseranges] [-m num_offpulse mean phase ranges] file path

//function for slicing a string. Taken from stack overflow
void slice_str(const char * str, char * buffer, size_t start, size_t end)
{
    size_t j = 0;
    for ( size_t i = start; i <= end; ++i ) {
        buffer[j++] = str[i];
    }
    buffer[j] = 0;
}

int main(int argc, char* argv[])
{	
	//taking inputs from the command line arguments using getopt. Usage is simply: [-p num_onpulse_phaseranges] [-s num_offpulse stats phaseranges] [-m num_offpulse mean phase ranges] file path. All arguments are mandatory
	int num_onpulse_phaseranges=0,num_offpulsestats_ranges=0, num_offpulsemean_ranges=0,opt;
	while ((opt = getopt(argc, argv, "p:m:s:f:")) != -1) {
		switch (opt){
			case'p':
				num_onpulse_phaseranges= atoi(optarg);
				break;
			case's':
				num_offpulsestats_ranges= atoi(optarg);
				break;
			case'm':
				num_offpulsemean_ranges= atoi(optarg);
				break;
			case'f':
				num_offpulsemean_ranges= atoi(optarg);
				break;
			default: /* '?' */
				fprintf(stderr, "Usage: [-p num_onpulse_phaseranges] [-s num_offpulse stats phaseranges] [-m num_offpulse mean phase ranges] file path\n");
				exit(EXIT_FAILURE);
		}
	}

	//printf("%s",argv[optind]);
	printf("%d %d\n",optind,argc);
	
	//sanity checks on the arguments
	if (optind >= argc) {
	fprintf(stderr, "path to pulses file required\n");
	exit(EXIT_FAILURE);
	}
	if(num_onpulse_phaseranges==0){
	fprintf(stderr, "number of on pulse phase ranges is required. Use -p\n");
	exit(EXIT_FAILURE);
	}
	if(num_offpulsestats_ranges==0){
	fprintf(stderr, "number of off pulse phase ranges for stats is required. Use -s\n");
	exit(EXIT_FAILURE);
	}
	if(num_offpulsemean_ranges==0){
	fprintf(stderr, "number of offpulse phase ranges for mean is required. Use -m\n");
	exit(EXIT_FAILURE);
	}
	
//manually putting the pulse and noise phase ranges for getting noise stats------------------------
	//asking for inputting the pulse phase range num_onpulse_phaseranges number of times
	//float onpulse_phaserangesmin[2]={.14, .34}, onpulse_phaserangesmax[2]={.27,.65}; //for J2145-0750
	//float onpulse_phaserangesmin[1]={.809}, onpulse_phaserangesmax[1]={.8905}; //for B0329+54
	//float onpulse_phaserangesmin[1]={.841}, onpulse_phaserangesmax[1]={.919}; //for B0740-28
	//float onpulse_phaserangesmin[2]={.368,.828}, onpulse_phaserangesmax[2]={.391,.865}; //for B1642-03
	float onpulse_phaserangesmin[num_onpulse_phaseranges],onpulse_phaserangesmax[num_onpulse_phaseranges];
	for(int i=0;i<num_onpulse_phaseranges;i++)
	{
		printf("\nEnter min and max for on pulse phase range %d: ",i+1);
		scanf("%f %f",&onpulse_phaserangesmin[i], &onpulse_phaserangesmax[i]);
	}

	//asking for inputting the noise phase range num_offpulse_phaseranges number of times
	//float offpulsestats_rangesmin[1]={0.7275}, offpulsestats_rangesmax[1]={0.809}; //for B0329+54
	//float offpulsemean_rangesmin[2]={0.0,0.8905}, offpulsemean_rangesmax[2]={.7275,1.0}; //for B0329+54
	//float offpulsestats_rangesmin[1]={0.763}, offpulsestats_rangesmax[1]={0.841}; //for B0740-28
	//float offpulsemean_rangesmin[2]={0.0,0.919}, offpulsemean_rangesmax[2]={.763,1.0}; //for B0740-28
	//float offpulsestats_rangesmin[1]={0.768}, offpulsestats_rangesmax[1]={0.828}; //for B1642-03
	//float offpulsemean_rangesmin[3]={0.0,0.391,.865}, offpulsemean_rangesmax[3]={.368,.768,1.0}; //for B1642-03
	//float offpulsestats_rangesmin[3]={0.0,0.27,0.65}, offpulsestats_rangesmax[3]={.14,.34,0.88}; //for J2145-
	//float offpulsemean_rangesmin[1]={.88}, offpulsemean_rangesmax[1]={1.0}; //for B1642-03

	float offpulsestats_rangesmin[num_offpulsestats_ranges],offpulsestats_rangesmax[num_offpulsestats_ranges];
	for(int i=0;i<num_offpulsestats_ranges;i++)
	{
		printf("\nEnter min and max for off pulse stats phase range %d: ",i+1);
		scanf("%f %f",&offpulsestats_rangesmin[i], &offpulsestats_rangesmax[i]);
	}

	float offpulsemean_rangesmin[num_offpulsemean_ranges],offpulsemean_rangesmax[num_offpulsemean_ranges];
	for(int i=0;i<num_offpulsemean_ranges;i++)
	{
		printf("\nEnter min and max for off pulse mean phase range %d: ",i+1);
		scanf("%f %f",&offpulsemean_rangesmin[i], &offpulsemean_rangesmax[i]);
	}

	//num_onpulse_phaseranges=2;
	//num_offpulsestats_ranges=3;
	//int num_offpulsemean_ranges=1;

//-----------------------------------------------------------------------------

	//reading the header from the pulses file (first 10 lines) 
	char *pulses_path= (char*) malloc(1* strlen(argv[optind]));
	strcpy(pulses_path,argv[optind]);
	
	char *line,*header_info[11],buffer[300];
	for(int i=0;i<11;i++) header_info[i]=malloc(sizeof(char)*150);
	FILE *f;
	size_t len;
	
	f= fopen(pulses_path,"rb");
	for(int i=0;i<11;i++) 
	{
		getline(&line, &len, f);
		strcpy(header_info[i],line);
	}

	//getting the relevant numbers from the header
	int num_pulses, num_pulsesamples;
	slice_str(header_info[7],buffer,35,strlen(header_info[7]));
	num_pulses= atoi(buffer);
	slice_str(header_info[8],buffer,35,strlen(header_info[8]));
	num_pulsesamples= atoi(buffer);
	
	//reading the data and calculating the sum of sample and noise values in each given phase range. Also calculating the average noise power in a single pulse
	float *onpulse_sum[num_pulses], phase, *offpulse_avg=(float*) malloc(sizeof(float)*num_pulses), *offpulsestats_sum=(float*) malloc(sizeof(float)* num_pulses);
	double val;

	int num_offpulsesamples_mean, num_offpulsesamples_stats; //no of samples considered to be noise in the whole pulse (for finding offpulse_avg energy in a pulse and for calculating off-pulse energy in a pulse respectively)  

	//this array has number of samples in each of the given onpulse phase ranges	
	int num_onpulsesamples[num_onpulse_phaseranges];


	for(int i=0;i<num_pulses;i++)
	{
		onpulse_sum[i]=(float*)malloc(sizeof(float)*num_onpulse_phaseranges);
		num_offpulsesamples_mean=0;
		num_offpulsesamples_stats=0;

		for(int k=0;k<num_onpulse_phaseranges;k++) num_onpulsesamples[k]= 0; 
		
		for(int j=0;j<num_pulsesamples;j++)
		{	
			fread(&val, sizeof(double), 1, f);
			phase = (j+1)/((float)num_pulsesamples);
			
			//this loop checks if the jth sample belongs to any of the given onpulse phase ranges
			for(int k=0;k<num_onpulse_phaseranges;k++)
			{
				if(phase>onpulse_phaserangesmin[k] && phase<onpulse_phaserangesmax[k]) 
				{
				onpulse_sum[i][k]+=(float) val;
				num_onpulsesamples[k]++;
				}			
			}

			//this loop checks if the jth sample belongs to any of the given offpulse phase ranges
			for(int k=0;k<num_offpulsestats_ranges;k++)
			{
				if(phase>offpulsestats_rangesmin[k] && phase<offpulsestats_rangesmax[k]) 
				{
				offpulsestats_sum[i]+=(float) val;
				num_offpulsesamples_stats++;
				}
			}

			for(int k=0;k<num_offpulsemean_ranges;k++)
			{
				if(phase>offpulsemean_rangesmin[k] && phase<offpulsemean_rangesmax[k]) 
				{
				offpulse_avg[i]+= (float) val;
				num_offpulsesamples_mean++;
				}			
			}
		}
		offpulse_avg[i]/=num_offpulsesamples_mean; //finding average noise level for the ith pulse	
		
	}
	//printf("%d %d\n",num_offpulsesamples_mean, num_offpulsesamples_stats);

	//printf("%d %d\n",num_onpulsesamples[0], num_onpulsesamples[1]);	

	//finally calculating the energy *only* due to pulsar in each onpulse phase range
	//and subtracting the offpulse_avg from offpulsestats_sum to get offpulse_energy array
	float *onpulse_energy[num_pulses], *offpulse_energy=(float*) malloc(sizeof(float)* num_pulses);
	
	float offpulseenergy_distributionmean=0.0,offpulseenergy_distributionsigma=0.0;
	for(int i=0;i<num_pulses;i++) 
	{	
		onpulse_energy[i]=(float*) malloc(sizeof(float)*num_onpulse_phaseranges);
		
		for(int k=0;k<num_onpulse_phaseranges;k++)
		{
		onpulse_energy[i][k] = onpulse_sum[i][k] - offpulse_avg[i]*num_onpulsesamples[k]; //subtracting the common noise level from the kth phase range in the ith pulse
		}

		offpulse_energy[i]=offpulsestats_sum[i] - offpulse_avg[i]*num_offpulsesamples_stats; //subtracting the common noise level from the summed noise samples. basically sets the distribution of noise samples to 0
		offpulseenergy_distributionmean+=offpulse_energy[i];
				
	}
	offpulseenergy_distributionmean/=num_pulses;
	
	for(int i=0;i<num_pulses;i++) 
	{
	offpulseenergy_distributionsigma+= (offpulseenergy_distributionmean- offpulse_energy[i])*(offpulseenergy_distributionmean- offpulse_energy[i]);
	}
	offpulseenergy_distributionsigma/=num_pulses;
	offpulseenergy_distributionsigma= sqrtf(offpulseenergy_distributionsigma);

	//FILE *noisestats= fopen("../noise_stats/noisestats.txt", "a");
	//fprintf(noisestats,"%d %f %f\n",num_offpulsesamples_stats,offpulseenergy_distributionmean,offpulseenergy_distributionsigma);
	

	//making the header for the output onpulse_energy.txt file
	char header[800];
	sprintf(header, "#Path                             :%s\n#no. of pulses                    :%d\n#no. of on-pulse phase ranges     :%d\n",pulses_path,num_pulses, num_onpulse_phaseranges);

	for(int i=0;i<num_onpulse_phaseranges;i++)
	{
	sprintf(buffer,"#on-pulse phase range %d min       :%f\n",i+1,onpulse_phaserangesmin[i]);
	strcat(header,buffer);
	sprintf(buffer,"#on-pulse phase range %d max       :%f\n",i+1,onpulse_phaserangesmax[i]);
	strcat(header,buffer);
	}

	for(int i=0;i<num_offpulsestats_ranges;i++)
	{
	sprintf(buffer,"#off-pulse stats phase range %d min:%f\n",i+1,offpulsestats_rangesmin[i]);
	strcat(header,buffer);
	sprintf(buffer,"#off-pulse stats phase range %d max:%f\n",i+1,offpulsestats_rangesmax[i]);
	strcat(header,buffer);
	}

	for(int i=0;i<num_offpulsemean_ranges;i++)
	{
	sprintf(buffer,"#off-pulse avg phase range %d min  :%f\n",i+1,offpulsemean_rangesmin[i]);
	strcat(header,buffer);
	sprintf(buffer,"#off-pulse avg phase range %d max  :%f\n",i+1,offpulsemean_rangesmax[i]);
	strcat(header,buffer);
	}	
	
	sprintf(buffer,"#Total no. of lines is %d(including this line) and the numeric values at each line start from the 36th character\n",3+2*(num_onpulse_phaseranges+num_offpulsestats_ranges+num_offpulsemean_ranges)+1);
	strcat(header,buffer);

	
	//making 3 txt files for pulse energy, noise energy and avg noise level. writing the header and energy values to the files
	FILE *onpulse_energy_file = fopen("onpulse_energy.txt", "w"), *offpulse_energy_file = fopen("offpulse_energy.txt", "w"), *offpulse_avgfile= fopen("offpulse_avg.txt", "w");
	fprintf(onpulse_energy_file, header); fprintf(offpulse_avgfile, header);

	sprintf(buffer, "#offpulsesamples forStats(9thline):%d\n",num_offpulsesamples_stats);
	strcat(header,buffer);	
	fprintf(offpulse_energy_file, header);

	for(int i=0;i<num_pulses;i++)
	{
		fprintf(onpulse_energy_file, "%d",i+1);fprintf(offpulse_energy_file, "%d",i+1); fprintf(offpulse_avgfile, "%d",i+1);

		for(int k=0;k<num_onpulse_phaseranges;k++)
		{
		fprintf(onpulse_energy_file, " %f", onpulse_energy[i][k]);
		}
		fprintf(onpulse_energy_file, "\n");

		fprintf(offpulse_energy_file, " %f\n", offpulse_energy[i]);
		fprintf(offpulse_avgfile, " %f\n",offpulse_avg[i]);
	}
	
	// automatically saving and opening the onpulse_energy plot as png using gnuplot
	FILE *onpulse_gnuplot = popen("gnuplot -persistent", "w"); 
	FILE *shotwell = popen("shotwell onpulse_energy.png","r");
	int numcommands= 5;
	
	char title[100] = "set title 'on-pulse energy - "; strcat(title,pulses_path);
	char *gnuplotcommands[numcommands]; 
	gnuplotcommands[0]=title; gnuplotcommands[1]="set xlabel 'pulse number'"; gnuplotcommands[2]="set terminal png size 1000,600"; gnuplotcommands[3]= "set output 'onpulse_energy.png'"; 
	gnuplotcommands[4]=(char*) malloc(sizeof(char)*200); strcpy(gnuplotcommands[4],"plot ");

	char string[50];
	for (int i=0; i< num_onpulse_phaseranges; i++) 
	{
	
	sprintf(string,"\"onpulse_energy.txt\" u 1:%d w lp title \"Phase %.4f - %.4f\",",2+i, onpulse_phaserangesmin[i], onpulse_phaserangesmax[i]);
	strcat(gnuplotcommands[4],string);
	}

	for(int i=0;i<numcommands ;i++)
	{
	fprintf(onpulse_gnuplot, "%s\n", gnuplotcommands[i]);
	}


	// automatically saving the noise energy plot as png using gnuplot
	FILE *offpulse_energy_gnuplot = popen("gnuplot -persistent", "w");
	numcommands= 5;
	strcpy(title,"set title 'off-pulse energy - "); strcat(title,pulses_path);
	sprintf(string,"  |  mean - %f | rms - %f",offpulseenergy_distributionmean,offpulseenergy_distributionsigma); strcat(title,string);	

	*gnuplotcommands[numcommands]; 
	gnuplotcommands[0]=title; gnuplotcommands[1]="set xlabel 'pulse number'"; gnuplotcommands[2]="set terminal png size 1000,600"; gnuplotcommands[3]= "set output 'offpulse_energy.png'"; gnuplotcommands[4]="plot \"offpulse_energy.txt\" u 1:2 w lp"; 

	for(int i=0;i<numcommands ;i++)
	{
	fprintf(offpulse_energy_gnuplot, "%s\n", gnuplotcommands[i]);
	}
		

	// automatically saving the average  average noise level plot as png using gnuplot
	FILE *offpulse_avg_gnuplot = popen("gnuplot -persistent", "w");
	numcommands= 5;
	strcpy(title,"set title 'off-pulse average energy level - "); strcat(title,pulses_path);
	
	*gnuplotcommands[numcommands]; 
	gnuplotcommands[0]=title; gnuplotcommands[1]="set xlabel 'pulse number'"; gnuplotcommands[2]="set terminal png size 1000,600"; gnuplotcommands[3]= "set output 'offpulse_avg.png'"; gnuplotcommands[4]="plot \"offpulse_avg.txt\" u 1:2 w lp"; 

	for(int i=0;i<numcommands ;i++)
	{
	fprintf(offpulse_avg_gnuplot, "%s\n", gnuplotcommands[i]);
	}
			
	return 0;
}		


	
