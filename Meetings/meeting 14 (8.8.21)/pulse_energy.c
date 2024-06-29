#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include<stdbool.h>

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
	//taking inputs from the command line arguments using getopt. Usage is simply -n [no. of phase ranges] [data file]. Both arguments are mandatory
	int num_phaseranges=0, opt;
	while ((opt = getopt(argc, argv, "n:")) != -1) {
		switch (opt){
			case'n':
				num_phaseranges= atoi(optarg);
				break;
			default: /* '?' */
				fprintf(stderr, "Usage: [-n num_phaseranges] filepath\n");
				exit(EXIT_FAILURE);
		}
	}
	
	if (optind >= argc) {
	fprintf(stderr, "path to pulses file required\n");
	exit(EXIT_FAILURE);
	}
	if(num_phaseranges==0){
	fprintf(stderr, "number of phase ranges is required. Use -n\n");
	exit(EXIT_FAILURE);
	}
	
	//asking for inputting the phase range num_phaseranges number of times
	float phaseranges_min[num_phaseranges], phaseranges_max[num_phaseranges];
	for(int i=0;i<num_phaseranges;i++)
	{
		printf("\nEnter min and max for phase range %d: ",i+1);
		scanf("%f %f",&phaseranges_min[i], &phaseranges_max[i]);
	}

	//reading the header from the pulses file (first 10 lines) 
	char *pulses_path= (char*) malloc(1* strlen(argv[optind]));
	strcpy(pulses_path,argv[optind]);
	
	char *line,*header_info[11],buffer[100];
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
	float *pulse_sum[num_pulses], phase, *pulse_avgnoise=(float*) malloc(sizeof(float)*num_pulses), *noise_sum=(float*) malloc(sizeof(float)* num_pulses);
	double val;
	int num_noisesamples, reqnum_noisesamples=0.2* num_pulsesamples; //no of samples considered to be noise in the whole pulse (for finding average noise power in a pulse and for calculating noise energy in a pulse respectively)  

	for(int i=0;i<num_pulses;i++)
	{
		pulse_sum[i]=(float*)malloc(sizeof(float)*num_phaseranges);
		num_noisesamples=0; 
		
		for(int j=0;j<num_pulsesamples;j++)
		{	
			fread(&val, sizeof(double), 1, f);
			bool sampleis_noise=true; //helps decide if the jth sample is a noise value or whether it falls inside the given phase ranges
			phase = j/((float)num_pulsesamples);
			
			//this loop checks if the jth sample belongs to any of the given phase ranges
			for(int k=0;k<num_phaseranges;k++)
			{
				if(phase>phaseranges_min[k] && phase<phaseranges_max[k]) 
				{
				pulse_sum[i][k]+=(float) val;
				sampleis_noise=false;
				}			
			}

			if(sampleis_noise) 
			{
			pulse_avgnoise[i]+= (float) val;
			num_noisesamples++;
			}
			if(num_noisesamples== reqnum_noisesamples) noise_sum[i]=pulse_avgnoise[i]; //if we have summed reqnum_noisesamples then we save the summed noise values before more noise samples are added to pulse_avgnoise[i] 
		}
		
		pulse_avgnoise[i]/=num_noisesamples; //finding average noise level for the ith pulse	
		
	}
	//printf("%d %d\n",num_noisesamples, reqnum_noisesamples);

	//this array has number of samples in each of the given phase ranges
	int samples_phaseranges[num_phaseranges];
	for(int k=0;k<num_phaseranges;k++) samples_phaseranges[k]=(phaseranges_max[k]-phaseranges_min[k])*num_pulsesamples;	
	
	//finally calculating the energy *only* due to pulsar in each phase range
	//and subtracting the avg_noise from noise_sum to get noise_energy array
	float *pulse_energy[num_pulses], *noise_energy=(float*) malloc(sizeof(float)* num_pulses);
	for(int i=0;i<num_pulses;i++) 
	{	
		pulse_energy[i]=(float*) malloc(sizeof(float)*num_phaseranges);
		for(int k=0;k<num_phaseranges;k++)
		{
		pulse_energy[i][k] = pulse_sum[i][k] - pulse_avgnoise[i]*samples_phaseranges[k]; //subtracting the common noise level from the kth phase range in the ith pulse
		}
		noise_energy[i]=noise_sum[i] - pulse_avgnoise[i]*reqnum_noisesamples;	//subtracting the common noise level from the summed noise samples. basically sets the distribution of noise samples to 0
	}
	
	//making the header for the output pulse_energy.txt file
	char header[500];
	sprintf(header, "#Path                             :%s\n#no. of pulses                    :%d\n#no. of phase ranges              :%d\n",pulses_path,num_pulses, num_phaseranges);

	for(int i=0;i<num_phaseranges;i++)
	{
	sprintf(buffer,"#Phase range %d min                :%f\n",i+1,phaseranges_min[i]);
	strcat(header,buffer);
	sprintf(buffer,"#Phase range %d max                :%f\n",i+1,phaseranges_max[i]);
	strcat(header,buffer);
	}
	
	sprintf(buffer,"#Total no. of lines is %d(including this line) and the numeric values at each line start from the 36th character\n",3+2*num_phaseranges+1);
	strcat(header,buffer);

	
	//making 3 txt files for pulse energy, noise energy and avg noise level. writing the header and energy values to the files
	FILE *pulse_energy_file = fopen("pulse_energy.txt", "w"), *noise_energy_file = fopen("noise_energy.txt", "w"), *avg_noise_file= fopen("avg_noiselevel.txt", "w");
	fprintf(pulse_energy_file, header); fprintf(avg_noise_file, header);

	sprintf(buffer, "#no.of noise sample added(9thline):%d\n",reqnum_noisesamples);
	strcat(header,buffer);	
	fprintf(noise_energy_file, header);

	for(int i=0;i<num_pulses;i++)
	{
		fprintf(pulse_energy_file, "%d",i+1);fprintf(noise_energy_file, "%d",i+1); fprintf(avg_noise_file, "%d",i+1);

		for(int k=0;k<num_phaseranges;k++)
		{
		fprintf(pulse_energy_file, " %f", pulse_energy[i][k]);
		}
		fprintf(pulse_energy_file, "\n");

		fprintf(noise_energy_file, " %f\n", noise_energy[i]);
		fprintf(avg_noise_file, " %f\n",pulse_avgnoise[i]);
	}
	
	// automatically saving and opening the pulse_energy plot as png using gnuplot
	FILE *pulse_gnuplot = popen("gnuplot -persistent", "w"); 
	FILE *shotwell = popen("shotwell pulse_energy.png","r");
	int numcommands= 5;
	
	char title[100] = "set title 'Pulse energy - "; strcat(title,pulses_path);
	char *gnuplotcommands[numcommands]; 
	gnuplotcommands[0]=title; gnuplotcommands[1]="set xlabel 'pulse number'"; gnuplotcommands[2]="set terminal png size 1000,600"; gnuplotcommands[3]= "set output 'pulse_energy.png'"; 
	gnuplotcommands[4]=(char*) malloc(sizeof(char)*200); strcpy(gnuplotcommands[4],"plot ");

	for (int i=0; i< num_phaseranges; i++) 
	{
	char string[50];
	sprintf(string,"\"pulse_energy.txt\" u 1:%d w lp title \"Phase %.3f - %.3f\",",2+i, phaseranges_min[i], phaseranges_max[i]);
	strcat(gnuplotcommands[4],string);
	}

	for(int i=0;i<numcommands ;i++)
	{
	fprintf(pulse_gnuplot, "%s\n", gnuplotcommands[i]);
	}


	// automatically saving (and opening?) the noise energy plot as png using gnuplot
	FILE *noise_energy_gnuplot = popen("gnuplot -persistent", "w");
	numcommands= 5;
	strcpy(title,"set title 'Noise Energy - "); strcat(title,pulses_path);
	
	*gnuplotcommands[numcommands]; 
	gnuplotcommands[0]=title; gnuplotcommands[1]="set xlabel 'pulse number'"; gnuplotcommands[2]="set terminal png size 1000,600"; gnuplotcommands[3]= "set output 'noise_energy.png'"; gnuplotcommands[4]="plot \"noise_energy.txt\" u 1:2 w lp"; 

	for(int i=0;i<numcommands ;i++)
	{
	fprintf(noise_energy_gnuplot, "%s\n", gnuplotcommands[i]);
	}
		

	// automatically saving (and opening?) the average  average noise level plot as png using gnuplot
	FILE *avg_noise_gnuplot = popen("gnuplot -persistent", "w");
	numcommands= 5;
	strcpy(title,"set title 'Average Noise level - "); strcat(title,pulses_path);
	
	*gnuplotcommands[numcommands]; 
	gnuplotcommands[0]=title; gnuplotcommands[1]="set xlabel 'pulse number'"; gnuplotcommands[2]="set terminal png size 1000,600"; gnuplotcommands[3]= "set output 'avg_noiselevel.png'"; gnuplotcommands[4]="plot \"avg_noiselevel.txt\" u 1:2 w lp"; 

	for(int i=0;i<numcommands ;i++)
	{
	fprintf(avg_noise_gnuplot, "%s\n", gnuplotcommands[i]);
	}
			
	return 0;
}		


	
