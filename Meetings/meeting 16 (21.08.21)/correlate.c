#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include<stdbool.h>

//the rationale behind the variable names is as follows: we are calculating the correlation of 2 functions(could be identical), called f and g which would appear as arrays having values of the pulse energy functions (belonging to different phase ranges). 
//Note: this code assumes that the given file as 2 columns. And one cant have f and g coming from two different files.

int num_pulses;

//function for slicing a string. Taken from stack overflow
void slice_str(const char * str, char * buffer, size_t start, size_t end)
{
    size_t j = 0;
    for ( size_t i = start; i <= end; ++i ) {
        buffer[j++] = str[i];
    }
    buffer[j] = 0;
}

//function for 'shifting' an array by shiftnum places and saving the modified array in a different array pointer. Inspired from http://www.jot.fm/issues/issue_2010_03/column2.pdf
void shift(float* array,float* shifted_array, int shiftnum) 
{	 
	for (int i = 0; i < num_pulses; i++) 
	{
		if ((shiftnum + i >= 0)&&(shiftnum+i < num_pulses)) shifted_array[i] = array[shiftnum + i];
		else shifted_array[i]=0;
	}
} 

//function for normalising a array ie shifting the values by their mean and diving them by their rms
void normalise(float* array,float* normalised_array)
{
	float mean=0.0, rms=0.0;
	for(int i=0;i<num_pulses;i++) mean+=array[i];
	mean/=num_pulses;

	for(int i=0;i<num_pulses;i++)	
	{
		rms+=(mean-array[i])*(mean-array[i]);
		normalised_array[i]=array[i]-mean;
	}
	rms/=num_pulses;
	rms=sqrtf(rms);

	for(int i=0;i<num_pulses;i++) normalised_array[i]/=rms;
}

//return the 'dot product' of the f function with the (shifted) g function in every step of the correlation loop.
float dot(float* f,float* shifted_g) 
{
	float sum = 0.0;
	for (int i = 0; i < num_pulses; i++) sum += f[i] * shifted_g[i];
	
	return sum;
} 

int main(int argc, char* argv[])
{	
	//taking inputs from the command line arguments using getopt. Usage is [-f column number of the data to be used as the f function] [-g column number of the data to be used as the g function] [path of file having data for both functions]. Both arguments are mandatory
	int totalcols=2,fcol=0,gcol=0, opt;
	while ((opt = getopt(argc, argv, "f:g:")) != -1) {
		switch (opt){
		
			case'f':
				fcol=atoi(optarg);
				break;
			case'g':
				gcol=atoi(optarg);
				break;
			default: /* '?' */
				fprintf(stderr, "Usage: [-p num_pulsephaseranges] [-n num_noisephaseranges] filepath\n");
				exit(EXIT_FAILURE);
		}
	}

	char *file_path=malloc(sizeof(char)*50);
	strcpy(file_path,argv[optind]);
	
//loading the complete header.
	int header_length=8; //no of lines in header
	char *line,*header_info[header_length],buffer[100];
	for(int i=0;i<header_length;i++) header_info[i]=(char*) malloc(sizeof(char)*100);
	FILE *file;
	size_t len;
	
	file= fopen(file_path,"r");
	for(int i=0;i<header_length;i++) 
	{
		getline(&line, &len, file);
		strcpy(header_info[i],line);
	}

//saving the number of total pulses from the header
	slice_str(header_info[1],buffer,35,strlen(header_info[1]));
	num_pulses= atoi(buffer);

	float f_phaseranges[2], g_phaseranges[2];
	for(int i=0;i<2;i++)
	{
		slice_str(header_info[3+(fcol-1)*2+i],buffer,35,strlen(header_info[3+(fcol-1)*2+i]));
		f_phaseranges[i]= atof(buffer);
	}
	for(int i=0;i<2;i++)
	{
		slice_str(header_info[3+(gcol-1)*2+i],buffer,35,strlen(header_info[3+(gcol-1)*2+i]));
		g_phaseranges[i]= atof(buffer);
	}


//loading data of both columns in data array
	float *data[totalcols];
	for(int i=0;i<totalcols;i++) data[i]=(float*) malloc(sizeof(float)*num_pulses);	
	int n=0;
	for(int i=0;i<num_pulses; i++) 
	{
		getline(&line, &len, file);
		sscanf(line,"%d %f %f\n",&n,&data[0][i],&data[1][i]);
	}

//making separate arrays for f and g function
	float *f=(float*) malloc(sizeof(float)*num_pulses),*g=(float*) malloc(sizeof(float)*num_pulses);
	for(int i=0;i<num_pulses;i++) 
	{
		f[i]=data[fcol-1][i];
		g[i]=data[gcol-1][i];
	}

//the actual correlation loop. saves the correlation function in 'result' array
	float *normalised_f= (float*) malloc(sizeof(float)* num_pulses);
	float *shifted_g= (float*) malloc(sizeof(float)* num_pulses), *normalised_g= (float*) malloc(sizeof(float)* num_pulses);
	float *result=(float*) malloc(sizeof(float) * (2*num_pulses - 1));

	normalise(f,normalised_f);
	//normalise(g,normalised_g);	
	for (int i = -num_pulses+1; i < num_pulses; i++)
	{
		shift(g,shifted_g, i);
		normalise(shifted_g,normalised_g);
		result[i+num_pulses-1]=dot(normalised_f,normalised_g)/num_pulses;
	}

//saving the result in a text file

	//making the header
	char filename[50],header[200];
	sprintf(header,"#file path: %s\n#number of pulses: %d\n#this is correlation function between pulse energy (function) of 2 phase ranges: (%.4f,%.4f) and (%.4f,%.4f)\n",file_path, num_pulses, f_phaseranges[0],f_phaseranges[1],g_phaseranges[0],g_phaseranges[1]);

	sprintf(filename,"correlated%d%d.txt",fcol,gcol);
	FILE *correlated=fopen(filename,"w");
	
	fprintf(correlated,header);	
	for(int i = -num_pulses+1; i < num_pulses; i++) fprintf(correlated,"%d %f\n",i,result[i+num_pulses-1]);

	return 0;

}
