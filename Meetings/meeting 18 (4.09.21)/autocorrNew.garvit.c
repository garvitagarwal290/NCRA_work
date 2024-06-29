#define _FILE_OFFSET_BITS 64
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "cpgplot.h"
#define lowFreq 500.00000000
#define upFreq 400.0000000
#define subPulseFreqThreshold 2
#define powerThreshold 3
#define spectraPlotThresh 5
#define microThreshold 10
#include<fftw3.h>

struct visualPulse
{
	int pulse;
	int microatlowfreq;
	int microathighfreq;
};

void plotone_single(char *outfile, int pulseNum, int polFlag) 			//this function is used by the single frequency portion for both plotting the folded profile as well as individual pulses
{
	int numSamples=0,i;							//Outfile is the file which contains the folded/single pulse data
	float t,signal,maxP=0,minP=0;
	FILE *fileID;
	fileID=fopen(outfile,"r");
	
	while (fscanf(fileID,"%f %f",&t,&signal)==2)
	{
		if (numSamples==1)
			minP=signal;
		numSamples++;
	}
	rewind(fileID);
	float times[numSamples];
	float power[numSamples];

	for(i=0;i<numSamples;i++)
	{
		fscanf(fileID,"%f %f",&t,&signal);
		times[i]=t;
		power[i]=signal;
		if (power[i]>maxP)
		{
			maxP=power[i];
		}
		else if (power[i]<minP)
		{
			minP=power[i];
		}
	}
	//cpgsubp(1,1);
	cpgenv(times[0],times[numSamples-1],minP,maxP,0,0);
	char title[50];
	if(polFlag==1)
		sprintf(title,"Pulse Number %d",pulseNum);
	else if(polFlag==2)
		sprintf(title,"Folded profile");
		
	cpglab("Phase","Intensity (arbitrary units)",title);

	cpgline(numSamples,times,power);
}

void plottwo_single(char *fold, char *pulse, int pulseNum)		//This is the code used by single frequency portion in ordr to plot the folded profile along with the single pulses
{
	int numSamples=0,i;
        float t,signal,maxP=0,minP=0;
        FILE *fileID;
        fileID=fopen(pulse,"r");

        while (fscanf(fileID,"%f %f",&t,&signal)==2)
        {
                if (numSamples==1)
                        minP=signal;
                numSamples++;
        }
        rewind(fileID);
        float times[numSamples];
        float power[numSamples];

        for(i=0;i<numSamples;i++)
        {
                fscanf(fileID,"%f %f",&t,&signal);
                times[i]=t;
                power[i]=signal;
                if (power[i]>maxP)
                {
                        maxP=power[i];
                }
                else if (power[i]<minP)
                {
                        minP=power[i];
                }
        }
	
	cpgsubp(1,2);
        cpgenv(times[0],times[numSamples-1],minP,maxP,0,0);
        char title[50];

	sprintf(title,"Pulse Number %d",pulseNum);
	
	cpglab("Phase","Intensity (arbitrary units)",title);
        cpgline(numSamples,times,power);
	
	minP=maxP;
	fclose(fileID);
	maxP=0;
	fileID=fopen(fold,"r");

	for(i=0;i<numSamples;i++)
        {
                fscanf(fileID,"%f %f",&t,&signal);
                times[i]=t;
                power[i]=signal;
                if (power[i]>maxP)
                {
                        maxP=power[i];
                }
                else if (power[i]<minP && power[i]!=0)
                {
                        minP=power[i];
                }
        }
	cpgpanl(1,1);
	cpgenv(times[0],times[numSamples-1],minP,maxP,0,0);
	cpglab("Phase","Intensity (arbitrary units)","Folded Profile");

	cpgline(numSamples,times,power);
}

float absgal(float num)
{
	if(num<0)
		return -num;
	else
		return num;
}

float maxMod(float *folda, int last)			//folda is the double equivalent of the float intensities.
{
        int i=0;
        float max=0;
        for(i=0;i<last;i++)
        {
                if (absgal(folda[i])>max)
                        max=folda[i];
        }

        return max;
}

double diffTime(int hour_1,int min_1, double sec_1, int hour_2, int min_2, double sec_2)		//calculating difference between the two timestamps		
{
        double secDiff=3600*(hour_2-hour_1)+60*(min_2-min_1)+sec_2-sec_1;
        return secDiff;
}

double calculateRelDelay(double startDiff, double DM,double timeOffset)					//total shift that needs to be corrected for
{
	double dispDelay=4.148808*pow((double)10,3)*(pow((double)lowFreq,-2)-pow((double)upFreq,-2))*DM;	
	printf("\nDispersion delay is %lf s for DM %lf pc/cc.\n",dispDelay, DM);
	double relDiscard=startDiff-dispDelay-timeOffset;
	return relDiscard;	
}

int readFloatToDouble(int numElements, float *readBuffer, double *store, FILE *readFile)
{
	if(fread(readBuffer,sizeof(float),numElements,readFile)!=numElements)
		return 0;
	
	int i;
	for (i=0;i<numElements;i++)
	{
		store[i]=(double)readBuffer[i];
	}
	return 1;
}

double meanBase(double *folda, double *phase, double pulseStart, double pulseEnd, int last)
{
	int i=0,count=0;
	double baseline=0, tempVal=0, tempPhase=0;//the baseline is the mean of the offpulse
	for (i=0;i<last;i++)
	{
		tempPhase=phase[i];
		if (tempPhase>=pulseStart && tempPhase<=pulseEnd)
			continue;
		else
		{
			tempVal=folda[i];
			baseline=(baseline*count+tempVal)/(count+1);
			count++;
		}
	}
	if (baseline==0)
		return -1;
	else
		return baseline;
}

double stdevBaseIntg(double *folda, double *phase, double pulseStart, double pulseEnd, int last, int intg)
{	// The main thing to notice is that the parameters used are pulseStart and pulseEnd which were defined right at the beginning of the code and are independent of the pulse window one has zoomed into.
	// The intg variation of stdevBase takes into account the integration of the pulses
	int i=0,count=0,newSize=last/intg;
	double tempStore=0,tempPhase=0;
	double *folda_use=(double*)malloc((last/intg)*sizeof(double));
	double *phase_use=(double*)malloc((last/intg)*sizeof(double));
	for(i=0;i<last;i++)					/*We are averaging out both the amplitude values and the phase values over the integration interval and then calculating the stdev over this integrated data*/ 
	{
		tempStore+=folda[i];
		tempPhase+=phase[i];
		count++;
		if(count%intg==0)
		{
			folda_use[count/intg-1]=tempStore/intg;
			phase_use[count/intg-1]=tempPhase/intg;
			tempStore=0;
			tempPhase=0;
		}
	}
	
	double baseline=meanBase(folda_use,phase_use,pulseStart,pulseEnd,newSize),sumSquare=0,stdev=0;
	tempPhase=0;
	count=0;
	
	for(i=0;i<newSize;i++)
	{
		tempPhase=phase_use[i];
		if (tempPhase>=pulseStart && tempPhase<=pulseEnd)
			continue;
		else	
		{
			sumSquare+=folda_use[i]*folda_use[i];
			count++;
		}
	}
	
	stdev=sqrt(sumSquare/count-baseline*baseline);
	free(folda_use);
	free(phase_use);
	return stdev;
}

double stdevBase(double *folda, double *phase, double pulseStart, double pulseEnd, int last) //THis has no integartion portion, calculating the stdev over the whole data
{
	int i=0,count=0;
	double baseline=meanBase(folda,phase,pulseStart,pulseEnd,last),sumSquare=0,tempPhase=0,stdev=0;
	
	for(i=0;i<last;i++)
	{
		tempPhase=phase[i];
		if (tempPhase>=pulseStart && tempPhase<=pulseEnd)
			continue;
		else	
		{
			sumSquare+=folda[i]*folda[i];
			count++;
		}
	}
	
	stdev=sqrt(sumSquare/count-baseline*baseline);
	return stdev;
}

double max(double *folda, int last)
{
	int i=0;
	double max=0;
	for(i=0;i<last;i++)
	{
		if (folda[i]>max)
			max=folda[i];
	}
	
	return max;
}

double min(double *folda, int last)
{
	int i=0;
	double min=folda[0];
	for(i=0;i<last;i++)
	{
		if (folda[i]<min)
			min=folda[i];
	}
	
	return min;
}

void smoothPulse(float *pulse, float *smoothPulse, int arraySize, int boxSize) //We are box car averagthe pulse data in order to calculate the smoothed pulse and then storing the smoothed pulse in smoothPulse
{
	int i,j;
	if (boxSize%2==0)
		boxSize++;
		
	int shiftSize=(boxSize-1)/2,countSamples=0;
	float tempVal=0;
	
	for(i=0;i<arraySize;i++)
	{
		countSamples=0;
		tempVal=0;
		for(j=i-shiftSize;j<=i+shiftSize;j++)
		{
			if(j<0 || j>=arraySize)
				continue;
			
			tempVal+=pulse[j];
			countSamples++;
		}
		smoothPulse[i]=tempVal/countSamples;
	}
}
void findWidths(float *gsbData, float *gwbData, int arraySize_gsb, int arraySize_gwb, float *widths, double t_gsb, double t_gwb, int intg_gsb, int intg_gwb)
{	//
	float lag=5, refWidth=0.5;			
	int autosize_gsb=(int)(lag/(t_gsb*intg_gsb)),autosize_gwb=(int)(lag/(t_gwb*intg_gwb)),i,j;
	float* autocorrGSB=(float*)malloc(autosize_gsb*sizeof(float));
	float* autocorrGWB=(float*)malloc(autosize_gwb*sizeof(float));
	float* lagsGSB=(float*)malloc(autosize_gsb*sizeof(float));
	float* lagsGWB=(float*)malloc(autosize_gwb*sizeof(float));
		
	float sum=0;
	float autoNorm=0;
	float slope=0;
	float refCorrGSB=0,refCorrGWB=0;
	
	for(i=0;i<autosize_gsb;i++)
	{
		sum=0;
		for(j=0;j<arraySize_gsb;j++)
		{
			if(j+i>=arraySize_gsb)
				sum+=gsbData[j]*gsbData[j+i-arraySize_gsb];
			else
				sum+=gsbData[j]*gsbData[j+i];
		}
		
		autocorrGSB[i]=sum;
		if(i==0)
			autoNorm=autocorrGSB[i];

		autocorrGSB[i]=autocorrGSB[i]/autoNorm;
		lagsGSB[i]=(float)(t_gsb*intg_gsb*i);
	}
	for(i=0;i<autosize_gwb;i++)
	{
		sum=0;
		for(j=0;j<arraySize_gwb;j++)
		{
			if(j+i>=arraySize_gwb)
				sum+=gwbData[j]*gwbData[j+i-arraySize_gwb];
			else
				sum+=gwbData[j]*gwbData[j+i];
		}
		
		autocorrGWB[i]=sum;
		if(i==0)
			autoNorm=autocorrGWB[i];

		autocorrGWB[i]=autocorrGWB[i]/autoNorm;
		lagsGWB[i]=(float)(t_gwb*intg_gwb*i);
	}
	printf("\n");
	for(i=0;i<autosize_gsb;i++)
	{
		if(autocorrGSB[i]<=autocorrGSB[i+1] && autocorrGSB[i+1]<=autocorrGSB[i+2] && autocorrGSB[i+2]<=autocorrGSB[i+3] && autocorrGSB[i]<=autocorrGSB[i-1] && autocorrGSB[i-1]<=autocorrGSB[i-2] && autocorrGSB[i-2]<=autocorrGSB[i-3])
		{
			printf("Minima found at GSB lag %f ms\n",lagsGSB[i]);
			refCorrGSB=autocorrGSB[0]-0.5*(autocorrGSB[0]-autocorrGSB[i]);
			for(j=0;j<i;j++)
			{
				if(autocorrGSB[j]<refCorrGSB)
				{
					slope=(autocorrGSB[j-1]-autocorrGSB[j])/(lagsGSB[j-1]-lagsGSB[j]);
					widths[0]=2*((refCorrGSB-autocorrGSB[j])/slope+lagsGSB[j]);
					break;
				}
			}
			break;
		}
		if(i==autosize_gsb-1)
			widths[0]=-1;
	}
	
	for(i=0;i<autosize_gwb;i++)
	{
		if(autocorrGWB[i]<=autocorrGWB[i+1] && autocorrGWB[i+1]<=autocorrGWB[i+2] && autocorrGWB[i+2]<=autocorrGWB[i+3] && autocorrGWB[i]<=autocorrGWB[i-1] && autocorrGWB[i-1]<=autocorrGWB[i-2] && autocorrGWB[i-2]<=autocorrGWB[i-3])
		{
			printf("Minima found at GWB lag %f ms\n",lagsGWB[i]);
			refCorrGWB=autocorrGWB[0]-0.5*(autocorrGWB[0]-autocorrGWB[i]);
			for(j=0;j<i;j++)
			{
				if(autocorrGWB[j]<refCorrGWB)
				{
					slope=(autocorrGWB[j-1]-autocorrGWB[j])/(lagsGWB[j-1]-lagsGWB[j]);
					widths[1]=2*((refCorrGWB-autocorrGWB[j])/slope+lagsGWB[j]);
					break;
				}
			}
			break;
		}
		if(i==autosize_gwb-1)
			widths[1]=-1;
	}
}

void findWidths_derSG(float *gsbData, float *gwbData, int arraySize_gsb, int arraySize_gwb, float *widths, double t_gsb, double t_gwb, int intg_gsb, int intg_gwb, double *sg_coeff, int sg_coeff_size)
{
	float lag=5, refWidth=0.5;
	int autosize_gsb=(int)(lag/(t_gsb*intg_gsb)),autosize_gwb=(int)(lag/(t_gwb*intg_gwb)),i,j;
	float* autocorrGSB=(float*)malloc(autosize_gsb*sizeof(float));
	float* autocorrGWB=(float*)malloc(autosize_gwb*sizeof(float));
	float* lagsGSB=(float*)malloc(autosize_gsb*sizeof(float));
	float* lagsGWB=(float*)malloc(autosize_gwb*sizeof(float));
		
	float sum=0;
	float autoNorm=0;
	float slope=0;
	float refCorrGSB=0,refCorrGWB=0;
	
	for(i=0;i<autosize_gsb;i++)
	{
		sum=0;
		for(j=0;j<arraySize_gsb;j++)
		{
			if(j+i>=arraySize_gsb)
				sum+=gsbData[j]*gsbData[j+i-arraySize_gsb];
			else
				sum+=gsbData[j]*gsbData[j+i];
		}
		
		autocorrGSB[i]=sum;
		if(i==0)
			autoNorm=autocorrGSB[i];

		autocorrGSB[i]=autocorrGSB[i]/autoNorm;
		lagsGSB[i]=(float)(t_gsb*intg_gsb*i);
	}
	for(i=0;i<autosize_gwb;i++)
	{
		sum=0;
		for(j=0;j<arraySize_gwb;j++)
		{
			if(j+i>=arraySize_gwb)
				sum+=gwbData[j]*gwbData[j+i-arraySize_gwb];
			else
				sum+=gwbData[j]*gwbData[j+i];
		}
		
		autocorrGWB[i]=sum;
		if(i==0)
			autoNorm=autocorrGWB[i];

		autocorrGWB[i]=autocorrGWB[i]/autoNorm;
		lagsGWB[i]=(float)(t_gwb*intg_gwb*i);
	}
	printf("\n");
	
	autocorrGWB[0]=autocorrGWB[1]+autocorrGWB[1]-autocorrGWB[2];
	autocorrGSB[0]=autocorrGSB[1]+autocorrGSB[1]-autocorrGSB[2];
	
	float *autocorrGSB_SGder=(float*)malloc(autosize_gsb*sizeof(float));
	float *autocorrGWB_SGder=(float*)malloc(autosize_gwb*sizeof(float));	
	
	float *autocorrGSB_SG=(float*)malloc(autosize_gsb*sizeof(float));
	float *autocorrGWB_SG=(float*)malloc(autosize_gwb*sizeof(float));
	
	for(i=0;i<autosize_gsb;i++)
	{
		autocorrGSB_SGder[i]=0;
		autocorrGSB_SG[i]=0;
		
		for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
		{
			if(j<0)
			{
				autocorrGSB_SGder[i]+=autocorrGSB[-j]*(j-i);
				
				if(j-i<0)
				{
					autocorrGSB_SG[i]+=autocorrGSB[-j]*sg_coeff[i-j];
				}
				else
				{
					autocorrGSB_SG[i]+=autocorrGSB[-j]*sg_coeff[j-i];
				}
			}	
			else
			{
				autocorrGSB_SGder[i]+=autocorrGSB[j]*(j-i);
				
				if(j-i<0)
				{
					autocorrGSB_SG[i]+=autocorrGSB[j]*sg_coeff[i-j];
				}
				else
				{
					autocorrGSB_SG[i]+=autocorrGSB[j]*sg_coeff[j-i];
				}
			}
		}
	}
	
	for(i=0;i<autosize_gwb;i++)
	{
		autocorrGWB_SGder[i]=0;
		autocorrGWB_SG[i]=0;
		
		for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
		{
			if(j<0)
			{
				autocorrGWB_SGder[i]+=autocorrGWB[-j]*(j-i);
				
				if(j-i<0)
				{
					autocorrGWB_SG[i]+=autocorrGWB[-j]*sg_coeff[i-j];
				}
				else
				{
					autocorrGWB_SG[i]+=autocorrGWB[-j]*sg_coeff[j-i];
				}
			}	
			else
			{
				autocorrGWB_SGder[i]+=autocorrGWB[j]*(j-i);
				
				if(j-i<0)
				{
					autocorrGWB_SG[i]+=autocorrGWB[j]*sg_coeff[i-j];
				}
				else
				{
					autocorrGWB_SG[i]+=autocorrGWB[j]*sg_coeff[j-i];
				}
			}
		}
	}
	
	for(i=0;i<autosize_gsb;i++)
	{
		autocorrGSB[i]=autocorrGSB_SG[i];
	}
	for(i=0;i<autosize_gwb;i++)
	{
		autocorrGWB[i]=autocorrGWB_SG[i];
	}
	
	
	printf("\nFinding minima from Savitzky Golay smoothed derivative\n");
	
	float slope_der,minima,autocorrminVal;
	for(i=0;i<autosize_gsb-1;i++)
	{
		if(autocorrGSB_SGder[i]<0 && autocorrGSB_SGder[i+1]>0)
		{
			slope_der=(autocorrGSB_SGder[i+1]-autocorrGSB_SGder[i])/(lagsGSB[i+1]-lagsGSB[i]);
			minima=lagsGSB[i]-autocorrGSB_SGder[i]/slope_der;
			slope=(autocorrGSB[i+1]-autocorrGSB[i])/(lagsGSB[i+1]-lagsGSB[i]);
			autocorrminVal=autocorrGSB[i]+slope*(minima-lagsGSB[i]);
			printf("Minima found at GSB lag %f ms with ACF GSB = %f\n",minima,autocorrminVal);
			refCorrGSB=autocorrGSB[0]-0.5*(autocorrGSB[0]-autocorrminVal);
			for(j=0;j<i;j++)
			{
				if(autocorrGSB[j]<refCorrGSB && j>0)
				{
					slope=(autocorrGSB[j-1]-autocorrGSB[j])/(lagsGSB[j-1]-lagsGSB[j]);
					widths[0]=2*((refCorrGSB-autocorrGSB[j])/slope+lagsGSB[j]);
					break;
				}
			}
			break;
		}
		if(i==autosize_gsb-2)
			widths[0]=-1;
	}
	for(i=0;i<autosize_gwb-1;i++)
	{
		if(autocorrGWB_SGder[i]<0 && autocorrGWB_SGder[i+1]>0)
		{
			slope_der=(autocorrGWB_SGder[i+1]-autocorrGWB_SGder[i])/(lagsGWB[i+1]-lagsGWB[i]);
			minima=lagsGWB[i]-autocorrGWB_SGder[i]/slope_der;
			slope=(autocorrGWB[i+1]-autocorrGWB[i])/(lagsGWB[i+1]-lagsGWB[i]);
			autocorrminVal=autocorrGWB[i]+slope*(minima-lagsGWB[i]);
			printf("Minima found at GWB lag %f ms with ACF GWB = %f\n",minima,autocorrminVal);
			refCorrGWB=autocorrGWB[0]-0.5*(autocorrGWB[0]-autocorrminVal);
			for(j=0;j<i;j++)
			{
				if(autocorrGWB[j]<refCorrGWB && j>0)
				{
					slope=(autocorrGWB[j-1]-autocorrGWB[j])/(lagsGWB[j-1]-lagsGWB[j]);
					widths[1]=2*((refCorrGWB-autocorrGWB[j])/slope+lagsGWB[j]);
					break;
				}
			}
			break;
		}
		if(i==autosize_gwb-2)
			widths[1]=-1;
	}
	
	free(autocorrGSB);
	free(autocorrGWB);
	free(autocorrGSB_SG);
	free(autocorrGWB_SG);
	free(autocorrGSB_SGder);
	free(autocorrGWB_SGder);
}
double filterStrong(double *folda, double *phase, double pulseStart, double pulseEnd, int last,int intg)	//THis function filters the pulse according to the threshol;d given based on the generic SNR calculation
{
	double baseline=meanBase(folda,phase,pulseStart,pulseEnd,last),stdev=stdevBase(folda,phase,pulseStart,pulseEnd,last),tempPhase=0,intSquare=0,intMean=0,SNR=0;
	int i,count=0;
	for(i=0;i<last;i++)
	{
		tempPhase=phase[i];
		if (tempPhase>=pulseStart && tempPhase<=pulseEnd)
		{
			intSquare+=folda[i];
			count++;
		}
		else
			continue;
	}
	intMean=intSquare/count;
	
	//SNR=((intMean-baseline)/stdev)*sqrt((count*1.0)/(intg*1.0));
	SNR=(intMean-baseline)/stdev;						//THis is the SNR calculation
	return SNR;
}

double filterStrongPeak(double *folda, double *phase, double pulseStart, double pulseEnd, int last) //The strong peak refers to the peakSNR calculation
{
	double baseline=meanBase(folda,phase,pulseStart,pulseEnd,last),stdev=stdevBase(folda,phase,pulseStart,pulseEnd,last),intMax=0,SNR=0;
	
	intMax=max(folda,last);
	
	SNR=(intMax-baseline)/stdev;
	return SNR;
}

void plottwoNoPar(char *file1, char *file2, int pulseNum)
{
	int numSamples_1=0,numSamples_2=0,i;
        float t,signal,maxP=0,minP=0;
        FILE *fileID;
        fileID=fopen(file1,"r");

        while (fscanf(fileID,"%f %f",&t,&signal)==2)
        {
                if (numSamples_1==1)
                        minP=signal;
                numSamples_1++;
        }
        rewind(fileID);
        float times_1[numSamples_1];
        float power_1[numSamples_1];

        for(i=0;i<numSamples_1;i++)
        {
                fscanf(fileID,"%f %f",&t,&signal);
                times_1[i]=t;
                power_1[i]=signal;
                if (power_1[i]>maxP)
                {
                        maxP=power_1[i];
                }
                else if (power_1[i]<minP && power_1[i]!=0)
                {
                        minP=power_1[i];
                }
        }
	
        cpgenv(times_1[0],times_1[numSamples_1-1],minP,maxP,0,0);
        char title[50];
	
	if(pulseNum<0)
		sprintf(title,"Folded profile at High Freq");
	else
		sprintf(title,"Pulse Number %d at High Freq",pulseNum);
	
	cpglab("Phase","Intensity",title);
        cpgline(numSamples_1,times_1,power_1);
		
	minP=maxP;
	fclose(fileID);
	maxP=0;
	fileID=fopen(file2,"r");
	
	while (fscanf(fileID,"%f %f",&t,&signal)==2)
        {
                if (numSamples_2==1)
                        minP=signal;
                numSamples_2++;
        }
        rewind(fileID);
	
        float times_2[numSamples_2];
        float power_2[numSamples_2];


	for(i=0;i<numSamples_2;i++)
        {
                fscanf(fileID,"%f %f",&t,&signal);
                times_2[i]=t;
                power_2[i]=signal;
                if (power_2[i]>maxP)
                {
                        maxP=power_2[i];
                }
                else if (power_2[i]<minP && power_2[i]!=0)
                {
                        minP=power_2[i];
                }
        }
	cpgpanl(1,1);
	cpgenv(times_2[0],times_2[numSamples_2-1],minP,maxP,0,0);

	if(pulseNum<0)
                sprintf(title,"Folded profile at Low Freq");
        else
                sprintf(title,"Pulse Number %d at Low Freq",pulseNum);

	cpglab("Phase","Intensity",title);

	cpgline(numSamples_2,times_2,power_2);
}

void plottwo(char *file1, char *file2, int pulseNum)		//This is used by the dual frequency part in order to plot individual pulses/folded pulses at the two different frequencies
{
	int numSamples_1=0,numSamples_2=0,i;
        float t,signal,maxP=0,minP=0;
        FILE *fileID;
        fileID=fopen(file1,"r");

        while (fscanf(fileID,"%f %f",&t,&signal)==2)
        {
                if (numSamples_1==1)
                        minP=signal;
                numSamples_1++;
        }
        rewind(fileID);
        float times_1[numSamples_1];
        float power_1[numSamples_1];

        for(i=0;i<numSamples_1;i++)
        {
                fscanf(fileID,"%f %f",&t,&signal);
                times_1[i]=t;
                power_1[i]=signal;
                if (power_1[i]>maxP)
                {
                        maxP=power_1[i];
                }
                else if (power_1[i]<minP && power_1[i]!=0)
                {
                        minP=power_1[i];
                }
        }
	
	cpgsubp(1,2);
        cpgenv(times_1[0],times_1[numSamples_1-1],minP,maxP,0,0);
        char title[50];
	
	if(pulseNum<0)
		sprintf(title,"Folded profile at High Freq");
	else
		sprintf(title,"Pulse Number %d at High Freq",pulseNum);
	
	cpglab("Phase","Intensity",title);
        cpgline(numSamples_1,times_1,power_1);
		
	minP=maxP;
	fclose(fileID);
	maxP=0;
	fileID=fopen(file2,"r");
	
	while (fscanf(fileID,"%f %f",&t,&signal)==2)
        {
                if (numSamples_2==1)
                        minP=signal;
                numSamples_2++;
        }
        rewind(fileID);
	
        float times_2[numSamples_2];
        float power_2[numSamples_2];


	for(i=0;i<numSamples_2;i++)
        {
                fscanf(fileID,"%f %f",&t,&signal);
                times_2[i]=t;
                power_2[i]=signal;
                if (power_2[i]>maxP)
                {
                        maxP=power_2[i];
                }
                else if (power_2[i]<minP && power_2[i]!=0)
                {
                        minP=power_2[i];
                }
        }
	cpgpanl(1,1);
	cpgenv(times_2[0],times_2[numSamples_2-1],minP,maxP,0,0);

	if(pulseNum<0)
                sprintf(title,"Folded profile at Low Freq");
        else
                sprintf(title,"Pulse Number %d at Low Freq",pulseNum);

	cpglab("Phase","Intensity",title);

	cpgline(numSamples_2,times_2,power_2);
}

void plotfour(char *file1, char *file2, char *file3, char *file4, int pulseNum)	//Used by dual frequency to plot the individual pulses at the two frequencies along with the folded profiles.
{
	int numSamples_1=0,numSamples_2=0,i;
        float t,signal,maxP=0,minP=0;
        FILE *fileID;
        fileID=fopen(file1,"r");

        while (fscanf(fileID,"%f %f",&t,&signal)==2)
        {
                if (numSamples_1==1)
                        minP=signal;
                numSamples_1++;
        }
        rewind(fileID);
        float times_1[numSamples_1];
        float power_1[numSamples_1];

        for(i=0;i<numSamples_1;i++)
        {
                fscanf(fileID,"%f %f",&t,&signal);
                times_1[i]=t;
                power_1[i]=signal;
                if (power_1[i]>maxP)
                {
                        maxP=power_1[i];
                }
                else if (power_1[i]<minP && power_1[i]!=0)
                {
                        minP=power_1[i];
                }
        }
	
	cpgsubp(1,4);
        cpgenv(times_1[0],times_1[numSamples_1-1],minP,maxP,0,0);
        char title[50];
	
	sprintf(title,"Pulse Number %d at High Freq",pulseNum);
	
	cpglab("Phase","Intensity",title);
        cpgline(numSamples_1,times_1,power_1);
		
	minP=maxP;
	fclose(fileID);
	maxP=0;
	fileID=fopen(file2,"r");
	
	for(i=0;i<numSamples_1;i++)
        {
                fscanf(fileID,"%f %f",&t,&signal);
                times_1[i]=t;
                power_1[i]=signal;
                if (power_1[i]>maxP)
                {
                        maxP=power_1[i];
                }
                else if (power_1[i]<minP && power_1[i]!=0)
                {
                        minP=power_1[i];
                }
        }
        
	cpgpanl(1,1);
        cpgenv(times_1[0],times_1[numSamples_1-1],minP,maxP,0,0);
	
	sprintf(title,"Folded profile at High Freq",pulseNum);
	
	cpglab("Phase","Intensity",title);
        cpgline(numSamples_1,times_1,power_1);
        
       	minP=maxP;
	fclose(fileID);
	maxP=0;
	fileID=fopen(file3,"r");
	
	while (fscanf(fileID,"%f %f",&t,&signal)==2)
        {
                if (numSamples_2==1)
                        minP=signal;
                numSamples_2++;
        }
        rewind(fileID);
	
        float times_2[numSamples_2];
        float power_2[numSamples_2];

	for(i=0;i<numSamples_2;i++)
        {
                fscanf(fileID,"%f %f",&t,&signal);
                times_2[i]=t;
                power_2[i]=signal;
                if (power_2[i]>maxP)
                {
                        maxP=power_2[i];
                }
                else if (power_2[i]<minP && power_2[i]!=0)
                {
                        minP=power_2[i];
                }
        }
	cpgpanl(1,2);
	cpgenv(times_2[0],times_2[numSamples_2-1],minP,maxP,0,0);

        sprintf(title,"Pulse Number %d at Low Freq",pulseNum);
	cpglab("Phase","Intensity",title);
	cpgline(numSamples_2,times_2,power_2);
	
	minP=maxP;
	fclose(fileID);
	maxP=0;
	fileID=fopen(file4,"r");
	
	for(i=0;i<numSamples_2;i++)
        {
                fscanf(fileID,"%f %f",&t,&signal);
                times_2[i]=t;
                power_2[i]=signal;
                if (power_2[i]>maxP)
                {
                        maxP=power_2[i];
                }
                else if (power_2[i]<minP && power_2[i]!=0)
                {
                        minP=power_2[i];
                }
        }
        
       	cpgpanl(1,3);
	cpgenv(times_2[0],times_2[numSamples_2-1],minP,maxP,0,0);

        sprintf(title,"Folded Profile at Low Freq",pulseNum);
	cpglab("Phase","Intensity",title);
	cpgline(numSamples_2,times_2,power_2);
}

void interpolate_float(float *shortSeries, float *longSeries, double shortRes, double longRes, float *interPolate, int numLongElements)
{	//interPolFlag = 1 means GWB series needs to be interpolated, and vice versa; LongRes is longer time series
	//As we donot have the same reolution for GSB and GWB, i cases of cross-correlation the data needs to be intrapolated so that both can be artificially made of the same resolution.
	float shortSlope;
	int i=0, shortPos;
	
	for(i=0;i<numLongElements;i++)
	{
		shortPos=(int)((i*longRes)/shortRes);
		shortSlope=(shortSeries[shortPos+1]-shortSeries[shortPos])/shortRes;
		interPolate[i]=shortSlope*(i*longRes-shortRes*shortPos)+shortSeries[shortPos];
	}
}

void interpolate(double *shortSeries, double *longSeries, double shortRes, double longRes, double *interPolate, int numLongElements)
{	//interPolFlag = 1 means GWB series needs to be interpolated, and vice versa; LongRes is longer time series
	double shortSlope;
	int i=0, shortPos;
	
	for(i=0;i<numLongElements;i++)
	{
		shortPos=(int)((i*longRes)/shortRes);
		shortSlope=(shortSeries[shortPos+1]-shortSeries[shortPos])/shortRes;
		interPolate[i]=shortSlope*(i*longRes-shortRes*shortPos)+shortSeries[shortPos];
	}
}

void detectQuasiPeriod(float *corrdata_gsb, float *corrdata_gwb, int corrsize_gsb, int corrsize_gwb, float *peaks, double t_gsb, double t_gwb, int intg_gsb, int intg_gwb, float threshFreq_gsb, float threshFreq_gwb)
{
	int minExpo_gsb=(int)floor((double)log((double)corrsize_gsb)/log(2))+1,minExpo_gwb=(int)floor((double)log((double)corrsize_gwb)/log(2))+1, fftsize_gsb=(int)pow((double)2,minExpo_gsb),fftsize_gwb=(int)pow((double)2,minExpo_gwb),i;
	float Nyqfreq_gsb=(float)1/(2*t_gsb*intg_gsb),Nyqfreq_gwb=(float)1/(2*t_gwb*intg_gwb),maxGSBfreq=0,maxGWBfreq=0,maxGSB=0,maxGWB=0;
	printf("\nGSB Sample size is %d. GSB FFT size is %d\nGWB Sample size is %d. GWB FFT size is %d",corrsize_gsb,fftsize_gsb,corrsize_gwb,fftsize_gwb);
	
	double *FFTUse_GSB=(double*)malloc(fftsize_gsb*sizeof(double));	
	double *FFTUse_GWB=(double*)malloc(fftsize_gwb*sizeof(double));	
	
	float rmsRes_gsb=0, rmsRes_gwb=0;
	
	for(i=0;i<corrsize_gsb;i++)
		rmsRes_gsb+=(corrdata_gsb[i]*corrdata_gsb[i]);
	for(i=0;i<corrsize_gwb;i++)
		rmsRes_gwb+=(corrdata_gwb[i]*corrdata_gwb[i]);		
		
	for(i=0;i<fftsize_gsb;i++)
	{
		if (i<corrsize_gsb)
			FFTUse_GSB[i]=corrdata_gsb[i];
		else
			FFTUse_GSB[i]=0;
	}
	
	for(i=0;i<fftsize_gwb;i++)
	{
		if (i<corrsize_gwb)
			FFTUse_GWB[i]=corrdata_gwb[i];
		else
			FFTUse_GWB[i]=0;
	}
	
	fftw_complex specGSB[fftsize_gsb/2+1];
	fftw_complex specGWB[fftsize_gwb/2+1];	
			
	fftw_plan planGSB, planGWB;
	planGSB=fftw_plan_dft_r2c_1d(fftsize_gsb,FFTUse_GSB,specGSB,FFTW_ESTIMATE);
	planGWB=fftw_plan_dft_r2c_1d(fftsize_gwb,FFTUse_GWB,specGWB,FFTW_ESTIMATE);
			
	fftw_execute_dft_r2c(planGSB,FFTUse_GSB,specGSB);
	fftw_execute_dft_r2c(planGWB,FFTUse_GWB,specGWB);	
	
	float *powGSB=(float*)malloc((fftsize_gsb/2+1)*sizeof(float));
	float *powGWB=(float*)malloc((fftsize_gwb/2+1)*sizeof(float));
	
	float *freqGSB=(float*)malloc((fftsize_gsb/2+1)*sizeof(float));
	float *freqGWB=(float*)malloc((fftsize_gwb/2+1)*sizeof(float));
	
	for(i=0;i<fftsize_gsb/2+1;i++)
	{
		powGSB[i]=(float)(specGSB[i][0]*specGSB[i][0]+specGSB[i][1]*specGSB[i][1])/(rmsRes_gsb*fftsize_gsb/2);
		freqGSB[i]=2.0*(float)i/fftsize_gsb*Nyqfreq_gsb;		
	}		
	for(i=0;i<fftsize_gwb/2+1;i++)
	{
		powGWB[i]=(float)(specGWB[i][0]*specGWB[i][0]+specGWB[i][1]*specGWB[i][1])/(rmsRes_gwb*fftsize_gwb/2);
		freqGWB[i]=2.0*(float)i/fftsize_gwb*Nyqfreq_gwb;		
	}		
	
	float minPow=powGSB[0],maxPow=0,totSqPowGSB=0,totSqPowGWB=0,totPowGSB=0,totPowGWB=0;
	int normRMSGSB=0,normRMSGWB=0;
	
	for(i=0;i<fftsize_gsb/2+1;i++)
	{
		if(powGSB[i]>maxPow)
			maxPow=powGSB[i];
		if(powGSB[i]<minPow)
			minPow=powGSB[i];
		
		if(freqGSB[i]>threshFreq_gsb)
		{
			totPowGSB+=powGSB[i];
			totSqPowGSB+=(powGSB[i]*powGSB[i]);
			normRMSGSB++;
		}
	}
	for(i=0;i<fftsize_gwb/2+1;i++)
	{
		if(powGWB[i]>maxPow)
			maxPow=powGWB[i];
		if(powGWB[i]<minPow)
			minPow=powGWB[i];
		
		if(freqGWB[i]>threshFreq_gwb)
		{
			totPowGWB+=powGWB[i];
			totSqPowGWB+=(powGWB[i]*powGWB[i]);
			normRMSGWB++;
		}
	}
	
	float rmsPowGSB=sqrt(totSqPowGSB/normRMSGSB-(totPowGSB/normRMSGSB)*(totPowGSB/normRMSGSB));
	float rmsPowGWB=sqrt(totSqPowGWB/normRMSGWB-(totPowGWB/normRMSGWB)*(totPowGWB/normRMSGWB));
	
	for(i=0;i<fftsize_gsb/2+1;i++)
	{
		if(powGSB[i]>powGSB[i-1] && powGSB[i-1]>powGSB[i-2] && powGSB[i]>powGSB[i+1] && powGSB[i+1]>powGSB[i+2] && freqGSB[i]>threshFreq_gsb)
		{
			if(powGSB[i]>powerThreshold*rmsPowGSB && powGSB[i]>maxGSB)
			{
				maxGSBfreq=freqGSB[i];
				maxGSB=powGSB[i];
			}
		}	
	}
	for(i=0;i<fftsize_gwb/2+1;i++)
	{
		if(powGWB[i]>powGWB[i-1] && powGWB[i-1]>powGWB[i-2] && powGWB[i]>powGWB[i+1] && powGWB[i+1]>powGWB[i+2] && freqGWB[i]>threshFreq_gwb)
		{
			if(powGWB[i]>powerThreshold*rmsPowGWB && powGWB[i]>maxGWB)
			{
				maxGWBfreq=freqGWB[i];
				maxGWB=powGWB[i];
			}
		}	
	}
	
	if(maxGSBfreq!=0)
		peaks[0]=1.0/maxGSBfreq;
	else
		peaks[0]=-1;
	
	if(maxGWBfreq!=0)
		peaks[1]=1.0/maxGWBfreq;
	else
		peaks[1]=-1;
}

void smoothSubPowerSpectrum(char *filename_GSB, char *filename_GWB, char *FFTName, double *phase_gsb, double *phase_gwb, double s, double *folda_gsb, double *folda_gwb, int last_gsb, int last_gwb,double pstart, double pend, double pulseStart, double pulseEnd, int intg_gsb, int intg_gwb, int pulseNum, double fold, double t_gsb, double t_gwb, double smoothDur_gsb, double smoothDur_gwb, double gsbDelay, double gwbDelay,float *readData_gsb,float *readData_gwb, double gsbWidth, double gwbWidth,int option)
//
{
	switch (option)
	{
	case 2:
	{
	FILE *data_GSB, *data_GWB, *FFTFile;
	data_GSB=fopen(filename_GSB,"rb");
	data_GWB=fopen(filename_GWB,"rb");
	FFTFile=fopen(FFTName,"w");
	int i,numshift_gsb, numshift_gwb,j;
	int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), skipnum_gwb=(int)((pulseNum-1)*(fold/t_gwb)), sampleDelay_gsb=gsbDelay/t_gsb*1000, sampleDelay_gwb=gwbDelay/t_gwb*1000;

	numshift_gsb=(int)(last_gsb*s);
	numshift_gwb=(int)(last_gwb*s);	
	fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);
	fseek(data_GWB,(skipnum_gwb+sampleDelay_gwb+numshift_gwb)*sizeof(float),SEEK_SET);	
	readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB);
	readFloatToDouble(last_gwb,readData_gwb,folda_gwb,data_GWB);
	
	int corrsize_gsb, corrsize_gwb, corrpos=0;
	corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
	corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *corrdata_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	float *phaseData_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *phaseData_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	
	double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),range_gsb=max(folda_gsb,last_gsb)-baseline_gsb,baseline_gwb=meanBase(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb),range_gwb=max(folda_gwb,last_gwb)-baseline_gwb;
	for(i=0;i<last_gsb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gsb && j<last_gsb;j++)
		{
			dpoint+=folda_gsb[j];
			ppoint+=phase_gsb[j];
		}
		if(j==last_gsb)
			break;
		ppoint=ppoint/intg_gsb;
		dpoint=dpoint/intg_gsb;
		dpoint=(dpoint-baseline_gsb)/range_gsb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
		{
			corrdata_gsb[corrpos]=dpoint;
			phaseData_gsb[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_gsb;
	}
	corrpos=0;
	for(i=0;i<last_gwb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gwb && j<last_gwb;j++)
		{
			dpoint+=folda_gwb[j];
			ppoint+=phase_gwb[j];
		}
		if(j==last_gwb)
			break;
		ppoint=ppoint/intg_gwb;
		dpoint=dpoint/intg_gwb;
		dpoint=(dpoint-baseline_gwb)/range_gwb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gwb)
		{
			corrdata_gwb[corrpos]=dpoint;
			phaseData_gwb[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_gwb;
	}	
	
	int smoothWindow_gsb, smoothWindow_gwb;	
	float* smoothGSB=(float*)malloc(corrsize_gsb*sizeof(float));
	float* smoothGWB=(float*)malloc(corrsize_gwb*sizeof(float));
	
	smoothWindow_gsb=(int)(smoothDur_gsb/(t_gsb*intg_gsb));
	smoothWindow_gwb=(int)(smoothDur_gwb/(t_gwb*intg_gwb));
	printf("\nSmoothing window size is %d bins for GSB; %d bins for GWB.\n",smoothWindow_gsb,smoothWindow_gwb);
	smoothPulse(corrdata_gsb,smoothGSB,corrsize_gsb,smoothWindow_gsb);
	smoothPulse(corrdata_gwb,smoothGWB,corrsize_gwb,smoothWindow_gwb);
	
	if(gsbWidth==0 || gwbWidth==0)
	{
		printf("\nError! Widths cannot be zero!\n");
		return;
	}
	
	float threshFreq_gsb=1.0/gsbWidth*subPulseFreqThreshold;
	float threshFreq_gwb=1.0/gwbWidth*subPulseFreqThreshold;
	printf("\nFrequency thresholds set at GSB %f kHz and GWB %f kHz\n",threshFreq_gsb,threshFreq_gwb);
	
	float maxSmooth=0,minSmooth=smoothGSB[0],maxRes=0,minRes=0,rmsRes_gsb=0,rmsRes_gwb=0;
	
	for(i=0;i<corrsize_gsb;i++)
	{
		corrdata_gsb[i]=corrdata_gsb[i]-smoothGSB[i];
		rmsRes_gsb+=corrdata_gsb[i]*corrdata_gsb[i];
		
		if(smoothGSB[i]>maxSmooth)
			maxSmooth=smoothGSB[i];

		if(smoothGSB[i]<minSmooth)
			minSmooth=smoothGSB[i];
		
		if(corrdata_gsb[i]>maxRes)
			maxRes=corrdata_gsb[i];
			
		if(corrdata_gsb[i]<minRes)
			minRes=corrdata_gsb[i];
	}
	for(i=0;i<corrsize_gwb;i++)
	{
		corrdata_gwb[i]=corrdata_gwb[i]-smoothGWB[i];
		rmsRes_gwb+=corrdata_gwb[i]*corrdata_gwb[i];
		
		if(smoothGWB[i]>maxSmooth)
			maxSmooth=smoothGWB[i];

		if(smoothGWB[i]<minSmooth)
			minSmooth=smoothGWB[i];
		
		if(corrdata_gwb[i]>maxRes)
			maxRes=corrdata_gwb[i];
			
		if(corrdata_gwb[i]<minRes)
			minRes=corrdata_gwb[i];
	}
	
	int minExpo_gsb=(int)floor((double)log((double)corrsize_gsb)/log(2))+1,minExpo_gwb=(int)floor((double)log((double)corrsize_gwb)/log(2))+1, fftsize_gsb=(int)pow((double)2,minExpo_gsb),fftsize_gwb=(int)pow((double)2,minExpo_gwb);
	float Nyqfreq_gsb=(float)1/(2*t_gsb*intg_gsb),Nyqfreq_gwb=(float)1/(2*t_gwb*intg_gwb),maxGSBfreq=0,maxGWBfreq=0,maxGSB=0,maxGWB=0;
	printf("\nGSB Sample size is %d. GSB FFT size is %d\nGWB Sample size is %d. GWB FFT size is %d",corrsize_gsb,fftsize_gsb,corrsize_gwb,fftsize_gwb);
	
	double *FFTUse_GSB=(double*)malloc(fftsize_gsb*sizeof(double));	
	double *FFTUse_GWB=(double*)malloc(fftsize_gwb*sizeof(double));	
	
	for(i=0;i<fftsize_gsb;i++)
	{
		if (i<corrsize_gsb)
			FFTUse_GSB[i]=corrdata_gsb[i];
		else
			FFTUse_GSB[i]=0;
	}
	
	for(i=0;i<fftsize_gwb;i++)
	{
		if (i<corrsize_gwb)
			FFTUse_GWB[i]=corrdata_gwb[i];
		else
			FFTUse_GWB[i]=0;
	}
	
	fftw_complex specGSB[fftsize_gsb/2+1];
	fftw_complex specGWB[fftsize_gwb/2+1];	
			
	fftw_plan planGSB, planGWB;
	planGSB=fftw_plan_dft_r2c_1d(fftsize_gsb,FFTUse_GSB,specGSB,FFTW_ESTIMATE);
	planGWB=fftw_plan_dft_r2c_1d(fftsize_gwb,FFTUse_GWB,specGWB,FFTW_ESTIMATE);
			
	fftw_execute_dft_r2c(planGSB,FFTUse_GSB,specGSB);
	fftw_execute_dft_r2c(planGWB,FFTUse_GWB,specGWB);	
	
	float *powGSB=(float*)malloc((fftsize_gsb/2+1)*sizeof(float));
	float *powGWB=(float*)malloc((fftsize_gwb/2+1)*sizeof(float));
	
	float *freqGSB=(float*)malloc((fftsize_gsb/2+1)*sizeof(float));
	float *freqGWB=(float*)malloc((fftsize_gwb/2+1)*sizeof(float));
	
	for(i=0;i<fftsize_gsb/2+1;i++)
	{
		powGSB[i]=(float)(specGSB[i][0]*specGSB[i][0]+specGSB[i][1]*specGSB[i][1])/(rmsRes_gsb*fftsize_gsb/2);
		freqGSB[i]=2.0*(float)i/fftsize_gsb*Nyqfreq_gsb;		
		
		fprintf(FFTFile,"%f\t%f\n",freqGSB[i],powGSB[i]);
	}		
	for(i=0;i<fftsize_gwb/2+1;i++)
	{
		powGWB[i]=(float)(specGWB[i][0]*specGWB[i][0]+specGWB[i][1]*specGWB[i][1])/(rmsRes_gwb*fftsize_gwb/2);
		freqGWB[i]=2.0*(float)i/fftsize_gwb*Nyqfreq_gwb;		
		
		fprintf(FFTFile,"%f\t%f\n",freqGWB[i],powGWB[i]);
	}		
	
	fclose(FFTFile);
	
	float minPow=powGSB[0],maxPow=0,totSqPowGSB=0,totSqPowGWB=0,totPowGSB=0,totPowGWB=0;
	int normRMSGSB=0,normRMSGWB=0;
	
	for(i=0;i<fftsize_gsb/2+1;i++)
	{
		if(powGSB[i]>maxPow)
			maxPow=powGSB[i];
		if(powGSB[i]<minPow)
			minPow=powGSB[i];
		
		if(freqGSB[i]>threshFreq_gsb)
		{
			totPowGSB+=powGSB[i];
			totSqPowGSB+=(powGSB[i]*powGSB[i]);
			normRMSGSB++;
		}
	}
	for(i=0;i<fftsize_gwb/2+1;i++)
	{
		if(powGWB[i]>maxPow)
			maxPow=powGWB[i];
		if(powGWB[i]<minPow)
			minPow=powGWB[i];
		
		if(freqGWB[i]>threshFreq_gwb)
		{
			totPowGWB+=powGWB[i];
			totSqPowGWB+=(powGWB[i]*powGWB[i]);
			normRMSGWB++;
		}
	}
	
	float rmsPowGSB=sqrt(totSqPowGSB/normRMSGSB-(totPowGSB/normRMSGSB)*(totPowGSB/normRMSGSB));
	float rmsPowGWB=sqrt(totSqPowGWB/normRMSGWB-(totPowGWB/normRMSGWB)*(totPowGWB/normRMSGWB));
	
	printf("\nRMSPOW_GSB=%f\tRMSPOW_GWB=%f\n",rmsPowGSB,rmsPowGWB);
	
	for(i=0;i<fftsize_gsb/2+1;i++)
	{
		if(powGSB[i]>powGSB[i-1] && powGSB[i-1]>powGSB[i-2] && powGSB[i]>powGSB[i+1] && powGSB[i+1]>powGSB[i+2] && freqGSB[i]>threshFreq_gsb)
		{
			if(powGSB[i]>powerThreshold*rmsPowGSB && powGSB[i]>maxGSB)
			{
				maxGSBfreq=freqGSB[i];
				maxGSB=powGSB[i];
			}
		}
		freqGSB[i]=(float)log10(2.0*(double)i/fftsize_gsb*Nyqfreq_gsb);		
	}
	for(i=0;i<fftsize_gwb/2+1;i++)
	{
		if(powGWB[i]>powGWB[i-1] && powGWB[i-1]>powGWB[i-2] && powGWB[i]>powGWB[i+1] && powGWB[i+1]>powGWB[i+2] && freqGWB[i]>threshFreq_gwb)
		{
			if(powGWB[i]>powerThreshold*rmsPowGWB && powGWB[i]>maxGWB)
			{
				maxGWBfreq=freqGWB[i];
				maxGWB=powGWB[i];
			}
		}
		freqGWB[i]=(float)log10(2.0*(double)i/fftsize_gwb*Nyqfreq_gwb);		
	}
	
	if(maxGSBfreq!=0)
		printf("\nPeriod GSB found at %f ms",1.0/maxGSBfreq);
	else
		printf("\nPeriod GSB not found");
	if(maxGWBfreq!=0)
		printf("\nPeriod GWB found at %f ms",1.0/maxGWBfreq);
	else
		printf("\nPeriod GWB not found");
			
	
	cpgsubp(1,3);
        cpgsci(1);
        cpgenv(pstart,pend,minRes,maxRes,0,0);
       	char title[50];
	sprintf(title,"Residuals for pulse number %d",pulseNum);
	cpglab("Phase","Residual Intensity",title);
	cpgsci(2);
        cpgline(corrsize_gsb,phaseData_gsb,corrdata_gsb);
        cpgsci(3);
        cpgline(corrsize_gwb,phaseData_gwb,corrdata_gwb);

        cpgpanl(1,1);
        cpgsci(1);
        cpgenv(pstart,pend,minSmooth,maxSmooth,0,0);
	sprintf(title,"Smoothed pulse for pulse number %d",pulseNum);
	cpglab("Phase","Intensity",title);
	cpgsci(2);
        cpgline(corrsize_gsb,phaseData_gsb,smoothGSB);
        cpgsci(3);
        cpgline(corrsize_gwb,phaseData_gwb,smoothGWB);
        
        cpgpanl(1,2);
        cpgsci(1);
        cpgenv(-1,freqGSB[fftsize_gsb/2],minPow,maxPow,0,10);
        sprintf(title,"Power spectrum of residuals for pulse number %d",pulseNum); 
	cpglab("Frequency (kHz)","Normalised Power",title);
	cpgsci(2);
        cpgline(fftsize_gsb/2,freqGSB,powGSB);
        cpgsci(3);
        cpgline(fftsize_gwb/2,freqGWB,powGWB);
        
	free(phaseData_gsb);
	free(phaseData_gwb);
	free(corrdata_gsb);
	free(corrdata_gwb);
	free(freqGSB);
	free(powGSB);
	free(freqGWB);
	free(powGWB);
	break;
	}
	case 1:
	{
	FILE *data_GSB, *FFTFile;
	data_GSB=fopen(filename_GSB,"rb");
	FFTFile=fopen(FFTName,"w");
	int i,numshift_gsb,j;
	int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), sampleDelay_gsb=gsbDelay/t_gsb*1000;

	numshift_gsb=(int)(last_gsb*s);
	fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);	
	readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB);
	
	int corrsize_gsb, corrpos=0;
	corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *phaseData_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	
	double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),range_gsb=max(folda_gsb,last_gsb)-baseline_gsb;
	for(i=0;i<last_gsb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gsb && j<last_gsb;j++)
		{
			dpoint+=folda_gsb[j];
			ppoint+=phase_gsb[j];
		}
		if(j==last_gsb)
			break;
		ppoint=ppoint/intg_gsb;
		dpoint=dpoint/intg_gsb;
		dpoint=(dpoint-baseline_gsb)/range_gsb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
		{
			corrdata_gsb[corrpos]=dpoint;
			phaseData_gsb[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_gsb;
	}
	corrpos=0;
	
	int smoothWindow_gsb;	
	float* smoothGSB=(float*)malloc(corrsize_gsb*sizeof(float));
	
	smoothWindow_gsb=(int)(smoothDur_gsb/(t_gsb*intg_gsb));
	printf("\nSmoothing window size is %d bins \n",smoothWindow_gsb);
	smoothPulse(corrdata_gsb,smoothGSB,corrsize_gsb,smoothWindow_gsb);
	
	if(gsbWidth==0)
	{
		printf("\nError! Widths cannot be zero!\n");
		return;
	}
	
	float threshFreq_gsb=1.0/gsbWidth*subPulseFreqThreshold;
	printf("\nFrequency thresholds set at GSB %f kHz \n",threshFreq_gsb);
	
	float maxSmooth=0,minSmooth=smoothGSB[0],maxRes=0,minRes=0,rmsRes_gsb=0;
	
	for(i=0;i<corrsize_gsb;i++)
	{
		corrdata_gsb[i]=corrdata_gsb[i]-smoothGSB[i];
		rmsRes_gsb+=corrdata_gsb[i]*corrdata_gsb[i];
		
		if(smoothGSB[i]>maxSmooth)
			maxSmooth=smoothGSB[i];

		if(smoothGSB[i]<minSmooth)
			minSmooth=smoothGSB[i];
		
		if(corrdata_gsb[i]>maxRes)
			maxRes=corrdata_gsb[i];
			
		if(corrdata_gsb[i]<minRes)
			minRes=corrdata_gsb[i];
	}
	
	int minExpo_gsb=(int)floor((double)log((double)corrsize_gsb)/log(2))+1, fftsize_gsb=(int)pow((double)2,minExpo_gsb);
	float Nyqfreq_gsb=(float)1/(2*t_gsb*intg_gsb),maxGSBfreq=0,maxGSB=0;
	printf("\nSample size is %d. FFT size is %d\n",corrsize_gsb,fftsize_gsb);
	
	double *FFTUse_GSB=(double*)malloc(fftsize_gsb*sizeof(double));		
	
	for(i=0;i<fftsize_gsb;i++)
	{
		if (i<corrsize_gsb)
			FFTUse_GSB[i]=corrdata_gsb[i];
		else
			FFTUse_GSB[i]=0;
	}
	
	fftw_complex specGSB[fftsize_gsb/2+1];
			
	fftw_plan planGSB;
	planGSB=fftw_plan_dft_r2c_1d(fftsize_gsb,FFTUse_GSB,specGSB,FFTW_ESTIMATE);
			
	fftw_execute_dft_r2c(planGSB,FFTUse_GSB,specGSB);	
	
	float *powGSB=(float*)malloc((fftsize_gsb/2+1)*sizeof(float));
	
	float *freqGSB=(float*)malloc((fftsize_gsb/2+1)*sizeof(float));
	
	for(i=0;i<fftsize_gsb/2+1;i++)
	{
		powGSB[i]=(float)(specGSB[i][0]*specGSB[i][0]+specGSB[i][1]*specGSB[i][1])/(rmsRes_gsb*fftsize_gsb/2);
		freqGSB[i]=2.0*(float)i/fftsize_gsb*Nyqfreq_gsb;		
		
		fprintf(FFTFile,"%f\t%f\n",freqGSB[i],powGSB[i]);
	}			
	
	fclose(FFTFile);
	
	float minPow=powGSB[0],maxPow=0,totSqPowGSB=0,totPowGSB=0;
	int normRMSGSB=0;
	
	for(i=0;i<fftsize_gsb/2+1;i++)
	{
		if(powGSB[i]>maxPow)
			maxPow=powGSB[i];
		if(powGSB[i]<minPow)
			minPow=powGSB[i];
		
		if(freqGSB[i]>threshFreq_gsb)
		{
			totPowGSB+=powGSB[i];
			totSqPowGSB+=(powGSB[i]*powGSB[i]);
			normRMSGSB++;
		}
	}
	
	float rmsPowGSB=sqrt(totSqPowGSB/normRMSGSB-(totPowGSB/normRMSGSB)*(totPowGSB/normRMSGSB));
	
	printf("\nRMSPOW_GSB=%f\n",rmsPowGSB);
	
	for(i=0;i<fftsize_gsb/2+1;i++)
	{
		if(powGSB[i]>powGSB[i-1] && powGSB[i-1]>powGSB[i-2] && powGSB[i]>powGSB[i+1] && powGSB[i+1]>powGSB[i+2] && freqGSB[i]>threshFreq_gsb)
		{
			if(powGSB[i]>powerThreshold*rmsPowGSB && powGSB[i]>maxGSB)
			{
				maxGSBfreq=freqGSB[i];
				maxGSB=powGSB[i];
			}
		}
		freqGSB[i]=(float)log10(2.0*(double)i/fftsize_gsb*Nyqfreq_gsb);		
	}
	
	if(maxGSBfreq!=0)
		printf("\nPeriod GSB found at %f ms",1.0/maxGSBfreq);
	else
		printf("\nPeriod GSB not found");
			
	
	cpgsubp(1,3);
        cpgsci(1);
        cpgenv(pstart,pend,minRes,maxRes,0,0);
       	char title[50];
	sprintf(title,"Residuals for pulse number %d",pulseNum);
	cpglab("Phase","Residual Intensity",title);
	cpgsci(2);
        cpgline(corrsize_gsb,phaseData_gsb,corrdata_gsb);

        cpgpanl(1,1);
        cpgsci(1);
        cpgenv(pstart,pend,minSmooth,maxSmooth,0,0);
	sprintf(title,"Smoothed pulse for pulse number %d",pulseNum);
	cpglab("Phase","Intensity",title);
	cpgsci(2);
        cpgline(corrsize_gsb,phaseData_gsb,smoothGSB);
        
        cpgpanl(1,2);
        cpgsci(1);
        cpgenv(-1,freqGSB[fftsize_gsb/2],minPow,maxPow,0,10);
        sprintf(title,"Power spectrum of residuals for pulse number %d",pulseNum); 
	cpglab("Frequency (kHz)","Normalised Power",title);
	cpgsci(2);
        cpgline(fftsize_gsb/2,freqGSB,powGSB);
        
	free(phaseData_gsb);
	free(corrdata_gsb);
	free(freqGSB);
	free(powGSB);
	}
	}
}
void giveRelStrengthSmoothed(float *corrdata_gsb, float *corrdata_gwb, float *smoothGSB, float *smoothGWB, int corrsize_gsb, int corrsize_gwb, float *relArray)
{
	int i,j;
	float totPowGSB=0, totPowGWB=0;
	float totSmoothGSB=0, totSmoothGWB=0;
	
	for(i=0;i<corrsize_gsb;i++)
	{	
		totPowGSB+=corrdata_gsb[i]*corrdata_gsb[i];
		totSmoothGSB+=smoothGSB[i]*smoothGSB[i];
	}
	for(i=0;i<corrsize_gwb;i++)
	{	
		totPowGWB+=corrdata_gwb[i]*corrdata_gwb[i];
		totSmoothGWB+=smoothGWB[i]*smoothGWB[i];
	}
	
	relArray[0]=totPowGSB/totSmoothGSB;
	relArray[1]=totPowGWB/totSmoothGWB;
	relArray[2]=totSmoothGSB;
	relArray[3]=totSmoothGWB;
	
	printf("\nRatios are GSB = %f \t GWB = %f\n",relArray[0],relArray[1]);
}

void smoothSubSGACF_intrapolated(char *filename_GSB, char *filename_GWB, char *autoName, char *SG_coeffs, double *phase_gsb, double *phase_gwb, double s, double *folda_gsb, double *folda_gwb, int last_gsb, int last_gwb,double pstart, double pend, double pulseStart, double pulseEnd, int intg_gsb, int intg_gwb, int pulseNum, double fold, double t_gsb, double t_gwb, double lag, double smoothDur_gsb, double smoothDur_gwb, double gsbDelay, double gwbDelay,float *readData_gsb,float *readData_gwb, int sg_coeff_size,int option,char *outfile)
{
	//Savitzky Golay filter coefficients
			FILE *out=fopen(outfile,"w");
			FILE *SG_coeff_file=fopen(SG_coeffs,"r");
			int sg_read[13];
			double sg_coeff[13];
			int sg_order, sg_norm,i;

			if(SG_coeff_file==NULL)
			{
				printf("\nCoefficients file not found!");
				return;
			}	

			while(1)
			{
				if(fscanf(SG_coeff_file,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",&sg_order,&sg_norm,sg_read,sg_read+1,sg_read+2,sg_read+3,sg_read+4, sg_read+5, sg_read+6, sg_read+7, sg_read+8,sg_read+9, sg_read+10,sg_read+11,sg_read+12)==15)
				{
					if(sg_order==2*sg_coeff_size+1)
					{
						printf("\nFound coefficients for smoothing window of %d points with normalisation of %d\n. Coefficients are:",sg_order,sg_norm);
						for(i=0;i<13;i++)
						{
							sg_coeff[i]=(double)sg_read[i]/sg_norm;
							printf("%f  ",sg_coeff[i]);
						}
						printf("\n");
						break;
					}
					else
						continue;
				}
				else
				{
					printf("\nsg Coeffs not found!\n");
					return;
				}

			}
			fclose(SG_coeff_file);	
	
			FILE *data_GSB, *data_GWB, *autoFile;
			data_GSB=fopen(filename_GSB,"rb");
			data_GWB=fopen(filename_GWB,"rb");
			autoFile=fopen(autoName,"w");
			int numshift_gsb, numshift_gwb,j;
			int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), skipnum_gwb=(int)((pulseNum-1)*(fold/t_gwb)), sampleDelay_gsb=gsbDelay/t_gsb*1000, sampleDelay_gwb=gwbDelay/t_gwb*1000;

			numshift_gsb=(int)(last_gsb*s);
			numshift_gwb=(int)(last_gwb*s);	
			fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);
			fseek(data_GWB,(skipnum_gwb+sampleDelay_gwb+numshift_gwb)*sizeof(float),SEEK_SET);	
			readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB);
			readFloatToDouble(last_gwb,readData_gwb,folda_gwb,data_GWB);
	
			int corrsize_gsb, corrsize_gwb, corrpos=0;
			corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
			corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
			float dpoint=0,ppoint=0;
			float *corrdata_gsb1=(float*)malloc(corrsize_gsb*sizeof(float));
			float *corrdata_gwb1=(float*)malloc(corrsize_gwb*sizeof(float));
			float *phaseData_gsb1=(float*)malloc(corrsize_gsb*sizeof(float));
			float *phaseData_gwb1=(float*)malloc(corrsize_gwb*sizeof(float));
	
			double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),range_gsb=max(folda_gsb,last_gsb)-baseline_gsb,baseline_gwb=meanBase(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb),range_gwb=max(folda_gwb,last_gwb)-baseline_gwb;
			for(i=0;i<last_gsb;)
			{
				dpoint=0;
				ppoint=0;
				for(j=i;j<i+intg_gsb && j<last_gsb;j++)
				{
					dpoint+=folda_gsb[j];
					ppoint+=phase_gsb[j];
				}
				if(j==last_gsb)
					break;
				ppoint=ppoint/intg_gsb;
				dpoint=dpoint/intg_gsb;
				dpoint=(dpoint-baseline_gsb)/range_gsb;
		
				if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
				{
					corrdata_gsb1[corrpos]=dpoint;
					phaseData_gsb1[corrpos]=ppoint;
					corrpos++;
				}
				i+=intg_gsb;
			}
			corrpos=0;
			for(i=0;i<last_gwb;)
			{
				dpoint=0;
				ppoint=0;
				for(j=i;j<i+intg_gwb && j<last_gwb;j++)
				{
					dpoint+=folda_gwb[j];
					ppoint+=phase_gwb[j];
				}
				if(j==last_gwb)
					break;
				ppoint=ppoint/intg_gwb;
				dpoint=dpoint/intg_gwb;
				dpoint=(dpoint-baseline_gwb)/range_gwb;
		
				if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gwb)
				{
					corrdata_gwb1[corrpos]=dpoint;
					phaseData_gwb1[corrpos]=ppoint;
					corrpos++;
				}
				i+=intg_gwb;
			}
			corrpos=0;

			float *crossArray1, *crossArray2,*phaseArray1,*phaseArray2;
			float t_corr; int corrsize;
			if(corrsize_gsb>corrsize_gwb)
			{
				crossArray1=corrdata_gsb1;
				crossArray2=(float*)malloc(corrsize_gsb*sizeof(float));
				phaseArray1=phaseData_gsb1;
				phaseArray2=(float*)malloc(corrsize_gsb*sizeof(float));
				interpolate_float(corrdata_gwb1,corrdata_gsb1,t_gwb*intg_gwb,t_gsb*intg_gsb,crossArray2,corrsize_gsb);
				interpolate_float(phaseData_gwb1,phaseData_gsb1,t_gwb*intg_gwb,t_gsb*intg_gsb,phaseArray2,corrsize_gsb);
				t_gwb=t_gsb*intg_gsb;
				corrsize_gwb=corrsize_gsb;
				printf("Yo");
			}
			else
			{
				crossArray2=corrdata_gwb1;
				crossArray1=(float*)malloc(corrsize_gwb*sizeof(float));
				phaseArray2=phaseData_gwb1;
				phaseArray1=(float*)malloc(corrsize_gwb*sizeof(float));
				interpolate_float(corrdata_gsb1,corrdata_gwb1,t_gsb*intg_gsb,t_gwb*intg_gwb,crossArray1,corrsize_gwb);
				interpolate_float(phaseData_gsb1,phaseData_gwb1,t_gsb*intg_gsb,t_gwb*intg_gwb,phaseArray1,corrsize_gwb);
				t_gsb=t_gwb*intg_gwb;
				corrsize_gsb=corrsize_gwb;
			} 

			float *corrdata_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
			float *corrdata_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
			float *phaseData_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
			float *phaseData_gwb=(float*)malloc(corrsize_gwb*sizeof(float));

			for (i=0;i<corrsize_gsb;i++)
			{
				corrdata_gsb[i]=crossArray1[i];
				//printf("%lf\n",corrdata_gsb[i]);
				//printf("%lf\n",crossArray1[i]);
				//printf("%lf\n",corrdata_gsb1[i]);
			}
			
			for (i=0;i<corrsize_gwb;i++)
			{
				corrdata_gwb[i]=crossArray2[i];
			}

			for (i=0;i<corrsize_gsb;i++)
			{
				phaseData_gsb[i]=phaseArray1[i];
			}
			
			for (i=0;i<corrsize_gwb;i++)
			{
				phaseData_gwb[i]=phaseArray2[i];
			}		
	
			int smoothWindow_gsb, smoothWindow_gwb;	
			float* smoothGSB=(float*)malloc(corrsize_gsb*sizeof(float));
			float* smoothGWB=(float*)malloc(corrsize_gwb*sizeof(float));
	
			smoothWindow_gsb=(int)(smoothDur_gsb/(t_gsb*intg_gsb));
			smoothWindow_gwb=(int)(smoothDur_gwb/(t_gwb*intg_gwb));
			printf("\nSmoothing window size is %d bins for GSB; %d bins for GWB.\n",smoothWindow_gsb,smoothWindow_gwb);
			smoothPulse(corrdata_gsb,smoothGSB,corrsize_gsb,smoothWindow_gsb);
			smoothPulse(corrdata_gwb,smoothGWB,corrsize_gwb,smoothWindow_gwb);
	
			float maxSmooth=0,minSmooth=smoothGSB[0],maxRes=0,minRes=0,rmsRes_gsb=0,rmsRes_gwb=0;
	
			for(i=0;i<corrsize_gsb;i++)
			{
				corrdata_gsb[i]=corrdata_gsb[i]-smoothGSB[i];
				rmsRes_gsb+=corrdata_gsb[i]*corrdata_gsb[i];
					
				if(smoothGSB[i]>maxSmooth)
					maxSmooth=smoothGSB[i];
		
				if(smoothGSB[i]<minSmooth)
					minSmooth=smoothGSB[i];
		
				if(corrdata_gsb[i]>maxRes)
					maxRes=corrdata_gsb[i];
					
				if(corrdata_gsb[i]<minRes)
					minRes=corrdata_gsb[i];
			}
			printf("Power at GSB = %f \n",rmsRes_gsb);
			for(i=0;i<corrsize_gwb;i++)
			{
				corrdata_gwb[i]=corrdata_gwb[i]-smoothGWB[i];
				rmsRes_gwb+=corrdata_gwb[i]*corrdata_gwb[i];
		
				if(smoothGWB[i]>maxSmooth)
					maxSmooth=smoothGWB[i];

				if(smoothGWB[i]<minSmooth)
					minSmooth=smoothGWB[i];
		
				if(corrdata_gwb[i]>maxRes)
					maxRes=corrdata_gwb[i];
			
				if(corrdata_gwb[i]<minRes)
					minRes=corrdata_gwb[i];
			}
			printf("Power at GWB = %f \n",rmsRes_gwb);
		
			for(i=0;i<last_gsb;i++)
				folda_gsb[i]=(folda_gsb[i]-baseline_gsb)/range_gsb;
			for(i=0;i<last_gwb;i++)
				folda_gwb[i]=(folda_gwb[i]-baseline_gwb)/range_gwb;
		
			float stdev_gsb=stdevBaseIntg(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb,intg_gsb),stdev_gwb=stdevBaseIntg(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb,intg_gwb);
	
			int normRes_gsb=corrsize_gsb, normRes_gwb=corrsize_gwb;
			rmsRes_gsb/=(normRes_gsb*stdev_gsb*stdev_gsb);
			rmsRes_gwb/=(normRes_gwb*stdev_gwb*stdev_gwb);
	
			rmsRes_gsb=sqrt((double)rmsRes_gsb);
			rmsRes_gwb=sqrt((double)rmsRes_gwb);

			float maxDev_gsb=maxMod(corrdata_gsb,corrsize_gsb),maxDev_gwb=maxMod(corrdata_gwb,corrsize_gwb);
			float maxDevRMS_gsb=maxDev_gsb/stdev_gsb,maxDevRMS_gwb=maxDev_gwb/stdev_gwb;
	
			printf("\nRMS deviation of GSB = %f\tGWB = %f for STDEV_GSB=%f\tSTDEV_GWB=%f",rmsRes_gsb,rmsRes_gwb,stdev_gsb,stdev_gwb);
			printf("\nMax_RMS_Dev of GSB = %f\tGWB = %f",maxDevRMS_gsb,maxDevRMS_gwb);
			float autoWidths[2];
			findWidths_derSG(corrdata_gsb,corrdata_gwb,corrsize_gsb,corrsize_gwb,autoWidths,t_gsb,t_gwb,intg_gsb,intg_gwb,sg_coeff,sg_coeff_size);
			printf("\nWidths detected GSB = %f ms\tGWB = %f ms\n",autoWidths[0],autoWidths[1]);	
	
			int autosize_gsb=(int)(lag/(t_gsb*intg_gsb));
			int autosize_gwb=(int)(lag/(t_gwb*intg_gwb));
			float* autocorrGSB=(float*)malloc(autosize_gsb*sizeof(float));
			float* autocorrGWB=(float*)malloc(autosize_gwb*sizeof(float));
			float* lagsGSB=(float*)malloc(autosize_gsb*sizeof(float));
			float* lagsGWB=(float*)malloc(autosize_gwb*sizeof(float));	
			
			float maxAuto=0,minAuto=1;
		
			float sum=0;
			float autoNorm=0;
	
			for(i=0;i<autosize_gsb;i++)
			{
				sum=0;
				for(j=0;j<corrsize_gsb;j++)
				{
					if(j+i>=corrsize_gsb)
					{
						sum+=corrdata_gsb[j]*corrdata_gsb[j+i-corrsize_gsb];
					}
					else
					{
						sum+=corrdata_gsb[j]*corrdata_gsb[j+i];
					}
				}
		
				autocorrGSB[i]=sum;
		
				if(i==0)
					autoNorm=autocorrGSB[i];

				autocorrGSB[i]=autocorrGSB[i]/autoNorm;
		
				if(autocorrGSB[i]>maxAuto)
				{
					maxAuto=autocorrGSB[i];
				}
				if(autocorrGSB[i]<minAuto)
				{
					minAuto=autocorrGSB[i];
				}
				
				lagsGSB[i]=(t_gsb*intg_gsb*i);
				fprintf(autoFile,"%f\t%f\n",lagsGSB[i],autocorrGSB[i]);
			}
			for(i=0;i<autosize_gwb;i++)
			{
				sum=0;
				for(j=0;j<corrsize_gwb;j++)
				{
					if(j+i>=corrsize_gwb)
					{
						sum+=corrdata_gwb[j]*corrdata_gwb[j+i-corrsize_gwb];
					}
					else
					{
						sum+=corrdata_gwb[j]*corrdata_gwb[j+i];
					}
				}
		
				autocorrGWB[i]=sum;
		
				if(i==0)
					autoNorm=autocorrGWB[i];

				autocorrGWB[i]=autocorrGWB[i]/autoNorm;
		
				if(autocorrGWB[i]>maxAuto)
				{
					maxAuto=autocorrGWB[i];
				}
				if(autocorrGWB[i]<minAuto)
				{
					minAuto=autocorrGWB[i];
				}
		
				lagsGWB[i]=(t_gwb*intg_gwb*i);
				fprintf(autoFile,"%f\t%f\n",lagsGWB[i],autocorrGWB[i]);
			}

	
	
			float *autocorrGSB_SG=(float*)malloc(autosize_gsb*sizeof(float));
			float *autocorrGWB_SG=(float*)malloc(autosize_gwb*sizeof(float));
	
			for(i=0;i<autosize_gsb;i++)
			{
				autocorrGSB_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							autocorrGSB_SG[i]+=autocorrGSB[-j]*sg_coeff[i-j];
						}
						else
						{
							autocorrGSB_SG[i]+=autocorrGSB[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							autocorrGSB_SG[i]+=autocorrGSB[j]*sg_coeff[i-j];
						}
						else
						{
							autocorrGSB_SG[i]+=autocorrGSB[j]*sg_coeff[j-i];
						}
					}
				}
			}

	
	
			for(i=0;i<autosize_gwb;i++)
			{
				autocorrGWB_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							autocorrGWB_SG[i]+=autocorrGWB[-j]*sg_coeff[i-j];
						}
						else
						{
							autocorrGWB_SG[i]+=autocorrGWB[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							autocorrGWB_SG[i]+=autocorrGWB[j]*sg_coeff[i-j];
						}
						else
						{
							autocorrGWB_SG[i]+=autocorrGWB[j]*sg_coeff[j-i];
						}
					}
				}
			}
	
			for(i=0;i<autosize_gsb;i++)
			{
				autocorrGSB[i]=autocorrGSB_SG[i];
			}
			for(i=0;i<autosize_gwb;i++)
			{
				autocorrGWB[i]=autocorrGWB_SG[i];
			}

        		fclose(autoFile);
        		fclose(data_GSB);
        		fclose(data_GWB);
        
        		cpgsubp(1,3);

        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       			char title[50];
			sprintf(title,"Residuals for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Phase","Residual Intensity",title);
			cpgsci(2);
        		cpgline(corrsize_gsb,phaseData_gsb,corrdata_gsb);
        		cpgsci(3);
        		cpgline(corrsize_gwb,phaseData_gwb,corrdata_gwb);


        		cpgpanl(1,1);
        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
			sprintf(title,"Smoothed pulse for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Phase","Intensity of smoothed part",title);
			cpgsci(2);
        		cpgline(corrsize_gsb,phaseData_gsb,smoothGSB);
        		cpgsci(3);
        		cpgline(corrsize_gwb,phaseData_gwb,smoothGWB);


        
        		cpgpanl(1,2);
        		cpgsci(1);
        		cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        		sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Intensity. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Time lag (ms)","Autocorrelation amplitude",title);
			cpgsci(2);
        		cpgline(autosize_gsb,lagsGSB,autocorrGSB);
        		cpgsci(3);
        		cpgline(autosize_gwb,lagsGWB,autocorrGWB);

	//Vincross 
			/*FILE *crfile=fopen("crossresidual.txt","w");

			for (i=0 ; i<corrsize_gsb ; i++)
			{
				corrdata_gsb[i]=corrdata_gsb[i]*50;
				//fprintf(garbage,"%.10lf ",corrdata_gsb[i]);
			}
			//fprintf(garbage,"\n\n\n\n");

			for (i=0 ; i<corrsize_gwb ; i++)
			{
				corrdata_gwb[i]=corrdata_gwb[i]*50;
				//fprintf(garbage,"%.10lf ",corrdata_gwb[i]);
			}

	

			double *crossArray1, *crossArray2;
			double t_corr; int corrsize;
			if(corrsize_gsb>corrsize_gwb)
			{
				crossArray1=corrdata_gsb;
				crossArray2=(double*)malloc(corrsize_gsb*sizeof(double));
				interpolate(corrdata_gwb,corrdata_gsb,t_gwb*intg_gwb,t_gsb*intg_gsb,crossArray2,corrsize_gsb);
				t_corr=t_gsb*intg_gsb;
				corrsize=corrsize_gsb;
			}
			else
			{
				crossArray2=corrdata_gwb;
				crossArray1=(double*)malloc(corrsize_gwb*sizeof(double));
				interpolate(corrdata_gsb,corrdata_gwb,t_gsb*intg_gsb,t_gwb*intg_gwb,crossArray1,corrsize_gwb);
				t_corr=t_gwb*intg_gwb;
				corrsize=corrsize_gwb;
			} 

			int crosssize=2*(int)((lag*1.0)/(t_corr*1.0));
			double *crosscorr=(double*)malloc(crosssize*sizeof(double));
			double sum1=0;
			int c=0;
			for(i=-crosssize/2;i<crosssize/2;i++)
			{
				sum1=0;
				for(j=0;j<corrsize;j++)
				{
					if(j+i>=corrsize)
						sum1+=crossArray1[j]*crossArray2[j+i-corrsize];
					else if(j+i<0)
						sum1+=crossArray1[j]*crossArray2[j+i+corrsize];
					else
						sum1+=crossArray1[j]*crossArray2[j+i];	
				}
				crosscorr[i+crosssize/2]=sum1;
				sum1=0;
			}
	
			float lags[crosssize];
			float crossfloat[crosssize];
			float minCross=crosscorr[crosssize-1],maxCross=0;
			for(i=-crosssize/2;i<crosssize/2;i++)
        		{	
        		        fprintf(out,"%.10lf\t%.10lf\n",(double)(t_corr*i),crosscorr[i+crosssize/2]);
        		}


        		fclose(out);
			int maxindex=-crosssize/2;
			for(i=-crosssize/2;i<crosssize/2;i++)
			{
				lags[i+crosssize/2]=(float)(t_corr*i);
				crossfloat[i+crosssize/2]=(float)crosscorr[i+crosssize/2];
				if(crossfloat[i+crosssize/2]>maxCross)
				{
					maxCross=crossfloat[i+crosssize/2];
				}
				else if(crossfloat[i+crosssize/2]<minCross)
					minCross=crossfloat[i+crosssize/2];
					
			}

			for (i=0;i<crosssize;i++)
			{
				crossfloat[i]=(crossfloat[i]*1.0)/(maxCross*1.0);
			}

			float *crosscorr_SG=(float*)malloc(crosssize*sizeof(float));
			
			for(i=0;i<crosssize;i++)
			{
				crosscorr_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							crosscorr_SG[i]+=crosscorr[-j]*sg_coeff[i-j];
						}
						else
						{
							crosscorr_SG[i]+=crosscorr[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							crosscorr_SG[i]+=crosscorr[j]*sg_coeff[i-j];
						}
						else
						{
							crosscorr_SG[i]+=crosscorr[j]*sg_coeff[j-i];
						}
					}
				}
			}

			for(i=0;i<crosssize;i++)
			{
				crosscorr[i]=crosscorr_SG[i];
			}
			for (i=0;i<crosssize;i++)
			{
				fprintf(crfile,"%.10lf",crossfloat[i]);
			} 

			cpgpanl(1,3);
			cpgsci(1);
			cpgenv(lags[0],lags[crosssize-1],0,1,0,0);
			sprintf(title,"Cross-correlation of residual part (HF is held stationary and LF is moved from +lag to -lag) for pulse number %d",pulseNum);		
			cpglab("Time lag (ms)","Cross-correlation amplitude",title);
        		cpgline(crosssize,lags,crossfloat);


			
			for (i=0 ; i<corrsize_gsb ; i++)
			{
				corrdata_gsb[i]=corrdata_gsb[i]/50;
			}

			for (i=0 ; i<corrsize_gwb ; i++)
			{
				corrdata_gwb[i]=corrdata_gwb[i]/50;
			}

			int microphase=0;
			printf("\nDo you want to change the phase range of the micropulse? 1=yes, 0=no\n");
			scanf("%d",&microphase);
			while (microphase==1)
			{
				printf ("\nenter phase range\n");
				printf("Start phase:");
				scanf ("%lf",&pstart);
				printf("End phase:");
				scanf ("%lf",&pend);

				int corrsize_gwb1=(int)(last_gwb/intg_gwb*(pend-pstart));
				int corrsize_gsb1=(int)(last_gsb/intg_gsb*(pend-pstart));

				int pos=0;

				float dpoint=0,ppoint=0;
				float *corrdata_gsb1=(float*)malloc(corrsize_gsb1*sizeof(float));
				float *corrdata_gwb1=(float*)malloc(corrsize_gwb1*sizeof(float));
				float *phaseData_gsb1=(float*)malloc(corrsize_gsb1*sizeof(float));
				float *phaseData_gwb1=(float*)malloc(corrsize_gwb1*sizeof(float));
				
				for(i=0;i<corrsize_gsb;i++)
				{
					if (phaseData_gsb[i]>=pstart && phaseData_gsb[i]<=pend)
					{
						corrdata_gsb1[pos]=corrdata_gsb[i]/50;
						phaseData_gsb1[pos]=phaseData_gsb[i];
					}
				}
				corrpos=0;
				for(i=0;i<corrsize_gwb;i++)
				{
					if (phaseData_gwb[i]>=pstart && phaseData_gwb[i]<=pend)
					{
						corrdata_gwb1[pos]=corrdata_gwb[i]/50;
						phaseData_gwb1[pos]=phaseData_gwb[i];
					}
				}
				corrpos=0;
				for(i=0;i<corrsize_gsb1;i++)
				{
					if(corrdata_gsb1[i]>maxRes)
						maxRes=corrdata_gsb1[i];
						
					if(corrdata_gsb1[i]<minRes)
						minRes=corrdata_gsb1[i];
					printf("\n%lf\t%lf",phaseData_gsb1[i],corrdata_gsb[i]);
				}	
				
				for(i=0;i<corrsize_gwb1;i++)
				{
					if(corrdata_gwb1[i]>maxRes)
						maxRes=corrdata_gwb1[i];
				
					if(corrdata_gwb1[i]<minRes)
						minRes=corrdata_gwb1[i];
				}

				cpgeras;

				cpgsubp(1,4);

				cpgsci(1);
        			cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       				char title[50];
				sprintf(title,"Residuals for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Phase","Residual Intensity",title);
				cpgsci(2);
        			cpgline(corrsize_gsb1,phaseData_gsb1,corrdata_gsb1);
        			cpgsci(3);
        			cpgline(corrsize_gwb1,phaseData_gwb1,corrdata_gwb1);

				cpgpanl(1,1);
        			cpgsci(1);
        			cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
				sprintf(title,"Smoothed pulse for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Phase","Intensity of smoothed part",title);
				cpgsci(2);
        			cpgline(corrsize_gsb,phaseData_gsb,smoothGSB);
        			cpgsci(3);
        			cpgline(corrsize_gwb,phaseData_gwb,smoothGWB);
	
	
        	
        			cpgpanl(1,2);
        			cpgsci(1);
        			cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        			sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Intensity. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Time lag (ms)","Autocorrelation amplitude",title);
				cpgsci(2);
        			cpgline(autosize_gsb,lagsGSB,autocorrGSB);
        			cpgsci(3);
        			cpgline(autosize_gwb,lagsGWB,autocorrGWB);

				cpgpanl(1,3);
				cpgsci(1);
				cpgenv(lags[0],lags[crosssize-1],0,1,0,0);
				sprintf(title,"Cross-correlation of residual part (HF is held stationary and LF is moved from +lag to -lag) for pulse number %d",pulseNum);		
				cpglab("Time lag (ms)","Cross-correlation amplitude",title);
	        		cpgline(crosssize,lags,crossfloat);

				printf("\nDo you want to re-size the phase window? 1=yes, 0=no\n");
				scanf("%d",&microphase);
			}*/			
				
                                                  

	
	//Vincross

			


			free(phaseData_gsb);
			free(phaseData_gwb);
			free(corrdata_gsb);
			free(corrdata_gwb);
			free(lagsGSB);
			free(lagsGWB);
			free(autocorrGSB);
			free(autocorrGWB);
			free(autocorrGSB_SG);
			free(autocorrGWB_SG);
}

void smoothSubSGACF(char *filename_GSB, char *filename_GWB, char *autoName, char *SG_coeffs, double *phase_gsb, double *phase_gwb, double s, double *folda_gsb, double *folda_gwb, int last_gsb, int last_gwb,double pstart, double pend, double pulseStart, double pulseEnd, int intg_gsb, int intg_gwb, int pulseNum, double fold, double t_gsb, double t_gwb, double lag, double smoothDur_gsb, double smoothDur_gwb, double gsbDelay, double gwbDelay,float *readData_gsb,float *readData_gwb, int sg_coeff_size,int option,char *outfile)
{
	switch (option)
	{
		case 2:
		{ 
			//Savitzky Golay filter coefficients
			FILE *out=fopen(outfile,"w");
			FILE *SG_coeff_file=fopen(SG_coeffs,"r");
			int sg_read[13];
			double sg_coeff[13];
			int sg_order, sg_norm,i;

			if(SG_coeff_file==NULL)
			{
				printf("\nCoefficients file not found!");
				return;
			}	

			while(1)
			{
				if(fscanf(SG_coeff_file,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",&sg_order,&sg_norm,sg_read,sg_read+1,sg_read+2,sg_read+3,sg_read+4, sg_read+5, sg_read+6, sg_read+7, sg_read+8,sg_read+9, sg_read+10,sg_read+11,sg_read+12)==15)
				{
					if(sg_order==2*sg_coeff_size+1)
					{
						printf("\nFound coefficients for smoothing window of %d points with normalisation of %d\n. Coefficients are:",sg_order,sg_norm);
						for(i=0;i<13;i++)
						{
							sg_coeff[i]=(double)sg_read[i]/sg_norm;
							printf("%f  ",sg_coeff[i]);
						}
						printf("\n");
						break;
					}
					else
						continue;
				}
				else
				{
					printf("\nsg Coeffs not found!\n");
					return;
				}

			}
			fclose(SG_coeff_file);	
	
			FILE *data_GSB, *data_GWB, *autoFile;
			data_GSB=fopen(filename_GSB,"rb");
			data_GWB=fopen(filename_GWB,"rb");
			autoFile=fopen(autoName,"w");
			int numshift_gsb, numshift_gwb,j;
			int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), skipnum_gwb=(int)((pulseNum-1)*(fold/t_gwb)), sampleDelay_gsb=gsbDelay/t_gsb*1000, sampleDelay_gwb=gwbDelay/t_gwb*1000;

			numshift_gsb=(int)(last_gsb*s);
			numshift_gwb=(int)(last_gwb*s);	
			fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);
			fseek(data_GWB,(skipnum_gwb+sampleDelay_gwb+numshift_gwb)*sizeof(float),SEEK_SET);	
			readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB);
			readFloatToDouble(last_gwb,readData_gwb,folda_gwb,data_GWB);
	
			int corrsize_gsb, corrsize_gwb, corrpos=0;
			corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
			corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
			float dpoint=0,ppoint=0;
			float *corrdata_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
			float *corrdata_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
			float *phaseData_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
			float *phaseData_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	
			double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),range_gsb=max(folda_gsb,last_gsb)-baseline_gsb,baseline_gwb=meanBase(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb),range_gwb=max(folda_gwb,last_gwb)-baseline_gwb;
			for(i=0;i<last_gsb;)
			{
				dpoint=0;
				ppoint=0;
				for(j=i;j<i+intg_gsb && j<last_gsb;j++)
				{
					dpoint+=folda_gsb[j];
					ppoint+=phase_gsb[j];
				}
				if(j==last_gsb)
					break;
				ppoint=ppoint/intg_gsb;
				dpoint=dpoint/intg_gsb;
				dpoint=(dpoint-baseline_gsb)/range_gsb;
		
				if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
				{
					corrdata_gsb[corrpos]=dpoint;
					phaseData_gsb[corrpos]=ppoint;
					corrpos++;
				}
				i+=intg_gsb;
			}
			corrpos=0;
			for(i=0;i<last_gwb;)
			{
				dpoint=0;
				ppoint=0;
				for(j=i;j<i+intg_gwb && j<last_gwb;j++)
				{
					dpoint+=folda_gwb[j];
					ppoint+=phase_gwb[j];
				}
				if(j==last_gwb)
					break;
				ppoint=ppoint/intg_gwb;
				dpoint=dpoint/intg_gwb;
				dpoint=(dpoint-baseline_gwb)/range_gwb;
		
				if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gwb)
				{
					corrdata_gwb[corrpos]=dpoint;
					phaseData_gwb[corrpos]=ppoint;
					corrpos++;
				}
				i+=intg_gwb;
			}
			corrpos=0;	
	
			int smoothWindow_gsb, smoothWindow_gwb;	
			float* smoothGSB=(float*)malloc(corrsize_gsb*sizeof(float));
			float* smoothGWB=(float*)malloc(corrsize_gwb*sizeof(float));
	
			smoothWindow_gsb=(int)(smoothDur_gsb/(t_gsb*intg_gsb));
			smoothWindow_gwb=(int)(smoothDur_gwb/(t_gwb*intg_gwb));
			printf("\nSmoothing window size is %d bins for GSB; %d bins for GWB.\n",smoothWindow_gsb,smoothWindow_gwb);
			smoothPulse(corrdata_gsb,smoothGSB,corrsize_gsb,smoothWindow_gsb);
			smoothPulse(corrdata_gwb,smoothGWB,corrsize_gwb,smoothWindow_gwb);
	
			float maxSmooth=0,minSmooth=smoothGSB[0],maxRes=0,minRes=0,rmsRes_gsb=0,rmsRes_gwb=0;
	
			for(i=0;i<corrsize_gsb;i++)
			{
				corrdata_gsb[i]=corrdata_gsb[i]-smoothGSB[i];
				rmsRes_gsb+=corrdata_gsb[i]*corrdata_gsb[i];
					
				if(smoothGSB[i]>maxSmooth)
					maxSmooth=smoothGSB[i];
		
				if(smoothGSB[i]<minSmooth)
					minSmooth=smoothGSB[i];
		
				if(corrdata_gsb[i]>maxRes)
					maxRes=corrdata_gsb[i];
					
				if(corrdata_gsb[i]<minRes)
					minRes=corrdata_gsb[i];
			}
			printf("Power at GSB = %f \n",rmsRes_gsb);
			for(i=0;i<corrsize_gwb;i++)
			{
				corrdata_gwb[i]=corrdata_gwb[i]-smoothGWB[i];
				rmsRes_gwb+=corrdata_gwb[i]*corrdata_gwb[i];
		
				if(smoothGWB[i]>maxSmooth)
					maxSmooth=smoothGWB[i];

				if(smoothGWB[i]<minSmooth)
					minSmooth=smoothGWB[i];
		
				if(corrdata_gwb[i]>maxRes)
					maxRes=corrdata_gwb[i];
			
				if(corrdata_gwb[i]<minRes)
					minRes=corrdata_gwb[i];
			}
			printf("Power at GWB = %f \n",rmsRes_gwb);
		
			for(i=0;i<last_gsb;i++)
				folda_gsb[i]=(folda_gsb[i]-baseline_gsb)/range_gsb;
			for(i=0;i<last_gwb;i++)
				folda_gwb[i]=(folda_gwb[i]-baseline_gwb)/range_gwb;
		
			float stdev_gsb=stdevBaseIntg(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb,intg_gsb),stdev_gwb=stdevBaseIntg(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb,intg_gwb);
	
			int normRes_gsb=corrsize_gsb, normRes_gwb=corrsize_gwb;
			rmsRes_gsb/=(normRes_gsb*stdev_gsb*stdev_gsb);
			rmsRes_gwb/=(normRes_gwb*stdev_gwb*stdev_gwb);
	
			rmsRes_gsb=sqrt((double)rmsRes_gsb);
			rmsRes_gwb=sqrt((double)rmsRes_gwb);

			float maxDev_gsb=maxMod(corrdata_gsb,corrsize_gsb),maxDev_gwb=maxMod(corrdata_gwb,corrsize_gwb);
			float maxDevRMS_gsb=maxDev_gsb/stdev_gsb,maxDevRMS_gwb=maxDev_gwb/stdev_gwb;
	
			printf("\nRMS deviation of GSB = %f\tGWB = %f for STDEV_GSB=%f\tSTDEV_GWB=%f",rmsRes_gsb,rmsRes_gwb,stdev_gsb,stdev_gwb);
			printf("\nMax_RMS_Dev of GSB = %f\tGWB = %f",maxDevRMS_gsb,maxDevRMS_gwb);
			float autoWidths[2];
			findWidths_derSG(corrdata_gsb,corrdata_gwb,corrsize_gsb,corrsize_gwb,autoWidths,t_gsb,t_gwb,intg_gsb,intg_gwb,sg_coeff,sg_coeff_size);
			printf("\nWidths detected GSB = %f ms\tGWB = %f ms\n",autoWidths[0],autoWidths[1]);	
	
			int autosize_gsb=(int)(lag/(t_gsb*intg_gsb));
			int autosize_gwb=(int)(lag/(t_gwb*intg_gwb));
			float* autocorrGSB=(float*)malloc(autosize_gsb*sizeof(float));
			float* autocorrGWB=(float*)malloc(autosize_gwb*sizeof(float));
			float* lagsGSB=(float*)malloc(autosize_gsb*sizeof(float));
			float* lagsGWB=(float*)malloc(autosize_gwb*sizeof(float));	
			
			float maxAuto=0,minAuto=1;
		
			float sum=0;
			float autoNorm=0;
	
			for(i=0;i<autosize_gsb;i++)
			{
				sum=0;
				for(j=0;j<corrsize_gsb;j++)
				{
					if(j+i>=corrsize_gsb)
					{
						sum+=corrdata_gsb[j]*corrdata_gsb[j+i-corrsize_gsb];
					}
					else
					{
						sum+=corrdata_gsb[j]*corrdata_gsb[j+i];
					}
				}
		
				autocorrGSB[i]=sum;
		
				if(i==0)
					autoNorm=autocorrGSB[i];

				autocorrGSB[i]=autocorrGSB[i]/autoNorm;
		
				if(autocorrGSB[i]>maxAuto)
				{
					maxAuto=autocorrGSB[i];
				}
				if(autocorrGSB[i]<minAuto)
				{
					minAuto=autocorrGSB[i];
				}
				
				lagsGSB[i]=(t_gsb*intg_gsb*i);
				fprintf(autoFile,"%f\t%f\n",lagsGSB[i],autocorrGSB[i]);
			}
			for(i=0;i<autosize_gwb;i++)
			{
				sum=0;
				for(j=0;j<corrsize_gwb;j++)
				{
					if(j+i>=corrsize_gwb)
					{
						sum+=corrdata_gwb[j]*corrdata_gwb[j+i-corrsize_gwb];
					}
					else
					{
						sum+=corrdata_gwb[j]*corrdata_gwb[j+i];
					}
				}
		
				autocorrGWB[i]=sum;
		
				if(i==0)
					autoNorm=autocorrGWB[i];

				autocorrGWB[i]=autocorrGWB[i]/autoNorm;
		
				if(autocorrGWB[i]>maxAuto)
				{
					maxAuto=autocorrGWB[i];
				}
				if(autocorrGWB[i]<minAuto)
				{
					minAuto=autocorrGWB[i];
				}
		
				lagsGWB[i]=(t_gwb*intg_gwb*i);
				fprintf(autoFile,"%f\t%f\n",lagsGWB[i],autocorrGWB[i]);
			}

	
	
			float *autocorrGSB_SG=(float*)malloc(autosize_gsb*sizeof(float));
			float *autocorrGWB_SG=(float*)malloc(autosize_gwb*sizeof(float));
	
			for(i=0;i<autosize_gsb;i++)
			{
				autocorrGSB_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							autocorrGSB_SG[i]+=autocorrGSB[-j]*sg_coeff[i-j];
						}
						else
						{
							autocorrGSB_SG[i]+=autocorrGSB[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							autocorrGSB_SG[i]+=autocorrGSB[j]*sg_coeff[i-j];
						}
						else
						{
							autocorrGSB_SG[i]+=autocorrGSB[j]*sg_coeff[j-i];
						}
					}
				}
			}

	
	
			for(i=0;i<autosize_gwb;i++)
			{
				autocorrGWB_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							autocorrGWB_SG[i]+=autocorrGWB[-j]*sg_coeff[i-j];
						}
						else
						{
							autocorrGWB_SG[i]+=autocorrGWB[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							autocorrGWB_SG[i]+=autocorrGWB[j]*sg_coeff[i-j];
						}
						else
						{
							autocorrGWB_SG[i]+=autocorrGWB[j]*sg_coeff[j-i];
						}
					}
				}
			}
	
			for(i=0;i<autosize_gsb;i++)
			{
				autocorrGSB[i]=autocorrGSB_SG[i];
			}
			for(i=0;i<autosize_gwb;i++)
			{
				autocorrGWB[i]=autocorrGWB_SG[i];
			}

        		fclose(autoFile);
        		fclose(data_GSB);
        		fclose(data_GWB);
        
        		cpgsubp(1,4);

        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       			char title[50];
			sprintf(title,"Residuals for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Phase","Residual Intensity",title);
			cpgsci(2);
        		cpgline(corrsize_gsb,phaseData_gsb,corrdata_gsb);
        		cpgsci(3);
        		cpgline(corrsize_gwb,phaseData_gwb,corrdata_gwb);


        		cpgpanl(1,1);
        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
			sprintf(title,"Smoothed pulse for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Phase","Intensity of smoothed part",title);
			cpgsci(2);
        		cpgline(corrsize_gsb,phaseData_gsb,smoothGSB);
        		cpgsci(3);
        		cpgline(corrsize_gwb,phaseData_gwb,smoothGWB);


        
        		cpgpanl(1,2);
        		cpgsci(1);
        		cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        		sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Intensity. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Time lag (ms)","Autocorrelation amplitude",title);
			cpgsci(2);
        		cpgline(autosize_gsb,lagsGSB,autocorrGSB);
        		cpgsci(3);
        		cpgline(autosize_gwb,lagsGWB,autocorrGWB);

	//Vincross 
			FILE *crfile=fopen("crossresidual.txt","w");

			/*for (i=0 ; i<corrsize_gsb ; i++)
			{
				corrdata_gsb[i]=corrdata_gsb[i]*50;
				//fprintf(garbage,"%.10lf ",corrdata_gsb[i]);
			}
			//fprintf(garbage,"\n\n\n\n");

			for (i=0 ; i<corrsize_gwb ; i++)
			{
				corrdata_gwb[i]=corrdata_gwb[i]*50;
				//fprintf(garbage,"%.10lf ",corrdata_gwb[i]);
			}*/

	

			float *crossArray1, *crossArray2;
			float t_corr; int corrsize;
			if(corrsize_gsb>corrsize_gwb)
			{
				crossArray1=corrdata_gsb;
				crossArray2=(float*)malloc(corrsize_gsb*sizeof(float));
				interpolate_float(corrdata_gwb,corrdata_gsb,t_gwb*intg_gwb,t_gsb*intg_gsb,crossArray2,corrsize_gsb);
				t_corr=t_gsb*intg_gsb;
				corrsize=corrsize_gsb;
			}
			else
			{
				crossArray2=corrdata_gwb;
				crossArray1=(float*)malloc(corrsize_gwb*sizeof(float));
				interpolate_float(corrdata_gsb,corrdata_gwb,t_gsb*intg_gsb,t_gwb*intg_gwb,crossArray1,corrsize_gwb);
				t_corr=t_gwb*intg_gwb;
				corrsize=corrsize_gwb;
			} 

			int crosssize=2*(int)((lag*1.0)/(t_corr*1.0));
			float *crosscorr=(float*)malloc(crosssize*sizeof(float));
			float sum1=0;
			int c=0;
			for(i=-crosssize/2;i<crosssize/2;i++)
			{
				sum1=0;
				for(j=0;j<corrsize;j++)
				{
					if(j+i>=corrsize)
						sum1+=crossArray1[j]*crossArray2[j+i-corrsize];
					else if(j+i<0)
						sum1+=crossArray1[j]*crossArray2[j+i+corrsize];
					else
						sum1+=crossArray1[j]*crossArray2[j+i];	
				}
				crosscorr[i+crosssize/2]=sum1;
				sum1=0;
			}
	
			float lags[crosssize];
			float crossfloat[crosssize];
			float minCross=crosscorr[crosssize-1],maxCross=0;
			for(i=-crosssize/2;i<crosssize/2;i++)
        		{	
        		        fprintf(out,"%.10lf\t%.10lf\n",(float)(t_corr*i),crosscorr[i+crosssize/2]);
        		}


        		fclose(out);
			int maxindex=-crosssize/2;
			for(i=-crosssize/2;i<crosssize/2;i++)
			{
				lags[i+crosssize/2]=(float)(t_corr*i);
				crossfloat[i+crosssize/2]=(float)crosscorr[i+crosssize/2];
				if(crossfloat[i+crosssize/2]>maxCross)
				{
					maxCross=crossfloat[i+crosssize/2];
				}
				else if(crossfloat[i+crosssize/2]<minCross)
					minCross=crossfloat[i+crosssize/2];
					
			}

			for (i=0;i<crosssize;i++)
			{
				crossfloat[i]=(crossfloat[i]*1.0)/(maxCross*1.0);
			}

			float *crosscorr_SG=(float*)malloc(crosssize*sizeof(float));
			
			for(i=0;i<crosssize;i++)
			{
				crosscorr_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							crosscorr_SG[i]+=crosscorr[-j]*sg_coeff[i-j];
						}
						else
						{
							crosscorr_SG[i]+=crosscorr[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							crosscorr_SG[i]+=crosscorr[j]*sg_coeff[i-j];
						}
						else
						{
							crosscorr_SG[i]+=crosscorr[j]*sg_coeff[j-i];
						}
					}
				}
			}
	
			for(i=0;i<crosssize;i++)
			{
				crosscorr[i]=crosscorr_SG[i];
				
			}

			float maxpos=0;
			float maxval=0;

			for(i=0;i<crosssize;i++)
			{
				if(crossfloat[i]>maxval)
				{
					maxval=crossfloat[i];
					maxpos=lags[i];
				}
			}				
			
			printf("\nCross-correlation peak at %f ms",maxpos);			
			
			for (i=0;i<crosssize;i++)
			{
				fprintf(crfile,"%.10lf",crossfloat[i]);
			} 

			cpgpanl(1,3);
			cpgsci(1);
			cpgenv(lags[0],lags[crosssize-1],0,1,0,0);
			sprintf(title,"Cross-correlation of residual part (HF is held stationary and LF is moved from +lag to -lag) for pulse number %d",pulseNum);		
			cpglab("Time lag (ms)","Cross-correlation amplitude",title);
        		cpgline(crosssize,lags,crossfloat);

			
			
			/*for (i=0 ; i<corrsize_gsb ; i++)
			{
				corrdata_gsb[i]=corrdata_gsb[i]/50;
			}

			for (i=0 ; i<corrsize_gwb ; i++)
			{
				corrdata_gwb[i]=corrdata_gwb[i]/50;
			}

			/*int microphase=0;
			printf("\nDo you want to change the phase range of the micropulse? 1=yes, 0=no\n");
			scanf("%d",&microphase);
			while (microphase==1)
			{
				printf ("\nenter phase range\n");
				printf("Start phase:");
				scanf ("%lf",&pstart);
				printf("End phase:");
				scanf ("%lf",&pend);

				int corrsize_gwb1=(int)(last_gwb/intg_gwb*(pend-pstart));
				int corrsize_gsb1=(int)(last_gsb/intg_gsb*(pend-pstart));

				int pos=0;

				float dpoint=0,ppoint=0;
				float *corrdata_gsb1=(float*)malloc(corrsize_gsb1*sizeof(float));
				float *corrdata_gwb1=(float*)malloc(corrsize_gwb1*sizeof(float));
				float *phaseData_gsb1=(float*)malloc(corrsize_gsb1*sizeof(float));
				float *phaseData_gwb1=(float*)malloc(corrsize_gwb1*sizeof(float));
				
				for(i=0;i<corrsize_gsb;i++)
				{
					if (phaseData_gsb[i]>=pstart && phaseData_gsb[i]<=pend)
					{
						corrdata_gsb1[pos]=corrdata_gsb[i]/50;
						phaseData_gsb1[pos]=phaseData_gsb[i];
					}
				}
				corrpos=0;
				for(i=0;i<corrsize_gwb;i++)
				{
					if (phaseData_gwb[i]>=pstart && phaseData_gwb[i]<=pend)
					{
						corrdata_gwb1[pos]=corrdata_gwb[i]/50;
						phaseData_gwb1[pos]=phaseData_gwb[i];
					}
				}
				corrpos=0;
				for(i=0;i<corrsize_gsb1;i++)
				{
					if(corrdata_gsb1[i]>maxRes)
						maxRes=corrdata_gsb1[i];
						
					if(corrdata_gsb1[i]<minRes)
						minRes=corrdata_gsb1[i];
					printf("\n%lf\t%lf",phaseData_gsb1[i],corrdata_gsb[i]);
				}	
				
				for(i=0;i<corrsize_gwb1;i++)
				{
					if(corrdata_gwb1[i]>maxRes)
						maxRes=corrdata_gwb1[i];
				
					if(corrdata_gwb1[i]<minRes)
						minRes=corrdata_gwb1[i];
				}

				cpgeras;

				cpgsubp(1,4);

				cpgsci(1);
        			cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       				char title[50];
				sprintf(title,"Residuals for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Phase","Residual Intensity",title);
				cpgsci(2);
        			cpgline(corrsize_gsb1,phaseData_gsb1,corrdata_gsb1);
        			cpgsci(3);
        			cpgline(corrsize_gwb1,phaseData_gwb1,corrdata_gwb1);

				cpgpanl(1,1);
        			cpgsci(1);
        			cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
				sprintf(title,"Smoothed pulse for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Phase","Intensity of smoothed part",title);
				cpgsci(2);
        			cpgline(corrsize_gsb,phaseData_gsb,smoothGSB);
        			cpgsci(3);
        			cpgline(corrsize_gwb,phaseData_gwb,smoothGWB);
	
	
        	
        			cpgpanl(1,2);
        			cpgsci(1);
        			cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        			sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Intensity. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Time lag (ms)","Autocorrelation amplitude",title);
				cpgsci(2);
        			cpgline(autosize_gsb,lagsGSB,autocorrGSB);
        			cpgsci(3);
        			cpgline(autosize_gwb,lagsGWB,autocorrGWB);

				cpgpanl(1,3);
				cpgsci(1);
				cpgenv(lags[0],lags[crosssize-1],0,1,0,0);
				sprintf(title,"Cross-correlation of residual part (HF is held stationary and LF is moved from +lag to -lag) for pulse number %d",pulseNum);		
				cpglab("Time lag (ms)","Cross-correlation amplitude",title);
	        		cpgline(crosssize,lags,crossfloat);

				printf("\nDo you want to re-size the phase window? 1=yes, 0=no\n");
				scanf("%d",&microphase);
			}*/			
				
                                                  

	
	//Vincross

			


			free(phaseData_gsb);
			free(phaseData_gwb);
			free(corrdata_gsb);
			free(corrdata_gwb);
			free(lagsGSB);
			free(lagsGWB);
			free(autocorrGSB);
			free(autocorrGWB);
			free(autocorrGSB_SG);
			free(autocorrGWB_SG);
			break;
		}
		case 1:
		{
			FILE *SG_coeff_file=fopen(SG_coeffs,"r");
			int sg_read[13];
			double sg_coeff[13];
			int sg_order, sg_norm,i;

			if(SG_coeff_file==NULL)
			{
				printf("\nCoefficients file not found!");
				return;
			}	

			while(1)
			{
				if(fscanf(SG_coeff_file,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",&sg_order,&sg_norm,sg_read,sg_read+1,sg_read+2,sg_read+3,sg_read+4, sg_read+5, sg_read+6, sg_read+7, sg_read+8,sg_read+9, sg_read+10,sg_read+11,sg_read+12)==15)
				{
					if(sg_order==2*sg_coeff_size+1)
					{
						printf("\nFound coefficients for smoothing window of %d points with normalisation of %d\n. Coefficients are:",sg_order,sg_norm);
						for(i=0;i<13;i++)
						{
							sg_coeff[i]=(double)sg_read[i]/sg_norm;
							printf("%f  ",sg_coeff[i]);
						}
						printf("\n");
						break;
					}
					else
						continue;
				}
				else
				{
					printf("\nsg Coeffs not found!\n");
					return;
				}

			}
			fclose(SG_coeff_file);	
	
			FILE *data_GSB, *autoFile;
			data_GSB=fopen(filename_GSB,"rb");
			autoFile=fopen(autoName,"w");
			int numshift_gsb,j;
			int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), sampleDelay_gsb=gsbDelay/t_gsb*1000;

			numshift_gsb=(int)(last_gsb*s);
			fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);	
			readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB);
	
			int corrsize_gsb, corrpos=0;
			corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
			float dpoint=0,ppoint=0;
			float *corrdata_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
			float *phaseData_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	
			double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),range_gsb=max(folda_gsb,last_gsb)-baseline_gsb;
			for(i=0;i<last_gsb;)
			{
				dpoint=0;
				ppoint=0;
				for(j=i;j<i+intg_gsb && j<last_gsb;j++)
				{
					dpoint+=folda_gsb[j];
					ppoint+=phase_gsb[j];
				}
				if(j==last_gsb)
					break;
				ppoint=ppoint/intg_gsb;
				dpoint=dpoint/intg_gsb;
				dpoint=(dpoint-baseline_gsb)/range_gsb;
		
				if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
				{
					corrdata_gsb[corrpos]=dpoint;
					phaseData_gsb[corrpos]=ppoint;
					corrpos++;
				}
				i+=intg_gsb;
			}
			corrpos=0;
	
			int smoothWindow_gsb, smoothWindow_gwb;	
			float* smoothGSB=(float*)malloc(corrsize_gsb*sizeof(float));
	
			smoothWindow_gsb=(int)(smoothDur_gsb/(t_gsb*intg_gsb));
			printf("\nSmoothing window size is %d bins.\n",smoothWindow_gsb);
			smoothPulse(corrdata_gsb,smoothGSB,corrsize_gsb,smoothWindow_gsb);
	
			float maxSmooth=0,minSmooth=smoothGSB[0],maxRes=0,minRes=0,rmsRes_gsb=0;
	
			for(i=0;i<corrsize_gsb;i++)
			{
				corrdata_gsb[i]=corrdata_gsb[i]-smoothGSB[i];
				rmsRes_gsb+=corrdata_gsb[i]*corrdata_gsb[i];
		
				if(smoothGSB[i]>maxSmooth)
					maxSmooth=smoothGSB[i];

				if(smoothGSB[i]<minSmooth)
					minSmooth=smoothGSB[i];
		
				if(corrdata_gsb[i]>maxRes)
					maxRes=corrdata_gsb[i];
			
				if(corrdata_gsb[i]<minRes)
					minRes=corrdata_gsb[i];
			}
			printf("Power = %f \n",rmsRes_gsb);
	
			for(i=0;i<last_gsb;i++)
				folda_gsb[i]=(folda_gsb[i]-baseline_gsb)/range_gsb;

			float stdev_gsb=stdevBaseIntg(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb,intg_gsb);
	
			int normRes_gsb=corrsize_gsb;
			rmsRes_gsb/=(normRes_gsb*stdev_gsb*stdev_gsb);
	
			rmsRes_gsb=sqrt((double)rmsRes_gsb);

			float maxDev_gsb=maxMod(corrdata_gsb,corrsize_gsb);
			float maxDevRMS_gsb=maxDev_gsb/stdev_gsb;
	
			printf("\nRMS deviation = %f\t for STDEV=%f\t",rmsRes_gsb,stdev_gsb);
			printf("\nMax_RMS_Dev = %f\t",maxDevRMS_gsb);
			float autoWidths[2];
			findWidths_derSG(corrdata_gsb,corrdata_gsb,corrsize_gsb,corrsize_gsb,autoWidths,t_gsb,t_gsb,intg_gsb,intg_gsb,sg_coeff,sg_coeff_size);
			printf("\nWidths detected = %f ms \n",autoWidths[0]);	
	
			int autosize_gsb=(int)(lag/(t_gsb*intg_gsb));
			float* autocorrGSB=(float*)malloc(autosize_gsb*sizeof(float));
			float* lagsGSB=(float*)malloc(autosize_gsb*sizeof(float));
	
			float maxAuto=0,minAuto=1;
		
			float sum=0;
			float autoNorm=0;
	
			for(i=0;i<autosize_gsb;i++)
			{
				sum=0;
				for(j=0;j<corrsize_gsb;j++)
				{
					if(j+i>=corrsize_gsb)
					{
						sum+=corrdata_gsb[j]*corrdata_gsb[j+i-corrsize_gsb];
					}
					else
					{
						sum+=corrdata_gsb[j]*corrdata_gsb[j+i];
					}
				}
		
				autocorrGSB[i]=sum;
		
				if(i==0)
					autoNorm=autocorrGSB[i];

				autocorrGSB[i]=autocorrGSB[i]/autoNorm;
		
				if(autocorrGSB[i]>maxAuto)
				{
					maxAuto=autocorrGSB[i];
				}
				if(autocorrGSB[i]<minAuto)
				{
					minAuto=autocorrGSB[i];
				}
		
				lagsGSB[i]=(t_gsb*intg_gsb*i);
				fprintf(autoFile,"%f\t%f\n",lagsGSB[i],autocorrGSB[i]);
			}
	
			float *autocorrGSB_SG=(float*)malloc(autosize_gsb*sizeof(float));
	
			for(i=0;i<autosize_gsb;i++)
			{
				autocorrGSB_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							autocorrGSB_SG[i]+=autocorrGSB[-j]*sg_coeff[i-j];
						}
						else
						{
							autocorrGSB_SG[i]+=autocorrGSB[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							autocorrGSB_SG[i]+=autocorrGSB[j]*sg_coeff[i-j];
						}
						else
						{
							autocorrGSB_SG[i]+=autocorrGSB[j]*sg_coeff[j-i];
						}
					}
				}
			}



			for(i=0;i<autosize_gsb;i++)
			{
				autocorrGSB[i]=autocorrGSB_SG[i];
			}
	
        		fclose(autoFile);
        		fclose(data_GSB);
        
        		cpgsubp(1,3);
        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       			char title[50];
			sprintf(title,"Residuals for pulse number %d",pulseNum);
			cpglab("Phase","Residual Intensity",title);
			cpgsci(2);
        		cpgline(corrsize_gsb,phaseData_gsb,corrdata_gsb);
        
        		cpgpanl(1,1);
        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
			sprintf(title,"Smoothed pulse for pulse number %d",pulseNum);
			cpglab("Phase","Intensity of smoothed part",title);
			cpgsci(2);
        		cpgline(corrsize_gsb,phaseData_gsb,smoothGSB);
        
        		cpgpanl(1,2);
        		cpgsci(1);
        		cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        		sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Intensity",pulseNum);
			cpglab("Time lag (ms)","Autocorrelation amplitude",title);
			cpgsci(2);
        		cpgline(autosize_gsb,lagsGSB,autocorrGSB);


			free(phaseData_gsb);
			free(corrdata_gsb);
			free(lagsGSB);
			free(autocorrGSB);
			free(autocorrGSB_SG);		
			break;
		}
	}
}


void smoothSubACF(char *filename_GSB, char *filename_GWB, char *autoName, double *phase_gsb, double *phase_gwb, double s, double *folda_gsb, double *folda_gwb, int last_gsb, int last_gwb,double pstart, double pend, double pulseStart, double pulseEnd, int intg_gsb, int intg_gwb, int pulseNum, double fold, double t_gsb, double t_gwb, double lag, double smoothDur_gsb, double smoothDur_gwb, double gsbDelay, double gwbDelay,float *readData_gsb,float *readData_gwb)
{
	FILE *data_GSB, *data_GWB, *autoFile;
	data_GSB=fopen(filename_GSB,"rb");
	data_GWB=fopen(filename_GWB,"rb");
	autoFile=fopen(autoName,"w");
	int i,numshift_gsb, numshift_gwb,j;
	int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), skipnum_gwb=(int)((pulseNum-1)*(fold/t_gwb)), sampleDelay_gsb=gsbDelay/t_gsb*1000, sampleDelay_gwb=gwbDelay/t_gwb*1000;

	numshift_gsb=(int)(last_gsb*s);
	numshift_gwb=(int)(last_gwb*s);	
	fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);
	fseek(data_GWB,(skipnum_gwb+sampleDelay_gwb+numshift_gwb)*sizeof(float),SEEK_SET);	
	readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB);
	readFloatToDouble(last_gwb,readData_gwb,folda_gwb,data_GWB);
	
	int corrsize_gsb, corrsize_gwb, corrpos=0;
	corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
	corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *corrdata_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	float *phaseData_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *phaseData_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	
	double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),range_gsb=max(folda_gsb,last_gsb)-baseline_gsb,baseline_gwb=meanBase(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb),range_gwb=max(folda_gwb,last_gwb)-baseline_gwb;
	for(i=0;i<last_gsb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gsb && j<last_gsb;j++)
		{
			dpoint+=folda_gsb[j];
			ppoint+=phase_gsb[j];
		}
		if(j==last_gsb)
			break;
		ppoint=ppoint/intg_gsb;
		dpoint=dpoint/intg_gsb;
		dpoint=(dpoint-baseline_gsb)/range_gsb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
		{
			corrdata_gsb[corrpos]=dpoint;
			phaseData_gsb[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_gsb;
	}
	corrpos=0;
	for(i=0;i<last_gwb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gwb && j<last_gwb;j++)
		{
			dpoint+=folda_gwb[j];
			ppoint+=phase_gwb[j];
		}
		if(j==last_gwb)
			break;
		ppoint=ppoint/intg_gwb;
		dpoint=dpoint/intg_gwb;
		dpoint=(dpoint-baseline_gwb)/range_gwb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gwb)
		{
			corrdata_gwb[corrpos]=dpoint;
			phaseData_gwb[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_gwb;
	}	
	
	int smoothWindow_gsb, smoothWindow_gwb;	
	float* smoothGSB=(float*)malloc(corrsize_gsb*sizeof(float));
	float* smoothGWB=(float*)malloc(corrsize_gwb*sizeof(float));
	
	smoothWindow_gsb=(int)(smoothDur_gsb/(t_gsb*intg_gsb));
	smoothWindow_gwb=(int)(smoothDur_gwb/(t_gwb*intg_gwb));
	printf("\nSmoothing window size is %d bins for GSB; %d bins for GWB.\n",smoothWindow_gsb,smoothWindow_gwb);
	smoothPulse(corrdata_gsb,smoothGSB,corrsize_gsb,smoothWindow_gsb);
	smoothPulse(corrdata_gwb,smoothGWB,corrsize_gwb,smoothWindow_gwb);
	
	float maxSmooth=0,minSmooth=smoothGSB[0],maxRes=0,minRes=0,rmsRes_gsb=0,rmsRes_gwb=0;
	
	for(i=0;i<corrsize_gsb;i++)
	{
		corrdata_gsb[i]=corrdata_gsb[i]-smoothGSB[i];
		rmsRes_gsb+=corrdata_gsb[i]*corrdata_gsb[i];
		
		if(smoothGSB[i]>maxSmooth)
			maxSmooth=smoothGSB[i];

		if(smoothGSB[i]<minSmooth)
			minSmooth=smoothGSB[i];
		
		if(corrdata_gsb[i]>maxRes)
			maxRes=corrdata_gsb[i];
			
		if(corrdata_gsb[i]<minRes)
			minRes=corrdata_gsb[i];
	}
	for(i=0;i<corrsize_gwb;i++)
	{
		corrdata_gwb[i]=corrdata_gwb[i]-smoothGWB[i];
		rmsRes_gwb+=corrdata_gwb[i]*corrdata_gwb[i];
		
		if(smoothGWB[i]>maxSmooth)
			maxSmooth=smoothGWB[i];

		if(smoothGWB[i]<minSmooth)
			minSmooth=smoothGWB[i];
		
		if(corrdata_gwb[i]>maxRes)
			maxRes=corrdata_gwb[i];
			
		if(corrdata_gwb[i]<minRes)
			minRes=corrdata_gwb[i];
	}
	
	for(i=0;i<last_gsb;i++)
		folda_gsb[i]=(folda_gsb[i]-baseline_gsb)/range_gsb;
	for(i=0;i<last_gwb;i++)
		folda_gwb[i]=(folda_gwb[i]-baseline_gwb)/range_gwb;
		
	float stdev_gsb=stdevBaseIntg(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb,intg_gsb),stdev_gwb=stdevBaseIntg(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb,intg_gwb);
	
	int normRes_gsb=corrsize_gsb, normRes_gwb=corrsize_gwb;
	rmsRes_gsb/=(normRes_gsb*stdev_gsb*stdev_gsb);
	rmsRes_gwb/=(normRes_gwb*stdev_gwb*stdev_gwb);
	
	rmsRes_gsb=sqrt((double)rmsRes_gsb);
	rmsRes_gwb=sqrt((double)rmsRes_gwb);
	
	printf("\nRMS deviation of GSB = %f\tGWB = %f for STDEV_GSB=%f\tSTDEV_GWB=%f",rmsRes_gsb,rmsRes_gwb,stdev_gsb,stdev_gwb);
	
	float autoWidths[2];
	findWidths(corrdata_gsb,corrdata_gwb,corrsize_gsb,corrsize_gwb,autoWidths,t_gsb,t_gwb,intg_gsb,intg_gwb);
	printf("\nWidths detected GSB = %f ms\tGWB = %f ms\n",autoWidths[0],autoWidths[1]);	
	
	int autosize_gsb=(int)(lag/(t_gsb*intg_gsb));
	int autosize_gwb=(int)(lag/(t_gwb*intg_gwb));
	float* autocorrGSB=(float*)malloc(autosize_gsb*sizeof(float));
	float* autocorrGWB=(float*)malloc(autosize_gwb*sizeof(float));
	float* lagsGSB=(float*)malloc(autosize_gsb*sizeof(float));
	float* lagsGWB=(float*)malloc(autosize_gwb*sizeof(float));	
	
	float maxAuto=0,minAuto=1;
		
	float sum=0;
	float autoNorm=0;
	
	for(i=0;i<autosize_gsb;i++)
	{
		sum=0;
		for(j=0;j<corrsize_gsb;j++)
		{
			if(j+i>=corrsize_gsb)
			{
				sum+=corrdata_gsb[j]*corrdata_gsb[j+i-corrsize_gsb];
			}
			else
			{
				sum+=corrdata_gsb[j]*corrdata_gsb[j+i];
			}
		}
		
		autocorrGSB[i]=sum;
		
		if(i==0)
			autoNorm=autocorrGSB[i];

		autocorrGSB[i]=autocorrGSB[i]/autoNorm;
		
		if(autocorrGSB[i]>maxAuto)
		{
			maxAuto=autocorrGSB[i];
		}
		if(autocorrGSB[i]<minAuto)
		{
			minAuto=autocorrGSB[i];
		}
		
		lagsGSB[i]=(t_gsb*intg_gsb*i);
		fprintf(autoFile,"%f\t%f\n",lagsGSB[i],autocorrGSB[i]);
	}
	for(i=0;i<autosize_gwb;i++)
	{
		sum=0;
		for(j=0;j<corrsize_gwb;j++)
		{
			if(j+i>=corrsize_gwb)
			{
				sum+=corrdata_gwb[j]*corrdata_gwb[j+i-corrsize_gwb];
			}
			else
			{
				sum+=corrdata_gwb[j]*corrdata_gwb[j+i];
			}
		}
		
		autocorrGWB[i]=sum;
		
		if(i==0)
			autoNorm=autocorrGWB[i];

		autocorrGWB[i]=autocorrGWB[i]/autoNorm;
		
		if(autocorrGWB[i]>maxAuto)
		{
			maxAuto=autocorrGWB[i];
		}
		if(autocorrGWB[i]<minAuto)
		{
			minAuto=autocorrGWB[i];
		}
		
		lagsGWB[i]=(t_gwb*intg_gwb*i);
		fprintf(autoFile,"%f\t%f\n",lagsGWB[i],autocorrGWB[i]);
	}

        fclose(autoFile);
        fclose(data_GSB);
        fclose(data_GWB);
        
        cpgsubp(1,3);
        cpgsci(1);
        cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       	char title[50];
	sprintf(title,"Residuals for pulse number %d",pulseNum);
	cpglab("Phase","Residual Stokes parameters",title);
	cpgsci(2);
        cpgline(corrsize_gsb,phaseData_gsb,corrdata_gsb);
        cpgsci(3);
        cpgline(corrsize_gwb,phaseData_gwb,corrdata_gwb);
        
        cpgpanl(1,1);
        cpgsci(1);
        cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
	sprintf(title,"Smoothed pulse for pulse number %d",pulseNum);
	cpglab("Phase","Stokes parameters",title);
	cpgsci(2);
        cpgline(corrsize_gsb,phaseData_gsb,smoothGSB);
        cpgsci(3);
        cpgline(corrsize_gwb,phaseData_gwb,smoothGWB);
        
        cpgpanl(1,2);
        cpgsci(1);
        cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Stokes parameters",pulseNum);
	cpglab("Time lag (ms)","Autocorrelation amplitude",title);
	cpgsci(2);
        cpgline(autosize_gsb,lagsGSB,autocorrGSB);
        cpgsci(3);
        cpgline(autosize_gwb,lagsGWB,autocorrGWB);
        
	free(phaseData_gsb);
	free(phaseData_gwb);
	free(corrdata_gsb);
	free(corrdata_gwb);
	free(lagsGSB);
	free(lagsGWB);
	free(autocorrGSB);
	free(autocorrGWB);
}

int microCandidates(char *filename_GSB, char *filename_GWB, FILE *microFile, double *phase_gsb, double *phase_gwb, double s, double *folda_gsb, double *folda_gwb, int last_gsb, int last_gwb,double pstart, double pend, double pulseStart, double pulseEnd, int intg_gsb, int intg_gwb, int pulseNum, double fold, double t_gsb, double t_gwb, double lag, double smoothDur_gsb, double smoothDur_gwb, double gsbDelay, double gwbDelay,float *readData_gsb,float *readData_gwb, float cutOffhighfreq,float cutOfflowfreq, float gsbWidth, float gwbWidth)
{
	FILE *data_GSB, *data_GWB;
	data_GWB=fopen(filename_GWB,"rb");
	data_GSB=fopen(filename_GSB,"rb");

	if(data_GSB==NULL)
	{
		printf("\nCould not open GSB file!");
		return 0;
	}
	if(data_GWB==NULL)
	{
		printf("\nCould not open GWB file!");
		return 0;
	}
	
	int i,numshift_gsb, numshift_gwb,j;
	int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), skipnum_gwb=(int)((pulseNum-1)*(fold/t_gwb)), sampleDelay_gsb=gsbDelay/t_gsb*1000, sampleDelay_gwb=gwbDelay/t_gwb*1000;

	numshift_gsb=(int)(last_gsb*s);
	numshift_gwb=(int)(last_gwb*s);	
	fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);
	fseek(data_GWB,(skipnum_gwb+sampleDelay_gwb+numshift_gwb)*sizeof(float),SEEK_SET);	
	
	if(readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB)!=1)
		return 0;
	if(readFloatToDouble(last_gwb,readData_gwb,folda_gwb,data_GWB)!=1)
		return 0;
	
	fclose(data_GSB);
	fclose(data_GWB);
	
	int corrsize_gsb, corrsize_gwb, corrpos=0;
	corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
	corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *corrdata_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	float *phaseData_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *phaseData_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	
	double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),range_gsb=max(folda_gsb,last_gsb)-baseline_gsb,baseline_gwb=meanBase(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb),range_gwb=max(folda_gwb,last_gwb)-baseline_gwb;
	for(i=0;i<last_gsb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gsb && j<last_gsb;j++)
		{
			dpoint+=folda_gsb[j];
			ppoint+=phase_gsb[j];
		}
		if(j==last_gsb)
			break;
		ppoint=ppoint/intg_gsb;
		dpoint=dpoint/intg_gsb;
		dpoint=(dpoint-baseline_gsb)/range_gsb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
		{
			corrdata_gsb[corrpos]=dpoint;
			phaseData_gsb[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_gsb;
	}
	corrpos=0;
	for(i=0;i<last_gwb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gwb && j<last_gwb;j++)
		{
			dpoint+=folda_gwb[j];
			ppoint+=phase_gwb[j];
		}
		if(j==last_gwb)
			break;
		ppoint=ppoint/intg_gwb;
		dpoint=dpoint/intg_gwb;
		dpoint=(dpoint-baseline_gwb)/range_gwb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gwb)
		{
			corrdata_gwb[corrpos]=dpoint;
			phaseData_gwb[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_gwb;
	}	
	
	int smoothWindow_gsb, smoothWindow_gwb;	
	float* smoothGSB=(float*)malloc(corrsize_gsb*sizeof(float));
	float* smoothGWB=(float*)malloc(corrsize_gwb*sizeof(float));
	
	smoothWindow_gsb=(int)(smoothDur_gsb/(t_gsb*intg_gsb));
	smoothWindow_gwb=(int)(smoothDur_gwb/(t_gwb*intg_gwb));
	smoothPulse(corrdata_gsb,smoothGSB,corrsize_gsb,smoothWindow_gsb);
	smoothPulse(corrdata_gwb,smoothGWB,corrsize_gwb,smoothWindow_gwb);
	
	if(gsbWidth==0 || gwbWidth==0)
	{
		printf("\nError! Widths cannot be zero!\n");
		return;
	}
	
	float threshFreq_gsb=1.0/gsbWidth*subPulseFreqThreshold;
	float threshFreq_gwb=1.0/gwbWidth*subPulseFreqThreshold;
	
	float maxSmooth=0,minSmooth=smoothGSB[0],maxRes=0,minRes=0,rmsRes_gsb=0,rmsRes_gwb=0;
	
	for(i=0;i<corrsize_gsb;i++)
	{
		corrdata_gsb[i]=corrdata_gsb[i]-smoothGSB[i];
		rmsRes_gsb+=corrdata_gsb[i]*corrdata_gsb[i];
		
		if(smoothGSB[i]>maxSmooth)
			maxSmooth=smoothGSB[i];

		if(smoothGSB[i]<minSmooth)
			minSmooth=smoothGSB[i];
		
		if(corrdata_gsb[i]>maxRes)
			maxRes=corrdata_gsb[i];
			
		if(corrdata_gsb[i]<minRes)
			minRes=corrdata_gsb[i];
	}
	for(i=0;i<corrsize_gwb;i++)
	{
		corrdata_gwb[i]=corrdata_gwb[i]-smoothGWB[i];
		rmsRes_gwb+=corrdata_gwb[i]*corrdata_gwb[i];
		
		if(smoothGWB[i]>maxSmooth)
			maxSmooth=smoothGWB[i];

		if(smoothGWB[i]<minSmooth)
			minSmooth=smoothGWB[i];
		
		if(corrdata_gwb[i]>maxRes)
			maxRes=corrdata_gwb[i];
			
		if(corrdata_gwb[i]<minRes)
			minRes=corrdata_gwb[i];
	}
	
	for(i=0;i<last_gsb;i++)
		folda_gsb[i]=(folda_gsb[i]-baseline_gsb)/range_gsb;
	for(i=0;i<last_gwb;i++)
		folda_gwb[i]=(folda_gwb[i]-baseline_gwb)/range_gwb;
		
	float stdev_gsb=stdevBaseIntg(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb,intg_gsb),stdev_gwb=stdevBaseIntg(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb,intg_gwb);
	
	int normRes_gsb=corrsize_gsb, normRes_gwb=corrsize_gwb;
	rmsRes_gsb/=(normRes_gsb*stdev_gsb*stdev_gsb);
	rmsRes_gwb/=(normRes_gwb*stdev_gwb*stdev_gwb);
	
	rmsRes_gsb=sqrt((double)rmsRes_gsb);
	rmsRes_gwb=sqrt((double)rmsRes_gwb);
	printf("\nProcessing pulse number %d. RMS deviation of GSB = %f\tGWB = %f",pulseNum,rmsRes_gsb,rmsRes_gwb);
	
	float *autoWidths=(float*)malloc(2*sizeof(float));
	float *quasiPeriods=(float*)malloc(2*sizeof(float));
	float *relArray=(float*)malloc(4*sizeof(float));
	
	if(rmsRes_gsb>=cutOfflowfreq || rmsRes_gwb>=cutOffhighfreq)
	{
		printf("\nExceeded microstructure cut-off.");
		detectQuasiPeriod(corrdata_gsb,corrdata_gwb,corrsize_gsb,corrsize_gwb,quasiPeriods,t_gsb,t_gwb,intg_gsb,intg_gwb,threshFreq_gsb,threshFreq_gwb);
		if(quasiPeriods[0]>0)
			printf("\nQuasiperiod detected GSB at %f ms",quasiPeriods[0]);
		if(quasiPeriods[1]>0)
			printf("\nQuasiperiod detected GWB at %f ms",quasiPeriods[1]);
		
		findWidths(corrdata_gsb,corrdata_gwb,corrsize_gsb,corrsize_gwb,autoWidths,t_gsb,t_gwb,intg_gsb,intg_gwb);

		printf("\nWidths detected GSB = %f\tGWB = %f",autoWidths[0],autoWidths[1]);
		
		giveRelStrengthSmoothed(corrdata_gsb,corrdata_gwb,smoothGSB,smoothGWB,corrsize_gsb,corrsize_gwb,relArray);
		
		printf("\n");
			
		fprintf(microFile,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",pulseNum,rmsRes_gsb,rmsRes_gwb,autoWidths[0],autoWidths[1],quasiPeriods[0],quasiPeriods[1],relArray[0],relArray[1],relArray[2],relArray[3]);
	}
	
	free(corrdata_gsb);
	free(corrdata_gwb);
	free(phaseData_gsb);
	free(phaseData_gwb);
	free(smoothGSB);
	free(smoothGWB);
	free(autoWidths);
	free(quasiPeriods);
	free(relArray);
	
	return 1;
}

int microCandidates_SG(char *filename_GSB, char *filename_GWB, FILE *microFile, char *SG_coeffs, double *phase_gsb, double *phase_gwb, double s, double *folda_gsb, double *folda_gwb, int last_gsb, int last_gwb,double pstart, double pend, double pulseStart, double pulseEnd, int intg_gsb, int intg_gwb, int pulseNum, double fold, double t_gsb, double t_gwb, double lag, double smoothDur_gsb, double smoothDur_gwb, double gsbDelay, double gwbDelay,float *readData_gsb,float *readData_gwb, float cutOffhighfreq,float cutOfflowfreq, float gsbWidth, float gwbWidth, int sg_coeff_size,int option)
{
	switch (option)
	{
	case 2:
	{
	//Savitzky Golay filter coefficients
        FILE *SG_coeff_file=fopen(SG_coeffs,"r");
        int *sg_read=(int*)malloc(13*sizeof(int));
        double *sg_coeff=(double*)malloc(13*sizeof(double));
        int sg_order, sg_norm,i;
                                                                                                                                                                                                                                                          
        if(SG_coeff_file==NULL)
        {
        	printf("\nCoefficients file not found!");
        	return;
        }	
                                                                                                                                                                                                                                                          
        while(1)
        {
        	if(fscanf(SG_coeff_file,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",&sg_order,&sg_norm,sg_read,sg_read+1,sg_read+2,sg_read+3,sg_read+4, sg_read+5, sg_read+6, sg_read+7, sg_read+8,sg_read+9, sg_read+10,sg_read+11,sg_read+12)==15)
        	{
        		if(sg_order==2*sg_coeff_size+1)
        		{
        			printf("\nFound coefficients for smoothing window of %d points with normalisation of %d\n. Coefficients are:",sg_order,sg_norm);
        			for(i=0;i<13;i++)
        			{
        				sg_coeff[i]=(double)sg_read[i]/sg_norm;
        				printf("%f  ",sg_coeff[i]);
        			}
        			printf("\n");
        			break;
        		}
        		else
        			continue;
        	}
        	else
        	{
        		printf("\nsg Coeffs not found!\n");
        		return;
        	}
                                                                                                                                                                                                                                                          
        }
	
	fclose(SG_coeff_file);
	FILE *data_GSB, *data_GWB;
	data_GWB=fopen(filename_GWB,"rb");
	data_GSB=fopen(filename_GSB,"rb");

	if(data_GSB==NULL)
	{
		printf("\nCould not open GSB file!");
		return 0;
	}
	if(data_GWB==NULL)
	{
		printf("\nCould not open GWB file!");
		return 0;
	}
	
	int numshift_gsb, numshift_gwb,j;
	int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), skipnum_gwb=(int)((pulseNum-1)*(fold/t_gwb)), sampleDelay_gsb=gsbDelay/t_gsb*1000, sampleDelay_gwb=gwbDelay/t_gwb*1000;

	numshift_gsb=(int)(last_gsb*s);
	numshift_gwb=(int)(last_gwb*s);	
	fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);
	fseek(data_GWB,(skipnum_gwb+sampleDelay_gwb+numshift_gwb)*sizeof(float),SEEK_SET);	
	
	if(readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB)!=1)
		return 0;
	if(readFloatToDouble(last_gwb,readData_gwb,folda_gwb,data_GWB)!=1)
		return 0;
	
	fclose(data_GSB);
	fclose(data_GWB);
	
	int corrsize_gsb, corrsize_gwb, corrpos=0;
	corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
	corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *corrdata_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	float *phaseData_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *phaseData_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	
	double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),range_gsb=max(folda_gsb,last_gsb)-baseline_gsb,baseline_gwb=meanBase(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb),range_gwb=max(folda_gwb,last_gwb)-baseline_gwb;
	for(i=0;i<last_gsb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gsb && j<last_gsb;j++)
		{
			dpoint+=folda_gsb[j];
			ppoint+=phase_gsb[j];
		}
		if(j==last_gsb)
			break;
		ppoint=ppoint/intg_gsb;
		dpoint=dpoint/intg_gsb;
		dpoint=(dpoint-baseline_gsb)/range_gsb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
		{
			corrdata_gsb[corrpos]=dpoint;
			phaseData_gsb[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_gsb;
	}
	corrpos=0;
	for(i=0;i<last_gwb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gwb && j<last_gwb;j++)
		{
			dpoint+=folda_gwb[j];
			ppoint+=phase_gwb[j];
		}
		if(j==last_gwb)
			break;
		ppoint=ppoint/intg_gwb;
		dpoint=dpoint/intg_gwb;
		dpoint=(dpoint-baseline_gwb)/range_gwb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gwb)
		{
			corrdata_gwb[corrpos]=dpoint;
			phaseData_gwb[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_gwb;
	}	
	
	int smoothWindow_gsb, smoothWindow_gwb;	
	float* smoothGSB=(float*)malloc(corrsize_gsb*sizeof(float));
	float* smoothGWB=(float*)malloc(corrsize_gwb*sizeof(float));
	
	smoothWindow_gsb=(int)(smoothDur_gsb/(t_gsb*intg_gsb));
	smoothWindow_gwb=(int)(smoothDur_gwb/(t_gwb*intg_gwb));
	smoothPulse(corrdata_gsb,smoothGSB,corrsize_gsb,smoothWindow_gsb);
	smoothPulse(corrdata_gwb,smoothGWB,corrsize_gwb,smoothWindow_gwb);
	
	if(gsbWidth==0 || gwbWidth==0)
	{
		printf("\nError! Widths cannot be zero!\n");
		return;
	}
	
	float threshFreq_gsb=1.0/gsbWidth*subPulseFreqThreshold;
	float threshFreq_gwb=1.0/gwbWidth*subPulseFreqThreshold;
	
	float maxSmooth=0,minSmooth=smoothGSB[0],maxRes=0,minRes=0,rmsRes_gsb=0,rmsRes_gwb=0,microstrPowergsb=0,microstrPowergwb=0;
	
	for(i=0;i<corrsize_gsb;i++)
	{
		corrdata_gsb[i]=corrdata_gsb[i]-smoothGSB[i];
		rmsRes_gsb+=corrdata_gsb[i]*corrdata_gsb[i];
		
		if(smoothGSB[i]>maxSmooth)
			maxSmooth=smoothGSB[i];

		if(smoothGSB[i]<minSmooth)
			minSmooth=smoothGSB[i];
		
		if(corrdata_gsb[i]>maxRes)
			maxRes=corrdata_gsb[i];
			
		if(corrdata_gsb[i]<minRes)
			minRes=corrdata_gsb[i];
	}
	microstrPowergsb=rmsRes_gsb;
	for(i=0;i<corrsize_gwb;i++)
	{
		corrdata_gwb[i]=corrdata_gwb[i]-smoothGWB[i];
		rmsRes_gwb+=corrdata_gwb[i]*corrdata_gwb[i];
		
		if(smoothGWB[i]>maxSmooth)
			maxSmooth=smoothGWB[i];

		if(smoothGWB[i]<minSmooth)
			minSmooth=smoothGWB[i];
		
		if(corrdata_gwb[i]>maxRes)
			maxRes=corrdata_gwb[i];
			
		if(corrdata_gwb[i]<minRes)
			minRes=corrdata_gwb[i];
	}
	microstrPowergwb=rmsRes_gwb;
	
	for(i=0;i<last_gsb;i++)
		folda_gsb[i]=(folda_gsb[i]-baseline_gsb)/range_gsb;
	for(i=0;i<last_gwb;i++)
		folda_gwb[i]=(folda_gwb[i]-baseline_gwb)/range_gwb;
		
	float stdev_gsb=stdevBaseIntg(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb,intg_gsb),stdev_gwb=stdevBaseIntg(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb,intg_gwb);
	
	int normRes_gsb=corrsize_gsb, normRes_gwb=corrsize_gwb;
	rmsRes_gsb/=(normRes_gsb*stdev_gsb*stdev_gsb);
	rmsRes_gwb/=(normRes_gwb*stdev_gwb*stdev_gwb);
	
	rmsRes_gsb=sqrt((double)rmsRes_gsb);
	rmsRes_gwb=sqrt((double)rmsRes_gwb);
	
	float maxDev_gsb=maxMod(corrdata_gsb,corrsize_gsb),maxDev_gwb=maxMod(corrdata_gwb,corrsize_gwb);
        float maxDevRMS_gsb=maxDev_gsb/stdev_gsb,maxDevRMS_gwb=maxDev_gwb/stdev_gwb;
	
	printf("\nProcessing pulse number %d. RMS deviation of GSB = %f\tGWB = %f. MaxDevRMS GSB = %f GWB = %f",pulseNum,rmsRes_gsb,rmsRes_gwb,maxDevRMS_gsb,maxDevRMS_gwb);
	
	float *autoWidths=(float*)malloc(2*sizeof(float));
	float *quasiPeriods=(float*)malloc(2*sizeof(float));
	float *relArray=(float*)malloc(4*sizeof(float));
	
	if(rmsRes_gsb>=cutOfflowfreq || rmsRes_gwb>=cutOffhighfreq || maxDevRMS_gsb>=microThreshold || maxDevRMS_gwb>=microThreshold)
	{
		printf("\nExceeded microstructure cut-off.");
		detectQuasiPeriod(corrdata_gsb,corrdata_gwb,corrsize_gsb,corrsize_gwb,quasiPeriods,t_gsb,t_gwb,intg_gsb,intg_gwb,threshFreq_gsb,threshFreq_gwb);
		if(quasiPeriods[0]>0)
			printf("\nQuasiperiod detected GSB at %f ms",quasiPeriods[0]);
		if(quasiPeriods[1]>0)
			printf("\nQuasiperiod detected GWB at %f ms",quasiPeriods[1]);
		
		findWidths_derSG(corrdata_gsb,corrdata_gwb,corrsize_gsb,corrsize_gwb,autoWidths,t_gsb,t_gwb,intg_gsb,intg_gwb,sg_coeff,sg_coeff_size);

		printf("\nWidths detected GSB = %f\tGWB = %f",autoWidths[0],autoWidths[1]);
		
		giveRelStrengthSmoothed(corrdata_gsb,corrdata_gwb,smoothGSB,smoothGWB,corrsize_gsb,corrsize_gwb,relArray);
		
		printf("\n");
			
		fprintf(microFile,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n",pulseNum,rmsRes_gsb,rmsRes_gwb,maxDevRMS_gsb,maxDevRMS_gwb,autoWidths[0],autoWidths[1],quasiPeriods[0],quasiPeriods[1],relArray[0],relArray[1],relArray[2],relArray[3],microstrPowergsb,microstrPowergwb);
	}
	
	free(corrdata_gsb);
	free(corrdata_gwb);
	free(phaseData_gsb);
	free(phaseData_gwb);
	free(smoothGSB);
	free(smoothGWB);
	free(autoWidths);
	free(quasiPeriods);
	free(relArray);
	
	return 1;
	}
	case 1:
	{
	FILE *SG_coeff_file=fopen(SG_coeffs,"r");
        int *sg_read=(int*)malloc(13*sizeof(int));
        double *sg_coeff=(double*)malloc(13*sizeof(double));
        int sg_order, sg_norm,i;
                                                                                                                                                                                                                                                          
        if(SG_coeff_file==NULL)
        {
        	printf("\nCoefficients file not found!");
        	return;
        }	
                                                                                                                                                                                                                                                          
        while(1)
        {
        	if(fscanf(SG_coeff_file,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",&sg_order,&sg_norm,sg_read,sg_read+1,sg_read+2,sg_read+3,sg_read+4, sg_read+5, sg_read+6, sg_read+7, sg_read+8,sg_read+9, sg_read+10,sg_read+11,sg_read+12)==15)
        	{
        		if(sg_order==2*sg_coeff_size+1)
        		{
        			printf("\nFound coefficients for smoothing window of %d points with normalisation of %d\n. Coefficients are:",sg_order,sg_norm);
        			for(i=0;i<13;i++)
        			{
        				sg_coeff[i]=(double)sg_read[i]/sg_norm;
        				printf("%f  ",sg_coeff[i]);
        			}
        			printf("\n");
        			break;
        		}
        		else
        			continue;
        	}
        	else
        	{
        		printf("\nsg Coeffs not found!\n");
        		return;
        	}
                                                                                                                                                                                                                                                          
        }
	
	fclose(SG_coeff_file);
	FILE *data_GSB;
	data_GSB=fopen(filename_GSB,"rb");

	if(data_GSB==NULL)
	{
		printf("\nCould not open GSB file!");
		return 0;
	}
		
	int numshift_gsb,j;
	int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), sampleDelay_gsb=gsbDelay/t_gsb*1000;

	numshift_gsb=(int)(last_gsb*s);
	fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);	
	
	if(readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB)!=1)
		return 0;
	
	fclose(data_GSB);
	
	int corrsize_gsb, corrpos=0;
	corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *phaseData_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	
	double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),range_gsb=max(folda_gsb,last_gsb)-baseline_gsb;
	for(i=0;i<last_gsb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gsb && j<last_gsb;j++)
		{
			dpoint+=folda_gsb[j];
			ppoint+=phase_gsb[j];
		}
		if(j==last_gsb)
			break;
		ppoint=ppoint/intg_gsb;
		dpoint=dpoint/intg_gsb;
		dpoint=(dpoint-baseline_gsb)/range_gsb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
		{
			corrdata_gsb[corrpos]=dpoint;
			phaseData_gsb[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_gsb;
	}
	corrpos=0;
	
	int smoothWindow_gsb, smoothWindow_gwb;	
	float* smoothGSB=(float*)malloc(corrsize_gsb*sizeof(float));
	
	smoothWindow_gsb=(int)(smoothDur_gsb/(t_gsb*intg_gsb));
	smoothPulse(corrdata_gsb,smoothGSB,corrsize_gsb,smoothWindow_gsb);
	
	if(gsbWidth==0)
	{
		printf("\nError! Widths cannot be zero!\n");
		return;
	}
	
	float threshFreq_gsb=1.0/gsbWidth*subPulseFreqThreshold;
	
	float maxSmooth=0,minSmooth=smoothGSB[0],maxRes=0,minRes=0,rmsRes_gsb=0,microstrPowergsb=0;
	
	for(i=0;i<corrsize_gsb;i++)
	{
		corrdata_gsb[i]=corrdata_gsb[i]-smoothGSB[i];
		rmsRes_gsb+=corrdata_gsb[i]*corrdata_gsb[i];
		
		if(smoothGSB[i]>maxSmooth)
			maxSmooth=smoothGSB[i];

		if(smoothGSB[i]<minSmooth)
			minSmooth=smoothGSB[i];
		
		if(corrdata_gsb[i]>maxRes)
			maxRes=corrdata_gsb[i];
			
		if(corrdata_gsb[i]<minRes)
			minRes=corrdata_gsb[i];
	}
	microstrPowergsb=rmsRes_gsb;
	
	for(i=0;i<last_gsb;i++)
		folda_gsb[i]=(folda_gsb[i]-baseline_gsb)/range_gsb;
		
	float stdev_gsb=stdevBaseIntg(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb,intg_gsb);

	int normRes_gsb=corrsize_gsb;
	rmsRes_gsb/=(normRes_gsb*stdev_gsb*stdev_gsb);
	
	rmsRes_gsb=sqrt((double)rmsRes_gsb);
	
	float maxDev_gsb=maxMod(corrdata_gsb,corrsize_gsb);
        float maxDevRMS_gsb=maxDev_gsb/stdev_gsb;
	
	printf("\nProcessing pulse number %d. RMS deviation = %f . MaxDevRMS GSB = %f",pulseNum,rmsRes_gsb,maxDevRMS_gsb);
	
	float *autoWidths=(float*)malloc(2*sizeof(float));
	float *quasiPeriods=(float*)malloc(2*sizeof(float));
	float *relArray=(float*)malloc(4*sizeof(float));
	
	if(rmsRes_gsb>=cutOfflowfreq || maxDevRMS_gsb>=microThreshold )
	{
		printf("\nExceeded microstructure cut-off.");
		detectQuasiPeriod(corrdata_gsb,corrdata_gsb,corrsize_gsb,corrsize_gsb,quasiPeriods,t_gsb,t_gsb,intg_gsb,intg_gsb,threshFreq_gsb,threshFreq_gsb);
		if(quasiPeriods[0]>0)
			printf("\nQuasiperiod detected GSB at %f ms",quasiPeriods[0]);
		
		findWidths_derSG(corrdata_gsb,corrdata_gsb,corrsize_gsb,corrsize_gsb,autoWidths,t_gsb,t_gsb,intg_gsb,intg_gsb,sg_coeff,sg_coeff_size);

		printf("\nWidths detected = %f ",autoWidths[0]);
		
		giveRelStrengthSmoothed(corrdata_gsb,corrdata_gsb,smoothGSB,smoothGSB,corrsize_gsb,corrsize_gsb,relArray);
		
		printf("\n");
			
		fprintf(microFile,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",pulseNum,rmsRes_gsb,maxDevRMS_gsb,autoWidths[0],quasiPeriods[0],relArray[0],relArray[2],microstrPowergsb);
	}
	
	free(corrdata_gsb);
	free(phaseData_gsb);
	free(smoothGSB);
	free(autoWidths);
	free(quasiPeriods);
	free(relArray);
	return 1;
	}
	}
}


void correlate(char *filename, char *outname, int pulseNum, int last, double s, double *folda,double pstart, double pend,double *phase,int intg, double fold, double t, double lag, double pulseStart, double pulseEnd, float *readData, double relDelay, int freqFlag)
{
	FILE *data,*out;
	data=fopen(filename,"rb");
	out=fopen(outname,"w");
	int i,numshift,j;
	int skipnum=(int)((pulseNum-1)*(fold/t)), sampleDelay=relDelay/t*1000;
	
	/*for(i=0;i<(skipnum/last);i++)
	{
		fread(folda,sizeof(double),last,data);
	}
	fread(folda,sizeof(double),(skipnum%last),data);*/

	numshift=(int)(last*s);	
	fseek(data,(skipnum+sampleDelay+numshift)*sizeof(float),SEEK_SET);
	readFloatToDouble(last,readData,folda,data);
	
	double baseline=meanBase(folda,phase,pulseStart,pulseEnd,last),range=max(folda,last)-baseline, dpoint=0,ppoint=0;
	int corrsize=(int)(last/intg*(pend-pstart)),corrpos=0;
	double *corrdata=(double*)malloc(corrsize*sizeof(double));
	
	for(i=0;i<last;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg && j<last;j++)
		{
			dpoint+=folda[j];
			ppoint+=phase[j];
		}
		if(j==last)
			break;
		ppoint=ppoint/intg;
		dpoint=dpoint/intg;
		dpoint=(dpoint-baseline)/range;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize)
		{
			corrdata[corrpos++]=dpoint;
		}
		i+=intg;
	}
	
	int autosize=(int)(lag/(t*intg));
	double *autocorr=(double*)malloc(autosize*sizeof(double));
	double sum=0, AutoNorm=0;
	for(i=0;i<autosize;i++)
	{
		sum=0;
		for(j=0;j<corrsize;j++)
		{
			if(j+i>=corrsize)
			{
				sum+=corrdata[j]*corrdata[j+i-corrsize];
			}
			else
			{
				sum+=corrdata[j]*corrdata[j+i];
			}
		}
		
		autocorr[i]=sum;
		if(i==0)
			AutoNorm=autocorr[i];
		
		autocorr[i]/=AutoNorm;
	}
	
	float lags[autosize];
	float autofloat[autosize];
	float minAuto=autocorr[autosize-1],maxAuto=0;
	for(i=0;i<autosize;i++)
        {	
                fprintf(out,"%.10lf\t%.10lf\n",(double)(t*intg*i),autocorr[i]);
        }


        fclose(out);
	
	for(i=0;i<autosize;i++)
	{
		lags[i]=(float)(t*intg*i);
		autofloat[i]=(float)autocorr[i];
		if(autofloat[i]>maxAuto)
			maxAuto=autofloat[i];
		else if(autofloat[i]<minAuto)
			minAuto=autofloat[i];
	}
	
	cpgsci(1);
	cpgenv(lags[0],lags[autosize-1],minAuto,maxAuto,0,0);
	char title[50];

	if(freqFlag==0)
		sprintf(title,"Autocorrelation plot for pulse number %d at High Freq",pulseNum);		
	else
		sprintf(title,"Autocorrelation plot for pulse number %d at Low Freq",pulseNum);

        cpglab("Time lag (ms)","Autocorrelation amplitude",title);
        cpgline(autosize,lags,autofloat);

  	fclose(data);
}

void createAverageSpectra(char *filename_GSB,char *filename_GWB, char *outname, char *foldFile_GSB, char *foldFile_GWB, int last_gsb, int last_gwb, double s, double *folda_gsb, double *folda_gwb ,double pstart, double pend, double pulseStart, double pulseEnd, double *phase_gsb, double *phase_gwb, int intg_gsb, int intg_gwb, double fold, double t_gsb, double t_gwb, float *readData_gsb, float *readData_gwb, double gsbDelay, double gwbDelay)
{
	printf("\nPlease wait ... \n");
	FILE *data_GSB, *data_GWB, *out;
	data_GSB=fopen(filename_GSB,"rb");
	data_GWB=fopen(filename_GWB,"rb");
	out=fopen(outname,"w");
	int i,numshift_gsb, numshift_gwb,j;
	int corrsize_gsb, corrsize_gwb, corrpos=0;
	corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
	corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
	double dpoint=0,ppoint=0;
	float *corrdata_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *corrdata_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	float *phaseData_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *phaseData_gwb=(float*)malloc(corrsize_gwb*sizeof(float));

	numshift_gsb=(int)(last_gsb*s);
	numshift_gwb=(int)(last_gwb*s);	
	int sampleDelay_gsb=gsbDelay/t_gsb*1000, sampleDelay_gwb=gwbDelay/t_gwb*1000;
	fseek(data_GSB,(sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);
	fseek(data_GWB,(sampleDelay_gwb+numshift_gwb)*sizeof(float),SEEK_SET);
	
	double baseline_gsb,baseline_gwb;
	float maxSpec=0,minSpec=0;
	
	float stdev_gsb, stdev_gwb;
	
	float *crossArray1, *crossArray2;
	double t_corr, crossArray1_stdev, crossArray2_stdev; 
	int corrsize,flag;
	
	if(corrsize_gsb>corrsize_gwb)
	{
		crossArray1=corrdata_gsb;
		crossArray2=(float*)malloc(corrsize_gsb*sizeof(float));
		t_corr=t_gsb*intg_gsb;
		corrsize=corrsize_gsb;
		flag=1;
	}
	else
	{
		crossArray1=corrdata_gwb;
		crossArray2=(float*)malloc(corrsize_gwb*sizeof(float));	
		t_corr=t_gwb*intg_gwb;
		corrsize=corrsize_gwb;
		flag=2;
	}
	
	float *specArray=(float*)malloc(corrsize*sizeof(float));
	float specSum=0,specMean=0;
	int numSpecs=0,b;
	long int foldSample=0;
	float *foldSpectra=(float*)malloc(corrsize*sizeof(float));
	int *countSpectra=(int*)malloc(corrsize*sizeof(int));
	memset(foldSpectra,0,corrsize*sizeof(float));
	memset(countSpectra,0,corrsize*sizeof(float));
		
	while(readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB)==1 && readFloatToDouble(last_gwb,readData_gwb,folda_gwb,data_GWB)==1)
	{
		baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb);
		baseline_gwb=meanBase(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb);
		
		corrpos=0;
		for(i=0;i<last_gsb;)
		{
			dpoint=0;
			ppoint=0;
			for(j=i;j<i+intg_gsb && j<last_gsb;j++)
			{
				dpoint+=folda_gsb[j];
				ppoint+=phase_gsb[j];
			}
			if(j==last_gsb)
				break;
			ppoint=ppoint/intg_gsb;
			dpoint=dpoint/intg_gsb;
			dpoint=(dpoint-baseline_gsb)/baseline_gsb;
		
			if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
			{
				corrdata_gsb[corrpos]=dpoint;
				phaseData_gsb[corrpos]=ppoint;	
				corrpos++;
			}
			i+=intg_gsb;
		
		}
		corrpos=0;
		for(i=0;i<last_gwb;)
		{
			dpoint=0;
			ppoint=0;
			for(j=i;j<i+intg_gwb && j<last_gwb;j++)
			{
				dpoint+=folda_gwb[j];
				ppoint+=phase_gwb[j];
			}
			if(j==last_gwb)
				break;
			ppoint=ppoint/intg_gwb;
			dpoint=dpoint/intg_gwb;
			dpoint=(dpoint-baseline_gwb)/baseline_gwb;
		
			if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gwb)
			{
				corrdata_gwb[corrpos]=dpoint;
				phaseData_gwb[corrpos]=ppoint;
				corrpos++;
			}
			i+=intg_gwb;
		}	
		for(i=0;i<last_gsb;i++)
			folda_gsb[i]=(folda_gsb[i]-baseline_gsb)/baseline_gsb;
		for(i=0;i<last_gwb;i++)
			folda_gwb[i]=(folda_gwb[i]-baseline_gwb)/baseline_gwb;
		
		stdev_gsb=stdevBaseIntg(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb,intg_gsb);
		stdev_gwb=stdevBaseIntg(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb,intg_gwb);

		if(corrsize_gsb>corrsize_gwb)
		{
			interpolate_float(corrdata_gwb,corrdata_gsb,t_gwb*intg_gwb,t_gsb*intg_gsb,crossArray2,corrsize_gsb);
			crossArray1_stdev=stdev_gsb;
			crossArray2_stdev=stdev_gwb;
		}
		else
		{
			interpolate_float(corrdata_gsb,corrdata_gwb,t_gsb*intg_gsb,t_gwb*intg_gwb,crossArray2,corrsize_gwb);
			crossArray1_stdev=stdev_gwb;
			crossArray2_stdev=stdev_gsb;
		}
		
		for(i=0;i<corrsize;i++)
		{
			if(crossArray1[i]>spectraPlotThresh*crossArray1_stdev && crossArray2[i]>spectraPlotThresh*crossArray2_stdev)
			{
				if(flag==1)
					specArray[i]=crossArray1[i]/crossArray2[i];
				else
					specArray[i]=crossArray2[i]/crossArray1[i];
					
				specSum+=specArray[i];
				numSpecs++;
			}
			else
				specArray[i]=-1;
		}
		for(i=0;i<corrsize;i++)
		{
			if(specArray[i]!=-1)
			{
				/*b=(int)floor((foldSample*t_corr-fold*(pend-pstart)*floor((double)foldSample*t_corr/(fold*(pend-pstart))))/t_corr);
				foldSample++;*/
				foldSpectra[i]+=(specArray[i]-foldSpectra[i])/(countSpectra[i]+1);				
				countSpectra[i]++;
			}
		}
	}
	
	maxSpec=0;
	for(i=0;i<corrsize;i++)
	{
		if(foldSpectra[i]!=0)
		{
			specSum+=foldSpectra[i];
			numSpecs++;
		}
	}
	specMean=specSum/numSpecs;
	for(i=0;i<corrsize;i++)
	{
		if(foldSpectra[i]==0)
			foldSpectra[i]=specMean;
	
		if(flag==1)
			fprintf(out,"%f\t%f\n",phaseData_gsb[i],foldSpectra[i]);
		else
			fprintf(out,"%f\t%f\n",phaseData_gwb[i],foldSpectra[i]);
		
		if(i==0)
			minSpec=foldSpectra[i];
		if(foldSpectra[i]>maxSpec)
			maxSpec=foldSpectra[i];
		if(foldSpectra[i]<minSpec)
			minSpec=foldSpectra[i];
	}

	fclose(out);

        char title[50];	
	
	cpgsch(2);
	cpgsubp(1,3);
	plottwoNoPar(foldFile_GWB,foldFile_GSB,-1);       

	cpgpanl(1,2);
	cpgenv(pstart,pend,minSpec,maxSpec,0,0);
	sprintf(title,"Average Pulse Spectrum");
	cpglab("Phase","Ratio",title);
	if(flag==1)
	        cpgline(corrsize,phaseData_gsb,foldSpectra);
        else
	        cpgline(corrsize,phaseData_gwb,foldSpectra);
	        
	fclose(data_GSB);
	fclose(data_GWB);
	fclose(out);
}

void plotSpectra(char *filename_GSB,char *filename_GWB, char *outname, int pulseNum, int last_gsb, int last_gwb, double s, double *folda_gsb, double *folda_gwb ,double pstart, double pend, double pulseStart, double pulseEnd, double *phase_gsb, double *phase_gwb, int intg_gsb, int intg_gwb, double fold, double t_gsb, double t_gwb, float *readData_gsb, float *readData_gwb, double gsbDelay, double gwbDelay)
{
	FILE *data_GSB, *data_GWB, *out;
	data_GSB=fopen(filename_GSB,"rb");
	data_GWB=fopen(filename_GWB,"rb");
	out=fopen(outname,"w");
	int i,numshift_gsb, numshift_gwb,j;
	int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), skipnum_gwb=(int)((pulseNum-1)*(fold/t_gwb)), sampleDelay_gsb=gsbDelay/t_gsb*1000, sampleDelay_gwb=gwbDelay/t_gwb*1000;
	
	numshift_gsb=(int)(last_gsb*s);
	numshift_gwb=(int)(last_gwb*s);	
	fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);
	fseek(data_GWB,(skipnum_gwb+sampleDelay_gwb+numshift_gwb)*sizeof(float),SEEK_SET);	
	readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB);
	readFloatToDouble(last_gwb,readData_gwb,folda_gwb,data_GWB);
	
	int corrsize_gsb, corrsize_gwb, corrpos=0;
	corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
	corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
	double dpoint=0,ppoint=0;
	float *corrdata_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *corrdata_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	float *phaseData_gsb=(float*)malloc(corrsize_gsb*sizeof(float));
	float *phaseData_gwb=(float*)malloc(corrsize_gwb*sizeof(float));
	
	double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),baseline_gwb=meanBase(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb);
	float maxP_gsb=0,minP_gsb=0,maxP_gwb=0,minP_gwb=0,maxSpec=0,minSpec=0;
	
	for(i=0;i<last_gsb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gsb && j<last_gsb;j++)
		{
			dpoint+=folda_gsb[j];
			ppoint+=phase_gsb[j];
		}
		if(j==last_gsb)
			break;
		ppoint=ppoint/intg_gsb;
		dpoint=dpoint/intg_gsb;
		dpoint=(dpoint-baseline_gsb)/baseline_gsb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
		{
			corrdata_gsb[corrpos]=dpoint;
			phaseData_gsb[corrpos]=ppoint;
			if(corrdata_gsb[corrpos]>maxP_gsb)
				maxP_gsb=corrdata_gsb[corrpos];
			if(corrdata_gsb[corrpos]<minP_gsb)
				minP_gsb=corrdata_gsb[corrpos];	
			corrpos++;
		}
		i+=intg_gsb;
		
	}
	corrpos=0;
	for(i=0;i<last_gwb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gwb && j<last_gwb;j++)
		{
			dpoint+=folda_gwb[j];
			ppoint+=phase_gwb[j];
		}
		if(j==last_gwb)
			break;
		ppoint=ppoint/intg_gwb;
		dpoint=dpoint/intg_gwb;
		dpoint=(dpoint-baseline_gwb)/baseline_gwb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gwb)
		{
			corrdata_gwb[corrpos]=dpoint;
			phaseData_gwb[corrpos]=ppoint;
			if(corrdata_gwb[corrpos]>maxP_gwb)
				maxP_gwb=corrdata_gwb[corrpos];
			if(corrdata_gwb[corrpos]<minP_gwb)
				minP_gwb=corrdata_gwb[corrpos];	
			corrpos++;
		}
		i+=intg_gwb;
	}	
	for(i=0;i<last_gsb;i++)
		folda_gsb[i]=(folda_gsb[i]-baseline_gsb)/baseline_gsb;
	for(i=0;i<last_gwb;i++)
		folda_gwb[i]=(folda_gwb[i]-baseline_gwb)/baseline_gwb;
		
	float stdev_gsb=stdevBaseIntg(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb,intg_gsb),stdev_gwb=stdevBaseIntg(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb,intg_gwb);
	
	float *crossArray1, *crossArray2;
	double t_corr, crossArray1_stdev, crossArray2_stdev; 
	int corrsize,flag;
	if(corrsize_gsb>corrsize_gwb)
	{
		crossArray1=corrdata_gsb;
		crossArray2=(float*)malloc(corrsize_gsb*sizeof(float));
		interpolate_float(corrdata_gwb,corrdata_gsb,t_gwb*intg_gwb,t_gsb*intg_gsb,crossArray2,corrsize_gsb);
		t_corr=t_gsb*intg_gsb;
		corrsize=corrsize_gsb;
		crossArray1_stdev=stdev_gsb;
		crossArray2_stdev=stdev_gwb;
		flag=1;
	}
	else
	{
		crossArray1=corrdata_gwb;
		crossArray2=(float*)malloc(corrsize_gwb*sizeof(float));
		interpolate_float(corrdata_gsb,corrdata_gwb,t_gsb*intg_gsb,t_gwb*intg_gwb,crossArray2,corrsize_gwb);
		t_corr=t_gwb*intg_gwb;
		corrsize=corrsize_gwb;
		crossArray1_stdev=stdev_gwb;
		crossArray2_stdev=stdev_gsb;
		flag=2;
	}
	
	float *specArray=(float*)malloc(corrsize*sizeof(float));
	float specSum=0,specMean=0;
	int numSpecs=0;
	for(i=0;i<corrsize;i++)
	{
		if(crossArray1[i]>spectraPlotThresh*crossArray1_stdev && crossArray2[i]>spectraPlotThresh*crossArray2_stdev)
		{
			if(flag==1)
				specArray[i]=crossArray1[i]/crossArray2[i];
			else
				specArray[i]=crossArray2[i]/crossArray1[i];
				
			specSum+=specArray[i];
			numSpecs++;
		}
		else
			specArray[i]=-1;
	}
	specMean=specSum/numSpecs;
	maxSpec=0;
	for(i=0;i<corrsize;i++)
	{
		if(specArray[i]==-1)
			specArray[i]=specMean;
		if(i==0)
			minSpec=specArray[i];
		
		if(specArray[i]>maxSpec)
			maxSpec=specArray[i];
		if(specArray[i]<minSpec)
			minSpec=specArray[i];
		if(flag==1)
			fprintf(out,"%f\t%f\n",phaseData_gsb[i],specArray[i]);
		else
			fprintf(out,"%f\t%f\n",phaseData_gwb[i],specArray[i]);
	}

        char title[50];	

	cpgsubp(1,3);
        cpgenv(pstart,pend,minP_gwb,maxP_gwb,0,0);
	sprintf(title,"Pulse Number %d at High Freq",pulseNum);
	cpglab("Phase","Intensity",title);
	if(flag==1)
	        cpgline(corrsize,phaseData_gsb,crossArray2);
        else
	        cpgline(corrsize,phaseData_gwb,crossArray1);        
	       
	cpgpanl(1,1);
	cpgenv(pstart,pend,minP_gsb,maxP_gsb,0,0);
	sprintf(title,"Pulse Number %d at Low Freq",pulseNum);
	cpglab("Phase","Intensity",title);
	if(flag==1)
	        cpgline(corrsize,phaseData_gsb,crossArray1);
        else
	        cpgline(corrsize,phaseData_gwb,crossArray2);        
	        
	cpgpanl(1,2);
	cpgenv(pstart,pend,minSpec,maxSpec,0,0);
	sprintf(title,"SPECTRA: Pulse Number %d",pulseNum);
	cpglab("Phase","Ratio",title);
	if(flag==1)
	        cpgline(corrsize,phaseData_gsb,specArray);
        else
	        cpgline(corrsize,phaseData_gwb,specArray);
	        
	fclose(data_GSB);
	fclose(data_GWB);
	fclose(out);
}

void crossCorrelate(char *filename_GSB,char *filename_GWB, char *outname, int pulseNum, int last_gsb, int last_gwb, double s, double *folda_gsb, double *folda_gwb ,double pstart, double pend, double pulseStart, double pulseEnd, double *phase_gsb, double *phase_gwb, int intg_gsb, int intg_gwb, double fold, double t_gsb, double t_gwb, double lag, float *readData_gsb, float *readData_gwb, double gsbDelay, double gwbDelay)
{
	FILE *data_GSB, *data_GWB, *out;
	data_GSB=fopen(filename_GSB,"rb");
	data_GWB=fopen(filename_GWB,"rb");
	out=fopen(outname,"w");
	int i,numshift_gsb, numshift_gwb,j;
	int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), skipnum_gwb=(int)((pulseNum-1)*(fold/t_gwb)), sampleDelay_gsb=gsbDelay/t_gsb*1000, sampleDelay_gwb=gwbDelay/t_gwb*1000;
	
	/*for(i=0;i<(skipnum/last);i++)
	{
		fread(folda,sizeof(double),last,data);
	}
	fread(folda,sizeof(double),(skipnum%last),data);*/

	numshift_gsb=(int)(last_gsb*s);
	numshift_gwb=(int)(last_gwb*s);	
	fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);
	fseek(data_GWB,(skipnum_gwb+sampleDelay_gwb+numshift_gwb)*sizeof(float),SEEK_SET);	
	readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB);
	readFloatToDouble(last_gwb,readData_gwb,folda_gwb,data_GWB);
	
	int corrsize_gsb, corrsize_gwb, corrpos=0;
	corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
	corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
		
	double dpoint=0,ppoint=0;
	double *corrdata_gsb=(double*)malloc(corrsize_gsb*sizeof(double));
	double *corrdata_gwb=(double*)malloc(corrsize_gwb*sizeof(double));
	double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),range_gsb=max(folda_gsb,last_gsb)-baseline_gsb,baseline_gwb=meanBase(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb),range_gwb=max(folda_gwb,last_gwb)-baseline_gwb;
	
	for(i=0;i<last_gsb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gsb && j<last_gsb;j++)
		{
			dpoint+=folda_gsb[j];
			ppoint+=phase_gsb[j];
		}
		if(j==last_gsb)
			break;
		ppoint=ppoint/intg_gsb;
		dpoint=dpoint/intg_gsb;
		dpoint=(dpoint-baseline_gsb)/range_gsb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
		{
			corrdata_gsb[corrpos++]=dpoint;
		}
		i+=intg_gsb;
	}
	corrpos=0;
	for(i=0;i<last_gwb;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_gwb && j<last_gwb;j++)
		{
			dpoint+=folda_gwb[j];
			ppoint+=phase_gwb[j];
		}
		if(j==last_gwb)
			break;
		ppoint=ppoint/intg_gwb;
		dpoint=dpoint/intg_gwb;
		dpoint=(dpoint-baseline_gwb)/range_gwb;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gwb)
		{
			corrdata_gwb[corrpos++]=dpoint;
		}
		i+=intg_gwb;
	}	
	
	double *crossArray1, *crossArray2;
	double t_corr; int corrsize;
	if(corrsize_gsb>corrsize_gwb)
	{
		crossArray1=corrdata_gsb;
		crossArray2=(double*)malloc(corrsize_gsb*sizeof(double));
		interpolate(corrdata_gwb,corrdata_gsb,t_gwb*intg_gwb,t_gsb*intg_gsb,crossArray2,corrsize_gsb);
		t_corr=t_gsb*intg_gsb;
		corrsize=corrsize_gsb;
	}
	else
	{
		crossArray1=corrdata_gwb;
		crossArray2=(double*)malloc(corrsize_gwb*sizeof(double));
		interpolate(corrdata_gsb,corrdata_gwb,t_gsb*intg_gsb,t_gwb*intg_gwb,crossArray2,corrsize_gwb);
		t_corr=t_gwb*intg_gwb;
		corrsize=corrsize_gwb;
	}
	
	int crosssize=2*(int)(lag/(t_corr));
	double *crosscorr=(double*)malloc(crosssize*sizeof(double));
	double sum=0;
	for(i=-crosssize/2;i<crosssize/2;i++)
	{
		sum=0;
		for(j=0;j<corrsize;j++)
		{
			if(j+i>=corrsize)
				sum+=crossArray1[j]*crossArray2[j+i-corrsize];
			else if(j+i<0)
				sum+=crossArray1[j]*crossArray2[j+i+corrsize];
			else
				sum+=crossArray1[j]*crossArray2[j+i];
		}
		crosscorr[i+crosssize/2]=sum;
	}
	
	float lags[crosssize];
	float crossfloat[crosssize];
	float minCross=crosscorr[crosssize-1],maxCross=0;
	for(i=-crosssize/2;i<crosssize/2;i++)
        {	
                fprintf(out,"%.10lf\t%.10lf\n",(double)(t_corr*i),crosscorr[i+crosssize/2]);
        }


        fclose(out);
	
	for(i=-crosssize/2;i<crosssize/2;i++)
	{
		lags[i+crosssize/2]=(float)(t_corr*i);
		crossfloat[i+crosssize/2]=(float)crosscorr[i+crosssize/2];
		if(crossfloat[i+crosssize/2]>maxCross)
			maxCross=crossfloat[i+crosssize/2];
		else if(crossfloat[i+crosssize/2]<minCross)
			minCross=crossfloat[i+crosssize/2];
	}
	
	cpgsci(1);
	cpgenv(lags[0],lags[crosssize-1],minCross,maxCross,0,0);
	char title[50];

	sprintf(title,"Cross-correlation plot for pulse number %d",pulseNum);		

        cpglab("Time lag (ms)","Cross-correlation amplitude",title);
        cpgline(crosssize,lags,crossfloat);

  	fclose(data_GSB);
  	fclose(data_GWB);
}

void produceCrossCorrelate(char *filename_GSB,char *filename_GWB, int pulseNum, int last_gsb, int last_gwb, double s, double *folda_gsb, double *folda_gwb ,double pstart, double pend, double pulseStart, double pulseEnd, double *phase_gsb, double *phase_gwb, int intg_gsb, int intg_gwb, double fold, double t_gsb, double t_gwb, double lag, float *readData_gsb, float *readData_gwb, double gsbDelay, double gwbDelay, double *crosscorr)
{
        FILE *data_GSB, *data_GWB, *out;
        data_GSB=fopen(filename_GSB,"rb");
        data_GWB=fopen(filename_GWB,"rb");
        int i,numshift_gsb, numshift_gwb,j;
        int skipnum_gsb=(int)((pulseNum-1)*(fold/t_gsb)), skipnum_gwb=(int)((pulseNum-1)*(fold/t_gwb)), sampleDelay_gsb=gsbDelay/t_gsb*1000, sampleDelay_gwb=gwbDelay/t_gwb*1000;

	numshift_gsb=(int)(last_gsb*s);
        numshift_gwb=(int)(last_gwb*s);
        fseek(data_GSB,(skipnum_gsb+sampleDelay_gsb+numshift_gsb)*sizeof(float),SEEK_SET);
        fseek(data_GWB,(skipnum_gwb+sampleDelay_gwb+numshift_gwb)*sizeof(float),SEEK_SET);
        readFloatToDouble(last_gsb,readData_gsb,folda_gsb,data_GSB);
        readFloatToDouble(last_gwb,readData_gwb,folda_gwb,data_GWB);

        int corrsize_gsb, corrsize_gwb, corrpos=0;
        corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
        corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));

        double dpoint=0,ppoint=0;
        double *corrdata_gsb=(double*)malloc(corrsize_gsb*sizeof(double));
        double *corrdata_gwb=(double*)malloc(corrsize_gwb*sizeof(double));
        double baseline_gsb=meanBase(folda_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb),range_gsb=max(folda_gsb,last_gsb)-baseline_gsb,baseline_gwb=meanBase(folda_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb),range_gwb=max(folda_gwb,last_gwb)-baseline_gwb;

        for(i=0;i<last_gsb;)
        {
                dpoint=0;
                ppoint=0;
                for(j=i;j<i+intg_gsb && j<last_gsb;j++)
                {
                        dpoint+=folda_gsb[j];
                        ppoint+=phase_gsb[j];
                }
                if(j==last_gsb)
                        break;
                ppoint=ppoint/intg_gsb;
                dpoint=dpoint/intg_gsb;
                dpoint=(dpoint-baseline_gsb)/range_gsb;

                if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gsb)
                {
                        corrdata_gsb[corrpos++]=dpoint;
                }
                i+=intg_gsb;
        }
        corrpos=0;
        for(i=0;i<last_gwb;)
        {
                dpoint=0;
                ppoint=0;
                for(j=i;j<i+intg_gwb && j<last_gwb;j++)
                {
                        dpoint+=folda_gwb[j];
                        ppoint+=phase_gwb[j];
                }
                if(j==last_gwb)
                        break;
                ppoint=ppoint/intg_gwb;
                dpoint=dpoint/intg_gwb;
                dpoint=(dpoint-baseline_gwb)/range_gwb;

                if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_gwb)
                {
                        corrdata_gwb[corrpos++]=dpoint;
                }
                i+=intg_gwb;
        }

        double *crossArray1, *crossArray2;
        double t_corr; int corrsize;
        if(corrsize_gsb>corrsize_gwb)
        {
                crossArray1=corrdata_gsb;
                crossArray2=(double*)malloc(corrsize_gsb*sizeof(double));
                interpolate(corrdata_gwb,corrdata_gsb,t_gwb*intg_gwb,t_gsb*intg_gsb,crossArray2,corrsize_gsb);
                t_corr=t_gsb*intg_gsb;
                corrsize=corrsize_gsb;
        }
        else
        {
                crossArray1=corrdata_gwb;
                crossArray2=(double*)malloc(corrsize_gwb*sizeof(double));
                interpolate(corrdata_gsb,corrdata_gwb,t_gsb*intg_gsb,t_gwb*intg_gwb,crossArray2,corrsize_gwb);
                t_corr=t_gwb*intg_gwb;
                corrsize=corrsize_gwb;
        }

        int crosssize=2*(int)(lag/(t_corr));
        double sum=0;
        for(i=-crosssize/2;i<crosssize/2;i++)
        {
                sum=0;
                for(j=0;j<corrsize;j++)
                {
                        if(j+i>=corrsize)
                                sum+=crossArray1[j]*crossArray2[j+i-corrsize];
                        else if(j+i<0)
                                sum+=crossArray1[j]*crossArray2[j+i+corrsize];
                        else
                                sum+=crossArray1[j]*crossArray2[j+i];
                }
                crosscorr[i+crosssize/2]=sum;
        }

        fclose(data_GSB);
        fclose(data_GWB);
}

void averageCrossCorrelate(char *filename_GSB,char *filename_GWB, int last_gsb, int last_gwb, double s, double *folda_gsb, double *folda_gwb ,double pstart, double pend, double pulseStart, double pulseEnd, double *phase_gsb, double *phase_gwb, int intg_gsb, int intg_gwb, double fold, double t_gsb, double t_gwb, double lag, float *readData_gsb, float *readData_gwb, double gsbDelay, double gwbDelay, char *strongPulse, char *filecorr)
{
        int corrsize_gsb, corrsize_gwb, corrpos=0;
        corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
        corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));
        
        double t_corr; int corrsize;
        if(corrsize_gsb>corrsize_gwb)
        {
                t_corr=t_gsb*intg_gsb;
                corrsize=corrsize_gsb;
        }
        else
        {
                t_corr=t_gwb*intg_gwb;
                corrsize=corrsize_gwb;
        }

        int crosssize=2*(int)(lag/(t_corr));
	int pulseNum,count=0,i=0;
	double SNR;
	double *crosscorr=(double*)malloc(crosssize*sizeof(double));
	double *avgcrosscorr=(double*)malloc(crosssize*sizeof(double));
	
	for(i=0;i<crosssize;i++)
	{
		avgcrosscorr[i]=0;
	}
	
	FILE *listPulse=fopen(strongPulse,"r");
	FILE *out=fopen(filecorr,"w");
	
	while(fscanf(listPulse,"%d %lf",&pulseNum,&SNR)==2)
	{
		printf("\nAdding pulse number %d with SNR %lf",pulseNum, SNR);
		produceCrossCorrelate(filename_GSB,filename_GWB,pulseNum,last_gsb,last_gwb,s,folda_gsb,folda_gwb,pstart,pend,pulseStart,pulseEnd,phase_gsb,phase_gwb,intg_gsb,intg_gwb,  fold,t_gsb,t_gwb,lag,readData_gsb,readData_gwb,gsbDelay,gwbDelay,crosscorr);
		for(i=0;i<crosssize;i++)
		{
			avgcrosscorr[i]+=crosscorr[i];
		}
		count++;
	}
	free(crosscorr);
	for(i=0;i<crosssize;i++)
	{
		avgcrosscorr[i]=avgcrosscorr[i]/count;
	}
	for(i=-crosssize/2;i<crosssize/2;i++)
        {	
                fprintf(out,"%.10lf\t%.10lf\n",(double)(t_corr*i),avgcrosscorr[i+crosssize/2]);
        }
	
	fclose(out);
	float lags[crosssize];
	float crossfloat[crosssize];
	float minCross=avgcrosscorr[0],maxCross=0;
	
	for(i=0;i<crosssize;i++)
	{
		lags[i]=(float)(t_corr*(i-crosssize/2));
		crossfloat[i]=(float)avgcrosscorr[i];
		if(crossfloat[i]>maxCross)
			maxCross=crossfloat[i];
		else if(crossfloat[i]<minCross)
			minCross=crossfloat[i];
	}
	
	cpgsci(1);
	cpgenv(lags[0],lags[crosssize-1],minCross,maxCross,0,0);
	char title[50];
	sprintf(title,"Average crosscorrelation plot");		

        cpglab("Time lag (ms)","Crosscorrelation amplitude",title);
        cpgline(crosssize,lags,crossfloat);
}

void produceCorrelate(char *filename, int pulseNum, int last, double s, double *folda,double pstart, double pend,double *phase,int intg, double fold, double t, double lag,double pulseStart, double pulseEnd, double *autocorr, float *readData, double relDelay)
{
	FILE *data;
	data=fopen(filename,"rb");
	int i,numshift,j;
	int skipnum=(int)((pulseNum-1)*(fold/t)), sampleDelay=relDelay/t*1000;
	
	/*for(i=0;i<(skipnum/last);i++)
	{
		fread(folda,sizeof(double),last,data);
	}
	fread(folda,sizeof(double),(skipnum%last),data);*/
	
	numshift=(int)(last*s);	
	fseek(data,(skipnum+sampleDelay+numshift)*sizeof(float),SEEK_SET);
	readFloatToDouble(last,readData,folda,data);
	
	double baseline=meanBase(folda,phase,pulseStart,pulseEnd,last),range=max(folda,last)-baseline, dpoint=0,ppoint=0;
	int corrsize=(int)(last/intg*(pend-pstart)),corrpos=0;
	double *corrdata=(double*)malloc(corrsize*sizeof(double));
	
	for(i=0;i<last;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg && j<last;j++)
		{
			dpoint+=folda[j];
			ppoint+=phase[j];
		}
		if(j==last)
			break;
		ppoint=ppoint/intg;
		dpoint=dpoint/intg;
		dpoint=(dpoint-baseline)/range;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize)
		{
			corrdata[corrpos++]=dpoint;
		}
		i+=intg;
	}
	
	int autosize=(int)(lag/(t*intg));
	
	double sum=0;
	for(i=0;i<autosize;i++)
	{
		sum=0;
		for(j=0;j<corrsize;j++)
		{
			if(j+i>=corrsize)
			{
				sum+=corrdata[j]*corrdata[j+i-corrsize];
			}
			else
			{
				sum+=corrdata[j]*corrdata[j+i];
			}
		}
		
		autocorr[i]=sum;
	}
	
  	fclose(data);
}

void averageAutoCorrelate(char *filename, int last, double s, double *folda,double pstart, double pend,double *phase,int intg, double fold, double t, double lag, double pulseStart, double pulseEnd, char *strongPulse,char *filecorr, float *readData, double relDelay, int freqFlag)
{
	int autosize=(int)(lag/(t*intg)),pulseNum,count=0,i=0;
	double SNR;
	double *autocorr=(double*)malloc(autosize*sizeof(double));
	double *avgcorr=(double*)malloc(autosize*sizeof(double));
	
	for(i=0;i<autosize;i++)
	{
		avgcorr[i]=0;
	}
	
	FILE *listPulse=fopen(strongPulse,"r");
	FILE *out=fopen(filecorr,"w");
	
	while(fscanf(listPulse,"%d %lf",&pulseNum,&SNR)==2)
	{
		printf("\nAdding pulse number %d with SNR %lf",pulseNum, SNR);
		produceCorrelate(filename,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,lag,pulseStart,pulseEnd,autocorr,readData,relDelay);
		for(i=0;i<autosize;i++)
		{
			avgcorr[i]+=autocorr[i];
		}
		count++;
	}
	free(autocorr);
	for(i=0;i<autosize;i++)
	{
		avgcorr[i]=avgcorr[i]/count;
	}
	for(i=0;i<autosize;i++)
        {	
                fprintf(out,"%.10lf\t%.10lf\n",(double)(t*intg*i),avgcorr[i]);
        }
	fclose(out);
	float lags[autosize];
	float autofloat[autosize];
	float minAuto=avgcorr[autosize-1],maxAuto=0;
	
	for(i=0;i<autosize;i++)
	{
		lags[i]=(float)(t*intg*i);
		autofloat[i]=(float)avgcorr[i];
		if(autofloat[i]>maxAuto)
			maxAuto=autofloat[i];
		else if(autofloat[i]<minAuto)
			minAuto=autofloat[i];
	}
	
	cpgsci(1);
	cpgenv(lags[0],lags[autosize-1],minAuto,maxAuto,0,0);
	char title[50];
	if(freqFlag==0)
		sprintf(title,"Average autocorrelation plot at High Freq");		
	else
		sprintf(title,"Average autocorrelation plot at Low Freq");		

        cpglab("Time lag (ms)","Autocorrelation amplitude",title);
        cpgline(autosize,lags,autofloat);
}

double getAvgHalfWidth(char *avgAutocorrFile, double t, float lag, int intg)
{
	FILE *avgRead=fopen(avgAutocorrFile,"r");
	int autosize=(int)(lag/(t*intg)),i;
	float *avgAutocorr=(float*)malloc(autosize*sizeof(float));
	float *lags=(float*)malloc(autosize*sizeof(float));
	float slope, refWidth=0.5, width=0;
	
	if(avgRead==NULL)
	{
		printf("\nERROR! AvgAutocorr file now found!\n");
		return;
	}
	
	for(i=0;i<autosize;i++)
		fscanf(avgRead,"%f %f",lags+i,avgAutocorr+i);
		
	float refCorr=refWidth*avgAutocorr[1];
	for(i=0;i<autosize;i++)
	{
		if(avgAutocorr[i]<refCorr)
		{
			slope=(avgAutocorr[i-1]-avgAutocorr[i])/(lags[i-1]-lags[i]);
			width=(refCorr-avgAutocorr[i])/slope+lags[i];
			break;
		}
	}

	if(width==0)
		printf("\nERROR: REFERENCE POINT NOT FOUND!!\n");
	else
		printf("\nReference width found at lag of %f ms for reference point of %f\n",width,refCorr);
		
	return (double)width;
}

void showStrongPulses(char *filename,int last,double s, double *folda, double *phase,int intg, double fold, double t, double pulseStart, double pulseEnd,double minSNR,char *outname, float *readData, double relDelay)
{
	FILE *data, *out;
	data=fopen(filename,"rb");
	out=fopen(outname,"w");
	int i,numshift=(int)(last*s),pulseNum=1,skipnum=0,sampleDelay;
	double SNR=0;
	sampleDelay=relDelay/t*1000;

	while(1)
	{
		skipnum=(int)((pulseNum-1)*(fold/t));
		fseek(data,(skipnum+sampleDelay+numshift)*sizeof(float),SEEK_SET);
		if(readFloatToDouble(last,readData,folda,data)!=1)
			break;
		
		SNR=filterStrong(folda,phase,pulseStart,pulseEnd,last,intg);
		if (SNR>=minSNR)
		{
			printf("\nPulse Number %d with SNR %f",pulseNum,SNR);
			fprintf(out,"%d\t%f\n",pulseNum,SNR);
		}
		pulseNum++;
	}
	fclose(out);
	fclose(data);
}

void showStrongPulsesPeak(char *filename,int last,double s, double *folda, double *phase,int intg, double fold, double t, double pulseStart, double pulseEnd,double minSNR,char *outname,float *readData,double relDelay)
{
	FILE *data, *out;
	data=fopen(filename,"rb");
	out=fopen(outname,"w");
	int i,numshift=(int)(last*s),pulseNum=1,skipnum=0,sampleDelay;
	double SNR=0;
	sampleDelay=relDelay/t*1000;

	while(1)
	{
		skipnum=(int)((pulseNum-1)*(fold/t));
		fseek(data,(skipnum+sampleDelay+numshift)*sizeof(float),SEEK_SET);
		if(readFloatToDouble(last,readData,folda,data)!=1)
			break;
		
		SNR=filterStrongPeak(folda,phase,pulseStart,pulseEnd,last);
		if (SNR>=minSNR)
		{
			printf("\nPulse Number %d with SNR %f",pulseNum,SNR);
			fprintf(out,"%d\t%f\n",pulseNum,SNR);
		}
		pulseNum++;
	}
	fclose(out);
	fclose(data);
}

/*void produceADP(double *ADPwindow, int numADP, double tsamp, char *outfile, int pulseNum, int panlflag)
{
	int minExpo=(int)floor(log((double)numADP)/log(2))+1,fftsize=(int)pow((double)2,minExpo),i=0;
	double Nyqfreq=(double)1/(2*tsamp);
	printf("\nSample size is %d. FFT size is %d",numADP,fftsize);
	double *ADPuse=(double*)malloc(fftsize*sizeof(double));
	for(i=0;i<fftsize;i++)
	{
		if (i<numADP)
			ADPuse[i]=ADPwindow[i];
		else
			ADPuse[i]=0;
	}
	
	fftw_complex spec[fftsize/2+1];
	fftw_plan p;
	p=fftw_plan_dft_r2c_1d(fftsize,ADPuse,spec,FFTW_ESTIMATE);
	fftw_execute_dft_r2c(p,ADPuse,spec);
	
	float *specInt=(float*)malloc((fftsize/2+1)*sizeof(float));
	float *freq=(float*)malloc((fftsize/2+1)*sizeof(float));
	float *ADP=(float*)malloc((fftsize/2+1)*sizeof(float));
	
	for(i=0;i<fftsize/2+1;i++)
	{
		specInt[i]=(float)(spec[i][0]*spec[i][0]+spec[i][1]*spec[i][1]);
		freq[i]=(float)2*i/fftsize*Nyqfreq;
		ADP[i]=freq[i]*freq[i]*specInt[i]*specInt[i];
	}
	
	FILE *out=fopen("outfile","w");
	for(i=0;i<fftsize/2+1;i++)
	{
		fprintf(out,"%f\t%f\n",freq[i],ADP[i]);
	}
	
	fclose(out);

	float minADP=ADP[0],maxADP=0;
	
	for(i=0;i<fftsize/2+1;i++)
	{
		if(ADP[i]>maxADP)
			maxADP=ADP[i];
		if(ADP[i]<minADP)
			minADP=ADP[i];
	}
	
	cpgpanl(1,panlflag);
	cpgenv(freq[0],freq[fftsize/2],minADP,maxADP,0,0);
	char title[50];
	if (panlflag==1)
		sprintf(title,"ADP for pulse number %d in L pol",pulseNum);		
	else
		sprintf(title,"ADP for pulse number %d in R pol",pulseNum);		

        cpglab("Frequency (kHz)","ADP",title);
        cpgline(fftsize/2,freq,ADP);
	fftw_destroy_plan(p);
}

void ADPcreatePulse(char *filename, char *outfile,int pulseNum, int last,double s,double *folda,double pstart, double pend,double *phase,int intg, double fold, double t, double pulseStart, double pulseEnd, int panlflag)
{
	FILE *data;
	data=fopen(filename,"rb");
	int i,numshift,j,pos=0;
	int skipnum=(int)((pulseNum-1)*(fold/t)),numADP=(int)((pend-pstart)*last/intg);
	
	fseek(data,skipnum*sizeof(double),SEEK_SET);
	numshift=(int)(last*s);
	fread(folda,sizeof(double),numshift,data);
	fread(folda,sizeof(double),last,data);

	double baseline=meanBase(folda,phase,pulseStart,pulseEnd,last),range=max(folda,last)-baseline, dpoint=0,ppoint=0;
	double *ADPwindow=(double*)malloc(numADP*sizeof(double));
	
	for(i=0;i<last;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg && j<last;j++)
		{
			dpoint+=folda[j];
			ppoint+=phase[j];
		}
		if(j==last)
			break;
		ppoint=ppoint/intg;
		dpoint=dpoint/intg;
		dpoint=(dpoint-baseline)/range;
		
		if(ppoint>=pstart && ppoint<=pend && pos<numADP)
		{
			ADPwindow[pos]=dpoint;
			pos++;
		}
		i+=intg;
	}
	
	produceADP(ADPwindow,numADP,t*intg,outfile,pulseNum,panlflag);
}*/

void writefile(char *filename, char *outfile,int pulseNum, int last,double s,double *folda,double pstart, double pend,double *phase,int intg, double fold, double t, double pulseStart, double pulseEnd, float *readData,double relDelay)
{
	FILE *data,*out;
	data=fopen(filename,"rb");
	out=fopen(outfile,"w");
	int i,numshift,j;
	int skipnum=(int)((pulseNum-1)*(fold/t)), sampleDelay=relDelay/t*1000;
	
	/*for(i=0;i<(skipnum/last);i++)
	{
		fread(folda,sizeof(double),last,data);
	}
	fread(folda,sizeof(double),(skipnum%last),data);*/
	
	numshift=(int)(last*s);
	fseek(data,(skipnum+sampleDelay+numshift)*sizeof(float),SEEK_SET);
	if(readFloatToDouble(last,readData,folda,data)!=1)
	{
		printf("\nReached end of file!\n");
		return;
	}
		

	double baseline=meanBase(folda,phase,pulseStart,pulseEnd,last),range=max(folda,last)-baseline, dpoint=0,ppoint=0;	
	
	for(i=0;i<last;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg && j<last;j++)
		{
			dpoint+=folda[j];
			ppoint+=phase[j];
		}
		if(j==last)
			break;
		ppoint=ppoint/intg;
		dpoint=dpoint/intg;
		dpoint=(dpoint-baseline)/range;
		
		if(ppoint>=pstart && ppoint<=pend)
		{
			fprintf(out,"%.10lf\t%.10lf\n",ppoint,dpoint);
		}
		i+=intg;
	}
	fclose(data);
	fclose(out);
}

void foldprofile(char* filename, char *foldname, double fold, double t, double *phase,double pstart,double pend, int numshift,double *profile, int option, int intg, double pulseStart, double pulseEnd, double relDelay)
{
	int i,last;
	if (floor((double)fold/t)==fold/t)
	{
		last=floor(fold/t);
	}
	else
	{
		last=floor(fold/t)+1;
	}
	if(option==1)
	{
		FILE *data,*output;
		data=fopen(filename,"rb");
		int sampleDelay=relDelay/t*1000;
		fseek(data,sampleDelay*sizeof(float),SEEK_SET);
		output=fopen(foldname,"w");
		int sample=0,b=0,sampleNum=32768;
	
		printf("\nFolding profile. Please wait.\n");
	
		double *folda=(double*)malloc(last*sizeof(double));

		int *count=(int*)malloc(last*sizeof(int));
		
		double store[sampleNum];
		float readBuffer[sampleNum];



		for(i=0;i<last;i++)
		{
			count[i]=0;
			folda[i]=0;
		}
		if (data!=NULL)
		{
			while(readFloatToDouble(sampleNum,readBuffer,store,data)==1)
			{	
				for(i=0;i<sampleNum;i++)
				{
					b=floor((sample*t-fold*floor((double)sample*t/fold))/t);
					sample++;
					folda[b]+=(store[i]-folda[b])/(count[b]+1);				
					count[b]++;
				}
			}
			
			for(i=0;i<last;i++)
			{	
				fprintf(output,"%.10lf\t%.10lf\n",phase[i],folda[i]);
			}
			memcpy(profile,folda,last*sizeof(double));
			fclose(output);
			fclose(data);
		}
		else
		{
			printf("\nFile not accessible\n");
		}	
	}
	else if(option==2)
	{
		FILE *output;
		int j;
		output=fopen(foldname,"w");		
		double dpoint=0,ppoint=0,foldBase=meanBase(profile,phase,pulseStart,pulseEnd,last),foldRange=max(profile,last)-foldBase;
		for(i=0;i<last;)
		{
			dpoint=0;
			ppoint=0;
			for(j=i;j<i+intg && j<last;j++)
			{
				dpoint+=profile[j];
				ppoint+=phase[j];
			}
			
			if(j==last)
				break;
				
			ppoint=ppoint/intg;
			dpoint=dpoint/intg;
			dpoint=(dpoint-foldBase)/foldRange;
			
			if(ppoint>=pstart && ppoint<=pend)
			{
				fprintf(output,"%.10lf\t%.10lf\n",ppoint,dpoint);
			}
			i+=intg;
		}
		fclose(output);
	}
	else if(option==3)
	{		
		if(numshift>=0)
		{	
			double *store=(double*)malloc(numshift*sizeof(double));
			for(i=0;i<numshift;i++)
			{
				store[i]=profile[i];
			}
			for(i=numshift;i<last;i++)
			{
				profile[i-numshift]=profile[i];
			}
			for(i=0;i<numshift;i++)
			{
				profile[i+last-numshift]=store[i];
			}
		}
		else
		{
			numshift=-numshift;
			double *store=(double*)malloc(numshift*sizeof(double));
			for(i=last-numshift;i<last;i++)
			{
				store[i-last+numshift]=profile[i];
			}
			for(i=last-1;i>=numshift;i--)
			{
				profile[i]=profile[i-numshift];
			}
			for(i=0;i<numshift;i++)
			{
				profile[i]=store[i];
			}
		}
	}
}

int recordMicro(char *filename_GSB,char *outfile_GSB, char *filename_GWB, char *outfile_GWB, int pulseNum,int last_gsb,int last_gwb, double s,double* folda_gsb,double *folda_gwb, double pstart,double pend,double* phase_gsb,double *phase_gwb,int intg_gsb,int intg_gwb, double fold,double t_gsb,double t_gwb, double pulseStart,double pulseEnd,float* readData_gsb,float *readData_gwb,double gsbDelay,double gwbDelay, FILE* recordFile,int pulseWindow,char *foldfile_GSB,char *foldfile_GWB,double *profile_gsb,double* profile_gwb,int numshift)
{
	char *input=(char*)malloc(100*sizeof(char));
	writefile(filename_GSB,outfile_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,readData_gsb,gsbDelay);
	writefile(filename_GWB,outfile_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,readData_gwb,gwbDelay);	

	foldprofile(filename_GSB,foldfile_GSB,fold,t_gsb,phase_gsb,pstart,pend,numshift,profile_gsb,2,intg_gsb,pulseStart,pulseEnd,gsbDelay);
	foldprofile(filename_GWB,foldfile_GWB,fold,t_gwb,phase_gwb,pstart,pend,numshift,profile_gwb,2,intg_gwb,pulseStart,pulseEnd,gwbDelay);
			
	cpgslct(pulseWindow);
	cpgsch(1.5);
	plottwo(outfile_GWB,outfile_GSB,pulseNum);
	printf("\nIs this pulse usable? y for yes, n for no, e for end check:");
        fgets(input,100,stdin);
	int valInput=0;

	while(valInput==0)
	{
		if(input[0]=='y')
		{
			fprintf(recordFile,"%d\t1\n",pulseNum);
			valInput=1;
		}
		else if(input[0]=='e')
		{
			return -1;
		}
		else if(input[0]=='n')
		{
			fprintf(recordFile,"%d\t-1\n",pulseNum);
			valInput=1;
		}
		
		else
		{
			printf("Invalid input. try again.\n");
		}
	}
	free(input);
	return 0;
}

int foldedCrossCorr(char *fold_GSB,char *fold_GWB,double pstart,double pend,double lag,double t_gsb,double t_gwb,int intg_gsb,int intg_gwb,double fold,int last_gsb,int last_gwb)
{
	int i,j;
	FILE *avgprofile_GSB=fopen(fold_GSB, "r");
	FILE *avgprofile_GWB=fopen(fold_GWB, "r");
	//printf("Yo");
	if (avgprofile_GSB==NULL)
	{
		printf("No file found for GSB");
		return 1;
	}
	if (avgprofile_GWB==NULL)
	{
		printf("No file found foe GWB");
		return 1;
	}
	//printf("%lf %lf %lf %lf",pstart,pend,t_gsb,t_gwb);
	int corrsize_gsb, corrsize_gwb, corrpos=0;
	corrsize_gwb=(int)(last_gwb/intg_gwb*(pend-pstart));
	corrsize_gsb=(int)(last_gsb/intg_gsb*(pend-pstart));	
	//printf("%d %d",num_gsb,num_gwb);
	double *corrdata_gsb=(double*)(malloc(corrsize_gsb*sizeof(double)));
	double *phase_gsb=(double*)(malloc(corrsize_gsb*sizeof(double)));
	double *corrdata_gwb=(double*)(malloc(corrsize_gwb*sizeof(double)));
	double *phase_gwb=(double*)(malloc(corrsize_gwb*sizeof(double)));

	char singleLine[100];
	int num=0;

	while (!feof(avgprofile_GSB) && (num<=corrsize_gsb))
		{
			fgets(singleLine,100,avgprofile_GSB);		// 
			sscanf(singleLine, "%lf\t%lf", &phase_gsb[num],&corrdata_gsb[num]);	
			num++;
		}

	num=0;
	
	while (!feof(avgprofile_GWB) && (num<=corrsize_gwb))
		{
			fgets(singleLine,100,avgprofile_GWB);		// 
			sscanf(singleLine, "%lf\t%lf", &phase_gwb[num],&corrdata_gwb[num]);	
			num++;
		}

	/*for (i=0;i<corrsize_gwb;i++)
	{
		printf("%lf\n",phase_gwb[i]);
	}*/

	double *crossArray1, *crossArray2;
	double t_corr; int corrsize;
	//Finding the frequency which has higher density of data points with/without integration so that the other frequency can be brought to a similar resolution
	if(corrsize_gsb>corrsize_gwb)
	{
		crossArray1=corrdata_gsb;
		crossArray2=(double*)malloc(corrsize_gsb*sizeof(double));
		interpolate(corrdata_gwb,corrdata_gsb,t_gwb*intg_gwb,t_gsb*intg_gsb,crossArray2,corrsize_gsb);
		t_corr=t_gsb*intg_gsb;
		corrsize=corrsize_gsb;
	}
	else
	{
		crossArray2=corrdata_gwb;
		crossArray1=(double*)malloc(corrsize_gwb*sizeof(double));
		interpolate(corrdata_gsb,corrdata_gwb,t_gsb*intg_gsb,t_gwb*intg_gwb,crossArray1,corrsize_gwb);
		t_corr=t_gwb*intg_gwb;
		corrsize=corrsize_gwb;
	}


	/*for(i=0;i<corrsize;i++)
	{
		printf("%lf\n",crossArray1[i]);
	} */
			
			//Now crossArray1 and 2 containf the identical resolution data which can now be cross correlated
	int crosssize=2*(int)((lag*1.0)/(t_corr*1.0));
	double *crosscorr=(double*)malloc(crosssize*sizeof(double));
	double sum1=0;
	int c=0;  //Performing crosscorrelation 
	for(i=-crosssize/2;i<crosssize/2;i++)
	{
		sum1=0;
		for(j=0;j<corrsize;j++)
		{
			if(j+i>=corrsize)
				sum1+=crossArray1[j]*crossArray2[j+i-corrsize];
			else if(j+i<0)
				sum1+=crossArray1[j]*crossArray2[j+i+corrsize];
			else
				sum1+=crossArray1[j]*crossArray2[j+i];	
		}
		crosscorr[i+crosssize/2]=sum1;
		sum1=0;
	}

	
	double lags[crosssize];
	double crossdouble[crosssize];
	double minCross=crosscorr[crosssize-1],maxCross=0;


	for(i=-crosssize/2;i<crosssize/2;i++)
	{
		lags[i+crosssize/2]=(float)(t_corr*i);
		crossdouble[i+crosssize/2]=(float)crosscorr[i+crosssize/2];
		if(crossdouble[i+crosssize/2]>maxCross)
		{
			maxCross=crossdouble[i+crosssize/2];
		}
		else if(crossdouble[i+crosssize/2]<minCross)
			minCross=crossdouble[i+crosssize/2];
	}

	float *crossfloat=(float*)(malloc(crosssize*sizeof(float)));
	float lagsfloat[crosssize];
	float minCrossfloat=(float)minCross,maxCrossfloat=(float)(maxCross);

	for (i=0;i<crosssize;i++)
	{
		crossfloat[i]=(float)crossdouble[i];
		lagsfloat[i]=(float)(lags[i]);
	}


	for (i=0;i<crosssize;i++)
	{
		crossfloat[i]=(crossfloat[i]*1.0)/(maxCrossfloat*1.0);
		//printf("%lf\n",crossfloat[i]);
	}

	float maxpos=0;
	float maxval=0;

	for(i=0;i<crosssize;i++)
	{
		if(crossfloat[i]>maxval)
		{
			maxval=crossfloat[i];
			maxpos=lagsfloat[i];
		}
	}				
			
	printf("\nCross-correlation peak at %f ms",maxpos);
	

	cpgsci(1);
	cpgenv(lagsfloat[0],lagsfloat[crosssize-1],0,1,0,0);
	char title[50];
	sprintf(title,"Folded-profile crosscorrelation plot");		

        cpglab("Time lag (ms)","Crosscorrelation amplitude",title);
        cpgline(crosssize,lagsfloat,crossfloat);
			

}




int main(int argc, char* argv[])
{

  int option;
  printf("\nSingle frequency=1 or dual frequency=2\n");
  scanf("%d",&option);
	switch(option)
	{
	case 2:
	{
        double t_gsb,t_gwb,fold,value,time,k,lag=0,DM,smoothDur_gsb=0,smoothDur_gwb=0;
	int last_gsb,last_gwb,i,sample=0,count=0,scanPulseNum=0;
	int b=0;
	char *input=(char*)malloc(100*sizeof(char));

        int pulsarnum;
        printf("B0329+54_1 ---- press 1\nB0329+54_2 ---- press 2\nB0628-28_1 ---- press 3\nB0628-28_2 ---- press 4\nB0740-28_1 ---- press 5\nB0950+08 ---- press 6\nB1642-03_1 ---- press 7\nB1642-03_2 ---- press 8\nB1749-28_2 ---- press 9\nB2016+28 ---- press 10\nB2016+28_1 ---- press 11\nB2045 ---- press 12\nB1642-03 half bands with adjusted DM ---- press 13\nB1642-03 Band3 ch.#64-192 (25 MHz pieces) ---- press 14\nB1642-03 Band3 ch.#192-320 (25 MHz pieces) ---- press 15\nB1642-03 Band3 ch.#320-448 (25 MHz pieces) ---- press 16\nB1642-03 25 MHz Sub 7&6 ---- press 17\nB1642-03 25 MHz Sub 7&5 ---- press 18\nB1642-03 25 MHz Sub 7&4 ---- press 19\nB1642-03 25 MHz Sub 7&3 ---- press 20\nB1642-03 25 MHz Sub 7&2 ---- press 21\nJ2145 ---- press 22\n");
        scanf ("%d", &pulsarnum);

        if (pulsarnum > 22 || pulsarnum < 1)
            {
              printf ("Wrong input\n");
              return 1;
            }

        

        int in;
        char gsbfile[1000],gwbfile[1000],timegsb[1000],timegwb[1000],file[1000];
	float gsbreso,gwbreso,pulseper,dm,timeOffset;
	FILE *parameter = fopen("/Data/ygupta/psr_microstr_work/recovered_code/dual_freq_analysis/Parameters2.txt", "r");
	
	if (parameter == NULL) 
        {
        perror("Error: Failed to open file.");
        return 1;
        }

	for (in = 0; in < pulsarnum; in++)
           {
            if (fscanf(parameter, "%*s %f %f %f %f %f %s %s %s %s",&timeOffset,&gsbreso,&gwbreso,&pulseper,&dm,gsbfile,gwbfile,timegsb,timegwb) == 1);
            }
        
	//printf("\nGSB file is %s\n",gsbfile);
	printf("\nGSB resolution is %f",gsbreso);	
        printf("\nGWB resolution %f",gwbreso);
        printf("\nPulse period is %f",pulseper);
        printf("\nDM is %f",dm);
        printf("\nTime Offset is %f\n",timeOffset);
        
	char singleLine[100];
	int index=-1;
	int pulseList[10000];
	int totalpulses=0;
        FILE *pulsesFile;
	int pulseFromFile=0,checkForFile=0;
	printf ("\nIf you want to read pulses from a file press 1,else press 0\n");
	scanf ("%d",&checkForFile);
	if (checkForFile==1)
           {
		pulseFromFile=1;
		printf("\nGive file path\n");
		scanf("%s",file);
		//printf("%s",file);
		pulsesFile=fopen(file,"r");
		while (!feof(pulsesFile))
		{
			fgets(singleLine,100,pulsesFile);
			sscanf(singleLine, "%d", &pulseList[totalpulses]);	
			totalpulses++;
		}
		totalpulses-=1;
		printf("\nTotal number of pulses in file = %d\n",totalpulses);
		//int i;
		//for (i=0;i<totalpulses;i++)
			//printf("\n%d",pulseList[i]);
	    }
         else if (checkForFile!=0) 
         {
		printf("\nWrong option chosen.\n");
		return 1;
	 }
	
	t_gsb=gsbreso;
	t_gwb=gwbreso;
	fold=pulseper;
	DM=dm;

	char *filename_GSB=gsbfile;
	char *filename_GWB=gwbfile;
	char *outfile_GSB="plot_GSB.txt";
	char *outfile_GWB="plot_GWB.txt";
	char *foldfile_GSB="fold_GSB.txt";
	char *foldfile_GWB="fold_GWB.txt";
	char *outcorr_GSB="corr_GSB.txt";
	char *outcorr_GWB="corr_GWB.txt";
	char *avgcorr_GSB="avgcorr_GSB.txt";
	char *avgcorr_GWB="avgcorr_GWB.txt";
	char *strong_GSB="strongpulse_GSB.txt";
	char *strong_GWB="strongpulse_GWB.txt";
	char *adp_GSB="ADP_GSB.txt";
	char *adp_GWB="ADP_GWB.txt";
	char *avg_crossCorr="avgCrossCorr.txt";
	char *crossCorr="crossCorr.txt";
	//Vincross
	char *resicrosscorr="Residual_crosscorr.txt";
	char *subAutoFile="subAuto.txt";
	char *FFTName="SpectrumFile.txt";
	char *microFileName="microStat.txt";
	char *recordFileName="recordStat.txt";
	char *specFile="specPlot.txt";
	char *avgPulseSpec="avgSpec.txt";
	char *SG_coeffs="SG_coeffs.txt";
	char *visualMicroFile="Vis_micropulse.txt";

	/*double t_gsb,t_gwb,fold,value,time,k,lag=0,DM,smoothDur_gsb=0,smoothDur_gwb=0;
	int last_gsb,last_gwb,i,sample=0,count=0,scanPulseNum=0;
	int b=0;
	char *input=(char*)malloc(100*sizeof(char));*/

        /*t_gsb=(double)atof(argv[3]);
	t_gwb=(double)atof(argv[4]);
	fold=(double)atof(argv[5]);
	DM=(double)atof(argv[6]);*/
	
	if (floor((double)fold/t_gsb)==fold/t_gsb)
		last_gsb=floor(fold/t_gsb);
	else
		last_gsb=floor(fold/t_gsb)+1;
	
	if (floor((double)fold/t_gwb)==fold/t_gwb)
                last_gwb=floor(fold/t_gwb);
        else
                last_gwb=floor(fold/t_gwb)+1;
        FILE *GWBStamp=fopen(timegwb,"r");
        FILE *GSBStamp=fopen(timegsb,"r");

	if(GWBStamp==NULL)
	{
		printf("\nGWB timestamp not found! Exit.\n");
		exit(-1);
	}
	if(GSBStamp==NULL)
	{
		printf("\nGSB timestamp not found! Exit.\n");
		exit(-1);
	}
	
	float cutOffhighfreq=0,cutOfflowfreq=0;
        char tempStore[100];
        int gwb_hour, gwb_min, gsb_hour, gsb_min, gsb_sec, temp;
        double gwb_sec, gsb_fracSec, gsb_fullSec;
        fscanf(GWBStamp,"%s %s %s %s\n",tempStore,tempStore,tempStore,tempStore);
        fscanf(GWBStamp,"IST Time: %d:%d:%lf",&gwb_hour,&gwb_min,&gwb_sec);
        fscanf(GSBStamp,"%d %d %d %d %d %d %lf",&temp,&temp,&temp,&gsb_hour,&gsb_min,&gsb_sec,&gsb_fracSec);
	gsb_fullSec=gsb_sec+gsb_fracSec;
	//double relDiff=diffTime(gwb_hour,gwb_min,gwb_sec,gsb_hour,gsb_min,gsb_fullSec);
	double relDiff=0;
	//double relDelay=calculateRelDelay(relDiff,DM,timeOffset);
	double relDelay=0;
	printf("\nGWB timestamp is %02d:%02d:%2.9lf.\nGSB timestamp is %02d:%02d:%2.9lf.\n",gwb_hour,gwb_min,gwb_sec,gsb_hour,gsb_min,gsb_fullSec);
	printf("\nTime difference between timestamps is %lf secs (GSB ahead of GWB)\nRelative delay to be adjusted is %lf s.\n",relDiff,relDelay);

	double *folda_gsb=(double*)malloc(last_gsb*sizeof(double));
	double *folda_gwb=(double*)malloc(last_gwb*sizeof(double));
	float *readData_gsb=(float*)malloc(last_gsb*sizeof(float));
	float *readData_gwb=(float*)malloc(last_gwb*sizeof(float));
	double *phase_gsb=(double*)malloc(last_gsb*sizeof(double));
	double *phase_gwb=(double*)malloc(last_gwb*sizeof(double));
	double *profile_gsb=(double*)malloc(last_gsb*sizeof(double));	
	double *profile_gwb=(double*)malloc(last_gwb*sizeof(double));
	
	//printf("\nPulse Viewer\nInstructions:\nUse 'n' to go to next pulse\nUse 'b' to go the earlier pulse\nUse 'p <start phase> <end phase>' to set a phase range (Default is full period)\nUse 's <phase shift>' to shift phase window (Phase shift should remain between 0 and 1 !)\nUse 'f' to show folded profile\nUse 'h' to hide folded profile\nUse 'i <Integration>' to integrate profile\nUse 'g <Pulse Number>' to get profile for a specified pulse number\nUse 'a <Time lag>' to get autocorrelation plot for the current pulse, integration and phase range\nUse 'c <MinSNR>' to get list of pulse numbers above a given threshold\nUse 't <MinSNR>' to get list of pulse numbers with peak above a given SNR\nUse 'v <Time lag(ms)> to get the average ACF for the strongest pulses\nUse 'y' to view the ADP for the specified window\nUse 'k' to get instructions\nUse 'q' to terminate program\n");
	int pulseNum=0,numshift=0,plotswitch=0,intg_gsb=1,intg_gwb=1,curshift=0,backOp,visualPulseindex=0,val=0;
	double s=0,pstart=0,pend=1,pulseStart=0,pulseEnd=1,minSNR=0;
	char *pch;
	struct  visualPulse visarr[1000];
	int autoflag=0, autoWindow, psWindow, foldWindow;	
	int sg_coeff_size=0;
	int pulseWindow=cpgopen("/xs");
	cpgask(0);
	
	double width_gsb, width_gwb;

	double gsbDelay, gwbDelay;
	if(relDelay>=0)
	{
		gsbDelay=0;
		gwbDelay=relDelay;
	}
	else
	{
		gsbDelay=-relDelay;
		gwbDelay=0;
	}
	
	for(i=0;i<last_gsb;i++)
	{
		phase_gsb[i]=(double)i/last_gsb;
	}
	for(i=0;i<last_gwb;i++)
	{
                phase_gwb[i]=(double)i/last_gwb;
	}
	
	foldprofile(filename_GSB,foldfile_GSB,fold,t_gsb,phase_gsb,pstart,pend,numshift,profile_gsb,1,intg_gsb,pulseStart,pulseEnd,gsbDelay);
	foldprofile(filename_GWB,foldfile_GWB,fold,t_gwb,phase_gwb,pstart,pend,numshift,profile_gwb,1,intg_gwb,pulseStart,pulseEnd,gwbDelay);
	
	printf("\nGSB default integration is %d samples. GWB default integration is %d samples\n",intg_gsb,intg_gwb);
	
	printf("\nPlease check visually and specify pulse on window for baseline and noise statistics. Showing the two folded profiles.\n");

	cpgslct(pulseWindow);
	cpgsch(1.5);
	plottwo(foldfile_GWB,foldfile_GSB,-1);
	
	printf("\nEnter pulse start:");
	scanf("%lf",&pulseStart);
	printf("\nEnter pulse end:");
	scanf("%lf",&pulseEnd);
	

	FILE* vismicro=fopen(visualMicroFile,"w");

	while(1)
	{
		printf("\nEnter option:");
		fgets(input,100,stdin);

		if (input[0]=='l')
		{

			if (visualPulseindex==0)
			{
				fprintf(vismicro,"pulseNum\thighfreq\tlowfreq\n");			
			}

			printf("\nDo you see micropulse at low frequency? (1=yes,-1=no,0=o)\n");
			scanf("%d",&val);
			if (val==1)
				visarr[visualPulseindex].microatlowfreq=1;
			else if (val==-1)
				visarr[visualPulseindex].microatlowfreq=-1;
			else if (val==0)
				visarr[visualPulseindex].microatlowfreq=0;
			else
			{
				printf("Wrong input");
				continue;
			}

			printf("\nDo you see micropulse at high frequency? (1=yes,-1=no,0=o)\n");
			scanf("%d",&val);
			if (val==1)
				visarr[visualPulseindex].microathighfreq=1;
			else if (val==-1)
				visarr[visualPulseindex].microathighfreq=-1;
			else if (val==0)
				visarr[visualPulseindex].microathighfreq=0;
			else
			{
				printf("Wrong input");
				continue;
			}

			visarr[visualPulseindex].pulse=pulseNum;

			fprintf(vismicro,"%d\t%d\t%d\n",pulseNum,visarr[visualPulseindex].microathighfreq,visarr[visualPulseindex].microatlowfreq);			

			visualPulseindex++;
		}
			
		
		if(input[0]=='n')
		{
			if ((checkForFile==1)&&(input[2]!='s'))
			{
				if (input[2]=='e')
				{
					if (pulseFromFile==1)
					{
						printf("\nMoving out of file mode\n");
						pulseFromFile=0;
					}
					else
					{
						printf("\nContinuing from the last read pulse\n");
						pulseFromFile=1;
					}
	
				}
				else if ((index<(totalpulses-1))&&(pulseFromFile==1))
				{
					index++;
					pulseNum=pulseList[index];
					
				}
				else if ((pulseFromFile==1)&&(index=totalpulses-1))
				{
					printf("\nReached end of file\n");
					pulseFromFile=0;
				}
				else
					pulseNum++;
			}
			else if(input[2]=='s')
			{
				if (checkForFile==1)
				{
					printf("\nTo jump to the beginning of the same file press c.\nIf you want to change to a different file press d.\n\n");
					fgets (input,100,stdin);
					if (input[0]=='c')
					{
						printf("\nBeginning at the start of the list\n");
						pulseFromFile=1;
						index=-1;
					}
					else if (input[0]=='d')
					{	
						totalpulses=0;
						printf("\nGive the file path\n");
						scanf("%s",file);
						pulsesFile=fopen(file,"r");
						while (!feof(pulsesFile))
						{
							fgets(singleLine,100,pulsesFile);
							sscanf(singleLine, "%d", &pulseList[totalpulses]);	
							totalpulses++;
						}
						totalpulses-=2;
						printf("\nTotal number of pulses in file = %d\n",totalpulses);
						pulseFromFile=1;
						index=-1;
					}
				}
				else
				{
					totalpulses=0;
					printf("\nGive the file path\n");
					scanf("%s",file);
					pulsesFile=fopen(file,"r");
					while (!feof(pulsesFile))
					{
						fgets(singleLine,100,pulsesFile);
						sscanf(singleLine, "%d", &pulseList[totalpulses]);	
						totalpulses++;
					}
					totalpulses-=2;
					printf("\nTotal number of pulses in file = %d\n",totalpulses);
					pulseFromFile=1;
					index=-1;
					checkForFile=1;
				}
			}
			else
				pulseNum++;
			

			writefile(filename_GSB,outfile_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,readData_gsb,gsbDelay);
			writefile(filename_GWB,outfile_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,readData_gwb,gwbDelay);	

			foldprofile(filename_GSB,foldfile_GSB,fold,t_gsb,phase_gsb,pstart,pend,numshift,profile_gsb,2,intg_gsb,pulseStart,pulseEnd,gsbDelay);
			foldprofile(filename_GWB,foldfile_GWB,fold,t_gwb,phase_gwb,pstart,pend,numshift,profile_gwb,2,intg_gwb,pulseStart,pulseEnd,gwbDelay);
			
			cpgslct(pulseWindow);
			if(plotswitch==0)		
			{
				cpgsch(1.5);
				plottwo(outfile_GWB,outfile_GSB,pulseNum);
			}
			else
			{
				cpgsch(2);
				plotfour(outfile_GWB,foldfile_GWB,outfile_GSB,foldfile_GSB,pulseNum);
			}
			
		}
		else if(input[0]=='b')
		{
			if(pulseFromFile==1)
			{
				if (index>0)
				{
					index--;
					pulseNum=pulseList[index];
				}
				else
				{
					printf("\nCant go below this while in file mode, moving out of mode\n");
					pulseFromFile=0;
				}
			}
			else
				pulseNum--;
			if(pulseNum<=0)
			{
				printf("\nInvalid range! Pulse Number is %d",pulseNum);
				pulseNum=0;
				continue;
			}
                        writefile(filename_GSB,outfile_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,readData_gsb,gsbDelay);
                        writefile(filename_GWB,outfile_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,readData_gwb,gwbDelay);      
                                
			foldprofile(filename_GSB,foldfile_GSB,fold,t_gsb,phase_gsb,pstart,pend,numshift,profile_gsb,2,intg_gsb,pulseStart,pulseEnd,gsbDelay);
                        foldprofile(filename_GWB,foldfile_GWB,fold,t_gwb,phase_gwb,pstart,pend,numshift,profile_gwb,2,intg_gwb,pulseStart,pulseEnd,gwbDelay);

                        cpgslct(pulseWindow);
			if(plotswitch==0)		
			{
				cpgsch(1.5);
				plottwo(outfile_GWB,outfile_GSB,pulseNum);
			}
			else
			{
				cpgsch(2);
				plotfour(outfile_GWB,foldfile_GWB,outfile_GSB,foldfile_GSB,pulseNum);
			}
		}
		else if(input[0]=='p')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");			
			if(pch==NULL)
			{
				printf("\np_start not specified!\n");
				continue;
			}
			pstart=atof(pch);
			pch=strtok(NULL," ");						
			if(pch==NULL)
			{
                                printf("\np_end not specified!\n");
                                continue;
                        }  
			pend=atof(pch);
			
			printf("\nNew Phase start is %lf. New phase end is %lf.\n",pstart,pend);
		
                        writefile(filename_GSB,outfile_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,readData_gsb,gsbDelay);
                        writefile(filename_GWB,outfile_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,readData_gwb,gwbDelay);      
                        foldprofile(filename_GSB,foldfile_GSB,fold,t_gsb,phase_gsb,pstart,pend,numshift,profile_gsb,2,intg_gsb,pulseStart,pulseEnd,gsbDelay);
                        foldprofile(filename_GWB,foldfile_GWB,fold,t_gwb,phase_gwb,pstart,pend,numshift,profile_gwb,2,intg_gwb,pulseStart,pulseEnd,gwbDelay);

                        cpgslct(pulseWindow);
			if(plotswitch==0)		
			{
				cpgsch(1.5);
				plottwo(outfile_GWB,outfile_GSB,pulseNum);
			}
			else
			{
				cpgsch(2);
				plotfour(outfile_GWB,foldfile_GWB,outfile_GSB,foldfile_GSB,pulseNum);
			}
		}
		else if(input[0]=='q')
		{
			printf("\nTerminating program!\n");
			cpgend();
			return 0;
		}
/*		else if(input[0]=='s')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			curshift=(int)(atof(pch)*last);
			s+=atof(pch);
			
			pulseStart=pulseStart-s;
			if(pulseStart<0)
				pulseStart=1+pulseStart;
			pulseEnd=pulseEnd-s;
			if(pulseEnd<0)
				pulseEnd=1+pulseEnd;
				
			printf("\nNew shift is %lf.New pulseStart is %lf. New pulseEnd is %lf\n",s,pulseStart,pulseEnd);
			
			writefile(filename,outfile,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,readData);
			foldprofile(filename,foldfile,fold,t,phase,pstart,pend,curshift,profile,3,intg,pulseStart,pulseEnd);
			foldprofile(filename,foldfile,fold,t,phase,pstart,pend,numshift,profile,2,intg,pulseStart,pulseEnd);
			
			cpgslct(pulseWindow);
			if(plotswitch==0)			
				plotone(outfile,pulseNum,1);
			else
				plottwo(foldfile,outfile,pulseNum);
		}*/
		else if(input[0]=='f')
		{
			plotswitch=1;
writefile(filename_GSB,outfile_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,readData_gsb,gsbDelay);
                        writefile(filename_GWB,outfile_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,readData_gwb,gwbDelay);      
                                
			foldprofile(filename_GSB,foldfile_GSB,fold,t_gsb,phase_gsb,pstart,pend,numshift,profile_gsb,2,intg_gsb,pulseStart,pulseEnd,gsbDelay);
                        foldprofile(filename_GWB,foldfile_GWB,fold,t_gwb,phase_gwb,pstart,pend,numshift,profile_gwb,2,intg_gwb,pulseStart,pulseEnd,gwbDelay);

                        cpgslct(pulseWindow);
			cpgsch(2);
			plotfour(outfile_GWB,foldfile_GWB,outfile_GSB,foldfile_GSB,pulseNum);

			if (input[2] == 's')
			{

				int i,last_gsb,last_gwb,count=0;
				if (floor((double)fold/t_gsb)==fold/t_gsb)
				{
					last_gsb=floor(fold/t_gsb);
				}
				else
				{
					last_gsb=floor(fold/t_gsb)+1;
				}

				if (floor((double)fold/t_gwb)==fold/t_gwb)
				{
					last_gwb=floor(fold/t_gwb);
				}
				else
				{
					last_gwb=floor(fold/t_gwb)+1;
				}


				
				double baseline_gsb=0,baseline_gwb=0,onpulse_mean_gsb=0,onpulse_mean_gwb=0;
				double stdev_off_gsb=0,stdev_off_gwb=0,tempPhase;
				double SNR_gsb,SNR_gwb;

				baseline_gwb=meanBase(profile_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb);
				stdev_off_gwb=stdevBase(profile_gwb,phase_gwb,pulseStart,pulseEnd,last_gwb);

				for(i=0;i<last_gwb;i++)
				{
					tempPhase=phase_gwb[i];
					if(tempPhase>=pulseStart && tempPhase<=pulseEnd)
					{
						onpulse_mean_gwb+=profile_gwb[i]-baseline_gwb;
						count++;
					}			
				}
				SNR_gwb=(onpulse_mean_gwb)/((stdev_off_gwb)*sqrt(count));
				printf("\nBest estimate for SNR of folded profile GWB:%lf",SNR_gwb);

				count=0;

				baseline_gsb=meanBase(profile_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb);
				stdev_off_gsb=stdevBase(profile_gsb,phase_gsb,pulseStart,pulseEnd,last_gsb);

				for(i=0;i<last_gsb;i++)
				{
					tempPhase=phase_gsb[i];
					if(tempPhase>=pulseStart && tempPhase<=pulseEnd)
					{
						onpulse_mean_gsb+=profile_gsb[i]-baseline_gsb;
						count++;
					}			
				}
				SNR_gsb=(onpulse_mean_gsb)/((stdev_off_gsb)*sqrt(count));
				printf("\nBest estimate for SNR of folded profile GSB:%lf\n",SNR_gsb);

				double ratio;
				ratio=SNR_gsb/SNR_gwb;
				printf("SNR of GSB / SNR of GWB (folded profile) = %lf\n",ratio);

			}
				
				
		}
		else if(input[0]=='h')
		{
			plotswitch=0;
writefile(filename_GSB,outfile_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,readData_gsb,gsbDelay);
                        writefile(filename_GWB,outfile_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,readData_gwb,gwbDelay);      
                                
			foldprofile(filename_GSB,foldfile_GSB,fold,t_gsb,phase_gsb,pstart,pend,numshift,profile_gsb,2,intg_gsb,pulseStart,pulseEnd,gsbDelay);
                        foldprofile(filename_GWB,foldfile_GWB,fold,t_gwb,phase_gwb,pstart,pend,numshift,profile_gwb,2,intg_gwb,pulseStart,pulseEnd,gwbDelay);

                        cpgslct(pulseWindow);
			cpgsch(1.5);
			plottwo(outfile_GWB,outfile_GSB,pulseNum);
		}
		else if(input[0]=='i')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");			
			if(pch==NULL)
			{	
				printf("\nBackend not specified!\n");
				continue;
			}
			backOp=atoi(pch);
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nIntegration not specified!\n");
				continue;
			}
			if(backOp==1)
			{
				intg_gwb=atoi(pch);
				printf("\nNew integration is %d samples for GWB\n", intg_gwb);
			}
			else if(backOp==2)
			{
				intg_gsb=atoi(pch);
				printf("\nNew integration is %d samples for GSB\n", intg_gsb);
			}
			else
			{
				printf("\nInvalid backend option!\n");
				continue;
			}
			
			writefile(filename_GSB,outfile_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,readData_gsb,gsbDelay);
                        writefile(filename_GWB,outfile_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,readData_gwb,gwbDelay);      
                        foldprofile(filename_GSB,foldfile_GSB,fold,t_gsb,phase_gsb,pstart,pend,numshift,profile_gsb,2,intg_gsb,pulseStart,pulseEnd,gsbDelay);
                        foldprofile(filename_GWB,foldfile_GWB,fold,t_gwb,phase_gwb,pstart,pend,numshift,profile_gwb,2,intg_gwb,pulseStart,pulseEnd,gwbDelay);

                        cpgslct(pulseWindow);
			if(plotswitch==0)		
			{
				cpgsch(1.5);
				plottwo(outfile_GWB,outfile_GSB,pulseNum);
			}
			else
			{
				cpgsch(2);
				plotfour(outfile_GWB,foldfile_GWB,outfile_GSB,foldfile_GSB,pulseNum);
			}
		}
		else if(input[0]=='k')
		{
			printf("\nInstructions:\nUse 'n' to go to next pulse\nUse 'b' to go the earlier pulse\nUse 'p <start phase> <end phase>' to set a phase range (Default is full period)\nUse 's <phase shift>' to shift phase window (Phase shift should remain between 0 and 1 !)\nUse 'f' to show folded profile\nUse 'h' to hide folded profile\nUse 'i <Integration>' to integrate profile\nUse 'g <Pulse Number>' to get profile for a specified pulse number\nUse 'a <Time lag>' to get autocorrelation plot for the current pulse, integration and phase range\nUse 'c <MinSNR>' to get list of pulse numbers above a given threshold\nUse ‘m <min SNR>’ to output list of pulses above some average SNR \nUse 't <MinSNR>' to get list of pulse numbers with peak above a given SNR\nUse 'v <Time lag(ms)> to get the average ACF for the strongest pulses\nUse 'y' to view the ADP for the specified window\nUse 'k' to get instructions\nUse 'q' to terminate program\nUse 'e <time lag>' to get average cross correlation\n\nSUBTRACTION FEATURE:\n\nUse 'd <smoothing (ms)> <microCuttOff(SNR)> <SG order>’ to run Analysis of Current Pulse, that smoothes the pulse, subtracts the smoothed version, and then computes widths from the subtracted ACF smoothed with a Savitzky Golay filter and computes periods from the power spectrum of the subtracted features\n\n");			
		}
		else if(input[0]=='a')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo limiting lag specified!\n");
			else
			{
				lag=atof(pch);
				if(autoflag==0)
				{
					autoflag=1;
					autoWindow=cpgopen("/xs");
					cpgask(0);
					cpgsubp(1,2);
				}
				cpgslct(autoWindow);
				cpgsubp(1,2);
				cpgsch(1.5);
				correlate(filename_GWB,outcorr_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,lag,pulseStart,pulseEnd,readData_gwb,gwbDelay,0);
				cpgpanl(1,1);
				correlate(filename_GSB,outcorr_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,lag,pulseStart,pulseEnd,readData_gsb,gsbDelay,1);
			}
		}
		else if(input[0]=='z')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo limiting lag specified!\n");
			else
			{
				lag=atof(pch);
				if(autoflag==0)
				{
					autoflag=1;
					autoWindow=cpgopen("/xs");
					cpgask(0);
					cpgsubp(1,2);
				}
				cpgslct(autoWindow);
				cpgsubp(1,3);
				cpgsch(2.5);
				correlate(filename_GWB,outcorr_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,lag,pulseStart,pulseEnd,readData_gwb,gwbDelay,0);
				cpgpanl(1,1);
				correlate(filename_GSB,outcorr_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,lag,pulseStart,pulseEnd,readData_gsb,gsbDelay,1);
				cpgpanl(1,2);
				crossCorrelate(filename_GSB,filename_GWB,crossCorr,pulseNum,last_gsb,last_gwb,s,folda_gsb,folda_gwb,pstart,pend,pulseStart,pulseEnd,phase_gsb,phase_gwb,intg_gsb,intg_gwb,fold,t_gsb,t_gwb,lag,readData_gsb,readData_gwb,gsbDelay,gwbDelay);
				
			}
		}
		else if(input[0]=='g')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");			
			if(pch==NULL)
			{
                                printf("\npulseNum not specified!\n");
                                continue;
		        }  
			pulseNum=atoi(pch);

			writefile(filename_GSB,outfile_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,readData_gsb,gsbDelay);
                        writefile(filename_GWB,outfile_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,readData_gwb,gwbDelay);      
                        foldprofile(filename_GSB,foldfile_GSB,fold,t_gsb,phase_gsb,pstart,pend,numshift,profile_gsb,2,intg_gsb,pulseStart,pulseEnd,gsbDelay);
                        foldprofile(filename_GWB,foldfile_GWB,fold,t_gwb,phase_gwb,pstart,pend,numshift,profile_gwb,2,intg_gwb,pulseStart,pulseEnd,gwbDelay);

                        cpgslct(pulseWindow);
                        if(plotswitch==0)		
			{
				cpgsch(1.5);
				plottwo(outfile_GWB,outfile_GSB,pulseNum);
			}
			else
			{
				cpgsch(2);
				plotfour(outfile_GWB,foldfile_GWB,outfile_GSB,foldfile_GSB,pulseNum);
			}
		}
		else if(input[0]=='m')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo minimum SNR specified!\n");
			else
			{
				minSNR=atof(pch);
				printf("\nStrong pulses for GWB:\n");
				showStrongPulses(filename_GWB,last_gwb,s,folda_gwb,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,minSNR,strong_GWB,readData_gwb,gwbDelay);
				printf("\nStrong pulses for GSB:\n");
				showStrongPulses(filename_GSB,last_gsb,s,folda_gsb,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,minSNR,strong_GSB,readData_gsb,gsbDelay);				
			}
		}
		else if(input[0]=='t')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo minimum SNR specified!\n");
			else
			{
				minSNR=atof(pch);
				printf("\nStrong PEAK pulses for GWB:\n");
				showStrongPulsesPeak(filename_GWB,last_gwb,s,folda_gwb,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,minSNR,strong_GWB,readData_gwb,gwbDelay);
				printf("\nStrong PEAK pulses for GSB:\n");
				showStrongPulsesPeak(filename_GSB,last_gsb,s,folda_gsb,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,minSNR,strong_GSB,readData_gsb,gsbDelay);				
			}
		}
		else if(input[0]=='v')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo limiting lag specified!\n");
			else
			{
				lag=atof(pch);
				if(autoflag==0)
				{
					autoflag=1;
					autoWindow=cpgopen("/xs");
					cpgask(0);
					cpgsubp(1,2);
				}
				cpgslct(autoWindow);
				cpgsubp(1,2);
				cpgsch(1.5);
				averageAutoCorrelate(filename_GWB,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,lag,pulseStart,pulseEnd,strong_GWB,avgcorr_GWB,readData_gwb,gwbDelay,0);
				cpgpanl(1,1);
				averageAutoCorrelate(filename_GSB,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,lag,pulseStart,pulseEnd,strong_GSB,avgcorr_GSB,readData_gsb,gsbDelay,1);
				
				width_gsb=getAvgHalfWidth(avgcorr_GSB,t_gsb,lag,intg_gsb);
				width_gwb=getAvgHalfWidth(avgcorr_GWB,t_gwb,lag,intg_gwb);
			}
		}
		else if(input[0]=='u')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo smoothing specified!\n");
				continue;
			}			
			if(pch[0]=='d')
			{
				if(width_gsb==0 || width_gwb==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur_gsb=width_gsb;
				smoothDur_gwb=width_gwb;
				printf("\nTaking default ACF half-width GSB = %f; ACF half-width GWB = %f\n",smoothDur_gsb,smoothDur_gwb);
			}
/*			else if(pch[0]=='s')
			{
				if(averageAutocorrHalfWidths[0]==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur=averageAutocorrHalfWidths[0];
				printf("\nTaking total intensity ACF half-width = %f ms as smooth duration\n",averageAutocorrHalfWidths[0]);
			}*/
			else		
			{	smoothDur_gsb=atof(pch);
				if(smoothDur_gsb<0)
				{
					printf("\nInvalid smoothing duration!\n");
					continue;
				}	
				smoothDur_gwb=smoothDur_gsb;
			}
			
			pch=strtok(NULL," ");						
			if(pch==NULL)
			{
				printf("\nNo lag specified!\n");
				continue;
			}
			lag=atof(pch);
			
			if(lag<0)
			{
				printf("\nInvalid lag!\n");
				continue;
			}
			
			if(autoflag==0)
			{
				autoWindow=cpgopen("/xs");
				cpgask(0);
				autoflag=1;
			}
			
			cpgslct(autoWindow);
			
/*			if(smoothDur==-1)
				printf("\nPlotting pulse smoothed to ACF half-width and autocorrelation upto lag %f ms.\n",lag);
			else*/
			printf("\nPlotting pulse smoothed to GSB: %f ms GWB: %f ms, and autocorrelation upto lag %f ms.\n",smoothDur_gsb,smoothDur_gwb,lag);
			
			cpgsch(2.5);
			smoothSubACF(filename_GSB,filename_GWB,subAutoFile,phase_gsb,phase_gwb,s,folda_gsb,folda_gwb,last_gsb,last_gwb,pstart,pend,pulseStart,pulseEnd,intg_gsb,intg_gwb,pulseNum,fold,t_gsb,t_gwb,lag,smoothDur_gsb,smoothDur_gwb,gsbDelay,gwbDelay,readData_gsb,readData_gwb);
		}
		else if(input[0]=='d')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo smoothing specified!\n");
				continue;
			}			
			if(pch[0]=='d')
			{
				if(width_gsb==0 || width_gwb==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur_gsb=width_gsb;
				smoothDur_gwb=width_gwb;
				printf("\nTaking default ACF half-width GSB = %f; ACF half-width GWB = %f\n",smoothDur_gsb,smoothDur_gwb);
			}
/*			else if(pch[0]=='s')
			{
				if(averageAutocorrHalfWidths[0]==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur=averageAutocorrHalfWidths[0];
				printf("\nTaking total intensity ACF half-width = %f ms as smooth duration\n",averageAutocorrHalfWidths[0]);
			}*/
			else		
			{	smoothDur_gsb=atof(pch);
				if(smoothDur_gsb<0)
				{
					printf("\nInvalid smoothing duration!\n");
					continue;
				}	
				smoothDur_gwb=smoothDur_gsb;
			}
			
			pch=strtok(NULL," ");						
			if(pch==NULL)
			{
				printf("\nNo lag specified!\n");
				continue;
			}
			lag=atof(pch);
			
			if(lag<0)
			{
				printf("\nInvalid lag!\n");
				continue;
			}
			
			pch=strtok(NULL," ");     
                        if(pch==NULL)
                        {
                                printf("\nNo SG order specified!\n");
                                continue;
  	                }
			sg_coeff_size=atoi(pch);
			
			if(autoflag==0)
			{
				autoWindow=cpgopen("/xs");
				cpgask(0);
				autoflag=1;
			}
			
			cpgslct(autoWindow);
			
/*			if(smoothDur==-1)
				printf("\nPlotting pulse smoothed to ACF half-width and autocorrelation upto lag %f ms.\n",lag);
			else*/
			printf("\nPlotting pulse smoothed to GSB: %f ms GWB: %f ms, and autocorrelation upto lag %f ms.\n",smoothDur_gsb,smoothDur_gwb,lag);
			
			cpgsch(2.5);
			smoothSubSGACF(filename_GSB,filename_GWB,subAutoFile,SG_coeffs,phase_gsb,phase_gwb,s,folda_gsb,folda_gwb,last_gsb,last_gwb,pstart,pend,pulseStart,pulseEnd,intg_gsb,intg_gwb,pulseNum,fold,t_gsb,t_gwb,lag,smoothDur_gsb,smoothDur_gwb,gsbDelay,gwbDelay,readData_gsb,readData_gwb,sg_coeff_size,option,resicrosscorr);
		}
		else if(input[0]=='s')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo smoothing specified!\n");
				continue;
			}			
			if(pch[0]=='d')
			{
				if(width_gsb==0 || width_gwb==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur_gsb=width_gsb;
				smoothDur_gwb=width_gwb;
				printf("\nTaking default ACF half-width GSB = %f; ACF half-width GWB = %f\n",smoothDur_gsb,smoothDur_gwb);
			}
		/*	else if(pch[0]=='s')
			{
				if(averageAutocorrHalfWidths[0]==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur=averageAutocorrHalfWidths[0];
				printf("\nTaking total intensity ACF half-width = %f ms as smooth duration\n",averageAutocorrHalfWidths[0]);
			}*/
			else		
			{	smoothDur_gsb=atof(pch);
				if(smoothDur_gsb<0)
				{
					printf("\nInvalid smoothing duration!\n");
					continue;
				}	
				smoothDur_gwb=smoothDur_gsb;
			}
			
			if(autoflag==0)
			{
				autoWindow=cpgopen("/xs");
				cpgask(0);
				autoflag=1;
			}
			
			cpgslct(autoWindow);
			
/*			if(smoothDur==-1)
				printf("\nPlotting pulse smoothed to ACF half-width and autocorrelation upto lag %f ms.\n",lag);
			else*/
			printf("\nPlotting pulse smoothed to GSB: %f ms GWB: %f ms, and power spectrum of residuals.\n",smoothDur_gsb,smoothDur_gwb);
			
			cpgsch(2.5);
			smoothSubPowerSpectrum(filename_GSB,filename_GWB,FFTName,phase_gsb,phase_gwb,s,folda_gsb,folda_gwb,last_gsb,last_gwb,pstart,pend,pulseStart,pulseEnd,intg_gsb,intg_gwb,pulseNum,fold,t_gsb,t_gwb,smoothDur_gsb,smoothDur_gwb,gsbDelay,gwbDelay,readData_gsb,readData_gwb,width_gsb,width_gwb,option);
		}
		else if (input[0]=='c')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo smoothing specified!\n");
				continue;
			}			
			if(pch[0]=='d')
			{
				if(width_gsb==0 || width_gwb==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur_gsb=width_gsb;
				smoothDur_gwb=width_gwb;
				printf("\nTaking default ACF half-width GSB = %f; ACF half-width GWB = %f\n",smoothDur_gsb,smoothDur_gwb);
			}
		/*	else if(pch[0]=='s')
			{
				if(averageAutocorrHalfWidths[0]==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur=averageAutocorrHalfWidths[0];
				printf("\nTaking total intensity ACF half-width = %f ms as smooth duration\n",averageAutocorrHalfWidths[0]);
			}*/
			else		
			{	smoothDur_gsb=atof(pch);
				if(smoothDur_gsb<0)
				{
					printf("\nInvalid smoothing duration!\n");
					continue;
				}	microFileName;
				smoothDur_gwb=smoothDur_gsb;
			}
			
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo high freq cutOff specified!\n");
				continue;
			}
			cutOffhighfreq=atof(pch);
			
			if(cutOffhighfreq<0)
			{
				printf("\nInvalid cut off!\n");
				continue;
			}

			if(pch==NULL)
			{
				printf("\nNo low freq cutOff specified!\n");
				continue;
			}
			cutOfflowfreq=atof(pch);
			
			if(cutOfflowfreq<0)
			{
				printf("\nInvalid cut off!\n");
				continue;
			}
			
/*			if(smoothDur==-1)
				printf("\nSmoothing window is ACF half-width. cutOff is %f. Frequency threshold set at %d times RMS\n",cutOff,powerThreshold);		
			else*/
			printf("\nSmoothing window is GSB %f ms; GWB %f ms. highfreq cutOff is %f. lowfreq cutOff is %f. Frequency threshold set at %d times RMS\n",smoothDur_gsb,smoothDur_gwb,cutOffhighfreq,cutOfflowfreq,powerThreshold);
			
			scanPulseNum=1;
			FILE *microFile=fopen(microFileName,"w");		
			fprintf(microFile,"pulseNum\trmsGSB\trmsGWB\twidthGSB\twidthGWB\tperiodGSB\tperiodGWB\trelStGSB\trelStGWB\ttotSmoothGSB\ttotSmoothGWB\n");			
			while(microCandidates(filename_GSB,filename_GWB,microFile,phase_gsb,phase_gwb,s,folda_gsb,folda_gwb,last_gsb,last_gwb,pstart,pend,pulseStart,pulseEnd,intg_gsb,intg_gwb,scanPulseNum,fold,t_gsb,t_gwb,lag,smoothDur_gsb,smoothDur_gwb,gsbDelay,gwbDelay,readData_gsb,readData_gwb,cutOffhighfreq,cutOfflowfreq,width_gsb,width_gwb)==1)
			{
				scanPulseNum++;
			}
			
			printf("\nReached end of file,\n");	
			fclose(microFile);
		}
		else if (input[0]=='y')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo smoothing specified!\n");
				continue;
			}			
			if(pch[0]=='d')
			{
				if(width_gsb==0 || width_gwb==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur_gsb=width_gsb;
				smoothDur_gwb=width_gwb;
				printf("\nTaking default ACF half-width GSB = %f; ACF half-width GWB = %f\n",smoothDur_gsb,smoothDur_gwb);
			}
		/*	else if(pch[0]=='s')
			{
				if(averageAutocorrHalfWidths[0]==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur=averageAutocorrHalfWidths[0];
				printf("\nTaking total intensity ACF half-width = %f ms as smooth duration\n",averageAutocorrHalfWidths[0]);
			}*/
			else		
			{	smoothDur_gsb=atof(pch);
				if(smoothDur_gsb<0)
				{
					printf("\nInvalid smoothing duration!\n");
					continue;
				}	
				smoothDur_gwb=smoothDur_gsb;
			}
			
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo highfreq cutOff specified!\n");
				continue;
			}
			cutOffhighfreq=atof(pch);

			if(cutOffhighfreq<0)
			{
				printf("\nInvalid cut off!\n");
				continue;
			}

			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo lowfreq cutOff specified!\n");
				continue;
			}
			cutOfflowfreq=atof(pch);
			
			if(cutOfflowfreq<0)
			{
				printf("\nInvalid cut off!\n");
				continue;
			}
			
			pch=strtok(NULL," ");	
			if(pch==NULL)
                        {
                                printf("\nNo SG order specified!\n");
                                continue;
  	                }
			sg_coeff_size=atoi(pch);
			
/*			if(smoothDur==-1)
				printf("\nSmoothing window is ACF half-width. cutOff is %f. Frequency threshold set at %d times RMS\n",cutOff,powerThreshold);		
			else*/
			printf("\nSmoothing window is GSB %f ms; GWB %f ms. cutOffhigh is %f.cutOfflow is %f. Frequency threshold set at %d times RMS. SG order is %d\n",smoothDur_gsb,smoothDur_gwb,cutOffhighfreq,cutOfflowfreq,powerThreshold,sg_coeff_size);
			
			scanPulseNum=1;
			FILE *microFile=fopen(microFileName,"w");		
			fprintf(microFile,"pulseNum\trmsGSB\trmsGWB\tmaxDevRNSGSB\tmaxDevRMSGWB\twidthGSB\twidthGWB\tperiodGSB\tperiodGWB\trelStGSB\trelStGWB\ttotSmoothGSB\ttotSmoothGWB\tMicrostrPowerGSB\tMicrostrPowerGWB\n");			
			while(microCandidates_SG(filename_GSB,filename_GWB,microFile,SG_coeffs,phase_gsb,phase_gwb,s,folda_gsb,folda_gwb,last_gsb,last_gwb,pstart,pend,pulseStart,pulseEnd,intg_gsb,intg_gwb,scanPulseNum,fold,t_gsb,t_gwb,lag,smoothDur_gsb,smoothDur_gwb,gsbDelay,gwbDelay,readData_gsb,readData_gwb,cutOffhighfreq,cutOfflowfreq,width_gsb,width_gwb,sg_coeff_size,option)==1)
			{
				scanPulseNum++;
			}
			
			printf("\nReached end of file,\n");	
			fclose(microFile);
		}
		else if(input[0]=='r')
		{
			FILE* microFile=fopen(microFileName,"r");
			FILE* recordFile=fopen(recordFileName,"w");
			int recPulseNumber;
			float usel;
			while(fscanf(microFile,"%d %f %f %f %f %f %f %f %f %f %f %f %f",&recPulseNumber,&usel,&usel,&usel,&usel,&usel,&usel,&usel,&usel,&usel,&usel,&usel,&usel)==13)
			{
				if(recordMicro(filename_GSB,outfile_GSB,filename_GWB,outfile_GWB,recPulseNumber,last_gsb,last_gwb,s,folda_gsb,folda_gwb,pstart,pend,phase_gsb,phase_gwb,intg_gsb,intg_gwb,fold,t_gsb,t_gwb,pulseStart,pulseEnd,readData_gsb,readData_gwb,gsbDelay,gwbDelay,recordFile,pulseWindow,foldfile_GSB,foldfile_GWB,profile_gsb,profile_gwb,numshift)==-1)
					break;
			}
			fclose(recordFile);
			fclose(microFile);

		}
		/*else if(input[0]=='x')
		{
			cpgslct(pulseWindow);
			plotSpectra(filename_GSB,filename_GWB,specFile,pulseNum,last_gsb,last_gwb,s,folda_gsb,folda_gwb,pstart,pend,pulseStart,pulseEnd,phase_gsb,phase_gwb,intg_gsb,intg_gwb,fold,t_gsb,t_gwb,readData_gsb,readData_gwb,gsbDelay,gwbDelay);
		}*/
		else if(input[0]=='x')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo smoothing specified!\n");
				continue;
			}			
			if(pch[0]=='d')
			{
				if(width_gsb==0 || width_gwb==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur_gsb=width_gsb;
				smoothDur_gwb=width_gwb;
				printf("\nTaking default ACF half-width GSB = %f; ACF half-width GWB = %f\n",smoothDur_gsb,smoothDur_gwb);
			}
/*			else if(pch[0]=='s')
			{
				if(averageAutocorrHalfWidths[0]==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur=averageAutocorrHalfWidths[0];
				printf("\nTaking total intensity ACF half-width = %f ms as smooth duration\n",averageAutocorrHalfWidths[0]);
			}*/
			else		
			{	smoothDur_gsb=atof(pch);
				if(smoothDur_gsb<0)
				{
					printf("\nInvalid smoothing duration!\n");
					continue;
				}	
				smoothDur_gwb=smoothDur_gsb;
			}
			
			pch=strtok(NULL," ");						
			if(pch==NULL)
			{
				printf("\nNo lag specified!\n");
				continue;
			}
			lag=atof(pch);
			
			if(lag<0)
			{
				printf("\nInvalid lag!\n");
				continue;
			}
			
			pch=strtok(NULL," ");     
                        if(pch==NULL)
                        {
                                printf("\nNo SG order specified!\n");
                                continue;
  	                }
			sg_coeff_size=atoi(pch);
			
			if(autoflag==0)
			{
				autoWindow=cpgopen("/xs");
				cpgask(0);
				autoflag=1;
			}
			
			cpgslct(autoWindow);
			
/*			if(smoothDur==-1)
				printf("\nPlotting pulse smoothed to ACF half-width and autocorrelation upto lag %f ms.\n",lag);
			else*/
			printf("\nPlotting pulse smoothed to GSB: %f ms GWB: %f ms, and autocorrelation upto lag %f ms.\n",smoothDur_gsb,smoothDur_gwb,lag);
			
			cpgsch(2.5);
			smoothSubSGACF_intrapolated(filename_GSB,filename_GWB,subAutoFile,SG_coeffs,phase_gsb,phase_gwb,s,folda_gsb,folda_gwb,last_gsb,last_gwb,pstart,pend,pulseStart,pulseEnd,intg_gsb,intg_gwb,pulseNum,fold,t_gsb,t_gwb,lag,smoothDur_gsb,smoothDur_gwb,gsbDelay,gwbDelay,readData_gsb,readData_gwb,sg_coeff_size,option,resicrosscorr);	
		}
		else if(input[0]=='j')
		{
			createAverageSpectra(filename_GSB,filename_GWB,avgPulseSpec,foldfile_GSB,foldfile_GWB,last_gsb,last_gwb,s,folda_gsb,folda_gwb,pstart,pend,pulseStart,pulseEnd,phase_gsb,phase_gwb,intg_gsb,intg_gwb,fold,t_gsb,t_gwb,readData_gsb,readData_gwb,gsbDelay,gwbDelay);
		}
/*		else if(input[0]=='y')
		{
			if(adpflag==0)
			{
				adpflag=1;
				adpWindow=cpgopen("/xs");
				cpgask(0);
				cpgsubp(1,2);
			}
			cpgslct(adpWindow);
			ADPcreatePulse(filename_R,adp_R,pulseNum,last,s,folda_R,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,2);
			ADPcreatePulse(filename_L,adp_L,pulseNum,last,s,folda_L,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,1);
		}*/
		else if (input[0]=='w')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo option specified!\n");
				continue;
			}
			if(pch[0]=='a')
			{
				pch=strtok(input," ");
				pch=strtok(NULL," ");
				if (pch==NULL)
					printf("\nNo limiting lag specified!\n");
				else
				{
					lag=atof(pch);
					psWindow=cpgopen("autoCorr.ps/CPS");
					cpgask(0);
					cpgslct(psWindow);					
					cpgsubp(1,2);

					cpgsch(1.5);
					correlate(filename_GWB,outcorr_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,lag,pulseStart,pulseEnd,readData_gwb,gwbDelay,0);
					cpgpanl(1,1);
					correlate(filename_GSB,outcorr_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,lag,pulseStart,pulseEnd,readData_gsb,gsbDelay,1);
					cpgend();
					foldWindow=cpgopen("/xs");
					cpgask(0);
					if(autoflag==1)
					{
						autoWindow=cpgopen("/xs");
						cpgask(0);
					}
				}
			}
			else if (pch[0]=='p')
			{
				printf("\nPrinting current pulse window.\n");
				psWindow=cpgopen("pulsePol.ps/CPS");
				cpgask(0);
				cpgslct(psWindow);
				
				writefile(filename_GSB,outfile_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,readData_gsb,gsbDelay);
				writefile(filename_GWB,outfile_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,readData_gwb,gwbDelay);	

				foldprofile(filename_GSB,foldfile_GSB,fold,t_gsb,phase_gsb,pstart,pend,numshift,profile_gsb,2,intg_gsb,pulseStart,pulseEnd,gsbDelay);
				foldprofile(filename_GWB,foldfile_GWB,fold,t_gwb,phase_gwb,pstart,pend,numshift,profile_gwb,2,intg_gwb,pulseStart,pulseEnd,gwbDelay);
			
				if(plotswitch==0)		
				{
					cpgsch(1.5);
					plottwo(outfile_GWB,outfile_GSB,pulseNum);
				}
				else
				{
					cpgsch(2);
					plotfour(outfile_GWB,foldfile_GWB,outfile_GSB,foldfile_GSB,pulseNum);
				}
					
				cpgend();
				foldWindow=cpgopen("/xs");
				cpgask(0);
				if(autoflag==1)
				{
					autoWindow=cpgopen("/xs");
					cpgask(0);
				}
			}
			else if (pch[0]=='f')
			{
				printf("\nPrinting folded profile.\n");
				psWindow=cpgopen("foldPol.ps/CPS");
				cpgask(0);
				cpgslct(psWindow);
				
				cpgsch(1.5);
				plottwo(foldfile_GWB,foldfile_GSB,-1);
					
				cpgend();
				foldWindow=cpgopen("/xs");
				cpgask(0);
				if(autoflag==1)
				{
					autoWindow=cpgopen("/xs");
					cpgask(0);
				}
			}
			else if (input[0]=='z')
			{
				
			}
		}
		else if (input[0]=='e')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo limiting lag specified!\n");
			else
			{
				lag=atof(pch);
				if(autoflag==0)
				{
					autoflag=1;
					autoWindow=cpgopen("/xs");
					cpgask(0);
				}
				cpgslct(autoWindow);
				averageCrossCorrelate(filename_GSB,filename_GWB,last_gsb,last_gwb,s,folda_gsb,folda_gwb,pstart,pend,pulseStart, pulseEnd,phase_gsb,phase_gwb,intg_gsb,intg_gwb,fold,t_gsb,t_gwb,lag,readData_gsb,readData_gwb,gsbDelay,gwbDelay,strong_GSB,avg_crossCorr);
			}
		}
		else if (input[0]=='o')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo limiting lag specified!\n");
			else
			{
				lag=atof(pch);
				if(autoflag==0)
				{
					autoflag=1;
					autoWindow=cpgopen("/xs");
					cpgask(0);
				}
			}
			cpgslct(autoWindow);
			foldedCrossCorr(foldfile_GSB,foldfile_GWB,pstart,pend,lag,t_gsb,t_gwb,intg_gsb,intg_gwb,fold,last_gsb,last_gwb);		

		}
		else
		{
			printf("\nInvalid option. Try again.\n");
		}
        
	}
	fclose (pulsesFile);
	fclose (vismicro);
	cpgend();
	return EXIT_SUCCESS;
	}
	case 1:
	{
	double t,fold,value,time,k,lag=0,DM,smoothDur=0;
	int last,i,sample=0,count=0,scanPulseNum=0;
	int b=0;
	char *input=(char*)malloc(100*sizeof(char));

        int in;
        char datafile[1000],file[1000];
	float reso,pulseper,dm;
	FILE *parameter = fopen("/Data/ygupta/psr_microstr_work/recovered_code/dual_freq_analysis/Parameters3.garvit.txt", "r");
	
	if (parameter == NULL) 
        {
        perror("Error: Failed to open file.");
        return 1;
        }

	char *parameters[50],line[300];
	for(i=0;i<50;i++) parameters[i]=(char*) malloc(sizeof(char)*300);
	int num_pulsars=0;	

	for(i=0;fgets(line,sizeof(line),parameter);i++)	
	{
		strcpy(parameters[i],line);
		printf("%s\n",parameters[i]);
	}
	num_pulsars=i;
	fclose(parameter);

        int pulsarnum;
	printf("Enter pulsar number ");
	scanf("%d", &pulsarnum);

	if (pulsarnum > num_pulsars || pulsarnum < 1)
            {
              printf ("Wrong input\n");
              return 1;
            }

	char pulsarchosen[100];
	int num;
        sscanf(parameters[pulsarnum-1], "%d)%s %f %f %f %s",&num,pulsarchosen,&reso,&pulseper,&dm,datafile);

	printf("\nPulsar is %s",pulsarchosen);        
        printf("\nResolution %f",reso);
        printf("\nPulse period is %f",pulseper);
        printf("\nDM is %f",dm);

	char singleLine[100];
	int index=-1;
	int pulseList[100000];
	int totalpulses=0;
        FILE *pulsesFile;
	int pulseFromFile=0,checkForFile=0;
	printf ("\nIf you want to read pulses from a file press 1,else press 0\n");
	scanf ("%d",&checkForFile);
	if (checkForFile==1)
           {
		pulseFromFile=1;
		printf("\nGive file path\n");
		scanf("%s",file);
		//printf("%s",file);
		pulsesFile=fopen(file,"r");
		while (!feof(pulsesFile))
		{
			fgets(singleLine,100,pulsesFile);
			sscanf(singleLine, "%d", &pulseList[totalpulses]);	
			totalpulses++;
		}
		totalpulses-=1;
		printf("\nTotal number of pulses in file = %d\n",totalpulses);
		//int i;
		//for (i=0;i<totalpulses;i++)
			//printf("\n%d",pulseList[i]);
	    }
         else if (checkForFile!=0) 
         {
		printf("\nWrong option chosen.\n");
		return 1;
	 }
	
	t=reso;
	fold=pulseper;
	DM=dm;

	char *filename=datafile;
	char *outfile="plot.txt";
	char *foldfile="fold.txt";
	char *outcorr="corr.txt";
	char *avgcorr="avgcorr.txt";
	char *strong="strongpulse.txt";
	char *adp="ADP.txt";
	char *avg_crossCorr="avgCrossCorr.txt";
	char *crossCorr="crossCorr.txt";
	char *subAutoFile="subAuto.txt";
	char *FFTName="SpectrumFile.txt";
	char *microFileName="microStat.txt";
	char *recordFileName="recordStat.txt";
	char *specFile="specPlot.txt";
	char *avgPulseSpec="avgSpec.txt";
	char *SG_coeffs="SG_coeffs.txt";
	char *visualMicroFile="Vis_micropulse.txt";

	/*double t_gsb,t_gwb,fold,value,time,k,lag=0,DM,smoothDur_gsb=0,smoothDur_gwb=0;
	int last_gsb,last_gwb,i,sample=0,count=0,scanPulseNum=0;
	int b=0;
	char *input=(char*)malloc(100*sizeof(char));*/

        /*t_gsb=(double)atof(argv[3]);
	t_gwb=(double)atof(argv[4]);
	fold=(double)atof(argv[5]);
	DM=(double)atof(argv[6]);*/
	
	if (floor((double)fold/t)==fold/t)
		last=floor(fold/t);
	else
		last=floor(fold/t)+1;
	
	
	float cutOff=0;
        
	double *folda=(double*)malloc(last*sizeof(double));
	float *readData=(float*)malloc(last*sizeof(float));
	double *phase=(double*)malloc(last*sizeof(double));
	double *profile=(double*)malloc(last*sizeof(double));	
	
	//printf("\nPulse Viewer\nInstructions:\nUse 'n' to go to next pulse\nUse 'b' to go the earlier pulse\nUse 'p <start phase> <end phase>' to set a phase range (Default is full period)\nUse 's <phase shift>' to shift phase window (Phase shift should remain between 0 and 1 !)\nUse 'f' to show folded profile\nUse 'h' to hide folded profile\nUse 'i <Integration>' to integrate profile\nUse 'g <Pulse Number>' to get profile for a specified pulse number\nUse 'a <Time lag>' to get autocorrelation plot for the current pulse, integration and phase range\nUse 'c <MinSNR>' to get list of pulse numbers above a given threshold\nUse 't <MinSNR>' to get list of pulse numbers with peak above a given SNR\nUse 'v <Time lag(ms)> to get the average ACF for the strongest pulses\nUse 'y' to view the ADP for the specified window\nUse 'k' to get instructions\nUse 'q' to terminate program\n");
	int pulseNum=0,numshift=0,plotswitch=0,intg=1,curshift=0,backOp,visualPulseindex=0,val=0;
	double s=0,pstart=0,pend=1,pulseStart=0,pulseEnd=1,minSNR=0;
	char *pch;
	struct  visualPulse visarr[1000];
	int autoflag=0, autoWindow, psWindow, foldWindow;	
	int sg_coeff_size=0;
	int pulseWindow=cpgopen("/xs");
	cpgask(0);
	
	double width;
	
	for(i=0;i<last;i++)
	{
		phase[i]=(double)i/last;
	}

	foldprofile(filename,foldfile,fold,t,phase,pstart,pend,numshift,profile,1,intg,pulseStart,pulseEnd,0);

	printf("\nDefault integration is %d samples.\n",intg);
	
	printf("\nPlease check visually and specify pulse on window for baseline and noise statistics. Showing the two folded profiles.\n");

	cpgslct(pulseWindow);
	cpgsch(1.5);
	plotone_single(foldfile,-1,2);
	
	printf("\nEnter pulse start:");
	scanf("%lf",&pulseStart);
	printf("\nEnter pulse end:");
	scanf("%lf",&pulseEnd);
	

	FILE* vismicro=fopen(visualMicroFile,"w");

	//to get past the \n character from stdin taken before
	int c;
	while ((c = getchar()) != '\n' && c != EOF) { }
		
	while(1)
	{
		printf("\nEnter option:\t\tUse 'k' to get instructions and 'q' to terminate program\n");
		fgets(input,100,stdin);

		if (input[0]=='l')
		{

			if (visualPulseindex==0)
			{
				fprintf(vismicro,"pulseNum\thighfreq\tlowfreq\n");			
			}

			printf("\nDo you see micropulse at low frequency? (1=yes,-1=no,0=o)\n");
			scanf("%d",&val);
			if (val==1)
				visarr[visualPulseindex].microatlowfreq=1;
			else if (val==-1)
				visarr[visualPulseindex].microatlowfreq=-1;
			else if (val==0)
				visarr[visualPulseindex].microatlowfreq=0;
			else
			{
				printf("Wrong input");
				continue;
			}

			printf("\nDo you see micropulse at high frequency? (1=yes,-1=no,0=o)\n");
			scanf("%d",&val);
			if (val==1)
				visarr[visualPulseindex].microathighfreq=1;
			else if (val==-1)
				visarr[visualPulseindex].microathighfreq=-1;
			else if (val==0)
				visarr[visualPulseindex].microathighfreq=0;
			else
			{
				printf("Wrong input");
				continue;
			}

			visarr[visualPulseindex].pulse=pulseNum;

			fprintf(vismicro,"%d\t%d\t%d\n",pulseNum,visarr[visualPulseindex].microathighfreq,visarr[visualPulseindex].microatlowfreq);			

			visualPulseindex++;
		}
			
		
		if(input[0]=='n')
		{
			if ((checkForFile==1)&&(input[2]!='s'))
			{
				if (input[2]=='e')
				{
					if (pulseFromFile==1)
					{
						printf("\nMoving out of file mode\n");
						pulseFromFile=0;
					}
					else
					{
						printf("\nContinuing from the last read pulse\n");
						pulseFromFile=1;
					}
	
				}
				else if ((index<(totalpulses-1))&&(pulseFromFile==1))
				{
					index++;
					pulseNum=pulseList[index];
					
				}
				else if ((pulseFromFile==1)&&(index=totalpulses-1))
				{
					printf("\nReached end of file\n");
					pulseFromFile=0;
				}
				else
					pulseNum++;
			}
			else if(input[2]=='s')
			{
				if (checkForFile==1)
				{
					printf("\nTo jump to the beginning of the same file press c.\nIf you want to change to a different file press d.\n\n");
					fgets (input,100,stdin);
					if (input[0]=='c')
					{
						printf("\nBeginning at the start of the list\n");
						pulseFromFile=1;
						index=-1;
					}
					else if (input[0]=='d')
					{	
						totalpulses=0;
						printf("\nGive the file path\n");
						scanf("%s",file);
						pulsesFile=fopen(file,"r");
						while (!feof(pulsesFile))
						{
							fgets(singleLine,100,pulsesFile);
							sscanf(singleLine, "%d", &pulseList[totalpulses]);	
							totalpulses++;
						}
						totalpulses-=2;
						printf("\nTotal number of pulses in file = %d\n",totalpulses);
						pulseFromFile=1;
						index=-1;
					}
				}
				else
				{
					totalpulses=0;
					printf("\nGive the file path\n");
					scanf("%s",file);
					pulsesFile=fopen(file,"r");
					while (!feof(pulsesFile))
					{
						fgets(singleLine,100,pulsesFile);
						sscanf(singleLine, "%d", &pulseList[totalpulses]);	
						totalpulses++;
					}
					totalpulses-=2;
					printf("\nTotal number of pulses in file = %d\n",totalpulses);
					pulseFromFile=1;
					index=-1;
					checkForFile=1;
				}
			}
			else
				pulseNum++;
			

			writefile(filename,outfile,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,readData,0);
			
			foldprofile(filename,foldfile,fold,t,phase,pstart,pend,numshift,profile,2,intg,pulseStart,pulseEnd,0);
			
			cpgslct(pulseWindow);
			if(plotswitch==0)		
			{
				cpgsch(1.5);
				plotone_single(outfile,pulseNum,1);
			}
			else
			{
				cpgsch(2);
				plottwo_single(foldfile,outfile,pulseNum);
			}
			
		}
		else if(input[0]=='b')
		{
			if(pulseFromFile==1)
			{
				if (index>0)
				{
					index--;
					pulseNum=pulseList[index];
				}
				else
				{
					printf("\nCant go below this while in file mode, moving out of mode\n");
					pulseFromFile=0;
				}
			}
			else
				pulseNum--;
			if(pulseNum<=0)
			{
				printf("\nInvalid range! Pulse Number is %d",pulseNum);
				pulseNum=0;
				continue;
			}
                        writefile(filename,outfile,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,readData,0);
                        
			foldprofile(filename,foldfile,fold,t,phase,pstart,pend,numshift,profile,2,intg,pulseStart,pulseEnd,0);
                       
			cpgslct(pulseWindow);
			if(plotswitch==0)		
			{
				cpgsch(1.5);
				plotone_single(outfile,pulseNum,1);
			}
			else
			{
				cpgsch(2);
				plottwo_single(foldfile,outfile,pulseNum);
			}
		}
		else if(input[0]=='p')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");			
			if(pch==NULL)
			{
				printf("\np_start not specified!\n");
				continue;
			}
			pstart=atof(pch);
			pch=strtok(NULL," ");						
			if(pch==NULL)
			{
                                printf("\np_end not specified!\n");
                                continue;
                        }  
			pend=atof(pch);
			
			printf("\nNew Phase start is %lf. New phase end is %lf.\n",pstart,pend);
		
                        writefile(filename,outfile,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,readData,0);
                        
			foldprofile(filename,foldfile,fold,t,phase,pstart,pend,numshift,profile,2,intg,pulseStart,pulseEnd,0);

                        cpgslct(pulseWindow);
			if(plotswitch==0)		
			{
				cpgsch(1.5);
				plotone_single(outfile,pulseNum,1);
			}
			else
			{
				cpgsch(2);
				plottwo_single(foldfile,outfile,pulseNum);
			}
		}
		else if(input[0]=='q')
		{
			printf("\nTerminating program!\n");
			cpgend();
			return 0;
		}
/*		else if(input[0]=='s')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			curshift=(int)(atof(pch)*last);
			s+=atof(pch);
			
			pulseStart=pulseStart-s;
			if(pulseStart<0)
				pulseStart=1+pulseStart;
			pulseEnd=pulseEnd-s;
			if(pulseEnd<0)
				pulseEnd=1+pulseEnd;
				
			printf("\nNew shift is %lf.New pulseStart is %lf. New pulseEnd is %lf\n",s,pulseStart,pulseEnd);
			
			writefile(filename,outfile,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,readData);
			foldprofile(filename,foldfile,fold,t,phase,pstart,pend,curshift,profile,3,intg,pulseStart,pulseEnd);
			foldprofile(filename,foldfile,fold,t,phase,pstart,pend,numshift,profile,2,intg,pulseStart,pulseEnd);
			
			cpgslct(pulseWindow);
			if(plotswitch==0)			
				plotone(outfile,pulseNum,1);
			else
				plottwo(foldfile,outfile,pulseNum);
		}*/
		else if(input[0]=='f')
		{
			plotswitch=1;
			writefile(filename,outfile,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,readData,0);
                            
			foldprofile(filename,foldfile,fold,t,phase,pstart,pend,numshift,profile,2,intg,pulseStart,pulseEnd,0);

                        cpgslct(pulseWindow);
			cpgsch(2);
			plottwo_single(foldfile,outfile,pulseNum);
		}
		else if(input[0]=='h')
		{
			plotswitch=0;
			writefile(filename,outfile,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,readData,0);
                               
			foldprofile(filename,foldfile,fold,t,phase,pstart,pend,numshift,profile,2,intg,pulseStart,pulseEnd,0);
                        
                        cpgslct(pulseWindow);
			cpgsch(1.5);
			plotone_single(outfile,pulseNum,1);
		}
		else if(input[0]=='i')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");			
			
			if(pch==NULL)
			{
				printf("\nIntegration not specified!\n");
				continue;
			}

			intg=atoi(pch);
			printf("\nNew integration is %d samples for GWB\n", intg);
			
			writefile(filename,outfile,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,readData,0);
                         
                        foldprofile(filename,foldfile,fold,t,phase,pstart,pend,numshift,profile,2,intg,pulseStart,pulseEnd,0);
                        
                        cpgslct(pulseWindow);
			if(plotswitch==0)		
			{
				cpgsch(1.5);
				plotone_single(outfile,pulseNum,1);
			}
			else
			{
				cpgsch(2);
				plottwo_single(foldfile,outfile,pulseNum);
			}
		}
		else if(input[0]=='k')
		{
			printf("\nInstructions:\nUse 'n' to go to next pulse\nUse 'b' to go the earlier pulse\nUse 'p <start phase> <end phase>' to set a phase range (Default is full period)\nUse 's <phase shift>' to shift phase window (Phase shift should remain between 0 and 1 !)\nUse 'f' to show folded profile\nUse 'h' to hide folded profile\nUse 'i <Integration>' to integrate profile\nUse 'g <Pulse Number>' to get profile for a specified pulse number\nUse ‘m <min SNR>’ to output list of pulses above some average SNR \nUse 'a <Time lag>' to get autocorrelation plot for the current pulse, integration and phase range\nUse 'c <MinSNR>' to get list of pulse numbers above a given threshold\nUse 't <MinSNR>' to get list of pulse numbers with peak above a given SNR\nUse 'v <Time lag(ms)> to get the average ACF for the strongest pulses\nUse 'y' to view the ADP for the specified window\nUse 'k' to get instructions\nUse 'q' to terminate program\nUse 'e <time lag>' to get average cross correlation\n\nSUBTRACTION FEATURE:\n\nUse 'd <smoothing (ms)> <microCuttOff(SNR)> <SG order>’ to run Analysis of Current Pulse, that smoothes the pulse, subtracts the smoothed version, and then computes widths from the subtracted ACF smoothed with a Savitzky Golay filter and computes periods from the power spectrum of the subtracted features\n\n");			
		}
		else if(input[0]=='a')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo limiting lag specified!\n");
			else
			{
				lag=atof(pch);
				if(autoflag==0)
				{
					autoflag=1;
					autoWindow=cpgopen("/xs");
					cpgask(0);
					cpgsubp(1,2);
				}
				cpgslct(autoWindow);
				cpgsubp(1,2);
				cpgsch(1.5);
				correlate(filename,outcorr,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,lag,pulseStart,pulseEnd,readData,0,0);
				cpgpanl(1,1);
				correlate(filename,outcorr,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,lag,pulseStart,pulseEnd,readData,0,0);
			}
		}
		/*else if(input[0]=='z')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo limiting lag specified!\n");
			else
			{
				lag=atof(pch);
				if(autoflag==0)
				{
					autoflag=1;
					autoWindow=cpgopen("/xs");
					cpgask(0);
					cpgsubp(1,2);
				}
				cpgslct(autoWindow);
				cpgsubp(1,3);
				cpgsch(2.5);
				correlate(filename_GWB,outcorr_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,lag,pulseStart,pulseEnd,readData_gwb,gwbDelay,0);
				cpgpanl(1,1);
				correlate(filename_GSB,outcorr_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,lag,pulseStart,pulseEnd,readData_gsb,gsbDelay,1);
				cpgpanl(1,2);
				crossCorrelate(filename_GSB,filename_GWB,crossCorr,pulseNum,last_gsb,last_gwb,s,folda_gsb,folda_gwb,pstart,pend,pulseStart,pulseEnd,phase_gsb,phase_gwb,intg_gsb,intg_gwb,fold,t_gsb,t_gwb,lag,readData_gsb,readData_gwb,gsbDelay,gwbDelay);
				
			}
		}*/
		else if(input[0]=='g')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");			
			if(pch==NULL)
			{
                                printf("\npulseNum not specified!\n");
                                continue;
		        }  
			pulseNum=atoi(pch);

			writefile(filename,outfile,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,readData,0);
                        
                        foldprofile(filename,foldfile,fold,t,phase,pstart,pend,numshift,profile,2,intg,pulseStart,pulseEnd,0);

                        cpgslct(pulseWindow);
                        if(plotswitch==0)		
			{
				cpgsch(1.5);
				plotone_single(outfile,pulseNum,1);
			}
			else
			{
				cpgsch(2);
				plottwo_single(foldfile,outfile,pulseNum);
			}
		}
		else if(input[0]=='m')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo minimum SNR specified!\n");
			else
			{
				minSNR=atof(pch);
				printf("\nStrong pulses are:\n");
				showStrongPulses(filename,last,s,folda,phase,intg,fold,t,pulseStart,pulseEnd,minSNR,strong,readData,0);		
			}
		}
		else if(input[0]=='t')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo minimum SNR specified!\n");
			else
			{
				minSNR=atof(pch);
				printf("\nStrong PEAK pulses are:\n");
				showStrongPulsesPeak(filename,last,s,folda,phase,intg,fold,t,pulseStart,pulseEnd,minSNR,strong,readData,0);				
			}
		}
		else if(input[0]=='v')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo limiting lag specified!\n");
			else
			{
				lag=atof(pch);
				if(autoflag==0)
				{
					autoflag=1;
					autoWindow=cpgopen("/xs");
					cpgask(0);
					cpgsubp(1,2);
				}
				cpgslct(autoWindow);
				cpgsubp(1,1);
				cpgsch(1.5);
				averageAutoCorrelate(filename,last,s,folda,pstart,pend,phase,intg,fold,t,lag,pulseStart,pulseEnd,strong,avgcorr,readData,0,0);
				
				width=getAvgHalfWidth(avgcorr,t,lag,intg);
			}
		}
		/*else if(input[0]=='u')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo smoothing specified!\n");
				continue;
			}			
			if(pch[0]=='d')
			{
				if(width_gsb==0 || width_gwb==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur_gsb=width_gsb;
				smoothDur_gwb=width_gwb;
				printf("\nTaking default ACF half-width GSB = %f; ACF half-width GWB = %f\n",smoothDur_gsb,smoothDur_gwb);
			}

			else		
			{	smoothDur_gsb=atof(pch);
				if(smoothDur_gsb<0)
				{
					printf("\nInvalid smoothing duration!\n");
					continue;
				}	
				smoothDur_gwb=smoothDur_gsb;
			}
			
			pch=strtok(NULL," ");						
			if(pch==NULL)
			{
				printf("\nNo lag specified!\n");
				continue;
			}
			lag=atof(pch);
			
			if(lag<0)
			{
				printf("\nInvalid lag!\n");
				continue;
			}
			
			if(autoflag==0)
			{
				autoWindow=cpgopen("/xs");
				cpgask(0);
				autoflag=1;
			}
			
			cpgslct(autoWindow);
			

			printf("\nPlotting pulse smoothed to GSB: %f ms GWB: %f ms, and autocorrelation upto lag %f ms.\n",smoothDur_gsb,smoothDur_gwb,lag);
			
			cpgsch(2.5);
			smoothSubACF(filename_GSB,filename_GWB,subAutoFile,phase_gsb,phase_gwb,s,folda_gsb,folda_gwb,last_gsb,last_gwb,pstart,pend,pulseStart,pulseEnd,intg_gsb,intg_gwb,pulseNum,fold,t_gsb,t_gwb,lag,smoothDur_gsb,smoothDur_gwb,gsbDelay,gwbDelay,readData_gsb,readData_gwb);
		}*/
		else if(input[0]=='d')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo smoothing specified!\n");
				continue;
			}			
			if(pch[0]=='d')
			{
				if(width==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur=width;
				printf("\nTaking default ACF half-width = %f\n",smoothDur);
			}
			else		
			{	smoothDur=atof(pch);
				if(smoothDur<0)
				{
					printf("\nInvalid smoothing duration!\n");
					continue;
				}	
			}
			
			pch=strtok(NULL," ");						
			if(pch==NULL)
			{
				printf("\nNo lag specified!\n");
				continue;
			}
			lag=atof(pch);
			
			if(lag<0)
			{
				printf("\nInvalid lag!\n");
				continue;
			}
			
			pch=strtok(NULL," ");     
                        if(pch==NULL)
                        {
                                printf("\nNo SG order specified!\n");
                                continue;
  	                }
			sg_coeff_size=atoi(pch);
			
			if(autoflag==0)
			{
				autoWindow=cpgopen("/xs");
				cpgask(0);
				autoflag=1;
			}
			
			cpgslct(autoWindow);
			
			printf("\nPlotting pulse smoothed to : %f ms and autocorrelation upto lag %f ms.\n",smoothDur,lag);
			
			cpgsch(2.5);
			smoothSubSGACF(filename,filename,subAutoFile,SG_coeffs,phase,phase,s,folda,folda,last,last,pstart,pend,pulseStart,pulseEnd,intg,intg,pulseNum,fold,t,t,lag,smoothDur,smoothDur,0,0,readData,readData,sg_coeff_size,option,NULL);
		}
		else if(input[0]=='s')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo smoothing specified!\n");
				continue;
			}			
			if(pch[0]=='d')
			{
				if(width==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur=width;
				printf("\nTaking default ACF half-width = %f\n",smoothDur);
			}
		/*	else if(pch[0]=='s')
			{
				if(averageAutocorrHalfWidths[0]==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur=averageAutocorrHalfWidths[0];
				printf("\nTaking total intensity ACF half-width = %f ms as smooth duration\n",averageAutocorrHalfWidths[0]);
			}*/
			else		
			{	smoothDur=atof(pch);
				if(smoothDur<0)
				{
					printf("\nInvalid smoothing duration!\n");
					continue;
				}	
			}
			
			if(autoflag==0)
			{
				autoWindow=cpgopen("/xs");
				cpgask(0);
				autoflag=1;
			}
			
			cpgslct(autoWindow);
			
/*			if(smoothDur==-1)
				printf("\nPlotting pulse smoothed to ACF half-width and autocorrelation upto lag %f ms.\n",lag);
			else*/
			printf("\nPlotting pulse smoothed : %f ms  and power spectrum of residual.\n",smoothDur);
			
			cpgsch(2.5);
			smoothSubPowerSpectrum(filename,filename,FFTName,phase,phase,s,folda,folda,last,last,pstart,pend,pulseStart,pulseEnd,intg,intg,pulseNum,fold,t,t,smoothDur,smoothDur,0,0,readData,readData,width,width,option);
		}
		/*else if (input[0]=='c')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo smoothing specified!\n");
				continue;
			}			
			if(pch[0]=='d')
			{
				if(width_gsb==0 || width_gwb==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur_gsb=width_gsb;
				smoothDur_gwb=width_gwb;
				printf("\nTaking default ACF half-width GSB = %f; ACF half-width GWB = %f\n",smoothDur_gsb,smoothDur_gwb);
			}

			else		
			{	smoothDur_gsb=atof(pch);
				if(smoothDur_gsb<0)
				{
					printf("\nInvalid smoothing duration!\n");
					continue;
				}	microFileName;
				smoothDur_gwb=smoothDur_gsb;
			}
			
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo high freq cutOff specified!\n");
				continue;
			}
			cutOffhighfreq=atof(pch);
			
			if(cutOffhighfreq<0)
			{
				printf("\nInvalid cut off!\n");
				continue;
			}

			if(pch==NULL)
			{
				printf("\nNo low freq cutOff specified!\n");
				continue;
			}
			cutOfflowfreq=atof(pch);
			
			if(cutOfflowfreq<0)
			{
				printf("\nInvalid cut off!\n");
				continue;
			}
			

			printf("\nSmoothing window is GSB %f ms; GWB %f ms. highfreq cutOff is %f. lowfreq cutOff is %f. Frequency threshold set at %d times RMS\n",smoothDur_gsb,smoothDur_gwb,cutOffhighfreq,cutOfflowfreq,powerThreshold);
			
			scanPulseNum=1;
			FILE *microFile=fopen(microFileName,"w");		
			fprintf(microFile,"pulseNum\trmsGSB\trmsGWB\twidthGSB\twidthGWB\tperiodGSB\tperiodGWB\trelStGSB\trelStGWB\ttotSmoothGSB\ttotSmoothGWB\n");			
			while(microCandidates(filename_GSB,filename_GWB,microFile,phase_gsb,phase_gwb,s,folda_gsb,folda_gwb,last_gsb,last_gwb,pstart,pend,pulseStart,pulseEnd,intg_gsb,intg_gwb,scanPulseNum,fold,t_gsb,t_gwb,lag,smoothDur_gsb,smoothDur_gwb,gsbDelay,gwbDelay,readData_gsb,readData_gwb,cutOffhighfreq,cutOfflowfreq,width_gsb,width_gwb)==1)
			{
				scanPulseNum++;
			}
			
			printf("\nReached end of file,\n");	
			fclose(microFile);
		}*/
		else if (input[0]=='y')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo smoothing specified!\n");
				continue;
			}			
			if(pch[0]=='d')
			{
				if(width==0)
				{
					printf("\nERROR: AvgWidth is zero!\n");
					continue;
				}
				smoothDur=width;
				printf("\nTaking default ACF half-width = %f\n",smoothDur);
			}

			else		
			{	smoothDur=atof(pch);
				if(smoothDur<0)
				{
					printf("\nInvalid smoothing duration!\n");
					continue;
				}	
			}
			
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo highfreq cutOff specified!\n");
				continue;
			}
			cutOff=atof(pch);

			if(cutOff<0)
			{
				printf("\nInvalid cut off!\n");
				continue;
			}
			
			pch=strtok(NULL," ");	
			if(pch==NULL)
                        {
                                printf("\nNo SG order specified!\n");
                                continue;
  	                }
			sg_coeff_size=atoi(pch);
			

			printf("\nSmoothing window is %f ms. CutOff is %f. Frequency threshold set at %d times RMS. SG order is %d\n",smoothDur,cutOff,powerThreshold,sg_coeff_size);
			
			scanPulseNum=1;
			FILE *microFile=fopen(microFileName,"w");		
			fprintf(microFile,"pulseNum\trms\tmaxDevRMS\twidth\tperiod\trelSt\ttotSmooth\tMicrostrPower\n");			
			while(microCandidates_SG(filename,filename,microFile,SG_coeffs,phase,phase,s,folda,folda,last,last,pstart,pend,pulseStart,pulseEnd,intg,intg,scanPulseNum,fold,t,t,lag,smoothDur,smoothDur,0,0,readData,readData,cutOff,cutOff,width,width,sg_coeff_size,option)==1)
			{
				scanPulseNum++;
			}
			
			printf("\nReached end of file,\n");	
			fclose(microFile);
		}
		/*else if(input[0]=='r')
		{
			FILE* microFile=fopen(microFileName,"r");
			FILE* recordFile=fopen(recordFileName,"w");
			int recPulseNumber;
			float usel;
			while(fscanf(microFile,"%d %f %f %f %f %f %f %f %f %f %f %f %f",&recPulseNumber,&usel,&usel,&usel,&usel,&usel,&usel,&usel,&usel,&usel,&usel,&usel,&usel)==13)
			{
				if(recordMicro(filename_GSB,outfile_GSB,filename_GWB,outfile_GWB,recPulseNumber,last_gsb,last_gwb,s,folda_gsb,folda_gwb,pstart,pend,phase_gsb,phase_gwb,intg_gsb,intg_gwb,fold,t_gsb,t_gwb,pulseStart,pulseEnd,readData_gsb,readData_gwb,gsbDelay,gwbDelay,recordFile,pulseWindow,foldfile_GSB,foldfile_GWB,profile_gsb,profile_gwb,numshift)==-1)
					break;
			}
			fclose(recordFile);
			fclose(microFile);

		}
		else if(input[0]=='x')
		{
			cpgslct(pulseWindow);
			plotSpectra(filename_GSB,filename_GWB,specFile,pulseNum,last_gsb,last_gwb,s,folda_gsb,folda_gwb,pstart,pend,pulseStart,pulseEnd,phase_gsb,phase_gwb,intg_gsb,intg_gwb,fold,t_gsb,t_gwb,readData_gsb,readData_gwb,gsbDelay,gwbDelay);
		}
		else if(input[0]=='j')
		{
			createAverageSpectra(filename_GSB,filename_GWB,avgPulseSpec,foldfile_GSB,foldfile_GWB,last_gsb,last_gwb,s,folda_gsb,folda_gwb,pstart,pend,pulseStart,pulseEnd,phase_gsb,phase_gwb,intg_gsb,intg_gwb,fold,t_gsb,t_gwb,readData_gsb,readData_gwb,gsbDelay,gwbDelay);
		}

		else if (input[0]=='w')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if(pch==NULL)
			{
				printf("\nNo option specified!\n");
				continue;
			}
			if(pch[0]=='a')
			{
				pch=strtok(input," ");
				pch=strtok(NULL," ");
				if (pch==NULL)
					printf("\nNo limiting lag specified!\n");
				else
				{
					lag=atof(pch);
					psWindow=cpgopen("autoCorr.ps/CPS");
					cpgask(0);
					cpgslct(psWindow);					
					cpgsubp(1,2);

					cpgsch(1.5);
					correlate(filename_GWB,outcorr_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,lag,pulseStart,pulseEnd,readData_gwb,gwbDelay,0);
					cpgpanl(1,1);
					correlate(filename_GSB,outcorr_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,lag,pulseStart,pulseEnd,readData_gsb,gsbDelay,1);
					cpgend();
					foldWindow=cpgopen("/xs");
					cpgask(0);
					if(autoflag==1)
					{
						autoWindow=cpgopen("/xs");
						cpgask(0);
					}
				}
			}
			else if (pch[0]=='p')
			{
				printf("\nPrinting current pulse window.\n");
				psWindow=cpgopen("pulsePol.ps/CPS");
				cpgask(0);
				cpgslct(psWindow);
				
				writefile(filename_GSB,outfile_GSB,pulseNum,last_gsb,s,folda_gsb,pstart,pend,phase_gsb,intg_gsb,fold,t_gsb,pulseStart,pulseEnd,readData_gsb,gsbDelay);
				writefile(filename_GWB,outfile_GWB,pulseNum,last_gwb,s,folda_gwb,pstart,pend,phase_gwb,intg_gwb,fold,t_gwb,pulseStart,pulseEnd,readData_gwb,gwbDelay);	

				foldprofile(filename_GSB,foldfile_GSB,fold,t_gsb,phase_gsb,pstart,pend,numshift,profile_gsb,2,intg_gsb,pulseStart,pulseEnd,gsbDelay);
				foldprofile(filename_GWB,foldfile_GWB,fold,t_gwb,phase_gwb,pstart,pend,numshift,profile_gwb,2,intg_gwb,pulseStart,pulseEnd,gwbDelay);
			
				if(plotswitch==0)		
				{
					cpgsch(1.5);
					plottwo(outfile_GWB,outfile_GSB,pulseNum);
				}
				else
				{
					cpgsch(2);
					plotfour(outfile_GWB,foldfile_GWB,outfile_GSB,foldfile_GSB,pulseNum);
				}
					
				cpgend();
				foldWindow=cpgopen("/xs");
				cpgask(0);
				if(autoflag==1)
				{
					autoWindow=cpgopen("/xs");
					cpgask(0);
				}
			}
			else if (pch[0]=='f')
			{
				printf("\nPrinting folded profile.\n");
				psWindow=cpgopen("foldPol.ps/CPS");
				cpgask(0);
				cpgslct(psWindow);
				
				cpgsch(1.5);
				plottwo(foldfile_GWB,foldfile_GSB,-1);
					
				cpgend();
				foldWindow=cpgopen("/xs");
				cpgask(0);
				if(autoflag==1)
				{
					autoWindow=cpgopen("/xs");
					cpgask(0);
				}
			}
			else if (input[0]=='z')
			{
				
			}
		}
		else if (input[0]=='e')
		{
			pch=strtok(input," ");
			pch=strtok(NULL," ");
			if (pch==NULL)
				printf("\nNo limiting lag specified!\n");
			else
			{
				lag=atof(pch);
				if(autoflag==0)
				{
					autoflag=1;
					autoWindow=cpgopen("/xs");
					cpgask(0);
				}
				cpgslct(autoWindow);
				averageCrossCorrelate(filename_GSB,filename_GWB,last_gsb,last_gwb,s,folda_gsb,folda_gwb,pstart,pend,pulseStart, pulseEnd,phase_gsb,phase_gwb,intg_gsb,intg_gwb,fold,t_gsb,t_gwb,lag,readData_gsb,readData_gwb,gsbDelay,gwbDelay,strong_GSB,avg_crossCorr);
			}
		}*/

		//else if(input[0]=='\n') continue;//printf("\nenter!");

		else
		{
			printf("\nInvalid option. Try again.\n\n");
		}
        
	}
	fclose (pulsesFile);
	fclose (vismicro);
	cpgend();
	return EXIT_SUCCESS;
	}
	default:
	{
		printf("\nWrong input.Bbye\n");
		return 0;
	}
	}
}
