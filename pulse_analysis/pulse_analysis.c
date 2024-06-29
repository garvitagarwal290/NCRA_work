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
#include<stdbool.h>

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
	cpgsubp(1,1);
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
void findWidths(float *rf1Data, float *rf2Data, int arraySize_rf1, int arraySize_rf2, float *widths, double t_rf1, double t_rf2, int intg_rf1, int intg_rf2)
{	//
	float lag=5, refWidth=0.5;			
	int autosize_rf1=(int)(lag/(t_rf1*intg_rf1)),autosize_rf2=(int)(lag/(t_rf2*intg_rf2)),i,j;
	float* autocorr_rf1=(float*)malloc(autosize_rf1*sizeof(float));
	float* autocorr_rf2=(float*)malloc(autosize_rf2*sizeof(float));
	float* lags_rf1=(float*)malloc(autosize_rf1*sizeof(float));
	float* lags_rf2=(float*)malloc(autosize_rf2*sizeof(float));
		
	float sum=0;
	float autoNorm=0;
	float slope=0;
	float refCorr_rf1=0,refCorr_rf2=0;
	
	for(i=0;i<autosize_rf1;i++)
	{
		sum=0;
		for(j=0;j<arraySize_rf1;j++)
		{
			if(j+i>=arraySize_rf1)
				sum+=rf1Data[j]*rf1Data[j+i-arraySize_rf1];
			else
				sum+=rf1Data[j]*rf1Data[j+i];
		}
		
		autocorr_rf1[i]=sum;
		if(i==0)
			autoNorm=autocorr_rf1[i];

		autocorr_rf1[i]=autocorr_rf1[i]/autoNorm;
		lags_rf1[i]=(float)(t_rf1*intg_rf1*i);
	}
	for(i=0;i<autosize_rf2;i++)
	{
		sum=0;
		for(j=0;j<arraySize_rf2;j++)
		{
			if(j+i>=arraySize_rf2)
				sum+=rf2Data[j]*rf2Data[j+i-arraySize_rf2];
			else
				sum+=rf2Data[j]*rf2Data[j+i];
		}
		
		autocorr_rf2[i]=sum;
		if(i==0)
			autoNorm=autocorr_rf2[i];

		autocorr_rf2[i]=autocorr_rf2[i]/autoNorm;
		lags_rf2[i]=(float)(t_rf2*intg_rf2*i);
	}
	printf("\n");
	for(i=0;i<autosize_rf1;i++)
	{
		if(autocorr_rf1[i]<=autocorr_rf1[i+1] && autocorr_rf1[i+1]<=autocorr_rf1[i+2] && autocorr_rf1[i+2]<=autocorr_rf1[i+3] && autocorr_rf1[i]<=autocorr_rf1[i-1] && autocorr_rf1[i-1]<=autocorr_rf1[i-2] && autocorr_rf1[i-2]<=autocorr_rf1[i-3])
		{
			printf("Minima found at _rf1 lag %f ms\n",lags_rf1[i]);
			refCorr_rf1=autocorr_rf1[0]-0.5*(autocorr_rf1[0]-autocorr_rf1[i]);
			for(j=0;j<i;j++)
			{
				if(autocorr_rf1[j]<refCorr_rf1)
				{
					slope=(autocorr_rf1[j-1]-autocorr_rf1[j])/(lags_rf1[j-1]-lags_rf1[j]);
					widths[0]=2*((refCorr_rf1-autocorr_rf1[j])/slope+lags_rf1[j]);
					break;
				}
			}
			break;
		}
		if(i==autosize_rf1-1)
			widths[0]=-1;
	}
	
	for(i=0;i<autosize_rf2;i++)
	{
		if(autocorr_rf2[i]<=autocorr_rf2[i+1] && autocorr_rf2[i+1]<=autocorr_rf2[i+2] && autocorr_rf2[i+2]<=autocorr_rf2[i+3] && autocorr_rf2[i]<=autocorr_rf2[i-1] && autocorr_rf2[i-1]<=autocorr_rf2[i-2] && autocorr_rf2[i-2]<=autocorr_rf2[i-3])
		{
			printf("Minima found at _rf2 lag %f ms\n",lags_rf2[i]);
			refCorr_rf2=autocorr_rf2[0]-0.5*(autocorr_rf2[0]-autocorr_rf2[i]);
			for(j=0;j<i;j++)
			{
				if(autocorr_rf2[j]<refCorr_rf2)
				{
					slope=(autocorr_rf2[j-1]-autocorr_rf2[j])/(lags_rf2[j-1]-lags_rf2[j]);
					widths[1]=2*((refCorr_rf2-autocorr_rf2[j])/slope+lags_rf2[j]);
					break;
				}
			}
			break;
		}
		if(i==autosize_rf2-1)
			widths[1]=-1;
	}
}

void findWidths_derSG(float *rf1Data, float *rf2Data, int arraySize_rf1, int arraySize_rf2, float *widths, double t_rf1, double t_rf2, int intg_rf1, int intg_rf2, double *sg_coeff, int sg_coeff_size)
{
	float lag=5, refWidth=0.5;
	int autosize_rf1=(int)(lag/(t_rf1*intg_rf1)),autosize_rf2=(int)(lag/(t_rf2*intg_rf2)),i,j;
	float* autocorr_rf1=(float*)malloc(autosize_rf1*sizeof(float));
	float* autocorr_rf2=(float*)malloc(autosize_rf2*sizeof(float));
	float* lags_rf1=(float*)malloc(autosize_rf1*sizeof(float));
	float* lags_rf2=(float*)malloc(autosize_rf2*sizeof(float));
		
	float sum=0;
	float autoNorm=0;
	float slope=0;
	float refCorr_rf1=0,refCorr_rf2=0;
	
	for(i=0;i<autosize_rf1;i++)
	{
		sum=0;
		for(j=0;j<arraySize_rf1;j++)
		{
			if(j+i>=arraySize_rf1)
				sum+=rf1Data[j]*rf1Data[j+i-arraySize_rf1];
			else
				sum+=rf1Data[j]*rf1Data[j+i];
		}
		
		autocorr_rf1[i]=sum;
		if(i==0)
			autoNorm=autocorr_rf1[i];

		autocorr_rf1[i]=autocorr_rf1[i]/autoNorm;
		lags_rf1[i]=(float)(t_rf1*intg_rf1*i);
	}
	for(i=0;i<autosize_rf2;i++)
	{
		sum=0;
		for(j=0;j<arraySize_rf2;j++)
		{
			if(j+i>=arraySize_rf2)
				sum+=rf2Data[j]*rf2Data[j+i-arraySize_rf2];
			else
				sum+=rf2Data[j]*rf2Data[j+i];
		}
		
		autocorr_rf2[i]=sum;
		if(i==0)
			autoNorm=autocorr_rf2[i];

		autocorr_rf2[i]=autocorr_rf2[i]/autoNorm;
		lags_rf2[i]=(float)(t_rf2*intg_rf2*i);
	}
	printf("\n");
	
	autocorr_rf2[0]=autocorr_rf2[1]+autocorr_rf2[1]-autocorr_rf2[2];
	autocorr_rf1[0]=autocorr_rf1[1]+autocorr_rf1[1]-autocorr_rf1[2];
	
	float *autocorr_rf1_SGder=(float*)malloc(autosize_rf1*sizeof(float));
	float *autocorr_rf2_SGder=(float*)malloc(autosize_rf2*sizeof(float));	
	
	float *autocorr_rf1_SG=(float*)malloc(autosize_rf1*sizeof(float));
	float *autocorr_rf2_SG=(float*)malloc(autosize_rf2*sizeof(float));
	
	for(i=0;i<autosize_rf1;i++)
	{
		autocorr_rf1_SGder[i]=0;
		autocorr_rf1_SG[i]=0;
		
		for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
		{
			if(j<0)
			{
				autocorr_rf1_SGder[i]+=autocorr_rf1[-j]*(j-i);
				
				if(j-i<0)
				{
					autocorr_rf1_SG[i]+=autocorr_rf1[-j]*sg_coeff[i-j];
				}
				else
				{
					autocorr_rf1_SG[i]+=autocorr_rf1[-j]*sg_coeff[j-i];
				}
			}	
			else
			{
				autocorr_rf1_SGder[i]+=autocorr_rf1[j]*(j-i);
				
				if(j-i<0)
				{
					autocorr_rf1_SG[i]+=autocorr_rf1[j]*sg_coeff[i-j];
				}
				else
				{
					autocorr_rf1_SG[i]+=autocorr_rf1[j]*sg_coeff[j-i];
				}
			}
		}
	}
	
	for(i=0;i<autosize_rf2;i++)
	{
		autocorr_rf2_SGder[i]=0;
		autocorr_rf2_SG[i]=0;
		
		for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
		{
			if(j<0)
			{
				autocorr_rf2_SGder[i]+=autocorr_rf2[-j]*(j-i);
				
				if(j-i<0)
				{
					autocorr_rf2_SG[i]+=autocorr_rf2[-j]*sg_coeff[i-j];
				}
				else
				{
					autocorr_rf2_SG[i]+=autocorr_rf2[-j]*sg_coeff[j-i];
				}
			}	
			else
			{
				autocorr_rf2_SGder[i]+=autocorr_rf2[j]*(j-i);
				
				if(j-i<0)
				{
					autocorr_rf2_SG[i]+=autocorr_rf2[j]*sg_coeff[i-j];
				}
				else
				{
					autocorr_rf2_SG[i]+=autocorr_rf2[j]*sg_coeff[j-i];
				}
			}
		}
	}
	
	for(i=0;i<autosize_rf1;i++)
	{
		autocorr_rf1[i]=autocorr_rf1_SG[i];
	}
	for(i=0;i<autosize_rf2;i++)
	{
		autocorr_rf2[i]=autocorr_rf2_SG[i];
	}
	
	
	printf("\nFinding minima from Savitzky Golay smoothed derivative\n");
	
	float slope_der,minima,autocorrminVal;
	for(i=0;i<autosize_rf1-1;i++)
	{
		if(autocorr_rf1_SGder[i]<0 && autocorr_rf1_SGder[i+1]>0)
		{
			slope_der=(autocorr_rf1_SGder[i+1]-autocorr_rf1_SGder[i])/(lags_rf1[i+1]-lags_rf1[i]);
			minima=lags_rf1[i]-autocorr_rf1_SGder[i]/slope_der;
			slope=(autocorr_rf1[i+1]-autocorr_rf1[i])/(lags_rf1[i+1]-lags_rf1[i]);
			autocorrminVal=autocorr_rf1[i]+slope*(minima-lags_rf1[i]);
			printf("Minima found at _rf1 lag %f ms with ACF _rf1 = %f\n",minima,autocorrminVal);
			refCorr_rf1=autocorr_rf1[0]-0.5*(autocorr_rf1[0]-autocorrminVal);
			for(j=0;j<i;j++)
			{
				if(autocorr_rf1[j]<refCorr_rf1 && j>0)
				{
					slope=(autocorr_rf1[j-1]-autocorr_rf1[j])/(lags_rf1[j-1]-lags_rf1[j]);
					widths[0]=2*((refCorr_rf1-autocorr_rf1[j])/slope+lags_rf1[j]);
					break;
				}
			}
			break;
		}
		if(i==autosize_rf1-2)
			widths[0]=-1;
	}
	for(i=0;i<autosize_rf2-1;i++)
	{
		if(autocorr_rf2_SGder[i]<0 && autocorr_rf2_SGder[i+1]>0)
		{
			slope_der=(autocorr_rf2_SGder[i+1]-autocorr_rf2_SGder[i])/(lags_rf2[i+1]-lags_rf2[i]);
			minima=lags_rf2[i]-autocorr_rf2_SGder[i]/slope_der;
			slope=(autocorr_rf2[i+1]-autocorr_rf2[i])/(lags_rf2[i+1]-lags_rf2[i]);
			autocorrminVal=autocorr_rf2[i]+slope*(minima-lags_rf2[i]);
			printf("Minima found at _rf2 lag %f ms with ACF _rf2 = %f\n",minima,autocorrminVal);
			refCorr_rf2=autocorr_rf2[0]-0.5*(autocorr_rf2[0]-autocorrminVal);
			for(j=0;j<i;j++)
			{
				if(autocorr_rf2[j]<refCorr_rf2 && j>0)
				{
					slope=(autocorr_rf2[j-1]-autocorr_rf2[j])/(lags_rf2[j-1]-lags_rf2[j]);
					widths[1]=2*((refCorr_rf2-autocorr_rf2[j])/slope+lags_rf2[j]);
					break;
				}
			}
			break;
		}
		if(i==autosize_rf2-2)
			widths[1]=-1;
	}
	
	free(autocorr_rf1);
	free(autocorr_rf2);
	free(autocorr_rf1_SG);
	free(autocorr_rf2_SG);
	free(autocorr_rf1_SGder);
	free(autocorr_rf2_SGder);
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
{	//interPolFlag = 1 means _rf2 series needs to be interpolated, and vice versa; LongRes is longer time series
	//As we donot have the same reolution for _rf1 and _rf2, i cases of cross-correlation the data needs to be intrapolated so that both can be artificially made of the same resolution.
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
{	//interPolFlag = 1 means _rf2 series needs to be interpolated, and vice versa; LongRes is longer time series
	double shortSlope;
	int i=0, shortPos;
	
	for(i=0;i<numLongElements;i++)
	{
		shortPos=(int)((i*longRes)/shortRes);
		shortSlope=(shortSeries[shortPos+1]-shortSeries[shortPos])/shortRes;
		interPolate[i]=shortSlope*(i*longRes-shortRes*shortPos)+shortSeries[shortPos];
	}
}

void detectQuasiPeriod(float *corrdata_rf1, float *corrdata_rf2, int corrsize_rf1, int corrsize_rf2, float *peaks, double t_rf1, double t_rf2, int intg_rf1, int intg_rf2, float threshFreq_rf1, float threshFreq_rf2)
{
	int minExpo_rf1=(int)floor((double)log((double)corrsize_rf1)/log(2))+1,minExpo_rf2=(int)floor((double)log((double)corrsize_rf2)/log(2))+1, fftsize_rf1=(int)pow((double)2,minExpo_rf1),fftsize_rf2=(int)pow((double)2,minExpo_rf2),i;
	float Nyqfreq_rf1=(float)1/(2*t_rf1*intg_rf1),Nyqfreq_rf2=(float)1/(2*t_rf2*intg_rf2),max_rf1freq=0,max_rf2freq=0,max_rf1=0,max_rf2=0;
	printf("\n_rf1 Sample size is %d. _rf1 FFT size is %d\n_rf2 Sample size is %d. _rf2 FFT size is %d",corrsize_rf1,fftsize_rf1,corrsize_rf2,fftsize_rf2);
	
	double *FFTUse_rf1=(double*)malloc(fftsize_rf1*sizeof(double));	
	double *FFTUse_rf2=(double*)malloc(fftsize_rf2*sizeof(double));	
	
	float rmsRes_rf1=0, rmsRes_rf2=0;
	
	for(i=0;i<corrsize_rf1;i++)
		rmsRes_rf1+=(corrdata_rf1[i]*corrdata_rf1[i]);
	for(i=0;i<corrsize_rf2;i++)
		rmsRes_rf2+=(corrdata_rf2[i]*corrdata_rf2[i]);		
		
	for(i=0;i<fftsize_rf1;i++)
	{
		if (i<corrsize_rf1)
			FFTUse_rf1[i]=corrdata_rf1[i];
		else
			FFTUse_rf1[i]=0;
	}
	
	for(i=0;i<fftsize_rf2;i++)
	{
		if (i<corrsize_rf2)
			FFTUse_rf2[i]=corrdata_rf2[i];
		else
			FFTUse_rf2[i]=0;
	}
	
	fftw_complex spec_rf1[fftsize_rf1/2+1];
	fftw_complex spec_rf2[fftsize_rf2/2+1];	
			
	fftw_plan plan_rf1, plan_rf2;
	plan_rf1=fftw_plan_dft_r2c_1d(fftsize_rf1,FFTUse_rf1,spec_rf1,FFTW_ESTIMATE);
	plan_rf2=fftw_plan_dft_r2c_1d(fftsize_rf2,FFTUse_rf2,spec_rf2,FFTW_ESTIMATE);
			
	fftw_execute_dft_r2c(plan_rf1,FFTUse_rf1,spec_rf1);
	fftw_execute_dft_r2c(plan_rf2,FFTUse_rf2,spec_rf2);	
	
	float *pow_rf1=(float*)malloc((fftsize_rf1/2+1)*sizeof(float));
	float *pow_rf2=(float*)malloc((fftsize_rf2/2+1)*sizeof(float));
	
	float *freq_rf1=(float*)malloc((fftsize_rf1/2+1)*sizeof(float));
	float *freq_rf2=(float*)malloc((fftsize_rf2/2+1)*sizeof(float));
	
	for(i=0;i<fftsize_rf1/2+1;i++)
	{
		pow_rf1[i]=(float)(spec_rf1[i][0]*spec_rf1[i][0]+spec_rf1[i][1]*spec_rf1[i][1])/(rmsRes_rf1*fftsize_rf1/2);
		freq_rf1[i]=2.0*(float)i/fftsize_rf1*Nyqfreq_rf1;		
	}		
	for(i=0;i<fftsize_rf2/2+1;i++)
	{
		pow_rf2[i]=(float)(spec_rf2[i][0]*spec_rf2[i][0]+spec_rf2[i][1]*spec_rf2[i][1])/(rmsRes_rf2*fftsize_rf2/2);
		freq_rf2[i]=2.0*(float)i/fftsize_rf2*Nyqfreq_rf2;		
	}		
	
	float minPow=pow_rf1[0],maxPow=0,totSqPow_rf1=0,totSqPow_rf2=0,totPow_rf1=0,totPow_rf2=0;
	int normRMS_rf1=0,normRMS_rf2=0;
	
	for(i=0;i<fftsize_rf1/2+1;i++)
	{
		if(pow_rf1[i]>maxPow)
			maxPow=pow_rf1[i];
		if(pow_rf1[i]<minPow)
			minPow=pow_rf1[i];
		
		if(freq_rf1[i]>threshFreq_rf1)
		{
			totPow_rf1+=pow_rf1[i];
			totSqPow_rf1+=(pow_rf1[i]*pow_rf1[i]);
			normRMS_rf1++;
		}
	}
	for(i=0;i<fftsize_rf2/2+1;i++)
	{
		if(pow_rf2[i]>maxPow)
			maxPow=pow_rf2[i];
		if(pow_rf2[i]<minPow)
			minPow=pow_rf2[i];
		
		if(freq_rf2[i]>threshFreq_rf2)
		{
			totPow_rf2+=pow_rf2[i];
			totSqPow_rf2+=(pow_rf2[i]*pow_rf2[i]);
			normRMS_rf2++;
		}
	}
	
	float rmsPow_rf1=sqrt(totSqPow_rf1/normRMS_rf1-(totPow_rf1/normRMS_rf1)*(totPow_rf1/normRMS_rf1));
	float rmsPow_rf2=sqrt(totSqPow_rf2/normRMS_rf2-(totPow_rf2/normRMS_rf2)*(totPow_rf2/normRMS_rf2));
	
	for(i=0;i<fftsize_rf1/2+1;i++)
	{
		if(pow_rf1[i]>pow_rf1[i-1] && pow_rf1[i-1]>pow_rf1[i-2] && pow_rf1[i]>pow_rf1[i+1] && pow_rf1[i+1]>pow_rf1[i+2] && freq_rf1[i]>threshFreq_rf1)
		{
			if(pow_rf1[i]>powerThreshold*rmsPow_rf1 && pow_rf1[i]>max_rf1)
			{
				max_rf1freq=freq_rf1[i];
				max_rf1=pow_rf1[i];
			}
		}	
	}
	for(i=0;i<fftsize_rf2/2+1;i++)
	{
		if(pow_rf2[i]>pow_rf2[i-1] && pow_rf2[i-1]>pow_rf2[i-2] && pow_rf2[i]>pow_rf2[i+1] && pow_rf2[i+1]>pow_rf2[i+2] && freq_rf2[i]>threshFreq_rf2)
		{
			if(pow_rf2[i]>powerThreshold*rmsPow_rf2 && pow_rf2[i]>max_rf2)
			{
				max_rf2freq=freq_rf2[i];
				max_rf2=pow_rf2[i];
			}
		}	
	}
	
	if(max_rf1freq!=0)
		peaks[0]=1.0/max_rf1freq;
	else
		peaks[0]=-1;
	
	if(max_rf2freq!=0)
		peaks[1]=1.0/max_rf2freq;
	else
		peaks[1]=-1;
}

void smoothSubPowerSpectrum(char *filename_rf1, char *filename_rf2, char *FFTName, double *phase_rf1, double *phase_rf2, double s, double *folda_rf1, double *folda_rf2, int last_rf1, int last_rf2,double pstart, double pend, double pulseStart, double pulseEnd, int intg_rf1, int intg_rf2, int pulseNum, double fold, double t_rf1, double t_rf2, double smoothDur_rf1, double smoothDur_rf2, double rf1Delay, double rf2Delay,float *readData_rf1,float *readData_rf2, double _rf1Width, double _rf2Width,int option)
//
{
	switch (option)
	{
	case 2:
	{
	FILE *data_rf1, *data_rf2, *FFTFile;
	data_rf1=fopen(filename_rf1,"rb");
	data_rf2=fopen(filename_rf2,"rb");
	FFTFile=fopen(FFTName,"w");
	int i,numshift_rf1, numshift_rf2,j;
	int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), skipnum_rf2=(int)((pulseNum-1)*(fold/t_rf2)), sampleDelay_rf1=rf1Delay/t_rf1*1000, sampleDelay_rf2=rf2Delay/t_rf2*1000;

	numshift_rf1=(int)(last_rf1*s);
	numshift_rf2=(int)(last_rf2*s);	
	fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);
	fseek(data_rf2,(skipnum_rf2+sampleDelay_rf2+numshift_rf2)*sizeof(float),SEEK_SET);	
	readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1);
	readFloatToDouble(last_rf2,readData_rf2,folda_rf2,data_rf2);
	
	int corrsize_rf1, corrsize_rf2, corrpos=0;
	corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
	corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *corrdata_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	float *phaseData_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *phaseData_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
	double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),range_rf1=max(folda_rf1,last_rf1)-baseline_rf1,baseline_rf2=meanBase(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2),range_rf2=max(folda_rf2,last_rf2)-baseline_rf2;
	for(i=0;i<last_rf1;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
		{
			dpoint+=folda_rf1[j];
			ppoint+=phase_rf1[j];
		}
		if(j==last_rf1)
			break;
		ppoint=ppoint/intg_rf1;
		dpoint=dpoint/intg_rf1;
		dpoint=(dpoint-baseline_rf1)/range_rf1;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
		{
			corrdata_rf1[corrpos]=dpoint;
			phaseData_rf1[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_rf1;
	}
	corrpos=0;
	for(i=0;i<last_rf2;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf2 && j<last_rf2;j++)
		{
			dpoint+=folda_rf2[j];
			ppoint+=phase_rf2[j];
		}
		if(j==last_rf2)
			break;
		ppoint=ppoint/intg_rf2;
		dpoint=dpoint/intg_rf2;
		dpoint=(dpoint-baseline_rf2)/range_rf2;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf2)
		{
			corrdata_rf2[corrpos]=dpoint;
			phaseData_rf2[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_rf2;
	}	
	
	int smoothWindow_rf1, smoothWindow_rf2;	
	float* smooth_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float* smooth_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
	smoothWindow_rf1=(int)(smoothDur_rf1/(t_rf1*intg_rf1));
	smoothWindow_rf2=(int)(smoothDur_rf2/(t_rf2*intg_rf2));
	printf("\nSmoothing window size is %d bins for _rf1; %d bins for _rf2.\n",smoothWindow_rf1,smoothWindow_rf2);
	smoothPulse(corrdata_rf1,smooth_rf1,corrsize_rf1,smoothWindow_rf1);
	smoothPulse(corrdata_rf2,smooth_rf2,corrsize_rf2,smoothWindow_rf2);
	
	if(_rf1Width==0 || _rf2Width==0)
	{
		printf("\nError! Widths cannot be zero!\n");
		return;
	}
	
	float threshFreq_rf1=1.0/_rf1Width*subPulseFreqThreshold;
	float threshFreq_rf2=1.0/_rf2Width*subPulseFreqThreshold;
	printf("\nFrequency thresholds set at _rf1 %f kHz and _rf2 %f kHz\n",threshFreq_rf1,threshFreq_rf2);
	
	float maxSmooth=0,minSmooth=smooth_rf1[0],maxRes=0,minRes=0,rmsRes_rf1=0,rmsRes_rf2=0;
	
	for(i=0;i<corrsize_rf1;i++)
	{
		corrdata_rf1[i]=corrdata_rf1[i]-smooth_rf1[i];
		rmsRes_rf1+=corrdata_rf1[i]*corrdata_rf1[i];
		
		if(smooth_rf1[i]>maxSmooth)
			maxSmooth=smooth_rf1[i];

		if(smooth_rf1[i]<minSmooth)
			minSmooth=smooth_rf1[i];
		
		if(corrdata_rf1[i]>maxRes)
			maxRes=corrdata_rf1[i];
			
		if(corrdata_rf1[i]<minRes)
			minRes=corrdata_rf1[i];
	}
	for(i=0;i<corrsize_rf2;i++)
	{
		corrdata_rf2[i]=corrdata_rf2[i]-smooth_rf2[i];
		rmsRes_rf2+=corrdata_rf2[i]*corrdata_rf2[i];
		
		if(smooth_rf2[i]>maxSmooth)
			maxSmooth=smooth_rf2[i];

		if(smooth_rf2[i]<minSmooth)
			minSmooth=smooth_rf2[i];
		
		if(corrdata_rf2[i]>maxRes)
			maxRes=corrdata_rf2[i];
			
		if(corrdata_rf2[i]<minRes)
			minRes=corrdata_rf2[i];
	}
	
	int minExpo_rf1=(int)floor((double)log((double)corrsize_rf1)/log(2))+1,minExpo_rf2=(int)floor((double)log((double)corrsize_rf2)/log(2))+1, fftsize_rf1=(int)pow((double)2,minExpo_rf1),fftsize_rf2=(int)pow((double)2,minExpo_rf2);
	float Nyqfreq_rf1=(float)1/(2*t_rf1*intg_rf1),Nyqfreq_rf2=(float)1/(2*t_rf2*intg_rf2),max_rf1freq=0,max_rf2freq=0,max_rf1=0,max_rf2=0;
	printf("\n_rf1 Sample size is %d. _rf1 FFT size is %d\n_rf2 Sample size is %d. _rf2 FFT size is %d",corrsize_rf1,fftsize_rf1,corrsize_rf2,fftsize_rf2);
	
	double *FFTUse_rf1=(double*)malloc(fftsize_rf1*sizeof(double));	
	double *FFTUse_rf2=(double*)malloc(fftsize_rf2*sizeof(double));	
	
	for(i=0;i<fftsize_rf1;i++)
	{
		if (i<corrsize_rf1)
			FFTUse_rf1[i]=corrdata_rf1[i];
		else
			FFTUse_rf1[i]=0;
	}
	
	for(i=0;i<fftsize_rf2;i++)
	{
		if (i<corrsize_rf2)
			FFTUse_rf2[i]=corrdata_rf2[i];
		else
			FFTUse_rf2[i]=0;
	}
	
	fftw_complex spec_rf1[fftsize_rf1/2+1];
	fftw_complex spec_rf2[fftsize_rf2/2+1];	
			
	fftw_plan plan_rf1, plan_rf2;
	plan_rf1=fftw_plan_dft_r2c_1d(fftsize_rf1,FFTUse_rf1,spec_rf1,FFTW_ESTIMATE);
	plan_rf2=fftw_plan_dft_r2c_1d(fftsize_rf2,FFTUse_rf2,spec_rf2,FFTW_ESTIMATE);
			
	fftw_execute_dft_r2c(plan_rf1,FFTUse_rf1,spec_rf1);
	fftw_execute_dft_r2c(plan_rf2,FFTUse_rf2,spec_rf2);	
	
	float *pow_rf1=(float*)malloc((fftsize_rf1/2+1)*sizeof(float));
	float *pow_rf2=(float*)malloc((fftsize_rf2/2+1)*sizeof(float));
	
	float *freq_rf1=(float*)malloc((fftsize_rf1/2+1)*sizeof(float));
	float *freq_rf2=(float*)malloc((fftsize_rf2/2+1)*sizeof(float));
	
	for(i=0;i<fftsize_rf1/2+1;i++)
	{
		pow_rf1[i]=(float)(spec_rf1[i][0]*spec_rf1[i][0]+spec_rf1[i][1]*spec_rf1[i][1])/(rmsRes_rf1*fftsize_rf1/2);
		freq_rf1[i]=2.0*(float)i/fftsize_rf1*Nyqfreq_rf1;		
		
		fprintf(FFTFile,"%f\t%f\n",freq_rf1[i],pow_rf1[i]);
	}		
	for(i=0;i<fftsize_rf2/2+1;i++)
	{
		pow_rf2[i]=(float)(spec_rf2[i][0]*spec_rf2[i][0]+spec_rf2[i][1]*spec_rf2[i][1])/(rmsRes_rf2*fftsize_rf2/2);
		freq_rf2[i]=2.0*(float)i/fftsize_rf2*Nyqfreq_rf2;		
		
		fprintf(FFTFile,"%f\t%f\n",freq_rf2[i],pow_rf2[i]);
	}		
	
	fclose(FFTFile);
	
	float minPow=pow_rf1[0],maxPow=0,totSqPow_rf1=0,totSqPow_rf2=0,totPow_rf1=0,totPow_rf2=0;
	int normRMS_rf1=0,normRMS_rf2=0;
	
	for(i=0;i<fftsize_rf1/2+1;i++)
	{
		if(pow_rf1[i]>maxPow)
			maxPow=pow_rf1[i];
		if(pow_rf1[i]<minPow)
			minPow=pow_rf1[i];
		
		if(freq_rf1[i]>threshFreq_rf1)
		{
			totPow_rf1+=pow_rf1[i];
			totSqPow_rf1+=(pow_rf1[i]*pow_rf1[i]);
			normRMS_rf1++;
		}
	}
	for(i=0;i<fftsize_rf2/2+1;i++)
	{
		if(pow_rf2[i]>maxPow)
			maxPow=pow_rf2[i];
		if(pow_rf2[i]<minPow)
			minPow=pow_rf2[i];
		
		if(freq_rf2[i]>threshFreq_rf2)
		{
			totPow_rf2+=pow_rf2[i];
			totSqPow_rf2+=(pow_rf2[i]*pow_rf2[i]);
			normRMS_rf2++;
		}
	}
	
	float rmsPow_rf1=sqrt(totSqPow_rf1/normRMS_rf1-(totPow_rf1/normRMS_rf1)*(totPow_rf1/normRMS_rf1));
	float rmsPow_rf2=sqrt(totSqPow_rf2/normRMS_rf2-(totPow_rf2/normRMS_rf2)*(totPow_rf2/normRMS_rf2));
	
	printf("\nRMSPOW_rf1=%f\tRMSPOW_rf2=%f\n",rmsPow_rf1,rmsPow_rf2);
	
	for(i=0;i<fftsize_rf1/2+1;i++)
	{
		if(pow_rf1[i]>pow_rf1[i-1] && pow_rf1[i-1]>pow_rf1[i-2] && pow_rf1[i]>pow_rf1[i+1] && pow_rf1[i+1]>pow_rf1[i+2] && freq_rf1[i]>threshFreq_rf1)
		{
			if(pow_rf1[i]>powerThreshold*rmsPow_rf1 && pow_rf1[i]>max_rf1)
			{
				max_rf1freq=freq_rf1[i];
				max_rf1=pow_rf1[i];
			}
		}
		freq_rf1[i]=(float)log10(2.0*(double)i/fftsize_rf1*Nyqfreq_rf1);		
	}
	for(i=0;i<fftsize_rf2/2+1;i++)
	{
		if(pow_rf2[i]>pow_rf2[i-1] && pow_rf2[i-1]>pow_rf2[i-2] && pow_rf2[i]>pow_rf2[i+1] && pow_rf2[i+1]>pow_rf2[i+2] && freq_rf2[i]>threshFreq_rf2)
		{
			if(pow_rf2[i]>powerThreshold*rmsPow_rf2 && pow_rf2[i]>max_rf2)
			{
				max_rf2freq=freq_rf2[i];
				max_rf2=pow_rf2[i];
			}
		}
		freq_rf2[i]=(float)log10(2.0*(double)i/fftsize_rf2*Nyqfreq_rf2);		
	}
	
	if(max_rf1freq!=0)
		printf("\nPeriod _rf1 found at %f ms",1.0/max_rf1freq);
	else
		printf("\nPeriod _rf1 not found");
	if(max_rf2freq!=0)
		printf("\nPeriod _rf2 found at %f ms",1.0/max_rf2freq);
	else
		printf("\nPeriod _rf2 not found");
			
	
	cpgsubp(1,3);
        cpgsci(1);
        cpgenv(pstart,pend,minRes,maxRes,0,0);
       	char title[50];
	sprintf(title,"Residuals for pulse number %d",pulseNum);
	cpglab("Phase","Residual Intensity",title);
	cpgsci(2);
        cpgline(corrsize_rf1,phaseData_rf1,corrdata_rf1);
        cpgsci(3);
        cpgline(corrsize_rf2,phaseData_rf2,corrdata_rf2);

        cpgpanl(1,1);
        cpgsci(1);
        cpgenv(pstart,pend,minSmooth,maxSmooth,0,0);
	sprintf(title,"Smoothed pulse for pulse number %d",pulseNum);
	cpglab("Phase","Intensity",title);
	cpgsci(2);
        cpgline(corrsize_rf1,phaseData_rf1,smooth_rf1);
        cpgsci(3);
        cpgline(corrsize_rf2,phaseData_rf2,smooth_rf2);
        
        cpgpanl(1,2);
        cpgsci(1);
        cpgenv(-1,freq_rf1[fftsize_rf1/2],minPow,maxPow,0,10);
        sprintf(title,"Power spectrum of residuals for pulse number %d",pulseNum); 
	cpglab("Frequency (kHz)","Normalised Power",title);
	cpgsci(2);
        cpgline(fftsize_rf1/2,freq_rf1,pow_rf1);
        cpgsci(3);
        cpgline(fftsize_rf2/2,freq_rf2,pow_rf2);
        
	free(phaseData_rf1);
	free(phaseData_rf2);
	free(corrdata_rf1);
	free(corrdata_rf2);
	free(freq_rf1);
	free(pow_rf1);
	free(freq_rf2);
	free(pow_rf2);
	break;
	}
	case 1:
	{
	FILE *data_rf1, *FFTFile;
	data_rf1=fopen(filename_rf1,"rb");
	FFTFile=fopen(FFTName,"w");
	int i,numshift_rf1,j;
	int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), sampleDelay_rf1=rf1Delay/t_rf1*1000;

	numshift_rf1=(int)(last_rf1*s);
	fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);	
	readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1);
	
	int corrsize_rf1, corrpos=0;
	corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *phaseData_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	
	double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),range_rf1=max(folda_rf1,last_rf1)-baseline_rf1;
	for(i=0;i<last_rf1;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
		{
			dpoint+=folda_rf1[j];
			ppoint+=phase_rf1[j];
		}
		if(j==last_rf1)
			break;
		ppoint=ppoint/intg_rf1;
		dpoint=dpoint/intg_rf1;
		dpoint=(dpoint-baseline_rf1)/range_rf1;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
		{
			corrdata_rf1[corrpos]=dpoint;
			phaseData_rf1[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_rf1;
	}
	corrpos=0;
	
	int smoothWindow_rf1;	
	float* smooth_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	
	smoothWindow_rf1=(int)(smoothDur_rf1/(t_rf1*intg_rf1));
	printf("\nSmoothing window size is %d bins \n",smoothWindow_rf1);
	smoothPulse(corrdata_rf1,smooth_rf1,corrsize_rf1,smoothWindow_rf1);
	
	if(_rf1Width==0)
	{
		printf("\nError! Widths cannot be zero!\n");
		return;
	}
	
	float threshFreq_rf1=1.0/_rf1Width*subPulseFreqThreshold;
	printf("\nFrequency thresholds set at _rf1 %f kHz \n",threshFreq_rf1);
	
	float maxSmooth=0,minSmooth=smooth_rf1[0],maxRes=0,minRes=0,rmsRes_rf1=0;
	
	for(i=0;i<corrsize_rf1;i++)
	{
		corrdata_rf1[i]=corrdata_rf1[i]-smooth_rf1[i];
		rmsRes_rf1+=corrdata_rf1[i]*corrdata_rf1[i];
		
		if(smooth_rf1[i]>maxSmooth)
			maxSmooth=smooth_rf1[i];

		if(smooth_rf1[i]<minSmooth)
			minSmooth=smooth_rf1[i];
		
		if(corrdata_rf1[i]>maxRes)
			maxRes=corrdata_rf1[i];
			
		if(corrdata_rf1[i]<minRes)
			minRes=corrdata_rf1[i];
	}
	
	int minExpo_rf1=(int)floor((double)log((double)corrsize_rf1)/log(2))+1, fftsize_rf1=(int)pow((double)2,minExpo_rf1);
	float Nyqfreq_rf1=(float)1/(2*t_rf1*intg_rf1),max_rf1freq=0,max_rf1=0;
	printf("\nSample size is %d. FFT size is %d\n",corrsize_rf1,fftsize_rf1);
	
	double *FFTUse_rf1=(double*)malloc(fftsize_rf1*sizeof(double));		
	
	for(i=0;i<fftsize_rf1;i++)
	{
		if (i<corrsize_rf1)
			FFTUse_rf1[i]=corrdata_rf1[i];
		else
			FFTUse_rf1[i]=0;
	}
	
	fftw_complex spec_rf1[fftsize_rf1/2+1];
			
	fftw_plan plan_rf1;
	plan_rf1=fftw_plan_dft_r2c_1d(fftsize_rf1,FFTUse_rf1,spec_rf1,FFTW_ESTIMATE);
			
	fftw_execute_dft_r2c(plan_rf1,FFTUse_rf1,spec_rf1);	
	
	float *pow_rf1=(float*)malloc((fftsize_rf1/2+1)*sizeof(float));
	
	float *freq_rf1=(float*)malloc((fftsize_rf1/2+1)*sizeof(float));
	
	for(i=0;i<fftsize_rf1/2+1;i++)
	{
		pow_rf1[i]=(float)(spec_rf1[i][0]*spec_rf1[i][0]+spec_rf1[i][1]*spec_rf1[i][1])/(rmsRes_rf1*fftsize_rf1/2);
		freq_rf1[i]=2.0*(float)i/fftsize_rf1*Nyqfreq_rf1;		
		
		fprintf(FFTFile,"%f\t%f\n",freq_rf1[i],pow_rf1[i]);
	}			
	
	fclose(FFTFile);
	
	float minPow=pow_rf1[0],maxPow=0,totSqPow_rf1=0,totPow_rf1=0;
	int normRMS_rf1=0;
	
	for(i=0;i<fftsize_rf1/2+1;i++)
	{
		if(pow_rf1[i]>maxPow)
			maxPow=pow_rf1[i];
		if(pow_rf1[i]<minPow)
			minPow=pow_rf1[i];
		
		if(freq_rf1[i]>threshFreq_rf1)
		{
			totPow_rf1+=pow_rf1[i];
			totSqPow_rf1+=(pow_rf1[i]*pow_rf1[i]);
			normRMS_rf1++;
		}
	}
	
	float rmsPow_rf1=sqrt(totSqPow_rf1/normRMS_rf1-(totPow_rf1/normRMS_rf1)*(totPow_rf1/normRMS_rf1));
	
	printf("\nRMSPOW_rf1=%f\n",rmsPow_rf1);
	
	for(i=0;i<fftsize_rf1/2+1;i++)
	{
		if(pow_rf1[i]>pow_rf1[i-1] && pow_rf1[i-1]>pow_rf1[i-2] && pow_rf1[i]>pow_rf1[i+1] && pow_rf1[i+1]>pow_rf1[i+2] && freq_rf1[i]>threshFreq_rf1)
		{
			if(pow_rf1[i]>powerThreshold*rmsPow_rf1 && pow_rf1[i]>max_rf1)
			{
				max_rf1freq=freq_rf1[i];
				max_rf1=pow_rf1[i];
			}
		}
		freq_rf1[i]=(float)log10(2.0*(double)i/fftsize_rf1*Nyqfreq_rf1);		
	}
	
	if(max_rf1freq!=0)
		printf("\nPeriod _rf1 found at %f ms",1.0/max_rf1freq);
	else
		printf("\nPeriod _rf1 not found");
			
	
	cpgsubp(1,3);
        cpgsci(1);
        cpgenv(pstart,pend,minRes,maxRes,0,0);
       	char title[50];
	sprintf(title,"Residuals for pulse number %d",pulseNum);
	cpglab("Phase","Residual Intensity",title);
	cpgsci(2);
        cpgline(corrsize_rf1,phaseData_rf1,corrdata_rf1);

        cpgpanl(1,1);
        cpgsci(1);
        cpgenv(pstart,pend,minSmooth,maxSmooth,0,0);
	sprintf(title,"Smoothed pulse for pulse number %d",pulseNum);
	cpglab("Phase","Intensity",title);
	cpgsci(2);
        cpgline(corrsize_rf1,phaseData_rf1,smooth_rf1);
        
        cpgpanl(1,2);
        cpgsci(1);
        cpgenv(-1,freq_rf1[fftsize_rf1/2],minPow,maxPow,0,10);
        sprintf(title,"Power spectrum of residuals for pulse number %d",pulseNum); 
	cpglab("Frequency (kHz)","Normalised Power",title);
	cpgsci(2);
        cpgline(fftsize_rf1/2,freq_rf1,pow_rf1);
        
	free(phaseData_rf1);
	free(corrdata_rf1);
	free(freq_rf1);
	free(pow_rf1);
	}
	}
}
void giveRelStrengthSmoothed(float *corrdata_rf1, float *corrdata_rf2, float *smooth_rf1, float *smooth_rf2, int corrsize_rf1, int corrsize_rf2, float *relArray)
{
	int i,j;
	float totPow_rf1=0, totPow_rf2=0;
	float totSmooth_rf1=0, totSmooth_rf2=0;
	
	for(i=0;i<corrsize_rf1;i++)
	{	
		totPow_rf1+=corrdata_rf1[i]*corrdata_rf1[i];
		totSmooth_rf1+=smooth_rf1[i]*smooth_rf1[i];
	}
	for(i=0;i<corrsize_rf2;i++)
	{	
		totPow_rf2+=corrdata_rf2[i]*corrdata_rf2[i];
		totSmooth_rf2+=smooth_rf2[i]*smooth_rf2[i];
	}
	
	relArray[0]=totPow_rf1/totSmooth_rf1;
	relArray[1]=totPow_rf2/totSmooth_rf2;
	relArray[2]=totSmooth_rf1;
	relArray[3]=totSmooth_rf2;
	
	printf("\nRatios are _rf1 = %f \t _rf2 = %f\n",relArray[0],relArray[1]);
}

void smoothSubSGACF_intrapolated(char *filename_rf1, char *filename_rf2, char *autoName, char *SG_coeffs, double *phase_rf1, double *phase_rf2, double s, double *folda_rf1, double *folda_rf2, int last_rf1, int last_rf2,double pstart, double pend, double pulseStart, double pulseEnd, int intg_rf1, int intg_rf2, int pulseNum, double fold, double t_rf1, double t_rf2, double lag, double smoothDur_rf1, double smoothDur_rf2, double rf1Delay, double rf2Delay,float *readData_rf1,float *readData_rf2, int sg_coeff_size,int option,char *outfile)
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
	
			FILE *data_rf1, *data_rf2, *autoFile;
			data_rf1=fopen(filename_rf1,"rb");
			data_rf2=fopen(filename_rf2,"rb");
			autoFile=fopen(autoName,"w");
			int numshift_rf1, numshift_rf2,j;
			int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), skipnum_rf2=(int)((pulseNum-1)*(fold/t_rf2)), sampleDelay_rf1=rf1Delay/t_rf1*1000, sampleDelay_rf2=rf2Delay/t_rf2*1000;

			numshift_rf1=(int)(last_rf1*s);
			numshift_rf2=(int)(last_rf2*s);	
			fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);
			fseek(data_rf2,(skipnum_rf2+sampleDelay_rf2+numshift_rf2)*sizeof(float),SEEK_SET);	
			readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1);
			readFloatToDouble(last_rf2,readData_rf2,folda_rf2,data_rf2);
	
			int corrsize_rf1, corrsize_rf2, corrpos=0;
			corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
			corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
			float dpoint=0,ppoint=0;
			float *corrdata_rf11=(float*)malloc(corrsize_rf1*sizeof(float));
			float *corrdata_rf21=(float*)malloc(corrsize_rf2*sizeof(float));
			float *phaseData_rf11=(float*)malloc(corrsize_rf1*sizeof(float));
			float *phaseData_rf21=(float*)malloc(corrsize_rf2*sizeof(float));
	
			double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),range_rf1=max(folda_rf1,last_rf1)-baseline_rf1,baseline_rf2=meanBase(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2),range_rf2=max(folda_rf2,last_rf2)-baseline_rf2;
			for(i=0;i<last_rf1;)
			{
				dpoint=0;
				ppoint=0;
				for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
				{
					dpoint+=folda_rf1[j];
					ppoint+=phase_rf1[j];
				}
				if(j==last_rf1)
					break;
				ppoint=ppoint/intg_rf1;
				dpoint=dpoint/intg_rf1;
				dpoint=(dpoint-baseline_rf1)/range_rf1;
		
				if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
				{
					corrdata_rf11[corrpos]=dpoint;
					phaseData_rf11[corrpos]=ppoint;
					corrpos++;
				}
				i+=intg_rf1;
			}
			corrpos=0;
			for(i=0;i<last_rf2;)
			{
				dpoint=0;
				ppoint=0;
				for(j=i;j<i+intg_rf2 && j<last_rf2;j++)
				{
					dpoint+=folda_rf2[j];
					ppoint+=phase_rf2[j];
				}
				if(j==last_rf2)
					break;
				ppoint=ppoint/intg_rf2;
				dpoint=dpoint/intg_rf2;
				dpoint=(dpoint-baseline_rf2)/range_rf2;
		
				if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf2)
				{
					corrdata_rf21[corrpos]=dpoint;
					phaseData_rf21[corrpos]=ppoint;
					corrpos++;
				}
				i+=intg_rf2;
			}
			corrpos=0;

			float *crossArray1, *crossArray2,*phaseArray1,*phaseArray2;
			float t_corr; int corrsize;
			if(corrsize_rf1>corrsize_rf2)
			{
				crossArray1=corrdata_rf11;
				crossArray2=(float*)malloc(corrsize_rf1*sizeof(float));
				phaseArray1=phaseData_rf11;
				phaseArray2=(float*)malloc(corrsize_rf1*sizeof(float));
				interpolate_float(corrdata_rf21,corrdata_rf11,t_rf2*intg_rf2,t_rf1*intg_rf1,crossArray2,corrsize_rf1);
				interpolate_float(phaseData_rf21,phaseData_rf11,t_rf2*intg_rf2,t_rf1*intg_rf1,phaseArray2,corrsize_rf1);
				t_rf2=t_rf1*intg_rf1;
				corrsize_rf2=corrsize_rf1;
				printf("Yo");
			}
			else
			{
				crossArray2=corrdata_rf21;
				crossArray1=(float*)malloc(corrsize_rf2*sizeof(float));
				phaseArray2=phaseData_rf21;
				phaseArray1=(float*)malloc(corrsize_rf2*sizeof(float));
				interpolate_float(corrdata_rf11,corrdata_rf21,t_rf1*intg_rf1,t_rf2*intg_rf2,crossArray1,corrsize_rf2);
				interpolate_float(phaseData_rf11,phaseData_rf21,t_rf1*intg_rf1,t_rf2*intg_rf2,phaseArray1,corrsize_rf2);
				t_rf1=t_rf2*intg_rf2;
				corrsize_rf1=corrsize_rf2;
			} 

			float *corrdata_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
			float *corrdata_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
			float *phaseData_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
			float *phaseData_rf2=(float*)malloc(corrsize_rf2*sizeof(float));

			for (i=0;i<corrsize_rf1;i++)
			{
				corrdata_rf1[i]=crossArray1[i];
				//printf("%lf\n",corrdata_rf1[i]);
				//printf("%lf\n",crossArray1[i]);
				//printf("%lf\n",corrdata_rf11[i]);
			}
			
			for (i=0;i<corrsize_rf2;i++)
			{
				corrdata_rf2[i]=crossArray2[i];
			}

			for (i=0;i<corrsize_rf1;i++)
			{
				phaseData_rf1[i]=phaseArray1[i];
			}
			
			for (i=0;i<corrsize_rf2;i++)
			{
				phaseData_rf2[i]=phaseArray2[i];
			}		
	
			int smoothWindow_rf1, smoothWindow_rf2;	
			float* smooth_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
			float* smooth_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
			smoothWindow_rf1=(int)(smoothDur_rf1/(t_rf1*intg_rf1));
			smoothWindow_rf2=(int)(smoothDur_rf2/(t_rf2*intg_rf2));
			printf("\nSmoothing window size is %d bins for _rf1; %d bins for _rf2.\n",smoothWindow_rf1,smoothWindow_rf2);
			smoothPulse(corrdata_rf1,smooth_rf1,corrsize_rf1,smoothWindow_rf1);
			smoothPulse(corrdata_rf2,smooth_rf2,corrsize_rf2,smoothWindow_rf2);
	
			float maxSmooth=0,minSmooth=smooth_rf1[0],maxRes=0,minRes=0,rmsRes_rf1=0,rmsRes_rf2=0;
	
			for(i=0;i<corrsize_rf1;i++)
			{
				corrdata_rf1[i]=corrdata_rf1[i]-smooth_rf1[i];
				rmsRes_rf1+=corrdata_rf1[i]*corrdata_rf1[i];
					
				if(smooth_rf1[i]>maxSmooth)
					maxSmooth=smooth_rf1[i];
		
				if(smooth_rf1[i]<minSmooth)
					minSmooth=smooth_rf1[i];
		
				if(corrdata_rf1[i]>maxRes)
					maxRes=corrdata_rf1[i];
					
				if(corrdata_rf1[i]<minRes)
					minRes=corrdata_rf1[i];
			}
			printf("Power at _rf1 = %f \n",rmsRes_rf1);
			for(i=0;i<corrsize_rf2;i++)
			{
				corrdata_rf2[i]=corrdata_rf2[i]-smooth_rf2[i];
				rmsRes_rf2+=corrdata_rf2[i]*corrdata_rf2[i];
		
				if(smooth_rf2[i]>maxSmooth)
					maxSmooth=smooth_rf2[i];

				if(smooth_rf2[i]<minSmooth)
					minSmooth=smooth_rf2[i];
		
				if(corrdata_rf2[i]>maxRes)
					maxRes=corrdata_rf2[i];
			
				if(corrdata_rf2[i]<minRes)
					minRes=corrdata_rf2[i];
			}
			printf("Power at _rf2 = %f \n",rmsRes_rf2);
		
			for(i=0;i<last_rf1;i++)
				folda_rf1[i]=(folda_rf1[i]-baseline_rf1)/range_rf1;
			for(i=0;i<last_rf2;i++)
				folda_rf2[i]=(folda_rf2[i]-baseline_rf2)/range_rf2;
		
			float stdev_rf1=stdevBaseIntg(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1,intg_rf1),stdev_rf2=stdevBaseIntg(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2,intg_rf2);
	
			int normRes_rf1=corrsize_rf1, normRes_rf2=corrsize_rf2;
			rmsRes_rf1/=(normRes_rf1*stdev_rf1*stdev_rf1);
			rmsRes_rf2/=(normRes_rf2*stdev_rf2*stdev_rf2);
	
			rmsRes_rf1=sqrt((double)rmsRes_rf1);
			rmsRes_rf2=sqrt((double)rmsRes_rf2);

			float maxDev_rf1=maxMod(corrdata_rf1,corrsize_rf1),maxDev_rf2=maxMod(corrdata_rf2,corrsize_rf2);
			float maxDevRMS_rf1=maxDev_rf1/stdev_rf1,maxDevRMS_rf2=maxDev_rf2/stdev_rf2;
	
			printf("\nRMS deviation of _rf1 = %f\t_rf2 = %f for STDEV_rf1=%f\tSTDEV_rf2=%f",rmsRes_rf1,rmsRes_rf2,stdev_rf1,stdev_rf2);
			printf("\nMax_RMS_Dev of _rf1 = %f\t_rf2 = %f",maxDevRMS_rf1,maxDevRMS_rf2);
			float autoWidths[2];
			findWidths_derSG(corrdata_rf1,corrdata_rf2,corrsize_rf1,corrsize_rf2,autoWidths,t_rf1,t_rf2,intg_rf1,intg_rf2,sg_coeff,sg_coeff_size);
			printf("\nWidths detected _rf1 = %f ms\t_rf2 = %f ms\n",autoWidths[0],autoWidths[1]);	
	
			int autosize_rf1=(int)(lag/(t_rf1*intg_rf1));
			int autosize_rf2=(int)(lag/(t_rf2*intg_rf2));
			float* autocorr_rf1=(float*)malloc(autosize_rf1*sizeof(float));
			float* autocorr_rf2=(float*)malloc(autosize_rf2*sizeof(float));
			float* lags_rf1=(float*)malloc(autosize_rf1*sizeof(float));
			float* lags_rf2=(float*)malloc(autosize_rf2*sizeof(float));	
			
			float maxAuto=0,minAuto=1;
		
			float sum=0;
			float autoNorm=0;
	
			for(i=0;i<autosize_rf1;i++)
			{
				sum=0;
				for(j=0;j<corrsize_rf1;j++)
				{
					if(j+i>=corrsize_rf1)
					{
						sum+=corrdata_rf1[j]*corrdata_rf1[j+i-corrsize_rf1];
					}
					else
					{
						sum+=corrdata_rf1[j]*corrdata_rf1[j+i];
					}
				}
		
				autocorr_rf1[i]=sum;
		
				if(i==0)
					autoNorm=autocorr_rf1[i];

				autocorr_rf1[i]=autocorr_rf1[i]/autoNorm;
		
				if(autocorr_rf1[i]>maxAuto)
				{
					maxAuto=autocorr_rf1[i];
				}
				if(autocorr_rf1[i]<minAuto)
				{
					minAuto=autocorr_rf1[i];
				}
				
				lags_rf1[i]=(t_rf1*intg_rf1*i);
				fprintf(autoFile,"%f\t%f\n",lags_rf1[i],autocorr_rf1[i]);
			}
			for(i=0;i<autosize_rf2;i++)
			{
				sum=0;
				for(j=0;j<corrsize_rf2;j++)
				{
					if(j+i>=corrsize_rf2)
					{
						sum+=corrdata_rf2[j]*corrdata_rf2[j+i-corrsize_rf2];
					}
					else
					{
						sum+=corrdata_rf2[j]*corrdata_rf2[j+i];
					}
				}
		
				autocorr_rf2[i]=sum;
		
				if(i==0)
					autoNorm=autocorr_rf2[i];

				autocorr_rf2[i]=autocorr_rf2[i]/autoNorm;
		
				if(autocorr_rf2[i]>maxAuto)
				{
					maxAuto=autocorr_rf2[i];
				}
				if(autocorr_rf2[i]<minAuto)
				{
					minAuto=autocorr_rf2[i];
				}
		
				lags_rf2[i]=(t_rf2*intg_rf2*i);
				fprintf(autoFile,"%f\t%f\n",lags_rf2[i],autocorr_rf2[i]);
			}

	
	
			float *autocorr_rf1_SG=(float*)malloc(autosize_rf1*sizeof(float));
			float *autocorr_rf2_SG=(float*)malloc(autosize_rf2*sizeof(float));
	
			for(i=0;i<autosize_rf1;i++)
			{
				autocorr_rf1_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[-j]*sg_coeff[i-j];
						}
						else
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[j]*sg_coeff[i-j];
						}
						else
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[j]*sg_coeff[j-i];
						}
					}
				}
			}

	
	
			for(i=0;i<autosize_rf2;i++)
			{
				autocorr_rf2_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							autocorr_rf2_SG[i]+=autocorr_rf2[-j]*sg_coeff[i-j];
						}
						else
						{
							autocorr_rf2_SG[i]+=autocorr_rf2[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							autocorr_rf2_SG[i]+=autocorr_rf2[j]*sg_coeff[i-j];
						}
						else
						{
							autocorr_rf2_SG[i]+=autocorr_rf2[j]*sg_coeff[j-i];
						}
					}
				}
			}
	
			for(i=0;i<autosize_rf1;i++)
			{
				autocorr_rf1[i]=autocorr_rf1_SG[i];
			}
			for(i=0;i<autosize_rf2;i++)
			{
				autocorr_rf2[i]=autocorr_rf2_SG[i];
			}

        		fclose(autoFile);
        		fclose(data_rf1);
        		fclose(data_rf2);
        
        		cpgsubp(1,3);

        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       			char title[50];
			sprintf(title,"Residuals for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Phase","Residual Intensity",title);
			cpgsci(2);
        		cpgline(corrsize_rf1,phaseData_rf1,corrdata_rf1);
        		cpgsci(3);
        		cpgline(corrsize_rf2,phaseData_rf2,corrdata_rf2);


        		cpgpanl(1,1);
        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
			sprintf(title,"Smoothed pulse for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Phase","Intensity of smoothed part",title);
			cpgsci(2);
        		cpgline(corrsize_rf1,phaseData_rf1,smooth_rf1);
        		cpgsci(3);
        		cpgline(corrsize_rf2,phaseData_rf2,smooth_rf2);


        
        		cpgpanl(1,2);
        		cpgsci(1);
        		cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        		sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Intensity. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Time lag (ms)","Autocorrelation amplitude",title);
			cpgsci(2);
        		cpgline(autosize_rf1,lags_rf1,autocorr_rf1);
        		cpgsci(3);
        		cpgline(autosize_rf2,lags_rf2,autocorr_rf2);

	//Vincross 
			/*FILE *crfile=fopen("crossresidual.txt","w");

			for (i=0 ; i<corrsize_rf1 ; i++)
			{
				corrdata_rf1[i]=corrdata_rf1[i]*50;
				//fprintf(garbage,"%.10lf ",corrdata_rf1[i]);
			}
			//fprintf(garbage,"\n\n\n\n");

			for (i=0 ; i<corrsize_rf2 ; i++)
			{
				corrdata_rf2[i]=corrdata_rf2[i]*50;
				//fprintf(garbage,"%.10lf ",corrdata_rf2[i]);
			}

	

			double *crossArray1, *crossArray2;
			double t_corr; int corrsize;
			if(corrsize_rf1>corrsize_rf2)
			{
				crossArray1=corrdata_rf1;
				crossArray2=(double*)malloc(corrsize_rf1*sizeof(double));
				interpolate(corrdata_rf2,corrdata_rf1,t_rf2*intg_rf2,t_rf1*intg_rf1,crossArray2,corrsize_rf1);
				t_corr=t_rf1*intg_rf1;
				corrsize=corrsize_rf1;
			}
			else
			{
				crossArray2=corrdata_rf2;
				crossArray1=(double*)malloc(corrsize_rf2*sizeof(double));
				interpolate(corrdata_rf1,corrdata_rf2,t_rf1*intg_rf1,t_rf2*intg_rf2,crossArray1,corrsize_rf2);
				t_corr=t_rf2*intg_rf2;
				corrsize=corrsize_rf2;
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


			
			for (i=0 ; i<corrsize_rf1 ; i++)
			{
				corrdata_rf1[i]=corrdata_rf1[i]/50;
			}

			for (i=0 ; i<corrsize_rf2 ; i++)
			{
				corrdata_rf2[i]=corrdata_rf2[i]/50;
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

				int corrsize_rf21=(int)(last_rf2/intg_rf2*(pend-pstart));
				int corrsize_rf11=(int)(last_rf1/intg_rf1*(pend-pstart));

				int pos=0;

				float dpoint=0,ppoint=0;
				float *corrdata_rf11=(float*)malloc(corrsize_rf11*sizeof(float));
				float *corrdata_rf21=(float*)malloc(corrsize_rf21*sizeof(float));
				float *phaseData_rf11=(float*)malloc(corrsize_rf11*sizeof(float));
				float *phaseData_rf21=(float*)malloc(corrsize_rf21*sizeof(float));
				
				for(i=0;i<corrsize_rf1;i++)
				{
					if (phaseData_rf1[i]>=pstart && phaseData_rf1[i]<=pend)
					{
						corrdata_rf11[pos]=corrdata_rf1[i]/50;
						phaseData_rf11[pos]=phaseData_rf1[i];
					}
				}
				corrpos=0;
				for(i=0;i<corrsize_rf2;i++)
				{
					if (phaseData_rf2[i]>=pstart && phaseData_rf2[i]<=pend)
					{
						corrdata_rf21[pos]=corrdata_rf2[i]/50;
						phaseData_rf21[pos]=phaseData_rf2[i];
					}
				}
				corrpos=0;
				for(i=0;i<corrsize_rf11;i++)
				{
					if(corrdata_rf11[i]>maxRes)
						maxRes=corrdata_rf11[i];
						
					if(corrdata_rf11[i]<minRes)
						minRes=corrdata_rf11[i];
					printf("\n%lf\t%lf",phaseData_rf11[i],corrdata_rf1[i]);
				}	
				
				for(i=0;i<corrsize_rf21;i++)
				{
					if(corrdata_rf21[i]>maxRes)
						maxRes=corrdata_rf21[i];
				
					if(corrdata_rf21[i]<minRes)
						minRes=corrdata_rf21[i];
				}

				cpgeras;

				cpgsubp(1,4);

				cpgsci(1);
        			cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       				char title[50];
				sprintf(title,"Residuals for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Phase","Residual Intensity",title);
				cpgsci(2);
        			cpgline(corrsize_rf11,phaseData_rf11,corrdata_rf11);
        			cpgsci(3);
        			cpgline(corrsize_rf21,phaseData_rf21,corrdata_rf21);

				cpgpanl(1,1);
        			cpgsci(1);
        			cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
				sprintf(title,"Smoothed pulse for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Phase","Intensity of smoothed part",title);
				cpgsci(2);
        			cpgline(corrsize_rf1,phaseData_rf1,smooth_rf1);
        			cpgsci(3);
        			cpgline(corrsize_rf2,phaseData_rf2,smooth_rf2);
	
	
        	
        			cpgpanl(1,2);
        			cpgsci(1);
        			cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        			sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Intensity. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Time lag (ms)","Autocorrelation amplitude",title);
				cpgsci(2);
        			cpgline(autosize_rf1,lags_rf1,autocorr_rf1);
        			cpgsci(3);
        			cpgline(autosize_rf2,lags_rf2,autocorr_rf2);

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

			


			free(phaseData_rf1);
			free(phaseData_rf2);
			free(corrdata_rf1);
			free(corrdata_rf2);
			free(lags_rf1);
			free(lags_rf2);
			free(autocorr_rf1);
			free(autocorr_rf2);
			free(autocorr_rf1_SG);
			free(autocorr_rf2_SG);
}

void smoothSubSGACF(char *filename_rf1, char *filename_rf2, char *autoName, char *SG_coeffs, double *phase_rf1, double *phase_rf2, double s, double *folda_rf1, double *folda_rf2, int last_rf1, int last_rf2,double pstart, double pend, double pulseStart, double pulseEnd, int intg_rf1, int intg_rf2, int pulseNum, double fold, double t_rf1, double t_rf2, double lag, double smoothDur_rf1, double smoothDur_rf2, double rf1Delay, double rf2Delay,float *readData_rf1,float *readData_rf2, int sg_coeff_size,int option,char *outfile)
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
	
			FILE *data_rf1, *data_rf2, *autoFile;
			data_rf1=fopen(filename_rf1,"rb");
			data_rf2=fopen(filename_rf2,"rb");
			autoFile=fopen(autoName,"w");
			int numshift_rf1, numshift_rf2,j;
			int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), skipnum_rf2=(int)((pulseNum-1)*(fold/t_rf2)), sampleDelay_rf1=rf1Delay/t_rf1*1000, sampleDelay_rf2=rf2Delay/t_rf2*1000;

			numshift_rf1=(int)(last_rf1*s);
			numshift_rf2=(int)(last_rf2*s);	
			fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);
			fseek(data_rf2,(skipnum_rf2+sampleDelay_rf2+numshift_rf2)*sizeof(float),SEEK_SET);	
			readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1);
			readFloatToDouble(last_rf2,readData_rf2,folda_rf2,data_rf2);
	
			int corrsize_rf1, corrsize_rf2, corrpos=0;
			corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
			corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
			float dpoint=0,ppoint=0;
			float *corrdata_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
			float *corrdata_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
			float *phaseData_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
			float *phaseData_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
			double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),range_rf1=max(folda_rf1,last_rf1)-baseline_rf1,baseline_rf2=meanBase(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2),range_rf2=max(folda_rf2,last_rf2)-baseline_rf2;
			for(i=0;i<last_rf1;)
			{
				dpoint=0;
				ppoint=0;
				for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
				{
					dpoint+=folda_rf1[j];
					ppoint+=phase_rf1[j];
				}
				if(j==last_rf1)
					break;
				ppoint=ppoint/intg_rf1;
				dpoint=dpoint/intg_rf1;
				dpoint=(dpoint-baseline_rf1)/range_rf1;
		
				if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
				{
					corrdata_rf1[corrpos]=dpoint;
					phaseData_rf1[corrpos]=ppoint;
					corrpos++;
				}
				i+=intg_rf1;
			}
			corrpos=0;
			for(i=0;i<last_rf2;)
			{
				dpoint=0;
				ppoint=0;
				for(j=i;j<i+intg_rf2 && j<last_rf2;j++)
				{
					dpoint+=folda_rf2[j];
					ppoint+=phase_rf2[j];
				}
				if(j==last_rf2)
					break;
				ppoint=ppoint/intg_rf2;
				dpoint=dpoint/intg_rf2;
				dpoint=(dpoint-baseline_rf2)/range_rf2;
		
				if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf2)
				{
					corrdata_rf2[corrpos]=dpoint;
					phaseData_rf2[corrpos]=ppoint;
					corrpos++;
				}
				i+=intg_rf2;
			}
			corrpos=0;	
	
			int smoothWindow_rf1, smoothWindow_rf2;	
			float* smooth_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
			float* smooth_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
			smoothWindow_rf1=(int)(smoothDur_rf1/(t_rf1*intg_rf1));
			smoothWindow_rf2=(int)(smoothDur_rf2/(t_rf2*intg_rf2));
			printf("\nSmoothing window size is %d bins for _rf1; %d bins for _rf2.\n",smoothWindow_rf1,smoothWindow_rf2);
			smoothPulse(corrdata_rf1,smooth_rf1,corrsize_rf1,smoothWindow_rf1);
			smoothPulse(corrdata_rf2,smooth_rf2,corrsize_rf2,smoothWindow_rf2);
	
			float maxSmooth=0,minSmooth=smooth_rf1[0],maxRes=0,minRes=0,rmsRes_rf1=0,rmsRes_rf2=0;
	
			for(i=0;i<corrsize_rf1;i++)
			{
				corrdata_rf1[i]=corrdata_rf1[i]-smooth_rf1[i];
				rmsRes_rf1+=corrdata_rf1[i]*corrdata_rf1[i];
					
				if(smooth_rf1[i]>maxSmooth)
					maxSmooth=smooth_rf1[i];
		
				if(smooth_rf1[i]<minSmooth)
					minSmooth=smooth_rf1[i];
		
				if(corrdata_rf1[i]>maxRes)
					maxRes=corrdata_rf1[i];
					
				if(corrdata_rf1[i]<minRes)
					minRes=corrdata_rf1[i];
			}
			printf("Power at _rf1 = %f \n",rmsRes_rf1);
			for(i=0;i<corrsize_rf2;i++)
			{
				corrdata_rf2[i]=corrdata_rf2[i]-smooth_rf2[i];
				rmsRes_rf2+=corrdata_rf2[i]*corrdata_rf2[i];
		
				if(smooth_rf2[i]>maxSmooth)
					maxSmooth=smooth_rf2[i];

				if(smooth_rf2[i]<minSmooth)
					minSmooth=smooth_rf2[i];
		
				if(corrdata_rf2[i]>maxRes)
					maxRes=corrdata_rf2[i];
			
				if(corrdata_rf2[i]<minRes)
					minRes=corrdata_rf2[i];
			}
			printf("Power at _rf2 = %f \n",rmsRes_rf2);
		
			for(i=0;i<last_rf1;i++)
				folda_rf1[i]=(folda_rf1[i]-baseline_rf1)/range_rf1;
			for(i=0;i<last_rf2;i++)
				folda_rf2[i]=(folda_rf2[i]-baseline_rf2)/range_rf2;
		
			float stdev_rf1=stdevBaseIntg(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1,intg_rf1),stdev_rf2=stdevBaseIntg(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2,intg_rf2);
	
			int normRes_rf1=corrsize_rf1, normRes_rf2=corrsize_rf2;
			rmsRes_rf1/=(normRes_rf1*stdev_rf1*stdev_rf1);
			rmsRes_rf2/=(normRes_rf2*stdev_rf2*stdev_rf2);
	
			rmsRes_rf1=sqrt((double)rmsRes_rf1);
			rmsRes_rf2=sqrt((double)rmsRes_rf2);

			float maxDev_rf1=maxMod(corrdata_rf1,corrsize_rf1),maxDev_rf2=maxMod(corrdata_rf2,corrsize_rf2);
			float maxDevRMS_rf1=maxDev_rf1/stdev_rf1,maxDevRMS_rf2=maxDev_rf2/stdev_rf2;
	
			printf("\nRMS deviation of _rf1 = %f\t_rf2 = %f for STDEV_rf1=%f\tSTDEV_rf2=%f",rmsRes_rf1,rmsRes_rf2,stdev_rf1,stdev_rf2);
			printf("\nMax_RMS_Dev of _rf1 = %f\t_rf2 = %f",maxDevRMS_rf1,maxDevRMS_rf2);
			float autoWidths[2];
			findWidths_derSG(corrdata_rf1,corrdata_rf2,corrsize_rf1,corrsize_rf2,autoWidths,t_rf1,t_rf2,intg_rf1,intg_rf2,sg_coeff,sg_coeff_size);
			printf("\nWidths detected _rf1 = %f ms\t_rf2 = %f ms\n",autoWidths[0],autoWidths[1]);	
	
			int autosize_rf1=(int)(lag/(t_rf1*intg_rf1));
			int autosize_rf2=(int)(lag/(t_rf2*intg_rf2));
			float* autocorr_rf1=(float*)malloc(autosize_rf1*sizeof(float));
			float* autocorr_rf2=(float*)malloc(autosize_rf2*sizeof(float));
			float* lags_rf1=(float*)malloc(autosize_rf1*sizeof(float));
			float* lags_rf2=(float*)malloc(autosize_rf2*sizeof(float));	
			
			float maxAuto=0,minAuto=1;
		
			float sum=0;
			float autoNorm=0;
	
			for(i=0;i<autosize_rf1;i++)
			{
				sum=0;
				for(j=0;j<corrsize_rf1;j++)
				{
					if(j+i>=corrsize_rf1)
					{
					sum+=corrdata_rf1[j]*corrdata_rf1[j+i-corrsize_rf1];
					}
					else
					{
						sum+=corrdata_rf1[j]*corrdata_rf1[j+i];
					}
				}
		
				autocorr_rf1[i]=sum;
		
				if(i==0)
					autoNorm=autocorr_rf1[i];

				autocorr_rf1[i]=autocorr_rf1[i]/autoNorm;
		
				if(autocorr_rf1[i]>maxAuto)
				{
					maxAuto=autocorr_rf1[i];
				}
				if(autocorr_rf1[i]<minAuto)
				{
					minAuto=autocorr_rf1[i];
				}
				
				lags_rf1[i]=(t_rf1*intg_rf1*i);
				fprintf(autoFile,"%f\t%f\n",lags_rf1[i],autocorr_rf1[i]);
			}
			for(i=0;i<autosize_rf2;i++)
			{
				sum=0;
				for(j=0;j<corrsize_rf2;j++)
				{
					if(j+i>=corrsize_rf2)
					{
						sum+=corrdata_rf2[j]*corrdata_rf2[j+i-corrsize_rf2];
					}
					else
					{
						sum+=corrdata_rf2[j]*corrdata_rf2[j+i];
					}
				}
		
				autocorr_rf2[i]=sum;
		
				if(i==0)
					autoNorm=autocorr_rf2[i];

				autocorr_rf2[i]=autocorr_rf2[i]/autoNorm;
		
				if(autocorr_rf2[i]>maxAuto)
				{
					maxAuto=autocorr_rf2[i];
				}
				if(autocorr_rf2[i]<minAuto)
				{
					minAuto=autocorr_rf2[i];
				}
		
				lags_rf2[i]=(t_rf2*intg_rf2*i);
				fprintf(autoFile,"%f\t%f\n",lags_rf2[i],autocorr_rf2[i]);
			}

	
	
			float *autocorr_rf1_SG=(float*)malloc(autosize_rf1*sizeof(float));
			float *autocorr_rf2_SG=(float*)malloc(autosize_rf2*sizeof(float));
	
			for(i=0;i<autosize_rf1;i++)
			{
				autocorr_rf1_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[-j]*sg_coeff[i-j];
						}
						else
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[j]*sg_coeff[i-j];
						}
						else
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[j]*sg_coeff[j-i];
						}
					}
				}
			}

	
	
			for(i=0;i<autosize_rf2;i++)
			{
				autocorr_rf2_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							autocorr_rf2_SG[i]+=autocorr_rf2[-j]*sg_coeff[i-j];
						}
						else
						{
							autocorr_rf2_SG[i]+=autocorr_rf2[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							autocorr_rf2_SG[i]+=autocorr_rf2[j]*sg_coeff[i-j];
						}
						else
						{
							autocorr_rf2_SG[i]+=autocorr_rf2[j]*sg_coeff[j-i];
						}
					}
				}
			}
	
			for(i=0;i<autosize_rf1;i++)
			{
				autocorr_rf1[i]=autocorr_rf1_SG[i];
			}
			for(i=0;i<autosize_rf2;i++)
			{
				autocorr_rf2[i]=autocorr_rf2_SG[i];
			}

        		fclose(autoFile);
        		fclose(data_rf1);
        		fclose(data_rf2);
        
        		cpgsubp(1,4);

        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       			char title[50];
			sprintf(title,"Residuals for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Phase","Residual Intensity",title);
			cpgsci(2);
        		cpgline(corrsize_rf1,phaseData_rf1,corrdata_rf1);
        		cpgsci(3);
        		cpgline(corrsize_rf2,phaseData_rf2,corrdata_rf2);


        		cpgpanl(1,1);
        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
			sprintf(title,"Smoothed pulse for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Phase","Intensity of smoothed part",title);
			cpgsci(2);
        		cpgline(corrsize_rf1,phaseData_rf1,smooth_rf1);
        		cpgsci(3);
        		cpgline(corrsize_rf2,phaseData_rf2,smooth_rf2);


        
        		cpgpanl(1,2);
        		cpgsci(1);
        		cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        		sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Intensity. Red = low frequency, Green = high frequency",pulseNum);
			cpglab("Time lag (ms)","Autocorrelation amplitude",title);
			cpgsci(2);
        		cpgline(autosize_rf1,lags_rf1,autocorr_rf1);
        		cpgsci(3);
        		cpgline(autosize_rf2,lags_rf2,autocorr_rf2);

	//Vincross 
			FILE *crfile=fopen("crossresidual.txt","w");

			/*for (i=0 ; i<corrsize_rf1 ; i++)
			{
				corrdata_rf1[i]=corrdata_rf1[i]*50;
				//fprintf(garbage,"%.10lf ",corrdata_rf1[i]);
			}
			//fprintf(garbage,"\n\n\n\n");

			for (i=0 ; i<corrsize_rf2 ; i++)
			{
				corrdata_rf2[i]=corrdata_rf2[i]*50;
				//fprintf(garbage,"%.10lf ",corrdata_rf2[i]);
			}*/

	

			float *crossArray1, *crossArray2;
			float t_corr; int corrsize;
			if(corrsize_rf1>corrsize_rf2)
			{
				crossArray1=corrdata_rf1;
				crossArray2=(float*)malloc(corrsize_rf1*sizeof(float));
				interpolate_float(corrdata_rf2,corrdata_rf1,t_rf2*intg_rf2,t_rf1*intg_rf1,crossArray2,corrsize_rf1);
				t_corr=t_rf1*intg_rf1;
				corrsize=corrsize_rf1;
			}
			else
			{
				crossArray2=corrdata_rf2;
				crossArray1=(float*)malloc(corrsize_rf2*sizeof(float));
				interpolate_float(corrdata_rf1,corrdata_rf2,t_rf1*intg_rf1,t_rf2*intg_rf2,crossArray1,corrsize_rf2);
				t_corr=t_rf2*intg_rf2;
				corrsize=corrsize_rf2;
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

			
			
			/*for (i=0 ; i<corrsize_rf1 ; i++)
			{
				corrdata_rf1[i]=corrdata_rf1[i]/50;
			}

			for (i=0 ; i<corrsize_rf2 ; i++)
			{
				corrdata_rf2[i]=corrdata_rf2[i]/50;
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

				int corrsize_rf21=(int)(last_rf2/intg_rf2*(pend-pstart));
				int corrsize_rf11=(int)(last_rf1/intg_rf1*(pend-pstart));

				int pos=0;

				float dpoint=0,ppoint=0;
				float *corrdata_rf11=(float*)malloc(corrsize_rf11*sizeof(float));
				float *corrdata_rf21=(float*)malloc(corrsize_rf21*sizeof(float));
				float *phaseData_rf11=(float*)malloc(corrsize_rf11*sizeof(float));
				float *phaseData_rf21=(float*)malloc(corrsize_rf21*sizeof(float));
				
				for(i=0;i<corrsize_rf1;i++)
				{
					if (phaseData_rf1[i]>=pstart && phaseData_rf1[i]<=pend)
					{
						corrdata_rf11[pos]=corrdata_rf1[i]/50;
						phaseData_rf11[pos]=phaseData_rf1[i];
					}
				}
				corrpos=0;
				for(i=0;i<corrsize_rf2;i++)
				{
					if (phaseData_rf2[i]>=pstart && phaseData_rf2[i]<=pend)
					{
						corrdata_rf21[pos]=corrdata_rf2[i]/50;
						phaseData_rf21[pos]=phaseData_rf2[i];
					}
				}
				corrpos=0;
				for(i=0;i<corrsize_rf11;i++)
				{
					if(corrdata_rf11[i]>maxRes)
						maxRes=corrdata_rf11[i];
						
					if(corrdata_rf11[i]<minRes)
						minRes=corrdata_rf11[i];
					printf("\n%lf\t%lf",phaseData_rf11[i],corrdata_rf1[i]);
				}	
				
				for(i=0;i<corrsize_rf21;i++)
				{
					if(corrdata_rf21[i]>maxRes)
						maxRes=corrdata_rf21[i];
				
					if(corrdata_rf21[i]<minRes)
						minRes=corrdata_rf21[i];
				}

				cpgeras;

				cpgsubp(1,4);

				cpgsci(1);
        			cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       				char title[50];
				sprintf(title,"Residuals for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Phase","Residual Intensity",title);
				cpgsci(2);
        			cpgline(corrsize_rf11,phaseData_rf11,corrdata_rf11);
        			cpgsci(3);
        			cpgline(corrsize_rf21,phaseData_rf21,corrdata_rf21);

				cpgpanl(1,1);
        			cpgsci(1);
        			cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
				sprintf(title,"Smoothed pulse for pulse number %d. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Phase","Intensity of smoothed part",title);
				cpgsci(2);
        			cpgline(corrsize_rf1,phaseData_rf1,smooth_rf1);
        			cpgsci(3);
        			cpgline(corrsize_rf2,phaseData_rf2,smooth_rf2);
	
	
        	
        			cpgpanl(1,2);
        			cpgsci(1);
        			cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        			sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Intensity. Red = low frequency, Green = high frequency",pulseNum);
				cpglab("Time lag (ms)","Autocorrelation amplitude",title);
				cpgsci(2);
        			cpgline(autosize_rf1,lags_rf1,autocorr_rf1);
        			cpgsci(3);
        			cpgline(autosize_rf2,lags_rf2,autocorr_rf2);

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

			


			free(phaseData_rf1);
			free(phaseData_rf2);
			free(corrdata_rf1);
			free(corrdata_rf2);
			free(lags_rf1);
			free(lags_rf2);
			free(autocorr_rf1);
			free(autocorr_rf2);
			free(autocorr_rf1_SG);
			free(autocorr_rf2_SG);
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
	
			FILE *data_rf1, *autoFile;
			data_rf1=fopen(filename_rf1,"rb");
			autoFile=fopen(autoName,"w");
			int numshift_rf1,j;
			int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), sampleDelay_rf1=rf1Delay/t_rf1*1000;

			numshift_rf1=(int)(last_rf1*s);
			fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);	
			readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1);
	
			int corrsize_rf1, corrpos=0;
			corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
			float dpoint=0,ppoint=0;
			float *corrdata_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
			float *phaseData_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	
			double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),range_rf1=max(folda_rf1,last_rf1)-baseline_rf1;
			for(i=0;i<last_rf1;)
			{
				dpoint=0;
				ppoint=0;
				for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
				{
					dpoint+=folda_rf1[j];
					ppoint+=phase_rf1[j];
				}
				if(j==last_rf1)
					break;
				ppoint=ppoint/intg_rf1;
				dpoint=dpoint/intg_rf1;
				dpoint=(dpoint-baseline_rf1)/range_rf1;
		
				if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
				{
					corrdata_rf1[corrpos]=dpoint;
					phaseData_rf1[corrpos]=ppoint;
					corrpos++;
				}
				i+=intg_rf1;
			}
			corrpos=0;
	
			int smoothWindow_rf1, smoothWindow_rf2;	
			float* smooth_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	
			smoothWindow_rf1=(int)(smoothDur_rf1/(t_rf1*intg_rf1));
			printf("\nSmoothing window size is %d bins.\n",smoothWindow_rf1);
			smoothPulse(corrdata_rf1,smooth_rf1,corrsize_rf1,smoothWindow_rf1);
	
			float maxSmooth=0,minSmooth=smooth_rf1[0],maxRes=0,minRes=0,rmsRes_rf1=0;
	
			for(i=0;i<corrsize_rf1;i++)
			{
				corrdata_rf1[i]=corrdata_rf1[i]-smooth_rf1[i];
				rmsRes_rf1+=corrdata_rf1[i]*corrdata_rf1[i];
		
				if(smooth_rf1[i]>maxSmooth)
					maxSmooth=smooth_rf1[i];

				if(smooth_rf1[i]<minSmooth)
					minSmooth=smooth_rf1[i];
		
				if(corrdata_rf1[i]>maxRes)
					maxRes=corrdata_rf1[i];
			
				if(corrdata_rf1[i]<minRes)
					minRes=corrdata_rf1[i];
			}
			printf("Power = %f \n",rmsRes_rf1);
	
			for(i=0;i<last_rf1;i++)
				folda_rf1[i]=(folda_rf1[i]-baseline_rf1)/range_rf1;

			float stdev_rf1=stdevBaseIntg(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1,intg_rf1);
	
			int normRes_rf1=corrsize_rf1;
			rmsRes_rf1/=(normRes_rf1*stdev_rf1*stdev_rf1);
	
			rmsRes_rf1=sqrt((double)rmsRes_rf1);

			float maxDev_rf1=maxMod(corrdata_rf1,corrsize_rf1);
			float maxDevRMS_rf1=maxDev_rf1/stdev_rf1;
	
			printf("\nRMS deviation = %f\t for STDEV=%f\t",rmsRes_rf1,stdev_rf1);
			printf("\nMax_RMS_Dev = %f\t",maxDevRMS_rf1);
			float autoWidths[2];
			findWidths_derSG(corrdata_rf1,corrdata_rf1,corrsize_rf1,corrsize_rf1,autoWidths,t_rf1,t_rf1,intg_rf1,intg_rf1,sg_coeff,sg_coeff_size);
			printf("\nWidths detected = %f ms \n",autoWidths[0]);	
	
			int autosize_rf1=(int)(lag/(t_rf1*intg_rf1));
			float* autocorr_rf1=(float*)malloc(autosize_rf1*sizeof(float));
			float* lags_rf1=(float*)malloc(autosize_rf1*sizeof(float));
	
			float maxAuto=0,minAuto=1;
		
			float sum=0;
			float autoNorm=0;
	
			for(i=0;i<autosize_rf1;i++)
			{
				sum=0;
				for(j=0;j<corrsize_rf1;j++)
				{
					if(j+i>=corrsize_rf1)
					{
						sum+=corrdata_rf1[j]*corrdata_rf1[j+i-corrsize_rf1];
					}
					else
					{
						sum+=corrdata_rf1[j]*corrdata_rf1[j+i];
					}
				}
		
				autocorr_rf1[i]=sum;
		
				if(i==0)
					autoNorm=autocorr_rf1[i];

				autocorr_rf1[i]=autocorr_rf1[i]/autoNorm;
		
				if(autocorr_rf1[i]>maxAuto)
				{
					maxAuto=autocorr_rf1[i];
				}
				if(autocorr_rf1[i]<minAuto)
				{
					minAuto=autocorr_rf1[i];
				}
		
				lags_rf1[i]=(t_rf1*intg_rf1*i);
				fprintf(autoFile,"%f\t%f\n",lags_rf1[i],autocorr_rf1[i]);
			}
	
			float *autocorr_rf1_SG=(float*)malloc(autosize_rf1*sizeof(float));
	
			for(i=0;i<autosize_rf1;i++)
			{
				autocorr_rf1_SG[i]=0;
		
				for(j=i-sg_coeff_size;j<=i+sg_coeff_size;j++)
				{
					if(j<0)
					{
				
						if(j-i<0)
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[-j]*sg_coeff[i-j];
						}
						else
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[-j]*sg_coeff[j-i];
						}
					}	
					else
					{
				
						if(j-i<0)
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[j]*sg_coeff[i-j];
						}
						else
						{
							autocorr_rf1_SG[i]+=autocorr_rf1[j]*sg_coeff[j-i];
						}
					}
				}
			}



			for(i=0;i<autosize_rf1;i++)
			{
				autocorr_rf1[i]=autocorr_rf1_SG[i];
			}
	
        		fclose(autoFile);
        		fclose(data_rf1);
        
        		cpgsubp(1,3);
        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       			char title[50];
			sprintf(title,"Residuals for pulse number %d",pulseNum);
			cpglab("Phase","Residual Intensity",title);
			cpgsci(2);
        		cpgline(corrsize_rf1,phaseData_rf1,corrdata_rf1);
        
        		cpgpanl(1,1);
        		cpgsci(1);
        		cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
			sprintf(title,"Smoothed pulse for pulse number %d",pulseNum);
			cpglab("Phase","Intensity of smoothed part",title);
			cpgsci(2);
        		cpgline(corrsize_rf1,phaseData_rf1,smooth_rf1);
        
        		cpgpanl(1,2);
        		cpgsci(1);
        		cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        		sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Intensity",pulseNum);
			cpglab("Time lag (ms)","Autocorrelation amplitude",title);
			cpgsci(2);
        		cpgline(autosize_rf1,lags_rf1,autocorr_rf1);


			free(phaseData_rf1);
			free(corrdata_rf1);
			free(lags_rf1);
			free(autocorr_rf1);
			free(autocorr_rf1_SG);		
			break;
		}
	}
}


void smoothSubACF(char *filename_rf1, char *filename_rf2, char *autoName, double *phase_rf1, double *phase_rf2, double s, double *folda_rf1, double *folda_rf2, int last_rf1, int last_rf2,double pstart, double pend, double pulseStart, double pulseEnd, int intg_rf1, int intg_rf2, int pulseNum, double fold, double t_rf1, double t_rf2, double lag, double smoothDur_rf1, double smoothDur_rf2, double rf1Delay, double rf2Delay,float *readData_rf1,float *readData_rf2)
{
	FILE *data_rf1, *data_rf2, *autoFile;
	data_rf1=fopen(filename_rf1,"rb");
	data_rf2=fopen(filename_rf2,"rb");
	autoFile=fopen(autoName,"w");
	int i,numshift_rf1, numshift_rf2,j;
	int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), skipnum_rf2=(int)((pulseNum-1)*(fold/t_rf2)), sampleDelay_rf1=rf1Delay/t_rf1*1000, sampleDelay_rf2=rf2Delay/t_rf2*1000;

	numshift_rf1=(int)(last_rf1*s);
	numshift_rf2=(int)(last_rf2*s);	
	fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);
	fseek(data_rf2,(skipnum_rf2+sampleDelay_rf2+numshift_rf2)*sizeof(float),SEEK_SET);	
	readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1);
	readFloatToDouble(last_rf2,readData_rf2,folda_rf2,data_rf2);
	
	int corrsize_rf1, corrsize_rf2, corrpos=0;
	corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
	corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *corrdata_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	float *phaseData_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *phaseData_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
	double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),range_rf1=max(folda_rf1,last_rf1)-baseline_rf1,baseline_rf2=meanBase(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2),range_rf2=max(folda_rf2,last_rf2)-baseline_rf2;
	for(i=0;i<last_rf1;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
		{
			dpoint+=folda_rf1[j];
			ppoint+=phase_rf1[j];
		}
		if(j==last_rf1)
			break;
		ppoint=ppoint/intg_rf1;
		dpoint=dpoint/intg_rf1;
		dpoint=(dpoint-baseline_rf1)/range_rf1;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
		{
			corrdata_rf1[corrpos]=dpoint;
			phaseData_rf1[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_rf1;
	}
	corrpos=0;
	for(i=0;i<last_rf2;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf2 && j<last_rf2;j++)
		{
			dpoint+=folda_rf2[j];
			ppoint+=phase_rf2[j];
		}
		if(j==last_rf2)
			break;
		ppoint=ppoint/intg_rf2;
		dpoint=dpoint/intg_rf2;
		dpoint=(dpoint-baseline_rf2)/range_rf2;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf2)
		{
			corrdata_rf2[corrpos]=dpoint;
			phaseData_rf2[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_rf2;
	}	
	
	int smoothWindow_rf1, smoothWindow_rf2;	
	float* smooth_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float* smooth_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
	smoothWindow_rf1=(int)(smoothDur_rf1/(t_rf1*intg_rf1));
	smoothWindow_rf2=(int)(smoothDur_rf2/(t_rf2*intg_rf2));
	printf("\nSmoothing window size is %d bins for _rf1; %d bins for _rf2.\n",smoothWindow_rf1,smoothWindow_rf2);
	smoothPulse(corrdata_rf1,smooth_rf1,corrsize_rf1,smoothWindow_rf1);
	smoothPulse(corrdata_rf2,smooth_rf2,corrsize_rf2,smoothWindow_rf2);
	
	float maxSmooth=0,minSmooth=smooth_rf1[0],maxRes=0,minRes=0,rmsRes_rf1=0,rmsRes_rf2=0;
	
	for(i=0;i<corrsize_rf1;i++)
	{
		corrdata_rf1[i]=corrdata_rf1[i]-smooth_rf1[i];
		rmsRes_rf1+=corrdata_rf1[i]*corrdata_rf1[i];
		
		if(smooth_rf1[i]>maxSmooth)
			maxSmooth=smooth_rf1[i];

		if(smooth_rf1[i]<minSmooth)
			minSmooth=smooth_rf1[i];
		
		if(corrdata_rf1[i]>maxRes)
			maxRes=corrdata_rf1[i];
			
		if(corrdata_rf1[i]<minRes)
			minRes=corrdata_rf1[i];
	}
	for(i=0;i<corrsize_rf2;i++)
	{
		corrdata_rf2[i]=corrdata_rf2[i]-smooth_rf2[i];
		rmsRes_rf2+=corrdata_rf2[i]*corrdata_rf2[i];
		
		if(smooth_rf2[i]>maxSmooth)
			maxSmooth=smooth_rf2[i];

		if(smooth_rf2[i]<minSmooth)
			minSmooth=smooth_rf2[i];
		
		if(corrdata_rf2[i]>maxRes)
			maxRes=corrdata_rf2[i];
			
		if(corrdata_rf2[i]<minRes)
			minRes=corrdata_rf2[i];
	}
	
	for(i=0;i<last_rf1;i++)
		folda_rf1[i]=(folda_rf1[i]-baseline_rf1)/range_rf1;
	for(i=0;i<last_rf2;i++)
		folda_rf2[i]=(folda_rf2[i]-baseline_rf2)/range_rf2;
		
	float stdev_rf1=stdevBaseIntg(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1,intg_rf1),stdev_rf2=stdevBaseIntg(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2,intg_rf2);
	
	int normRes_rf1=corrsize_rf1, normRes_rf2=corrsize_rf2;
	rmsRes_rf1/=(normRes_rf1*stdev_rf1*stdev_rf1);
	rmsRes_rf2/=(normRes_rf2*stdev_rf2*stdev_rf2);
	
	rmsRes_rf1=sqrt((double)rmsRes_rf1);
	rmsRes_rf2=sqrt((double)rmsRes_rf2);
	
	printf("\nRMS deviation of _rf1 = %f\t_rf2 = %f for STDEV_rf1=%f\tSTDEV_rf2=%f",rmsRes_rf1,rmsRes_rf2,stdev_rf1,stdev_rf2);
	
	float autoWidths[2];
	findWidths(corrdata_rf1,corrdata_rf2,corrsize_rf1,corrsize_rf2,autoWidths,t_rf1,t_rf2,intg_rf1,intg_rf2);
	printf("\nWidths detected _rf1 = %f ms\t_rf2 = %f ms\n",autoWidths[0],autoWidths[1]);	
	
	int autosize_rf1=(int)(lag/(t_rf1*intg_rf1));
	int autosize_rf2=(int)(lag/(t_rf2*intg_rf2));
	float* autocorr_rf1=(float*)malloc(autosize_rf1*sizeof(float));
	float* autocorr_rf2=(float*)malloc(autosize_rf2*sizeof(float));
	float* lags_rf1=(float*)malloc(autosize_rf1*sizeof(float));
	float* lags_rf2=(float*)malloc(autosize_rf2*sizeof(float));	
	float maxAuto=0,minAuto=1;
		
	float sum=0;
	float autoNorm=0;
	
	for(i=0;i<autosize_rf1;i++)
	{
		sum=0;
		for(j=0;j<corrsize_rf1;j++)
		{
			if(j+i>=corrsize_rf1)
			{
				sum+=corrdata_rf1[j]*corrdata_rf1[j+i-corrsize_rf1];
			}
			else
			{
				sum+=corrdata_rf1[j]*corrdata_rf1[j+i];
			}
		}
		
		autocorr_rf1[i]=sum;
		
		if(i==0)
			autoNorm=autocorr_rf1[i];

		autocorr_rf1[i]=autocorr_rf1[i]/autoNorm;
		
		if(autocorr_rf1[i]>maxAuto)
		{
			maxAuto=autocorr_rf1[i];
		}
		if(autocorr_rf1[i]<minAuto)
		{
			minAuto=autocorr_rf1[i];
		}
		
		lags_rf1[i]=(t_rf1*intg_rf1*i);
		fprintf(autoFile,"%f\t%f\n",lags_rf1[i],autocorr_rf1[i]);
	}
	for(i=0;i<autosize_rf2;i++)
	{
		sum=0;
		for(j=0;j<corrsize_rf2;j++)
		{
			if(j+i>=corrsize_rf2)
			{
				sum+=corrdata_rf2[j]*corrdata_rf2[j+i-corrsize_rf2];
			}
			else
			{
				sum+=corrdata_rf2[j]*corrdata_rf2[j+i];
			}
		}
		
		autocorr_rf2[i]=sum;
		
		if(i==0)
			autoNorm=autocorr_rf2[i];

		autocorr_rf2[i]=autocorr_rf2[i]/autoNorm;
		
		if(autocorr_rf2[i]>maxAuto)
		{
			maxAuto=autocorr_rf2[i];
		}
		if(autocorr_rf2[i]<minAuto)
		{
			minAuto=autocorr_rf2[i];
		}
		
		lags_rf2[i]=(t_rf2*intg_rf2*i);
		fprintf(autoFile,"%f\t%f\n",lags_rf2[i],autocorr_rf2[i]);
	}

        fclose(autoFile);
        fclose(data_rf1);
        fclose(data_rf2);
        
        cpgsubp(1,3);
        cpgsci(1);
        cpgenv((float)pstart,(float)pend,minRes,maxRes,0,0);
       	char title[50];
	sprintf(title,"Residuals for pulse number %d",pulseNum);
	cpglab("Phase","Residual Stokes parameters",title);
	cpgsci(2);
        cpgline(corrsize_rf1,phaseData_rf1,corrdata_rf1);
        cpgsci(3);
        cpgline(corrsize_rf2,phaseData_rf2,corrdata_rf2);
        
        cpgpanl(1,1);
        cpgsci(1);
        cpgenv((float)pstart,(float)pend,minSmooth,maxSmooth,0,0);
	sprintf(title,"Smoothed pulse for pulse number %d",pulseNum);
	cpglab("Phase","Stokes parameters",title);
	cpgsci(2);
        cpgline(corrsize_rf1,phaseData_rf1,smooth_rf1);
        cpgsci(3);
        cpgline(corrsize_rf2,phaseData_rf2,smooth_rf2);
        
        cpgpanl(1,2);
        cpgsci(1);
        cpgenv(0,(float)lag,minAuto,maxAuto,0,0);
        sprintf(title,"Autocorrelation plot of residuals for pulse number %d in Stokes parameters",pulseNum);
	cpglab("Time lag (ms)","Autocorrelation amplitude",title);
	cpgsci(2);
        cpgline(autosize_rf1,lags_rf1,autocorr_rf1);
        cpgsci(3);
        cpgline(autosize_rf2,lags_rf2,autocorr_rf2);
        
	free(phaseData_rf1);
	free(phaseData_rf2);
	free(corrdata_rf1);
	free(corrdata_rf2);
	//printf("hello3\n");	
	//free(lags_rf1);
	//free(lags_rf2);
	//printf("hello2\n");	
	free(autocorr_rf1);
	free(autocorr_rf2);
}

int microCandidates(char *filename_rf1, char *filename_rf2, FILE *microFile, double *phase_rf1, double *phase_rf2, double s, double *folda_rf1, double *folda_rf2, int last_rf1, int last_rf2,double pstart, double pend, double pulseStart, double pulseEnd, int intg_rf1, int intg_rf2, int pulseNum, double fold, double t_rf1, double t_rf2, double lag, double smoothDur_rf1, double smoothDur_rf2, double rf1Delay, double rf2Delay,float *readData_rf1,float *readData_rf2, float cutOffhighfreq,float cutOfflowfreq, float _rf1Width, float _rf2Width)
{
	FILE *data_rf1, *data_rf2;
	data_rf2=fopen(filename_rf2,"rb");
	data_rf1=fopen(filename_rf1,"rb");

	if(data_rf1==NULL)
	{
		printf("\nCould not open _rf1 file!");
		return 0;
	}
	if(data_rf2==NULL)
	{
		printf("\nCould not open _rf2 file!");
		return 0;
	}
	
	int i,numshift_rf1, numshift_rf2,j;
	int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), skipnum_rf2=(int)((pulseNum-1)*(fold/t_rf2)), sampleDelay_rf1=rf1Delay/t_rf1*1000, sampleDelay_rf2=rf2Delay/t_rf2*1000;

	numshift_rf1=(int)(last_rf1*s);
	numshift_rf2=(int)(last_rf2*s);	
	fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);
	fseek(data_rf2,(skipnum_rf2+sampleDelay_rf2+numshift_rf2)*sizeof(float),SEEK_SET);	
	
	if(readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1)!=1)
		return 0;
	if(readFloatToDouble(last_rf2,readData_rf2,folda_rf2,data_rf2)!=1)
		return 0;
	
	fclose(data_rf1);
	fclose(data_rf2);
	
	int corrsize_rf1, corrsize_rf2, corrpos=0;
	corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
	corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *corrdata_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	float *phaseData_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *phaseData_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
	double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),range_rf1=max(folda_rf1,last_rf1)-baseline_rf1,baseline_rf2=meanBase(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2),range_rf2=max(folda_rf2,last_rf2)-baseline_rf2;
	for(i=0;i<last_rf1;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
		{
			dpoint+=folda_rf1[j];
			ppoint+=phase_rf1[j];
		}
		if(j==last_rf1)
			break;
		ppoint=ppoint/intg_rf1;
		dpoint=dpoint/intg_rf1;
		dpoint=(dpoint-baseline_rf1)/range_rf1;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
		{
			corrdata_rf1[corrpos]=dpoint;
			phaseData_rf1[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_rf1;
	}
	corrpos=0;
	for(i=0;i<last_rf2;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf2 && j<last_rf2;j++)
		{
			dpoint+=folda_rf2[j];
			ppoint+=phase_rf2[j];
		}
		if(j==last_rf2)
			break;
		ppoint=ppoint/intg_rf2;
		dpoint=dpoint/intg_rf2;
		dpoint=(dpoint-baseline_rf2)/range_rf2;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf2)
		{
			corrdata_rf2[corrpos]=dpoint;
			phaseData_rf2[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_rf2;
	}	
	
	int smoothWindow_rf1, smoothWindow_rf2;	
	float* smooth_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float* smooth_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
	smoothWindow_rf1=(int)(smoothDur_rf1/(t_rf1*intg_rf1));
	smoothWindow_rf2=(int)(smoothDur_rf2/(t_rf2*intg_rf2));
	smoothPulse(corrdata_rf1,smooth_rf1,corrsize_rf1,smoothWindow_rf1);
	smoothPulse(corrdata_rf2,smooth_rf2,corrsize_rf2,smoothWindow_rf2);
	
	if(_rf1Width==0 || _rf2Width==0)
	{
		printf("\nError! Widths cannot be zero!\n");
		return;
	}
	
	float threshFreq_rf1=1.0/_rf1Width*subPulseFreqThreshold;
	float threshFreq_rf2=1.0/_rf2Width*subPulseFreqThreshold;
	
	float maxSmooth=0,minSmooth=smooth_rf1[0],maxRes=0,minRes=0,rmsRes_rf1=0,rmsRes_rf2=0;
	
	for(i=0;i<corrsize_rf1;i++)
	{
		corrdata_rf1[i]=corrdata_rf1[i]-smooth_rf1[i];
		rmsRes_rf1+=corrdata_rf1[i]*corrdata_rf1[i];
		
		if(smooth_rf1[i]>maxSmooth)
			maxSmooth=smooth_rf1[i];

		if(smooth_rf1[i]<minSmooth)
			minSmooth=smooth_rf1[i];
		
		if(corrdata_rf1[i]>maxRes)
			maxRes=corrdata_rf1[i];
			
		if(corrdata_rf1[i]<minRes)
			minRes=corrdata_rf1[i];
	}
	for(i=0;i<corrsize_rf2;i++)
	{
		corrdata_rf2[i]=corrdata_rf2[i]-smooth_rf2[i];
		rmsRes_rf2+=corrdata_rf2[i]*corrdata_rf2[i];
		
		if(smooth_rf2[i]>maxSmooth)
			maxSmooth=smooth_rf2[i];

		if(smooth_rf2[i]<minSmooth)
			minSmooth=smooth_rf2[i];
		
		if(corrdata_rf2[i]>maxRes)
			maxRes=corrdata_rf2[i];
			
		if(corrdata_rf2[i]<minRes)
			minRes=corrdata_rf2[i];
	}
	
	for(i=0;i<last_rf1;i++)
		folda_rf1[i]=(folda_rf1[i]-baseline_rf1)/range_rf1;
	for(i=0;i<last_rf2;i++)
		folda_rf2[i]=(folda_rf2[i]-baseline_rf2)/range_rf2;
		
	float stdev_rf1=stdevBaseIntg(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1,intg_rf1),stdev_rf2=stdevBaseIntg(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2,intg_rf2);
	
	int normRes_rf1=corrsize_rf1, normRes_rf2=corrsize_rf2;
	rmsRes_rf1/=(normRes_rf1*stdev_rf1*stdev_rf1);
	rmsRes_rf2/=(normRes_rf2*stdev_rf2*stdev_rf2);
	
	rmsRes_rf1=sqrt((double)rmsRes_rf1);
	rmsRes_rf2=sqrt((double)rmsRes_rf2);
	printf("\nProcessing pulse number %d. RMS deviation of _rf1 = %f\t_rf2 = %f",pulseNum,rmsRes_rf1,rmsRes_rf2);
	
	float *autoWidths=(float*)malloc(2*sizeof(float));
	float *quasiPeriods=(float*)malloc(2*sizeof(float));
	float *relArray=(float*)malloc(4*sizeof(float));
	
	if(rmsRes_rf1>=cutOfflowfreq || rmsRes_rf2>=cutOffhighfreq)
	{
		printf("\nExceeded microstructure cut-off.");
		detectQuasiPeriod(corrdata_rf1,corrdata_rf2,corrsize_rf1,corrsize_rf2,quasiPeriods,t_rf1,t_rf2,intg_rf1,intg_rf2,threshFreq_rf1,threshFreq_rf2);
		if(quasiPeriods[0]>0)
			printf("\nQuasiperiod detected _rf1 at %f ms",quasiPeriods[0]);
		if(quasiPeriods[1]>0)
			printf("\nQuasiperiod detected _rf2 at %f ms",quasiPeriods[1]);
		
		findWidths(corrdata_rf1,corrdata_rf2,corrsize_rf1,corrsize_rf2,autoWidths,t_rf1,t_rf2,intg_rf1,intg_rf2);

		printf("\nWidths detected _rf1 = %f\t_rf2 = %f",autoWidths[0],autoWidths[1]);
		
		giveRelStrengthSmoothed(corrdata_rf1,corrdata_rf2,smooth_rf1,smooth_rf2,corrsize_rf1,corrsize_rf2,relArray);
		
		printf("\n");
			
		fprintf(microFile,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",pulseNum,rmsRes_rf1,rmsRes_rf2,autoWidths[0],autoWidths[1],quasiPeriods[0],quasiPeriods[1],relArray[0],relArray[1],relArray[2],relArray[3]);
	}
	
	free(corrdata_rf1);
	free(corrdata_rf2);
	free(phaseData_rf1);
	free(phaseData_rf2);
	free(smooth_rf1);
	free(smooth_rf2);
	free(autoWidths);
	free(quasiPeriods);
	free(relArray);
	
	return 1;
}

int microCandidates_SG(char *filename_rf1, char *filename_rf2, FILE *microFile, char *SG_coeffs, double *phase_rf1, double *phase_rf2, double s, double *folda_rf1, double *folda_rf2, int last_rf1, int last_rf2,double pstart, double pend, double pulseStart, double pulseEnd, int intg_rf1, int intg_rf2, int pulseNum, double fold, double t_rf1, double t_rf2, double lag, double smoothDur_rf1, double smoothDur_rf2, double rf1Delay, double rf2Delay,float *readData_rf1,float *readData_rf2, float cutOffhighfreq,float cutOfflowfreq, float _rf1Width, float _rf2Width, int sg_coeff_size,int option)
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
	FILE *data_rf1, *data_rf2;
	data_rf2=fopen(filename_rf2,"rb");
	data_rf1=fopen(filename_rf1,"rb");

	if(data_rf1==NULL)
	{
		printf("\nCould not open _rf1 file!");
		return 0;
	}
	if(data_rf2==NULL)
	{
		printf("\nCould not open _rf2 file!");
		return 0;
	}
	
	int numshift_rf1, numshift_rf2,j;
	int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), skipnum_rf2=(int)((pulseNum-1)*(fold/t_rf2)), sampleDelay_rf1=rf1Delay/t_rf1*1000, sampleDelay_rf2=rf2Delay/t_rf2*1000;

	numshift_rf1=(int)(last_rf1*s);
	numshift_rf2=(int)(last_rf2*s);	
	fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);
	fseek(data_rf2,(skipnum_rf2+sampleDelay_rf2+numshift_rf2)*sizeof(float),SEEK_SET);	
	
	if(readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1)!=1)
		return 0;
	if(readFloatToDouble(last_rf2,readData_rf2,folda_rf2,data_rf2)!=1)
		return 0;
	
	fclose(data_rf1);
	fclose(data_rf2);
	
	int corrsize_rf1, corrsize_rf2, corrpos=0;
	corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
	corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *corrdata_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	float *phaseData_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *phaseData_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
	double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),range_rf1=max(folda_rf1,last_rf1)-baseline_rf1,baseline_rf2=meanBase(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2),range_rf2=max(folda_rf2,last_rf2)-baseline_rf2;
	for(i=0;i<last_rf1;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
		{
			dpoint+=folda_rf1[j];
			ppoint+=phase_rf1[j];
		}
		if(j==last_rf1)
			break;
		ppoint=ppoint/intg_rf1;
		dpoint=dpoint/intg_rf1;
		dpoint=(dpoint-baseline_rf1)/range_rf1;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
		{
			corrdata_rf1[corrpos]=dpoint;
			phaseData_rf1[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_rf1;
	}
	corrpos=0;
	for(i=0;i<last_rf2;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf2 && j<last_rf2;j++)
		{
			dpoint+=folda_rf2[j];
			ppoint+=phase_rf2[j];
		}
		if(j==last_rf2)
			break;
		ppoint=ppoint/intg_rf2;
		dpoint=dpoint/intg_rf2;
		dpoint=(dpoint-baseline_rf2)/range_rf2;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf2)
		{
			corrdata_rf2[corrpos]=dpoint;
			phaseData_rf2[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_rf2;
	}	
	
	int smoothWindow_rf1, smoothWindow_rf2;	
	float* smooth_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float* smooth_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
	smoothWindow_rf1=(int)(smoothDur_rf1/(t_rf1*intg_rf1));
	smoothWindow_rf2=(int)(smoothDur_rf2/(t_rf2*intg_rf2));
	smoothPulse(corrdata_rf1,smooth_rf1,corrsize_rf1,smoothWindow_rf1);
	smoothPulse(corrdata_rf2,smooth_rf2,corrsize_rf2,smoothWindow_rf2);
	
	if(_rf1Width==0 || _rf2Width==0)
	{
		printf("\nError! Widths cannot be zero!\n");
		return;
	}
	
	float threshFreq_rf1=1.0/_rf1Width*subPulseFreqThreshold;
	float threshFreq_rf2=1.0/_rf2Width*subPulseFreqThreshold;
	
	float maxSmooth=0,minSmooth=smooth_rf1[0],maxRes=0,minRes=0,rmsRes_rf1=0,rmsRes_rf2=0,microstrPower_rf1=0,microstrPower_rf2=0;
	
	for(i=0;i<corrsize_rf1;i++)
	{
		corrdata_rf1[i]=corrdata_rf1[i]-smooth_rf1[i];
		rmsRes_rf1+=corrdata_rf1[i]*corrdata_rf1[i];
		
		if(smooth_rf1[i]>maxSmooth)
			maxSmooth=smooth_rf1[i];

		if(smooth_rf1[i]<minSmooth)
			minSmooth=smooth_rf1[i];
		
		if(corrdata_rf1[i]>maxRes)
			maxRes=corrdata_rf1[i];
			
		if(corrdata_rf1[i]<minRes)
			minRes=corrdata_rf1[i];
	}
	microstrPower_rf1=rmsRes_rf1;
	for(i=0;i<corrsize_rf2;i++)
	{
		corrdata_rf2[i]=corrdata_rf2[i]-smooth_rf2[i];
		rmsRes_rf2+=corrdata_rf2[i]*corrdata_rf2[i];
		
		if(smooth_rf2[i]>maxSmooth)
			maxSmooth=smooth_rf2[i];

		if(smooth_rf2[i]<minSmooth)
			minSmooth=smooth_rf2[i];
		
		if(corrdata_rf2[i]>maxRes)
			maxRes=corrdata_rf2[i];
			
		if(corrdata_rf2[i]<minRes)
			minRes=corrdata_rf2[i];
	}
	microstrPower_rf2=rmsRes_rf2;
	
	for(i=0;i<last_rf1;i++)
		folda_rf1[i]=(folda_rf1[i]-baseline_rf1)/range_rf1;
	for(i=0;i<last_rf2;i++)
		folda_rf2[i]=(folda_rf2[i]-baseline_rf2)/range_rf2;
		
	float stdev_rf1=stdevBaseIntg(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1,intg_rf1),stdev_rf2=stdevBaseIntg(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2,intg_rf2);
	
	int normRes_rf1=corrsize_rf1, normRes_rf2=corrsize_rf2;
	rmsRes_rf1/=(normRes_rf1*stdev_rf1*stdev_rf1);
	rmsRes_rf2/=(normRes_rf2*stdev_rf2*stdev_rf2);
	
	rmsRes_rf1=sqrt((double)rmsRes_rf1);
	rmsRes_rf2=sqrt((double)rmsRes_rf2);
	
	float maxDev_rf1=maxMod(corrdata_rf1,corrsize_rf1),maxDev_rf2=maxMod(corrdata_rf2,corrsize_rf2);
        float maxDevRMS_rf1=maxDev_rf1/stdev_rf1,maxDevRMS_rf2=maxDev_rf2/stdev_rf2;
	
	printf("\nProcessing pulse number %d. RMS deviation of _rf1 = %f\t_rf2 = %f. MaxDevRMS _rf1 = %f _rf2 = %f",pulseNum,rmsRes_rf1,rmsRes_rf2,maxDevRMS_rf1,maxDevRMS_rf2);
	
	float *autoWidths=(float*)malloc(2*sizeof(float));
	float *quasiPeriods=(float*)malloc(2*sizeof(float));
	float *relArray=(float*)malloc(4*sizeof(float));
	
	if(rmsRes_rf1>=cutOfflowfreq || rmsRes_rf2>=cutOffhighfreq || maxDevRMS_rf1>=microThreshold || maxDevRMS_rf2>=microThreshold)
	{
		printf("\nExceeded microstructure cut-off.");
		detectQuasiPeriod(corrdata_rf1,corrdata_rf2,corrsize_rf1,corrsize_rf2,quasiPeriods,t_rf1,t_rf2,intg_rf1,intg_rf2,threshFreq_rf1,threshFreq_rf2);
		if(quasiPeriods[0]>0)
			printf("\nQuasiperiod detected _rf1 at %f ms",quasiPeriods[0]);
		if(quasiPeriods[1]>0)
			printf("\nQuasiperiod detected _rf2 at %f ms",quasiPeriods[1]);
		
		findWidths_derSG(corrdata_rf1,corrdata_rf2,corrsize_rf1,corrsize_rf2,autoWidths,t_rf1,t_rf2,intg_rf1,intg_rf2,sg_coeff,sg_coeff_size);

		printf("\nWidths detected _rf1 = %f\t_rf2 = %f",autoWidths[0],autoWidths[1]);
		
		giveRelStrengthSmoothed(corrdata_rf1,corrdata_rf2,smooth_rf1,smooth_rf2,corrsize_rf1,corrsize_rf2,relArray);
		
		printf("\n");
			
		fprintf(microFile,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n",pulseNum,rmsRes_rf1,rmsRes_rf2,maxDevRMS_rf1,maxDevRMS_rf2,autoWidths[0],autoWidths[1],quasiPeriods[0],quasiPeriods[1],relArray[0],relArray[1],relArray[2],relArray[3],microstrPower_rf1,microstrPower_rf2);
	}
	
	free(corrdata_rf1);
	free(corrdata_rf2);
	free(phaseData_rf1);
	free(phaseData_rf2);
	free(smooth_rf1);
	free(smooth_rf2);
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
	FILE *data_rf1;
	data_rf1=fopen(filename_rf1,"rb");

	if(data_rf1==NULL)
	{
		printf("\nCould not open _rf1 file!");
		return 0;
	}
		
	int numshift_rf1,j;
	int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), sampleDelay_rf1=rf1Delay/t_rf1*1000;

	numshift_rf1=(int)(last_rf1*s);
	fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);	
	
	if(readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1)!=1)
		return 0;
	
	fclose(data_rf1);
	
	int corrsize_rf1, corrpos=0;
	corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
	float dpoint=0,ppoint=0;
	float *corrdata_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *phaseData_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	
	double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),range_rf1=max(folda_rf1,last_rf1)-baseline_rf1;
	for(i=0;i<last_rf1;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
		{
			dpoint+=folda_rf1[j];
			ppoint+=phase_rf1[j];
		}
		if(j==last_rf1)
			break;
		ppoint=ppoint/intg_rf1;
		dpoint=dpoint/intg_rf1;
		dpoint=(dpoint-baseline_rf1)/range_rf1;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
		{
			corrdata_rf1[corrpos]=dpoint;
			phaseData_rf1[corrpos]=ppoint;
			corrpos++;
		}
		i+=intg_rf1;
	}
	corrpos=0;
	
	int smoothWindow_rf1, smoothWindow_rf2;	
	float* smooth_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	
	smoothWindow_rf1=(int)(smoothDur_rf1/(t_rf1*intg_rf1));
	smoothPulse(corrdata_rf1,smooth_rf1,corrsize_rf1,smoothWindow_rf1);
	
	if(_rf1Width==0)
	{
		printf("\nError! Widths cannot be zero!\n");
		return;
	}
	
	float threshFreq_rf1=1.0/_rf1Width*subPulseFreqThreshold;
	
	float maxSmooth=0,minSmooth=smooth_rf1[0],maxRes=0,minRes=0,rmsRes_rf1=0,microstrPower_rf1=0;
	
	for(i=0;i<corrsize_rf1;i++)
	{
		corrdata_rf1[i]=corrdata_rf1[i]-smooth_rf1[i];
		rmsRes_rf1+=corrdata_rf1[i]*corrdata_rf1[i];
		
		if(smooth_rf1[i]>maxSmooth)
			maxSmooth=smooth_rf1[i];

		if(smooth_rf1[i]<minSmooth)
			minSmooth=smooth_rf1[i];
		
		if(corrdata_rf1[i]>maxRes)
			maxRes=corrdata_rf1[i];
			
		if(corrdata_rf1[i]<minRes)
			minRes=corrdata_rf1[i];
	}
	microstrPower_rf1=rmsRes_rf1;
	
	for(i=0;i<last_rf1;i++)
		folda_rf1[i]=(folda_rf1[i]-baseline_rf1)/range_rf1;
		
	float stdev_rf1=stdevBaseIntg(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1,intg_rf1);

	int normRes_rf1=corrsize_rf1;
	rmsRes_rf1/=(normRes_rf1*stdev_rf1*stdev_rf1);
	
	rmsRes_rf1=sqrt((double)rmsRes_rf1);
	
	float maxDev_rf1=maxMod(corrdata_rf1,corrsize_rf1);
        float maxDevRMS_rf1=maxDev_rf1/stdev_rf1;
	
	printf("\nProcessing pulse number %d. RMS deviation = %f . MaxDevRMS _rf1 = %f",pulseNum,rmsRes_rf1,maxDevRMS_rf1);
	
	float *autoWidths=(float*)malloc(2*sizeof(float));
	float *quasiPeriods=(float*)malloc(2*sizeof(float));
	float *relArray=(float*)malloc(4*sizeof(float));
	
	if(rmsRes_rf1>=cutOfflowfreq || maxDevRMS_rf1>=microThreshold )
	{
		printf("\nExceeded microstructure cut-off.");
		detectQuasiPeriod(corrdata_rf1,corrdata_rf1,corrsize_rf1,corrsize_rf1,quasiPeriods,t_rf1,t_rf1,intg_rf1,intg_rf1,threshFreq_rf1,threshFreq_rf1);
		if(quasiPeriods[0]>0)
			printf("\nQuasiperiod detected _rf1 at %f ms",quasiPeriods[0]);
		
		findWidths_derSG(corrdata_rf1,corrdata_rf1,corrsize_rf1,corrsize_rf1,autoWidths,t_rf1,t_rf1,intg_rf1,intg_rf1,sg_coeff,sg_coeff_size);

		printf("\nWidths detected = %f ",autoWidths[0]);
		
		giveRelStrengthSmoothed(corrdata_rf1,corrdata_rf1,smooth_rf1,smooth_rf1,corrsize_rf1,corrsize_rf1,relArray);
		
		printf("\n");
			
		fprintf(microFile,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",pulseNum,rmsRes_rf1,maxDevRMS_rf1,autoWidths[0],quasiPeriods[0],relArray[0],relArray[2],microstrPower_rf1);
	}
	
	free(corrdata_rf1);
	free(phaseData_rf1);
	free(smooth_rf1);
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

void createAverageSpectra(char *filename_rf1,char *filename_rf2, char *outname, char *foldFile_rf1, char *foldFile_rf2, int last_rf1, int last_rf2, double s, double *folda_rf1, double *folda_rf2 ,double pstart, double pend, double pulseStart, double pulseEnd, double *phase_rf1, double *phase_rf2, int intg_rf1, int intg_rf2, double fold, double t_rf1, double t_rf2, float *readData_rf1, float *readData_rf2, double rf1Delay, double rf2Delay)
{
	printf("\nPlease wait ... \n");
	FILE *data_rf1, *data_rf2, *out;
	data_rf1=fopen(filename_rf1,"rb");
	data_rf2=fopen(filename_rf2,"rb");
	out=fopen(outname,"w");
	int i,numshift_rf1, numshift_rf2,j;
	int corrsize_rf1, corrsize_rf2, corrpos=0;
	corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
	corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
	double dpoint=0,ppoint=0;
	float *corrdata_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *corrdata_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	float *phaseData_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *phaseData_rf2=(float*)malloc(corrsize_rf2*sizeof(float));

	numshift_rf1=(int)(last_rf1*s);
	numshift_rf2=(int)(last_rf2*s);	
	int sampleDelay_rf1=rf1Delay/t_rf1*1000, sampleDelay_rf2=rf2Delay/t_rf2*1000;
	fseek(data_rf1,(sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);
	fseek(data_rf2,(sampleDelay_rf2+numshift_rf2)*sizeof(float),SEEK_SET);
	
	double baseline_rf1,baseline_rf2;
	float maxSpec=0,minSpec=0;
	
	float stdev_rf1, stdev_rf2;
	
	float *crossArray1, *crossArray2;
	double t_corr, crossArray1_stdev, crossArray2_stdev; 
	int corrsize,flag;
	
	if(corrsize_rf1>corrsize_rf2)
	{
		crossArray1=corrdata_rf1;
		crossArray2=(float*)malloc(corrsize_rf1*sizeof(float));
		t_corr=t_rf1*intg_rf1;
		corrsize=corrsize_rf1;
		flag=1;
	}
	else
	{
		crossArray1=corrdata_rf2;
		crossArray2=(float*)malloc(corrsize_rf2*sizeof(float));	
		t_corr=t_rf2*intg_rf2;
		corrsize=corrsize_rf2;
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
		
	while(readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1)==1 && readFloatToDouble(last_rf2,readData_rf2,folda_rf2,data_rf2)==1)
	{
		baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1);
		baseline_rf2=meanBase(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2);
		
		corrpos=0;
		for(i=0;i<last_rf1;)
		{
			dpoint=0;
			ppoint=0;
			for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
			{
				dpoint+=folda_rf1[j];
				ppoint+=phase_rf1[j];
			}
			if(j==last_rf1)
				break;
			ppoint=ppoint/intg_rf1;
			dpoint=dpoint/intg_rf1;
			dpoint=(dpoint-baseline_rf1)/baseline_rf1;
		
			if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
			{
				corrdata_rf1[corrpos]=dpoint;
				phaseData_rf1[corrpos]=ppoint;	
				corrpos++;
			}
			i+=intg_rf1;
		
		}
		corrpos=0;
		for(i=0;i<last_rf2;)
		{
			dpoint=0;
			ppoint=0;
			for(j=i;j<i+intg_rf2 && j<last_rf2;j++)
			{
				dpoint+=folda_rf2[j];
				ppoint+=phase_rf2[j];
			}
			if(j==last_rf2)
				break;
			ppoint=ppoint/intg_rf2;
			dpoint=dpoint/intg_rf2;
			dpoint=(dpoint-baseline_rf2)/baseline_rf2;
		
			if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf2)
			{
				corrdata_rf2[corrpos]=dpoint;
				phaseData_rf2[corrpos]=ppoint;
				corrpos++;
			}
			i+=intg_rf2;
		}	
		for(i=0;i<last_rf1;i++)
			folda_rf1[i]=(folda_rf1[i]-baseline_rf1)/baseline_rf1;
		for(i=0;i<last_rf2;i++)
			folda_rf2[i]=(folda_rf2[i]-baseline_rf2)/baseline_rf2;
		
		stdev_rf1=stdevBaseIntg(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1,intg_rf1);
		stdev_rf2=stdevBaseIntg(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2,intg_rf2);

		if(corrsize_rf1>corrsize_rf2)
		{
			interpolate_float(corrdata_rf2,corrdata_rf1,t_rf2*intg_rf2,t_rf1*intg_rf1,crossArray2,corrsize_rf1);
			crossArray1_stdev=stdev_rf1;
			crossArray2_stdev=stdev_rf2;
		}
		else
		{
			interpolate_float(corrdata_rf1,corrdata_rf2,t_rf1*intg_rf1,t_rf2*intg_rf2,crossArray2,corrsize_rf2);
			crossArray1_stdev=stdev_rf2;
			crossArray2_stdev=stdev_rf1;
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
			fprintf(out,"%f\t%f\n",phaseData_rf1[i],foldSpectra[i]);
		else
			fprintf(out,"%f\t%f\n",phaseData_rf2[i],foldSpectra[i]);
		
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
	plottwoNoPar(foldFile_rf2,foldFile_rf1,-1);       

	cpgpanl(1,2);
	cpgenv(pstart,pend,minSpec,maxSpec,0,0);
	sprintf(title,"Average Pulse Spectrum");
	cpglab("Phase","Ratio",title);
	if(flag==1)
	        cpgline(corrsize,phaseData_rf1,foldSpectra);
        else
	        cpgline(corrsize,phaseData_rf2,foldSpectra);
	        
	fclose(data_rf1);
	fclose(data_rf2);
	fclose(out);
}

void plotSpectra(char *filename_rf1,char *filename_rf2, char *outname, int pulseNum, int last_rf1, int last_rf2, double s, double *folda_rf1, double *folda_rf2 ,double pstart, double pend, double pulseStart, double pulseEnd, double *phase_rf1, double *phase_rf2, int intg_rf1, int intg_rf2, double fold, double t_rf1, double t_rf2, float *readData_rf1, float *readData_rf2, double rf1Delay, double rf2Delay)
{
	FILE *data_rf1, *data_rf2, *out;
	data_rf1=fopen(filename_rf1,"rb");
	data_rf2=fopen(filename_rf2,"rb");
	out=fopen(outname,"w");
	int i,numshift_rf1, numshift_rf2,j;
	int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), skipnum_rf2=(int)((pulseNum-1)*(fold/t_rf2)), sampleDelay_rf1=rf1Delay/t_rf1*1000, sampleDelay_rf2=rf2Delay/t_rf2*1000;
	
	numshift_rf1=(int)(last_rf1*s);
	numshift_rf2=(int)(last_rf2*s);	
	fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);
	fseek(data_rf2,(skipnum_rf2+sampleDelay_rf2+numshift_rf2)*sizeof(float),SEEK_SET);	
	readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1);
	readFloatToDouble(last_rf2,readData_rf2,folda_rf2,data_rf2);
	
	int corrsize_rf1, corrsize_rf2, corrpos=0;
	corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
	corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
	double dpoint=0,ppoint=0;
	float *corrdata_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *corrdata_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	float *phaseData_rf1=(float*)malloc(corrsize_rf1*sizeof(float));
	float *phaseData_rf2=(float*)malloc(corrsize_rf2*sizeof(float));
	
	double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),baseline_rf2=meanBase(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2);
	float maxP_rf1=0,minP_rf1=0,maxP_rf2=0,minP_rf2=0,maxSpec=0,minSpec=0;
	
	for(i=0;i<last_rf1;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
		{
			dpoint+=folda_rf1[j];
			ppoint+=phase_rf1[j];
		}
		if(j==last_rf1)
			break;
		ppoint=ppoint/intg_rf1;
		dpoint=dpoint/intg_rf1;
		dpoint=(dpoint-baseline_rf1)/baseline_rf1;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
		{
			corrdata_rf1[corrpos]=dpoint;
			phaseData_rf1[corrpos]=ppoint;
			if(corrdata_rf1[corrpos]>maxP_rf1)
				maxP_rf1=corrdata_rf1[corrpos];
			if(corrdata_rf1[corrpos]<minP_rf1)
				minP_rf1=corrdata_rf1[corrpos];	
			corrpos++;
		}
		i+=intg_rf1;
		
	}
	corrpos=0;
	for(i=0;i<last_rf2;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf2 && j<last_rf2;j++)
		{
			dpoint+=folda_rf2[j];
			ppoint+=phase_rf2[j];
		}
		if(j==last_rf2)
			break;
		ppoint=ppoint/intg_rf2;
		dpoint=dpoint/intg_rf2;
		dpoint=(dpoint-baseline_rf2)/baseline_rf2;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf2)
		{
			corrdata_rf2[corrpos]=dpoint;
			phaseData_rf2[corrpos]=ppoint;
			if(corrdata_rf2[corrpos]>maxP_rf2)
				maxP_rf2=corrdata_rf2[corrpos];
			if(corrdata_rf2[corrpos]<minP_rf2)
				minP_rf2=corrdata_rf2[corrpos];	
			corrpos++;
		}
		i+=intg_rf2;
	}	
	for(i=0;i<last_rf1;i++)
		folda_rf1[i]=(folda_rf1[i]-baseline_rf1)/baseline_rf1;
	for(i=0;i<last_rf2;i++)
		folda_rf2[i]=(folda_rf2[i]-baseline_rf2)/baseline_rf2;
		
	float stdev_rf1=stdevBaseIntg(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1,intg_rf1),stdev_rf2=stdevBaseIntg(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2,intg_rf2);
	
	float *crossArray1, *crossArray2;
	double t_corr, crossArray1_stdev, crossArray2_stdev; 
	int corrsize,flag;
	if(corrsize_rf1>corrsize_rf2)
	{
		crossArray1=corrdata_rf1;
		crossArray2=(float*)malloc(corrsize_rf1*sizeof(float));
		interpolate_float(corrdata_rf2,corrdata_rf1,t_rf2*intg_rf2,t_rf1*intg_rf1,crossArray2,corrsize_rf1);
		t_corr=t_rf1*intg_rf1;
		corrsize=corrsize_rf1;
		crossArray1_stdev=stdev_rf1;
		crossArray2_stdev=stdev_rf2;
		flag=1;
	}
	else
	{
		crossArray1=corrdata_rf2;
		crossArray2=(float*)malloc(corrsize_rf2*sizeof(float));
		interpolate_float(corrdata_rf1,corrdata_rf2,t_rf1*intg_rf1,t_rf2*intg_rf2,crossArray2,corrsize_rf2);
		t_corr=t_rf2*intg_rf2;
		corrsize=corrsize_rf2;
		crossArray1_stdev=stdev_rf2;
		crossArray2_stdev=stdev_rf1;
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
			fprintf(out,"%f\t%f\n",phaseData_rf1[i],specArray[i]);
		else
			fprintf(out,"%f\t%f\n",phaseData_rf2[i],specArray[i]);
	}

        char title[50];	

	cpgsubp(1,3);
        cpgenv(pstart,pend,minP_rf2,maxP_rf2,0,0);
	sprintf(title,"Pulse Number %d at High Freq",pulseNum);
	cpglab("Phase","Intensity",title);
	if(flag==1)
	        cpgline(corrsize,phaseData_rf1,crossArray2);
        else
	        cpgline(corrsize,phaseData_rf2,crossArray1);        
	       
	cpgpanl(1,1);
	cpgenv(pstart,pend,minP_rf1,maxP_rf1,0,0);
	sprintf(title,"Pulse Number %d at Low Freq",pulseNum);
	cpglab("Phase","Intensity",title);
	if(flag==1)
	        cpgline(corrsize,phaseData_rf1,crossArray1);
        else
	        cpgline(corrsize,phaseData_rf2,crossArray2);        
	        
	cpgpanl(1,2);
	cpgenv(pstart,pend,minSpec,maxSpec,0,0);
	sprintf(title,"SPECTRA: Pulse Number %d",pulseNum);
	cpglab("Phase","Ratio",title);
	if(flag==1)
	        cpgline(corrsize,phaseData_rf1,specArray);
        else
	        cpgline(corrsize,phaseData_rf2,specArray);
	        
	fclose(data_rf1);
	fclose(data_rf2);
	fclose(out);
}

void crossCorrelate(char *filename_rf1,char *filename_rf2, char *outname, int pulseNum, int last_rf1, int last_rf2, double s, double *folda_rf1, double *folda_rf2 ,double pstart, double pend, double pulseStart, double pulseEnd, double *phase_rf1, double *phase_rf2, int intg_rf1, int intg_rf2, double fold, double t_rf1, double t_rf2, double lag, float *readData_rf1, float *readData_rf2, double rf1Delay, double rf2Delay)
{
	FILE *data_rf1, *data_rf2, *out;
	data_rf1=fopen(filename_rf1,"rb");
	data_rf2=fopen(filename_rf2,"rb");
	out=fopen(outname,"w");
	int i,numshift_rf1, numshift_rf2,j;
	int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), skipnum_rf2=(int)((pulseNum-1)*(fold/t_rf2)), sampleDelay_rf1=rf1Delay/t_rf1*1000, sampleDelay_rf2=rf2Delay/t_rf2*1000;
	
	/*for(i=0;i<(skipnum/last);i++)
	{
		fread(folda,sizeof(double),last,data);
	}
	fread(folda,sizeof(double),(skipnum%last),data);*/

	numshift_rf1=(int)(last_rf1*s);
	numshift_rf2=(int)(last_rf2*s);	
	fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);
	fseek(data_rf2,(skipnum_rf2+sampleDelay_rf2+numshift_rf2)*sizeof(float),SEEK_SET);	
	readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1);
	readFloatToDouble(last_rf2,readData_rf2,folda_rf2,data_rf2);
	
	int corrsize_rf1, corrsize_rf2, corrpos=0;
	corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
	corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
		
	double dpoint=0,ppoint=0;
	double *corrdata_rf1=(double*)malloc(corrsize_rf1*sizeof(double));
	double *corrdata_rf2=(double*)malloc(corrsize_rf2*sizeof(double));
	double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),range_rf1=max(folda_rf1,last_rf1)-baseline_rf1,baseline_rf2=meanBase(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2),range_rf2=max(folda_rf2,last_rf2)-baseline_rf2;
	
	for(i=0;i<last_rf1;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
		{
			dpoint+=folda_rf1[j];
			ppoint+=phase_rf1[j];
		}
		if(j==last_rf1)
			break;
		ppoint=ppoint/intg_rf1;
		dpoint=dpoint/intg_rf1;
		dpoint=(dpoint-baseline_rf1)/range_rf1;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
		{
			corrdata_rf1[corrpos++]=dpoint;
		}
		i+=intg_rf1;
	}
	corrpos=0;
	for(i=0;i<last_rf2;)
	{
		dpoint=0;
		ppoint=0;
		for(j=i;j<i+intg_rf2 && j<last_rf2;j++)
		{
			dpoint+=folda_rf2[j];
			ppoint+=phase_rf2[j];
		}
		if(j==last_rf2)
			break;
		ppoint=ppoint/intg_rf2;
		dpoint=dpoint/intg_rf2;
		dpoint=(dpoint-baseline_rf2)/range_rf2;
		
		if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf2)
		{
			corrdata_rf2[corrpos++]=dpoint;
		}
		i+=intg_rf2;
	}	
	
	double *crossArray1, *crossArray2;
	double t_corr; int corrsize;
	if(corrsize_rf1>corrsize_rf2)
	{
		crossArray1=corrdata_rf1;
		crossArray2=(double*)malloc(corrsize_rf1*sizeof(double));
		interpolate(corrdata_rf2,corrdata_rf1,t_rf2*intg_rf2,t_rf1*intg_rf1,crossArray2,corrsize_rf1);
		t_corr=t_rf1*intg_rf1;
		corrsize=corrsize_rf1;
	}
	else
	{
		crossArray1=corrdata_rf2;
		crossArray2=(double*)malloc(corrsize_rf2*sizeof(double));
		interpolate(corrdata_rf1,corrdata_rf2,t_rf1*intg_rf1,t_rf2*intg_rf2,crossArray2,corrsize_rf2);
		t_corr=t_rf2*intg_rf2;
		corrsize=corrsize_rf2;
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

  	fclose(data_rf1);
  	fclose(data_rf2);
}

void produceCrossCorrelate(char *filename_rf1,char *filename_rf2, int pulseNum, int last_rf1, int last_rf2, double s, double *folda_rf1, double *folda_rf2 ,double pstart, double pend, double pulseStart, double pulseEnd, double *phase_rf1, double *phase_rf2, int intg_rf1, int intg_rf2, double fold, double t_rf1, double t_rf2, double lag, float *readData_rf1, float *readData_rf2, double rf1Delay, double rf2Delay, double *crosscorr)
{
        FILE *data_rf1, *data_rf2, *out;
        data_rf1=fopen(filename_rf1,"rb");
        data_rf2=fopen(filename_rf2,"rb");
        int i,numshift_rf1, numshift_rf2,j;
        int skipnum_rf1=(int)((pulseNum-1)*(fold/t_rf1)), skipnum_rf2=(int)((pulseNum-1)*(fold/t_rf2)), sampleDelay_rf1=rf1Delay/t_rf1*1000, sampleDelay_rf2=rf2Delay/t_rf2*1000;

	numshift_rf1=(int)(last_rf1*s);
        numshift_rf2=(int)(last_rf2*s);
        fseek(data_rf1,(skipnum_rf1+sampleDelay_rf1+numshift_rf1)*sizeof(float),SEEK_SET);
        fseek(data_rf2,(skipnum_rf2+sampleDelay_rf2+numshift_rf2)*sizeof(float),SEEK_SET);
        readFloatToDouble(last_rf1,readData_rf1,folda_rf1,data_rf1);
        readFloatToDouble(last_rf2,readData_rf2,folda_rf2,data_rf2);

        int corrsize_rf1, corrsize_rf2, corrpos=0;
        corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
        corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));

        double dpoint=0,ppoint=0;
        double *corrdata_rf1=(double*)malloc(corrsize_rf1*sizeof(double));
        double *corrdata_rf2=(double*)malloc(corrsize_rf2*sizeof(double));
        double baseline_rf1=meanBase(folda_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1),range_rf1=max(folda_rf1,last_rf1)-baseline_rf1,baseline_rf2=meanBase(folda_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2),range_rf2=max(folda_rf2,last_rf2)-baseline_rf2;

        for(i=0;i<last_rf1;)
        {
                dpoint=0;
                ppoint=0;
                for(j=i;j<i+intg_rf1 && j<last_rf1;j++)
                {
                        dpoint+=folda_rf1[j];
                        ppoint+=phase_rf1[j];
                }
                if(j==last_rf1)
                        break;
                ppoint=ppoint/intg_rf1;
                dpoint=dpoint/intg_rf1;
                dpoint=(dpoint-baseline_rf1)/range_rf1;

                if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf1)
                {
                        corrdata_rf1[corrpos++]=dpoint;
                }
                i+=intg_rf1;
        }
        corrpos=0;
        for(i=0;i<last_rf2;)
        {
                dpoint=0;
                ppoint=0;
                for(j=i;j<i+intg_rf2 && j<last_rf2;j++)
                {
                        dpoint+=folda_rf2[j];
                        ppoint+=phase_rf2[j];
                }
                if(j==last_rf2)
                        break;
                ppoint=ppoint/intg_rf2;
                dpoint=dpoint/intg_rf2;
                dpoint=(dpoint-baseline_rf2)/range_rf2;

                if(ppoint>=pstart && ppoint<=pend && corrpos<corrsize_rf2)
                {
                        corrdata_rf2[corrpos++]=dpoint;
                }
                i+=intg_rf2;
        }

        double *crossArray1, *crossArray2;
        double t_corr; int corrsize;
        if(corrsize_rf1>corrsize_rf2)
        {
                crossArray1=corrdata_rf1;
                crossArray2=(double*)malloc(corrsize_rf1*sizeof(double));
                interpolate(corrdata_rf2,corrdata_rf1,t_rf2*intg_rf2,t_rf1*intg_rf1,crossArray2,corrsize_rf1);
                t_corr=t_rf1*intg_rf1;
                corrsize=corrsize_rf1;
        }
        else
        {
                crossArray1=corrdata_rf2;
                crossArray2=(double*)malloc(corrsize_rf2*sizeof(double));
                interpolate(corrdata_rf1,corrdata_rf2,t_rf1*intg_rf1,t_rf2*intg_rf2,crossArray2,corrsize_rf2);
                t_corr=t_rf2*intg_rf2;
                corrsize=corrsize_rf2;
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

        fclose(data_rf1);
        fclose(data_rf2);
}

void averageCrossCorrelate(char *filename_rf1,char *filename_rf2, int last_rf1, int last_rf2, double s, double *folda_rf1, double *folda_rf2 ,double pstart, double pend, double pulseStart, double pulseEnd, double *phase_rf1, double *phase_rf2, int intg_rf1, int intg_rf2, double fold, double t_rf1, double t_rf2, double lag, float *readData_rf1, float *readData_rf2, double rf1Delay, double rf2Delay, char *strongPulse, char *filecorr)
{
        int corrsize_rf1, corrsize_rf2, corrpos=0;
        corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
        corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));
        
        double t_corr; int corrsize;
        if(corrsize_rf1>corrsize_rf2)
        {
                t_corr=t_rf1*intg_rf1;
                corrsize=corrsize_rf1;
        }
        else
        {
                t_corr=t_rf2*intg_rf2;
                corrsize=corrsize_rf2;
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
		produceCrossCorrelate(filename_rf1,filename_rf2,pulseNum,last_rf1,last_rf2,s,folda_rf1,folda_rf2,pstart,pend,pulseStart,pulseEnd,phase_rf1,phase_rf2,intg_rf1,intg_rf2,  fold,t_rf1,t_rf2,lag,readData_rf1,readData_rf2,rf1Delay,rf2Delay,crosscorr);
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

void produceADP(double *ADPwindow, int numADP, double tsamp, char *outfile, int pulseNum, int panlflag)
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
}

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

int recordMicro(char *filename_rf1,char *outfile_rf1, char *filename_rf2, char *outfile_rf2, int pulseNum,int last_rf1,int last_rf2, double s,double* folda_rf1,double *folda_rf2, double pstart,double pend,double* phase_rf1,double *phase_rf2,int intg_rf1,int intg_rf2, double fold,double t_rf1,double t_rf2, double pulseStart,double pulseEnd,float* readData_rf1,float *readData_rf2,double rf1Delay,double rf2Delay, FILE* recordFile,int pulseWindow,char *foldfile_rf1,char *foldfile_rf2,double *profile_rf1,double* profile_rf2,int numshift)
{
	char *input=(char*)malloc(100*sizeof(char));
	writefile(filename_rf1,outfile_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,readData_rf1,rf1Delay);
	writefile(filename_rf2,outfile_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,readData_rf2,rf2Delay);	

	foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,2,intg_rf1,pulseStart,pulseEnd,rf1Delay);
	foldprofile(filename_rf2,foldfile_rf2,fold,t_rf2,phase_rf2,pstart,pend,numshift,profile_rf2,2,intg_rf2,pulseStart,pulseEnd,rf2Delay);
			
	cpgslct(pulseWindow);
	cpgsch(1.5);
	plottwo(outfile_rf2,outfile_rf1,pulseNum);
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

int foldedCrossCorr(char *fold_rf1,char *fold_rf2,double pstart,double pend,double lag,double t_rf1,double t_rf2,int intg_rf1,int intg_rf2,double fold,int last_rf1,int last_rf2)
{
	int i,j;
	FILE *avgprofile_rf1=fopen(fold_rf1, "r");
	FILE *avgprofile_rf2=fopen(fold_rf2, "r");
	//printf("Yo");
	if (avgprofile_rf1==NULL)
	{
		printf("No file found for _rf1");
		return 1;
	}
	if (avgprofile_rf2==NULL)
	{
		printf("No file found foe _rf2");
		return 1;
	}
	//printf("%lf %lf %lf %lf",pstart,pend,t_rf1,t_rf2);
	int corrsize_rf1, corrsize_rf2, corrpos=0;
	corrsize_rf2=(int)(last_rf2/intg_rf2*(pend-pstart));
	corrsize_rf1=(int)(last_rf1/intg_rf1*(pend-pstart));	
	//printf("%d %d",num_rf1,num_rf2);
	double *corrdata_rf1=(double*)(malloc(corrsize_rf1*sizeof(double)));
	double *phase_rf1=(double*)(malloc(corrsize_rf1*sizeof(double)));
	double *corrdata_rf2=(double*)(malloc(corrsize_rf2*sizeof(double)));
	double *phase_rf2=(double*)(malloc(corrsize_rf2*sizeof(double)));

	char singleLine[100];
	int num=0;

	while (!feof(avgprofile_rf1) && (num<=corrsize_rf1))
		{
			fgets(singleLine,100,avgprofile_rf1);		// 
			sscanf(singleLine, "%lf\t%lf", &phase_rf1[num],&corrdata_rf1[num]);	
			num++;
		}

	num=0;
	
	while (!feof(avgprofile_rf2) && (num<=corrsize_rf2))
		{
			fgets(singleLine,100,avgprofile_rf2);		// 
			sscanf(singleLine, "%lf\t%lf", &phase_rf2[num],&corrdata_rf2[num]);	
			num++;
		}

	/*for (i=0;i<corrsize_rf2;i++)
	{
		printf("%lf\n",phase_rf2[i]);
	}*/

	double *crossArray1, *crossArray2;
	double t_corr; int corrsize;
	//Finding the frequency which has higher density of data points with/without integration so that the other frequency can be brought to a similar resolution
	if(corrsize_rf1>corrsize_rf2)
	{
		crossArray1=corrdata_rf1;
		crossArray2=(double*)malloc(corrsize_rf1*sizeof(double));
		interpolate(corrdata_rf2,corrdata_rf1,t_rf2*intg_rf2,t_rf1*intg_rf1,crossArray2,corrsize_rf1);
		t_corr=t_rf1*intg_rf1;
		corrsize=corrsize_rf1;
	}
	else
	{
		crossArray2=corrdata_rf2;
		crossArray1=(double*)malloc(corrsize_rf2*sizeof(double));
		interpolate(corrdata_rf1,corrdata_rf2,t_rf1*intg_rf1,t_rf2*intg_rf2,crossArray1,corrsize_rf2);
		t_corr=t_rf2*intg_rf2;
		corrsize=corrsize_rf2;
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


void plotDrifts(char *filename, double *folda,double *phase, double pulseStart,double pulseEnd,float *readData,int nbins,int last, double s,double fold, double t, char *foldfile, double *profile)
{	
	FILE *data;
	data=fopen(filename,"rb");
	int pulseNum=1,skipnum=0,NumofPulses,i,j,k;
	double tempPhase;
	printf("\nThis may take some time...\n");
	while(1)
	{
		skipnum=(int)((pulseNum-1)*(fold/t));
		fseek(data,skipnum*sizeof(float),SEEK_SET);

		if(readFloatToDouble(last,readData,folda,data)!=1)
			break;
		pulseNum++;
	}
	NumofPulses = pulseNum-1;
	
	double stepsize=(pulseEnd-pulseStart)/nbins;

	if (NumofPulses>5000)
		NumofPulses=5000;
/*
	float **Int2D;
	Int2D = malloc(NumofPulses * sizeof(float *));
	if(Int2D == NULL)
	{
		printf("out of memory");
		return;
	}
	for(i=0;i<NumofPulses; i++)
	{
		Int2D[i] = malloc(nbins * sizeof(float *));
		if(Int2D[i]==NULL)
		{
			printf("out of memory");
			return;
		}
	}
*/
	float **Int2D_2;
	Int2D_2 = malloc(NumofPulses * sizeof(float *));
	if(Int2D_2 == NULL)
	{
		printf("out of memory");
		return;
	}
	for(i=0;i<NumofPulses; i++)
	{
		Int2D_2[i] = malloc(nbins * sizeof(float *));
		if(Int2D_2[i]==NULL)
		{
			printf("out of memory");
			return;
		}
	}


	float Int2D[NumofPulses][nbins];
	//float Int2D_2[NumofPulses][nbins];

	for (i=0;i<NumofPulses;i++)
	{
		for(j=0;j<nbins;j++)
		{	
			Int2D[i][j]=0;	
			Int2D_2[i][j]=0;		
		}
			
	}
	
	float trans[6];
	float count;
	for (i=0;i<6;i++)
	{	if (i==1 || i==5)
			trans[i]=1;
		else
			trans[i]=0;
	}
	
	for (i=0;i<NumofPulses;i++)
	{
		for(j=0;j<nbins;j++)
		{	skipnum=(int)((i)*(fold/t));
			fseek(data,skipnum*sizeof(float),SEEK_SET);
	
			if(readFloatToDouble(last,readData,folda,data)!=1)
				break;
			count=0;
			for(k=0;k<last;k++)
			{	tempPhase = phase[k];
				if(tempPhase>pulseStart+j*stepsize && tempPhase<pulseStart + (j+1)*stepsize)
				{	
					Int2D[i][j]+=(float)folda[k];
					count++;
				}		
			}
			Int2D[i][j]/=count;
		}
	}

	int collapse,pulseNumStart,pulseNumEnd;
	char option='n';
	char foldopt='f';
	int pulseWindow=cpgopen("/xs");
	cpgslct(pulseWindow);
	cpgpap(4.75,1.5);
	int NumofPulses_buff=NumofPulses;
	float Int2D_buff[NumofPulses][nbins];
	
	for (i=0;i<NumofPulses;i++)
	{
		for(j=0;j<nbins;j++)
		{	
			Int2D_buff[i][j]=Int2D[i][j];		
		}
			
	}
	
	while(option=='n')
	{	NumofPulses=NumofPulses_buff;
		for (i=0;i<NumofPulses;i++)
		{
			for(j=0;j<nbins;j++)
			{	
				Int2D[i][j]=Int2D_buff[i][j];
				Int2D_2[i][j]=0;		
			}
			
		}
		printf("\nEnter number of pulses to be collapsed:\t");
		scanf("%d",&collapse);
	
		if(collapse>NumofPulses)
		{	
			printf("\nError! Not enough pulses to collapse.");
			return;	
		}	

		if(collapse!=1)
		{
			int rem = NumofPulses%collapse;

			for(i=0;i<NumofPulses-rem;i+=collapse)
			{
				for(j=0;j<nbins;j++)
				{	
					for(k=0;k<collapse;k++)
					{
						Int2D_2[i][j]+=Int2D[i+k][j];
					}

						Int2D_2[i][j]/=collapse;
			
					for(k=0;k<collapse;k++)
					{
						Int2D[i+k][j]=Int2D_2[i][j];
					}
				}
			}
		}

		float max1=0,min1=Int2D[0][0];
		for (i=0;i<NumofPulses;i++)
		{
			for(j=0;j<nbins;j++)
			{	
				if(Int2D[i][j]>max1)
					max1 = Int2D[i][j];
				if(Int2D[i][j]<min1)
					min1 = Int2D[i][j];
			}	
		}
		printf("Max = %f\tMin = %f",max1,min1);		
		
		int plotNum;		
		printf("\nHow many pulses do you wish to examine?\n(Press 0 to examine all pulses)\t");
		scanf("%d",&plotNum);
		
		if(plotNum>NumofPulses)
		{		
			printf("\nError! Not enough pulses!");
			return;
		}
		else if (plotNum==0)
		{
			pulseNumStart=1;
			printf("\nGenerating greyscale plot for all (%d) pulses...\n",NumofPulses);
			pulseNumEnd=NumofPulses;
		}
		else 
		{
			NumofPulses=plotNum;
			
			printf("Enter starting pulse number (1 or above)");
			scanf("%d",&pulseNumStart);		
			printf("\nGenerating greyscale plot for %d pulses starting from pulse number %d\n",NumofPulses,pulseNumStart);
			pulseNumEnd=pulseNumStart+NumofPulses-1;
		}
		
	X1:
		cpgsubp(1,1);
		cpgenv(1,nbins,pulseNumStart,pulseNumEnd,0,0);
		cpglab("Bins","Pulse Number","On-pulse energy");

		cpggray((float*)Int2D,nbins,NumofPulses_buff,1,nbins,pulseNumStart,pulseNumEnd,min1+(max1-min1)*0.03,max1-(max1-min1)*0.35,(float*)trans);
		
		printf("\nPress 'f' to see folded profile. Otherwise, press any key.\t");
		scanf(" %c",&foldopt);

		if(foldopt=='f')	
		{
		cpgsubp(1,2);
			foldprofile(filename,foldfile,fold,t,phase,pulseStart,pulseEnd,0,profile,2,1,pulseStart,pulseEnd,0);	
			plotone_single(foldfile,-1,2);
		cpgpanl(1,1);
			cpgenv(1,nbins,pulseNumStart,pulseNumEnd,0,0);
			cpglab("Bins","Pulse Number","On-pulse energy");
			cpggray((float*)Int2D,nbins,NumofPulses_buff,1,nbins,pulseNumStart,pulseNumEnd,min1+(max1-min1)*0.03,max1-(max1-min1)*0.35,(float*)trans);
		}
	X2:
		printf("\nTo continue analysis, press 'n'\nTo go back to previous plot, press 'p'\nTo quit this routine, press 'q'\t");
		scanf(" %c",&option);
		if (option=='n')
			continue;
		else if(option=='p')
			goto X1;
		else if(option=='q')
		{	
			printf("\nExiting routine...");
			cpgpap(8.75,0.8);			
			return;
		}
		else 
		{		
			printf("\nInvalid option. Try again");
			goto X2;
		}
	}
	
}




int main(int argc, char* argv[])
{

  int option;
  while(1)
  {
	  bool rerun=false;	
	
	  printf("\nSingle frequency=1 or dual frequency=2\n");
	  scanf("%d",&option);
		switch(option)
		{
		case 2:
		{
		double t_rf1,t_rf2,fold,value,time,k,lag=0,DM,smoothDur_rf1=0,smoothDur_rf2=0;
		int last_rf1,last_rf2,i,sample=0,count=0,scanPulseNum=0;
		int b=0;
		char *input=(char*)malloc(100*sizeof(char));

//-----------------------------------------------------------

		int in;
		char rf1file[1000],rf2file[10000],file[1000], timerf1[1000],timerf2[1000];
		float reso,rf1reso,rf2reso,pulseper,dm,timeOffset;
		FILE *parameter = fopen("/Data/ygupta/psr_microstr_work/recovered_code/dual_freq_analysis/parameters.pulse_analysis", "r");
		
		if (parameter == NULL) 
		{
		perror("Error: Failed to open file.");
		return 1;
		}

		char *parameters[50],line[600], pulsar[15];
		for(i=0;i<50;i++) parameters[i]=(char*) malloc(sizeof(char)*600);
		int num_pulsars=0, num=0;	

		for(i=0;fgets(line,sizeof(line),parameter);i++)	
		{
			strcpy(parameters[i],line);
			sscanf(parameters[i], "%d)%s %*s", &num,pulsar);
			printf("%d) %s\n",num,pulsar);
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
		sscanf(parameters[pulsarnum-1], "%*d)%s %f %f %f %f %f %s %s %s %s", pulsarchosen,&timeOffset,&rf1reso,&rf2reso,&pulseper,&dm,rf1file,rf2file,timerf1,timerf2);
		
		//printf("\n_rf1 file is %s\n",_rf1file);
		printf("\nrf1 resolution is %f",rf1reso);	
		printf("\nrf2 resolution %f",rf2reso);
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
		
		t_rf1=rf1reso;
		t_rf2=rf2reso;
		fold=pulseper;
		DM=dm;

		char *filename_rf1=rf1file;
		char *filename_rf2=rf2file;
		char *outfile_rf1="plot_rf1.txt";
		char *outfile_rf2="plot_rf2.txt";
		char *foldfile_rf1="fold_rf1.txt";
		char *foldfile_rf2="fold_rf2.txt";
		char *outcorr_rf1="corr_rf1.txt";
		char *outcorr_rf2="corr_rf2.txt";
		char *avgcorr_rf1="avgcorr_rf1.txt";
		char *avgcorr_rf2="avgcorr_rf2.txt";
		char *strong_rf1="strongpulse_rf1.txt";
		char *strong_rf2="strongpulse_rf2.txt";
		char *adp_rf1="ADP_rf1.txt";
		char *adp_rf2="ADP_rf2.txt";
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

		/*double t_rf1,t_rf2,fold,value,time,k,lag=0,DM,smoothDur_rf1=0,smoothDur_rf2=0;
		int last_rf1,last_rf2,i,sample=0,count=0,scanPulseNum=0;
		int b=0;
		char *input=(char*)malloc(100*sizeof(char));*/

		/*t_rf1=(double)atof(argv[3]);
		t_rf2=(double)atof(argv[4]);
		fold=(double)atof(argv[5]);
		DM=(double)atof(argv[6]);*/
		
		if (floor((double)fold/t_rf1)==fold/t_rf1)
			last_rf1=floor(fold/t_rf1);
		else
			last_rf1=floor(fold/t_rf1)+1;
		
		if (floor((double)fold/t_rf2)==fold/t_rf2)
		        last_rf2=floor(fold/t_rf2);
		else
		        last_rf2=floor(fold/t_rf2)+1;

		float cutOffhighfreq=0,cutOfflowfreq=0;
		char tempStore[100];
		int rf2_hour, rf2_min, rf1_hour, rf1_min, temp;
		double rf2_sec, rf1_sec,_rf1_fracSec, _rf1_fullSec;

		FILE *rf2_Stamp=fopen(timerf2,"r");

		if(rf2_Stamp==NULL)
		{
			printf("\nrf2 timestamp not found! Exit.\n");
			exit(-1);
		}

		fscanf(rf2_Stamp,"%s %s %s %s\n",tempStore,tempStore,tempStore,tempStore);
		fscanf(rf2_Stamp,"IST Time: %d:%d:%lf",&rf2_hour,&rf2_min,&rf2_sec);

		FILE *rf1_Stamp=fopen(timerf1,"r");
		if(rf1_Stamp==NULL)
		{
			printf("\nrf1 timestamp not found! Exit.\n");
			exit(-1);
		}
		
		fscanf(rf1_Stamp,"%s %s %s %s\n",tempStore,tempStore,tempStore,tempStore);
		fscanf(rf1_Stamp,"IST Time: %d:%d:%lf",&rf1_hour,&rf1_min,&rf1_sec);
		
		//double relDiff=diffTime(rf2_hour,rf2_min,rf2_sec,rf1_hour,rf1_min,rf1_sec);
		double relDiff=0;
		//double relDelay=calculateRelDelay(relDiff,DM,timeOffset);
		double relDelay=0;
		printf("\nrf2 timestamp is %02d:%02d:%2.9lf.\nrf1 timestamp is %02d:%02d:%2.9lf.\n",rf2_hour,rf2_min,rf2_sec,rf1_hour,rf1_min,rf1_sec);
		printf("\nTime difference between timestamps is %lf secs (rf1 ahead of rf2)\nRelative delay to be adjusted is %lf s.\n",relDiff,relDelay);

		double *folda_rf1=(double*)malloc(last_rf1*sizeof(double));
		double *folda_rf2=(double*)malloc(last_rf2*sizeof(double));
		float *readData_rf1=(float*)malloc(last_rf1*sizeof(float));
		float *readData_rf2=(float*)malloc(last_rf2*sizeof(float));
		double *phase_rf1=(double*)malloc(last_rf1*sizeof(double));
		double *phase_rf2=(double*)malloc(last_rf2*sizeof(double));
		double *profile_rf1=(double*)malloc(last_rf1*sizeof(double));	
		double *profile_rf2=(double*)malloc(last_rf2*sizeof(double));
		
		//printf("\nPulse Viewer\nInstructions:\nUse 'n' to go to next pulse\nUse 'b' to go the earlier pulse\nUse 'p <start phase> <end phase>' to set a phase range (Default is full period)\nUse 's <phase shift>' to shift phase window (Phase shift should remain between 0 and 1 !)\nUse 'f' to show folded profile\nUse 'h' to hide folded profile\nUse 'i <Integration>' to integrate profile\nUse 'g <Pulse Number>' to get profile for a specified pulse number\nUse 'a <Time lag>' to get autocorrelation plot for the current pulse, integration and phase range\nUse 'c <MinSNR>' to get list of pulse numbers above a given threshold\nUse 't <MinSNR>' to get list of pulse numbers with peak above a given SNR\nUse 'v <Time lag(ms)> to get the average ACF for the strongest pulses\nUse 'y' to view the ADP for the specified window\nUse 'k' to get instructions\nUse 'q' to terminate program\n");
		int pulseNum=0,numshift=0,plotswitch=0,intg_rf1=1,intg_rf2=1,curshift_rf1=0,curshift_rf2=0,backOp,visualPulseindex=0,val=0;
		double s=0,pstart=0,pend=1,pulseStart=0,pulseEnd=1,minSNR=0;
		char *pch;
		struct  visualPulse visarr[1000];
		int autoflag=0, autoWindow, psWindow, foldWindow;	
		int sg_coeff_size=0;
		int pulseWindow=cpgopen("/xs");
		cpgask(0);
		
		double width_rf1, width_rf2;

		double rf1Delay, rf2Delay;
		if(relDelay>=0)
		{
			rf1Delay=0;
			rf2Delay=relDelay;
		}
		else
		{
			rf1Delay=-relDelay;
			rf2Delay=0;
		}
		
		for(i=0;i<last_rf1;i++)
		{
			phase_rf1[i]=(double)i/last_rf1;
		}
		for(i=0;i<last_rf2;i++)
		{
		        phase_rf2[i]=(double)i/last_rf2;
		}
		
		foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,1,intg_rf1,pulseStart,pulseEnd,rf1Delay);
		foldprofile(filename_rf2,foldfile_rf2,fold,t_rf2,phase_rf2,pstart,pend,numshift,profile_rf2,1,intg_rf2,pulseStart,pulseEnd,rf2Delay);
		
		printf("\nrf1 default integration is %d samples. rf2 default integration is %d samples\n",intg_rf1,intg_rf2);
		
		printf("\nPlease check visually and specify pulse on window for baseline and noise statistics. Showing the two folded profiles.\n");

		cpgslct(pulseWindow);
		cpgsch(1.5);
		plottwo(foldfile_rf2,foldfile_rf1,-1);
		
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
				

				writefile(filename_rf1,outfile_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,readData_rf1,rf1Delay);
				writefile(filename_rf2,outfile_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,readData_rf2,rf2Delay);	

				foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,2,intg_rf1,pulseStart,pulseEnd,rf1Delay);
				foldprofile(filename_rf2,foldfile_rf2,fold,t_rf2,phase_rf2,pstart,pend,numshift,profile_rf2,2,intg_rf2,pulseStart,pulseEnd,rf2Delay);
				
				cpgslct(pulseWindow);
				if(plotswitch==0)		
				{
					cpgsch(1.5);
					plottwo(outfile_rf2,outfile_rf1,pulseNum);
				}
				else
				{
					cpgsch(2);
					plotfour(outfile_rf2,foldfile_rf2,outfile_rf1,foldfile_rf1,pulseNum);
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
		                writefile(filename_rf1,outfile_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,readData_rf1,rf1Delay);
		                writefile(filename_rf2,outfile_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,readData_rf2,rf2Delay);      
		                        
				foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,2,intg_rf1,pulseStart,pulseEnd,rf1Delay);
		                foldprofile(filename_rf2,foldfile_rf2,fold,t_rf2,phase_rf2,pstart,pend,numshift,profile_rf2,2,intg_rf2,pulseStart,pulseEnd,rf2Delay);

		                cpgslct(pulseWindow);
				if(plotswitch==0)		
				{
					cpgsch(1.5);
					plottwo(outfile_rf2,outfile_rf1,pulseNum);
				}
				else
				{
					cpgsch(2);
					plotfour(outfile_rf2,foldfile_rf2,outfile_rf1,foldfile_rf1,pulseNum);
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
			
		                writefile(filename_rf1,outfile_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,readData_rf1,rf1Delay);
		                writefile(filename_rf2,outfile_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,readData_rf2,rf2Delay);      
		                foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,2,intg_rf1,pulseStart,pulseEnd,rf1Delay);
		                foldprofile(filename_rf2,foldfile_rf2,fold,t_rf2,phase_rf2,pstart,pend,numshift,profile_rf2,2,intg_rf2,pulseStart,pulseEnd,rf2Delay);

		                cpgslct(pulseWindow);
				if(plotswitch==0)		
				{
					cpgsch(1.5);
					plottwo(outfile_rf2,outfile_rf1,pulseNum);
				}
				else
				{
					cpgsch(2);
					plotfour(outfile_rf2,foldfile_rf2,outfile_rf1,foldfile_rf1,pulseNum);
				}
			}
			else if(input[0]=='q')
			{
				printf("\nTerminating program!\n");
				cpgend();
				return 0;
			}
			else if(input[0]=='s')
			{
				pch=strtok(input," ");
				pch=strtok(NULL," ");
				curshift_rf1=(int)(atof(pch)*last_rf1);
				curshift_rf2=(int)(atof(pch)*last_rf2);
				s+=atof(pch);
				
				pulseStart=pulseStart-s;
				if(pulseStart<0)
					pulseStart=1+pulseStart;
				pulseEnd=pulseEnd-s;
				if(pulseEnd<0)
					pulseEnd=1+pulseEnd;
					
				printf("\nNew shift is %lf.New pulseStart is %lf. New pulseEnd is %lf\n",s,pulseStart,pulseEnd);
				
writefile(filename_rf1,outfile_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,readData_rf1,rf1Delay);
writefile(filename_rf2,outfile_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,readData_rf2,rf2Delay);		
		foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,curshift_rf1,profile_rf1,3,intg_rf1,pulseStart,pulseEnd,rf1Delay);
				foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,2,intg_rf1,pulseStart,pulseEnd,rf1Delay);
				

				cpgslct(pulseWindow);
				if(plotswitch==0)		
				{
					cpgsch(1.5);
					plottwo(outfile_rf2,outfile_rf1,pulseNum);
				}
				else
				{
					cpgsch(2);
					plotfour(outfile_rf2,foldfile_rf2,outfile_rf1,foldfile_rf1,pulseNum);
				}


			}
			else if(input[0]=='f')
			{
				plotswitch=1;
	writefile(filename_rf1,outfile_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,readData_rf1,rf1Delay);
		                writefile(filename_rf2,outfile_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,readData_rf2,rf2Delay);      
		                        
				foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,2,intg_rf1,pulseStart,pulseEnd,rf1Delay);
		                foldprofile(filename_rf2,foldfile_rf2,fold,t_rf2,phase_rf2,pstart,pend,numshift,profile_rf2,2,intg_rf2,pulseStart,pulseEnd,rf2Delay);

		                cpgslct(pulseWindow);
				cpgsch(2);
				plotfour(outfile_rf2,foldfile_rf2,outfile_rf1,foldfile_rf1,pulseNum);

				if (input[2] == 's')
				{

					int i,last_rf1,last_rf2,count=0;
					if (floor((double)fold/t_rf1)==fold/t_rf1)
					{
						last_rf1=floor(fold/t_rf1);
					}
					else
					{
						last_rf1=floor(fold/t_rf1)+1;
					}

					if (floor((double)fold/t_rf2)==fold/t_rf2)
					{
						last_rf2=floor(fold/t_rf2);
					}
					else
					{
						last_rf2=floor(fold/t_rf2)+1;
					}


					
					double baseline_rf1=0,baseline_rf2=0,onpulse_mean_rf1=0,onpulse_mean_rf2=0;
					double stdev_off_rf1=0,stdev_off_rf2=0,tempPhase;
					double SNR_rf1,SNR_rf2;

					baseline_rf2=meanBase(profile_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2);
					stdev_off_rf2=stdevBase(profile_rf2,phase_rf2,pulseStart,pulseEnd,last_rf2);

					for(i=0;i<last_rf2;i++)
					{
						tempPhase=phase_rf2[i];
						if(tempPhase>=pulseStart && tempPhase<=pulseEnd)
						{
							onpulse_mean_rf2+=profile_rf2[i]-baseline_rf2;
							count++;
						}			
					}
					SNR_rf2=(onpulse_mean_rf2)/((stdev_off_rf2)*sqrt(count));
					printf("\nBest estimate for SNR of folded profile _rf2:%lf",SNR_rf2);

					count=0;

					baseline_rf1=meanBase(profile_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1);
					stdev_off_rf1=stdevBase(profile_rf1,phase_rf1,pulseStart,pulseEnd,last_rf1);

					for(i=0;i<last_rf1;i++)
					{
						tempPhase=phase_rf1[i];
						if(tempPhase>=pulseStart && tempPhase<=pulseEnd)
						{
							onpulse_mean_rf1+=profile_rf1[i]-baseline_rf1;
							count++;
						}			
					}
					SNR_rf1=(onpulse_mean_rf1)/((stdev_off_rf1)*sqrt(count));
					printf("\nBest estimate for SNR of folded profile _rf1:%lf\n",SNR_rf1);

					double ratio;
					ratio=SNR_rf1/SNR_rf2;
					printf("SNR of _rf1 / SNR of _rf2 (folded profile) = %lf\n",ratio);

				}
					
					
			}
			else if(input[0]=='h')
			{
				plotswitch=0;
	writefile(filename_rf1,outfile_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,readData_rf1,rf1Delay);
		                writefile(filename_rf2,outfile_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,readData_rf2,rf2Delay);      
		                        
				foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,2,intg_rf1,pulseStart,pulseEnd,rf1Delay);
		                foldprofile(filename_rf2,foldfile_rf2,fold,t_rf2,phase_rf2,pstart,pend,numshift,profile_rf2,2,intg_rf2,pulseStart,pulseEnd,rf2Delay);

		                cpgslct(pulseWindow);
				cpgsch(1.5);
				plottwo(outfile_rf2,outfile_rf1,pulseNum);
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
					intg_rf2=atoi(pch);
					printf("\nNew integration is %d samples for rf2\n", intg_rf2);
				}
				else if(backOp==2)
				{
					intg_rf1=atoi(pch);
					printf("\nNew integration is %d samples for rf1\n", intg_rf1);
				}
				else
				{
					printf("\nInvalid backend option!\n");
					continue;
				}
				
				writefile(filename_rf1,outfile_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,readData_rf1,rf1Delay);
		                writefile(filename_rf2,outfile_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,readData_rf2,rf2Delay);      
		                foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,2,intg_rf1,pulseStart,pulseEnd,rf1Delay);
		                foldprofile(filename_rf2,foldfile_rf2,fold,t_rf2,phase_rf2,pstart,pend,numshift,profile_rf2,2,intg_rf2,pulseStart,pulseEnd,rf2Delay);

		                cpgslct(pulseWindow);
				if(plotswitch==0)		
				{
					cpgsch(1.5);
					plottwo(outfile_rf2,outfile_rf1,pulseNum);
				}
				else
				{
					cpgsch(2);
					plotfour(outfile_rf2,foldfile_rf2,outfile_rf1,foldfile_rf1,pulseNum);
				}
			}
			else if(input[0]=='k')
			{
				printf("\nInstructions:\nUse 'n' to go to next pulse\nUse 'b' to go the earlier pulse\nUse 'p <start phase> <end phase>' to set a phase range (Default is full period)\nUse 's <phase shift>' to shift phase window (Phase shift should remain between 0 and 1 !)\nUse 'f' to show folded profile\nUse 'h' to hide folded profile\nUse 'i <RF number 1/2> <Integration>' to integrate profile\nUse 'g <Pulse Number>' to get profile for a specified pulse number\nUse 'a <Time lag>' to get autocorrelation plot for the current pulse, integration and phase range\nUse ‘m <min SNR>’ to output list of pulses above some average SNR \nUse 't <MinSNR>' to get list of pulse numbers with peak above a given SNR\nUse 'v <Time lag(ms)> to get the average ACF for the strongest pulses\nUse 'y' to view the ADP for the specified window\nUse 'k' to get instructions\nUse 'q' to terminate program\nUse 'e <time lag>' to get average cross correlation\n\nSUBTRACTION FEATURE:\n\nUse 'd <smoothing (ms)> <microCuttOff(SNR)> <SG order>’ to run Analysis of Current Pulse, that smoothes the pulse, subtracts the smoothed version, and then computes widths from the subtracted ACF smoothed with a Savitzky Golay filter and computes periods from the power spectrum of the subtracted features\n\n");			
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
					correlate(filename_rf2,outcorr_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,lag,pulseStart,pulseEnd,readData_rf2,rf2Delay,0);
					cpgpanl(1,1);
					correlate(filename_rf1,outcorr_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,lag,pulseStart,pulseEnd,readData_rf1,rf1Delay,1);
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
					correlate(filename_rf2,outcorr_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,lag,pulseStart,pulseEnd,readData_rf2,rf2Delay,0);
					cpgpanl(1,1);
					correlate(filename_rf1,outcorr_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,lag,pulseStart,pulseEnd,readData_rf1,rf1Delay,1);
					cpgpanl(1,2);
					crossCorrelate(filename_rf1,filename_rf2,crossCorr,pulseNum,last_rf1,last_rf2,s,folda_rf1,folda_rf2,pstart,pend,pulseStart,pulseEnd,phase_rf1,phase_rf2,intg_rf1,intg_rf2,fold,t_rf1,t_rf2,lag,readData_rf1,readData_rf2,rf1Delay,rf2Delay);
					
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

				writefile(filename_rf1,outfile_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,readData_rf1,rf1Delay);
		                writefile(filename_rf2,outfile_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,readData_rf2,rf2Delay);      
		                foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,2,intg_rf1,pulseStart,pulseEnd,rf1Delay);
		                foldprofile(filename_rf2,foldfile_rf2,fold,t_rf2,phase_rf2,pstart,pend,numshift,profile_rf2,2,intg_rf2,pulseStart,pulseEnd,rf2Delay);

		                cpgslct(pulseWindow);
		                if(plotswitch==0)		
				{
					cpgsch(1.5);
					plottwo(outfile_rf2,outfile_rf1,pulseNum);
				}
				else
				{
					cpgsch(2);
					plotfour(outfile_rf2,foldfile_rf2,outfile_rf1,foldfile_rf1,pulseNum);
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
					printf("\nStrong pulses for _rf2:\n");
					showStrongPulses(filename_rf2,last_rf2,s,folda_rf2,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,minSNR,strong_rf2,readData_rf2,rf2Delay);
					printf("\nStrong pulses for _rf1:\n");
					showStrongPulses(filename_rf1,last_rf1,s,folda_rf1,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,minSNR,strong_rf1,readData_rf1,rf1Delay);				
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
					printf("\nStrong PEAK pulses for _rf2:\n");
					showStrongPulsesPeak(filename_rf2,last_rf2,s,folda_rf2,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,minSNR,strong_rf2,readData_rf2,rf2Delay);
					printf("\nStrong PEAK pulses for _rf1:\n");
					showStrongPulsesPeak(filename_rf1,last_rf1,s,folda_rf1,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,minSNR,strong_rf1,readData_rf1,rf1Delay);				
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
					averageAutoCorrelate(filename_rf2,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,lag,pulseStart,pulseEnd,strong_rf2,avgcorr_rf2,readData_rf2,rf2Delay,0);
					cpgpanl(1,1);
					averageAutoCorrelate(filename_rf1,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,lag,pulseStart,pulseEnd,strong_rf1,avgcorr_rf1,readData_rf1,rf1Delay,1);
					
					width_rf1=getAvgHalfWidth(avgcorr_rf1,t_rf1,lag,intg_rf1);
					width_rf2=getAvgHalfWidth(avgcorr_rf2,t_rf2,lag,intg_rf2);
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
					if(width_rf1==0 || width_rf2==0)
					{
						printf("\nERROR: AvgWidth is zero!\n");
						continue;
					}
					smoothDur_rf1=width_rf1;
					smoothDur_rf2=width_rf2;
					printf("\nTaking default ACF half-width _rf1 = %f; ACF half-width _rf2 = %f\n",smoothDur_rf1,smoothDur_rf2);
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
				{	smoothDur_rf1=atof(pch);
					if(smoothDur_rf1<0)
					{
						printf("\nInvalid smoothing duration!\n");
						continue;
					}	
					smoothDur_rf2=smoothDur_rf1;
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
				printf("\nPlotting pulse smoothed to _rf1: %f ms _rf2: %f ms, and autocorrelation upto lag %f ms.\n",smoothDur_rf1,smoothDur_rf2,lag);
				
				cpgsch(2.5);
				smoothSubACF(filename_rf1,filename_rf2,subAutoFile,phase_rf1,phase_rf2,s,folda_rf1,folda_rf2,last_rf1,last_rf2,pstart,pend,pulseStart,pulseEnd,intg_rf1,intg_rf2,pulseNum,fold,t_rf1,t_rf2,lag,smoothDur_rf1,smoothDur_rf2,rf1Delay,rf2Delay,readData_rf1,readData_rf2);
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
					if(width_rf1==0 || width_rf2==0)
					{
						printf("\nERROR: AvgWidth is zero!\n");
						continue;
					}
					smoothDur_rf1=width_rf1;
					smoothDur_rf2=width_rf2;
					printf("\nTaking default ACF half-width _rf1 = %f; ACF half-width _rf2 = %f\n",smoothDur_rf1,smoothDur_rf2);
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
				{	smoothDur_rf1=atof(pch);
					if(smoothDur_rf1<0)
					{
						printf("\nInvalid smoothing duration!\n");
						continue;
					}	
					smoothDur_rf2=smoothDur_rf1;
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
				printf("\nPlotting pulse smoothed to _rf1: %f ms _rf2: %f ms, and autocorrelation upto lag %f ms.\n",smoothDur_rf1,smoothDur_rf2,lag);
				
				cpgsch(2.5);
				smoothSubSGACF(filename_rf1,filename_rf2,subAutoFile,SG_coeffs,phase_rf1,phase_rf2,s,folda_rf1,folda_rf2,last_rf1,last_rf2,pstart,pend,pulseStart,pulseEnd,intg_rf1,intg_rf2,pulseNum,fold,t_rf1,t_rf2,lag,smoothDur_rf1,smoothDur_rf2,rf1Delay,rf2Delay,readData_rf1,readData_rf2,sg_coeff_size,option,resicrosscorr);
			}
			else if(input[0]==',')
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
					if(width_rf1==0 || width_rf2==0)
					{
						printf("\nERROR: AvgWidth is zero!\n");
						continue;
					}
					smoothDur_rf1=width_rf1;
					smoothDur_rf2=width_rf2;
					printf("\nTaking default ACF half-width _rf1 = %f; ACF half-width _rf2 = %f\n",smoothDur_rf1,smoothDur_rf2);
				}
				/*else if(pch[0]=='s')
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
				{	smoothDur_rf1=atof(pch);
					if(smoothDur_rf1<0)
					{
						printf("\nInvalid smoothing duration!\n");
						continue;
					}	
					smoothDur_rf2=smoothDur_rf1;
				}
				
				if(autoflag==0)
				{
					autoWindow=cpgopen("/xs");
					cpgask(0);
					autoflag=1;
				}
				
				cpgslct(autoWindow);
				
				if(smoothDur_rf1==-1 || smoothDur_rf2==-1)
					printf("\nPlotting pulse smoothed to ACF half-width and autocorrelation upto lag %f ms.\n",lag);
				else
				printf("\nPlotting pulse smoothed to _rf1: %f ms _rf2: %f ms, and power spectrum of residuals.\n",smoothDur_rf1,smoothDur_rf2);
				
				cpgsch(2.5);
				smoothSubPowerSpectrum(filename_rf1,filename_rf2,FFTName,phase_rf1,phase_rf2,s,folda_rf1,folda_rf2,last_rf1,last_rf2,pstart,pend,pulseStart,pulseEnd,intg_rf1,intg_rf2,pulseNum,fold,t_rf1,t_rf2,smoothDur_rf1,smoothDur_rf2,rf1Delay,rf2Delay,readData_rf1,readData_rf2,width_rf1,width_rf2,option);
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
					if(width_rf1==0 || width_rf2==0)
					{
						printf("\nERROR: AvgWidth is zero!\n");
						continue;
					}
					smoothDur_rf1=width_rf1;
					smoothDur_rf2=width_rf2;
					printf("\nTaking default ACF half-width rf1 = %f; ACF half-width rf2 = %f\n",smoothDur_rf1,smoothDur_rf2);
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
				{	smoothDur_rf1=atof(pch);
					if(smoothDur_rf1<0)
					{
						printf("\nInvalid smoothing duration!\n");
						continue;
					}	microFileName;
					smoothDur_rf2=smoothDur_rf1;
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
				printf("\nSmoothing window is rf1 %f ms; rf2 %f ms. highfreq cutOff is %f. lowfreq cutOff is %f. Frequency threshold set at %d times RMS\n",smoothDur_rf1,smoothDur_rf2,cutOffhighfreq,cutOfflowfreq,powerThreshold);
				
				scanPulseNum=1;
				FILE *microFile=fopen(microFileName,"w");		
				fprintf(microFile,"pulseNum\trms_rf1\trms_rf2\twidth_rf1\twidth_rf2\tperiod_rf1\tperiod_rf2\trelSt_rf1\trelSt_rf2\ttotSmooth_rf1\ttotSmooth_rf2\n");			
				while(microCandidates(filename_rf1,filename_rf2,microFile,phase_rf1,phase_rf2,s,folda_rf1,folda_rf2,last_rf1,last_rf2,pstart,pend,pulseStart,pulseEnd,intg_rf1,intg_rf2,scanPulseNum,fold,t_rf1,t_rf2,lag,smoothDur_rf1,smoothDur_rf2,rf1Delay,rf2Delay,readData_rf1,readData_rf2,cutOffhighfreq,cutOfflowfreq,width_rf1,width_rf2)==1)
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
					if(width_rf1==0 || width_rf2==0)
					{
						printf("\nERROR: AvgWidth is zero!\n");
						continue;
					}
					smoothDur_rf1=width_rf1;
					smoothDur_rf2=width_rf2;
					printf("\nTaking default ACF half-width _rf1 = %f; ACF half-width _rf2 = %f\n",smoothDur_rf1,smoothDur_rf2);
				}
				/*else if(pch[0]=='s')
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
				{	smoothDur_rf1=atof(pch);
					if(smoothDur_rf1<0)
					{
						printf("\nInvalid smoothing duration!\n");
						continue;
					}	
					smoothDur_rf2=smoothDur_rf1;
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
				
				/*if(smoothDur==-1)
					printf("\nSmoothing window is ACF half-width. cutOff is %f. Frequency threshold set at %d times RMS\n",cutOff,powerThreshold);		
				else*/
				printf("\nSmoothing window is _rf1 %f ms; _rf2 %f ms. cutOffhigh is %f.cutOfflow is %f. Frequency threshold set at %d times RMS. SG order is %d\n",smoothDur_rf1,smoothDur_rf2,cutOffhighfreq,cutOfflowfreq,powerThreshold,sg_coeff_size);
				
				scanPulseNum=1;
				FILE *microFile=fopen(microFileName,"w");		
				fprintf(microFile,"pulseNum\trms_rf1\trms_rf2\tmaxDevRNS_rf1\tmaxDevRMS_rf2\twidth_rf1\twidth_rf2\tperiod_rf1\tperiod_rf2\trelSt_rf1\trelSt_rf2\ttotSmooth_rf1\ttotSmooth_rf2\tMicrostrPower_rf1\tMicrostrPower_rf2\n");			
				while(microCandidates_SG(filename_rf1,filename_rf2,microFile,SG_coeffs,phase_rf1,phase_rf2,s,folda_rf1,folda_rf2,last_rf1,last_rf2,pstart,pend,pulseStart,pulseEnd,intg_rf1,intg_rf2,scanPulseNum,fold,t_rf1,t_rf2,lag,smoothDur_rf1,smoothDur_rf2,rf1Delay,rf2Delay,readData_rf1,readData_rf2,cutOffhighfreq,cutOfflowfreq,width_rf1,width_rf2,sg_coeff_size,option)==1)
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
					if(recordMicro(filename_rf1,outfile_rf1,filename_rf2,outfile_rf2,recPulseNumber,last_rf1,last_rf2,s,folda_rf1,folda_rf2,pstart,pend,phase_rf1,phase_rf2,intg_rf1,intg_rf2,fold,t_rf1,t_rf2,pulseStart,pulseEnd,readData_rf1,readData_rf2,rf1Delay,rf2Delay,recordFile,pulseWindow,foldfile_rf1,foldfile_rf2,profile_rf1,profile_rf2,numshift)==-1)
						break;
				}
				fclose(recordFile);
				fclose(microFile);

			}*/
			else if(input[0]=='>')
			{
				cpgslct(pulseWindow);
				plotSpectra(filename_rf1,filename_rf2,specFile,pulseNum,last_rf1,last_rf2,s,folda_rf1,folda_rf2,pstart,pend,pulseStart,pulseEnd,phase_rf1,phase_rf2,intg_rf1,intg_rf2,fold,t_rf1,t_rf2,readData_rf1,readData_rf2,rf1Delay,rf2Delay);
			}
			/*else if(input[0]=='x')
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
					if(width_rf1==0 || width_rf2==0)
					{
						printf("\nERROR: AvgWidth is zero!\n");
						continue;
					}
					smoothDur_rf1=width_rf1;
					smoothDur_rf2=width_rf2;
					printf("\nTaking default ACF half-width _rf1 = %f; ACF half-width _rf2 = %f\n",smoothDur_rf1,smoothDur_rf2);
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
			/*	else		
				{	smoothDur_rf1=atof(pch);
					if(smoothDur_rf1<0)
					{
						printf("\nInvalid smoothing duration!\n");
						continue;
					}	
					smoothDur_rf2=smoothDur_rf1;
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
				
				cpgslct(autoWindow);*/
				
	/*			if(smoothDur==-1)
					printf("\nPlotting pulse smoothed to ACF half-width and autocorrelation upto lag %f ms.\n",lag);
				else*/
				/*printf("\nPlotting pulse smoothed to _rf1: %f ms _rf2: %f ms, and autocorrelation upto lag %f ms.\n",smoothDur_rf1,smoothDur_rf2,lag);
				
				cpgsch(2.5);
				smoothSubSGACF_intrapolated(filename_rf1,filename_rf2,subAutoFile,SG_coeffs,phase_rf1,phase_rf2,s,folda_rf1,folda_rf2,last_rf1,last_rf2,pstart,pend,pulseStart,pulseEnd,intg_rf1,intg_rf2,pulseNum,fold,t_rf1,t_rf2,lag,smoothDur_rf1,smoothDur_rf2,rf1Delay,rf2Delay,readData_rf1,readData_rf2,sg_coeff_size,option,resicrosscorr);	
			}*/
			else if(input[0]=='j')
			{
				createAverageSpectra(filename_rf1,filename_rf2,avgPulseSpec,foldfile_rf1,foldfile_rf2,last_rf1,last_rf2,s,folda_rf1,folda_rf2,pstart,pend,pulseStart,pulseEnd,phase_rf1,phase_rf2,intg_rf1,intg_rf2,fold,t_rf1,t_rf2,readData_rf1,readData_rf2,rf1Delay,rf2Delay);
			}
			else if(input[0]=='/')
			{
				int adpflag, adpWindow;
		
				if(adpflag==0)
				{
					adpflag=1;
					adpWindow=cpgopen("/xs");
					cpgask(0);
					cpgsubp(1,2);
				}
				cpgslct(adpWindow);
				ADPcreatePulse(filename_rf1,adp_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,2);
				ADPcreatePulse(filename_rf2,adp_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,1);
			}
			/*else if (input[0]=='w')
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
					//pch=strtok(NULL," ");
					//pch=strtok(NULL," ");
					//printf("%s\n",pch);
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
						correlate(filename_rf2,outcorr_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,lag,pulseStart,pulseEnd,readData_rf2,rf2Delay,0);
						cpgpanl(1,1);
						correlate(filename_rf1,outcorr_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,lag,pulseStart,pulseEnd,readData_rf1,rf1Delay,1);
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
					
					writefile(filename_rf1,outfile_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,readData_rf1,rf1Delay);
					writefile(filename_rf2,outfile_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,readData_rf2,rf2Delay);	

					foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,2,intg_rf1,pulseStart,pulseEnd,rf1Delay);
					foldprofile(filename_rf2,foldfile_rf2,fold,t_rf2,phase_rf2,pstart,pend,numshift,profile_rf2,2,intg_rf2,pulseStart,pulseEnd,rf2Delay);
				
					if(plotswitch==0)		
					{
						cpgsch(1.5);
						plottwo(outfile_rf2,outfile_rf1,pulseNum);
					}
					else
					{
						cpgsch(2);
						plotfour(outfile_rf2,foldfile_rf2,outfile_rf1,foldfile_rf1,pulseNum);
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
					plottwo(foldfile_rf2,foldfile_rf1,-1);
						
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
			}*/
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
					averageCrossCorrelate(filename_rf1,filename_rf2,last_rf1,last_rf2,s,folda_rf1,folda_rf2,pstart,pend,pulseStart, pulseEnd,phase_rf1,phase_rf2,intg_rf1,intg_rf2,fold,t_rf1,t_rf2,lag,readData_rf1,readData_rf2,rf1Delay,rf2Delay,strong_rf1,avg_crossCorr);
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
				foldedCrossCorr(foldfile_rf1,foldfile_rf2,pstart,pend,lag,t_rf1,t_rf2,intg_rf1,intg_rf2,fold,last_rf1,last_rf2);		

			}

			else if(input[0]==';')
			{	
				double *folda,*phase;
				int last, rf;

				pch=strtok(input," ");
				pch=strtok(NULL," ");
		
				if(pch==NULL)
				{
					printf("\nRF number not specified \n");
					continue;
				}
				rf=atoi(pch);
				if (rf==1)
				{
					folda=folda_rf1;
					phase=phase_rf1;
					last=last_rf1;
				}
				else if (rf==2)
				{
					folda=folda_rf2;
					phase=phase_rf2;
					last=last_rf2;
				}
				else 
				{
					printf("Invalid RF number. Enter either 1 or 2\n");
					continue;
				}

 
				int i;
				double offmean, IntSignal=0, stdev,count=0,tempPhase,SNR;
				
				offmean=meanBase(folda,phase,pulseStart,pulseEnd,last);
				
				for(i=0;i<last;i++)
				{
					tempPhase=phase[i];
					if(tempPhase>=pulseStart && tempPhase<=pulseEnd)
					{
						IntSignal+=folda[i]-offmean;
						count++;
					}			
				}
				
				stdev=stdevBase(folda,phase,pulseStart,pulseEnd,last);
				SNR=IntSignal/(stdev*sqrt(count));
				printf("\nStandard deviation: %lf",stdev);
				printf("\nBest estimate for SNR: %lf",SNR);
				printf("\nEnergy of pulse: %lf",IntSignal);

				if(SNR>3)
					printf("\nSignal present");
				else
					printf("\nNulling");	
			}
			
			else if(input[0]==':')
			{
				double *folda,*phase, t;
				int last, rf;
				char* filename;
				float* readData;

				pch=strtok(input," ");
				pch=strtok(NULL," ");
		
				if(pch==NULL)
				{
					printf("\nRF number not specified \n");
					continue;
				}
				rf=atoi(pch);
				if (rf==1)
				{
					folda=folda_rf1;
					phase=phase_rf1;
					last=last_rf1;
					t=t_rf1;
					filename=filename_rf1;
					readData=readData_rf1;
				}
				else if (rf==2)
				{
					folda=folda_rf2;
					phase=phase_rf2;
					last=last_rf2;
					t=t_rf2;
					filename=filename_rf2;
					readData=readData_rf1;
				}
				else 
				{
					printf("Invalid RF number. Enter either 1 or 2\n");
					continue;
				}


				FILE *data;
				data=fopen(filename,"rb");
				int pulseNum=1,skipnum=0;
				float nulls=0;

				numshift=(int)(last*s); //garvit 15/11/21

				while(1)
				{
					skipnum=(int)((pulseNum-1)*(fold/t));
					fseek(data,(skipnum)*sizeof(float),SEEK_SET);
			
					readFloatToDouble(numshift,readData,folda,data); //garvit 15/11/21
					if(readFloatToDouble(last,readData,folda,data)!=1)
						break;
			
					int i;
					double offmean, IntSignal=0, stdev,count=0,tempPhase,SNR;
			
					offmean=meanBase(folda,phase,pulseStart,pulseEnd,last);
			
					for(i=0;i<last;i++)
					{
						tempPhase=phase[i];
						if(tempPhase>=pulseStart && tempPhase<=pulseEnd)
						{
							IntSignal+=folda[i]-offmean;							
							count++;
						}			
					}
				
			
					stdev=stdevBase(folda,phase,pulseStart,pulseEnd,last);
					SNR= IntSignal/(stdev*sqrt(count));

					if (SNR<=3)
					{
						printf("\nNulling in pulse number %d with SNR %f",pulseNum,SNR);
						nulls++;
					}
					pulseNum++;
				}

				fclose(data);
				printf("\nNumber of pulses = %d",pulseNum-1);
				printf("\nNulling fraction is %f %",(nulls*100)/(pulseNum-1));
			}

			else if(input[0]=='w')
			{
				double *folda,*phase, t;
				int last, rf;
				char* filename;
				float* readData;

				pch=strtok(input," ");
				pch=strtok(NULL," ");
		
				if(pch==NULL)
				{
					printf("\nRF number not specified \n");
					continue;
				}
				rf=atoi(pch);
				if (rf==1)
				{
					folda=folda_rf1;
					phase=phase_rf1;
					last=last_rf1;
					t=t_rf1;
					filename=filename_rf1;
					readData=readData_rf1;
				}
				else if (rf==2)
				{
					folda=folda_rf2;
					phase=phase_rf2;
					last=last_rf2;
					t=t_rf2;
					filename=filename_rf2;
					readData=readData_rf1;
				}
				else 
				{
					printf("Invalid RF number. Enter either 1 or 2\n");
					continue;
				}

				
				char *energyhist="energies.txt";
				FILE *data,*energy;
				data=fopen(filename,"rb");
				energy=fopen(energyhist,"w");
				int pulseNum=1,skipnum=0,NumofPulses;
				
				while(1)
				{
					skipnum=(int)((pulseNum-1)*(fold/t));
					fseek(data,(skipnum)*sizeof(float),SEEK_SET);
			
					if(readFloatToDouble(last,readData,folda,data)!=1)
						break;
					pulseNum++;
				}
				NumofPulses = pulseNum-1;	
				
				float *onmeanVal=(float*)malloc(NumofPulses*sizeof(float));
				float *offmeanVal=(float*)malloc(NumofPulses*sizeof(float));
				
				pulseNum=1;
				skipnum=0;	

				numshift=(int)(last*s); //garvit 15/11/21

				float offpulse_start, offpulse_end;
				printf("Enter off-pulse start: ");
				scanf("%f",&offpulse_start);
				printf("Enter off-pulse end: ");
				scanf("%f",&offpulse_end);

				while(1)
				{
					skipnum=(int)((pulseNum-1)*(fold/t));
					fseek(data,skipnum*sizeof(float),SEEK_SET);
				
					readFloatToDouble(numshift,readData,folda,data); //garvit 15/11/21
					if(readFloatToDouble(last,readData,folda,data)!=1)
						break;
					
					int i;
					double offmean=0,onmean=0,count=0,offcount=0,tempPhase;
					double baseline=meanBase(folda,phase,pulseStart,pulseEnd,last);
					double pulseWindow=pulseEnd-pulseStart;
					
					for(i=0;i<last;i++)
					{
						tempPhase=phase[i];
						if(tempPhase>=offpulse_start && tempPhase<=offpulse_end)
						{
							offmean+=folda[i];
							offcount++;
						}					


						/*if(pulseEnd+0.1+pulseWindow<1)
						{	tempPhase=phase[i];
							if(tempPhase>=pulseEnd+0.1 && tempPhase<=pulseEnd+0.1+pulseWindow)
							{
								offmean+=folda[i];
								offcount++;
							}
						}
						else
						{	tempPhase=phase[i];
							if((tempPhase>=pulseEnd+0.1 && tempPhase<=1) || (tempPhase>=0 && tempPhase<=pulseWindow-(1-pulseEnd-0.1)))
							{
								offmean+=folda[i];
								offcount++;
							}
						} */
					}
					
					for(i=0;i<last;i++)
					{
						tempPhase=phase[i];
						if(tempPhase>=pulseStart && tempPhase<=pulseEnd)
						{
							onmean+=folda[i];
							count++;
							
						}			
					}
					
					onmeanVal[pulseNum-1]=onmean/count-baseline;
					offmeanVal[pulseNum-1]=offmean/offcount-baseline;
					pulseNum++;
				}
				
				printf("Enter number of pulses to plot histogram for (max %d): \t",NumofPulses);
				int plotpulses;
				scanf("%d",&plotpulses);

				double onmeanmin=onmeanVal[0],onmeanmax=0,offmeanmin=offmeanVal[0],offmeanmax=0;
				for(i=0;i<plotpulses;i++)
				{
					if(onmeanVal[i]<onmeanmin)
						onmeanmin=onmeanVal[i];
					if(onmeanVal[i]>onmeanmax)
						onmeanmax=onmeanVal[i];
					if(offmeanVal[i]<offmeanmin)
						offmeanmin=offmeanVal[i];
					if(offmeanVal[i]>offmeanmax)
						offmeanmax=offmeanVal[i];

					fprintf(energy,"%f\t%f\n",(float)onmeanVal[i],(float)offmeanVal[i]);
				}
				
				printf("\nOn pulse mean:\nMin: %lf\nMax: %lf\nOff pulse mean:\nMin: %lf\nMax: %lf\n",onmeanmin,onmeanmax,offmeanmin,offmeanmax);
				const float *onmeanVal2 = onmeanVal;
				const float *offmeanVal2 = offmeanVal;
				
				cpgslct(pulseWindow);
				cpgpap(10.0,0.8);
				cpgsubp(1,2);
				//min2(onmeanmin,offmeanmin),max2(onmeanmax,offmeanmax)	
				
				cpgsci(1);

				printf("Enter number of bins for the on-pulse energy histogram (max 200): \t");
				int nbins_on;
				scanf("%d",&nbins_on);

				printf("Enter number of bins for the off-pulse energy histogram (max 200): \t");
				int nbins_off;
				scanf("%d",&nbins_off);

				cpgpanl(1,2);
				cpgenv(onmeanmin,onmeanmax,0,0.05*plotpulses,0,0);
				char title_on[50];
				sprintf(title_on,"On-pulse energy Histogram");
				cpglab("Energy","Frequency",title_on);		
				cpghist(plotpulses,onmeanVal2,onmeanmin,onmeanmax,nbins_on,3);
				//cpgsci(2);

				
				cpgpanl(1,1);
				cpgenv(offmeanmin,offmeanmax,0,0.05*plotpulses,0,0);
				char title_off[50];
				sprintf(title_off,"Off-pulse energy Histogram");
				cpglab("Energy","Frequency",title_off);
				cpghist(plotpulses,offmeanVal2,offmeanmin,offmeanmax,nbins_off,3);		
				
			}

			else if(input[0]=='x')
			{	

				double *folda,*phase, t, *profile;
				int last, rf;
				char* filename, *foldfile;
				float* readData;

				pch=strtok(input," ");
				pch=strtok(NULL," ");
		
				if(pch==NULL)
				{
					printf("\nRF number not specified \n");
					continue;
				}
				rf=atoi(pch);
				if (rf==1)
				{
					folda=folda_rf1;
					phase=phase_rf1;
					last=last_rf1;
					t=t_rf1;
					filename=filename_rf1;
					readData=readData_rf1;
					foldfile=foldfile_rf1;
					profile=profile_rf1;
				}
				else if (rf==2)
				{
					folda=folda_rf2;
					phase=phase_rf2;
					last=last_rf2;
					t=t_rf2;
					filename=filename_rf2;
					readData=readData_rf2;
					foldfile=foldfile_rf2;
					profile=profile_rf2;
				}
				else 
				{
					printf("Invalid RF number. Enter either 1 or 2\n");
					continue;
				}

				float phaseStart_gray, phaseEnd_gray;
				printf("Enter phase start:\t");
				scanf("%f",&phaseStart_gray);
				printf("Enter phase end:\t");
				scanf("%f",&phaseEnd_gray);

				int onPulsebins,collapsePulses;
				printf("\nEnter number of bins for the phase range %f to %f:\t",phaseStart_gray,phaseEnd_gray);
				scanf("%d",&onPulsebins);
				plotDrifts(filename, folda, phase, phaseStart_gray, phaseEnd_gray, readData, onPulsebins, last, s, fold, t, foldfile, profile);
							
			}

			

			else if(input[0]=='r') 
			{
				rerun=true;
				break;		
			}
			else
			{
				printf("\nInvalid option. Try again.\n");
			}
		
		}

		if (rerun==true) continue;

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
		char rf1file[1000],rf2file[10000],file[1000];
		float reso,rf1reso,rf2reso,pulseper,dm,timeOffset;
		FILE *parameter = fopen("/Data/ygupta/psr_microstr_work/recovered_code/dual_freq_analysis/parameters.pulse_analysis", "r");
		
		if (parameter == NULL) 
		{
		perror("Error: Failed to open file.");
		return 1;
		}

		char *parameters[50],line[600], pulsar[15];
		for(i=0;i<50;i++) parameters[i]=(char*) malloc(sizeof(char)*600);
		int num_pulsars=0, num=0;	

		for(i=0;fgets(line,sizeof(line),parameter);i++)	
		{
			strcpy(parameters[i],line);
			sscanf(parameters[i], "%d)%s %*s", &num,pulsar);
			printf("%d) %s\n",num,pulsar);
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
		sscanf(parameters[pulsarnum-1], "%*d)%s %f %f %f %f %f %s %s %*s %*s",pulsarchosen,&timeOffset,&rf1reso,&rf2reso,&pulseper,&dm,rf1file,rf2file);
	//"%d)%s %f %f %f %s",&num,pulsarchosen,&reso,&pulseper,&dm,datafile

		printf("\nPress 1 for first frequency data \nPress 2 for second frequency data: ");
		int freq;
		scanf("%d",&freq);

		char* filename;
		if (freq==1)
		{
			reso=rf1reso;
			filename=(char*)(rf1file);
		}
	
		else if (freq==2)
		{	
			reso=rf2reso;
			filename=(char*)(rf2file);
		}

		else 
		{
			printf("\nWrong Input");
			return 1;
		}


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

		//char *filename=datafile;
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

		/*double t_rf1,t_rf2,fold,value,time,k,lag=0,DM,smoothDur_rf1=0,smoothDur_rf2=0;
		int last_rf1,last_rf2,i,sample=0,count=0,scanPulseNum=0;
		int b=0;
		char *input=(char*)malloc(100*sizeof(char));*/

		/*t_rf1=(double)atof(argv[3]);
		t_rf2=(double)atof(argv[4]);
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
			printf("\nEnter option:\t\tUse 'k' to get instructions & 'q' to terminate program\n");
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
			else if(input[0]=='s')
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
				
				writefile(filename,outfile,pulseNum,last,s,folda,pstart,pend,phase,intg,fold,t,pulseStart,pulseEnd,readData,0);
				foldprofile(filename,foldfile,fold,t,phase,pstart,pend,curshift,profile,3,intg,pulseStart,pulseEnd,0);
				foldprofile(filename,foldfile,fold,t,phase,pstart,pend,numshift,profile,2,intg,pulseStart,pulseEnd,0);
				
				cpgslct(pulseWindow);
				if(plotswitch==0)			
					plotone_single(outfile,pulseNum,1);
				else
					plottwo_single(foldfile,outfile,pulseNum);
			}
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
				printf("\nNew integration is %d samples for _rf2\n", intg);
				
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
				printf("\nInstructions:\nUse 'n' to go to next pulse\nUse 'b' to go the earlier pulse\nUse 'p <start phase> <end phase>' to set a phase range (Default is full period)\nUse 's <phase shift>' to shift phase window (Phase shift should remain between 0 and 1 !)\nUse 'f' to show folded profile\nUse 'h' to hide folded profile\nUse 'i <Integration>' to integrate profile\nUse 'g <Pulse Number>' to get profile for a specified pulse number\nUse ‘m <min SNR>’ to output list of pulses above some average SNR \nUse 'a <Time lag>' to get autocorrelation plot for the current pulse, integration and phase range\nUse 't <MinSNR>' to get list of pulse numbers with peak above a given SNR\nUse 'v <Time lag(ms)> to get the average ACF for the strongest pulses\nUse 'y' to view the ADP for the specified window\nUse 'k' to get instructions\nUse 'q' to terminate program\nUse 'e <time lag>' to get average cross correlation\nUse 'r' to restart the program\n\nSUBTRACTION FEATURE:\n\nUse 'd <smoothing (ms)> <microCuttOff(SNR)> <SG order>’ to run Analysis of Current Pulse, that smoothes the pulse, subtracts the smoothed version, and then computes widths from the subtracted ACF smoothed with a Savitzky Golay filter and computes periods from the power spectrum of the subtracted features\n\n");			
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
					correlate(filename_rf2,outcorr_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,lag,pulseStart,pulseEnd,readData_rf2,rf2Delay,0);
					cpgpanl(1,1);
					correlate(filename_rf1,outcorr_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,lag,pulseStart,pulseEnd,readData_rf1,rf1Delay,1);
					cpgpanl(1,2);
					crossCorrelate(filename_rf1,filename_rf2,crossCorr,pulseNum,last_rf1,last_rf2,s,folda_rf1,folda_rf2,pstart,pend,pulseStart,pulseEnd,phase_rf1,phase_rf2,intg_rf1,intg_rf2,fold,t_rf1,t_rf2,lag,readData_rf1,readData_rf2,rf1Delay,rf2Delay);
					
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
					if(width ==0)
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
				
				if(autoflag==0)
				{
					autoWindow=cpgopen("/xs");
					cpgask(0);
					autoflag=1;
				}
				
				cpgslct(autoWindow);
				

				printf("\nPlotting pulse smoothed to: %f ms and autocorrelation upto lag %f ms.\n",smoothDur,lag);
				
				cpgsch(2.5);
				smoothSubACF(filename,filename,subAutoFile,phase,phase,s,folda,folda,last,last,pstart,pend,pulseStart,pulseEnd,intg,intg,pulseNum,fold,t,t,lag,smoothDur,smoothDur,0,0,readData,readData);
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
			else if(input[0]==',')
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
					if(width==0 )
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
					}	microFileName;
				}
				
				pch=strtok(NULL," ");
				if(pch==NULL)
				{
					printf("\nNo high freq cutOff specified!\n");
					continue;
				}
				cutOff=atof(pch);
				
				if(cutOff<0)
				{
					printf("\nInvalid cut off!\n");
					continue;
				}

				/*if(pch==NULL)
				{
					printf("\nNo low freq cutOff specified!\n");
					continue;
				}
				cutOfflowfreq=atof(pch);
				
				if(cutOfflowfreq<0)
				{
					printf("\nInvalid cut off!\n");
					continue;
				}*/
				

				printf("\nSmoothing window is %f ms. cutOff is %f. Frequency threshold set at %d times RMS\n",smoothDur,cutOff,powerThreshold);
				
				scanPulseNum=1;
				FILE *microFile=fopen(microFileName,"w");		
				fprintf(microFile,"pulseNum\trms\trms\twidth\twidth\tperiod\tperiod\trelSt\trelSt\ttotSmooth\ttotSmooth\n");			
				while(microCandidates(filename,filename,microFile,phase,phase,s,folda,folda,last,last,pstart,pend,pulseStart,pulseEnd,intg,intg,scanPulseNum,fold,t,t,lag,smoothDur,smoothDur,0,0,readData,readData,cutOff,cutOff,width,width)==1)
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
					if(recordMicro(filename_rf1,outfile_rf1,filename_rf2,outfile_rf2,recPulseNumber,last_rf1,last_rf2,s,folda_rf1,folda_rf2,pstart,pend,phase_rf1,phase_rf2,intg_rf1,intg_rf2,fold,t_rf1,t_rf2,pulseStart,pulseEnd,readData_rf1,readData_rf2,rf1Delay,rf2Delay,recordFile,pulseWindow,foldfile_rf1,foldfile_rf2,profile_rf1,profile_rf2,numshift)==-1)
						break;
				}
				fclose(recordFile);
				fclose(microFile);

			}*/
			/*else if(input[0]=='x')
			{
				cpgslct(pulseWindow);
				plotSpectra(filename,filename,specFile,pulseNum,last,last,s,folda,folda,pstart,pend,pulseStart,pulseEnd,phase,phase,intg,intg,fold,t,t,readData,readData,0,0);
			}*/
			else if(input[0]=='j')
			{
				createAverageSpectra(filename,filename,avgPulseSpec,foldfile,foldfile,last,last,s,folda,folda,pstart,pend,pulseStart,pulseEnd,phase,phase,intg,intg,fold,t,t,readData,readData,0,0);
			}

			/*else if (input[0]=='w')
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
						correlate(filename_rf2,outcorr_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,lag,pulseStart,pulseEnd,readData_rf2,rf2Delay,0);
						cpgpanl(1,1);
						correlate(filename_rf1,outcorr_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,lag,pulseStart,pulseEnd,readData_rf1,rf1Delay,1);
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
					
					writefile(filename_rf1,outfile_rf1,pulseNum,last_rf1,s,folda_rf1,pstart,pend,phase_rf1,intg_rf1,fold,t_rf1,pulseStart,pulseEnd,readData_rf1,rf1Delay);
					writefile(filename_rf2,outfile_rf2,pulseNum,last_rf2,s,folda_rf2,pstart,pend,phase_rf2,intg_rf2,fold,t_rf2,pulseStart,pulseEnd,readData_rf2,rf2Delay);	

					foldprofile(filename_rf1,foldfile_rf1,fold,t_rf1,phase_rf1,pstart,pend,numshift,profile_rf1,2,intg_rf1,pulseStart,pulseEnd,rf1Delay);
					foldprofile(filename_rf2,foldfile_rf2,fold,t_rf2,phase_rf2,pstart,pend,numshift,profile_rf2,2,intg_rf2,pulseStart,pulseEnd,rf2Delay);
				
					if(plotswitch==0)		
					{
						cpgsch(1.5);
						plottwo(outfile_rf2,outfile_rf1,pulseNum);
					}
					else
					{
						cpgsch(2);
						plotfour(outfile_rf2,foldfile_rf2,outfile_rf1,foldfile_rf1,pulseNum);
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
					plottwo(foldfile_rf2,foldfile_rf1,-1);
						
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
			}*/
			/*else if (input[0]=='e')
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
					averageCrossCorrelate(filename,filename,last,last,s,folda,folda,pstart,pend,pulseStart, pulseEnd,phase,phase,intg,intg,fold,t,t,lag,readData,readData,0,0,strong,avg_crossCorr);
				}
			}*/
		
			else if(input[0]==';')
			{	 
				int i;
				double offmean, IntSignal=0, stdev,count=0,tempPhase,SNR;
				
				offmean=meanBase(folda,phase,pulseStart,pulseEnd,last);
				
				for(i=0;i<last;i++)
				{
					tempPhase=phase[i];
					if(tempPhase>=pulseStart && tempPhase<=pulseEnd)
					{
						IntSignal+=folda[i]-offmean;
						count++;
					}			
				}
				
				stdev=stdevBase(folda,phase,pulseStart,pulseEnd,last);
				SNR=IntSignal/(stdev*sqrt(count));
				printf("\nStandard deviation: %lf",stdev);
				printf("\nBest estimate for SNR: %lf",SNR);
				printf("\nEnergy of pulse: %lf",IntSignal);

				if(SNR>3)
					printf("\nSignal present");
				else
					printf("\nNulling");	
			}
			
			else if(input[0]==':')
			{
				FILE *data;
				data=fopen(filename,"rb");
				int pulseNum=1,skipnum=0;
				float nulls=0;

				numshift=(int)(last*s); //garvit 15/11/21

				while(1)
				{
					skipnum=(int)((pulseNum-1)*(fold/t));
					fseek(data,(skipnum)*sizeof(float),SEEK_SET);
			
					readFloatToDouble(numshift,readData,folda,data); //garvit 15/11/21
					if(readFloatToDouble(last,readData,folda,data)!=1)
						break;
			
					int i;
					double offmean, IntSignal=0, stdev,count=0,tempPhase,SNR;
			
					offmean=meanBase(folda,phase,pulseStart,pulseEnd,last);
			
					for(i=0;i<last;i++)
					{
						tempPhase=phase[i];
						if(tempPhase>=pulseStart && tempPhase<=pulseEnd)
						{
							IntSignal+=folda[i]-offmean;							
							count++;
						}			
					}
				
			
					stdev=stdevBase(folda,phase,pulseStart,pulseEnd,last);
					SNR= IntSignal/(stdev*sqrt(count));

					if (SNR<=3)
					{
						printf("\nNulling in pulse number %d with SNR %f",pulseNum,SNR);
						nulls++;
					}
					pulseNum++;
				}

				fclose(data);
				printf("\nNumber of pulses = %d",pulseNum-1);
				printf("\nNulling fraction is %f %",(nulls*100)/(pulseNum-1));
			}

			else if(input[0]=='w')
			{
				char *energyhist="energies.txt";
				FILE *data,*energy;
				data=fopen(filename,"rb");
				energy=fopen(energyhist,"w");
				int pulseNum=1,skipnum=0,NumofPulses;
				
				while(1)
				{
					skipnum=(int)((pulseNum-1)*(fold/t));
					fseek(data,(skipnum)*sizeof(float),SEEK_SET);
			
					if(readFloatToDouble(last,readData,folda,data)!=1)
						break;
					pulseNum++;
				}
				NumofPulses = pulseNum-1;	
				
				float *onmeanVal=(float*)malloc(NumofPulses*sizeof(float));
				float *offmeanVal=(float*)malloc(NumofPulses*sizeof(float));
				
				pulseNum=1;
				skipnum=0;	

				numshift=(int)(last*s); //garvit 15/11/21

				float offpulse_start, offpulse_end;
				printf("Enter off-pulse start: ");
				scanf("%f",&offpulse_start);
				printf("Enter off-pulse end: ");
				scanf("%f",&offpulse_end);

				while(1)
				{
					skipnum=(int)((pulseNum-1)*(fold/t));
					fseek(data,skipnum*sizeof(float),SEEK_SET);
				
					readFloatToDouble(numshift,readData,folda,data); //garvit 15/11/21
					if(readFloatToDouble(last,readData,folda,data)!=1)
						break;
					
					int i;
					double offmean=0,onmean=0,count=0,offcount=0,tempPhase;
					double baseline=meanBase(folda,phase,pulseStart,pulseEnd,last);
					double pulseWindow=pulseEnd-pulseStart;
					
					for(i=0;i<last;i++)
					{
						tempPhase=phase[i];
						if(tempPhase>=offpulse_start && tempPhase<=offpulse_end)
						{
							offmean+=folda[i];
							offcount++;
						}					


						/*if(pulseEnd+0.1+pulseWindow<1)
						{	tempPhase=phase[i];
							if(tempPhase>=pulseEnd+0.1 && tempPhase<=pulseEnd+0.1+pulseWindow)
							{
								offmean+=folda[i];
								offcount++;
							}
						}
						else
						{	tempPhase=phase[i];
							if((tempPhase>=pulseEnd+0.1 && tempPhase<=1) || (tempPhase>=0 && tempPhase<=pulseWindow-(1-pulseEnd-0.1)))
							{
								offmean+=folda[i];
								offcount++;
							}
						} */
					}
					
					for(i=0;i<last;i++)
					{
						tempPhase=phase[i];
						if(tempPhase>=pulseStart && tempPhase<=pulseEnd)
						{
							onmean+=folda[i];
							count++;
							
						}			
					}
					
					onmeanVal[pulseNum-1]=onmean/count-baseline;
					offmeanVal[pulseNum-1]=offmean/offcount-baseline;
					pulseNum++;
				}
				
				printf("Enter number of pulses to plot histogram for (max %d): \t",NumofPulses);
				int plotpulses;
				scanf("%d",&plotpulses);

				double onmeanmin=onmeanVal[0],onmeanmax=0,offmeanmin=offmeanVal[0],offmeanmax=0;
				for(i=0;i<plotpulses;i++)
				{
					if(onmeanVal[i]<onmeanmin)
						onmeanmin=onmeanVal[i];
					if(onmeanVal[i]>onmeanmax)
						onmeanmax=onmeanVal[i];
					if(offmeanVal[i]<offmeanmin)
						offmeanmin=offmeanVal[i];
					if(offmeanVal[i]>offmeanmax)
						offmeanmax=offmeanVal[i];

					fprintf(energy,"%f\t%f\n",(float)onmeanVal[i],(float)offmeanVal[i]);
				}
				
				printf("\nOn pulse mean:\nMin: %lf\nMax: %lf\nOff pulse mean:\nMin: %lf\nMax: %lf\n",onmeanmin,onmeanmax,offmeanmin,offmeanmax);
				const float *onmeanVal2 = onmeanVal;
				const float *offmeanVal2 = offmeanVal;
				
				cpgslct(pulseWindow);
				cpgpap(10.0,0.8);
				cpgsubp(1,2);
				//min2(onmeanmin,offmeanmin),max2(onmeanmax,offmeanmax)	
				
				cpgsci(1);

				printf("Enter number of bins for the on-pulse energy histogram (max 200): \t");
				int nbins_on;
				scanf("%d",&nbins_on);

				printf("Enter number of bins for the off-pulse energy histogram (max 200): \t");
				int nbins_off;
				scanf("%d",&nbins_off);

				cpgpanl(1,2);
				cpgenv(onmeanmin,onmeanmax,0,0.05*plotpulses,0,0);
				char title_on[50];
				sprintf(title_on,"On-pulse energy Histogram");
				cpglab("Energy","Frequency",title_on);		
				cpghist(plotpulses,onmeanVal2,onmeanmin,onmeanmax,nbins_on,3);
				//cpgsci(2);

				
				cpgpanl(1,1);
				cpgenv(offmeanmin,offmeanmax,0,0.05*plotpulses,0,0);
				char title_off[50];
				sprintf(title_off,"Off-pulse energy Histogram");
				cpglab("Energy","Frequency",title_off);
				cpghist(plotpulses,offmeanVal2,offmeanmin,offmeanmax,nbins_off,3);		
				
			}

			else if(input[0]=='x')
			{	
				float phaseStart_gray, phaseEnd_gray;
				printf("Enter phase start:\t");
				scanf("%f",&phaseStart_gray);
				printf("Enter phase end:\t");
				scanf("%f",&phaseEnd_gray);

				int onPulsebins,collapsePulses;
				printf("\nEnter number of bins for the phase range %f to %f:\t",phaseStart_gray,phaseEnd_gray);
				scanf("%d",&onPulsebins);
				plotDrifts(filename, folda, phase, phaseStart_gray, phaseEnd_gray, readData, onPulsebins, last, s, fold, t, foldfile, profile);
							
			}

			else if(input[0]=='r') 
			{
				rerun=true;
				break;		
			}
			else
			{
				printf("\nInvalid option. Try again.\n\n");
			}
		
		}

		if (rerun==true) continue;

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
}
