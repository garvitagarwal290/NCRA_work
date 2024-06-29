void BasicAnalysis::subtractZeroDM(char* freqFlags,float centralTendency)
{
	float*	ptrZeroDM;
	float*	ptrCorrelationBandshape;
	float*	ptrMeanToRmsBandshape;
	float*  ptrNormalizedBandshape;
	float* 	ptrRawData;
	float 	count=0;					//Stores the number of channels added to get each time sample	
	int 	startChannel=info.startChannel;
	int 	nChan= info.stopChannel-startChannel;		//Number of channels to use
	int 	endExclude=info.noOfChannels-info.stopChannel;	//Number of channels to exclude from the end of the band
	int 	l= blockLength;
	float 	zeroDMMean,zeroDMRMS;
	for(int i=0;i<nChan;i++)
		if(!freqFlags[i])
			count++;
	ptrZeroDM=zeroDM;
	zeroDMMean=0.0;	
	//zeroDMRMS=0.0;
	for(int i=0;i<l;i++,ptrZeroDM++)
	{
		zeroDMMean+=*ptrZeroDM;
		//zeroDMRMS+=(*ptrZeroDM)*(*ptrZeroDM);
	}
	zeroDMMean/=l;
	//zeroDMRMS=zeroDMRMS/l-zeroDMMean*zeroDMMean;
	ptrZeroDM=zeroDM;
	ptrRawData=rawData;
	ptrCorrelationBandshape=correlationBandshape;
	for(int j=0;j<info.noOfChannels;j++,ptrCorrelationBandshape++)		
		*ptrCorrelationBandshape=0;
	for(int i=0;i<l;i++,ptrZeroDM++)
	{	
		ptrRawData+=startChannel;			//startChannel number of channels skipped at the start of the band
		ptrCorrelationBandshape=correlationBandshape+startChannel;
		ptrNormalizedBandshape=normalizedBandshape+startChannel;
		for(int j=0;j<nChan;j++,ptrRawData++,ptrCorrelationBandshape++,ptrNormalizedBandshape++)		
			*ptrCorrelationBandshape+=(*ptrZeroDM-zeroDMMean)*(*ptrRawData-*ptrNormalizedBandshape);
		ptrRawData+=endExclude;				//endExclude number of channels skipped at the end of the band
	}
	ptrZeroDM=zeroDM;
	ptrRawData=rawData;
	for(int i=0;i<l;i++,ptrZeroDM++)
	{	
		ptrRawData+=startChannel;			//startChannel number of channels skipped at the start of the band
		ptrCorrelationBandshape=correlationBandshape+startChannel;
		ptrNormalizedBandshape=normalizedBandshape+startChannel;
		//ptrMeanToRmsBandshape=meanToRmsBandshape+startChannel;
		for(int j=0;j<nChan;j++,ptrRawData++,ptrNormalizedBandshape++,ptrCorrelationBandshape++)//,ptrMeanToRmsBandshape++)		
			*ptrRawData=*ptrRawData-(*ptrZeroDM-zeroDMMean)*(*ptrCorrelationBandshape);
		ptrRawData+=endExclude;				//endExclude number of channels skipped at the end of the band
	}
}
