#include<AquireData.h>

#include<custom_cout.h> //7th Jan 2022 Garvit Agarwal. to use the custom cout for gbmon log

//Declaring static variables:
Information 	AquireData::info;
long int 	AquireData::eof;
DasHdrType*	AquireData::dataHdr;
DataBufType*	AquireData::dataBuffer;
DataTabType*	AquireData::dataTab;
unsigned short* AquireData::zeros;
int		AquireData::recNum=0;
int		AquireData::remainingData=0;
int		AquireData::currentReadBlock=0;
long int	AquireData::curPos=0;
double		AquireData::totalError=0.0;
float		AquireData::buffSizeSec;
int		AquireData::nbuff;		
struct timeval*	AquireData::startTimeStamp;


/*******************************************************************
*CONSTRUCTOR: AquireData::AquireData(Information _info)
*Information _info : contains all parameters.
*Invoked on first creation of an object of this type.
********************************************************************/
AquireData::AquireData(Information _info)
{
	 info=_info;
	 blockIndex=0;	 
	 hasReachedEof=0;
	 totalError=0.0;
	 remainingData=0;
	 //finding length of data file or initializing SHM according to the mode of operation
	 if(info.doReadFromFile)
	 {
	 	ifstream datafile;
	 	datafile.open(info.filepath,ios::binary);
		if(!datafile.is_open())
  		{
    			custom_cout::cout<<"Raw data file not found!"<<endl;
    			exit(1);
 		}
	 	datafile.seekg(0,ios::end);
	 	eof=datafile.tellg();
	 	datafile.close();
		if(eof==0)
		{
			custom_cout::cout<<"DATA FILE EMPTY"<<endl;
			exit(1);
		}
	 }
	// else
	//	initializeSHM();
	
	 
}
/*******************************************************************
*CONSTRUCTOR: AquireData::AquireData(int _blockIndex)
*int _blockIndex : current window number.
********************************************************************/
AquireData::AquireData(int _blockIndex)
{
	blockIndex=_blockIndex; 
	hasReachedEof=0;
	rawData=NULL;
	rawDataFloat=NULL;
	rawDataPolar=NULL;
	splittedRawData=NULL;
	nBuffTaken=0;
	if(info.isInline)
		headerInfo=new char[4096*nbuff];
	headerInfo=NULL;
}
/*******************************************************************
*DESTRUCTOR: AquireData::~AquireData()
*Free all memory
********************************************************************/
AquireData::~AquireData()
{	
	switch(info.sampleSizeBytes)
	{
		case 1:
			if(!info.doPolarMode)
			{
				delete[] rawDataChar;
			}
			else
				delete[] rawDataCharPolar;
			break;
		case 2:
			if(!info.doPolarMode)
				delete[] rawData;
			else
				delete[] rawDataPolar;
			break;
		case 4:		
			delete[] rawDataFloat;
			break;
	}
	
	//delete[] splittedRawData;
	
}
/*******************************************************************
*FUNCTION: void BasicAnalysis::readData()
*Wrapper to decide which method to call to aquire data
*******************************************************************/
void AquireData::readData()
{
	double er=info.blockSizeSec/info.samplingInterval;
	er=er-(long)er;
	totalError+=er;	
	blockLength = info.blockSizeSamples;	
  	if(!info.isInline && totalError>=1.0)
  	{
  		blockLength++;
  		totalError--;
  	}	
	if(info.doReadFromFile)
		readDataFromFile();
	else
		readFromSHM();
}
/*******************************************************************
*FUNCTION: void AquireData::initializeSHM()
*Initializes collect_psr shared memory of GSB/GWB
*******************************************************************/
void AquireData::initializeSHM()
{
  	int iddataHdr,idDataBuffer;
	switch(info.shmID)
	{
		case 1: //standard correlator shm
			iddataHdr = shmget(DAS_H_KEY,sizeof(DasHdrType),SHM_RDONLY);
			idDataBuffer = shmget(DAS_D_KEY,sizeof(DataBufType),SHM_RDONLY);
		break;
		case 2: //file shm
			iddataHdr = shmget(DAS_H_KEY_GPTOOL,sizeof(DasHdrType),SHM_RDONLY);
			idDataBuffer = shmget(DAS_D_KEY_GPTOOL,sizeof(DataBufType),SHM_RDONLY);
		break;
		case 3: //inline gptool shm
			iddataHdr = shmget(DAS_H_KEY_GPTOOL_INLINE,sizeof(DasHdrType),SHM_RDONLY);
			idDataBuffer = shmget(DAS_D_KEY_GPTOOL_INLINE,sizeof(DataBufType),SHM_RDONLY);
		break;
	}
	
	custom_cout::cout<<"iddataHdr="<<iddataHdr<<endl<<"idDataBuffer="<<idDataBuffer<<endl;
	
	if(iddataHdr < 0 || iddataHdr < 0)
	{
		exit(1);
	}
	dataHdr = (DasHdrType*)shmat(iddataHdr,0,SHM_RDONLY);
  	dataBuffer = (DataBufType*)shmat(idDataBuffer,0,SHM_RDONLY);
	if((dataBuffer)==(DataBufType*)-1)
	{
		custom_cout::cout<<"Cannot attach to BUFFER 3 shared memory!"<<endl;
		exit(1);
	}
	custom_cout::cout<<"Attached to shared memory:"<<dataHdr<<","<<dataBuffer<<endl;
	dataTab = dataBuffer-> dtab;
	custom_cout::cout<<"Max no of blocks="<<dataBuffer->maxblocks<<endl;
	zeros=new unsigned short[dataBuffer->blocksize-4096];
	if(info.isInline)
	{	
		buffSizeSec=(info.samplingInterval*(dataBuffer->blocksize-4096)/(info.sampleSizeBytes*info.noOfChannels));
		nbuff=(int)floor(info.blockSizeSec/buffSizeSec);
		info.blockSizeSec=nbuff*buffSizeSec;
		info.blockSizeSamples=(int)round(info.blockSizeSec/info.samplingInterval);

	}
	 /*   find a block, two blocks before the current block of the shm for reading data */
	//if(dataBuffer->cur_rec > (dataBuffer->maxblocks)/2)
	recNum = (dataBuffer->cur_rec+MaxDataBuf)%MaxDataBuf;
	currentReadBlock = dataTab[recNum].seqnum;
	
	int timeOff=sizeof(double)+sizeof(int)+sizeof(int)+sizeof(int);
	int blk_nano_off=timeOff+sizeof(struct timeval)+sizeof(int);
	startTimeStamp=(struct timeval*)(dataBuffer->buf+dataTab[recNum].rec*(dataBuffer->blocksize)+timeOff);
	
	int usec=0;
  	double nano_usec;
	struct tm *local_t = localtime(&startTimeStamp->tv_sec);
//  usec = 1000000*(timestamp_gps->tv_usec);
	double *blk_nano = (double *)(dataBuffer->buf+dataTab[recNum].rec*(dataBuffer->blocksize) + blk_nano_off);
  	usec = (startTimeStamp->tv_usec);
 	nano_usec = local_t->tm_sec; // putting int value in double:CAUTION
  	nano_usec += (startTimeStamp->tv_usec)/1e6;
  	nano_usec += (*blk_nano)/1e6;
	FILE* FATM=fopen("timestamp.gpt","w");
        char time_string[40];
        fprintf(FATM,"#Start time and date\n");
        fprintf(FATM,"IST Time: %02d:%02d:%012.9lf\n",local_t->tm_hour,local_t->tm_min, nano_usec);
        strftime (time_string, sizeof (time_string), "%d:%m:%Y", local_t);
        fprintf(FATM,"Date: %s\n",time_string);
        fclose(FATM);

}

/*******************************************************************
*FUNCTION: AquireData::readFromSHM()
*Reads from collect_psr shared memory of GSB/GWB
*******************************************************************/
int AquireData::readFromSHM()
{  	
	int DataOff=4096;
  	int bufferBlockLength=(dataBuffer->blocksize-DataOff)/info.sampleSizeBytes;  		
	long samplesToTake=info.noOfChannels*info.noOfPol*blockLength* info.sampleSizeBytes;	
	rawData=new unsigned  short int[blockLength*info.noOfPol*info.noOfChannels];	
	long int fetched=0;
	ofstream meanFile;
	meanFile.open("realTimeWarning.gpt",ios::app);
	ofstream recordFile;
	recordFile.open("blockLossRecord.gpt",ios::app);


////9 Jan 2022 Garvit Agarwal. Made several changes to 1) reduce the verbosity 2) clear the statements that are printed before the next stdout (for both console and log file) for compactness
	long int pos;
	pos=custom_cout::gbmon_log.tellp();

	while(fetched<samplesToTake)
	{
  		timeWaitTime+=omp_get_wtime();
		int flag=0;
		while((dataHdr->status == DAS_START) && (dataTab[recNum].flag &BufReady) == 0)
		{
			/*if(flag==0)
			{
				custom_cout::gbmon_log.seekp(pos);
				custom_cout::cout<<"\nWaiting for data block.."<<std::flush;
				flag=1;
			}*/
			usleep(2000);
		}
		/*if(flag==1){
			custom_cout::gbmon_log.seekp(pos);
			custom_cout::cout<<"Ready\n\n"<<std::flush;
		}*/		
		timeWaitTime-=omp_get_wtime();
		if(dataHdr->status != DAS_START)
		{
			/*if ((dataTab[recNum].flag & BufReady) == 0)
			{
				custom_cout::cout<<"DAS not in START mode!!"<<"\r"<<std::flush;
				return -1;
			}*/
			usleep(2000);
		}
		currentReadBlock = dataTab[recNum].seqnum;
		if(!info.isInline)
		{			
			/*custom_cout::cout<<endl<<"recNum "<<recNum<<endl;
			custom_cout::cout<<"dataBuffer->cur_rec "<<dataBuffer->cur_rec<<endl;
			custom_cout::cout<<"MaxDataBuf "<<MaxDataBuf<<endl;
			custom_cout::cout<<"dataBuffer->cur_block "<<dataBuffer->cur_block<<endl;
			custom_cout::cout<<"currentReadBlock "<<currentReadBlock<<endl;
			custom_cout::cout<<"dataBuffer->maxblocks "<<dataBuffer->maxblocks<<endl<<endl;
			*/

			custom_cout::gbmon_log.seekp(pos);
			custom_cout::cout<<"Current-scan block being	written: "<<dataBuffer->cur_block<<"	"<<"read: "<<currentReadBlock<<"\r"<<std::flush;		

		}
		//memcpy(&currentReadBlock, &dataTab[recNum].seqnum, sizeof(int));
		//dataBuffer->maxblocks-1
		if(dataBuffer->cur_block - currentReadBlock >=dataBuffer->maxblocks-1)
		{
			meanFile<<"Lag in block blockIndex="<<blockIndex<<endl;
			meanFile<<"recNum = "<<recNum<<", Reading Sequence: "<<currentReadBlock<<", Collect's Sequence: "<<(dataBuffer->cur_block-1)<<endl;
			
			//custom_cout::gbmon_log.seekp(pos);
			custom_cout::cout<<"Processing lagged behind..."<<std::flush;
			//recordFile<<blockIndex
			/*while(dataBuffer->cur_block - currentReadBlock >=dataBuffer->maxblocks-1)
			{
				if(samplesToTake-fetched>dataBuffer->blocksize-DataOff-remainingData)
				{
					memcpy(rawData+fetched/sizeof(unsigned short), zeros+remainingData, dataBuffer->blocksize-DataOff-remainingData);
  				fetched+=(dataBuffer->blocksize-DataOff-remainingData);
  				recNum=(recNum+1)%MaxDataBuf;
				remainingData=0;
  				}
  				else
  				{
					memcpy(rawData+fetched/sizeof(unsigned short), zeros+remainingData,samplesToTake-fetched);				
					remainingData=(samplesToTake-fetched);
					curPos+=samplesToTake;
					meanFile.close();
  					return -1;
  				}
				currentReadBlock = dataTab[(recNum)].seqnum;
			}*/

			//custom_cout::gbmon_log.seekp(pos);
			custom_cout::cout<<"recNum = "<<recNum<<", Reading Sequence: "<<currentReadBlock<<", Collect's Sequence: "<<(dataBuffer->cur_block-1)<<" blockIndex = "<<blockIndex<<"  Realiging..."<<std::flush;
			recNum = (dataBuffer->cur_rec-1-2+MaxDataBuf)%MaxDataBuf;
			currentReadBlock = dataTab[(recNum)].seqnum;
		}

		custom_cout::gbmon_log.seekp(pos);

		//custom_cout::cout<<"recNum = "<<recNum<<", Reading Sequence: "<<currentReadBlock<<", Collect's Sequence: "<<dataBuffer->cur_block-1<<endl;
		//custom_cout::cout<<"Blocksize="<<dataBuffer->blocksize<<endl;
		if(samplesToTake-fetched>dataBuffer->blocksize-DataOff-remainingData)
		{	
			if(info.isInline)
			{
				memcpy(headerInfo+nBuffTaken*DataOff, dataBuffer->buf+dataTab[recNum].rec*(dataBuffer->blocksize), DataOff);
				nBuffTaken++;
			}
  			memcpy(rawData+fetched/sizeof(short), dataBuffer->buf+dataTab[recNum].rec*(dataBuffer->blocksize)+DataOff+remainingData, dataBuffer->blocksize-DataOff-remainingData);
  			fetched+=(dataBuffer->blocksize-DataOff-remainingData);
  			recNum=(recNum+1)%MaxDataBuf;
			remainingData=0;
  		}
  		else
  		{
  			memcpy(rawData+fetched/sizeof(short), dataBuffer->buf+dataTab[recNum].rec*(dataBuffer->blocksize)+DataOff+remainingData, samplesToTake-fetched);
			remainingData=(samplesToTake-fetched);
  			break;
  		}
  	}
  	curPos+=samplesToTake;
	meanFile.close();
	return 1;
}
/*******************************************************************
*FUNCTION: void AquireData::readDataFromFile()
*In offline mode reads raw data from a file
*******************************************************************/
void AquireData::readDataFromFile()
{
	int 	 c,i;
	long int blockSizeBytes= info.noOfChannels*info.noOfPol* blockLength* info.sampleSizeBytes; //Number of bytes to read in eac block
							//Number of bytes that have already been read
	
	ifstream datafile;
	datafile.open(info.filepath,ios::binary);	
	datafile.seekg(curPos,ios::beg);
	//logic to handle reading last block
	if(curPos+blockSizeBytes> eof)
	{
		blockSizeBytes=eof-curPos;
		blockLength=blockSizeBytes/(info.sampleSizeBytes*info.noOfChannels*info.noOfPol);		//The number of time samples 
		hasReachedEof=1;
	}

	
	/*******************************************************************
	*Handles different data types. GMRT data is mostly of type short while 
	*certain processed data maybe floating point.
	*******************************************************************/
	switch(info.sampleSizeBytes)
	{
		case 1:
			if(!info.doPolarMode)
			{
				rawDataChar=new unsigned char[blockSizeBytes/info.sampleSizeBytes];
				datafile.read((char*)rawDataChar,blockSizeBytes);
			}	
			else
			{
				rawDataCharPolar=new char[blockSizeBytes/info.sampleSizeBytes];
				datafile.read((char*)rawDataCharPolar,blockSizeBytes);
			
			}
			break;
		case 2:
			if(!info.doPolarMode)
			{
				rawData=new unsigned short[blockSizeBytes/info.sampleSizeBytes];
				datafile.read((char*)rawData,blockSizeBytes);			
			}
			else
			{
				rawDataPolar=new short[blockSizeBytes/info.sampleSizeBytes];
				datafile.read((char*)rawDataPolar,blockSizeBytes);	
			}

			break;
		case 4:
			rawDataFloat=new float[blockSizeBytes/info.sampleSizeBytes];
			datafile.read((char*)rawDataFloat,blockSizeBytes);
			break;
	}
	datafile.close();		
	curPos+=blockSizeBytes;
}
float u16tofloat(short x)
{
    union { float f; int i; } u; u.f = 0.00f; u.i |= x;
    return u.f ; 
}
/**********************************************************************
*FUNCTION: void AquireData::splitRawData()
*If on polar mode, splitting of raw data into four polarization channels
*is done here. GMRT polarization data is read out in the following fomat:
T1_C1_P T1_C1_Q T1_C1_R T1_C1_S T1_C2_P T1_C2_Q T1_C2_R T1_C2_S...
T1_CN_S T2_C1_P ...........TN_CN_S
where P,Q,R,S are the four polarizations, Cn denotes the nth channel 
and Tn the nth time sample. So T15_C25_Q means the Q-polarization data
of channel 25 of the 15th time sample.
**********************************************************************/
void AquireData::splitRawData()
{
	
	splittedRawData=new float*[info.noOfPol];
	float **ptrSplittedRawData=new float*[info.noOfPol];
	long int length=blockLength*info.noOfChannels;	
	for(int k=0;k<info.noOfPol;k++)
	{
		splittedRawData[k]=new float[length];
		ptrSplittedRawData[k]=splittedRawData[k];
	}
	/*******************************************************************
	*Handles different data types. GMRT data is mostly of type short while 
	*certain processed data maybe floating point.
	*******************************************************************/
	
	float* ptrRawDataFloat=rawDataFloat;
	switch(info.sampleSizeBytes)
	{
		case 1:		
			if(!info.doPolarMode)
			{
				unsigned char* ptrRawData=rawDataChar;	
				for(int j=0;j<blockLength;j++)		//Refer to the GMRT polarization data format in Function description.
				{		
					for(int i=0;i<info.noOfChannels;i++,ptrRawData++,ptrSplittedRawData[0]++)
					{
							*(ptrSplittedRawData[0])=(*ptrRawData);
					}
				}
			}
			else
			{
				char* ptrRawData=rawDataCharPolar;	
				for(int j=0;j<blockLength;j++)		//Refer to the GMRT polarization data format in Function description.
				{		
					for(int i=0;i<info.noOfChannels;i++)
					{
							*(ptrSplittedRawData[0]++)=(*ptrRawData++);
							*(ptrSplittedRawData[1]++)=(*ptrRawData++);
							*(ptrSplittedRawData[2]++)=(*ptrRawData++);
							*(ptrSplittedRawData[3]++)=(*ptrRawData++);

					}
				}
			}
			break;
		case 2:		
			if(!info.doPolarMode)
			{
				unsigned short int* ptrRawData=rawData;	
				for(int j=0;j<blockLength;j++)		//Refer to the GMRT polarization data format in Function description.
				{		
					for(int i=0;i<info.noOfChannels;i++,ptrRawData++,ptrSplittedRawData[0]++)
					{
							*(ptrSplittedRawData[0])=(*ptrRawData);
					}
				}
			}
			else
			{
				short int* ptrRawData=rawDataPolar;	
				for(int j=0;j<blockLength;j++)		//Refer to the GMRT polarization data format in Function description.
				{		
					for(int i=0;i<info.noOfChannels;i++)
					{
							*(ptrSplittedRawData[0]++)=(*ptrRawData++);
							*(ptrSplittedRawData[1]++)=(*ptrRawData++);
							*(ptrSplittedRawData[2]++)=(*ptrRawData++);
							*(ptrSplittedRawData[3]++)=(*ptrRawData++);

					}
				}
			}
			break;
		case 4:					
			for(int j=0;j<blockLength;j++)		//Refer to the GMRT polarization data format in Function description.
			{		
				for(int i=0;i<info.noOfChannels;i++)
				{
					for(int k=0;k<info.noOfPol;k++,ptrRawDataFloat++)
					{				
						(*ptrSplittedRawData[k])=(*ptrRawDataFloat);	
						ptrSplittedRawData[k]++;
					}	
				}
			}
			break;
	}
	delete[] ptrSplittedRawData;
}
//End of AquireData implementation.
