/*************************************************************
shm library
This library is aimed at interfacing between various
kinds of shared memory structures used by different
pulsar codes at the GMRT.
				-Aditya Chowdhury, 14th March 2016
************************************************************/
#include "SHM.h"

#include<custom_cout.h> //7th Jan 2022 Garvit Agarwal. to use the custom cout for gbmon log

using namespace std;



DasHdrType*	Correlator::dataHdrWrite;
DataBufType*	Correlator::dataBufferWrite;
DataTabType*	Correlator::dataTabWrite;

DasHdrType*	Correlator::dataHdrRead;
DataBufType*	Correlator::dataBufferRead;
DataTabType*	Correlator::dataTabRead;

int		Correlator::recNumRead=0;
int		Correlator::recNumWrite=0;
long int	Correlator::currentReadBlock=0;


Correlator::Correlator(DasHdrType *_dataHdrRead ,DataBufType *_dataBufferRead)
{
	DataOff=4096;
	debug=1;
	dataHdrRead=_dataHdrRead;
	dataBufferRead=_dataBufferRead;
}
Correlator::Correlator(int _nchan,float _sampling)
{
	DataOff=4096;
	debug=1;
	nchan=_nchan;
	sampling=_sampling;
	dataHdrRead=NULL;
	dataBufferRead=NULL;
}
/*******************************************************************
*FUNCTION: void AquireData::initializeSHM()
*Initializes collect_psr shared memory of GSB/GWB
*******************************************************************/
void Correlator::initializeReadSHM()
{
  	int iddataHdr,idDataBuffer;
	iddataHdr = shmget(DAS_H_KEY,sizeof(DasHdrType),644);
	idDataBuffer = shmget(DAS_D_KEY,sizeof(DataBufType),0);
	custom_cout::cout<<"iddataHdr="<<iddataHdr<<endl<<"idDataBuffer="<<idDataBuffer<<endl;
	
	if(iddataHdr < 0 || iddataHdr < 0)
	{
		exit(1);
	}
	dataHdrRead = (DasHdrType*)shmat(iddataHdr,0,SHM_RDONLY);
  	dataBufferRead = (DataBufType*)shmat(idDataBuffer,0,SHM_RDONLY);
	if((dataBufferRead)==(DataBufType*)-1)
	{
		custom_cout::cout<<"Cannot attach to BUFFER1 shared memory!"<<endl;
		exit(1);
	}
	custom_cout::cout<<"Attached to shared memory:"<<dataHdrRead<<","<<dataBufferRead<<endl;
	dataTabRead = dataBufferRead-> dtab;
	custom_cout::cout<<"Max no of blocks="<<dataBufferRead->maxblocks<<endl;

	 /*   find a block, two blocks before the current block of the shm for reading data */
	//if(dataBuffer->cur_rec > (dataBuffer->maxblocks)/2)
	recNumRead = (dataBufferRead->cur_rec-2+MaxDataBuf)%MaxDataBuf;
	currentReadBlock = dataTabRead[recNumRead].seqnum;  	
}
void Correlator::copyHeaderInfo()
{
	dataHdrWrite->active=dataHdrRead->active;
	dataHdrWrite->status=dataHdrRead->status;
	dataHdrWrite->scan=dataHdrRead->scan;
	dataHdrWrite->scan_off=dataHdrRead->scan_off;
	dataHdrWrite->corr=dataHdrRead->corr;
	dataHdrWrite->model=dataHdrRead->model;
	dataHdrWrite->BeamHeader=dataHdrRead->BeamHeader;
 	dataBufferWrite->blocksize=dataBufferRead->blocksize;
	dataBufferWrite->maxblocks=dataBufferRead->maxblocks;
}
/*******************************************************************
*FUNCTION: void AquireData::initializeSHM()
*Initializes collect_psr shared memory of GSB/GWB
*******************************************************************/
int Correlator::initializeWriteSHM()
{
  	int iddataHdr,idDataBuffer;
	if(dataHdrRead==NULL)
	{
		iddataHdr = shmget(DAS_H_KEY_GPTOOL,sizeof(DasHdrType),IPC_CREAT | 0666 );
		idDataBuffer = shmget(DAS_D_KEY_GPTOOL,sizeof(DataBufType),IPC_CREAT | 0666 );
	}
	else
	{
		iddataHdr = shmget(DAS_H_KEY_GPTOOL_INLINE,sizeof(DasHdrType),IPC_CREAT | 0666 );
		idDataBuffer = shmget(DAS_D_KEY_GPTOOL_INLINE,sizeof(DataBufType),IPC_CREAT | 0666 );
	}
	custom_cout::cout<<"iddataHdr="<<iddataHdr<<endl<<"idDataBuffer="<<idDataBuffer<<endl;
	
	if(iddataHdr < 0 || iddataHdr < 0)
	{
		custom_cout::cout<<"Error creating shared memory"<<endl;
		exit(1);
	}
	dataHdrWrite = (DasHdrType*)shmat(iddataHdr,0,0);
  	dataBufferWrite = (DataBufType*)shmat(idDataBuffer,0,0);
	if((dataBufferWrite)==(DataBufType*)-1)
	{
		custom_cout::cout<<"Cannot attach to BUFFER2 shared memory!"<<endl;
		exit(1);
	}
	custom_cout::cout<<"Attached to write shared memory:"<<dataHdrWrite<<","<<dataBufferWrite<<endl;

	if(dataHdrRead==NULL)
	{	
		dataBufferWrite->blocksize=2*nchan*int(BLOCKTIME/sampling)+DataOff;
		dataBufferWrite->maxblocks=int(DAS_BUFSIZE/(dataBufferWrite->blocksize));
	}
	else
	{
		copyHeaderInfo();
		custom_cout::cout<<"header info copied"<<endl;	
	}
	dataTabWrite = dataBufferWrite-> dtab;
	dataBufferWrite->cur_rec=0;
	dataBufferWrite->cur_block=0;
	recNumWrite = (dataBufferWrite->cur_rec)%MaxDataBuf;
	dataTabWrite[recNumWrite].seqnum=0;
	dataTabWrite[recNumWrite].rec=0;
	for(int i=0;i<MaxDataBuf;i++)
		dataTabWrite[i].flag=0;
	custom_cout::cout<<"Max no of blocks="<<dataBufferWrite->maxblocks<<endl;
	dataHdrWrite->status = DAS_START;
	return dataBufferWrite->blocksize-DataOff;
}
void Correlator::writeToSHM(unsigned short int* rawData)
{  	
 	long int fetched=0;
	ofstream warnFile;
	dataTabWrite[recNumWrite].seqnum=dataBufferWrite->cur_block;
	dataTabWrite[recNumWrite].rec=(dataBufferWrite->cur_block)%(dataBufferWrite->maxblocks);
	
	
	if(debug)
	{
		custom_cout::cout<<endl<<"recNum "<<recNumWrite<<endl;
		custom_cout::cout<<"dataBuffer->cur_rec "<<dataBufferWrite->cur_rec<<endl;
		custom_cout::cout<<"MaxDataBuf "<<MaxDataBuf<<endl;
		custom_cout::cout<<"dataBuffer->cur_block "<<dataBufferWrite->cur_block<<endl;
		custom_cout::cout<<"dataBuffer->maxblocks "<<dataBufferWrite->maxblocks<<endl;
		custom_cout::cout<<"dataTabWrite[recNumWrite].rec "<<dataTabWrite[recNumWrite].rec<<endl;
		custom_cout::cout<<"dataTabWrite[recNumWrite].flag "<<dataTabWrite[recNumWrite].flag<<endl;
		custom_cout::cout<<"dataTabWrite[recNumWrite].seqnum"<<dataTabWrite[recNumWrite].seqnum<<endl;
		custom_cout::cout<<"dataHdrWrite->status "<<dataHdrWrite->status<<endl;
		custom_cout::cout<<dataBufferWrite->blocksize<<endl;
	}	
		
	memcpy(dataBufferWrite->buf+dataTabWrite[recNumWrite].rec*(dataBufferWrite->blocksize)+DataOff,rawData, dataBufferWrite->blocksize-DataOff);
	
  		
	dataBufferWrite->cur_rec=(recNumWrite+1)%MaxDataBuf;
	dataBufferWrite->cur_block+=1;
	warnFile.close();
	dataTabWrite[(recNumWrite+1)%MaxDataBuf].flag=0;
	dataTabWrite[(recNumWrite)].flag=BufReady;
	recNumWrite=(recNumWrite+1)%MaxDataBuf;  

}
void Correlator::writeToSHM(unsigned short int* rawData,char* header)
{  	
	dataTabWrite[recNumWrite].seqnum=dataBufferWrite->cur_block;
	dataTabWrite[recNumWrite].rec=(dataBufferWrite->cur_block)%(dataBufferWrite->maxblocks);
	if(debug)
	{
		custom_cout::cout<<endl<<"recNum "<<recNumWrite<<endl;
		custom_cout::cout<<"dataBuffer->cur_rec "<<dataBufferWrite->cur_rec<<endl;
		custom_cout::cout<<"MaxDataBuf "<<MaxDataBuf<<endl;
		custom_cout::cout<<"dataBuffer->cur_block "<<dataBufferWrite->cur_block<<endl;
		custom_cout::cout<<"dataBuffer->maxblocks "<<dataBufferWrite->maxblocks<<endl;
		custom_cout::cout<<"dataTabWrite[recNumWrite].rec "<<dataTabWrite[recNumWrite].rec<<endl;
		custom_cout::cout<<"dataTabWrite[recNumWrite].flag "<<dataTabWrite[recNumWrite].flag<<endl;
		custom_cout::cout<<"dataTabWrite[recNumWrite].seqnum"<<dataTabWrite[recNumWrite].seqnum<<endl;
		custom_cout::cout<<"dataHdrWrite->status "<<dataHdrWrite->status<<endl;
		custom_cout::cout<<dataBufferWrite->blocksize<<endl;
	}	
	custom_cout::cout<<"start memcopy rawdata"<<endl;	
	memcpy(dataBufferWrite->buf+dataTabWrite[recNumWrite].rec*(dataBufferWrite->blocksize)+DataOff,rawData, dataBufferWrite->blocksize-DataOff);
	custom_cout::cout<<"memcpy done"<<endl;
	custom_cout::cout<<"start memcopy header"<<endl;
	memcpy(dataBufferWrite->buf+dataTabWrite[recNumWrite].rec*(dataBufferWrite->blocksize),header,DataOff);	
 	custom_cout::cout<<"memcpy done"<<endl; 		
	dataBufferWrite->cur_rec=(recNumWrite+1)%MaxDataBuf;
	dataBufferWrite->cur_block+=1;
	dataTabWrite[(recNumWrite+1)%MaxDataBuf].flag=0;
	dataTabWrite[(recNumWrite)].flag=BufReady;
	recNumWrite=(recNumWrite+1)%MaxDataBuf;  

}
/*******************************************************************
*FUNCTION: AquireData::readFromSHM()
*Reads from collect_psr shared memory of GSB/GWB
*******************************************************************/
void Correlator::readFromSHM(unsigned short int* rawData)
{  	

	ofstream warnFile;
	warnFile.open("realTimeWarning.gpt",ios::app);
	int flag=0;
	while((dataHdrRead->status == DAS_START) && (dataTabRead[recNumRead].flag &BufReady) == 0)
	{
		usleep(2000);
		if(flag==0)
		{
			custom_cout::cout<<"Waiting"<<endl;
			flag=1;
		}
	}
	if(flag==1)
		custom_cout::cout<<"Ready"<<endl;
		
	if(dataHdrRead->status != DAS_START)
	{
		if ((dataTabRead[recNumRead].flag & BufReady) == 0)
		{
			custom_cout::cout<<"DAS not in START mode!!"<<endl;
			exit(0);
		}
	}
	currentReadBlock = dataTabRead[recNumRead].seqnum;
	if(debug)
	{
		custom_cout::cout<<endl<<"recNum "<<recNumRead<<endl;
		custom_cout::cout<<"dataBuffer->cur_rec "<<dataBufferRead->cur_rec<<endl;
		custom_cout::cout<<"MaxDataBuf "<<MaxDataBuf<<endl;
		custom_cout::cout<<"dataBuffer->cur_block "<<dataBufferRead->cur_block<<endl;
		custom_cout::cout<<"currentReadBlock "<<currentReadBlock<<endl;
		custom_cout::cout<<"dataBuffer->maxblocks "<<dataBufferRead->maxblocks<<endl<<endl;
	}
		
	if(dataBufferRead->cur_block - currentReadBlock >=dataBufferRead->maxblocks-1)
	{
		warnFile<<"recNum = "<<recNumRead<<", Reading Sequence: "<<currentReadBlock<<", Collect's Sequence: "<<(dataBufferRead->cur_block-1)<<endl;
		warnFile<<"Processing lagged behind..."<<endl;
			
		custom_cout::cout<<"recNum = "<<recNumRead<<", Reading Sequence: "<<currentReadBlock<<", Collect's Sequence: "<<(dataBufferRead->cur_block-1)<<"\nRealiging...\n";
		recNumRead = (dataBufferRead->cur_rec-1-2+MaxDataBuf)%MaxDataBuf;
		currentReadBlock = dataTabRead[(recNumRead)].seqnum;
	}
		
	memcpy(rawData, dataBufferRead->buf+dataTabRead[recNumRead].rec*(dataBufferRead->blocksize)+DataOff, dataBufferRead->blocksize-DataOff);
  	recNumRead=(recNumRead+1)%MaxDataBuf;
  	
	warnFile.close();
}
