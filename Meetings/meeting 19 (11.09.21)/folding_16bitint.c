
	float mean=0.0,rms=0.0;
	for(i=0; i<N ;i++)
	{
		fread(&val, sizeof(float), 1, f);
		time_series[i]=val;

		mean+=val;
	}
	fclose(f);

	mean/=N;	
	//calculating data rms	
	for(i=0; i<N ;i++)rms+=(mean-time_series[i])*(mean-time_series[i]);
	rms/=N;
	rms=sqrtf(rms);

	//defining a range for the data and changing outlier data points to border values
	int n=3;
	float float_range[2]={mean-n*rms,mean+n*rms};
	for(i=0; i<N ;i++)
	{
		if(time_series[i]>float_range[1]) time_series[i]=float_range[1];
		else if(time_series[i]<float_range[0]) time_series[i]=float_range[0];
	}

	//finally calculating the converted values and reconverting them to float
	short result[N];
	for(i=0; i<N ;i++) 
	{
		result[i]=(short)roundf((time_series[i]-float_range[0])/(2*n*rms) * 65536.0 - 32768.0);
		time_series[i]=(result[i]/32768.0)*n*rms+mean;
	}


