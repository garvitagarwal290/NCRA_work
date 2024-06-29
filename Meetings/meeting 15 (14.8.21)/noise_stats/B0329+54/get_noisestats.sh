
#this is to be run inside the pulse_energy directory of the pulsar

gcc ../../C_code/pulse_energy/pulse_energy.c -o ../../C_code/pulse_energy/pulse_energy -lm -std=gnu99 

f=`echo .1 | bc`
for i in {1..11..1}
do

	./../../C_code/pulse_energy/pulse_energy ../stacked_profiles/pulses.bin -p 1 -n 2 -f $f
	
	gnuplot -e "f=`echo $f` ;load '../noise_stats/gaussian_fit.gnu';set print '../noise_stats/fitstats.txt' append ;print sprintf(\"%d %f %f %f\",n,mu,abs(sigma),FIT_STDFIT**2)"

	rm hist.temp fit.log	

	f=`echo $f+.07 | bc`
done


