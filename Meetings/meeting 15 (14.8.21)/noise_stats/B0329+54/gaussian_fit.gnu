
#values to change as required: f in this case is a command line argument
binwidth= .000075 #.1
pulsar= "B0329+54"

#normal histogram routine----------------------------------------
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
n=round(f*4361)

bin(x,width) = width*floor(x/width) + width/2.0
set boxwidth binwidth* 1.0 #0.7
set style fill solid 1.0

set title sprintf("%s: Noise energy distribution (about the mean)\n no. of noise values summed - %d",pulsar,n)
set xlabel "noise values"
#----------------------------------------------------------------

#makes data file for histogram (later deleted automatically) to fit the gaussian on
set table 'hist.temp'
plot "noise_energy.txt" u (bin($2, binwidth)):(1.0) smooth freq with boxes
unset table

#the fit commands
gauss(x)=a/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))
mu=0.00001
a=.03
sigma=0.0003
fit gauss(x) "hist.temp" u 1:2 via a, sigma, mu

#plots both the histogram and the fitted gaussian function together
set terminal png size 1000,600
set output sprintf("../noise_stats/gaussianfit_%d.png",n)
set style line 2 lw 100
plot "noise_energy.txt" u (bin($2, binwidth)):(1.0) smooth freq with boxes lc rgb "red", gauss(x) w lp lw 2 lc rgb "green"





