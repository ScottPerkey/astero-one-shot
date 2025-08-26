#!/bin/bash
# copies over the relevant info from tess --- including history.data friles and ad_approx_rad_err.dat files
# creates a file called /Users/jzinn/nonad/src/nonad/data/count.index, which translates the history_##.data and ad_approx_rad_err_##.dat ## into a particular point in the grid (mass, feh, alpha)

set -e
masses=( 1.0 1.2 )
# masses=( 1.4 1.8 2.2 )
#0.8 0.9 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.3 2.4 2.6 2.8 3.0 3.2 3.4 )
fehs=( -1.0d0  -4.0d-1  -2.0d-1  0.0d0  2.0d-1  4.0d-1 )
zs=( 1.7d-3 6.8d-3  1.1d-2  1.7d-2  2.7d-2  4.3d-2 )
alphas=( 1.7 2.0 )
count=0 # use count = 0 for the 1.0 and 1.2 masses, but count=24 for 1.4 1.8 and 2.2.
cd /Users/jzinn/nonad/src/nonad/data/edd
if [[ count == 0 ]]; then
   echo "count mass feh alpha" > count.index
fi
   
for mass in ${masses[*]}; do
    for i in ${!fehs[@]}; do
	for j in ${!alphas[@]}; do
	    echo $mass
	    ~/t2m.sh "/home/zinn/mesa/nonad/nonad_${mass}msun_${fehs[$i]}feh_mass_loss_${alphas[$j]}alpha/LOGS//history.data" history_${count}.data
	    ~/t2m.sh "/home/zinn/mesa/nonad/nonad_${mass}msun_${fehs[$i]}feh_mass_loss_${alphas[$j]}alpha/LOGS//profiles.index" profiles_${count}.index
	    ~/t2m.sh "/home/zinn/mesa/nonad/nonad_${mass}msun_${fehs[$i]}feh_mass_loss_${alphas[$j]}alpha/LOGS/edd/ad_approx_rad_err.dat" ad_approx_rad_err_${count}.dat
	    
	    echo $count ${mass} ${fehs[$i]} ${alphas[$j]} >> count.index
	    count=`echo $count + 1 | bc -l`
	    echo $count

	done
    done
done

echo "copy_nonad.sh DONE"
    
