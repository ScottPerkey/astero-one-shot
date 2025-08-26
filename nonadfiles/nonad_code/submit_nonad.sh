#!/bin/bash
# will run to the tams for the paxton rc inlist for the case including mass loss
masses=( 1.0 1.2 )
#0.8 0.9 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.3 2.4 2.6 2.8 3.0 3.2 3.4 )
fehs=( -1.0d0  -4.0d-1  -2.0d-1  0.0d0  2.0d-1  4.0d-1 )
zs=( 1.7d-3 6.8d-3  1.1d-2  1.7d-2  2.7d-2  4.3d-2 )
alphas=( 1.7 2.0 )

# these have been sent to run
fehs=( 0.0d0 )
zs=( 1.7d-2 )
masses=( 1.2 )

# these have been sent to run
fehs=( -1.0d0  -4.0d-1  -2.0d-1 2.0d-1  4.0d-1 )
zs=( 1.7d-3 6.8d-3  1.1d-2 2.7d-2  4.3d-2 )
masses=( 1.2 )


# these have been sent to run
fehs=( -1.0d0  -4.0d-1 0.0d0 -2.0d-1 2.0d-1  4.0d-1 )
zs=( 1.7d-3 6.8d-3 1.7d-2 1.1d-2 2.7d-2  4.3d-2 )
masses=( 1.0 )

# sent to run
fehs=( 4.0d-1 )
zs=( 4.3d-2 )
masses=( 1.2 )

# need to rerun this because run into timestep problem.... also an abundance problem.
fehs=( -1.0d0 )
zs=( 1.7d-3 )
masses=( 1.2 )
alphas=( 2.0 )

# need to rerun this because run into timestep problem.... also an abundance problem.
fehs=( -4.0d-1 )
zs=( 6.8d-3 )
masses=( 1.2 )
alphas=( 1.7 )

for mass in ${masses[*]}; do
    for i in ${!fehs[@]}; do
	for j in ${!alphas[@]}; do
	    cd /home/zinn/mesa/nonad/
	    mkdir nonad_${mass}msun_${fehs[$i]}feh_mass_loss_${alphas[$j]}alpha
	    cd nonad_${mass}msun_${fehs[$i]}feh_mass_loss_${alphas[$j]}alpha
	    pwd
	    cp -f /home/zinn/mesa/nonad/inlists/* ./
	    # for mass loss, i've already made the correct choices for mass loss elsewhere so will override the normal
	    # inlist_common_nogold and inlist_common
	    
	    cp -rf $MESA_DIR/star/work/rn ./
	    cp -rf $MESA_DIR/star/work/src ./
	    cp -rf $MESA_DIR/star/work/mk ./
	    cp -rf $MESA_DIR/star/work/make ./
	    
	    # JCZ 221119
	    # change mass to desired mass
	    sed -i "s|initial_mass = 1.0|initial_mass = ${mass}|" inlist_base_lmesa_gyre
	    sed -i "s|Zbase = 0.02|Zbase = ${zs[$i]}|" inlist_base_lmesa_gyre
	    sed -i "s|initial_z = 0.02d0|initial_z = ${zs[$i]}|" inlist_base_lmesa_gyre
	    sed -i "s|max_years_for_timestep = 1d7| max_years_for_timestep = 1d6|" inlist_base_lmesa_gyre
	    sed -i "s|! mixing_length_alpha = 2.0| mixing_length_alpha = ${alphas[$j]}|" inlist_base_lmesa_gyre
	    sed -i "s|create_pre_main_sequence_model = .false.| create_pre_main_sequence_model = .true.|" inlist_base_lmesa_gyre
	    nohup bash -c "./mk && nice -n 19 ./rn" &
	done
    done
done

    
