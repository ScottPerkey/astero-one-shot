import pandas as pd
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import string
import shutil
import csv

from pysyd import plots
from pysyd.target import Target
from pysyd.utils import Parameters
from pysyd import PACKAGEDIR
from pysyd.plots import plot_1d_ed
from pysyd import utils


class pySYD_oneshot:
	def __init__(self,**kwargs):
		MPLSTYLE = os.path.join(PACKAGEDIR,'data','pysyd.mplstyle')
		plt.style.use(MPLSTYLE)
		#lets use this to call the values of the directories cause that shit is annoying to creat it each time
		#
		#			
		#
		#		Returns/Attributes
		#		
		#		"""
	

	def power_spec_filename_correction(self):
		#we need to do this other wise pySYD cannot read the values in the file correctly
		base_directory='/home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha'
		for dir_name in os.listdir(base_directory):
			dir_path = os.path.join(base_directory, dir_name)
			#make sure that the directory even exists
			if os.path.isdir(dir_path):
				string_number=dir_name.replace('profile','').strip()
				taken_number=string_number.split('.')[0]
				csv_file_path=os.path.join(dir_path,f'profile{taken_number}freq_amp.csv')
				# we do below because some of them are not created and/or are unfilled data
				if os.path.isfile(csv_file_path):
					#create the new file name that pySYD needs for the expected power spectrum

					pySYD_newname=f'{taken_number}_PS.txt'
					pySYD_newpath=os.path.join(dir_path,pySYD_newname)

					with open(csv_file_path,'r') as csv_file:
						read_original_txt=csv_file.read().replace(',',' ')
					with open(pySYD_newpath,'w') as new_file:
						new_file.write(read_original_txt)
					print(f'Correctly changed the name of the old literature created data of {pySYD_newpath} for pySYD')
		#this makes the 100_ps_txt format, it has to happen for the power spectrum because pysyd usually starts with light curve but we already have the power spectrum


	def process_txt_file(self,file_path):
        	data=[]
        	with open(file_path,'r') as file:
        	        for line in file:
        	                split_line=line.split()
        	                col1,col2 = map(float,line.split())
        	                data.append((col1,col2))
        	return max(data,key=lambda x: x[1])[0]


	def echelle_plot(self,name,star):
		fig, ax=plt.subplots(dpi=300)
		try:
			params = star.params['plotting']['estimates']
			use_dnu = params['use_dnu']
			interpolation='nearest'
			ax.imshow(params['ed'], extent=params['extent'], interpolation=interpolation, aspect='auto', origin='lower', cmap=plt.get_cmap(star.params['cmap']))
			ax.axvline([use_dnu], color='white', linestyle='--', linewidth=1.5, dashes=(5, 5))
			ax.set_title(r'$\rm \grave{E}chelle \,\, diagram$')
			ax.set_xlabel(r'$\rm \nu \,\, mod \,\, %.2f \,\, [\mu Hz]$' % use_dnu)
			ax.set_ylabel(r'$\rm \nu \,\, [\mu Hz]$')
			ax.set_xlim([params['extent'][0], params['extent'][1]])
			ax.set_ylim([params['extent'][2], params['extent'][3]])
			#probably change this later sorta
			plt.savefig(f'/home/032272043/projects/dnuGridCode/temp_move/ed_pysyd_oneshot/ed_{name}.png', dpi=300)
		except Exception as e:
			print(f"profile{name} cant ED : {e} for who cares cause frequency and echelle diagram done messed up I guess")
                        #plt.clf()       
                        #plot_1d_ed(star,filename=f'/home/032272043/projects/dnuGridCode/collapse_ed_pysyd/profile{name}_col_ed.png')


	def full_pysyd(self):
		os.mkdir('/home/032272043/projects/dnuGridCode/temp_move/ed_pysyd_oneshot')
		directory_path ='/home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha'
		file_list=os.listdir(directory_path)
		filtered_files=[file for file in file_list if file.endswith('.data.GYRE')]
		numbers = [int(re.search(r'\d+', file).group()) for file in filtered_files]


		csv_file='1M0_0D01.7ALPHA_SYDDNU_test4one_shot.csv'
		file_exists=os.path.isfile(csv_file)
		if not file_exists:
                	with open(csv_file, 'w') as file:

                        	writer=csv.writer(file)
                        	header=['ID','dnu_syd']
                        	writer.writerow(header)

		with open(csv_file,mode='a',newline='') as file:
			writer=csv.writer(file)
			for name in numbers:

				params=Parameters()
				params.add_targets(stars=name)
                        	#changed below usually 113 is name
	                       	#params.params[113]['show'], params.params[113]['verbose']=True,True 
				params.params[name]['show'],params.params[name]['verbose']=True,True

				star=Target(name, params)
                        	# temporary below just do 113 cause we know what its finna gonna look like lowkey
				star.params['inpdir']=f'/home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha/profile{name}.data.GYRE/'
				star.params['outdir']='/home/032272043/projects/dnuGridCode/temp_move'
                        	# temporary below just do 113 cause we know what its finna gonna look like lowkey
				file_path=f'/home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha/profile{name}.data.GYRE/{name}_PS.txt'
				col1_for_max_col2=self.process_txt_file(file_path)
				star.params['obs_numax']=col1_for_max_col2
				star.params['exp_dnu']=utils.delta_nu(star.params['obs_numax'])
				star.params['lower_ps']=star.params['obs_numax']-3*star.params['exp_dnu']
				star.params['upper_ps']=star.params['obs_numax']+3*star.params['exp_dnu']
				star.params['numax_smoo']=star.params['obs_numax']
				star.params['ps_mask']=([star.params["lower_ps"],star.params["upper_ps"]])
				star.load_data()
				star.params['results']["estimates"]={  'dnu':[]    }
				star.params['plotting']={'estimates':{  'obs_dnu':[] , 'zoom_freq':[]  ,'parameters':[] }}
				star.module='estimates'
                        	#remember that I changed the path for the saving of both the plot_1d_ed, and the regular echelle diagram
                        	#/research/CNSM-JZinn/anaconda3/lib/python3.11/site-packages/pysyd/target.py    
				star.bg_corr=star.power
				star.params['smooth_width']=0.005
				star.params['numax']=col1_for_max_col2
				star.params['resolution']=star.frequency[1]-star.frequency[0]
				star.i = 0
				star.params['oversampling_factor']=5
				star.params['smooth_ps']=0.0
				star.compute_acf()
				#try:
				#remember that we need to change if you update python to: 

				#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
				#from scipy.optimize import make_strictly_feasible
				#x0 = make_strictly_feasible(x0, lb, ub)
				#this is because the lowest value of dnu is not inbetween lb, ub that is calcualted for the frequency space
				#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

				#try:
				star.frequency_spacing()
				#except Exception as e:
				#	print(f"profile{name} doesnt work because:  {e} for frequency spacing")
				try:
					star.echelle_diagram()
				except Exception as e:
					print(f"profile{name} cant plot this in echelle because: {e} for echelle diragram")


				best_dnu=star.params['best_lag']
				#print(best_dnu)
				print(f'This was the best_dnu found by pySYD for profile{name}: {best_dnu}')
				row=[name,best_dnu]
				writer.writerow(row)
				self.echelle_plot(name,star)

	def delete_repeats(self):
		df=pd.read_csv('1M0_0D01.7ALPHA_SYDDNU_test4one_shot.csv')
		df_no_duplicates=df.drop_duplicates(keep=False)
		df_no_duplicates.to_csv('1M0_0D01.7ALPHA_SYDDNU_cleaned.csv',index=False)
		

	def sync_pysyd_csv(self):
		#remember for all explanations for this section
		column_names_added=['dnu_inf','nu_max','mass','radius','luminosity','mixing_length','metalicity']
		added_radius_tofile=pd.read_csv('1M0_0D01.7ALPHA_SYDDNU_cleaned.csv')
		csv_toextract_radius=pd.read_csv('old_lit_nonad_1M0_0D01.7ALPHA.csv')
		csv_toextract_radius['ID']=csv_toextract_radius['PROFID'].str.replace('profile','').astype(int)
		added_radius_tofile['ID']=added_radius_tofile['ID'].astype(int)
		for x in column_names_added:
			profile_radius_map=dict(zip(csv_toextract_radius['ID'],csv_toextract_radius[x]))
			added_radius_tofile[x]=added_radius_tofile['ID'].map(profile_radius_map)
		added_radius_tofile.to_csv('1M0_0D01.7ALPHA_SYDDNU_added_synced_cols.csv',index=False)
	
	# We should do everything at first for the old lit. Like make all of the same columns in pysyd aswell I think that makes the most sense.

# check this again but I am pretty sure this is everything that you need 
#gawk = pySYD_oneshot()
#gawk.power_spec_filename_correction()
#gawk.full_pysyd()
#gawk.echelle_plot()
#gawk.delete_repeats()
#gawk.sync_pysyd_csv()


