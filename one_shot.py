import os
import re
import numpy as np
import matplotlib.pyplot as plt
import csv
import subprocess
import math
import glob
import random as rand
import pandas as pd
# import pySYD_one_shot
from scipy.interpolate import interp1d
from astropy.io import ascii
from scipy.misc import derivative
from matplotlib.pyplot import figure
from nonadfiles.nonad_code.fileio import read_gyre
from scipy.stats import chi2
from nonadfiles.nonad_code.plot_modes import dnu_theory
from nonadfiles.nonad_code.plot_modes import return_modes
from nonadfiles.nonad_code.scaling_relation import logg_MR, loggteff2numax, numax, dnu_MR, freq_dyn, numax_sun, dnu_sun, L_sun, R_sun, M_sun
from nonadfiles.nonad_code.dnu_util import lorentzian




class full_create:
	def __init__(self,**kwargs):
		self.wowee=5
#		"""docu string meme
#		Parameters
#
#			
#
#		Returns/Attributes
#		
#		"""

#TRY THIS FOR THE SELF VALUSE 
 		
 		#get fdnu we calculaute already local from this code as of right now 2024-2-15, dnu_theory is the delta nu normal, the delta nu infinite is the mass an    d radius from return modes.
 		# bam gets us the delta nu local from the fake power specturm whih is the l-0 modes, which is the new, better way of doing it I think it takes the auto     correlation function of the power spectrum at numax (local), 1?  
 		
	def process(self,ID,read_path):
		font = {'family' : 'DejaVu Sans', 'weight' : 'normal','size'   : 16}
		plt.rc('font', **font)
		t = 4*math.pi*10
		nm = 40
                # based on what we did before gamma should be 0.12 so width should be changed to 0.0036, dont know where the original width = nm *.005 comes from? 
                # based on this : https://iopscience.iop.org/article/10.3847/1538-4357/835/2/172/pdf
		width = 0.0036
                #width = nm*.005
		dF = 1/(2*t)
	
		parts=read_path.split('/')
		relevant_path=parts[4]
		full_path=read_path+'/'
		df = return_modes(full_path+ID, 'nad', visibility_thresh=0.0,).reset_index(drop=True)
		df["freq"] = pd.to_numeric(df["freq"], downcast="float")
		df["ratio"] = pd.to_numeric(df["ratio"], downcast="float")
		df0 = df[df.l == 0]
		df1 = df[df.l == 1]
		df2 = df[df.l == 2]

		#multiply by delta nu sun, microhertz multiply by right side  
		numax = (df['M_star']/M_sun) * (df['R_star']/R_sun) ** (-7./4.) * (df['L_star']/L_sun) **(-1./8.) * numax_sun
		dnutheory = dnu_theory(df["l"],df["freq"],numax.loc[0])
		#dnutheory is the delta nu local, as the current literature does it, the weighted sum of the l=0 modes
		#below is the deltanu infinite
		dnuinf=math.sqrt(    (df['M_star'].iloc[0]/M_sun)/((df['R_star'].iloc[0]/R_sun)**3) )*dnu_sun

		fdnuStdLit= dnuinf/dnutheory
		#below is fdnu which is the sqrt of the density (infinite delta nu) divided by the dnu local (as calculated by the weights of the power spectru    m around the nu local), or how the curren literature does it.  fdnuStdLit= dnuinf/dnutheory

		#f1=[]
		#f2 = []
		#for i in np.arange(len(df0['freq'])-1):
		#	try: 
		#		flow= df0['freq'].iloc[i]
		#		fhi = df0['freq'].iloc[i + 1]
		#		sub1 = df1.loc[(df1['freq'] > flow) & (df1['freq'] < fhi)]
		#		fl1 = sub1.iloc[np.argmax(sub1['ratio'])]['freq']
		#		f1.append(fl1)
		#		print('AHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH')
		#	except ValueError: 
		#		if not input: 
		#			raise ValueError('no data ')
		#	else:           
		#		continue
		#f1=pd.to_numeric(f1,downcast="float")
		#print(f1)
		for i in np.arange(len(df0['freq'])-1):
			flow= df0['freq'].iloc[i]
			fhi = df0['freq'].iloc[i + 1]
			sub1 = df1.loc[(df1['freq'] > flow) & (df1['freq'] < fhi)]
			#print(sub1)
			#jfl1 = sub1.iloc[np.argmax(sub1['ratio'])]['freq']
			#f1.append(fl1)
		f = np.arange(0.0,df0["freq"].iloc[-1].real, dF)

		ampl0 = f*0
		ampl2 = f*0
		ampl1 = f*0
		a_l0 = 1
		a_l1 = 0.5
		a_l2= 0.75

		for loc in df0['freq']:
			ampl0 += lorentzian(f,loc.real,a_l0,width)

		amp_total = ampl0 + ampl1 + ampl2
		b = np.exp(1.05*np.log(numax.iloc[0]) -1.91)*1.0
		b=0.25*numax.iloc[0]
		ampnm = np.exp( - (f - numax.iloc[0])**2/2.0/b**2)
		ofac = 1
		noise = (chi2.rvs(2*ofac, size=len(ampnm)) * amp_total*ampnm )/2./ofac

		#only reals when dealing with nonad 
		a_csv = amp_total*ampnm
		profile_name= ID
		total_parent_dir="/home/032272043/projects/dnuGridCode/temp_move/test_runs"


		if os.path.exists(total_parent_dir):
			pass
		else:
			os.mkdir(total_parent_dir)

		parent_dir_results= os.path.join(total_parent_dir,relevant_path)
		if os.path.exists(parent_dir_results):
			pass
		else:
			os.mkdir(parent_dir_results)
			path=os.path.join(parent_dir_results,profile_name)
		data_creation_path=os.path.join(parent_dir_results,profile_name)
		if os.path.exists(data_creation_path):
			pass
		else:
			os.mkdir(data_creation_path)

		data = {'freq':f,'amp':a_csv}
		name_prof_saves=ID.split('.')[0]
		df_csv = pd.DataFrame(data)
		df_csv.to_csv(data_creation_path+'/'+f'{name_prof_saves}freq_amp.csv',header=False, index=False)

		#Frequency X Frequency/delta nu Mod  ------ ALL MODES
		#freq_ech0 = df0['freq'].div(dnutheory) % 1
		#plt.clf()
		#plt.figure(figsize=(10,6))
		#plt.scatter(freq_ech0, df0['freq'].real, color="black", label="l = 0 modes")
		#plt.legend()
		#plt.title("Echelle Plot")
		#plt.xlabel('Frequency Modulo {} $\mu$Hz'.format(dnutheory))
		#plt.ylabel('Frequency (uHz)')
		#plt.savefig(path+"zeromodes.png")

		# Lorentzian Model
		plt.clf()
		plt.figure(figsize=(10,6))
		plt.plot(f, ampl0, color="black", label="l = 0 modes")
		plt.legend()
		plt.title("Lorentzian Model")
		plt.xlabel('Frequency (uHz)')
		plt.ylabel('Amplitude')
		plt.xlim(numax[0]-dnutheory*3,numax[0]+dnutheory*3)
		plt.savefig(data_creation_path+'/'+f"{name_prof_saves}FreqAmp.png")

		#Background corrected power spectrum
		plt.clf()
		plt.figure(figsize=(10,6))
		plt.plot(f, ampnm, color="green")
		plt.plot(f, ampnm*amp_total, color="blue")
		plt.title("Background corrected power spectrum")
		plt.xlabel('Frequency (uHz)')
		plt.ylabel('Amplitude') #Signal to noise ratio
		plt.xlim(numax[0]-dnutheory*3,numax[0]+dnutheory*3)
		plt.savefig(data_creation_path+'/'+f"{name_prof_saves}Fgaussian.png")
 
		#Power Spectrum Plot
		plt.clf()
		plt.figure(figsize=(10,6))
		plt.plot(f, noise, color="black")
		plt.title("Power Spectrum")
		plt.xlabel('Frequency (uHz)')
		plt.ylabel('Amplitude')
		plt.xlim(numax[0]-dnutheory*3,numax[0]+dnutheory*3)
		plt.savefig(data_creation_path+'/'+f"{name_prof_saves}Fpowerspec.png")
 
		#return the list that goes in the csv
 
		comp=[name_prof_saves,fdnuStdLit,dnutheory, dnuinf,numax.iloc[0],df['M_star'].iloc[0],df['R_star'].iloc[0]/R_sun,df['L_star'].iloc[0],1.7, '0.0    d0']
		return comp

	def pat_match_files(self):
        	directory_path = '/research/CNSM-JZinn/nonad/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha/LOGS'
        	files=os.listdir(directory_path)
        	pattern=r'^profile(\d+)\.data\.GYRE\.gyre_nad.eigval.h5'
        	matching_files= [f for f in files if re.match(pattern,f)]
        	return matching_files, directory_path 

	def make_oldlit_csv(self):

		resulting_f,d_path=self.pat_match_files()
		suffix_to_remove='.gyre_nad.eigval.h5'
		trimmed=[file.replace(suffix_to_remove,'') for file in resulting_f]
		bad_files=[]
		filename='old_lit_nonad_1M0_0D01.7ALPHA.csv'

		header=['PROFID','fdnustdlit','dnu_model','dnu_inf','nu_max','mass','radius','luminosity', 'mixing_length', 'metalicity']
		try:
			with open(filename,'w') as csvfile:
				csvwriter=csv.writer(csvfile)
				csvwriter.writerow(header)
		except FileExistsError:
			print(f'{filename} already exists, no header or file creation')

		with open(filename,'a',newline='') as csvfile:
			csvwriter=csv.writer(csvfile)
			for i in trimmed:
			#       try: 
				data_CSV=self.process(i,d_path)
				csvwriter.writerow(data_CSV)
				#       except Exception as e:
				#               bad_files.append(i)
				#               continue
			#print(bad_files)

	def create_BAM_data(self):
		# all of this needs to be read after we do 
		original_oldlist_datapath='/home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha'
		
		integers=[]
		pattern=re.compile(r'profile(\d+)\.data\.GYRE')
		for item in os.listdir(original_oldlist_datapath):
		        match=pattern.match(item)
		        if match:
		                integer=int(match.group(1))
		                integers.append(integer)
		
		command1='mkdir -p /home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha/list_of_bam_runs'
		new_direc='/home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha/list_of_bam_runs'
		
		result1=subprocess.run(command1,shell=True,capture_output=True,text=True)
		# we now need to get the csv from the individual profiles for further plotting 
		os.chdir(new_direc)
		
		for i in integers :
		        command3=f'cp /home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha/profile{i}.data.GYRE/profile{i}freq_amp.csv  /home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha/list_of_bam_runs  '
		
		        result3=subprocess.run(command3,shell=True)
		#now we have all of the csv's in the right place, now
		#create the txt files for pySYD/BAM
		for filename in os.listdir('.'):
        		if filename.endswith('.csv'):
                		new_filename=filename[:-4]+'.txt'
                		os.rename(filename,new_filename)

		# Define the directory containing the text files
		directory = '/home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha/list_of_bam_runs/'
		# Get a list of all .txt files in the directory
		txt_files = glob.glob(os.path.join(directory, '*.txt'))

		# Loop through each file
		#we need to take replace all ',' with ' ' cause BAM needs that	
		for file_path in txt_files:
		    # Open the file, read its contents, and replace commas with spaces
		    with open(file_path, 'r') as file:
		        content = file.read()
		    
		    # Replace commas with spaces
		    modified_content = content.replace(',', ' ')
		    
		    # Write the modified content back to the file
		    with open(file_path, 'w') as file:
		        file.write(modified_content)
		
		print("All txt files have replaced ',' with ' '")
		# we need to but run_this_w_preprocess in our working directory
		

		#MAKE THE DEFINITION OF THE REPLACING COMMAS WITH SPACED THAT IS REQUIRED IN THIS HOE
		#USE SELF IN THIS HOE
		#REMEBER THE HIPASS RUNPREPROCESS SETTINGS THOSE ARE EXTREMELY IMPORTANT IN THIS HOE

		#command_copy_preprocess=f'cp $BAMHOME/run_this_w_preprocess /home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha/list_of_bam_runs' 
		#subprocess.run(command_copy_preprocess,shell=True)
		
		# read in all txt files
		# replace instances of ',' with ' '


		# we literally dont create the dnu like do the dnu I don't know how we are not doing it but like oops

				
		for k in integers:
		        command4=f'ls profile{k}freq_amp.txt > profile{k}_list_run'
		        result4=subprocess.run(command4,shell=True)
		#we need to reset the working directory
		for j in integers:
		        command5=f'$BAMHOME/preprocess.sh run_this_w_preprocess profile{j}_list_run'
		        result5=subprocess.run(command5,shell=True)
			


	def create_bam_csv(self):
		
		directory_main_path='/home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha'
		directory_get_txt='/home/032272043/projects/dnuGridCode/temp_move/test_runs/nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha/list_of_bam_runs'
		# create the csv once, just with the headers (column names)
		if not os.path.exists('1M0Z_BAM.csv'):
			pd.DataFrame(columns=['ID','dnu_BAM','dnu_unc_BAM']).to_csv('1M0Z_BAM.csv',index=False)

		# get all profile numbers here
		pattern=re.compile(r'profile(\d+)\.data\.GYRE')
		integers=[]
		for item in os.listdir(directory_main_path):
			match=pattern.match(item)
			if match:
				integers.append(int(match.group(1)))
		# I feel like we can just do this 
		# f'profile{k}_list_run.globalpars'
		# the global pars has the numac and the dnu
		# so like just use this it doesn't necesarily matter in this context

		#list_run.globalpar(this will set up some default behavior in how it fits dnu)

#>>> $BAMHOME/fitdnu.sh FILE.bgcorr fitbgsp_syd_dnu NUMAX

#Where NUMAX would be a number equal to numax for the star (e.g., 50.12).

#This will create a file called FILE.bgcorr.dnu, which, at time of writing, has two duplicate rows, which have the dnu in the first column and the uncertainty on dnu in the second.
		# first we need a NUMAX which is in global pars


		return
#		for j in integers:
#		        command5=f'$BAMHOME/fitdnu.sh profile{j}freq_amp.txt.clean.hipass.fill.psd.bgcorr fitbgsp_syd_dnu  {NUMAX}'
#		        result5=subprocess.run(command5,shell=True)
			



		for k in integers:
			try: 
				file_path=f'{directory_get_txt}/profile{k}freq_amp.txt.dnu'
				df=pd.read_csv(file_path,delim_whitespace=True,header=None,names=['dnu','dnu_uncertainty'])
				dnu_bam=df['dnu'].iloc[-1]
				dnu_unc_bam=df['dnu_uncertainty'].iloc[-1]

				pd.DataFrame([[f'profile{k}',dnu_bam_dnu_unc_bam]],
						coluns=['ID','dnu_BAM','dnu_unc_BAM'])\
					.to_csv('1M0Z_BAM.csv',mode='a',header=False,index=False)					
				#old
				#df_main=pd.DataFrame([[f'profile{k}',dnu_bam,dnu_unc_bam]],columns=[f'profile{k}','dnuBam','dnuUnc'])
				#df_main.to_csv('1M0Z_BAM.csv',mode='a',header=False,index=False)

			except FileNotFoundError:
				print(f'BAM file not found for profile{k}')
			except pd.errors.EmptyDataError:
				print(f'Empty BAM file for profile{k}')	
			except Exception as e:
				print(f'Error processing file profile{k} with error: {str(e)}')



	def asf_grid(self):
		#run it here
		# asf_create_txttorun.py
		# asf_run_txt.py
		# asf_csv_create.py
		import_pysyd_csv='/home/032272043/projects/dnuGridCode/1M0_0D01.7ALPHA_SYDDNU_update_radius.csv'
		sydadd2asf=pd.read_csv(import_pysyd_csv)
		headers=['evstate','logz','teff','mini','mass','logg']
		df=pd.DataFrame(columns=headers)
		def append_rows(row_data,outputname):
		        df=pd.DataFrame(row_data,columns=headers)
		        df.to_csv(outputname,mode='a',sep=' ',index=False,header=False)
		radius_solar=6.958e10
		mass_solar=1.989e33
		gravconst=6.67259e-8
		ID=sydadd2asf['ID']
		for i in range(len(sydadd2asf['ID'])):
		        rad_scale=sydadd2asf['radius'][i]
		        outputname=(f'{ID[i]}.txt')
		        df.to_csv(outputname,index=False,sep=' ')
		        grav_solar_log=np.log10((gravconst*mass_solar/(radius_solar*rad_scale)**2))
		#       for i in range(len(sydadd2asf['radius'])):
		        rowin=[[1,-1.721,5770,1.0,1.0,grav_solar_log]]
		
		        append_rows(rowin,outputname)
	
		
#	def sync_bam_lit(self):
		# the data is messy and we want it to be clean in the sense that all of the column names are the same, and that there is a consistent data sheet as in: 
#	def clean_asf_grid(self):
	# do the same thing but with asf grid 

exp_test=full_create()
#exp_test.make_oldlit_csv()
# we can do pysyd right here
# whatever=whatever the class is called then call from here
#exp_test.create_BAM_data()
exp_test.create_bam_csv()
#in order for the asf grid:
# asf_create_txttorun.py
# asf_run_txt.py
# asf_csv_create.py 
# I think we have to do an added radius to everything
# how else can we have the radius added? I feel like this is best for now unless we completely re do it. 
# lets just do what we do now, we can make it better later when we do all the changes and all that
# we are going to make a clean all data files at once before we plot
# we want to now create asfgrid data and all that 
# make_oldlist needs to happen first


# WHEN DO WE DO PYSYD?



# MAYBE MAIN OUR OWN INTEGER GET METHOD? ACTUALLY PROBABLY JUST DO IT IN INIT
# cause we can just have the integer whole class wide so that it can be upgraded either way
