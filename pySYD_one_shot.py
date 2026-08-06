import pandas as pd
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import string
import shutil
import csv

from multiprocessing.pool import ThreadPool
from threading import Lock
from pysyd import plots
from pysyd.target import Target
from pysyd.utils import Parameters
from pysyd import PACKAGEDIR
from pysyd.plots import plot_1d_ed
from pysyd import utils


class pySYD_oneshot:
	"""
    	A class to perform one-shot pySYD analysis on GYRE output files.
    
    	This class handles processing of stellar oscillation data, generating
    	power spectra, computing asteroseismic parameters, and creating echelle diagrams using the pySYD pipeline. See -> https://pysyd.readthedocs.io/en/latest/
    	
    	Attributes:
    	    MPLSTYLE (str): Path to matplotlib style file for pySYD plots.
    	"""

	def __init__(self,**kwargs):
		"""
        	Initialize the pySYD_oneshot class with matplotlib styling.
        
        	Args:
            	**kwargs: Additional keyword arguments (currently unused but 
                     allows for future flexibility).
   		"""

		MPLSTYLE = os.path.join(PACKAGEDIR,'data','pysyd.mplstyle')
		plt.style.use(MPLSTYLE)
		self.config={
			'dpi': 300,
			'smooth_width': 0.005,
			'oversampling_factor':5,
			'smooth_ps': 0.0,
			'show': True,
			'verbose':True,
			'save_echelle': True,
			**kwargs # overide defaults if needed
		}

	@staticmethod
	def li_spectrum_id(profile, mass, radius, metalicity):
		profile_match = re.search(r'\d+', str(profile))
		profile_num = profile_match.group(0) if profile_match else str(profile)
		return (
			f"profile{profile_num}_mass{float(mass):.4f}"
			f"_radius{float(radius):.4f}_metallicity{float(metalicity):.4f}_l0"
		)

	@classmethod
	def li_spectrum_id_from_row(cls, row):
		return cls.li_spectrum_id(row['PROFID'], row['mass'], row['radius'], row['metalicity'])
	

	def power_spec_filename_correction_nested(self,nonad_direc):
		"""
        	Convert CSV power spectrum files to pySYD-compatible TXT format.
        
        	Processes GYRE output files by converting comma-separated frequency-amplitude
        	data to space-separated format that pySYD can read.

		Supports two file naming patterns:
		1: profile{number}freq_amp.csv
		2: profile{number}_mass{...}.csv
        
        	Args:
            		nonad_direc (str): Directory containing the profiles.
			structure_type (srs): nested string meaning one csv file per directory, or flat, which is all of the files for one directory.
        	"""

		# we need to do this otherwise pySYD cannot read the values in the file correctly

		if not os.path.exists(nonad_direc):
			raise FileNotFoundError(f'Directory {nonad_direc} does not exist')

		for dir_name in os.listdir(nonad_direc):
			dir_path = os.path.join(nonad_direc, dir_name)
			# make sure that the directory exists and is created properly
			if os.path.isdir(dir_path):
				continue
			# start as below
			try:
				match=re.match(r'profile(\d+)',dir_name	)
				if not match:
					continue
				taken_number=match.group(1)
				csv_file_path = os.path.join(dir_path, f'profile{taken_number}freq_amp.csv')
				if not os.path.isfile(csv_file_path):
					continue
				pySYD_newpath = os.path.join(dir_path, f'{taken_number}_PS.txt')
				with open(csv_file_path, 'r') as csv_file:
					content = csv_file.read()

				lines=content.strip().split('\n')
				converted_lines=[]
				for line in lines:
					if ',' in line:
						converted_lines.append(line.replace(',',' '))
					else: 
						converted_lines.append(line)
				with open(pySYD_newpath, 'w') as new_file:
					new_file.write('\n'.join(converted_lines))
			except Exception as e: 
				print(f'Error processing {dir_name}: {e}')
				continue

	def power_spec_filename_correction_flat(self,csv_directory):
		if not os.path.exists(csv_directory):
			raise FileNotFoundError(f'Directory {csv_directory} does not exist')
		for file_name in os.listdir(csv_directory):
			file_path=os.path.join(csv_directory,file_name)
			if not os.path.isfile(file_path) or not file_name.endswith('_l0.csv'):
				continue		
			
			try: 
				source_id=os.path.splitext(file_name)[0]
				new_filename=f"{source_id}_PS.txt"
				pySYD_newpath=os.path.join(csv_directory,new_filename)

				with open(file_path, 'r') as csv_file:
					content=csv_file.read()

				lines=content.strip().split('\n')
				converted_lines=[]
				for line in lines: 
					if ',' in line: 
						converted_lines.append(line.replace(',', ' '))
					else: 
						converted_lines.append(line)
				with open(pySYD_newpath,'w') as new_file:
					new_file.write('\n'.join(converted_lines))

			except Exception as e:
				print(f'Error with {file_name}: {e}')
				continue
				

	def process_txt_file(self,file_path):
		"""
        	Process a text file to find the frequency with maximum amplitude.
        
        	Reads a two-column text file and returns the frequency value from the first
        	column that corresponds to the maximum amplitude in the second column.
        
        	Args:
            		file_path (str): Path to the text file containing frequency-amplitude data.
            
        	Returns:
            		float: The frequency value with the maximum amplitude.
        	"""
		data=[]
		with open(file_path,'r') as file:
			for line in file:
				split_line=line.split()
				col1,col2 = map(float,line.split())
				data.append((col1,col2))
			return max(data,key=lambda x: x[1])[0]

	def echelle_plot(self,bse_direc,name,star):
		"""
        	Generate and save an echelle diagram for a given star.
        	
        	Creates an echelle diagram using pySYD's plotting functionality and
        	saves it to the specified directory.
        	
        	Args:
        		bse_direc (str): Base directory where the plot will be saved.
        		name (int): Identifier for the star/profile.
        		star (Target): pySYD Target object containing seismic parameters.
        	"""

		fig, ax=plt.subplots(dpi=self.config['dpi'])
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
			#change below name 
			plt.savefig(f'{bse_direc}/ed_pysyd_oneshot_l0/ed_{name}.png', dpi=self.config['dpi'])
		except Exception as e:
			print(f"Cannot create profile{name} because of: {e}")
		finally:
			if fig:
				plt.close(fig)


	def delete_repeats(self,output_file,cleaned_file):
		"""
       		Remove duplicate entries from the results CSV file.
       		
       		Reads the output CSV file, removes any duplicate rows based on all columns,
       		and saves the cleaned data to a new file.
       		"""
		df=pd.read_csv(output_file) 
		if 'ID' in df.columns:
			df_no_duplicates=df.drop_duplicates(subset=['ID'],keep='last')
		else:
			df_no_duplicates=df.drop_duplicates(keep='last')
		df_no_duplicates.to_csv(cleaned_file,index=False)

	def full_pysyd_nested(self,bse_direc,nonad_direc,outputcsv):
		"""
        	Perform full pySYD analysis on all GYRE profiles in the directory.
        	
        	Processes each GYRE profile, computes seismic parameters (including Δν),
        	generates echelle diagrams, and saves results to a CSV file.
        	
        	Args:
        		bse_direc (str): Base directory for output file(s).
        		nonad_direc (str): Directory containing non-adiabatic profiles.
        	"""

		os.makedirs(f'{bse_direc}/ed_pysyd_oneshot',exist_ok=True)
		file_list=os.listdir(nonad_direc)
		filtered_files=[file for file in file_list if file.endswith('.data.GYRE')]
		numbers = [int(re.search(r'\d+', file).group()) for file in filtered_files]

		file_exists=os.path.isfile(outputcsv)
		if not file_exists:
                	with open(outputcsv, 'w') as file:

                        	writer=csv.writer(file)
                        	header=['ID','dnu_syd']
                        	writer.writerow(header)

		with open(outputcsv,mode='a',newline='') as file:
			writer=csv.writer(file)
			for name in numbers:

				params=Parameters()
				params.add_targets(stars=name)
				params.params[name]['show'],params.params[name]['verbose']=self.config['show'],self.config['verbose']

				star=Target(name, params)
				star.params['inpdir']=f'{nonad_direc}/profile{name}.data.GYRE/'
				star.params['outdir']=f'{bse_direc}'
				file_path=f'{nonad_direc}/profile{name}.data.GYRE/{name}_PS.txt'
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
				star.params['smooth_width']=self.config['smooth_width']
				star.params['numax']=col1_for_max_col2
				star.params['resolution']=star.frequency[1]-star.frequency[0]
				star.i = 0
				star.params['oversampling_factor']=self.config['oversampling_factor']
				star.params['smooth_ps']=self.config['smooth_ps']
				star.compute_acf()
				if not self.config['save_echelle']:
					star.echelle_diagram=lambda *args, **kwargs: None

				#remember that we need to change (if you update python to python3) to: 

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
					if self.config['save_echelle']:
						star.echelle_diagram()
				except Exception as e:
					print(f"profile{name} cant plot this in echelle because: {e} for echelle diragram")

				best_dnu=star.params['best_lag']
				print(f'This was the best_dnu found by pySYD for profile{name}: {best_dnu}')
				row=[name,best_dnu]
				writer.writerow(row)
				if self.config['save_echelle']:
					self.echelle_plot(bse_direc,name,star)

	def full_pysyd_flat(self,bse_direc,nonad_direc,outputcsv_flat,cores_mt):
		"""
        	Perform full pySYD analysis on all profiles in the directory.
        	
        	Processes each profile, computes seismic parameters (including Δν),
        	generates echelle diagrams, and saves results to a CSV file.
        	
        	Args:
        		bse_direc (str): Base directory for output file(s).
        		nonad_direc (str): Directory containing non-adiabatic profiles.
			cores_mt (int): Number of cores for multithreading.
        	"""
		os.makedirs(f'{bse_direc}/ed_pysyd_oneshot_l0',exist_ok=True)

		ps_files=[f for f in os.listdir(nonad_direc) if f.endswith('_PS.txt')] 
		unique_ps_pattern=re.compile(
			r'profile\d+_mass[-\d.]+_radius[-\d.]+_metallicity[-\d.]+_l0_PS\.txt$'
		)
		unique_ps_files=[f for f in ps_files if unique_ps_pattern.match(f)]
		if unique_ps_files:
			ps_files=unique_ps_files
		file_prefixes=[f[:-len('_PS.txt')] for f in ps_files]
		numbers = sorted(set(file_prefixes))

		if not numbers:
			print(f'No pySYD power spectra found in {nonad_direc}')
			return []

		with open(outputcsv_flat, 'w', newline='') as file:
			writer=csv.writer(file)
			header=['ID','dnu_syd']
			writer.writerow(header)



		# PUT threading here I think
		file_lock=Lock()
		def process_single_profile(full_name):		
			try:		
				params=Parameters()
				params.add_targets(stars=full_name)
				params.params[full_name]['show'],params.params[full_name]['verbose']=self.config['show'],self.config['verbose']

				star=Target(full_name, params)
				star.params['inpdir']=f'{nonad_direc}/'
				star.params['outdir']=f'{bse_direc}'
				print(star.params['inpdir'])
				print(f'this is the profile name: {full_name}')	
				
				
				file_path=os.path.join(nonad_direc,f'{full_name}_PS.txt')
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
				star.params['smooth_width']=self.config['smooth_width']
				star.params['numax']=col1_for_max_col2
				star.params['resolution']=star.frequency[1]-star.frequency[0]
				star.i = 0
				star.params['oversampling_factor']=self.config['oversampling_factor']
				star.params['smooth_ps']=self.config['smooth_ps']
				star.compute_acf()
				if not self.config['save_echelle']:
					star.echelle_diagram=lambda *args, **kwargs: None

				#remember that we need to change (if you update python to python3) to: 

				#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
				# below we want to put the thiing below on 
				# /home/032272043/projects/dnuGridCode/astero-one-shot/local/python-3.9.18/lib/python3.9/site-packages/scipy/optimize/_lsq/least_squares.p
				#from scipy.optimize import make_strictly_feasible
				#x0 = make_strictly_feasible(x0, lb, ub)
				#this is because the lowest value of dnu is not inbetween lb, ub that is calcualted for the frequency space
				#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

				try:
					star.frequency_spacing()
				except Exception as e:
					print(f"profile{full_name} doesnt work because:  {e} for frequency spacing")
				try:
					if self.config['save_echelle']:
						star.echelle_diagram()
				except Exception as e:
					print(f"profile{full_name} cant plot this in echelle because: {e} for echelle diragram")

				best_dnu=star.params['best_lag']
				print(f'This was the best_dnu found by pySYD for profile{full_name}: {best_dnu}')


				with file_lock:
					with open(outputcsv_flat, mode='a', newline='') as file:
						writer=csv.writer(file)
						row=[full_name,best_dnu]
						writer.writerow(row)
				if self.config['save_echelle']:
					self.echelle_plot(bse_direc,full_name,star)
				return full_name,best_dnu
			except Exception as e: 
				print(f'failed because of the fft: {e} this is profile {full_name}')
		with ThreadPool(processes=cores_mt) as pool:
			results=pool.map(process_single_profile, numbers)
		return results

	def sync_pysyd_csv_nested(self):
		"""
        	Synchronize pySYD results with additional stellar parameters.
        	
        	Merges the pySYD output with additional stellar parameters from the
        	literature CSV file, adding columns for dnu_inf, nu_max, mass, radius,
        	luminosity, mixing_length, and metallicity.
        	"""

		column_names_added=['dnu_inf','nu_max','mass','radius','luminosity','mixing_length','metalicity']
		added_radius_tofile=pd.read_csv(self.config['cleaned_csv'])
		csv_toextract_radius=pd.read_csv('old_lit_nonad_1M0_0D01.7ALPHA.csv')
		csv_toextract_radius['ID']=csv_toextract_radius['PROFID'].str.replace('profile','').astype(int)
		added_radius_tofile['ID']=added_radius_tofile['ID'].astype(int)
		for x in column_names_added:
			profile_radius_map=dict(zip(csv_toextract_radius['ID'],csv_toextract_radius[x]))
			added_radius_tofile[x]=added_radius_tofile['ID'].map(profile_radius_map)
		added_radius_tofile.to_csv(self.config['final_csv'],index=False)
			#'final_csv': '1M0_0D01.7ALPHA_SYDDNU_added_synced_cols.csv'

	def sync_pysyd_csv_flat(self,original_csv,cleaned_csv,synced_csv):
		"""
        	Synchronize pySYD results with additional stellar parameters.
        	
        	Merges the pySYD output with additional stellar parameters from the
        	literature CSV file, adding columns for dnu_inf, nu_max, mass, radius,
        	luminosity, mixing_length, and metallicity.
        	"""
		# we should really check later if this is appropriate/ the best approach to this

		df_tosync=pd.read_csv(cleaned_csv)
		df_read4sync=pd.read_csv(original_csv)

		df_tosync['pysyd_id']=df_tosync['ID'].astype(str)
		if 'pysyd_id' not in df_read4sync.columns:
			required_columns={'PROFID','mass','radius','metalicity'}
			if required_columns.issubset(df_read4sync.columns):
				df_read4sync['pysyd_id']=df_read4sync.apply(self.li_spectrum_id_from_row,axis=1)

		if 'pysyd_id' in df_read4sync.columns:
			df_read4sync=df_read4sync.drop_duplicates(subset=['pysyd_id'],keep='last')
			merged_df=pd.merge(
				df_tosync,
				df_read4sync,
				on='pysyd_id',
				how='inner',
				suffixes=('_pysyd','_theory')
			)
			if len(merged_df)>0:
				merged_df.to_csv(synced_csv,index=False)
				return
			print('No exact pySYD ID matches found; falling back to legacy profile/mass sync')

		def parse_id(id_value):
			id_str=str(id_value)
			if len(id_str)>=4:
				return int(id_str[:-4]), id_str[-4:]
			return None,None	

		df_tosync[['PROFID','mass_prefix']]=df_tosync['ID'].apply(
			lambda x: pd.Series(parse_id(x))
		)
		df_read4sync['mass_str']=df_read4sync['mass'].apply(lambda x: f"{x:.2f}")
		merged_df=pd.merge(
			df_tosync,
			df_read4sync,
			left_on=['PROFID','mass_prefix'],
			right_on=['PROFID','mass_str'],
			how='inner',
			suffixes=('_pysyd','_theory')
		)
		if len(merged_df)>0:
			merged_df.to_csv(synced_csv,index=False)
		else: 
			print("no matches found, need to so something different")
