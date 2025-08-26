AAAimport nonad
import glob
import subprocess
from pudb import set_trace as pause
gyre_command = '$GYRE_DIR/bin/gyre'

b = []

dirs = glob.glob('/Users/sarahmedina/PMS/mesa/nonad/*')
print(dirs)
dfs = []
for d in dirs:
    print(d)
    profiles = glob.glob(d+'/LOGS/*GYRE*')
    
    mass = float(d.split('msun')[0][-3:])
    feh = (d.split('feh')[0].split('_')[2])
    alpha = float(d.split('alpha')[0][-3:])
    
    for profile in profiles:
        b.append(profile.split('.GYRE')[0])
    for profile in b :
        gyre_infile_tmp = 'gyre_tmp.in'
        gyre_infile_template = '/Users/sarahmedina/PMS/codes/nonad/src/nonad/gyre_pms_template.in'
        with open(gyre_infile_tmp,'w') as outfile:
            with open(gyre_infile_template, 'r') as infile:
                outfile.write("".join(infile.readlines()).format(profile=profile))

        subprocess.call(gyre_command + ' ' + gyre_infile_tmp,shell=True)
    print(profile)
    for profile in b:
        try:
            nonad.plot_modes.main(profile=profile, plot=False)
        except:
            pass
        
    import pandas as pd
    from os import listdir
    filepaths = [d+'/LOGS/'+f for f in listdir(d+'/LOGS/') if f.endswith('.astero.in')]
    # this is the file that has the frequency information.
    #ad='ad'
    #eigval_file = profile + '.gyre_'+ad+'.eigval.h5'
    #print('reading in')
    #print(eigval_file)
    
    from functools import partial
    mapfunc = partial(pd.read_csv,  sep= '\s+')
    if len(filepaths) > 0:
        df= pd.concat(map(mapfunc, filepaths))
    
        df['mass'] = mass
        df['feh'] = feh
        df['alpha'] = alpha
        dfs.append(df)
#    except:

    

#print (df)


# when I am ready to create mult plots
import pickle
#submit_run= { "masses": [0.8 ,1.0, 1.2], "fehs": [0.0, 0.2], "alphas": [1.7 , 2.0]}
#submit_pickle_run=pickle.dump(submit_run)

#df.pickle('submit_run.pickle')
with open('data.pickle', 'wb') as f:
   
    pickle.dump(dfs, f)

#spliting files
#FILE_PATH= '~/PMS/mesa/nonad/'
#for profile in profiles:
 #   new = profile.split('.GYRE')[0]
  #  if not os.path.isdir(profile):
   #     os.mkdir(profile)
   # os.rename(os.path.join(profile, new))
