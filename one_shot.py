import time
import os
import re
import numpy as np
import csv
import subprocess
import math
import glob
import shutil
#import random as rand
import pandas as pd
try:
    import pySYD_one_shot
except ModuleNotFoundError:
    pySYD_one_shot = None

from multiprocessing.pool import ThreadPool
#from scipy.interpolate import LinearNDInterpolator
#from scipy.interpolate import interp1d
try:
    from astropy.io import ascii
except ModuleNotFoundError:
    ascii = None
#from scipy.interpolate import RegularGridInterpolator
#from scipy.interpolate import NearestNDInterpolator
from pathlib import Path
#from matplotlib.pyplot import figure
#from nonadfiles.nonad_code.fileio import read_gyre
from scipy.stats import chi2
try:
    from nonadfiles.nonad_code.plot_modes import dnu_theory
    from nonadfiles.nonad_code.plot_modes import return_modes
    from nonadfiles.nonad_code.scaling_relation import logg_MR, loggteff2numax, numax, dnu_MR, freq_dyn, numax_sun, dnu_sun, L_sun, R_sun, M_sun
    from nonadfiles.nonad_code.dnu_util import lorentzian
except ModuleNotFoundError:
    dnu_theory = None
    return_modes = None
    logg_MR = None
    loggteff2numax = None
    numax = None
    dnu_MR = None
    freq_dyn = None
    numax_sun = 3090.0
    dnu_sun = 135.1
    L_sun = 3.828e33
    R_sun = 6.957e10
    M_sun = 1.98847e33
    lorentzian = None

SOLAR_TEFF_K = 5777.0
SOLAR_LOGG_CGS = 4.43796037457


def _read_apokasc_table(path):
    path = Path(path)
    if ascii is not None:
        try:
            return ascii.read(path).to_pandas()
        except Exception:
            pass
    data = pd.read_csv(path, sep=r'\s+', header=None)
    data.columns = [f'col{i}' for i in range(1, len(data.columns) + 1)]
    return data


def _li_spectrum_id(profile, mass, radius, metalicity):
    profile_match = re.search(r'\d+', str(profile))
    profile_num = profile_match.group(0) if profile_match else str(profile)
    return (
        f"profile{profile_num}_mass{float(mass):.4f}"
        f"_radius{float(radius):.4f}_metallicity{float(metalicity):.4f}_l0"
    )


class full_create:
    """Workflow driver for one-shot pySYD, BAM, ASFGrid, and Li runs."""

    def __init__(self, **kwargs):

        self.bse_numax = []
        requested_wrk_drc = kwargs.get('wrk_drc')
        legacy_wrk_drc = Path('/home/032272043/astero-one-shot')
        if requested_wrk_drc is not None:
            self.wrk_drc = str(Path(requested_wrk_drc).expanduser())
        elif legacy_wrk_drc.exists():
            self.wrk_drc = str(legacy_wrk_drc)
        else:
            self.wrk_drc = str(Path(__file__).resolve().parent)
        self.test_runs_dir = Path(self.wrk_drc)/'test_runs'
        self.nonad_create_dir = Path(self.test_runs_dir)/'nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha'
        self.bam_runs = Path(self.nonad_create_dir)/'list_of_bam_runs'
        self.parquet_runs_dir = Path(self.wrk_drc)/'li_sfc_runs'
        self.parquet_runs_dir_l0 = Path(self.wrk_drc)/'li_sfc_runs_l0'
        self.parquet_files = Path(self.wrk_drc)/'parquet_files'
        self.apoxe_radius_interp = Path(self.wrk_drc)/'apokasc.txt'
        self.apox_data = _read_apokasc_table(self.apoxe_radius_interp)
        self.apox_data_radius = self.apox_data['col19']
        self.apox_data_err_radius = self.apox_data['col20']
        self.fdnu_li = []
        self.std_header = []
        self.output_std = '1M0_0D01.7ALPHA_SYDDNU_one_shot.csv'
        self.cleaned_std = '1M0_0D01.7ALPHA_SYDDNU_cleaned.csv'
        self.synced_std= '1M0_0D01.7ALPHA_SYDDNU_added_synced_cols.csv'

        self.output_li = '1M_li_one_shot.csv'
        self.cleaned_li = '1M_li_cleaned_one_shot.csv'
        self.synced_cleaned_li= '1M_li_added_synced_cols.csv'

        #above has not been fully adjusted yet 12.2.2025
        self.output_li_l0 = '1M_li_one_shot_l0.csv'
        self.cleaned_li_l0 = '1M_li_cleaned_one_shot_l0.csv'
        self.synced_cleaned_li_l0= '1M_li_added_synced_cols_l0.csv'
        #get fdnu we calculaute already local from this code as of right now 2022-2-15, dnu_theory is the delta nu normal, the delta nu infinite is the mass an    d radius from return modes.
        # bam gets us the delta nu local from the fake power specturm whih is the l-0 modes, which is the new, better way of doing it I think it takes the auto        correlation function of the power spectrum at numax (local), 1?

        self.output_bam_corrected = '1M0Z_BAM_corrected.csv'
        self.output_bam_manifest = 'BAM_run_manifest_standard.csv'
        self.output_li_bam_l0_corrected = '1M_li_BAM_l0_corrected.csv'
        self.output_li_bam_l0_manifest = 'BAM_run_manifest_li_l0.csv'
        self.asfgrid_dir = Path(self.wrk_drc)/'asfgrid_code'/'asfgrid_v0.0.6'
        self.output_li_asfgrid_l0 = '1M_li_asfgrid_l0.csv'
        self.output_li_asfgrid_l0_manifest = 'ASFGrid_run_manifest_li_l0.csv'
        self.output_li_asfgrid_l0_input = 'li_asfgrid_l0_input.txt'
        self.path_li_asfgrid_fdnu_plot = Path(self.wrk_drc)/'plots_directory'/'li_asfgrid_fdnu.png'

    def read_power_spectrum(self, spectrum_path):
        data = pd.read_csv(
            spectrum_path,
            sep=r'[\s,]+',
            header=None,
            engine='python',
            comment='#'
        )
        data = data.iloc[:, :2].copy()
        data.columns = ['frequency', 'amplitude']
        data['frequency'] = pd.to_numeric(data['frequency'], errors='coerce')
        data['amplitude'] = pd.to_numeric(data['amplitude'], errors='coerce')
        data = data.dropna(subset=['frequency', 'amplitude'])
        if data.empty:
            raise ValueError(f'No numeric frequency/amplitude rows in {spectrum_path}')
        return data

    def parse_spectrum_line(self, line):
        if isinstance(line, bytes):
            line = line.decode(errors='ignore')
        parts = re.split(r'[\s,]+', line.strip())
        if len(parts) < 2:
            return None
        try:
            return float(parts[0]), float(parts[1])
        except ValueError:
            return None

    def quick_spectrum_summary(self, spectrum_path, numax_guess):
        spectrum_path = Path(spectrum_path)
        first_pair = None
        second_pair = None
        with open(spectrum_path, 'rb') as handle:
            for raw_line in handle:
                pair = self.parse_spectrum_line(raw_line)
                if pair is None:
                    continue
                if first_pair is None:
                    first_pair = pair
                else:
                    second_pair = pair
                    break

        last_pair = None
        with open(spectrum_path, 'rb') as handle:
            handle.seek(0, os.SEEK_END)
            position = handle.tell()
            buffer = b''
            while position > 0 and last_pair is None:
                read_size = min(8192, position)
                position -= read_size
                handle.seek(position)
                buffer = handle.read(read_size) + buffer
                for raw_line in reversed(buffer.splitlines()):
                    pair = self.parse_spectrum_line(raw_line)
                    if pair is not None:
                        last_pair = pair
                        break

        if first_pair is None or last_pair is None:
            return self.spectrum_summary(spectrum_path)
        if second_pair is not None:
            resolution = abs(second_pair[0] - first_pair[0])
        else:
            resolution = np.nan
        if not np.isfinite(resolution) or resolution <= 0:
            resolution = np.nan
        return {
            'numax_guess': float(numax_guess),
            'max_amplitude': np.nan,
            'freq_min': min(first_pair[0], last_pair[0]),
            'freq_max': max(first_pair[0], last_pair[0]),
            'n_points': np.nan,
            'resolution': resolution
        }

    def spectrum_summary(self, spectrum_path):
        data = self.read_power_spectrum(spectrum_path)
        max_idx = data['amplitude'].idxmax()
        frequency = data['frequency'].to_numpy(dtype=float)
        diffs = np.diff(np.sort(np.unique(frequency)))
        positive_diffs = diffs[diffs > 0]
        resolution = float(np.median(positive_diffs)) if len(positive_diffs) else np.nan
        return {
            'numax_guess': float(data.loc[max_idx, 'frequency']),
            'max_amplitude': float(data.loc[max_idx, 'amplitude']),
            'freq_min': float(data['frequency'].min()),
            'freq_max': float(data['frequency'].max()),
            'n_points': int(len(data)),
            'resolution': resolution
        }

    def scaling_dnu_from_numax(self, numax_value):
        numax_value = max(float(numax_value), 1e-6)
        return float(dnu_sun * (numax_value / numax_sun) ** 0.75)

    def bam_bounds(self, summary, expected_dnu):
        expected_dnu = float(expected_dnu) if pd.notna(expected_dnu) and expected_dnu > 0 else self.scaling_dnu_from_numax(summary['numax_guess'])
        resolution = summary['resolution']
        if not np.isfinite(resolution) or resolution <= 0:
            resolution = max(expected_dnu / 20.0, 1e-4)

        freq_min = float(summary['freq_min'])
        freq_max = float(summary['freq_max'])
        numax_guess = float(summary['numax_guess'])
        ps_half_width = max(6.0 * expected_dnu, 0.25 * numax_guess, 5.0 * resolution)
        lower_freq = max(freq_min + resolution, numax_guess - ps_half_width)
        upper_freq = min(freq_max, numax_guess + ps_half_width)
        if upper_freq <= lower_freq:
            lower_freq = max(freq_min + resolution, freq_min)
            upper_freq = freq_max

        noise_high = max(freq_min + resolution, freq_max)
        noise_width = max(10.0 * expected_dnu, 0.1 * max(freq_max - freq_min, resolution), 5.0 * resolution)
        noise_low = max(freq_min, noise_high - noise_width)
        if noise_low >= noise_high:
            noise_low = max(freq_min, noise_high - resolution)

        max_lag_span = max((upper_freq - lower_freq) * 0.9, resolution)
        lag_lower = max(resolution, expected_dnu * 0.2)
        lag_upper = min(max(expected_dnu * 3.0, lag_lower + resolution), max_lag_span)
        if lag_upper <= lag_lower:
            lag_upper = lag_lower + resolution

        return {
            'lower_freq': lower_freq,
            'upper_freq': upper_freq,
            'noise_low': noise_low,
            'noise_high': noise_high,
            'lag_lower': lag_lower,
            'lag_upper': lag_upper,
            'smooth_ps': min(0.5, max(resolution, expected_dnu * 0.2)),
            'expected_dnu': expected_dnu,
            'resolution': resolution
        }

    def bam_parameter_text(self, bounds, mc_iter=500, enable_plots=False):
        plot_flag = 1 if enable_plots else 0
        return f""";;;; general parameters
{bounds['lower_freq']:.8g}		; low frequency limit
{bounds['upper_freq']:.8g}		; upper frequency limit
{int(mc_iter)}		; number of Monte-Carlo iterations (set 1 for no uncertainties)
''      	; name of file (in DnupipeV40/) containing frq range to be excluded for individual stars 3-col:[KICID fmin fmax]; one line per range
5.		; oversampling factor of input spectrum
;;;; background fit
2		; number of Harvey components in background fit
0   	    	; try to force all Harvey components to be in the final fit? (0-no,1-yes)
{bounds['noise_low']:.8g}		; lower limit for noise estimate
{bounds['noise_high']:.8g}		; upper limit for noise estimate
10		; smoothing width (# of frequency points)
30  	    	; width for independent bins for smoothing/fitting
1		; fix white noise component in fit?
0		; fix frequency intervals for components?
0   		; times*numax is fmax for fit 1st harvey model in 1st iteration; set to 0 if all components should be fit at once
;;;; deltanu determination
{bounds['lag_lower']:.8g}		; lower limit for lag
{bounds['lag_upper']:.8g}		; upper limit for lag
10.		; number of peaks for fitting AC
{bounds['smooth_ps']:.8g}		; smoothing of PS (muHz)
;;;; misc parameters
{plot_flag}		; enable plotting?
0  	    	; check results step by step (0-no,1-yes)?
{plot_flag}		; save screenshots?
screenshots_new/	;
0		; fit using python bayes module?
0		; additionally add a gaussian bump in the bayesian fitting?
1		; replace output files instead of appending stale rows
0		; replot only?
0		; bayesian dnu calculation?
"""

    def find_bam_executable(self, bam_executable=None):
        candidates = []
        if bam_executable:
            candidates.append(Path(bam_executable).expanduser())
        bam_home = os.environ.get('BAMHOME')
        if bam_home:
            candidates.append(Path(bam_home).expanduser() / 'fitdnu.sh')
        which_fitdnu = shutil.which('fitdnu.sh')
        if which_fitdnu:
            candidates.append(Path(which_fitdnu))
        candidates.append(Path('/research/CNSM-JZinn/BAM/fitdnu.sh'))

        for candidate in candidates:
            if candidate.exists():
                if os.access(candidate, os.X_OK):
                    return [str(candidate)]
                return ['bash', str(candidate)]
        checked = ', '.join(str(path) for path in candidates)
        raise FileNotFoundError(f'Could not find BAM fitdnu.sh. Checked: {checked}')

    def build_bam_manifest_from_spectra(self, spectra, metadata, output_csv, id_from_path, metadata_id_column, param_dir, mc_iter=500):
        param_dir = Path(param_dir)
        param_dir.mkdir(parents=True, exist_ok=True)
        metadata = metadata.copy()
        metadata[metadata_id_column] = metadata[metadata_id_column].astype(str)
        metadata_by_id = metadata.set_index(metadata_id_column, drop=False)
        rows = []

        for spectrum_path in spectra:
            spectrum_path = Path(spectrum_path)
            object_id = id_from_path(spectrum_path)
            if object_id not in metadata_by_id.index:
                continue
            meta_row = metadata_by_id.loc[object_id]
            if isinstance(meta_row, pd.DataFrame):
                meta_row = meta_row.iloc[-1]
            metadata_numax = meta_row.get('numax_raw', meta_row.get('nu_max', np.nan))
            if pd.notna(metadata_numax) and float(metadata_numax) > 0:
                summary = self.quick_spectrum_summary(spectrum_path, float(metadata_numax))
            else:
                summary = self.spectrum_summary(spectrum_path)
            expected_dnu = meta_row['dnu_theory'] if 'dnu_theory' in meta_row else meta_row.get('dnu_model', meta_row.get('dnu_inf', np.nan))
            bounds = self.bam_bounds(summary, expected_dnu)
            param_path = param_dir / f'{object_id}_fitbgsp_syd_dnu'
            old_dnu_path = Path(f'{spectrum_path}.dnu')
            old_dnu = np.nan
            old_dnu_unc = np.nan
            old_status = 'missing'
            if old_dnu_path.exists():
                try:
                    old_df = pd.read_csv(old_dnu_path, sep=r'\s+', header=None, names=['dnu', 'uncertainty'], engine='python')
                    if old_df.empty:
                        old_status = 'empty'
                    else:
                        old_dnu = float(old_df['dnu'].iloc[-1])
                        old_dnu_unc = float(old_df['uncertainty'].iloc[-1])
                        old_status = 'ok' if old_dnu > 0 and old_dnu_unc >= 0 else 'failure'
                except Exception:
                    old_status = 'parse_error'

            rows.append({
                'ID': object_id,
                'spectrum_path': str(spectrum_path),
                'param_path': str(param_path),
                'numax_guess': summary['numax_guess'],
                'max_amplitude': summary['max_amplitude'],
                'freq_min': summary['freq_min'],
                'freq_max': summary['freq_max'],
                'n_points': summary['n_points'],
                'resolution': bounds['resolution'],
                'expected_dnu': bounds['expected_dnu'],
                'lower_freq': bounds['lower_freq'],
                'upper_freq': bounds['upper_freq'],
                'noise_low': bounds['noise_low'],
                'noise_high': bounds['noise_high'],
                'lag_lower': bounds['lag_lower'],
                'lag_upper': bounds['lag_upper'],
                'smooth_ps': bounds['smooth_ps'],
                'mc_iter': int(mc_iter),
                'old_dnu_path': str(old_dnu_path),
                'old_dnu_BAM': old_dnu,
                'old_dnu_unc_BAM': old_dnu_unc,
                'old_bam_status': old_status,
                'dnu_inf': meta_row.get('dnu_inf', np.nan),
                'nu_max': meta_row.get('nu_max', meta_row.get('numax_raw', np.nan)),
                'mass': meta_row.get('mass', np.nan),
                'radius': meta_row.get('radius', np.nan),
                'luminosity': meta_row.get('luminosity', np.nan),
                'mixing_length': meta_row.get('mixing_length', np.nan),
                'metalicity': meta_row.get('metalicity', np.nan)
            })

        manifest = pd.DataFrame(rows)
        manifest.to_csv(output_csv, index=False)
        return manifest

    def write_bam_parameter_file(self, manifest_row, enable_plots=False):
        bounds = {
            'lower_freq': manifest_row['lower_freq'],
            'upper_freq': manifest_row['upper_freq'],
            'noise_low': manifest_row['noise_low'],
            'noise_high': manifest_row['noise_high'],
            'lag_lower': manifest_row['lag_lower'],
            'lag_upper': manifest_row['lag_upper'],
            'smooth_ps': manifest_row['smooth_ps']
        }
        param_path = Path(manifest_row['param_path'])
        param_path.parent.mkdir(parents=True, exist_ok=True)
        param_path.write_text(self.bam_parameter_text(bounds, mc_iter=manifest_row.get('mc_iter', 500), enable_plots=enable_plots))
        return param_path

    def collect_bam_csv_from_manifest(self, manifest_csv, output_csv):
        manifest = pd.read_csv(manifest_csv)
        rows = []
        for _, row in manifest.iterrows():
            dnu_path = Path(f"{row['spectrum_path']}.dnu")
            dnu_bam = np.nan
            dnu_unc = np.nan
            status = 'missing'
            if dnu_path.exists():
                try:
                    dnu_df = pd.read_csv(dnu_path, sep=r'\s+', header=None, names=['dnu', 'uncertainty'], engine='python')
                    if dnu_df.empty:
                        status = 'empty'
                    else:
                        dnu_bam = float(dnu_df['dnu'].iloc[-1])
                        dnu_unc = float(dnu_df['uncertainty'].iloc[-1])
                        if dnu_bam <= 0:
                            status = 'failure_code'
                        elif dnu_unc < 0:
                            status = 'negative_uncertainty'
                        else:
                            status = 'ok'
                except Exception:
                    status = 'parse_error'

            fdnu_bam = row['dnu_inf'] / dnu_bam if pd.notna(row.get('dnu_inf')) and dnu_bam and dnu_bam > 0 else np.nan
            rows.append({
                'ID': row['ID'],
                'dnu_BAM': dnu_bam,
                'dnu_unc_BAM': dnu_unc,
                'fdnu_BAM': fdnu_bam,
                'bam_status': status,
                'spectrum_path': row['spectrum_path'],
                'param_path': row['param_path'],
                'numax_guess': row['numax_guess'],
                'expected_dnu': row['expected_dnu'],
                'dnu_inf': row.get('dnu_inf', np.nan),
                'nu_max': row.get('nu_max', np.nan),
                'mass': row.get('mass', np.nan),
                'radius': row.get('radius', np.nan),
                'luminosity': row.get('luminosity', np.nan),
                'mixing_length': row.get('mixing_length', np.nan),
                'metalicity': row.get('metalicity', np.nan)
            })

        output = pd.DataFrame(rows)
        output.to_csv(output_csv, index=False)
        return output

    def run_bam_manifest(self, manifest_csv, output_csv, cores_mt=1, bam_executable=None, run_limit=None, enable_plots=False):
        command_prefix = self.find_bam_executable(bam_executable)
        manifest = pd.read_csv(manifest_csv)
        if run_limit is not None:
            manifest = manifest.head(int(run_limit))

        def run_one(row_dict):
            param_path = self.write_bam_parameter_file(pd.Series(row_dict), enable_plots=enable_plots)
            cmd = command_prefix + [row_dict['spectrum_path'], str(param_path), str(row_dict['numax_guess'])]
            result = subprocess.run(cmd, capture_output=True, text=True)
            return {
                'ID': row_dict['ID'],
                'returncode': result.returncode,
                'stdout': result.stdout[-1000:],
                'stderr': result.stderr[-1000:]
            }

        records = manifest.to_dict('records')
        if cores_mt and cores_mt > 1:
            with ThreadPool(processes=cores_mt) as pool:
                run_log = pool.map(run_one, records)
        else:
            run_log = [run_one(row) for row in records]
        pd.DataFrame(run_log).to_csv(Path(output_csv).with_suffix('.runlog.csv'), index=False)
        return self.collect_bam_csv_from_manifest(manifest_csv, output_csv)

    def build_standard_bam_manifest(self, output_csv=None, mc_iter=500):
        output_csv = output_csv or self.output_bam_manifest
        metadata = pd.read_csv('old_lit_nonad_1M0_0D01.7ALPHA.csv')
        metadata['ID'] = metadata['PROFID'].str.replace('profile', '', regex=False)
        spectra = sorted(Path(self.bam_runs).glob('profile*freq_amp.txt'))
        return self.build_bam_manifest_from_spectra(
            spectra,
            metadata,
            output_csv,
            id_from_path=lambda path: re.search(r'profile(\d+)freq_amp\.txt$', path.name).group(1),
            metadata_id_column='ID',
            param_dir=Path(self.bam_runs) / 'bam_params_corrected',
            mc_iter=mc_iter
        )

    def build_li_bam_l0_manifest(self, output_csv=None, mc_iter=500):
        output_csv = output_csv or self.output_li_bam_l0_manifest
        metadata = pd.read_csv('li_surfacecorrection_mass_const_l0.csv')
        metadata['ID'] = metadata.apply(
            lambda row: _li_spectrum_id(row['PROFID'], row['mass'], row['radius'], row['metalicity']),
            axis=1
        )
        spectra = sorted(Path(self.parquet_runs_dir_l0).glob('profile*_l0_PS.txt'))
        return self.build_bam_manifest_from_spectra(
            spectra,
            metadata,
            output_csv,
            id_from_path=lambda path: path.name[:-len('_PS.txt')],
            metadata_id_column='ID',
            param_dir=Path(self.parquet_runs_dir_l0) / 'bam_params_corrected',
            mc_iter=mc_iter
        )

    def process(self, ID, read_path):
        font = {'family': 'DejaVu Sans', 'weight': 'normal', 'size': 16}
        plt.rc('font', **font)
        t = 4*math.pi*10
        nm = 40
            # based on what we did before gamma should be 0.12 so width should be changed to 0.0036, dont know where the original width = nm *.005 comes from?
            # based on this : https://iopscience.iop.org/article/10.3847/1538-4357/835/2/172/pdf
        width = 0.0036
            #width = nm*.005
        dF = 1/(2*t)

        parts = read_path.split('/')
        relevant_path = parts[4]
        full_path = read_path+'/'
        df = return_modes(full_path+ID, ad='nad', visibility_thresh=0.0,).reset_index(drop=True)
        df["freq"] = pd.to_numeric(df["freq"], downcast="float")
        df["ratio"] = pd.to_numeric(df["ratio"], downcast="float")
        df0 = df[df.l == 0]
        df1 = df[df.l == 1]
        ef2 = df[df.l == 2]

        #multiply by delta nu sun, microhertz multiply by right side
        numax = (df['M_star']/M_sun) * (df['R_star']/R_sun) ** (-7./4.) * (df['L_star']/L_sun) ** (-1./8.) * numax_sun
        dnutheory = dnu_theory(df["l"], df["freq"], numax.loc[0])
        #dnutheory is the delta nu local, as the current literature does it, the weighted sum of the l=0 modes
        #below is the deltanu infinite
        dnuinf = math.sqrt((df['M_star'].iloc[0]/M_sun)/((df['R_star'].iloc[0]/R_sun)**3)) * dnu_sun

        fdnuStdLit = dnuinf/dnutheory
        #below is fdnu which is the sqrt of the density (infinite delta nu) divided by the dnu local (as calculated by the weights of the power spectru    m around the nu local), or how the curren literature does it.  fdnuStdLit= dnuinf/dnutheory

        #f1=[]
        #f2 = []
        #for i in np.arange(len(df0['freq'])-1):
        #    try:
        #        flow= df0['freq'].iloc[i]
        #        fhi = df0['freq'].iloc[i + 1]
        #        sub1 = df1.loc[(df1['freq'] > flow) & (df1['freq'] < fhi)]
        #        fl1 = sub1.iloc[np.argmax(sub1['ratio'])]['freq']
        #        f1.append(fl1)
        #    except ValueError:
        #        if not input:
        #            raise ValueError('no data ')
        #    else:
        #        continue
        #f1=pd.to_numeric(f1,downcast="float")
        #print(f1)
        for i in np.arange(len(df0['freq'])-1):
            flow = df0['freq'].iloc[i]
            fhi = df0['freq'].iloc[i + 1]
            sub1 = df1.loc[(df1['freq'] > flow) & (df1['freq'] < fhi)]
            #print(sub1)
            #jfl1 = sub1.iloc[np.argmax(sub1['ratio'])]['freq']
            #f1.append(fl1)
        f = np.arange(0.0, df0["freq"].iloc[-1].real, dF)

        ampl0 = f*0
        ampl2 = f*0
        ampl1 = f*0
        a_l0 = 1
        a_l1 = 0.5
        a_l2= 0.75

        for loc in df0['freq']:
            ampl0 += lorentzian(f, loc.real, a_l0, width)

        amp_total = ampl0 + ampl1 + ampl2
        b = np.exp(1.05*np.log(numax.iloc[0]) - 1.91)*1.0
        b = 0.25*numax.iloc[0]
        ampnm = np.exp(- (f - numax.iloc[0])**2/2.0/b**2)
        ofac = 1
        noise = (chi2.rvs(2*ofac, size=len(ampnm)) * amp_total*ampnm)/2./ofac

        #only reals when dealing with nonad
        a_csv = amp_total*ampnm
        profile_name = ID

        if os.path.exists(self.test_runs_dir):
            pass
        else:
            os.mkdir(self.test_runs_dir)

        parent_dir_results = os.path.join(self.test_runs_dir, relevant_path)
        if os.path.exists(parent_dir_results):
            pass
        else:
            os.mkdir(parent_dir_results)
            path = os.path.join(parent_dir_results, profile_name)
        data_creation_path = os.path.join(parent_dir_results, profile_name)
        if os.path.exists(data_creation_path):
            pass
        else:
            os.mkdir(data_creation_path)

        data = {'freq': f, 'amp': a_csv}
        csvwriter=csv.writer(csvfile)
        csvwriter.writerow(self.std_header)
        #except FileExistsError:
        #    print(f'{filename} already exists, no header or file creation')

        with open(filename,'a',newline='') as csvfile:
            csvwriter=csv.writer(csvfile)
            for i in trimmed:
            #        try:
                data_CSV=self.process(i,d_path)
                csvwriter.writerow(data_CSV)
                #        except Exception as e:
                #                bad_files.append(i)
                #                continue
            #print(bad_files)

    def create_BAM_data(self, cores_mt=1, run_limit=None, bam_executable=None, mc_iter=500):
        manifest = self.build_standard_bam_manifest(mc_iter=mc_iter)
        print(f'Prepared {len(manifest)} standard BAM spectra in {self.output_bam_manifest}')
        return self.run_bam_manifest(
            self.output_bam_manifest,
            self.output_bam_corrected,
            cores_mt=cores_mt,
            bam_executable=bam_executable,
            run_limit=run_limit,
            enable_plots=False
        )

    def create_bam_csv(self):
        if not Path(self.output_bam_manifest).exists():
            self.build_standard_bam_manifest()
        return self.collect_bam_csv_from_manifest(self.output_bam_manifest, self.output_bam_corrected)

    def sync_bam_csvformat(self):

        column_names_added=['dnu_inf','nu_max','mass','radius','luminosity','mixing_length','metalicity']
        added_radius_tofile=pd.read_csv('1M0Z_BAM.csv')
        csv_toextract_radius=pd.read_csv('old_lit_nonad_1M0_0D01.7ALPHA.csv')
        csv_toextract_radius['ID']=csv_toextract_radius['PROFID'].str.replace('profile','').astype(int)
        added_radius_tofile['ID']=added_radius_tofile['ID'].str.replace('profile','').astype(int)
        for x in column_names_added:
            profile_radius_map=dict(zip(csv_toextract_radius['ID'],csv_toextract_radius[x]))
            added_radius_tofile[x]=added_radius_tofile['ID'].map(profile_radius_map)
            added_radius_tofile.to_csv('1M0Z_BAM_added_synced_cols.csv',index=False)

    def asf_grid_datacreate(self,data_created=False):
        #create new directory if you want, this is just for the current directory when beign ran on server
        asf_data_direc=f'{self.wrk_drc}/asfgrid_code/asfgrid_v0.0.6/'
        os.chdir(asf_data_direc)
        ## were doing this based on pysyd data as below
        sydadd2asf=pd.read_csv(f'{self.wrk_drc}/1M0_0D01.7ALPHA_SYDDNU_added_synced_cols.csv')
        headers=['evstate','logz','teff','mini','mass','logg']
        df=pd.DataFrame(columns=headers)
        ID=sydadd2asf['ID']

        def append_rows(row_data,outputname):
            df=pd.DataFrame(row_data,columns=headers)
            df.to_csv(outputname,mode='a',sep=' ',index=False,header=False)
        #this is just so we dont hvae to re do this if we already have done this, really just for developer purposes
        radius_solar=6.958e10

        mass_solar=1.989e33
        gravconst=6.67259e-8
        #this just for developer testing, not needed, make it **kwargs soon or something so that it is optional maybe?
        if not (data_created):

            for i in range(len(sydadd2asf['ID'])):
                rad_scale=sydadd2asf['radius'][i]
                outputname=(f'{ID[i]}.txt')
                df.to_csv(outputname,index=False,sep=' ')
                grav_solar_log=np.log10((gravconst*mass_solar/(radius_solar*rad_scale)**2))
                #for i in range(len(sydadd2asf['radius'])):
                rowin=[[1,-1.721,5770,1.0,1.0,grav_solar_log]]
                append_rows(rowin,outputname)

            for i in range(len(sydadd2asf['ID'])):
                runstring=f'./asfgrid.py {ID[i]}.txt'
                subprocess.run(runstring,shell=True)

        # WE NEED TO ADD THE HEADER, DO IT HERE SO WE DONT HAVE TO CREATE IT IN IN SYNC
        header_output=['ID','radius','asf_dnu','asf_fdnu']
        rowin=[]
        # use the self.something such that we place it in our regular working directory
        asf_csv=f'{self.wrk_drc}/1M0_0D01.7ALPHA_asfgrid.csv'
        for i in range (0,len(ID)):
            string_txt=f'{ID[i]}.txt.out'
            df=pd.read_csv(string_txt)
            recov_radius=np.sqrt(gravconst*mass_solar/(10**float(df[' logg'])))/radius_solar
            rowin.append([f'{ID[i]}',recov_radius,float(df[' dnu']),float(df[' fdnu'])])
        asfout=pd.DataFrame(rowin,columns=header_output)
        asfout.to_csv(asf_csv,mode='w',index=False,header=True)

        #change back to working directory
        os.chdir(self.wrk_drc)

        # upload the data so update gitignore

    def asf_csv_sync_cols(self):
        column_names_added=['dnu_inf','nu_max','mass','radius','luminosity','mixing_length','metalicity']
        added_radius_tofile=pd.read_csv('1M0_0D01.7ALPHA_asfgrid.csv')
        csv_toextract_radius=pd.read_csv('old_lit_nonad_1M0_0D01.7ALPHA.csv')
        csv_toextract_radius['ID']=csv_toextract_radius['PROFID'].str.replace('profile','').astype(int)
        added_radius_tofile['ID']=added_radius_tofile['ID'].astype(int)
        for x in column_names_added:
            profile_radius_map=dict(zip(csv_toextract_radius['ID'],csv_toextract_radius[x]))
            added_radius_tofile[x]=added_radius_tofile['ID'].map(profile_radius_map)
            added_radius_tofile.to_csv('1M0_0D01.7ALPHA_asfgrid_added_synced_cols.csv',index=False)

    def asfgrid_python_env(self):
        """Return an environment where ASFGrid can import locally installed ebfpy."""
        env=os.environ.copy()
        site_packages=sorted((Path(self.wrk_drc)/'.venv_asfgrid'/'lib').glob('python*/site-packages'))
        if site_packages:
            py_path=str(site_packages[0])
            existing=env.get('PYTHONPATH')
            env['PYTHONPATH']=py_path if not existing else py_path+os.pathsep+existing
        return env

    def li_asfgrid_prepared_data(self, source_csv='li_surfacecorrection_mass_const_l0.csv', evstate=1, run_limit=None):
        """Prepare Li l=0 rows for ASFGrid script-mode input."""
        source_path=Path(self.wrk_drc)/source_csv
        li_data=pd.read_csv(source_path)
        required_columns=['PROFID','mass','radius','luminosity','metalicity','dnu_inf','nu_max','mixing_length','age','fdnu_li_massconst']
        missing=[column for column in required_columns if column not in li_data.columns]
        if missing:
            raise ValueError(f'{source_path} is missing required Li columns: {missing}')

        for column in required_columns:
            li_data[column]=pd.to_numeric(li_data[column],errors='coerce')
        li_data=li_data.replace([np.inf,-np.inf],np.nan).dropna(subset=required_columns).copy()
        li_data=li_data[
            (li_data['mass']>0) &
            (li_data['radius']>0) &
            (li_data['luminosity']>0)
        ].copy()
        if run_limit is not None:
            li_data=li_data.head(int(run_limit)).copy()

        li_data['ID']=li_data.apply(
            lambda row: _li_spectrum_id(row['PROFID'], row['mass'], row['radius'], row['metalicity']),
            axis=1
        )
        li_data['row_index']=np.arange(len(li_data),dtype=int)
        li_data['teff_asf']=SOLAR_TEFF_K*np.power(li_data['luminosity']/np.power(li_data['radius'],2),0.25)
        li_data['logg_asf']=SOLAR_LOGG_CGS+np.log10(li_data['mass']/np.power(li_data['radius'],2))
        li_data['mini_asf']=li_data['mass']
        li_data['evstate_asf']=int(evstate)
        return li_data

    def build_li_asfgrid_manifest(self, source_csv='li_surfacecorrection_mass_const_l0.csv', evstate=1, run_limit=None):
        """Write one ASFGrid script-mode input file plus a manifest for Li l=0 rows."""
        if not self.asfgrid_dir.exists():
            raise FileNotFoundError(f'ASFGrid directory not found: {self.asfgrid_dir}')

        li_data=self.li_asfgrid_prepared_data(source_csv=source_csv, evstate=evstate, run_limit=run_limit)
        input_path=self.asfgrid_dir/self.output_li_asfgrid_l0_input
        output_path=Path(str(input_path)+'.out')
        asf_input=pd.DataFrame({
            'row_index':li_data['row_index'],
            'evstate':li_data['evstate_asf'],
            'feh':li_data['metalicity'],
            'teff':li_data['teff_asf'],
            'mini':li_data['mini_asf'],
            'mass':li_data['mass'],
            'logg':li_data['logg_asf']
        })
        asf_input.to_csv(input_path,sep=' ',index=False)

        manifest=li_data[[
            'row_index','ID','PROFID','radius','mass','luminosity','metalicity',
            'mixing_length','age','dnu_inf','nu_max','fdnu_li_massconst',
            'teff_asf','logg_asf','mini_asf','evstate_asf'
        ]].copy()
        manifest['asf_input_file']=input_path.name
        manifest['asf_output_file']=output_path.name
        manifest_path=Path(self.wrk_drc)/self.output_li_asfgrid_l0_manifest
        manifest.to_csv(manifest_path,index=False)
        print(f'Prepared {len(manifest)} Li ASFGrid rows in {input_path}')
        return manifest_path,input_path,output_path

    def collect_li_asfgrid_csv(self, manifest_path=None, output_path=None):
        """Merge ASFGrid script output back onto the Li l=0 manifest."""
        if manifest_path is None:
            manifest_path=Path(self.wrk_drc)/self.output_li_asfgrid_l0_manifest
        else:
            manifest_path=Path(manifest_path)
        if output_path is None:
            output_path=self.asfgrid_dir/self.output_li_asfgrid_l0_input
            output_path=Path(str(output_path)+'.out')
        else:
            output_path=Path(output_path)

        manifest=pd.read_csv(manifest_path)
        asf_output=pd.read_csv(output_path,sep=r',\s*',engine='python')
        asf_output.columns=[column.strip() for column in asf_output.columns]
        asf_output['row_index']=pd.to_numeric(asf_output['row_index'],errors='coerce').astype('Int64')
        for column in ['dnu','numax','fdnu']:
            asf_output[column]=pd.to_numeric(asf_output[column],errors='coerce')

        result=manifest.merge(
            asf_output[['row_index','dnu','numax','fdnu']],
            on='row_index',
            how='left'
        )
        result=result.rename(columns={
            'dnu':'asf_dnu',
            'numax':'asf_numax',
            'fdnu':'asf_fdnu'
        })
        ok_mask=(
            np.isfinite(result['asf_dnu']) &
            np.isfinite(result['asf_numax']) &
            np.isfinite(result['asf_fdnu']) &
            (result['asf_dnu']>0) &
            (result['asf_fdnu']>0)
        )
        result['asf_status']=np.where(ok_mask,'ok','missing_or_invalid')
        output_csv=Path(self.wrk_drc)/self.output_li_asfgrid_l0
        result.to_csv(output_csv,index=False)
        print(f'Wrote {len(result)} Li ASFGrid rows to {output_csv}')
        print(result['asf_status'].value_counts(dropna=False).to_string())
        return output_csv

    def run_li_asfgrid_l0(self, source_csv='li_surfacecorrection_mass_const_l0.csv', evstate=1, run_limit=None):
        """Run Li l=0 rows through ASFGrid and collect a synchronized CSV."""
        manifest_path,input_path,output_path=self.build_li_asfgrid_manifest(
            source_csv=source_csv,
            evstate=evstate,
            run_limit=run_limit
        )
        result=subprocess.run(
            ['python','asfgrid.py',input_path.name],
            cwd=self.asfgrid_dir,
            env=self.asfgrid_python_env(),
            capture_output=True,
            text=True
        )
        if result.stdout:
            print(result.stdout[-4000:])
        if result.stderr:
            print(result.stderr[-4000:])
        if result.returncode != 0:
            raise RuntimeError(f'ASFGrid failed with exit code {result.returncode}')
        return self.collect_li_asfgrid_csv(manifest_path,output_path)

    def plot_li_asfgrid_fdnu(self, asf_csv=None):
        """Plot Li-through-ASFGrid f_Delta nu versus radius."""
        import matplotlib.pyplot as plt
        if asf_csv is None:
            asf_csv=Path(self.wrk_drc)/self.output_li_asfgrid_l0
        else:
            asf_csv=Path(asf_csv)
        data=pd.read_csv(asf_csv).replace([np.inf,-np.inf],np.nan)
        data=data.dropna(subset=['radius','asf_fdnu'])
        data=data[
            (data['radius']>0) &
            (data['asf_fdnu']>0) &
            (data['asf_status']=='ok')
        ].copy()
        data['asf_fdnu_li_convention']=1.0/data['asf_fdnu']

        self.path_li_asfgrid_fdnu_plot.parent.mkdir(parents=True,exist_ok=True)
        fig,ax=plt.subplots(figsize=(10.5,6.3),dpi=300)

        def binned_summary(frame, value_column, bin_count=45, min_count=8):
            clean=frame[['radius',value_column]].replace([np.inf,-np.inf],np.nan).dropna().copy()
            clean=clean[(clean['radius']>0) & (clean[value_column]>0)]
            if clean.empty:
                return pd.DataFrame()
            bins=np.geomspace(clean['radius'].min(),clean['radius'].max(),bin_count+1)
            clean['radius_bin']=pd.cut(clean['radius'],bins=bins,include_lowest=True)
            summary=clean.groupby('radius_bin',observed=True).agg(
                radius_median=('radius','median'),
                fdnu_median=(value_column,'median'),
                fdnu_q16=(value_column,lambda values: values.quantile(0.16)),
                fdnu_q84=(value_column,lambda values: values.quantile(0.84)),
                count=(value_column,'size')
            ).reset_index()
            return summary[summary['count']>=min_count].copy()

        if 'fdnu_li_massconst' in data.columns:
            li_reference=data.dropna(subset=['fdnu_li_massconst'])
            ax.scatter(
                li_reference['radius'],
                li_reference['fdnu_li_massconst'],
                color='#9aa4b2',
                alpha=0.16,
                s=4,
                edgecolors='none',
                label='Li reference rows'
            )
            li_summary=binned_summary(li_reference,'fdnu_li_massconst')
            if not li_summary.empty:
                ax.plot(li_summary['radius_median'],li_summary['fdnu_median'],color='#6a8cff',linewidth=1.9,label='Li binned median')
                ax.fill_between(
                    li_summary['radius_median'].to_numpy(dtype=float),
                    li_summary['fdnu_q16'].to_numpy(dtype=float),
                    li_summary['fdnu_q84'].to_numpy(dtype=float),
                    color='#6a8cff',
                    alpha=0.12,
                    linewidth=0
                )
        ax.scatter(
            data['radius'],
            data['asf_fdnu_li_convention'],
            color='#31a354',
            alpha=0.16,
            s=4,
            edgecolors='none',
            label=r'ASFGrid rows, inverted convention'
        )
        asf_summary=binned_summary(data,'asf_fdnu_li_convention')
        if not asf_summary.empty:
            ax.plot(asf_summary['radius_median'],asf_summary['fdnu_median'],color='#238b45',linewidth=1.9,label='ASFGrid binned median')
            ax.fill_between(
                asf_summary['radius_median'].to_numpy(dtype=float),
                asf_summary['fdnu_q16'].to_numpy(dtype=float),
                asf_summary['fdnu_q84'].to_numpy(dtype=float),
                color='#31a354',
                alpha=0.12,
                linewidth=0
            )
        ax.grid(True,color='#d7dce2',linestyle='--',linewidth=0.7,alpha=0.85)
        ax.set_xscale('log')
        ax.set_xlabel(r'Radius ($R_\odot$)',fontsize=12)
        ax.set_ylabel(r'$f_{\Delta \nu}$',fontsize=12)
        ax.set_title('Li l=0 sample processed through ASFGrid',fontsize=14,pad=10)
        ax.text(
            0.03,
            0.97,
            f'N = {len(data)}\nASFGrid ok rows\ninverted convention',
            transform=ax.transAxes,
            ha='left',
            va='top',
            fontsize=8.7,
            bbox=dict(boxstyle='round',facecolor='white',edgecolor='#9aa4b2',alpha=0.92)
        )
        ax.legend(loc='best',framealpha=0.9,fontsize=8)
        ax.tick_params(axis='both',which='major',labelsize=10)
        plt.tight_layout()
        plt.savefig(self.path_li_asfgrid_fdnu_plot,bbox_inches='tight',facecolor='white',dpi=300)
        plt.close()
        print(f'Wrote Li ASFGrid fdnu plot to {self.path_li_asfgrid_fdnu_plot}')
        return self.path_li_asfgrid_fdnu_plot

    def pySYD_calls(self,standard=False,li=False,ells=0,cores_mt=1):
        if pySYD_one_shot is None:
            raise ImportError('pySYD_one_shot could not be imported in this Python environment.')

        pysyd_obj=pySYD_one_shot.pySYD_oneshot()

        if standard == False:
            # remember to change belwo fia nythign happens

            pysyd_obj.power_spec_filename_correction_nested(self.nonad_create_dir)
            pysyd_obj.full_pysyd_nested(self.wrk_drc,self.nonad_create_dir,self.output_std)
            # below is just example string it does not have to be the csv as placed
            pysyd_obj.delete_repeats('li_surfacecorrection_mass_const_l0.csv','li_surfacecorrection_mass_const_l0_cleaned.csv')
            pysyd_obj.sync_pysyd_csv_nested()

        if li == False and ells==0:

            pysyd_obj.power_spec_filename_correction_flat(self.parquet_runs_dir_l0)

            pysyd_obj.full_pysyd_flat(self.wrk_drc,self.parquet_runs_dir_l0,self.output_li_l0,cores_mt)
            #how can we go about this, there is alot that we need to do like what is the problem 
            pysyd_obj.delete_repeats(self.output_li_l0,self.cleaned_li_l0)
            pysyd_obj.sync_pysyd_csv_flat('li_surfacecorrection_mass_const_l0.csv',self.cleaned_li_l0,self.synced_cleaned_li_l0)

    def li_BAM_chunk(self,sng_file):
        summary = self.spectrum_summary(sng_file)
        print('Initial NUMAX guess:', summary['numax_guess'])
        return summary['numax_guess']


    def create_BAM_data_li(self,num_chunks=1, run_limit=None, bam_executable=None, mc_iter=500):
        manifest = self.build_li_bam_l0_manifest(mc_iter=mc_iter)
        print(f'Prepared {len(manifest)} Li l0 BAM spectra in {self.output_li_bam_l0_manifest}')
        return self.run_bam_manifest(
            self.output_li_bam_l0_manifest,
            self.output_li_bam_l0_corrected,
            cores_mt=num_chunks,
            bam_executable=bam_executable,
            run_limit=run_limit,
            enable_plots=False
        )



    def surface_correction_li_mass_const_chunk(self,chunk):

        # CHANGE BELOW WHEN DOING MULTIPLE DIFFERENT ELLS OR JUST ONE AT A TIME
        #ells=0

        base_mask= (
            # only solar mettalicity stars
            (chunk['star_mass'] > 0.99) & (chunk['star_mass'] < 1.005) &

            # need them to exist properly
            ~chunk['Dnu_freq'].isna() &
            ~chunk['numax'].isna() &
            ~chunk['star_mass'].isna() &
            ~chunk['radius'].isna() &
            ~chunk['luminosity'].isna() &
            ~chunk['model_number'].isna() &
            ~chunk['massinit'].isna() &
            ~chunk['amlt'].isna() &
            ~chunk['fehinit'].isna() &
            ~chunk['star_age'].isna() &
            ~chunk['E_norm'].isna() &

            chunk['l'].apply(lambda x: isinstance(x,(list,np.ndarray))and len(x)>0) &
            chunk['fc'].apply(lambda x: isinstance(x,(list,np.ndarray)) and len(x)>0) &


            # make sure it is all positive (it should be)
            (chunk['Dnu_freq']) > 0  &
            (chunk['numax'] > 0) &
            (chunk['star_mass'] > 0) &
            (chunk['radius'] > 0) &
            (chunk['luminosity'] > 0) &
            (chunk['massinit'] > 0) &
            (chunk['amlt'] > 0) &
            (chunk['star_age'] > 0)
            # fehinit can be negative (metallicity)

        )

        filtered_chunk=chunk[base_mask].copy()

        if len(filtered_chunk)==0:
            return pd.DataFrame()

        l_single_list=[]
        fc_single_list=[]
        enorm_single_list=[]
        rows_to_keep=[]

        for idx, row in filtered_chunk.iterrows():

            l_array=row['l']
            fc_array=row['fc']
            enorm_array=row['E_norm']

            if len(l_array)!=len(fc_array):
                continue

            l0_values=[]
            fc0_values=[]
            enorm0_values=[]
            # how do we do the enorm0 as well?
            for l_val, fc_val, enorm_val in zip (l_array,fc_array,enorm_array):
                if l_val==0:
                    l0_values.append(l_val) #this should always be 0, it is just the single l=0 modes that we want
                    fc0_values.append(fc_val)
                    enorm0_values.append(enorm_val)

            if l0_values:
                rows_to_keep.append(row)
                l_single_list.append(l0_values)
                fc_single_list.append(fc0_values)
                enorm_single_list.append(enorm0_values)


        if not rows_to_keep:
            return pd.DataFrame()


        valid_chunk=pd.DataFrame(rows_to_keep).reset_index(drop=True)

        valid_chunk['l_single']=l_single_list
        valid_chunk['fc_single']=fc_single_list
        valid_chunk['enorm_single']=enorm_single_list

        #print(valid_chunk['l_single'])
        #print(valid_chunk['fc_single'])
        #return
        #print(valid_chunk['star_age'])
        valid_chunk['numax_raw']=(valid_chunk['star_mass'] *valid_chunk['radius']**(-7./4.)*
                    valid_chunk['luminosity']**(-1./8.) * numax_sun)

        valid_chunk['dnu_inf']= (np.sqrt(valid_chunk['star_mass']/
                    (valid_chunk['radius']**3)) * dnu_sun)
        # retry the dnu_theory from nonad, make sure that it works or something
        # plot the fdnu for the solar case
        # SOLAR MASS CASE SOLAR MASS CASR SOLAR MASS CASE SOLAR MASS CASE DO ALL OF THAT
        # compare fdnus from here, and then the dnu_theory


        #valid_chunk['l_series'] = valid_chunk['l'].apply(lambda x: pd.Series(x) if x is not None else pd.Series([],dtype=float))
        #valid_chunk['fc_series'] = valid_chunk['fc'].apply(lambda x: pd.Series(x) if x is not None else pd.Series([],dtype=float))


        valid_chunk['dnu_thry_li']=valid_chunk.apply(lambda row: dnu_theory(row['l_single'],row['fc_single'],row['numax_raw']),axis=1)
        valid_chunk['pysyd_id']=valid_chunk.apply(
            lambda row: _li_spectrum_id(
                row['model_number'],
                row['massinit'],
                row['radius'],
                row['fehinit']
            ),
            axis=1
        )

        chunk_results = pd.DataFrame({
            'pysyd_id': valid_chunk['pysyd_id'],
            'PROFID': valid_chunk['model_number'],
            'dnu_theory': valid_chunk['dnu_thry_li'],
            'numax_raw': valid_chunk['numax_raw'],
            'dnu_inf': valid_chunk['dnu_inf'],
            'nu_max': valid_chunk['numax'],
            'mass': valid_chunk['massinit'],
            'radius': valid_chunk['radius'],
            'luminosity': valid_chunk['luminosity'],
            'mixing_length': valid_chunk['amlt'],
            'metalicity': valid_chunk['fehinit'],
            'age': valid_chunk['star_age']
        })
        t = 4*math.pi*10
        nm = 40
        # based on what we did before gamma should be 0.12 so width should be changed to 0.0036, dont know where.
        # based on this : https://iopscience.iop.org/article/10.3847/1538-4357/835/2/172/pdf
        width = 0.0036
        #width = nm*.005
        dF = 1/(2*t)

        for idx, row in valid_chunk.iterrows():
            fc_values=np.array(row['fc_single'])
            E_norm_values=np.array(row['enorm_single'])

            f1=[]


            for i in np.arange(len(fc_values)-1):
                flow= fc_values[i]
                fhi = fc_values[i+1]

                mask_range=(fc_values > flow) & (fc_values<fhi)
                if np.any(mask_range):
                    sub_fc=fc_values[mask_range]
                    sub_E_norm=E_norm_values[mask_range]
                    jfl1=sub_fc[np.argmax(1/sub_E_norm)]
                    f1.append(jfl1)

            f = np.arange(0.0,fc_values[-1].real, dF)

            ampl0 = f*0
            ampl2 = f*0
            ampl1 = f*0
            a_l0 = 1
            a_l1 = 0.5
            a_l2 = 0.75

            #only reals when dealing with nonad
            for loc in fc_values:
                ampl0 += lorentzian(f, loc.real, a_l0, width)

            current_numax = row['numax']

            amp_total = ampl0 + ampl1 + ampl2
            b = 0.25*current_numax
            ampnm = np.exp(- (f - current_numax)**2/(2.0*b**2))
            ofac = 1
            noise = (chi2.rvs(2*ofac, size=len(ampnm)) * amp_total*ampnm)/2.0/ofac

            profile_indiv_data = pd.DataFrame({
                    'frequency': f,
                    'amplitude': amp_total*ampnm
            })

            if not os.path.exists(self.parquet_runs_dir_l0):
                os.makedirs(self.parquet_runs_dir_l0)

            # BELOW IS STRUCTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # we can make the strings such that it is (f'profile{valid_}_mass{number}_radius{number}i_metallicity{number}')
            filename = f"{row['pysyd_id']}.csv"
            filepath = Path(self.parquet_runs_dir_l0)/filename

            # do not want the headers or else pySYD cant read each header
            profile_indiv_data.to_csv(filepath,index=False,header=False)

        return chunk_results


    def surface_correction_li_mass_const(self,num_chunks=15,created_csv_and_data=True):

        # data as taken from special saved disk of https://ui.adsabs.harvard.edu/abs/2025AAS...24535303C/abstract
        if created_csv_and_data == False:
            arrs=pd.read_parquet(f'{self.parquet_files}/reimers_index_0_8191_rg_modes.parquet')

            #print(arrs.columns.tolist())

            realized_chunks=np.array_split(arrs,num_chunks)
            with ThreadPool(processes=num_chunks) as pool:
                results_chunks=pool.map(self.surface_correction_li_mass_const_chunk,realized_chunks)

            all_results=pd.concat(results_chunks,ignore_index=True)
            all_results['fdnu_li_massconst']=all_results['dnu_inf']/all_results['dnu_theory']
            all_results.to_csv('li_surfacecorrection_mass_const_l0.csv',index=False)



    def li_surf_interp(self):
        data=pd.read_csv(f'{self.wrk_drc}/li_data_rawcalculated.csv')
        dat_len=len(data)
        #interp_check=interp1d(data['metalicity'],data['radius'])
        #interp=RegularGridInterpolator((data['mass'],data['age'],data['metalicity'],data['mixing_length']),data['radius'],method='linear',fill_value=None,bounds_error=False)
        #print(interp)
        interp_check=LinearNDInterpolator(list(zip(data['mass'], data['age']/10E9,data['metalicity'],data['mixing_length'])),data['radius'])

        #possible_mixing=np.linspace(0.1,4.0,dat_len)
        #fixed_mass=np.ones(dat_len)
        #fixed_age=4*np.ones(dat_len)
        #fixed_metal=np.zeros(dat_len)

        #radius_interp_tone=interp_check(list(zip(fixed_mass,fixed_age,fixed_metal,possible_mixing)))
        #print(radius_interp_tone)
        #points=np.column_stack([
        #    data['mass'],
        #    data['age']/1e9,
        #    data['metalicity'],
        #    data['mixing_length']
        #    ])
        #
        #interpolator=NearestNDInterpolator(points,data['radius'])
        #
        #new_mixing= np.arange(0.1,4,0.001)
        #fixed_mass=1.0
        #fixed_age=4.0
        #fixed_metal=0.0
        #target_radius=1.0


        #query_points=np.column_stack([
        #    np.full_like(new_mixing, fixed_mass),
        #    np.full_like(new_mixing, fixed_age),
        #    np.full_like(new_mixing, fixed_metal),
        #    new_mixing
        #])

        #li_sfc_fdnu_interp=interpolator(query_points)
        #idx_closest=np.argmin(np.abs(li_sfc_fdnu_interp-target_radius))
        #best_mixing=new_mixing[idx_closest]
        #achieved_radius=li_sfc_fdnu_interp[idx_closest]



        #print(type(li_sfc_fdnu_interp))
        #with np.printoptions(threshold=np.inf, linewidth=200):
        #    print(best_mixing, achieved_radius)


        # above calcualtes the best mixing to be 2.3040



        #interpolator = LinearNDInterpolator(np.column_stack([data['mass'],data['age']/1e9,data['metalicity'],data['mixing_length']]),data['fdnu_li'],rescale=True)
        ##li_sfc_fdnu_interp=interpolator(1, 4E9, 0,new_mixing)
        #print(li_sfc_fdnu_interp)
        # given these interpolatoins, what fdnu?
        # we can also do that what radius can we receive back from it like in a reverse way?


        #li_sfc_fdnu_interp=interpolator(whatever numax, mass, metalicity, whatever)


        #interpolator = LinearNDInterpolator([data['nu_max'],data['metalicity'],data['mass'],data['mixing_length']],data['radius'])
        # so the logic is something like [M,t,[FeH],\alpha] -> Radius, find this such that we recreate the radius with varying the mixing length, keeping everything else the same I think. we want to get one. use 4 billion years old. Then once we do that, the alpha mixing length that we find from that we then want to keep that metalicity the same, but vary the others I think? ask joel


    def li_surf_interp_pt2(self):
        data=pd.read_csv(f'{self.wrk_drc}/li_data_rawcalculated.csv')
        dat_len=len(data)
        #interp_check=interp1d(data['metalicity'],data['radius'])
        #interp=RegularGridInterpolator((data['mass'],data['age'],data['metalicity'],data['mixing_length']),data['radius'],method='linear',fill_value=None,bounds_error=False)
        #print(interp)
        interp_check=LinearNDInterpolator(list(zip(data['mass'], data['age']/10E9,data['metalicity'],data['mixing_length'])),data['radius'])

        #possible_mixing=np.linspace(0.1,4.0,dat_len)
        #fixed_mass=np.ones(dat_len)
        #fixed_age=4*np.ones(dat_len)
        #fixed_metal=np.zeros(dat_len)

        #radius_interp_tone=interp_check(list(zip(fixed_mass,fixed_age,fixed_metal,possible_mixing)))
        #print(radius_interp_tone)
        #points=np.column_stack([
        #    data['mass'],
        #    data['age']/1e9,
        #    data['metalicity'],
        #    data['mixing_length']
        #    ])
        #
        #interpolator=NearestNDInterpolator(points,data['radius'])
        #
        #new_mixing= np.arange(0.1,4,0.001)
        #fixed_mass=1.0
        #fixed_age=4.0
        #fixed_metal=0.0
        #target_radius=1.0


        #query_points=np.column_stack([
        #    np.full_like(new_mixing, fixed_mass),
        #    np.full_like(new_mixing, fixed_age),
        #    np.full_like(new_mixing, fixed_metal),
        #    new_mixing
        #])

        #li_sfc_fdnu_interp=interpolator(query_points)
        #idx_closest=np.argmin(np.abs(li_sfc_fdnu_interp-target_radius))
        #best_mixing=new_mixing[idx_closest]
        #achieved_radius=li_sfc_fdnu_interp[idx_closest]



        #print(type(li_sfc_fdnu_interp))
        #with np.printoptions(threshold=np.inf, linewidth=200):
        #    print(best_mixing, achieved_radius)


        # above calcualtes the best mixing to be 2.3040



        #interpolator = LinearNDInterpolator(np.column_stack([data['mass'],data['age']/1e9,data['metalicity'],data['mixing_length']]),data['fdnu_li'],rescale=True)
        ##li_sfc_fdnu_interp=interpolator(1, 4E9, 0,new_mixing)
        #print(li_sfc_fdnu_interp)
        # given these interpolatoins, what fdnu?
        # we can also do that what radius can we receive back from it like in a reverse way?


        #li_sfc_fdnu_interp=interpolator(whatever numax, mass, metalicity, whatever)


        #interpolator = LinearNDInterpolator([data['nu_max'],data['metalicity'],data['mass'],data['mixing_length']],data['radius'])
        # so the logic is something like [M,t,[FeH],\alpha] -> Radius, find this such that we recreate the radius with varying the mixing length, keeping everything else the same I think. we want to get one. use 4 billion years old. Then once we do that, the alpha mixing length that we find from that we then want to keep that metalicity the same, but vary the others I think? ask joel

def main():
    start_time = time.time()
    true_oneshot = full_create()

    #true_oneshot.make_oldlit_csv()
    #pysyd has to happen first for data synchronization
    #true_oneshot.pysyd_calls(False,True)
    #true_oneshot.create_BAM_data()
    #true_oneshot.create_bam_csv()
    #true_oneshot.sync_bam_csvformat()
    #true_oneshot.asf_grid_datacreate(data_created=False)
    #true_oneshot.asf_csv_sync_cols()



    # as of 11.27.2025 below is true, we have the data the way we want it
    # make and ell=0 variable and separate out which ones we want ot look at in below method
    # for now it is just ells=0
    # FILES MUST END IN _PS.txt
    #true_oneshot.surface_correction_li_mass_const(16, False)
    #true_oneshot.pySYD_calls(True, False,0, 16)
    #true_oneshot.pySYD_calls(False,True,0, 16)

    #now lets get BAM to work
    #true_oneshot.create_BAM_data_li(16)
    end_time = time.time()
    print(f'The script took {end_time-start_time} seconds to run.')


# MAYBE MAIN OUR OWN INTEGER GET METHOD? ACTUALLY PROBABLY JUST DO IT IN INIT
# the preprocessw is required at all timesONLY IF DOING FULL BAM
# need fitbhsp_syd_dnu only in bam_runs because that is where we set up the enviornment



# 4.2.26 -> let's do the original fdnu of just li, well see about the numax
# also I think that the NUMAX is what is changing them
# for the fdnu for li lets make sure we filter by only solar mass, solar metallilcity
if __name__ == '__main__':
    main()
