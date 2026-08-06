import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from matplotlib.ticker import FormatStrFormatter

try:
	from astropy.io import ascii
except ModuleNotFoundError:
	ascii = None


PLOT_DPI=300
PLOT_FIGSIZE=(10.5,6.3)
PLOT_COLORS={
	'model':'#2f2f2f',
	'pysyd':'#2c7fb8',
	'bam':'#7b3294',
	'asfgrid':'#31a354',
	'li':'#6a8cff',
	'apokasc':'#b44cc2',
	'warning':'#c81d4f'
}


def _read_optional_csv(path, **kwargs):
	path = Path(path).expanduser()
	if path.exists():
		return pd.read_csv(path, **kwargs)
	return pd.DataFrame()


def _read_optional_parquet(path, **kwargs):
	path = Path(path).expanduser()
	if path.exists():
		return pd.read_parquet(path, **kwargs)
	return pd.DataFrame()


def _read_apokasc_data():
	for path in (Path('apoxeradiusinterp.txt'), Path('apokasc.txt')):
		if not path.exists():
			continue
		if ascii is not None:
			try:
				return ascii.read(path).to_pandas()
			except Exception:
				pass
		data = pd.read_csv(path, sep=r'\s+', header=None)
		data.columns = [f'col{i}' for i in range(1, len(data.columns) + 1)]
		return data
	raise FileNotFoundError('Expected apoxeradiusinterp.txt or apokasc.txt for APOKASC radii.')


def _finite_plot_data(radius, value, radius_name='radius', value_name='value', positive_radius=True):
	data=pd.DataFrame({
		radius_name:radius,
		value_name:value
	}).replace([np.inf,-np.inf],np.nan).dropna()
	if positive_radius:
		data=data[data[radius_name]>0].copy()
	data=data[np.isfinite(data[value_name])].copy()
	return data.sort_values(radius_name).reset_index(drop=True)


def _binned_radius_summary(radius, value, bins=38, min_count=3):
	data=_finite_plot_data(radius,value)
	if data.empty:
		return pd.DataFrame()
	radius_min=data['radius'].min()
	radius_max=data['radius'].max()
	if not np.isfinite(radius_min) or not np.isfinite(radius_max) or radius_min<=0 or radius_min==radius_max:
		return pd.DataFrame()
	data['radius_bin']=pd.cut(
		data['radius'],
		bins=np.geomspace(radius_min,radius_max,bins+1),
		include_lowest=True
	)
	summary=data.groupby('radius_bin',observed=True).agg(
		radius_median=('radius','median'),
		value_median=('value','median'),
		value_q16=('value',lambda values: values.quantile(0.16)),
		value_q84=('value',lambda values: values.quantile(0.84)),
		count=('value','size')
	).reset_index()
	return summary[summary['count']>=min_count].copy()


def _plot_binned_radius_summary(ax, radius, value, color, label, bins=38, min_count=3, linewidth=1.8, alpha=0.18):
	summary=_binned_radius_summary(radius,value,bins=bins,min_count=min_count)
	if summary.empty:
		return summary
	x=summary['radius_median'].to_numpy(dtype=float)
	y=summary['value_median'].to_numpy(dtype=float)
	q16=summary['value_q16'].to_numpy(dtype=float)
	q84=summary['value_q84'].to_numpy(dtype=float)
	ax.plot(x,y,color=color,linewidth=linewidth,label=label,zorder=4)
	ax.fill_between(x,q16,q84,color=color,alpha=alpha,linewidth=0,zorder=2)
	return summary


def _style_radius_plot(ax, ylabel=r'$f_{\Delta\nu}$', title=None, xlabel=r'Radius ($R_\odot$)', xscale='log'):
	if xscale:
		ax.set_xscale(xscale)
	ax.set_xlabel(xlabel,fontsize=12)
	ax.set_ylabel(ylabel,fontsize=12)
	if title:
		ax.set_title(title,fontsize=14,pad=10)
	ax.grid(True,color='#d7dce2',linestyle='--',linewidth=0.7,alpha=0.85)
	ax.tick_params(axis='both',which='major',labelsize=10)
	for spine in ['top','right']:
		ax.spines[spine].set_alpha(0.35)


def _add_stats_box(ax, lines, loc='upper left'):
	positions={
		'upper left':(0.03,0.97,'left','top'),
		'upper right':(0.97,0.97,'right','top'),
		'lower left':(0.03,0.03,'left','bottom'),
		'lower right':(0.97,0.03,'right','bottom')
	}
	x,y,ha,va=positions.get(loc,positions['upper left'])
	ax.text(
		x,
		y,
		'\n'.join(lines),
		transform=ax.transAxes,
		ha=ha,
		va=va,
		fontsize=8.7,
		bbox=dict(boxstyle='round',facecolor='white',edgecolor='#9aa4b2',alpha=0.92)
	)


class all_plots:
	def __init__(self,**kwargs):
		self.repo_dir=Path(kwargs.get('wrk_drc', Path(__file__).resolve().parent)).expanduser()
		self.plot_dir=self.repo_dir/'plots_directory'
		li_raw_path=self.repo_dir/'parquet_files'/'reimers_index_0_8191_rg_modes.parquet'
		li_raw_columns=[
			'model_number','star_mass','massinit','radius','numax','Dnu_freq','Dnu_int',
			'l','fc','E_norm','amlt','fehinit','star_age','luminosity'
		]
		li_raw=_read_optional_parquet(li_raw_path, columns=li_raw_columns) if kwargs.get('load_li_raw', False) else pd.DataFrame()
		self.data={
			'old_lit':pd.read_csv('old_lit_nonad_1M0_0D01.7ALPHA.csv'),
				'pySYD':pd.read_csv('1M0_0D01.7ALPHA_SYDDNU_added_synced_cols.csv'),
				'BAM':pd.read_csv('1M0Z_BAM_added_synced_cols.csv'),
				#change below to the nonad calculation
				'asf_raw_2022':_read_optional_csv(self.repo_dir/'asfgrid_code'/'asfgrid_v0.0.6'/'asfgrid_data_aug2022'/'dnu_grid_evstate.txt',sep=', ',engine='python'),
				'asfgrid':pd.read_csv('1M0_0D01.7ALPHA_asfgrid_added_synced_cols.csv'),
				'apokasc_data':_read_apokasc_data(),
				'li_raw': li_raw,
				'li_data_nonad': pd.read_csv('1M_li_added_synced_cols.csv'),
				'li_data_l0_nonad':pd.read_csv('1M_li_added_synced_cols_l0.csv'),
				'li_asfgrid_l0':_read_optional_csv(self.repo_dir/'1M_li_asfgrid_l0.csv')
				}
		if not os.path.exists(self.plot_dir):
			os.mkdir(self.plot_dir)
		else: 	
			print(f'Directory already created, named {self.plot_dir}')
			pass



	@property
	def radius_asf_raw(self):
		return self.data['asf_raw_2022']['sradius']
	@property
	def fdnu_asf_raw(self):
		return self.data['asf_raw_2022']['dnu_sc']/self.data['asf_raw_2022']['dnu_frq']
	@property
	def mass_asf_raw(self):
		return self.data['asf_raw_2022']['mass']


	@property
	def radius_li_nonad(self):
		return self.data['li_data_nonad']['radius']
	@property
	def fdnu_li_nonad(self):
		return self.data['li_data_nonad']['fdnu_li_massconst']
	@property
	def mass_li_nonad(self):
		return self.data['li_data_nonad']['mass']
	@property
	def amlt_li_nonad(self):
		return self.data['li_data_nonad']['mixing_length']

	@property
	def fdnu_li_l0(self):
		return self.data['li_data_l0_nonad']['fdnu_li_massconst']
	@property
	def fdnu_li_pySYD_l0(self):
		return self.data['li_data_l0_nonad']['dnu_inf']/self.data['li_data_l0_nonad']['dnu_syd']
	@property
	def radius_li_l0(self):
		return self.data['li_data_l0_nonad']['radius']
	@property
	def radius_li_log_l0(self):
		return np.log(self.data['li_data_l0_nonad']['radius'])

	@property
	def fdnu_li(self):
		return self.data['li_data_nonad']['fdnu_li_massconst']
	@property
	def fdnu_li_pySYD(self):
		return self.data['li_data_nonad']['dnu_inf']/self.data['li_data_nonad']['dnu_syd']
	@property
	def radius_li(self):
		return self.data['li_data_nonad']['radius']
	@property
	def radius_li_log(self):
		return np.log(self.data['li_data_nonad']['radius'])


	@property
	def old_lit_fdnu(self):
		return self.data['old_lit']['fdnustdlit']
	@property
	def old_lit_radius(self):
		return self.data['old_lit']['radius']
	@property
	def old_lit_radius_solar(self):
		return np.log(self.data['old_lit']['radius'])




	@property
	def fdnu_pySYD(self):
		return self.data['pySYD']['dnu_inf'] / self.data['pySYD']['dnu_syd']
	@property
	def radius_pySYD(self):
		return self.data['pySYD']['radius']
	@property
	def radius_solar_log_pySYD(self):
		return np.log(self.data['pySYD']['radius'])



	@property
	def fdnu_BAM(self):
		return self.data['BAM']['dnu_inf'] / self.data['BAM']['dnu_BAM']
	@property
	def radius_BAM(self):
		return self.data['BAM']['radius']
	@property
	def radius_solar_log_BAM(self):
		return np.log(self.data['BAM']['radius'])



	@property
	def fdnu_asfgrid(self):
		return self.data['asfgrid']['asf_fdnu'] 
	@property
	def radius_asfgrid(self):
		return self.data['asfgrid']['radius']
	@property
	def radius_solar_log_asfgrid(self):
		return np.log(self.data['asfgrid']['radius'])



	@property
	def radius_apokasc(self):
		return self.data['apokasc_data']['col19']
		#the radius of real_radius is based on the (stddev/rad)^2=..._(2*fdnudev/fdnu)
	@property
	def radius_err_apokasc(self):
		return self.data['apokasc_data']['col20']
	@property
	def fdnu_apokasc(self):
		return self.data['apokasc_data']['col14']
	@property
	def fdnu_apokasc_err(self):
		return self.data['apokasc_data']['col15']
	@property
	def mass_apokasc(self):
		return self.data['apokasc_data']['col17']
	@property
	def mass_apokasc_err(self):
		return self.data['apokasc_data']['col18']


	@property
	def path_oldlit(self):
		return self.plot_dir/'1M0_0D01.7ALPHA_stdlit.png'
	@property
	def path_pySYD_plot(self):
		return self.plot_dir/'1M0_0D01.7ALPHA_pySYD.png'
	@property
	def path_BAM_plot(self):
		return self.plot_dir/'1M0_0D01.7ALPHA_BAM.png'
	@property
	def path_asfgrid_plot(self):
		return self.plot_dir/'1M0_0D01.7ALPHA_asfgrid.png'
	@property
	def path_pySYDinterp_plot(self):
		return self.plot_dir/'Interpolated_pySYD_fdnu.png'
	@property
	def path_asfgridinterp_plot(self):
		return self.plot_dir/'Interpolated_asfgrid_fdnu.png'
	@property
	def path_pySYDdnuvsBAMdnu_plot(self):
		return self.plot_dir/'1M0_0D01.7ALPHA_syddnuvsbamdnu_frac.png'
	@property
	def path_li_fdnu_vs_pySYD_plot(self):
		return self.plot_dir/'fdnu_v_rad_li_pySYD.png'
	@property
	def path_li_fdnu(self):
		return self.plot_dir/'fdnu_li.png'
	@property
	def path_li_fdnu_l0(self):
		return self.plot_dir/'fdnu_li_l0.png'
	@property
	def path_li_fdnu_vs_pySYD_l0_plot(self):
		return self.plot_dir/'fdnu_v_rad_li_pySYD_l0.png'
	@property
	def path_li_fdnu_fractional_comp(self):
		return self.plot_dir/'fdnu_frac_compar_l0123.png'
	@property
	def path_li_vs_asfgrid(self):
		return self.plot_dir/'li_vs_asfgrid_fdnu.png'
	@property
	def path_li_vs_apokasc(self):
		return self.plot_dir/'li_vs_apokasc_fdnu.png'
	@property
	def path_li_nonad_color(self):
		return self.plot_dir/'li_nonad_amlt_colored.png'
	@property
	def path_unified_fdnu_jump_summary(self):
		return self.plot_dir/'unified_fdnu_jump_summary.png'
	@property
	def path_unified_fdnu_jump_table(self):
		return self.plot_dir/'unified_fdnu_jump_diagnostics.csv'
	def path_unified_fdnu_jump_case(self, profile_id):
		return self.plot_dir/f'unified_fdnu_jump_profile{int(profile_id)}.png'
	@property
	def path_nonad_pysyd_internal_table(self):
		return self.plot_dir/'nonad_pysyd_internal_diagnostics.csv'
	@property
	def path_nonad_pysyd_internal_summary(self):
		return self.plot_dir/'nonad_pysyd_internal_summary.png'
	@property
	def path_nonad_reference_long_table(self):
		return self.plot_dir/'nonad_reference_comparison_long.csv'
	@property
	def path_nonad_reference_comparison_table(self):
		return self.plot_dir/'nonad_reference_comparison_binned.csv'
	@property
	def path_nonad_reference_comparison_summary(self):
		return self.plot_dir/'nonad_reference_comparison_summary.png'
	def path_nonad_pysyd_internal_case(self, profile_id):
		return self.plot_dir/f'nonad_pysyd_internal_profile{int(profile_id)}.png'
	def path_pysyd_cutoff(self, profile_id):
		return self.plot_dir/f'pysyd_cutoff_profile{int(profile_id)}.png'

	def apokasc_radius_error_data(self):
		apokasc_data=pd.DataFrame({
			'radius':self.radius_apokasc,
			'radius_err':self.radius_err_apokasc
		}).replace([np.inf,-np.inf],np.nan).dropna()
		apokasc_data=apokasc_data[
			(apokasc_data['radius']>0) &
			(apokasc_data['radius_err']>=0)
		].copy()
		apokasc_data['frac_radius_err']=apokasc_data['radius_err']/apokasc_data['radius']
		return apokasc_data

	def add_interpolated_fdnu_uncertainty(self, apokasc_data, interp_func):
		lower_radius=np.maximum(
			apokasc_data['radius']-apokasc_data['radius_err'],
			np.finfo(float).tiny
		)
		upper_radius=apokasc_data['radius']+apokasc_data['radius_err']
		apokasc_data['fdnu_interp_lower']=interp_func(lower_radius)
		apokasc_data['fdnu_interp_upper']=interp_func(upper_radius)
		fdnu_envelope=apokasc_data[['fdnu_interp_lower','fdnu_interp','fdnu_interp_upper']]
		apokasc_data['fdnu_err_low']=apokasc_data['fdnu_interp']-fdnu_envelope.min(axis=1)
		apokasc_data['fdnu_err_high']=fdnu_envelope.max(axis=1)-apokasc_data['fdnu_interp']
		apokasc_data['fdnu_err_centered']=apokasc_data[['fdnu_err_low','fdnu_err_high']].max(axis=1)
		apokasc_data['fdnu_err_max']=apokasc_data[['fdnu_err_low','fdnu_err_high']].max(axis=1)
		return apokasc_data

	def interpolated_fdnu_errorbar_mask(self, apokasc_data, base_mask, upper_quantile=0.95):
		fdnu_err=apokasc_data.loc[base_mask,'fdnu_err_max']
		if fdnu_err.empty:
			return base_mask
		cutoff=fdnu_err.quantile(upper_quantile)
		if not np.isfinite(cutoff):
			return base_mask
		return base_mask & (apokasc_data['fdnu_err_max']<=cutoff)




	def plot_stdlit(self):
		data=_finite_plot_data(self.old_lit_radius,self.old_lit_fdnu)
		fig,ax=plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.plot(data['radius'],data['value'],color=PLOT_COLORS['model'],linewidth=1.3,alpha=0.75,label='model sequence')
		ax.scatter(data['radius'],data['value'],color=PLOT_COLORS['model'],alpha=0.68,s=20,edgecolors='white',linewidths=0.35)
		_style_radius_plot(ax,title='Original nonadiabatic model correction')
		_add_stats_box(ax,[
			f'N = {len(data)}',
			f'R = {data["radius"].min():.2g}-{data["radius"].max():.2g} R_sun',
			f'f range = {data["value"].min():.3f}-{data["value"].max():.3f}'
		])
		ax.legend(loc='lower left',framealpha=0.92,fontsize=8)
		plt.tight_layout()
		plt.savefig(self.path_oldlit,dpi=300,bbox_inches='tight',facecolor='white')
		plt.close()
		return self.path_oldlit

	def plot_pySYD(self):
		data=_finite_plot_data(self.radius_pySYD,self.fdnu_pySYD)
		fig,ax=plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.plot(data['radius'],data['value'],color=PLOT_COLORS['pysyd'],linewidth=1.1,alpha=0.45,label='radius-ordered sequence')
		ax.scatter(data['radius'],data['value'],color=PLOT_COLORS['pysyd'],alpha=0.78,s=24,edgecolors='white',linewidths=0.35,label='pySYD')
		_style_radius_plot(ax,title='pySYD correction for original sequence')
		_add_stats_box(ax,[
			f'N = {len(data)}',
			f'median f = {data["value"].median():.3f}',
			f'max f = {data["value"].max():.3f}'
		])
		ax.legend(loc='upper left',framealpha=0.92,fontsize=8)
		plt.tight_layout()
		plt.savefig(self.path_pySYD_plot,dpi=300,bbox_inches='tight',facecolor='white')
		plt.close()
		return self.path_pySYD_plot
		
	def plot_BAM(self):
		data=self.data['BAM'].replace([np.inf,-np.inf],np.nan).dropna(
			subset=['radius','dnu_inf','dnu_BAM','dnu_unc_BAM']
		).copy()
		data=data[(data['radius']>0) & (data['dnu_BAM']>0)].copy()
		data['fdnu_BAM_calc']=data['dnu_inf']/data['dnu_BAM']
		data['fdnu_BAM_unc']=np.abs(data['fdnu_BAM_calc']*data['dnu_unc_BAM']/data['dnu_BAM'])
		before=len(data)
		data=data[(data['radius']<100) & (data['fdnu_BAM_calc'].between(0.5,1.5))].sort_values('radius')
		fig, ax =plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.errorbar(
			data['radius'],
			data['fdnu_BAM_calc'],
			yerr=data['fdnu_BAM_unc'],
			fmt='o',
			color=PLOT_COLORS['bam'],
			ecolor='#7b3294',
			alpha=0.62,
			capsize=2.2,
			ms=3.5,
			elinewidth=0.7,
			capthick=0.7,
			label='BAM'
		)
		_plot_binned_radius_summary(ax,data['radius'],data['fdnu_BAM_calc'],PLOT_COLORS['model'],'median trend',bins=20,min_count=2,alpha=0.10)
		_style_radius_plot(ax,title='BAM correction for original sequence')
		_add_stats_box(ax,[
			f'N = {len(data)}',
			f'filtered = {before-len(data)}',
			f'median f = {data["fdnu_BAM_calc"].median():.3f}'
		])
		ax.legend(loc='best',framealpha=0.92,fontsize=8)
		plt.tight_layout()
		plt.savefig(self.path_BAM_plot,dpi=300,bbox_inches='tight',facecolor='white')
		plt.close()
		return self.path_BAM_plot



	def plot_asfgrid(self):
		data=_finite_plot_data(self.radius_asfgrid,self.fdnu_asfgrid)
		data['inverted']=1.0/data['value']
		fig, ax =plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.plot(data['radius'],data['value'],color=PLOT_COLORS['asfgrid'],linewidth=1.6,alpha=0.75,label='ASFGrid raw')
		ax.scatter(data['radius'],data['value'],color=PLOT_COLORS['asfgrid'],alpha=0.64,s=18,edgecolors='white',linewidths=0.3)
		ax.plot(data['radius'],data['inverted'],color='#006d2c',linewidth=1.1,alpha=0.65,linestyle='--',label='inverted convention')
		_style_radius_plot(ax,title='ASFGrid correction for original sequence')
		_add_stats_box(ax,[
			f'N = {len(data)}',
			f'raw median = {data["value"].median():.4f}',
			f'inv. median = {data["inverted"].median():.4f}'
		],loc='lower right')
		ax.legend(loc='best',framealpha=0.92,fontsize=8)
		plt.tight_layout()
		plt.savefig(self.path_asfgrid_plot,dpi=300,bbox_inches='tight',facecolor='white')

		plt.close()
		return self.path_asfgrid_plot

	def plot_pysyd_asf_interp(self):
		pysyd_data=pd.DataFrame({
			'radius':self.radius_pySYD,
			'fdnu':self.fdnu_pySYD
		}).replace([np.inf,-np.inf],np.nan).dropna().sort_values('radius')
		apokasc_data=self.apokasc_radius_error_data()

		interp_func=interp1d(
			pysyd_data['radius'],
			pysyd_data['fdnu'],
			kind='linear',
			bounds_error=False,
			fill_value=(pysyd_data['fdnu'].iloc[0],pysyd_data['fdnu'].iloc[-1])
		)
		apokasc_data['fdnu_interp']=interp_func(apokasc_data['radius'])
		apokasc_data=self.add_interpolated_fdnu_uncertainty(apokasc_data,interp_func)
		high_radius_unc=apokasc_data['frac_radius_err']>=0.15
		fdnu_errorbars=self.interpolated_fdnu_errorbar_mask(apokasc_data,high_radius_unc)

		fig, ax = plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.plot(pysyd_data['radius'],pysyd_data['fdnu'],color=PLOT_COLORS['model'],alpha=0.48,linewidth=1.2,label='pySYD model sequence')
		ax.scatter(
			apokasc_data['radius'],
			apokasc_data['fdnu_interp'],
			color=PLOT_COLORS['pysyd'],
			alpha=0.36,
			s=5,
			edgecolors='none',
			label='APOKASC radii mapped onto pySYD'
		)
		if fdnu_errorbars.any():
			ax.errorbar(
				apokasc_data.loc[fdnu_errorbars,'radius'],
				apokasc_data.loc[fdnu_errorbars,'fdnu_interp'],
				yerr=apokasc_data.loc[fdnu_errorbars,'fdnu_err_centered'].to_numpy(),
				fmt='none',
				ecolor=PLOT_COLORS['warning'],
				alpha=0.45,
				capsize=1.5,
				elinewidth=0.6,
				capthick=0.6,
				label=r'propagated radius uncertainty'
			)
		_style_radius_plot(ax,ylabel=r'Interpolated $f_{\Delta\nu}$',title='APOKASC radii interpolated onto pySYD sequence')
		_add_stats_box(ax,[
			f'APOKASC N = {len(apokasc_data)}',
			f'error bars shown = {int(fdnu_errorbars.sum())}',
			f'high R-err rows = {int(high_radius_unc.sum())}'
		],loc='lower left')
		ax.legend(loc='best', framealpha=0.9, fontsize=8)

		plt.tight_layout()
		plt.savefig(self.path_pySYDinterp_plot, bbox_inches='tight', dpi=300)
		plt.close()
		return self.path_pySYDinterp_plot



	def plot_asf_self_interp(self):
		asf_data=pd.DataFrame({
			'radius':self.radius_asfgrid,
			'fdnu':self.fdnu_asfgrid
		}).replace([np.inf,-np.inf],np.nan).dropna().sort_values('radius')
		apokasc_data=self.apokasc_radius_error_data()

		interp_func=interp1d(
			asf_data['radius'],
			asf_data['fdnu'],
			kind='linear',
			bounds_error=False,
			fill_value=(asf_data['fdnu'].iloc[0],asf_data['fdnu'].iloc[-1])
		)
		apokasc_data['fdnu_interp']=interp_func(apokasc_data['radius'])
		apokasc_data=self.add_interpolated_fdnu_uncertainty(apokasc_data,interp_func)
		high_radius_unc=apokasc_data['frac_radius_err']>=0.15
		fdnu_errorbars=self.interpolated_fdnu_errorbar_mask(apokasc_data,high_radius_unc)

		fig, ax = plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.plot(asf_data['radius'],asf_data['fdnu'],color=PLOT_COLORS['model'],alpha=0.5,linewidth=1.2,label='ASFGrid sequence')
		ax.scatter(
			apokasc_data['radius'],
			apokasc_data['fdnu_interp'],
			color=PLOT_COLORS['asfgrid'],
			alpha=0.36,
			s=5,
			edgecolors='none',
			label='APOKASC radii mapped onto ASFGrid'
		)
		if fdnu_errorbars.any():
			ax.errorbar(
				apokasc_data.loc[fdnu_errorbars,'radius'],
				apokasc_data.loc[fdnu_errorbars,'fdnu_interp'],
				yerr=apokasc_data.loc[fdnu_errorbars,'fdnu_err_centered'].to_numpy(),
				fmt='none',
				ecolor=PLOT_COLORS['warning'],
				alpha=0.45,
				capsize=1.5,
				elinewidth=0.6,
				capthick=0.6,
				label=r'propagated radius uncertainty'
			)

		_style_radius_plot(ax,ylabel=r'Interpolated $f_{\Delta\nu}$',title='APOKASC radii interpolated onto ASFGrid sequence')
		_add_stats_box(ax,[
			f'APOKASC N = {len(apokasc_data)}',
			f'error bars shown = {int(fdnu_errorbars.sum())}',
			f'high R-err rows = {int(high_radius_unc.sum())}'
		],loc='lower left')
		ax.legend(loc='best', framealpha=0.9, fontsize=8)

		plt.tight_layout()
		plt.savefig(self.path_asfgridinterp_plot, bbox_inches='tight', dpi=300)
		plt.close()
		return self.path_asfgridinterp_plot



	def plot_pySYD_BAM_fractional_uncertainty(self):
		merged_pySYD_BAM=pd.merge(self.data['pySYD'],self.data['BAM'],on=['ID','radius'],how='inner',suffixes=('_pySYD','BAM'))
		filtered=merged_pySYD_BAM.replace([np.inf,-np.inf],np.nan).dropna(
			subset=['radius','dnu_syd','dnu_BAM','dnu_unc_BAM','dnu_infBAM']
		).copy()
		filtered['fdnu_BAM_calc']=filtered['dnu_infBAM']/filtered['dnu_BAM']
		valid_mask=(
			(filtered['radius']<100) &
			(filtered['dnu_syd']>0) &
			(filtered['dnu_BAM']>0) &
			(filtered['dnu_unc_BAM']>=0) &
			(filtered['fdnu_BAM_calc'].between(0.5,1.5))
		)
		filtered=filtered[valid_mask].copy()
		frac_diff=(filtered['dnu_syd']-filtered['dnu_BAM'])/filtered['dnu_syd']
		frac_diff_unc=(filtered['dnu_unc_BAM']/filtered['dnu_syd']).abs()
		removed=len(merged_pySYD_BAM)-len(filtered)
		fig, ax = plt.subplots(figsize=PLOT_FIGSIZE, dpi=PLOT_DPI)
		ax.errorbar(filtered['radius'], frac_diff, yerr=frac_diff_unc,fmt='o', ecolor='#4d4d4d', color=PLOT_COLORS['bam'], alpha=0.64,capsize=2.4, ms=3.7, elinewidth=0.75, capthick=0.75,label='pySYD-BAM')
		ax.axhline(y=0, color=PLOT_COLORS['warning'], linestyle='--', alpha=0.75, linewidth=1.1, label='zero difference')
		_style_radius_plot(
			ax,
			ylabel=r'$(\Delta\nu_{\rm syd}-\Delta\nu_{\rm BAM})/\Delta\nu_{\rm syd}$',
			title='Fractional spacing difference: pySYD versus BAM'
		)
		ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
		stats_text = f'N = {len(filtered)}\nRemoved = {removed}\nMean = {frac_diff.mean():.4f}\nStd = {frac_diff.std():.4f}'
		_add_stats_box(ax,stats_text.split('\n'),loc='lower left')
		ax.legend(loc='best', framealpha=0.9, fontsize=8)
		plt.tight_layout()
		plt.savefig(self.path_pySYDdnuvsBAMdnu_plot,bbox_inches='tight',dpi=300)
		plt.close()
		return self.path_pySYDdnuvsBAMdnu_plot






	def plot_lifdnu(self):
		data=_finite_plot_data(self.radius_li,self.fdnu_li)
		fig,ax=plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.scatter(data['radius'],data['value'],color=PLOT_COLORS['li'],alpha=0.12,s=4,edgecolors='none',label='Li rows')
		_plot_binned_radius_summary(ax,data['radius'],data['value'],PLOT_COLORS['li'],'binned median',bins=45,min_count=8,alpha=0.16)
		_style_radius_plot(ax,title=r'Li reference $f_{\Delta\nu}$, all available modes')
		_add_stats_box(ax,[
			f'N = {len(data)}',
			f'median f = {data["value"].median():.3f}',
			f'R max = {data["radius"].max():.1f} R_sun'
		])
		ax.legend(loc='best', framealpha=0.9,fontsize=8)
		plt.tight_layout()
		plt.savefig(self.path_li_fdnu,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return self.path_li_fdnu








	def plot_lifdnu_l0(self):
		data=_finite_plot_data(self.radius_li_l0,self.fdnu_li_l0)
		fig,ax=plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.scatter(data['radius'],data['value'],color=PLOT_COLORS['asfgrid'],alpha=0.11,s=4,edgecolors='none',label=r'Li $\ell=0$ rows')
		_plot_binned_radius_summary(ax,data['radius'],data['value'],PLOT_COLORS['asfgrid'],'binned median',bins=45,min_count=8,alpha=0.16)
		_style_radius_plot(ax,title=r'Li reference $f_{\Delta\nu}$ from $\ell=0$ modes')
		_add_stats_box(ax,[
			f'N = {len(data)}',
			f'median f = {data["value"].median():.3f}',
			f'R max = {data["radius"].max():.1f} R_sun'
		])
		ax.legend(loc='best', framealpha=0.9,fontsize=8)
		plt.tight_layout()
		plt.savefig(self.path_li_fdnu_l0,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return self.path_li_fdnu_l0



	def plot_lifdnu_fractional(self):
		l0=self.data['li_data_l0_nonad'][['PROFID','radius','fdnu_li_massconst']].rename(columns={'fdnu_li_massconst':'fdnu_l0'}).copy()
		all_modes=self.data['li_data_nonad'][['PROFID','fdnu_li_massconst']].rename(columns={'fdnu_li_massconst':'fdnu_all'}).copy()
		all_modes=all_modes.groupby('PROFID',as_index=False).agg(fdnu_all=('fdnu_all','median'))
		matched=l0.merge(all_modes,on='PROFID',how='inner').replace([np.inf,-np.inf],np.nan).dropna()
		matched=matched[(matched['radius']>0) & (matched['fdnu_l0']>0) & (matched['fdnu_all']>0)].copy()
		matched['ratio']=matched['fdnu_l0']/matched['fdnu_all']
		fig,ax=plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.axhline(1.0,color=PLOT_COLORS['model'],linestyle=':',linewidth=1.1,label='equal correction')
		ax.scatter(matched['radius'],matched['ratio'],color=PLOT_COLORS['bam'],alpha=0.15,s=5,edgecolors='none',label='matched Li profiles')
		_plot_binned_radius_summary(ax,matched['radius'],matched['ratio'],PLOT_COLORS['bam'],'binned median',bins=42,min_count=5,alpha=0.14)
		_style_radius_plot(
			ax,
			ylabel=r'$f_{\Delta\nu}(\ell=0) / f_{\Delta\nu}(\ell=0,1,2,3)$',
			title=r'Li radial-mode correction relative to all-mode correction'
		)
		_add_stats_box(ax,[
			f'N = {len(matched)}',
			f'median ratio = {matched["ratio"].median():.3f}',
			f'16-84% = {matched["ratio"].quantile(0.16):.3f}-{matched["ratio"].quantile(0.84):.3f}'
		],loc='lower left')
		ax.legend(loc='best', framealpha=0.9,fontsize=8)
		plt.tight_layout()
		plt.savefig(self.path_li_fdnu_fractional_comp,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return self.path_li_fdnu_fractional_comp

	def plot_lifdnu_vs_pySYDfdnu(self):
		li_data=_finite_plot_data(self.radius_li,self.fdnu_li)
		pysyd_data=_finite_plot_data(self.radius_li,self.fdnu_li_pySYD)
		fig,ax=plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.scatter(li_data['radius'],li_data['value'],color=PLOT_COLORS['li'],alpha=0.12,s=4,edgecolors='none',label='Li reference rows')
		ax.scatter(pysyd_data['radius'],pysyd_data['value'],color=PLOT_COLORS['pysyd'],alpha=0.14,s=4,edgecolors='none',label='pySYD rows')
		_plot_binned_radius_summary(ax,li_data['radius'],li_data['value'],PLOT_COLORS['li'],'Li binned median',bins=45,min_count=8,alpha=0.12)
		_plot_binned_radius_summary(ax,pysyd_data['radius'],pysyd_data['value'],PLOT_COLORS['pysyd'],'pySYD binned median',bins=45,min_count=8,alpha=0.12)
		_style_radius_plot(ax,title='Li reference versus pySYD on Li spectra')
		_add_stats_box(ax,[
			f'N = {len(li_data)}',
			f'median Li = {li_data["value"].median():.3f}',
			f'median pySYD = {pysyd_data["value"].median():.3f}'
		])
		ax.legend(loc='best', framealpha=0.9,fontsize=8)
		plt.tight_layout()
		plt.savefig(self.path_li_fdnu_vs_pySYD_plot,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return self.path_li_fdnu_vs_pySYD_plot

	def plot_lifdnu_vs_pySYDfdnu_l0(self):
		li_data=_finite_plot_data(self.radius_li_l0,self.fdnu_li_l0)
		pysyd_data=_finite_plot_data(self.radius_li_l0,self.fdnu_li_pySYD_l0)
		fig,ax=plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.scatter(li_data['radius'],li_data['value'],color=PLOT_COLORS['li'],alpha=0.12,s=4,edgecolors='none',label=r'Li $\ell=0$ rows')
		ax.scatter(pysyd_data['radius'],pysyd_data['value'],color=PLOT_COLORS['pysyd'],alpha=0.14,s=4,edgecolors='none',label=r'pySYD $\ell=0$ rows')
		_plot_binned_radius_summary(ax,li_data['radius'],li_data['value'],PLOT_COLORS['li'],'Li binned median',bins=45,min_count=8,alpha=0.12)
		_plot_binned_radius_summary(ax,pysyd_data['radius'],pysyd_data['value'],PLOT_COLORS['pysyd'],'pySYD binned median',bins=45,min_count=8,alpha=0.12)
		_style_radius_plot(ax,title=r'Li $\ell=0$ reference versus pySYD')
		_add_stats_box(ax,[
			f'N = {len(li_data)}',
			f'median Li = {li_data["value"].median():.3f}',
			f'median pySYD = {pysyd_data["value"].median():.3f}'
		])
		ax.legend(loc='best', framealpha=0.9,fontsize=8)
		plt.tight_layout()
		plt.savefig(self.path_li_fdnu_vs_pySYD_l0_plot,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return self.path_li_fdnu_vs_pySYD_l0_plot


	def li_mode_table(self, row):
		l_values=np.asarray(row['l'])
		frequency_values=np.asarray(row['fc'],dtype=float)
		inertia_values=np.asarray(row['E_norm'],dtype=float)
		mode_count=min(len(l_values),len(frequency_values),len(inertia_values))
		mode_data=pd.DataFrame({
			'ell':l_values[:mode_count].astype(int),
			'frequency':frequency_values[:mode_count],
			'inertia':inertia_values[:mode_count]
		}).replace([np.inf,-np.inf],np.nan).dropna()
		mode_data=mode_data[(mode_data['frequency']>0) & (mode_data['inertia']>0)].copy()
		return mode_data

	def li_pysyd_candidates(self):
		candidate_tables=[]
		for data_key,source in (
			('li_data_l0_nonad','Li l=0 pySYD'),
			('li_data_nonad','Li pySYD')
		):
			data=self.data.get(data_key,pd.DataFrame())
			if data.empty or not {'PROFID','dnu_syd'}.issubset(data.columns):
				continue
			candidate=data.copy()
			candidate['pysyd_source']=source
			candidate['PROFID']=pd.to_numeric(candidate['PROFID'],errors='coerce')
			for column in ['dnu_syd','dnu_inf','nu_max','mass','radius','luminosity','mixing_length','metalicity','age']:
				if column in candidate.columns:
					candidate[column]=pd.to_numeric(candidate[column],errors='coerce')
			candidate=candidate.replace([np.inf,-np.inf],np.nan).dropna(subset=['PROFID','dnu_syd'])
			candidate=candidate[candidate['dnu_syd']>0].copy()
			candidate_tables.append(candidate)
		if not candidate_tables:
			return pd.DataFrame()
		return pd.concat(candidate_tables,ignore_index=True,sort=False)

	def match_li_pysyd_row(self, row):
		candidates=self.li_pysyd_candidates()
		if candidates.empty:
			return None

		model_number=int(row['model_number'])
		exact=candidates[candidates['PROFID'].astype(int)==model_number].copy()
		if exact.empty:
			return None

		raw_radius=float(row['radius'])
		raw_numax=float(row['numax'])
		raw_mass=float(row['massinit']) if 'massinit' in row.index and pd.notna(row['massinit']) else float(row['star_mass'])
		exact['match_distance']=0.0
		if 'radius' in exact.columns:
			exact['match_distance']+=((exact['radius']-raw_radius).abs()/max(abs(raw_radius),1.0)).fillna(0.0)
		if 'mass' in exact.columns:
			exact['match_distance']+=((exact['mass']-raw_mass).abs()/max(abs(raw_mass),1.0)).fillna(0.0)
		if 'nu_max' in exact.columns:
			exact['match_distance']+=((exact['nu_max']-raw_numax).abs()/max(abs(raw_numax),1.0)).fillna(0.0)
		return exact.sort_values('match_distance').iloc[0]

	def overlay_dnu_syd_comb(self, ax, plot_data, numax, dnu_syd, x_limits):
		if not np.isfinite(dnu_syd) or dnu_syd <= 0:
			return None
		l0_data=plot_data[plot_data['ell']==0].copy()
		if l0_data.empty:
			return None

		anchor_index=(l0_data['frequency']-numax).abs().idxmin()
		anchor_row=l0_data.loc[anchor_index]
		anchor_frequency=float(anchor_row['frequency'])
		anchor_inertia=float(anchor_row['inertia'])
		if x_limits is None:
			frequency_min=float(plot_data['frequency'].min())
			frequency_max=float(plot_data['frequency'].max())
		else:
			frequency_min,frequency_max=x_limits

		n_min=int(np.floor((frequency_min-anchor_frequency)/dnu_syd))
		n_max=int(np.ceil((frequency_max-anchor_frequency)/dnu_syd))
		comb_frequencies=[
			anchor_frequency+(n*dnu_syd)
			for n in range(n_min,n_max+1)
			if frequency_min <= anchor_frequency+(n*dnu_syd) <= frequency_max
		]
		for index,frequency in enumerate(comb_frequencies):
			ax.axvline(
				frequency,
				color='crimson',
				linestyle=':',
				linewidth=1.0,
				alpha=0.8,
				label=r'pySYD $\Delta\nu_{\rm syd}$ comb' if index == 0 else None
			)
		ax.scatter(
			[anchor_frequency],
			[anchor_inertia],
			facecolors='none',
			edgecolors='crimson',
			linewidths=2.1,
			s=125,
			zorder=6,
			label=r'comb anchor near $\nu_{\max}$'
		)

		next_frequency=anchor_frequency+dnu_syd
		if frequency_min <= next_frequency <= frequency_max:
			log_inertia=np.log10(plot_data['inertia'])
			y_position=10**(log_inertia.min()+0.12*(log_inertia.max()-log_inertia.min()))
			ax.annotate(
				'',
				xy=(next_frequency,y_position),
				xytext=(anchor_frequency,y_position),
				arrowprops=dict(arrowstyle='<->',color='crimson',linewidth=1.0)
			)
			ax.text(
				(anchor_frequency+next_frequency)/2.0,
				y_position*1.08,
				r'$\Delta\nu_{\rm syd}$',
				color='crimson',
				fontsize=8,
				horizontalalignment='center',
				verticalalignment='bottom'
			)
		return anchor_frequency

	def plot_inertia_li(self, targets=None, numax_window=6.0, mode_degrees=(0,1,2,3)):
		if self.data['li_raw'].empty:
			raise ValueError('Li raw parquet is not loaded. Create all_plots(load_li_raw=True) before plotting inertia.')

		required_columns=['model_number','star_mass','radius','numax','Dnu_freq','l','fc','E_norm']
		missing=[column for column in required_columns if column not in self.data['li_raw'].columns]
		if missing:
			raise ValueError(f'Li raw data is missing required inertia columns: {missing}')

		if targets is None:
			targets=[
				{'mass':1.9630752900409345,'radius':4.990136960972681},
				{'mass':1.963073493746886,'radius':5.368039943664774},
				{'mass':0.9950781856529303,'radius':153.6250562372917},
				{'mass':1.000080291584161,'radius':148.4494420049967},
				{'mass':1.0025484631740134,'radius':145.8381050829521}
			]

		scalar_columns=['model_number','star_mass','radius','numax','Dnu_freq']
		li_raw=self.data['li_raw'].copy()
		li_raw[scalar_columns]=li_raw[scalar_columns].replace([np.inf,-np.inf],np.nan)
		li_raw=li_raw.dropna(subset=scalar_columns).copy()
		li_raw=li_raw[
			(li_raw['star_mass']>0) &
			(li_raw['radius']>0) &
			(li_raw['numax']>0) &
			(li_raw['Dnu_freq']>0) &
			li_raw['l'].apply(lambda value: isinstance(value,(list,np.ndarray)) and len(value)>0) &
			li_raw['fc'].apply(lambda value: isinstance(value,(list,np.ndarray)) and len(value)>0) &
			li_raw['E_norm'].apply(lambda value: isinstance(value,(list,np.ndarray)) and len(value)>0)
		].copy()

		selected_rows=[]
		seen_models=set()
		for target in targets:
			mass_target=float(target['mass'])
			radius_target=float(target['radius'])
			distance=(
				((li_raw['star_mass']-mass_target).abs()/max(abs(mass_target),1.0))+
				((li_raw['radius']-radius_target).abs()/max(abs(radius_target),1.0))
			)
			row=li_raw.loc[distance.idxmin()]
			model_key=int(row['model_number'])
			if model_key in seen_models:
				continue
			seen_models.add(model_key)
			selected_rows.append(row)

		output_paths=[]
		degree_styles={
			0:('tab:blue','o'),
			1:('tab:orange','s'),
			2:('tab:green','^'),
			3:('tab:red','D')
		}
		self.plot_dir.mkdir(parents=True,exist_ok=True)

		for row in selected_rows:
			mode_data=self.li_mode_table(row)
			mode_data=mode_data[mode_data['ell'].isin(mode_degrees)].copy()
			if mode_data.empty:
				continue

			numax=float(row['numax'])
			dnu=float(row['Dnu_freq'])
			pysyd_row=self.match_li_pysyd_row(row)
			dnu_syd=np.nan if pysyd_row is None else float(pysyd_row['dnu_syd'])
			frequency_min=max(0.0,numax-(numax_window*dnu))
			frequency_max=numax+(numax_window*dnu)
			windowed=mode_data[
				(mode_data['frequency']>=frequency_min) &
				(mode_data['frequency']<=frequency_max)
			].copy()
			if len(windowed) >= 4:
				plot_data=windowed
				x_limits=(frequency_min,frequency_max)
			else:
				plot_data=mode_data
				x_limits=None

			fig,ax=plt.subplots(figsize=(10,6),dpi=300)
			for ell_value in sorted(plot_data['ell'].unique()):
				color,marker=degree_styles.get(int(ell_value),('black','o'))
				ell_data=plot_data[plot_data['ell']==ell_value]
				ax.scatter(
					ell_data['frequency'],
					ell_data['inertia'],
					color=color,
					marker=marker,
					alpha=0.75,
					s=18,
					edgecolors='none',
					label=rf'$\ell={int(ell_value)}$'
				)

			anchor_frequency=self.overlay_dnu_syd_comb(ax,plot_data,numax,dnu_syd,x_limits)
			ax.axvline(numax,color='black',linestyle='--',linewidth=1.0,alpha=0.75,label=r'$\nu_{\max}$')
			ax.axvspan(
				max(0.0,numax-dnu),
				numax+dnu,
				color='lightgray',
				alpha=0.25,
				label=r'$\nu_{\max}\pm\Delta\nu_{\rm freq}$'
			)
			if x_limits is not None:
				ax.set_xlim(*x_limits)
			ax.set_yscale('log')
			ax.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
			ax.set_xlabel(r'Frequency ($\mu$Hz)',fontsize=12)
			ax.set_ylabel(r'Normalized mode inertia $E_{\mathrm{norm}}$',fontsize=12)
			if pysyd_row is None:
				pysyd_text='dnu_syd = unavailable'
				source_text='pySYD match = none'
			else:
				pysyd_text=f'dnu_syd = {dnu_syd:.5g} uHz'
				source_text=f'pySYD match = {pysyd_row["pysyd_source"]}'
			anchor_text='comb anchor = unavailable' if anchor_frequency is None else f'comb anchor = {anchor_frequency:.5g} uHz'
			dnu_int=float(row['Dnu_int']) if 'Dnu_int' in row.index and pd.notna(row['Dnu_int']) else np.nan
			mass_initial=float(row['massinit']) if 'massinit' in row.index and pd.notna(row['massinit']) else np.nan
			luminosity=float(row['luminosity']) if 'luminosity' in row.index and pd.notna(row['luminosity']) else np.nan
			amlt=float(row['amlt']) if 'amlt' in row.index and pd.notna(row['amlt']) else np.nan
			feh=float(row['fehinit']) if 'fehinit' in row.index and pd.notna(row['fehinit']) else np.nan
			age_gyr=float(row['star_age'])/1e9 if 'star_age' in row.index and pd.notna(row['star_age']) else np.nan
			stats_text=(
				f'Model {int(row["model_number"])}\n'
				f'M = {row["star_mass"]:.4g} Msun, Minit = {mass_initial:.4g} Msun\n'
				f'R = {row["radius"]:.5g} Rsun, L = {luminosity:.5g} Lsun\n'
				f'age = {age_gyr:.4g} Gyr, [Fe/H] = {feh:.4g}, alpha = {amlt:.4g}\n'
				f'nu_max = {numax:.5g} uHz\n'
				f'Dnu_freq = {dnu:.5g} uHz, Dnu_int = {dnu_int:.5g} uHz\n'
				f'{pysyd_text}\n'
				f'{anchor_text}\n'
				f'{source_text}'
			)
			ax.text(
				0.02,
				0.98,
				stats_text,
				transform=ax.transAxes,
				verticalalignment='top',
				horizontalalignment='left',
				fontsize=8.0,
				bbox=dict(boxstyle='round',facecolor='white',edgecolor='gray',alpha=0.92)
			)
			ax.set_title(
				rf'Model {int(row["model_number"])}: {row["star_mass"]:.2f}$M_\odot$, '
				rf'{row["radius"]:.2f}$R_\odot$'
			)
			ax.legend(loc='best',framealpha=0.9,fontsize=8)
			ax.tick_params(axis='both',which='major',labelsize=10)
			plt.tight_layout()

			output_path=self.plot_dir/f'inertia_plot_model{int(row["model_number"])}_mass{row["star_mass"]:.4f}_radius{row["radius"]:.4f}.png'
			plt.savefig(output_path,bbox_inches='tight',facecolor='white',dpi=300)
			plt.close()
			output_paths.append(output_path)

		print(f'Wrote {len(output_paths)} Li inertia plots.')
		for output_path in output_paths:
			print(output_path)
		return output_paths


	def li_l0_fdnu_table(self):
		data=self.data['li_data_l0_nonad'].copy()
		numeric_columns=[
			'PROFID','dnu_syd','dnu_theory','numax_raw','dnu_inf','nu_max',
			'mass','radius','luminosity','mixing_length','metalicity','age',
			'fdnu_li_massconst'
		]
		for column in numeric_columns:
			if column in data.columns:
				data[column]=pd.to_numeric(data[column],errors='coerce')
		data=data.replace([np.inf,-np.inf],np.nan).dropna(
			subset=['ID','PROFID','radius','mass','metalicity','mixing_length','dnu_syd','dnu_inf','fdnu_li_massconst']
		).copy()
		data=data[
			(data['radius']>0) &
			(data['dnu_syd']>0) &
			(data['dnu_inf']>0) &
			(data['fdnu_li_massconst']>0)
		].copy()
		data['fdnu_li']=data['fdnu_li_massconst']
		data['fdnu_pysyd']=data['dnu_inf']/data['dnu_syd']
		data['track_mass']=data['mass'].round(4)
		data['track_metalicity']=data['metalicity'].round(4)
		data['track_mixing_length']=data['mixing_length'].round(4)
		data['track_id']=(
			data['track_mass'].astype(str)+'|'+
			data['track_metalicity'].astype(str)+'|'+
			data['track_mixing_length'].astype(str)
		)

		asfgrid=self.data.get('li_asfgrid_l0',pd.DataFrame()).copy()
		if not asfgrid.empty and {'ID','asf_fdnu'}.issubset(asfgrid.columns):
			asfgrid_columns=['ID','asf_dnu','asf_numax','asf_fdnu','asf_status']
			asfgrid=asfgrid[[column for column in asfgrid_columns if column in asfgrid.columns]].copy()
			for column in ['asf_dnu','asf_numax','asf_fdnu']:
				if column in asfgrid.columns:
					asfgrid[column]=pd.to_numeric(asfgrid[column],errors='coerce')
			asfgrid=asfgrid.replace([np.inf,-np.inf],np.nan).drop_duplicates(subset=['ID'])
			data=data.merge(asfgrid,on='ID',how='left')
		else:
			data['asf_dnu']=np.nan
			data['asf_numax']=np.nan
			data['asf_fdnu']=np.nan
			data['asf_status']='missing'
		data['fdnu_asfgrid_inverted']=1.0/data['asf_fdnu']
		data.loc[~np.isfinite(data['fdnu_asfgrid_inverted']),'fdnu_asfgrid_inverted']=np.nan

		data=data.sort_values(['track_id','radius']).reset_index(drop=True)
		data['track_order']=data.groupby('track_id').cumcount()
		data['previous_PROFID']=data.groupby('track_id')['PROFID'].shift(1)
		data['previous_radius']=data.groupby('track_id')['radius'].shift(1)
		for column in ['fdnu_pysyd','fdnu_li','fdnu_asfgrid_inverted']:
			data[f'{column}_previous']=data.groupby('track_id')[column].shift(1)
			data[f'{column}_jump_signed']=data[column]-data[f'{column}_previous']
			data[f'{column}_jump']=data[f'{column}_jump_signed'].abs()

		data['pysyd_li_disagreement']=(data['fdnu_pysyd']-data['fdnu_li']).abs()
		data['pysyd_asf_disagreement']=(data['fdnu_pysyd']-data['fdnu_asfgrid_inverted']).abs()
		data['pipeline_disagreement_max']=data[['pysyd_li_disagreement','pysyd_asf_disagreement']].max(axis=1)
		score_columns=['fdnu_pysyd_jump','pysyd_li_disagreement','pysyd_asf_disagreement']
		data['jump_score']=0.0
		for column in score_columns:
			data['jump_score']+=data[column].rank(pct=True).fillna(0.0)
		return data

	def li_l0_spectrum_path(self, row):
		id_text=str(row['ID'] if isinstance(row,pd.Series) else row)
		candidates=[
			self.repo_dir/'li_sfc_runs_l0'/f'{id_text}_PS.txt',
			self.repo_dir/'li_sfc_runs'/f'{id_text}_PS_l0.txt',
			self.repo_dir/'li_sfc_runs'/f'{id_text}_PS.txt'
		]
		for candidate in candidates:
			if candidate.exists():
				return candidate
		raise FileNotFoundError(f'No Li l=0 pySYD spectrum found for {id_text}: {candidates}')

	def read_li_l0_spectrum(self, row):
		spectrum_path=self.li_l0_spectrum_path(row)
		spectrum=pd.read_csv(
			spectrum_path,
			sep=r'[\s,]+',
			header=None,
			engine='python',
			comment='#'
		).iloc[:,:2].copy()
		spectrum.columns=['frequency','power']
		spectrum['frequency']=pd.to_numeric(spectrum['frequency'],errors='coerce')
		spectrum['power']=pd.to_numeric(spectrum['power'],errors='coerce')
		spectrum=spectrum.replace([np.inf,-np.inf],np.nan).dropna()
		spectrum=spectrum[(spectrum['frequency']>=0) & (spectrum['power']>=0)].copy()
		if spectrum.empty:
			raise ValueError(f'No usable Li l=0 spectrum rows for {row["ID"]}')
		power_max=spectrum['power'].max()
		if power_max <= 0:
			raise ValueError(f'Li l=0 spectrum for {row["ID"]} has no positive power values.')
		spectrum['power_norm']=spectrum['power']/power_max
		return spectrum_path,spectrum

	def li_l0_cutoff_diagnostics(self, row, window_factor=3.0, plot_factor=5.0):
		diagnostics={
			'cutoff_spectrum_found':False,
			'cutoff_lower_window':np.nan,
			'cutoff_upper_window':np.nan,
			'cutoff_peak_count':np.nan,
			'cutoff_peaks_inside':np.nan,
			'cutoff_peaks_left_near':np.nan,
			'cutoff_peaks_right_near':np.nan,
			'cutoff_near_edge_peak_count':np.nan,
			'cutoff_left_edge_gap_dnu':np.nan,
			'cutoff_right_edge_gap_dnu':np.nan,
			'cutoff_min_edge_gap_dnu':np.nan,
			'cutoff_clipped_side':'missing',
			'cutoff_spectrum_path':''
		}
		try:
			spectrum_path,spectrum=self.read_li_l0_spectrum(row)
			center=float(row['nu_max'])
			dnu=float(row['dnu_syd'])
			local,peaks,lower_window,upper_window,_,_=self.cutoff_peak_table(
				spectrum,
				center,
				dnu,
				window_factor=window_factor,
				plot_factor=plot_factor
			)
		except (FileNotFoundError,ValueError,KeyError) as exc:
			diagnostics['cutoff_error']=str(exc)
			return diagnostics

		if local.empty:
			diagnostics['cutoff_error']='empty local spectrum'
			return diagnostics
		left_near=pd.DataFrame()
		right_near=pd.DataFrame()
		left_gap=np.nan
		right_gap=np.nan
		if not peaks.empty:
			inside=peaks[peaks['inside_window']]
			left_outside=peaks[peaks['frequency']<lower_window]
			right_outside=peaks[peaks['frequency']>upper_window]
			left_near=left_outside[left_outside['frequency']>=lower_window-dnu]
			right_near=right_outside[right_outside['frequency']<=upper_window+dnu]
			if not left_outside.empty:
				left_gap=(lower_window-left_outside['frequency'].max())/dnu
			if not right_outside.empty:
				right_gap=(right_outside['frequency'].min()-upper_window)/dnu
			peaks_inside=len(inside)
		else:
			peaks_inside=0

		near_count=len(left_near)+len(right_near)
		if len(left_near)>0 and len(right_near)>0:
			clipped_side='both'
		elif len(left_near)>0:
			clipped_side='left'
		elif len(right_near)>0:
			clipped_side='right'
		else:
			clipped_side='none'
		edge_gaps=[gap for gap in [left_gap,right_gap] if np.isfinite(gap)]
		min_edge_gap=min(edge_gaps) if edge_gaps else np.nan
		diagnostics.update({
			'cutoff_spectrum_found':True,
			'cutoff_lower_window':lower_window,
			'cutoff_upper_window':upper_window,
			'cutoff_peak_count':len(peaks),
			'cutoff_peaks_inside':peaks_inside,
			'cutoff_peaks_left_near':len(left_near),
			'cutoff_peaks_right_near':len(right_near),
			'cutoff_near_edge_peak_count':near_count,
			'cutoff_left_edge_gap_dnu':left_gap,
			'cutoff_right_edge_gap_dnu':right_gap,
			'cutoff_min_edge_gap_dnu':min_edge_gap,
			'cutoff_clipped_side':clipped_side,
			'cutoff_spectrum_path':str(spectrum_path)
		})
		return diagnostics

	def match_li_raw_row_from_l0(self, row):
		if self.data['li_raw'].empty:
			return None
		required_columns=['model_number','massinit','radius','numax','Dnu_freq','l','fc','E_norm']
		if any(column not in self.data['li_raw'].columns for column in required_columns):
			return None
		model_number=int(row['PROFID'])
		li_raw=self.data['li_raw'].copy()
		model_rows=li_raw[li_raw['model_number'].astype(int)==model_number].copy()
		if model_rows.empty:
			return None
		scalar_columns=['massinit','radius','numax']
		model_rows[scalar_columns]=model_rows[scalar_columns].replace([np.inf,-np.inf],np.nan)
		model_rows=model_rows.dropna(subset=scalar_columns)
		if model_rows.empty:
			return None
		target_mass=float(row['mass'])
		target_radius=float(row['radius'])
		distance=(
			((model_rows['massinit']-target_mass).abs()/max(abs(target_mass),1.0))+
			((model_rows['radius']-target_radius).abs()/max(abs(target_radius),1.0))
		)
		return model_rows.loc[distance.idxmin()]

	def li_l0_inertia_diagnostics(self, row, window_factor=3.0, context_factor=6.0):
		diagnostics={
			'inertia_found':False,
			'inertia_raw_Dnu_freq':np.nan,
			'inertia_raw_Dnu_int':np.nan,
			'inertia_l0_modes_inside_pysyd_window':np.nan,
			'inertia_l0_context_count':np.nan,
			'inertia_anchor_frequency':np.nan,
			'inertia_local_l0_spacing':np.nan,
			'inertia_dnu_syd_l0_frac_offset':np.nan,
			'inertia_comb_mismatch_median':np.nan,
			'inertia_comb_mismatch_max':np.nan
		}
		raw_row=self.match_li_raw_row_from_l0(row)
		if raw_row is None:
			return diagnostics
		mode_data=self.li_mode_table(raw_row)
		l0_data=mode_data[mode_data['ell']==0].sort_values('frequency').copy()
		if l0_data.empty:
			return diagnostics

		numax=float(row['nu_max'])
		dnu_syd=float(row['dnu_syd'])
		raw_dnu=float(raw_row['Dnu_freq']) if 'Dnu_freq' in raw_row.index and pd.notna(raw_row['Dnu_freq']) else np.nan
		raw_dnu_int=float(raw_row['Dnu_int']) if 'Dnu_int' in raw_row.index and pd.notna(raw_row['Dnu_int']) else np.nan
		pysyd_lower=max(0.0,numax-(window_factor*dnu_syd))
		pysyd_upper=numax+(window_factor*dnu_syd)
		l0_inside=l0_data[l0_data['frequency'].between(pysyd_lower,pysyd_upper)]
		context_dnu=raw_dnu if np.isfinite(raw_dnu) and raw_dnu > 0 else dnu_syd
		context_lower=max(0.0,numax-(context_factor*context_dnu))
		context_upper=numax+(context_factor*context_dnu)
		l0_context=l0_data[l0_data['frequency'].between(context_lower,context_upper)].copy()
		if l0_context.empty:
			l0_context=l0_data.copy()
		frequencies=l0_context['frequency'].to_numpy(dtype=float)
		nearest_position=int(np.argmin(np.abs(frequencies-numax)))
		anchor_frequency=float(frequencies[nearest_position])
		spacings=[]
		if nearest_position > 0:
			spacings.append(frequencies[nearest_position]-frequencies[nearest_position-1])
		if nearest_position < len(frequencies)-1:
			spacings.append(frequencies[nearest_position+1]-frequencies[nearest_position])
		local_spacing=float(np.nanmedian(spacings)) if spacings else np.nan
		if np.isfinite(dnu_syd) and dnu_syd > 0:
			residuals=np.abs(((frequencies-anchor_frequency+(0.5*dnu_syd)) % dnu_syd)-(0.5*dnu_syd))/dnu_syd
			comb_mismatch_median=float(np.nanmedian(residuals))
			comb_mismatch_max=float(np.nanmax(residuals))
		else:
			comb_mismatch_median=np.nan
			comb_mismatch_max=np.nan
		frac_offset=(
			abs(dnu_syd-local_spacing)/local_spacing
			if np.isfinite(local_spacing) and local_spacing > 0 else np.nan
		)
		diagnostics.update({
			'inertia_found':True,
			'inertia_raw_Dnu_freq':raw_dnu,
			'inertia_raw_Dnu_int':raw_dnu_int,
			'inertia_l0_modes_inside_pysyd_window':len(l0_inside),
			'inertia_l0_context_count':len(l0_context),
			'inertia_anchor_frequency':anchor_frequency,
			'inertia_local_l0_spacing':local_spacing,
			'inertia_dnu_syd_l0_frac_offset':frac_offset,
			'inertia_comb_mismatch_median':comb_mismatch_median,
			'inertia_comb_mismatch_max':comb_mismatch_max
		})
		return diagnostics

	def select_unified_fdnu_jump_cases(self, table, count=6):
		candidates=table.dropna(subset=['fdnu_pysyd_jump']).copy()
		candidates=candidates[candidates['track_order']>0].copy()
		if candidates.empty:
			return candidates
		candidates=candidates.sort_values('jump_score',ascending=False)
		return candidates.drop_duplicates(subset=['track_id']).head(count).copy()

	def plot_unified_fdnu_jump_summary(self, table, diagnostics, cases):
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		fig,axes=plt.subplots(2,2,figsize=(13,10),dpi=300)
		ax_fdnu,ax_disagreement,ax_cutoff,ax_inertia=axes.ravel()

		ax_fdnu.scatter(table['radius'],table['fdnu_li'],color='lightgray',s=4,alpha=0.45,edgecolors='none',label='Li theory')
		ax_fdnu.scatter(table['radius'],table['fdnu_pysyd'],color='tab:blue',s=4,alpha=0.25,edgecolors='none',label='pySYD')
		asf_mask=table['fdnu_asfgrid_inverted'].notna()
		if asf_mask.any():
			ax_fdnu.scatter(
				table.loc[asf_mask,'radius'],
				table.loc[asf_mask,'fdnu_asfgrid_inverted'],
				color='tab:green',
				s=3,
				alpha=0.18,
				edgecolors='none',
				label='ASFGrid inverted'
			)
		if not cases.empty:
			ax_fdnu.scatter(
				cases['radius'],
				cases['fdnu_pysyd'],
				facecolors='none',
				edgecolors='crimson',
				linewidths=1.8,
				s=110,
				label='selected jump cases'
			)
		ax_fdnu.set_xscale('log')
		ax_fdnu.set_xlabel(r'Radius ($R_\odot$)')
		ax_fdnu.set_ylabel(r'$f_{\Delta\nu}$')
		ax_fdnu.set_title('Li l=0 pipeline comparison')
		ax_fdnu.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_fdnu.legend(loc='best',fontsize=8,framealpha=0.9)

		ax_disagreement.scatter(
			table['radius'],
			table['pysyd_li_disagreement'],
			color='lightgray',
			s=4,
			alpha=0.35,
			edgecolors='none',
			label='all Li l=0 rows'
		)
		if not diagnostics.empty:
			colors=diagnostics['cutoff_peaks_inside'].fillna(0)
			scatter=ax_disagreement.scatter(
				diagnostics['radius'],
				diagnostics['pysyd_li_disagreement'],
				c=colors,
				cmap='magma',
				s=34,
				alpha=0.85,
				edgecolors='black',
				linewidths=0.25,
				label='diagnosed jump candidates'
			)
			fig.colorbar(scatter,ax=ax_disagreement,label='peaks inside pySYD window')
		ax_disagreement.set_xscale('log')
		ax_disagreement.set_xlabel(r'Radius ($R_\odot$)')
		ax_disagreement.set_ylabel(r'$|f_{\Delta\nu,\rm pySYD}-f_{\Delta\nu,\rm Li}|$')
		ax_disagreement.set_title('Pipeline disagreement')
		ax_disagreement.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)

		if not diagnostics.empty:
			ax_cutoff.scatter(
				diagnostics['cutoff_peaks_inside'],
				diagnostics['fdnu_pysyd_jump'],
				c=diagnostics['pysyd_li_disagreement'],
				cmap='viridis',
				s=42,
				alpha=0.85,
				edgecolors='black',
				linewidths=0.25
			)
		ax_cutoff.set_xlabel('peaks inside pySYD window')
		ax_cutoff.set_ylabel(r'local pySYD $f_{\Delta\nu}$ jump')
		ax_cutoff.set_title('Window content vs jump size')
		ax_cutoff.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)

		if not diagnostics.empty:
			inertia_mask=diagnostics['inertia_found'].astype(bool)
			if inertia_mask.any():
				scatter=ax_inertia.scatter(
					diagnostics.loc[inertia_mask,'inertia_dnu_syd_l0_frac_offset'],
					diagnostics.loc[inertia_mask,'pysyd_li_disagreement'],
					c=diagnostics.loc[inertia_mask,'inertia_comb_mismatch_median'],
					cmap='plasma',
					s=44,
					alpha=0.88,
					edgecolors='black',
					linewidths=0.25
				)
				fig.colorbar(scatter,ax=ax_inertia,label='median comb mismatch')
		ax_inertia.set_xlabel(r'$|\Delta\nu_{\rm syd}-\Delta\nu_{\ell=0}|/\Delta\nu_{\ell=0}$')
		ax_inertia.set_ylabel(r'$|f_{\Delta\nu,\rm pySYD}-f_{\Delta\nu,\rm Li}|$')
		ax_inertia.set_title('Inertia/mode-spacing diagnostic')
		ax_inertia.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)

		plt.tight_layout()
		plt.savefig(self.path_unified_fdnu_jump_summary,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return self.path_unified_fdnu_jump_summary

	def plot_unified_fdnu_jump_case(self, row, table, window_factor=3.0, plot_factor=5.0):
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		fig,axes=plt.subplots(2,2,figsize=(13,10),dpi=300)
		ax_context,ax_spectrum,ax_inertia,ax_bar=axes.ravel()
		track=table[table['track_id']==row['track_id']].sort_values('radius')
		ax_context.plot(track['radius'],track['fdnu_li'],color='black',linewidth=1.1,label='Li theory')
		ax_context.plot(track['radius'],track['fdnu_pysyd'],color='tab:blue',linewidth=1.1,label='pySYD')
		if track['fdnu_asfgrid_inverted'].notna().any():
			ax_context.plot(track['radius'],track['fdnu_asfgrid_inverted'],color='tab:green',linewidth=1.0,label='ASFGrid inverted')
		ax_context.scatter([row['radius']],[row['fdnu_pysyd']],facecolors='none',edgecolors='crimson',s=180,linewidths=2.2,zorder=5)
		ax_context.set_xscale('log')
		ax_context.set_xlabel(r'Radius ($R_\odot$)')
		ax_context.set_ylabel(r'$f_{\Delta\nu}$')
		ax_context.set_title(f'Profile {int(row["PROFID"])} in its mass/metallicity track')
		ax_context.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_context.legend(loc='best',fontsize=8,framealpha=0.9)

		try:
			spectrum_path,spectrum=self.read_li_l0_spectrum(row)
			local,peaks,lower_window,upper_window,lower_plot,upper_plot=self.cutoff_peak_table(
				spectrum,
				float(row['nu_max']),
				float(row['dnu_syd']),
				window_factor=window_factor,
				plot_factor=plot_factor
			)
			ax_spectrum.plot(local['frequency'],local['power_norm'],color='black',linewidth=1.0,label='normalized power')
			ax_spectrum.axvspan(lower_window,upper_window,color='tab:blue',alpha=0.14,label='pySYD window')
			ax_spectrum.axvline(float(row['nu_max']),color='tab:blue',linestyle='-',linewidth=1.2,label=r'$\nu_{\max}$')
			ax_spectrum.axvline(lower_window,color='tab:blue',linestyle='--',linewidth=1.0)
			ax_spectrum.axvline(upper_window,color='tab:blue',linestyle='--',linewidth=1.0)
			if not peaks.empty:
				inside=peaks[peaks['inside_window']]
				outside=peaks[~peaks['inside_window']]
				if not inside.empty:
					ax_spectrum.scatter(inside['frequency'],inside['power_norm'],color='tab:green',s=26,zorder=5,label='inside peak')
				if not outside.empty:
					ax_spectrum.scatter(outside['frequency'],outside['power_norm'],facecolors='none',edgecolors='crimson',s=80,linewidths=1.5,zorder=6,label='outside peak')
			ax_spectrum.set_xlim(lower_plot,upper_plot)
			ax_spectrum.set_ylim(bottom=0)
			ax_spectrum.set_title('pySYD cutoff window')
			ax_spectrum.legend(loc='best',fontsize=8,framealpha=0.9)
		except (FileNotFoundError,ValueError,KeyError) as exc:
			ax_spectrum.text(0.5,0.5,f'No cutoff spectrum\n{exc}',ha='center',va='center',transform=ax_spectrum.transAxes)
		ax_spectrum.set_xlabel(r'Frequency ($\mu$Hz)')
		ax_spectrum.set_ylabel('Normalized power')
		ax_spectrum.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)

		raw_row=self.match_li_raw_row_from_l0(row)
		if raw_row is not None:
			mode_data=self.li_mode_table(raw_row)
			raw_dnu=float(raw_row['Dnu_freq']) if pd.notna(raw_row['Dnu_freq']) else float(row['dnu_syd'])
			context_lower=max(0.0,float(row['nu_max'])-(6.0*raw_dnu))
			context_upper=float(row['nu_max'])+(6.0*raw_dnu)
			plot_data=mode_data[mode_data['frequency'].between(context_lower,context_upper)].copy()
			if plot_data.empty:
				plot_data=mode_data.copy()
				x_limits=None
			else:
				x_limits=(context_lower,context_upper)
			degree_styles={0:('tab:blue','o'),1:('tab:orange','s'),2:('tab:green','^'),3:('tab:red','D')}
			for ell_value in sorted(plot_data['ell'].unique()):
				color,marker=degree_styles.get(int(ell_value),('black','o'))
				ell_data=plot_data[plot_data['ell']==ell_value]
				ax_inertia.scatter(ell_data['frequency'],ell_data['inertia'],color=color,marker=marker,alpha=0.75,s=16,edgecolors='none',label=rf'$\ell={int(ell_value)}$')
			self.overlay_dnu_syd_comb(ax_inertia,plot_data,float(row['nu_max']),float(row['dnu_syd']),x_limits)
			ax_inertia.axvline(float(row['nu_max']),color='black',linestyle='--',linewidth=1.0,alpha=0.75,label=r'$\nu_{\max}$')
			if x_limits is not None:
				ax_inertia.set_xlim(*x_limits)
			ax_inertia.set_yscale('log')
			ax_inertia.legend(loc='best',fontsize=7,framealpha=0.9)
		else:
			ax_inertia.text(0.5,0.5,'No matched raw Li inertia row',ha='center',va='center',transform=ax_inertia.transAxes)
		ax_inertia.set_xlabel(r'Frequency ($\mu$Hz)')
		ax_inertia.set_ylabel(r'Normalized mode inertia $E_{\mathrm{norm}}$')
		ax_inertia.set_title(r'Inertia modes with pySYD $\Delta\nu$ comb')
		ax_inertia.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)

		bar_labels=['Li theory','pySYD']
		bar_values=[float(row['fdnu_li']),float(row['fdnu_pysyd'])]
		bar_colors=['black','tab:blue']
		if pd.notna(row.get('fdnu_asfgrid_inverted',np.nan)):
			bar_labels.append('ASFGrid inv.')
			bar_values.append(float(row['fdnu_asfgrid_inverted']))
			bar_colors.append('tab:green')
		ax_bar.bar(bar_labels,bar_values,color=bar_colors,alpha=0.78)
		ax_bar.axhline(float(row['fdnu_li']),color='black',linestyle=':',linewidth=1.0)
		ax_bar.set_ylabel(r'$f_{\Delta\nu}$')
		ax_bar.set_title('Same model, different measurements')
		stats_text=(
			f'PROFID {int(row["PROFID"])}\n'
			f'R = {row["radius"]:.4g} Rsun, Minit = {row["mass"]:.4g} Msun\n'
			f'nu_max = {row["nu_max"]:.4g} uHz, dnu_syd = {row["dnu_syd"]:.4g} uHz\n'
			f'pySYD jump = {row["fdnu_pysyd_jump"]:.4g}\n'
			f'pySYD-Li diff = {row["pysyd_li_disagreement"]:.4g}'
		)
		ax_bar.text(
			0.03,
			0.97,
			stats_text,
			transform=ax_bar.transAxes,
			verticalalignment='top',
			horizontalalignment='left',
			fontsize=8.2,
			bbox=dict(boxstyle='round',facecolor='white',edgecolor='gray',alpha=0.92)
		)
		ax_bar.grid(True,axis='y',color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)

		plt.tight_layout()
		output_path=self.path_unified_fdnu_jump_case(row['PROFID'])
		plt.savefig(output_path,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return output_path

	def unified_fdnu_jump_diagnostics(self, profile_count=6, diagnostics_count=80):
		table=self.li_l0_fdnu_table()
		candidate_rows=self.select_unified_fdnu_jump_cases(table,count=diagnostics_count)
		records=[]
		for _,row in candidate_rows.iterrows():
			record=row.to_dict()
			record.update(self.li_l0_cutoff_diagnostics(row))
			record.update(self.li_l0_inertia_diagnostics(row))
			records.append(record)
		diagnostics=pd.DataFrame(records)
		cases=diagnostics.head(profile_count).copy() if not diagnostics.empty else diagnostics
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		diagnostics.to_csv(self.path_unified_fdnu_jump_table,index=False)
		summary_path=self.plot_unified_fdnu_jump_summary(table,diagnostics,cases)
		case_paths=[]
		for _,case_row in cases.iterrows():
			case_paths.append(self.plot_unified_fdnu_jump_case(case_row,table))
		print(f'Wrote unified f_Delta nu diagnostics table to {self.path_unified_fdnu_jump_table}')
		print(f'Wrote unified f_Delta nu summary plot to {summary_path}')
		for case_path in case_paths:
			print(case_path)
		return {
			'table':self.path_unified_fdnu_jump_table,
			'summary':summary_path,
			'cases':case_paths
		}


	def li_vs_asfgrid(self):
		base_mask=(
			(self.mass_asf_raw>=0.96) & (self.mass_asf_raw<=1.2) &
			(self.radius_asf_raw>=1) & (self.radius_asf_raw< 500) &
			(self.fdnu_asf_raw<=1.4)
		)

		li_data=_finite_plot_data(self.radius_li_nonad,self.fdnu_li_nonad)
		asf_data=_finite_plot_data(self.radius_asf_raw[base_mask],self.fdnu_asf_raw[base_mask])
		fig,ax=plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.scatter(li_data['radius'],li_data['value'],color=PLOT_COLORS['li'],alpha=0.08,s=4,edgecolors='none',label='Li reference rows')
		ax.scatter(asf_data['radius'],asf_data['value'],color=PLOT_COLORS['asfgrid'],alpha=0.08,s=4,edgecolors='none',label='ASFGrid grid rows')
		_plot_binned_radius_summary(ax,li_data['radius'],li_data['value'],PLOT_COLORS['li'],'Li binned median',bins=48,min_count=10,alpha=0.12)
		_plot_binned_radius_summary(ax,asf_data['radius'],asf_data['value'],PLOT_COLORS['asfgrid'],'ASFGrid binned median',bins=48,min_count=10,alpha=0.12)
		_style_radius_plot(ax,title='Li reference compared with ASFGrid grid')
		_add_stats_box(ax,[
			f'Li N = {len(li_data)}',
			f'ASFGrid N = {len(asf_data)}',
			r'mass window: 0.96-1.20 $M_\odot$'
		])
		ax.legend(loc='best', framealpha=0.9,fontsize=8)

		plt.tight_layout()
		plt.savefig(self.path_li_vs_asfgrid,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return self.path_li_vs_asfgrid

	def li_vs_apokasc(self):

		base_mask=(
			(self.mass_apokasc>=0.96) & 
			(self.mass_apokasc<=1.2) &
			(self.radius_apokasc>=1) & 
			(self.radius_apokasc< 200) &
			(self.fdnu_apokasc<=2.0)
		)

		li_data=_finite_plot_data(self.radius_li_nonad,self.fdnu_li_nonad)
		apokasc_data=_finite_plot_data(self.radius_apokasc[base_mask],self.fdnu_apokasc[base_mask])
		fig,ax=plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		ax.scatter(li_data['radius'],li_data['value'],color=PLOT_COLORS['li'],alpha=0.08,s=4,edgecolors='none',label='Li reference rows')
		ax.scatter(apokasc_data['radius'],apokasc_data['value'],color=PLOT_COLORS['apokasc'],alpha=0.16,s=6,edgecolors='none',label='APOKASC-3 rows')
		_plot_binned_radius_summary(ax,li_data['radius'],li_data['value'],PLOT_COLORS['li'],'Li binned median',bins=45,min_count=10,alpha=0.12)
		_plot_binned_radius_summary(ax,apokasc_data['radius'],apokasc_data['value'],PLOT_COLORS['apokasc'],'APOKASC binned median',bins=32,min_count=4,alpha=0.16)
		_style_radius_plot(ax,title='Li reference compared with APOKASC-3')
		_add_stats_box(ax,[
			f'Li N = {len(li_data)}',
			f'APOKASC N = {len(apokasc_data)}',
			r'mass window: 0.96-1.20 $M_\odot$'
		])
		ax.legend(loc='best', framealpha=0.9,fontsize=8)

		plt.tight_layout()
		plt.savefig(self.path_li_vs_apokasc,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return self.path_li_vs_apokasc


	def plot_lifdnu_amlt_color(self):
		data=self.data['li_data_nonad'].replace([np.inf,-np.inf],np.nan).dropna(
			subset=['radius','fdnu_li_massconst','mixing_length','mass']
		).copy()
		data=data[
			(data['radius']>0) &
			(data['fdnu_li_massconst']>0) &
			(data['mixing_length']>0) &
			(data['mass']>=0.99) &
			(data['mass']<=1.1)
		].copy()

		fig,ax=plt.subplots(figsize=PLOT_FIGSIZE,dpi=PLOT_DPI)
		scatter=ax.scatter(
			data['radius'],
			data['fdnu_li_massconst'],
			c=data['mixing_length'],
			cmap='viridis',
			s=8,
			alpha=0.72,
			edgecolors='none'
		)
		colorbar=fig.colorbar(scatter,ax=ax)
		colorbar.set_label(r'Mixing length $\alpha_{\rm MLT}$')
		_style_radius_plot(ax,title=r'Li reference correction colored by mixing length')
		_add_stats_box(ax,[
			f'N = {len(data)}',
			r'mass window: 0.99-1.10 $M_\odot$',
			f'alpha range = {data["mixing_length"].min():.2f}-{data["mixing_length"].max():.2f}'
		],loc='lower left')

		plt.tight_layout()
		plt.savefig(self.path_li_nonad_color,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return self.path_li_nonad_color

	def plot_legacy_comparison_suite(self):
		plot_methods=[
			self.plot_stdlit,
			self.plot_pySYD,
			self.plot_BAM,
			self.plot_asfgrid,
			self.plot_pysyd_asf_interp,
			self.plot_asf_self_interp,
			self.plot_pySYD_BAM_fractional_uncertainty,
			self.plot_lifdnu,
			self.plot_lifdnu_l0,
			self.plot_lifdnu_fractional,
			self.plot_lifdnu_vs_pySYDfdnu,
			self.plot_lifdnu_vs_pySYDfdnu_l0,
			self.li_vs_asfgrid,
			self.li_vs_apokasc,
			self.plot_lifdnu_amlt_color
		]
		output_paths=[]
		for plot_method in plot_methods:
			output_paths.append(plot_method())
		return output_paths

	def standard_pysyd_fdnu_data(self):
		data=self.data['pySYD'].replace([np.inf,-np.inf],np.nan).dropna(
			subset=['ID','radius','nu_max','dnu_syd','dnu_inf']
		).copy()
		for column in ['ID','radius','nu_max','dnu_syd','dnu_inf']:
			data[column]=pd.to_numeric(data[column],errors='coerce')
		data=data.dropna(subset=['ID','radius','nu_max','dnu_syd','dnu_inf']).copy()
		data=data[(data['radius']>0) & (data['nu_max']>0) & (data['dnu_syd']>0)].copy()
		data['fdnu_pysyd']=data['dnu_inf']/data['dnu_syd']
		return data.sort_values('radius')

	def read_apokasc3_benchmark(self, path=None):
		"""Read the filtered APOKASC-3/Gaia benchmark table.

		The table is downloaded from VizieR by the user and intentionally remains
		an ignored data artifact. APOKASC-3 reports inverse Gaia radii derived
		from Gaia luminosity and spectroscopic temperature; this method converts
		that quantity to radius and propagates its tabulated uncertainty.
		"""
		benchmark_path=Path(path).expanduser() if path is not None else self.repo_dir/'apokasc3_gaia_benchmark.tsv'
		if not benchmark_path.exists():
			raise FileNotFoundError(
				f'Missing {benchmark_path}. Download the filtered APOKASC-3 table '
				'from the VizieR query documented in apokasc3_gaia_benchmark_summary.txt.'
			)
		data=pd.read_csv(benchmark_path,sep='\t',comment='#')
		data.columns=[str(column).strip() for column in data.columns]
		data=data.rename(columns={
			'[Fe/H]':'feh',
			'e_[Fe/H]':'feh_err',
			'InvRGaia':'inv_radius_gaia',
			'e_InvRGaia':'inv_radius_gaia_err',
			'Radius':'radius_seismic',
			'e_Radius':'radius_seismic_err',
			'Mass':'mass',
			'e_Mass':'mass_err',
			'fDNu':'fdnu_observed',
			'e_fDNu':'fdnu_observed_err'
		})
		numeric_columns=[
			'KIC','Numax','e_Numax','DNu','e_DNu','fdnu_observed',
			'fdnu_observed_err','mass','mass_err','radius_seismic',
			'radius_seismic_err','Teff','e_Teff','feh','feh_err',
			'inv_radius_gaia','inv_radius_gaia_err','GaiaDR3'
		]
		for column in numeric_columns:
			if column in data.columns:
				data[column]=pd.to_numeric(data[column],errors='coerce')
		data['EvolSt']=data['EvolSt'].astype(str).str.strip()
		data=data[data['EvolSt'].str.upper().eq('RGB')].copy()
		data=data[
			(data['mass']>0.0) &
			(data['radius_seismic']>0.0) &
			(data['inv_radius_gaia']>0.0) &
			(data['fdnu_observed']>0.0)
		].copy()
		data['radius_gaia']=1.0/data['inv_radius_gaia']
		data['radius_gaia_err']=data['inv_radius_gaia_err']/(data['inv_radius_gaia']**2)
		data['radius_residual']=(data['radius_seismic']-data['radius_gaia'])/data['radius_gaia']
		return data.replace([np.inf,-np.inf],np.nan).dropna(
			subset=['radius_gaia','radius_gaia_err','radius_residual','fdnu_observed']
		).sort_values('radius_gaia').reset_index(drop=True)

	def _apokasc3_model_tracks(self):
		pysyd=self.standard_pysyd_fdnu_data().sort_values('radius')
		literature=pd.read_csv(self.repo_dir/'old_lit_nonad_1M0_0D01.7ALPHA.csv').copy()
		for column in ['radius','fdnustdlit']:
			literature[column]=pd.to_numeric(literature[column],errors='coerce')
		literature=literature.replace([np.inf,-np.inf],np.nan).dropna(subset=['radius','fdnustdlit']).sort_values('radius')
		return pysyd,literature

	def plot_apokasc3_gaia_benchmark(self, path=None):
		"""Compare Gaia-based and seismic radii and place f_Dnu in context."""
		data=self.read_apokasc3_benchmark(path=path)
		pysyd,literature=self._apokasc3_model_tracks()
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		output_path=self.plot_dir/'apokasc3_gaia_benchmark.png'
		fig,axes=plt.subplots(1,3,figsize=(17.0,5.7),dpi=PLOT_DPI)
		mass_values=data['mass'].to_numpy()
		scatter=axes[0].errorbar(
			data['radius_gaia'],data['radius_seismic'],
			xerr=data['radius_gaia_err'],yerr=data['radius_seismic_err'],
			fmt='o',markersize=3.5,alpha=0.38,color=PLOT_COLORS['apokasc'],
			ecolor='#b9a2c4',elinewidth=0.35,capsize=0
		)
		x_limits=(max(2.5,data['radius_gaia'].min()*0.85),data['radius_gaia'].max()*1.08)
		axes[0].plot(x_limits,x_limits,color='#343a40',linestyle='--',linewidth=1.1,label='1:1')
		axes[0].set_xscale('log'); axes[0].set_yscale('log'); axes[0].set_xlim(*x_limits); axes[0].set_ylim(*x_limits)
		axes[0].set_xlabel(r'Gaia-based radius ($R_\odot$)')
		axes[0].set_ylabel(r'APOKASC-3 seismic radius ($R_\odot$)')
		axes[0].set_title('Radius benchmark')
		axes[0].legend(loc='upper left',fontsize=8)

		axes[1].scatter(data['radius_gaia'],data['radius_residual'],c=mass_values,cmap='viridis',s=14,alpha=0.55,edgecolors='none')
		axes[1].axhline(0.0,color='#343a40',linestyle='--',linewidth=1.0)
		axes[1].axvline(50.0,color='#c81d4f',linestyle=':',linewidth=1.1,label=r'$50R_\odot$')
		axes[1].set_xscale('log'); axes[1].set_xlabel(r'Gaia-based radius ($R_\odot$)'); axes[1].set_ylabel(r'$(R_{\rm seis}-R_{\rm Gaia})/R_{\rm Gaia}$')
		axes[1].set_title('Radius residual')
		axes[1].legend(loc='best',fontsize=8)
		axes[1].grid(True,color='#d7dce2',linestyle=':',linewidth=0.6)

		axes[2].errorbar(data['radius_gaia'],data['fdnu_observed'],yerr=data['fdnu_observed_err'],fmt='o',markersize=3.2,alpha=0.38,color=PLOT_COLORS['apokasc'],ecolor='#b9a2c4',elinewidth=0.35,capsize=0,label='APOKASC-3 observed')
		axes[2].plot(pysyd['radius'],pysyd['fdnu_pysyd'],color=PLOT_COLORS['pysyd'],linewidth=1.5,label='nonad pySYD model')
		axes[2].plot(literature['radius'],literature['fdnustdlit'],color=PLOT_COLORS['model'],linewidth=1.4,label='standard literature model')
		axes[2].set_xscale('log'); axes[2].set_xlabel(r'Radius ($R_\odot$)'); axes[2].set_ylabel(r'$f_{\Delta\nu}$')
		axes[2].set_title(r'$f_{\Delta\nu}$ context')
		axes[2].legend(loc='best',fontsize=7.8)
		axes[2].grid(True,color='#d7dce2',linestyle=':',linewidth=0.6)
		_add_stats_box(axes[0],[f'N = {len(data):,}',r'$0.8<M/M_\odot<1.2$',r'$-0.3<[Fe/H]<0.3$', 'RGB only'],loc='lower right')
		fig.suptitle('APOKASC-3/Gaia benchmark for the nonadiabatic sequence',fontsize=15)
		fig.tight_layout(rect=[0,0,1,0.94])
		fig.savefig(output_path,bbox_inches='tight',facecolor='white',dpi=PLOT_DPI)
		plt.close(fig)
		return output_path

	def plot_apokasc3_high_radius_fdnu(self, path=None, radius_cut=50.0):
		"""Show the high-radius observational and model correction behavior."""
		data=self.read_apokasc3_benchmark(path=path)
		pysyd,literature=self._apokasc3_model_tracks()
		data=data[data['radius_gaia']>=radius_cut].copy()
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		output_path=self.plot_dir/'apokasc3_high_radius_fdnu.png'
		fig,ax=plt.subplots(figsize=(9.2,6.2),dpi=PLOT_DPI)
		ax.errorbar(data['radius_gaia'],data['fdnu_observed'],yerr=data['fdnu_observed_err'],fmt='o',markersize=4,alpha=0.48,color=PLOT_COLORS['apokasc'],ecolor='#b9a2c4',elinewidth=0.4,capsize=0,label='APOKASC-3 RGB')
		pysyd_high=pysyd[pysyd['radius']>=radius_cut]
		lit_high=literature[literature['radius']>=radius_cut]
		ax.plot(pysyd_high['radius'],pysyd_high['fdnu_pysyd'],color=PLOT_COLORS['pysyd'],linewidth=1.7,label='nonad pySYD model')
		ax.plot(lit_high['radius'],lit_high['fdnustdlit'],color=PLOT_COLORS['model'],linewidth=1.5,label='standard literature model')
		ax.axhline(1.0,color='#343a40',linestyle='--',linewidth=1.0,label=r'$f_{\Delta\nu}=1$')
		ax.set_xscale('log'); ax.set_xlabel(r'Gaia-based radius ($R_\odot$)'); ax.set_ylabel(r'$f_{\Delta\nu}$')
		ax.set_title(rf'High-radius $f_{{\Delta\nu}}$ benchmark ($R\geq{radius_cut:g}R_\odot$)')
		ax.grid(True,color='#d7dce2',linestyle=':',linewidth=0.7)
		ax.legend(loc='best',fontsize=9)
		_add_stats_box(ax,[f'N = {len(data):,}', 'RGB only', 'mass: 0.8-1.2 $M_\\odot$', 'metallicity: -0.3 to +0.3 dex', 'observed f_Dnu uses APOKASC-3/Mosser convention'],loc='upper right')
		fig.tight_layout()
		fig.savefig(output_path,bbox_inches='tight',facecolor='white',dpi=PLOT_DPI)
		plt.close(fig)
		return output_path

	def plot_comprehensive_fdnu_benchmark(self, path=None, radius_cut=30.0):
		"""Compare model, Li-pipeline, ASFGrid, and Gaia/APOKASC corrections."""
		li_table=self.li_l0_fdnu_table().copy()
		gaia=self.read_apokasc3_benchmark(path=path)
		pysyd,literature=self._apokasc3_model_tracks()
		asf=self.data['asfgrid'].copy()
		for column in ['radius','asf_fdnu']:
			asf[column]=pd.to_numeric(asf[column],errors='coerce')
		asf=asf.replace([np.inf,-np.inf],np.nan).dropna(subset=['radius','asf_fdnu']).sort_values('radius')
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		output_path=self.plot_dir/'comprehensive_fdnu_benchmark.png'
		fig,axes=plt.subplots(2,2,figsize=(14.5,10.5),dpi=PLOT_DPI)

		series=[
			('fdnu_li','Li reference',PLOT_COLORS['li']),
			('fdnu_pysyd','Li through pySYD',PLOT_COLORS['pysyd']),
			('fdnu_asfgrid_inverted','Li through ASFGrid (inverted)',PLOT_COLORS['asfgrid'])
		]
		for column,label,color in series:
			_plot_binned_radius_summary(axes[0,0],li_table['radius'],li_table[column],color,label,bins=48,min_count=20,linewidth=2.0,alpha=0.13)
		axes[0,0].plot(pysyd['radius'],pysyd['fdnu_pysyd'],color='#1769aa',linewidth=1.2,linestyle='--',label='original nonad pySYD')
		axes[0,0].plot(literature['radius'],literature['fdnustdlit'],color='#343a40',linewidth=1.2,linestyle='--',label='standard literature')
		axes[0,0].plot(asf['radius'],asf['asf_fdnu'],color='#238b45',linewidth=1.2,linestyle=':',label='original ASFGrid')
		axes[0,0].scatter(gaia['radius_gaia'],gaia['fdnu_observed'],color=PLOT_COLORS['apokasc'],s=8,alpha=0.22,edgecolors='none',label='APOKASC-3 observed')
		axes[0,0].set_xscale('log'); axes[0,0].set_xlabel(r'Radius ($R_\odot$)'); axes[0,0].set_ylabel(r'$f_{\Delta\nu}$'); axes[0,0].set_title('All correction tracks')
		axes[0,0].legend(loc='best',fontsize=7.4,ncol=2); axes[0,0].grid(True,color='#d7dce2',linestyle=':',linewidth=0.6)

		li_high=li_table[li_table['radius']>=radius_cut]
		for column,label,color in series:
			_plot_binned_radius_summary(axes[0,1],li_high['radius'],li_high[column],color,label,bins=30,min_count=10,linewidth=2.0,alpha=0.13)
		pysyd_high=pysyd[pysyd['radius']>=radius_cut]
		lit_high=literature[literature['radius']>=radius_cut]
		asf_high=asf[asf['radius']>=radius_cut]
		gaia_high=gaia[gaia['radius_gaia']>=radius_cut]
		axes[0,1].plot(pysyd_high['radius'],pysyd_high['fdnu_pysyd'],color='#1769aa',linewidth=1.2,linestyle='--',label='original nonad pySYD')
		axes[0,1].plot(lit_high['radius'],lit_high['fdnustdlit'],color='#343a40',linewidth=1.2,linestyle='--',label='standard literature')
		axes[0,1].plot(asf_high['radius'],asf_high['asf_fdnu'],color='#238b45',linewidth=1.2,linestyle=':',label='original ASFGrid')
		axes[0,1].errorbar(gaia_high['radius_gaia'],gaia_high['fdnu_observed'],yerr=gaia_high['fdnu_observed_err'],fmt='o',color=PLOT_COLORS['apokasc'],ecolor='#b9a2c4',markersize=3.2,alpha=0.55,elinewidth=0.4,capsize=0,label='APOKASC-3 observed')
		axes[0,1].set_xscale('log'); axes[0,1].set_xlabel(r'Radius ($R_\odot$)'); axes[0,1].set_ylabel(r'$f_{\Delta\nu}$'); axes[0,1].set_title(rf'High-radius comparison ($R\geq{radius_cut:g}R_\odot$)')
		axes[0,1].set_ylim(0.8,1.7); axes[0,1].legend(loc='best',fontsize=7.4); axes[0,1].grid(True,color='#d7dce2',linestyle=':',linewidth=0.6)

		li_table['pysyd_minus_li']=li_table['fdnu_pysyd']-li_table['fdnu_li']
		li_table['asf_minus_li']=li_table['fdnu_asfgrid_inverted']-li_table['fdnu_li']
		_plot_binned_radius_summary(axes[1,0],li_table['radius'],li_table['pysyd_minus_li'],PLOT_COLORS['pysyd'],'Li pySYD - Li reference',bins=48,min_count=20,linewidth=2.0,alpha=0.13)
		_plot_binned_radius_summary(axes[1,0],li_table['radius'],li_table['asf_minus_li'],PLOT_COLORS['asfgrid'],'Li ASFGrid - Li reference',bins=48,min_count=20,linewidth=2.0,alpha=0.13)
		axes[1,0].axhline(0.0,color='#343a40',linestyle='--',linewidth=1.0)
		axes[1,0].axvline(radius_cut,color='#c81d4f',linestyle=':',linewidth=1.0,label=rf'$R={radius_cut:g}R_\odot$')
		axes[1,0].set_xscale('log'); axes[1,0].set_xlabel(r'Radius ($R_\odot$)'); axes[1,0].set_ylabel(r'$\Delta f_{\Delta\nu}$ relative to Li reference'); axes[1,0].set_title('Pipeline/reference offsets')
		axes[1,0].legend(loc='best',fontsize=8); axes[1,0].grid(True,color='#d7dce2',linestyle=':',linewidth=0.6)

		axes[1,1].scatter(gaia['radius_gaia'],gaia['fdnu_observed'],color=PLOT_COLORS['apokasc'],s=10,alpha=0.28,edgecolors='none',label='APOKASC-3 observed')
		axes[1,1].plot(pysyd['radius'],pysyd['fdnu_pysyd'],color=PLOT_COLORS['pysyd'],linewidth=1.5,label='original nonad pySYD')
		axes[1,1].plot(literature['radius'],literature['fdnustdlit'],color='#343a40',linewidth=1.4,label='standard literature')
		axes[1,1].plot(asf['radius'],asf['asf_fdnu'],color=PLOT_COLORS['asfgrid'],linewidth=1.3,label='original ASFGrid')
		axes[1,1].axvline(radius_cut,color='#c81d4f',linestyle=':',linewidth=1.0)
		axes[1,1].set_xscale('log'); axes[1,1].set_xlabel(r'Radius ($R_\odot$)'); axes[1,1].set_ylabel(r'$f_{\Delta\nu}$'); axes[1,1].set_title('External Gaia/APOKASC context')
		axes[1,1].legend(loc='best',fontsize=8); axes[1,1].grid(True,color='#d7dce2',linestyle=':',linewidth=0.6)
		_add_stats_box(axes[1,1],[f'Li rows = {len(li_table):,}',f'Gaia/APOKASC rows = {len(gaia):,}','Li curves are binned medians','Gaia f_Dnu uses APOKASC-3 convention'],loc='upper right')
		fig.suptitle('Unified $f_{\\Delta\\nu}$ comparison across models, pipelines, and Gaia',fontsize=16)
		fig.tight_layout(rect=[0,0,1,0.95])
		fig.savefig(output_path,bbox_inches='tight',facecolor='white',dpi=PLOT_DPI)
		plt.close(fig)
		return output_path

	def apokasc3_gaia_required_fdnu(self, path=None, numax_sun=3090.0, dnu_sun=135.1, teff_sun=5772.0):
		"""Calculate f_Dnu required for the raw scaling radius to equal Gaia radius."""
		data=self.read_apokasc3_benchmark(path=path)
		data['radius_raw_scaling']=(
			(data['Numax']/numax_sun)*
			(data['DNu']/dnu_sun)**-2*
			np.sqrt(data['Teff']/teff_sun)
		)
		data['fdnu_gaia_required']=np.sqrt(data['radius_raw_scaling']/data['radius_gaia'])
		return data.replace([np.inf,-np.inf],np.nan).dropna(subset=['radius_raw_scaling','fdnu_gaia_required']).copy()

	def plot_gaia_required_fdnu_comparison(self, path=None, radius_cut=50.0):
		"""Rank correction tracks by the radius residual they produce against Gaia."""
		gaia=self.apokasc3_gaia_required_fdnu(path=path)
		li_table=self.li_l0_fdnu_table().copy()
		pysyd,literature=self._apokasc3_model_tracks()
		asf=self.data['asfgrid'].copy()
		for column in ['radius','asf_fdnu']:
			asf[column]=pd.to_numeric(asf[column],errors='coerce')
		asf=asf.replace([np.inf,-np.inf],np.nan).dropna(subset=['radius','asf_fdnu']).sort_values('radius')
		li_curves={
			'Li reference':('fdnu_li',PLOT_COLORS['li']),
			'Li through pySYD':('fdnu_pysyd',PLOT_COLORS['pysyd']),
			'Li through ASFGrid':('fdnu_asfgrid_inverted',PLOT_COLORS['asfgrid'])
		}
		candidate_sources={}
		for name,(column,color) in li_curves.items():
			binned=_binned_radius_summary(li_table['radius'],li_table[column],bins=60,min_count=20)
			candidate_sources[name]=(binned['radius_median'].to_numpy(),binned['value_median'].to_numpy(),color)
		candidate_sources['original nonad pySYD']=(pysyd['radius'].to_numpy(),pysyd['fdnu_pysyd'].to_numpy(),'#1769aa')
		candidate_sources['standard literature']=(literature['radius'].to_numpy(),literature['fdnustdlit'].to_numpy(),'#343a40')
		candidate_sources['original ASFGrid']=(asf['radius'].to_numpy(),asf['asf_fdnu'].to_numpy(),'#238b45')
		metrics=[]
		for name,(track_radius,track_fdnu,color) in candidate_sources.items():
			order=np.argsort(track_radius)
			track_radius=track_radius[order]; track_fdnu=track_fdnu[order]
			model_fdnu=np.interp(np.log(gaia['radius_gaia']),np.log(track_radius),track_fdnu,left=np.nan,right=np.nan)
			residual=gaia['radius_raw_scaling']/(model_fdnu**2)/gaia['radius_gaia']-1.0
			gaia[name]=model_fdnu
			gaia[f'{name}_radius_residual']=residual
			for cut_label,subset in [('all',np.ones(len(gaia),dtype=bool)),('high',gaia['radius_gaia']>=radius_cut)]:
				values=pd.Series(residual[subset]).replace([np.inf,-np.inf],np.nan).dropna().to_numpy()
				metrics.append({'method':name,'subset':cut_label,'mae':np.mean(np.abs(values)),'rmse':np.sqrt(np.mean(values**2)),'median':np.median(values),'n':len(values)})
		metrics=pd.DataFrame(metrics)
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		output_path=self.plot_dir/'gaia_required_fdnu_comparison.png'
		fig,axes=plt.subplots(2,2,figsize=(14.5,10.5),dpi=PLOT_DPI)
		axes[0,0].scatter(gaia['radius_gaia'],gaia['fdnu_gaia_required'],s=10,alpha=0.24,color=PLOT_COLORS['apokasc'],edgecolors='none',label='Gaia-required correction')
		for name,(track_radius,track_fdnu,color) in candidate_sources.items():
			axes[0,0].plot(track_radius,track_fdnu,color=color,linewidth=1.45,label=name)
		axes[0,0].set_xscale('log'); axes[0,0].set_xlabel(r'Gaia-based radius ($R_\odot$)'); axes[0,0].set_ylabel(r'$f_{\Delta\nu}$'); axes[0,0].set_title('Correction required by Gaia radius')
		axes[0,0].set_ylim(0.8,1.7); axes[0,0].legend(loc='best',fontsize=7.2,ncol=2); axes[0,0].grid(True,color='#d7dce2',linestyle=':',linewidth=0.6)

		for name,(track_radius,track_fdnu,color) in candidate_sources.items():
			_plot_binned_radius_summary(axes[0,1],gaia['radius_gaia'],gaia[f'{name}_radius_residual'],color,name,bins=38,min_count=12,linewidth=1.8,alpha=0.10)
		axes[0,1].axhline(0.0,color='#343a40',linestyle='--',linewidth=1.0)
		axes[0,1].axvline(radius_cut,color='#c81d4f',linestyle=':',linewidth=1.0,label=rf'$R={radius_cut:g}R_\odot$')
		axes[0,1].set_xscale('log'); axes[0,1].set_xlabel(r'Gaia-based radius ($R_\odot$)'); axes[0,1].set_ylabel(r'Predicted radius residual $(R_{\rm model}-R_{\rm Gaia})/R_{\rm Gaia}$'); axes[0,1].set_title('Radius residual from each correction')
		axes[0,1].legend(loc='best',fontsize=7.2,ncol=2); axes[0,1].grid(True,color='#d7dce2',linestyle=':',linewidth=0.6)

		method_order=list(candidate_sources)
		x=np.arange(len(method_order)); width=0.38
		all_mae=[metrics[(metrics.method==name)&(metrics.subset=='all')].iloc[0].mae for name in method_order]
		high_mae=[metrics[(metrics.method==name)&(metrics.subset=='high')].iloc[0].mae for name in method_order]
		axes[1,0].bar(x-width/2,all_mae,width,label='All RGB',color='#9ecae1')
		axes[1,0].bar(x+width/2,high_mae,width,label=rf'$R\geq{radius_cut:g}R_\odot$',color='#2171b5')
		axes[1,0].set_xticks(x,method_order,rotation=35,ha='right',fontsize=8); axes[1,0].set_ylabel('Mean absolute radius residual'); axes[1,0].set_title('Correction ranking by Gaia residual'); axes[1,0].legend(fontsize=8); axes[1,0].grid(True,axis='y',color='#d7dce2',linestyle=':',linewidth=0.6)

		high=gaia['radius_gaia']>=radius_cut
		box_values=[]
		box_labels=[]
		for name in method_order:
			values=gaia.loc[high,f'{name}_radius_residual'].replace([np.inf,-np.inf],np.nan).dropna().to_numpy()
			box_values.append(values); box_labels.append(name)
		box=axes[1,1].boxplot(box_values,patch_artist=True,showfliers=False,medianprops={'color':'#111827','linewidth':1.2})
		for patch,(name,(_,_,color)) in zip(box['boxes'],candidate_sources.items()):
			patch.set_facecolor(color); patch.set_alpha(0.55)
		axes[1,1].axhline(0.0,color='#343a40',linestyle='--',linewidth=1.0)
		axes[1,1].set_xticks(np.arange(1,len(box_labels)+1),box_labels,rotation=35,ha='right',fontsize=8); axes[1,1].set_ylabel('Radius residual'); axes[1,1].set_title(rf'High-radius residual distributions ($R\geq{radius_cut:g}R_\odot$)'); axes[1,1].grid(True,axis='y',color='#d7dce2',linestyle=':',linewidth=0.6)
		_add_stats_box(axes[1,1],[f'N = {int(high.sum())}',r'$R_{\rm raw}= (\nu_{\max}/\nu_{\max,\odot})(\Delta\nu/\Delta\nu_\odot)^{-2}(T_{\rm eff}/T_{\rm eff,\odot})^{1/2}$',r'then $f_{\Delta\nu,\rm Gaia}=\sqrt{R_{\rm raw}/R_{\rm Gaia}}$'],loc='upper right')
		fig.suptitle('Gaia-required $f_{\\Delta\\nu}$ and correction ranking',fontsize=16)
		fig.tight_layout(rect=[0,0,1,0.95])
		fig.savefig(output_path,bbox_inches='tight',facecolor='white',dpi=PLOT_DPI)
		plt.close(fig)
		return output_path,metrics

	def standard_pysyd_spectrum_path(self, profile_id):
		profile_id=int(profile_id)
		base_dir=self.repo_dir/'test_runs'/'nonad_1.0msun_0.0d0feh_mass_loss_1.7alpha'
		candidates=[
			base_dir/f'profile{profile_id}.data.GYRE'/f'{profile_id}_PS.txt',
			base_dir/'list_of_bam_runs'/f'profile{profile_id}freq_amp.txt'
		]
		for candidate in candidates:
			if candidate.exists():
				return candidate
		raise FileNotFoundError(f'No pySYD cutoff spectrum found for profile {profile_id}: {candidates}')

	def read_standard_pysyd_spectrum(self, profile_id):
		spectrum_path=self.standard_pysyd_spectrum_path(profile_id)
		spectrum=pd.read_csv(
			spectrum_path,
			sep=r'[\s,]+',
			header=None,
			engine='python',
			comment='#'
		).iloc[:,:2].copy()
		spectrum.columns=['frequency','power']
		spectrum['frequency']=pd.to_numeric(spectrum['frequency'],errors='coerce')
		spectrum['power']=pd.to_numeric(spectrum['power'],errors='coerce')
		spectrum=spectrum.replace([np.inf,-np.inf],np.nan).dropna()
		spectrum=spectrum[(spectrum['frequency']>=0) & (spectrum['power']>=0)].copy()
		if spectrum.empty:
			raise ValueError(f'No usable frequency/power rows for profile {profile_id}')
		power_max=spectrum['power'].max()
		if power_max <= 0:
			raise ValueError(f'Power spectrum for profile {profile_id} has no positive power values.')
		spectrum['power_norm']=spectrum['power']/power_max
		return spectrum_path,spectrum

	def plot_pysyd_echelle(self, profile_id, frequency_half_width=6.0):
		"""Plot a pySYD power spectrum folded on the measured dnu_syd.

		The source spectrum and synchronized pySYD row are kept separate so the
		figure documents the measured spacing used for the fold rather than a
		theoretical spacing chosen after inspecting the result.
		"""
		profile_id=int(profile_id)
		fdnu_data=self.standard_pysyd_fdnu_data()
		profile_rows=fdnu_data[fdnu_data['ID'].astype(int)==profile_id]
		if profile_rows.empty:
			raise ValueError(f'Profile {profile_id} is not in the pySYD synchronized table.')
		row=profile_rows.iloc[0]
		spectrum_path,spectrum=self.read_standard_pysyd_spectrum(profile_id)
		center=float(row['nu_max'])
		dnu=float(row['dnu_syd'])
		lower=max(0.0,center-(frequency_half_width*dnu))
		upper=center+(frequency_half_width*dnu)
		local=spectrum[
			(spectrum['frequency']>=lower) &
			(spectrum['frequency']<=upper)
		].copy()
		if local.empty:
			raise ValueError(f'Profile {profile_id} has no spectrum data near nu_max={center}.')

		local['folded_frequency']=np.mod(local['frequency']-center,dnu)
		local['log_power']=np.log10(np.maximum(local['power_norm'],1e-12))
		mass=float(row['mass'])/1.98847e33 if np.isfinite(float(row['mass'])) and float(row['mass'])>1e30 else float(row['mass'])
		metallicity=str(row['metalicity']).replace('    ',' ')
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		path=self.plot_dir/f'pysyd_echelle_profile{profile_id}.png'
		fig,ax=plt.subplots(figsize=(8.8,6.5),dpi=PLOT_DPI)
		scatter=ax.scatter(
			local['folded_frequency'],
			local['frequency'],
			c=local['log_power'],
			cmap='magma',
			s=10,
			alpha=0.82,
			edgecolors='none',
			rasterized=True
		)
		ax.axvline(0.0,color='#2c7fb8',linewidth=1.4,linestyle='--',label=r'fold origin: $\nu_{\max}$')
		ax.set_xlim(0.0,dnu)
		ax.set_ylim(lower,upper)
		ax.set_xlabel(r'Frequency modulo $\Delta\nu_{\rm syd}$')
		ax.set_ylabel(r'Frequency ($\mu$Hz)')
		ax.set_title(f'pySYD echelle diagram: profile {profile_id}',pad=10)
		ax.grid(True,color='#d7dce2',linestyle=':',linewidth=0.6,alpha=0.7)
		ax.legend(loc='lower right',fontsize=8,framealpha=0.9)
		colorbar=fig.colorbar(scatter,ax=ax,pad=0.02)
		colorbar.set_label('log10 normalized power')
		_add_stats_box(ax,[
			f'$R$ = {float(row["radius"]):.2f} $R_\\odot$',
			f'$M$ = {mass:.3f} $M_\\odot$',
			f'$\\nu_{{\\max}}$ = {center:.4g} $\\mu$Hz',
			f'$\\Delta\\nu_{{\\rm syd}}$ = {dnu:.4g} $\\mu$Hz',
			f'$f_{{\\Delta\\nu}}^{{\\rm syd}}$ = {float(row["fdnu_pysyd"]):.4g}',
			f'[Fe/H] = {metallicity}'
		],loc='upper right')
		fig.text(0.01,0.01,f'Source: {spectrum_path.name}; fold spacing is the measured pySYD dnu_syd',fontsize=7.5,color='#4b5563')
		fig.tight_layout(rect=[0,0.025,1,1])
		fig.savefig(path,bbox_inches='tight',facecolor='white',dpi=PLOT_DPI)
		plt.close(fig)
		return path

	def plot_pysyd_echelle_suite(self, profile_ids=(224,251,91)):
		"""Generate a small reference/intermediate/late pySYD echelle set."""
		return [self.plot_pysyd_echelle(profile_id) for profile_id in profile_ids]

	def cutoff_peak_table(self, spectrum, center, dnu, window_factor=3.0, plot_factor=5.0):
		lower_window=center-(window_factor*dnu)
		upper_window=center+(window_factor*dnu)
		lower_plot=max(0.0,center-(plot_factor*dnu))
		upper_plot=center+(plot_factor*dnu)
		local=spectrum[
			(spectrum['frequency']>=lower_plot) &
			(spectrum['frequency']<=upper_plot)
		].copy()
		if local.empty:
			return local,pd.DataFrame(),lower_window,upper_window,lower_plot,upper_plot

		frequency_values=local['frequency'].to_numpy()
		power_values=local['power_norm'].to_numpy()
		frequency_step=np.median(np.diff(frequency_values)) if len(frequency_values)>1 else dnu
		min_distance=max(1,int(0.45*dnu/max(frequency_step,np.finfo(float).eps)))
		peak_indices,_=find_peaks(
			power_values,
			prominence=0.02*np.nanmax(power_values),
			distance=min_distance
		)
		peaks=pd.DataFrame({
			'frequency':frequency_values[peak_indices],
			'power_norm':power_values[peak_indices]
		})
		if not peaks.empty:
			peaks['inside_window']=peaks['frequency'].between(lower_window,upper_window)
		return local,peaks,lower_window,upper_window,lower_plot,upper_plot

	def plot_pysyd_cutoff_profile(self, profile_id, window_factor=3.0, plot_factor=5.0):
		fdnu_data=self.standard_pysyd_fdnu_data()
		profile_id=int(profile_id)
		profile_rows=fdnu_data[fdnu_data['ID'].astype(int)==profile_id]
		if profile_rows.empty:
			raise ValueError(f'Profile {profile_id} is not in the pySYD synchronized table.')
		row=profile_rows.iloc[0]
		spectrum_path,spectrum=self.read_standard_pysyd_spectrum(profile_id)

		center=float(row['nu_max'])
		dnu=float(row['dnu_syd'])
		local,peaks,lower_window,upper_window,lower_plot,upper_plot=self.cutoff_peak_table(
			spectrum,
			center,
			dnu,
			window_factor=window_factor,
			plot_factor=plot_factor
		)
		if local.empty:
			raise ValueError(f'Profile {profile_id} has no spectrum data near nu_max={center}.')

		self.plot_dir.mkdir(parents=True,exist_ok=True)
		fig,(ax_context,ax_spectrum)=plt.subplots(
			2,
			1,
			figsize=(11,9),
			dpi=300,
			gridspec_kw={'height_ratios':[1.0,1.35]}
		)

		ax_context.scatter(
			fdnu_data['radius'],
			fdnu_data['fdnu_pysyd'],
			color='lightgray',
			s=18,
			alpha=0.75,
			edgecolors='none',
			label='standard pySYD grid'
		)
		ax_context.scatter(
			[row['radius']],
			[row['fdnu_pysyd']],
			color='crimson',
			s=42,
			zorder=5,
			label=f'profile {profile_id}'
		)
		ax_context.scatter(
			[row['radius']],
			[row['fdnu_pysyd']],
			facecolors='none',
			edgecolors='crimson',
			linewidths=2.8,
			s=620,
			zorder=4
		)
		ax_context.set_xscale('log')
		ax_context.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_context.set_xlabel(r'Radius ($R_\odot$)',fontsize=11)
		ax_context.set_ylabel(r'$f_{\Delta\nu}$ from pySYD',fontsize=11)
		ax_context.legend(loc='best',framealpha=0.9,fontsize=8)

		ax_spectrum.plot(
			local['frequency'],
			local['power_norm'],
			color='black',
			linewidth=1.2,
			alpha=0.9,
			label='normalized power'
		)
		ax_spectrum.axvspan(
			lower_window,
			upper_window,
			color='tab:blue',
			alpha=0.14,
			label=rf'pySYD window: $\nu_{{\max}}\pm {window_factor:g}\Delta\nu_{{\rm syd}}$'
		)
		ax_spectrum.axvline(center,color='tab:blue',linestyle='-',linewidth=1.4,label=r'$\nu_{\max}$')
		ax_spectrum.axvline(lower_window,color='tab:blue',linestyle='--',linewidth=1.2)
		ax_spectrum.axvline(upper_window,color='tab:blue',linestyle='--',linewidth=1.2)

		if not peaks.empty:
			inside=peaks[peaks['inside_window']]
			outside=peaks[~peaks['inside_window']]
			if not inside.empty:
				ax_spectrum.scatter(
					inside['frequency'],
					inside['power_norm'],
					color='tab:green',
					s=30,
					zorder=5,
					label='mode peak inside window'
				)
			if not outside.empty:
				ax_spectrum.scatter(
					outside['frequency'],
					outside['power_norm'],
					facecolors='none',
					edgecolors='crimson',
					linewidths=2.0,
					s=95,
					zorder=6,
					label='mode peak outside window'
				)
				for _,peak_row in outside.iterrows():
					ax_spectrum.axvline(
						peak_row['frequency'],
						color='crimson',
						linestyle=':',
						linewidth=1.0,
						alpha=0.75
					)

		stats_text=(
			f'Profile {profile_id}\n'
			f'Radius = {row["radius"]:.3f} R_sun\n'
			f'Mass = {row["mass"]/1.98847e33:.3f} M_sun\n'
			f'nu_max = {center:.4f} microHz\n'
			f'dnu_syd = {dnu:.4f} microHz\n'
			f'Window = [{lower_window:.4f}, {upper_window:.4f}] microHz'
		)
		ax_spectrum.text(
			0.02,
			0.96,
			stats_text,
			transform=ax_spectrum.transAxes,
			verticalalignment='top',
			horizontalalignment='left',
			fontsize=8.5,
			bbox=dict(boxstyle='round',facecolor='white',edgecolor='gray',alpha=0.92)
		)
		ax_spectrum.set_xlim(lower_plot,upper_plot)
		ax_spectrum.set_ylim(bottom=0)
		ax_spectrum.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_spectrum.set_xlabel(r'Frequency ($\mu$Hz)',fontsize=11)
		ax_spectrum.set_ylabel('Normalized power',fontsize=11)
		ax_spectrum.legend(loc='upper right',framealpha=0.9,fontsize=8)
		ax_spectrum.set_title(f'pySYD cutoff diagnostic for profile {profile_id}',fontsize=13)

		plt.tight_layout()
		output_path=self.path_pysyd_cutoff(profile_id)
		plt.savefig(output_path,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		print(f'Wrote pySYD cutoff plot for profile {profile_id} from {spectrum_path} to {output_path}')
		return output_path

	def pysyd_cutoff(self, profile_ids=(224,225,226,230,235,241,250,251)):
		output_paths=[]
		for profile_id in profile_ids:
			output_paths.append(self.plot_pysyd_cutoff_profile(profile_id))
		return output_paths

	def original_nonad_fdnu_table(self):
		data=self.standard_pysyd_fdnu_data().copy()
		data['ID']=data['ID'].astype(int)
		data=data.rename(columns={'fdnu_pysyd':'fdnu_nonad_pysyd'})

		old_lit=self.data['old_lit'].copy()
		old_lit['ID']=old_lit['PROFID'].astype(str).str.extract(r'(\d+)').astype(float)
		old_columns=['ID','fdnustdlit','dnu_model']
		old_lit=old_lit[[column for column in old_columns if column in old_lit.columns]].copy()
		for column in old_lit.columns:
			old_lit[column]=pd.to_numeric(old_lit[column],errors='coerce')
		old_lit=old_lit.dropna(subset=['ID']).copy()
		old_lit['ID']=old_lit['ID'].astype(int)
		old_lit=old_lit.rename(columns={
			'fdnustdlit':'fdnu_nonad_model',
			'dnu_model':'dnu_nonad_model'
		})

		asfgrid=self.data['asfgrid'].copy()
		asf_columns=['ID','asf_dnu','asf_fdnu']
		asfgrid=asfgrid[[column for column in asf_columns if column in asfgrid.columns]].copy()
		for column in asfgrid.columns:
			asfgrid[column]=pd.to_numeric(asfgrid[column],errors='coerce')
		asfgrid=asfgrid.dropna(subset=['ID']).copy()
		asfgrid['ID']=asfgrid['ID'].astype(int)
		asfgrid=asfgrid.rename(columns={'asf_dnu':'dnu_asfgrid','asf_fdnu':'fdnu_asfgrid_raw'})

		data=data.merge(old_lit,on='ID',how='left').merge(asfgrid,on='ID',how='left')
		data['fdnu_asfgrid_inverted']=1.0/data['fdnu_asfgrid_raw']
		data.loc[~np.isfinite(data['fdnu_asfgrid_inverted']),'fdnu_asfgrid_inverted']=np.nan
		data=data.sort_values('radius').reset_index(drop=True)
		data['previous_ID']=data['ID'].shift(1)
		data['previous_radius']=data['radius'].shift(1)
		data['fdnu_nonad_pysyd_previous']=data['fdnu_nonad_pysyd'].shift(1)
		data['fdnu_nonad_pysyd_jump_signed']=data['fdnu_nonad_pysyd']-data['fdnu_nonad_pysyd_previous']
		data['fdnu_nonad_pysyd_jump']=data['fdnu_nonad_pysyd_jump_signed'].abs()
		data['nonad_pysyd_model_offset']=data['fdnu_nonad_pysyd']-data['fdnu_nonad_model']
		data['nonad_pysyd_model_frac_offset']=data['nonad_pysyd_model_offset']/data['fdnu_nonad_model']
		data['nonad_pysyd_asfgrid_offset']=data['fdnu_nonad_pysyd']-data['fdnu_asfgrid_inverted']
		data['nonad_pysyd_asfgrid_frac_offset']=data['nonad_pysyd_asfgrid_offset']/data['fdnu_asfgrid_inverted']
		score_columns=['fdnu_nonad_pysyd_jump','nonad_pysyd_model_offset','nonad_pysyd_asfgrid_offset']
		data['internal_score']=0.0
		for column in score_columns:
			data['internal_score']+=data[column].abs().rank(pct=True).fillna(0.0)
		return data

	def original_nonad_cutoff_diagnostics(self, row, window_factor=3.0, plot_factor=5.0):
		diagnostics={
			'cutoff_spectrum_found':False,
			'cutoff_lower_window':np.nan,
			'cutoff_upper_window':np.nan,
			'cutoff_peak_count':np.nan,
			'cutoff_peaks_inside':np.nan,
			'cutoff_peaks_left_near':np.nan,
			'cutoff_peaks_right_near':np.nan,
			'cutoff_near_edge_peak_count':np.nan,
			'cutoff_left_edge_gap_dnu':np.nan,
			'cutoff_right_edge_gap_dnu':np.nan,
			'cutoff_min_edge_gap_dnu':np.nan,
			'cutoff_clipped_side':'missing',
			'cutoff_spectrum_path':''
		}
		try:
			spectrum_path,spectrum=self.read_standard_pysyd_spectrum(int(row['ID']))
			local,peaks,lower_window,upper_window,_,_=self.cutoff_peak_table(
				spectrum,
				float(row['nu_max']),
				float(row['dnu_syd']),
				window_factor=window_factor,
				plot_factor=plot_factor
			)
		except (FileNotFoundError,ValueError,KeyError) as exc:
			diagnostics['cutoff_error']=str(exc)
			return diagnostics
		if local.empty:
			diagnostics['cutoff_error']='empty local spectrum'
			return diagnostics

		left_gap=np.nan
		right_gap=np.nan
		if not peaks.empty:
			inside=peaks[peaks['inside_window']]
			left_outside=peaks[peaks['frequency']<lower_window]
			right_outside=peaks[peaks['frequency']>upper_window]
			left_near=left_outside[left_outside['frequency']>=lower_window-float(row['dnu_syd'])]
			right_near=right_outside[right_outside['frequency']<=upper_window+float(row['dnu_syd'])]
			if not left_outside.empty:
				left_gap=(lower_window-left_outside['frequency'].max())/float(row['dnu_syd'])
			if not right_outside.empty:
				right_gap=(right_outside['frequency'].min()-upper_window)/float(row['dnu_syd'])
		else:
			inside=pd.DataFrame()
			left_near=pd.DataFrame()
			right_near=pd.DataFrame()

		near_count=len(left_near)+len(right_near)
		if len(left_near)>0 and len(right_near)>0:
			clipped_side='both'
		elif len(left_near)>0:
			clipped_side='left'
		elif len(right_near)>0:
			clipped_side='right'
		else:
			clipped_side='none'
		edge_gaps=[gap for gap in [left_gap,right_gap] if np.isfinite(gap)]
		diagnostics.update({
			'cutoff_spectrum_found':True,
			'cutoff_lower_window':lower_window,
			'cutoff_upper_window':upper_window,
			'cutoff_peak_count':len(peaks),
			'cutoff_peaks_inside':len(inside),
			'cutoff_peaks_left_near':len(left_near),
			'cutoff_peaks_right_near':len(right_near),
			'cutoff_near_edge_peak_count':near_count,
			'cutoff_left_edge_gap_dnu':left_gap,
			'cutoff_right_edge_gap_dnu':right_gap,
			'cutoff_min_edge_gap_dnu':min(edge_gaps) if edge_gaps else np.nan,
			'cutoff_clipped_side':clipped_side,
			'cutoff_spectrum_path':str(spectrum_path)
		})
		return diagnostics

	def original_nonad_internal_diagnostics_table(self):
		table=self.original_nonad_fdnu_table()
		records=[]
		for _,row in table.iterrows():
			record=row.to_dict()
			record.update(self.original_nonad_cutoff_diagnostics(row))
			records.append(record)
		diagnostics=pd.DataFrame(records)
		if not diagnostics.empty:
			diagnostics['internal_score']+=diagnostics['cutoff_peaks_inside'].rank(pct=True,ascending=False).fillna(0.0)
		return diagnostics.sort_values('radius').reset_index(drop=True)

	def select_original_nonad_cases(self, table, count=6):
		candidates=table.dropna(subset=['fdnu_nonad_pysyd_jump']).copy()
		if candidates.empty:
			return candidates
		return candidates.sort_values('internal_score',ascending=False).head(count).copy()

	def plot_original_nonad_internal_summary(self, table, cases):
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		fig,axes=plt.subplots(2,2,figsize=(13,10),dpi=300)
		ax_fdnu,ax_jump,ax_offset,ax_window=axes.ravel()

		ax_fdnu.scatter(table['radius'],table['fdnu_nonad_model'],color='black',s=12,alpha=0.55,edgecolors='none',label='original nonad model')
		ax_fdnu.scatter(table['radius'],table['fdnu_nonad_pysyd'],color='tab:blue',s=18,alpha=0.65,edgecolors='none',label='pySYD on nonad spectra')
		ax_fdnu.scatter(table['radius'],table['fdnu_asfgrid_inverted'],color='tab:green',s=16,alpha=0.55,edgecolors='none',label='ASFGrid inverted')
		if not cases.empty:
			ax_fdnu.scatter(cases['radius'],cases['fdnu_nonad_pysyd'],facecolors='none',edgecolors='crimson',s=130,linewidths=2.0,label='selected cases')
		ax_fdnu.set_xscale('log')
		ax_fdnu.set_xlabel(r'Radius ($R_\odot$)')
		ax_fdnu.set_ylabel(r'$f_{\Delta\nu}$')
		ax_fdnu.set_title('Original nonad sequence')
		ax_fdnu.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_fdnu.legend(loc='best',fontsize=8,framealpha=0.9)

		scatter=ax_jump.scatter(
			table['radius'],
			table['fdnu_nonad_pysyd_jump'],
			c=table['cutoff_peaks_inside'],
			cmap='viridis',
			s=34,
			alpha=0.85,
			edgecolors='black',
			linewidths=0.2
		)
		fig.colorbar(scatter,ax=ax_jump,label='peaks inside pySYD window')
		ax_jump.set_xscale('log')
		ax_jump.set_xlabel(r'Radius ($R_\odot$)')
		ax_jump.set_ylabel(r'local pySYD $f_{\Delta\nu}$ jump')
		ax_jump.set_title('Internal pySYD jump size')
		ax_jump.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)

		ax_offset.axhline(0.0,color='black',linestyle=':',linewidth=1.0)
		ax_offset.scatter(table['radius'],table['nonad_pysyd_model_frac_offset'],color='tab:purple',s=18,alpha=0.65,edgecolors='none',label='vs nonad model')
		ax_offset.scatter(table['radius'],table['nonad_pysyd_asfgrid_frac_offset'],color='tab:green',s=18,alpha=0.55,edgecolors='none',label='vs ASFGrid inverted')
		if not cases.empty:
			ax_offset.scatter(cases['radius'],cases['nonad_pysyd_model_frac_offset'],facecolors='none',edgecolors='crimson',s=110,linewidths=1.7)
		ax_offset.set_xscale('log')
		ax_offset.set_xlabel(r'Radius ($R_\odot$)')
		ax_offset.set_ylabel('fractional offset of pySYD nonad')
		ax_offset.set_title('Offsets within the original sequence')
		ax_offset.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_offset.legend(loc='best',fontsize=8,framealpha=0.9)

		window_scatter=ax_window.scatter(
			table['cutoff_peaks_inside'],
			table['nonad_pysyd_model_frac_offset'].abs(),
			c=table['radius'],
			cmap='plasma',
			s=34,
			alpha=0.85,
			edgecolors='black',
			linewidths=0.2
		)
		fig.colorbar(window_scatter,ax=ax_window,label=r'Radius ($R_\odot$)')
		ax_window.set_xlabel('peaks inside pySYD window')
		ax_window.set_ylabel(r'$|f_{\Delta\nu,\rm pySYD}-f_{\Delta\nu,\rm model}|/f_{\Delta\nu,\rm model}$')
		ax_window.set_title('Window content vs pySYD offset')
		ax_window.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)

		plt.tight_layout()
		plt.savefig(self.path_nonad_pysyd_internal_summary,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return self.path_nonad_pysyd_internal_summary

	def plot_original_nonad_internal_case(self, row, table, window_factor=3.0, plot_factor=5.0):
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		fig,axes=plt.subplots(2,2,figsize=(13,10),dpi=300)
		ax_context,ax_spectrum,ax_zoom,ax_bar=axes.ravel()

		ax_context.scatter(table['radius'],table['fdnu_nonad_model'],color='black',s=12,alpha=0.45,edgecolors='none',label='original nonad model')
		ax_context.scatter(table['radius'],table['fdnu_nonad_pysyd'],color='tab:blue',s=16,alpha=0.6,edgecolors='none',label='pySYD on nonad spectra')
		ax_context.scatter(table['radius'],table['fdnu_asfgrid_inverted'],color='tab:green',s=14,alpha=0.5,edgecolors='none',label='ASFGrid inverted')
		ax_context.scatter([row['radius']],[row['fdnu_nonad_pysyd']],facecolors='none',edgecolors='crimson',s=160,linewidths=2.2,zorder=6)
		ax_context.set_xscale('log')
		ax_context.set_xlabel(r'Radius ($R_\odot$)')
		ax_context.set_ylabel(r'$f_{\Delta\nu}$')
		ax_context.set_title(f'Nonad pySYD case profile {int(row["ID"])}')
		ax_context.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_context.legend(loc='best',fontsize=8,framealpha=0.9)

		try:
			spectrum_path,spectrum=self.read_standard_pysyd_spectrum(int(row['ID']))
			local,peaks,lower_window,upper_window,lower_plot,upper_plot=self.cutoff_peak_table(
				spectrum,
				float(row['nu_max']),
				float(row['dnu_syd']),
				window_factor=window_factor,
				plot_factor=plot_factor
			)
			ax_spectrum.plot(local['frequency'],local['power_norm'],color='black',linewidth=1.1,label='normalized power')
			ax_spectrum.axvspan(lower_window,upper_window,color='tab:blue',alpha=0.14,label='pySYD window')
			ax_spectrum.axvline(float(row['nu_max']),color='tab:blue',linestyle='-',linewidth=1.2,label=r'$\nu_{\max}$')
			ax_spectrum.axvline(lower_window,color='tab:blue',linestyle='--',linewidth=1.0)
			ax_spectrum.axvline(upper_window,color='tab:blue',linestyle='--',linewidth=1.0)
			if not peaks.empty:
				inside=peaks[peaks['inside_window']]
				outside=peaks[~peaks['inside_window']]
				if not inside.empty:
					ax_spectrum.scatter(inside['frequency'],inside['power_norm'],color='tab:green',s=28,zorder=5,label='inside peak')
				if not outside.empty:
					ax_spectrum.scatter(outside['frequency'],outside['power_norm'],facecolors='none',edgecolors='crimson',s=80,linewidths=1.5,zorder=6,label='outside peak')
			ax_spectrum.set_xlim(lower_plot,upper_plot)
			ax_spectrum.set_ylim(bottom=0)
		except (FileNotFoundError,ValueError,KeyError) as exc:
			ax_spectrum.text(0.5,0.5,f'No spectrum\n{exc}',ha='center',va='center',transform=ax_spectrum.transAxes)
		ax_spectrum.set_xlabel(r'Frequency ($\mu$Hz)')
		ax_spectrum.set_ylabel('Normalized power')
		ax_spectrum.set_title('pySYD window on nonad spectrum')
		ax_spectrum.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_spectrum.legend(loc='best',fontsize=8,framealpha=0.9)

		row_index=int(table.index[table['ID']==row['ID']][0])
		zoom=table.iloc[max(0,row_index-5):min(len(table),row_index+6)].copy()
		ax_zoom.plot(zoom['radius'],zoom['fdnu_nonad_model'],color='black',marker='o',markersize=3,linewidth=1.0,label='original nonad model')
		ax_zoom.plot(zoom['radius'],zoom['fdnu_nonad_pysyd'],color='tab:blue',marker='o',markersize=3,linewidth=1.0,label='pySYD on nonad spectra')
		ax_zoom.plot(zoom['radius'],zoom['fdnu_asfgrid_inverted'],color='tab:green',marker='o',markersize=3,linewidth=1.0,label='ASFGrid inverted')
		ax_zoom.scatter([row['radius']],[row['fdnu_nonad_pysyd']],facecolors='none',edgecolors='crimson',s=140,linewidths=2.0,zorder=6)
		ax_zoom.set_xlabel(r'Radius ($R_\odot$)')
		ax_zoom.set_ylabel(r'$f_{\Delta\nu}$')
		ax_zoom.set_title('Local radius sequence')
		ax_zoom.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_zoom.legend(loc='best',fontsize=8,framealpha=0.9)

		labels=['nonad model','pySYD nonad','ASFGrid inv.']
		values=[row['fdnu_nonad_model'],row['fdnu_nonad_pysyd'],row['fdnu_asfgrid_inverted']]
		colors=['black','tab:blue','tab:green']
		ax_bar.bar(labels,values,color=colors,alpha=0.78)
		ax_bar.axhline(row['fdnu_nonad_model'],color='black',linestyle=':',linewidth=1.0)
		stats_text=(
			f'Profile {int(row["ID"])}\n'
			f'R = {row["radius"]:.4g} Rsun\n'
			f'nu_max = {row["nu_max"]:.4g} uHz\n'
			f'dnu_syd = {row["dnu_syd"]:.4g} uHz\n'
			f'pySYD jump = {row["fdnu_nonad_pysyd_jump"]:.4g}\n'
			f'window peaks = {int(row["cutoff_peaks_inside"])}'
		)
		ax_bar.text(
			0.03,
			0.97,
			stats_text,
			transform=ax_bar.transAxes,
			verticalalignment='top',
			horizontalalignment='left',
			fontsize=8.5,
			bbox=dict(boxstyle='round',facecolor='white',edgecolor='gray',alpha=0.92)
		)
		ax_bar.set_ylabel(r'$f_{\Delta\nu}$')
		ax_bar.set_title('Same nonad model, different measurements')
		ax_bar.grid(True,axis='y',color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)

		plt.tight_layout()
		output_path=self.path_nonad_pysyd_internal_case(row['ID'])
		plt.savefig(output_path,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return output_path

	def nonad_reference_long_table(self, original_table=None):
		if original_table is None:
			original_table=self.original_nonad_fdnu_table()
		series=[]
		original_sources=[
			('Original nonad model','original_nonad','fdnu_nonad_model'),
			('pySYD on nonad spectra','original_nonad','fdnu_nonad_pysyd'),
			('Original ASFGrid inverted','original_reference','fdnu_asfgrid_inverted')
		]
		for source,family,column in original_sources:
			if column in original_table.columns:
				part=original_table[['radius',column]].copy()
				part=part.rename(columns={column:'fdnu'})
				part['source']=source
				part['family']=family
				series.append(part)

		li_table=self.li_l0_fdnu_table()
		li_sources=[
			('Li reference','li_reference','fdnu_li'),
			('pySYD on Li spectra','li_reference','fdnu_pysyd'),
			('Li ASFGrid inverted','li_reference','fdnu_asfgrid_inverted')
		]
		for source,family,column in li_sources:
			if column in li_table.columns:
				part=li_table[['radius',column]].copy()
				part=part.rename(columns={column:'fdnu'})
				part['source']=source
				part['family']=family
				series.append(part)
		long_table=pd.concat(series,ignore_index=True)
		long_table=long_table.replace([np.inf,-np.inf],np.nan).dropna(subset=['radius','fdnu'])
		long_table=long_table[(long_table['radius']>0) & (long_table['fdnu']>0)].copy()
		return long_table

	def binned_nonad_reference_comparison(self, long_table, bin_count=24):
		radius_min=max(long_table['radius'].min(),np.finfo(float).tiny)
		radius_max=long_table['radius'].max()
		bins=np.geomspace(radius_min,radius_max,bin_count+1)
		long_table=long_table.copy()
		long_table['radius_bin']=pd.cut(long_table['radius'],bins=bins,include_lowest=True)
		binned=long_table.groupby(['source','family','radius_bin'],observed=True).agg(
			radius_median=('radius','median'),
			fdnu_median=('fdnu','median'),
			fdnu_q16=('fdnu',lambda value: value.quantile(0.16)),
			fdnu_q84=('fdnu',lambda value: value.quantile(0.84)),
			count=('fdnu','size')
		).reset_index()
		return binned[binned['count']>0].copy()

	def plot_nonad_reference_comparison(self, original_table, long_table, binned):
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		fig,axes=plt.subplots(2,2,figsize=(13,10),dpi=300)
		ax_raw,ax_binned,ax_offsets,ax_zoom=axes.ravel()
		colors={
			'Original nonad model':'black',
			'pySYD on nonad spectra':'tab:blue',
			'Original ASFGrid inverted':'tab:green',
			'Li reference':'gray',
			'pySYD on Li spectra':'tab:orange',
			'Li ASFGrid inverted':'mediumseagreen'
		}
		for source,source_data in long_table.groupby('source'):
			ax_raw.scatter(
				source_data['radius'],
				source_data['fdnu'],
				color=colors.get(source,'gray'),
				s=7 if 'Li' in source else 18,
				alpha=0.20 if 'Li' in source else 0.62,
				edgecolors='none',
				label=source
			)
		ax_raw.set_xscale('log')
		ax_raw.set_xlabel(r'Radius ($R_\odot$)')
		ax_raw.set_ylabel(r'$f_{\Delta\nu}$')
		ax_raw.set_title('Cross-dataset overlay, not model-matched')
		ax_raw.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_raw.legend(loc='best',fontsize=7.5,framealpha=0.9)

		for source,source_data in binned.groupby('source'):
			source_data=source_data.sort_values('radius_median')
			ax_binned.plot(
				source_data['radius_median'],
				source_data['fdnu_median'],
				color=colors.get(source,'gray'),
				linewidth=1.5,
				label=source
			)
			ax_binned.fill_between(
				source_data['radius_median'],
				source_data['fdnu_q16'],
				source_data['fdnu_q84'],
				color=colors.get(source,'gray'),
				alpha=0.12
			)
		ax_binned.set_xscale('log')
		ax_binned.set_xlabel(r'Radius ($R_\odot$)')
		ax_binned.set_ylabel(r'median $f_{\Delta\nu}$')
		ax_binned.set_title('Binned comparison')
		ax_binned.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_binned.legend(loc='best',fontsize=7.5,framealpha=0.9)

		ax_offsets.axhline(0.0,color='black',linestyle=':',linewidth=1.0)
		ax_offsets.scatter(
			original_table['radius'],
			original_table['nonad_pysyd_model_frac_offset'],
			color='tab:blue',
			s=18,
			alpha=0.65,
			edgecolors='none',
			label='nonad pySYD vs nonad model'
		)
		ax_offsets.scatter(
			original_table['radius'],
			original_table['nonad_pysyd_asfgrid_frac_offset'],
			color='tab:green',
			s=18,
			alpha=0.55,
			edgecolors='none',
			label='nonad pySYD vs ASFGrid inv.'
		)
		li_table=self.li_l0_fdnu_table()
		li_table['li_pysyd_frac_offset']=(li_table['fdnu_pysyd']-li_table['fdnu_li'])/li_table['fdnu_li']
		li_table['li_asfgrid_frac_offset']=(li_table['fdnu_asfgrid_inverted']-li_table['fdnu_li'])/li_table['fdnu_li']
		ax_offsets.scatter(
			li_table['radius'],
			li_table['li_pysyd_frac_offset'],
			color='tab:orange',
			s=4,
			alpha=0.10,
			edgecolors='none',
			label='Li pySYD vs Li reference'
		)
		ax_offsets.scatter(
			li_table['radius'],
			li_table['li_asfgrid_frac_offset'],
			color='mediumseagreen',
			s=4,
			alpha=0.10,
			edgecolors='none',
			label='Li ASFGrid inv. vs Li reference'
		)
		ax_offsets.set_xscale('log')
		ax_offsets.set_xlabel(r'Radius ($R_\odot$)')
		ax_offsets.set_ylabel('fractional offset')
		ax_offsets.set_title('Within-source offsets')
		ax_offsets.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_offsets.legend(loc='best',fontsize=7.2,framealpha=0.9)

		zoom_binned=binned[binned['radius_median']>=30].copy()
		for source,source_data in zoom_binned.groupby('source'):
			source_data=source_data.sort_values('radius_median')
			ax_zoom.plot(
				source_data['radius_median'],
				source_data['fdnu_median'],
				color=colors.get(source,'gray'),
				linewidth=1.5,
				label=source
			)
		ax_zoom.set_xscale('log')
		ax_zoom.set_xlabel(r'Radius ($R_\odot$)')
		ax_zoom.set_ylabel(r'median $f_{\Delta\nu}$')
		ax_zoom.set_title(r'Luminous giant zoom ($R \geq 30 R_\odot$)')
		ax_zoom.grid(True,color='lightgray',linestyle='--',linewidth=0.7,alpha=0.8)
		ax_zoom.legend(loc='best',fontsize=7.2,framealpha=0.9)

		plt.tight_layout()
		plt.savefig(self.path_nonad_reference_comparison_summary,bbox_inches='tight',facecolor='white',dpi=300)
		plt.close()
		return self.path_nonad_reference_comparison_summary

	def nonad_vs_reference_diagnostics(self, case_count=6):
		internal_table=self.original_nonad_internal_diagnostics_table()
		cases=self.select_original_nonad_cases(internal_table,count=case_count)
		self.plot_dir.mkdir(parents=True,exist_ok=True)
		internal_table.to_csv(self.path_nonad_pysyd_internal_table,index=False)
		internal_summary=self.plot_original_nonad_internal_summary(internal_table,cases)
		case_paths=[]
		for _,case_row in cases.iterrows():
			case_paths.append(self.plot_original_nonad_internal_case(case_row,internal_table))

		long_table=self.nonad_reference_long_table(internal_table)
		long_table.to_csv(self.path_nonad_reference_long_table,index=False)
		binned=self.binned_nonad_reference_comparison(long_table)
		binned.to_csv(self.path_nonad_reference_comparison_table,index=False)
		comparison_summary=self.plot_nonad_reference_comparison(internal_table,long_table,binned)
		print(f'Wrote original nonad pySYD diagnostics to {self.path_nonad_pysyd_internal_table}')
		print(f'Wrote original nonad pySYD summary plot to {internal_summary}')
		for case_path in case_paths:
			print(case_path)
		print(f'Wrote nonad/reference long comparison to {self.path_nonad_reference_long_table}')
		print(f'Wrote nonad/reference binned comparison to {self.path_nonad_reference_comparison_table}')
		print(f'Wrote nonad/reference comparison plot to {comparison_summary}')
		return {
			'internal_table':self.path_nonad_pysyd_internal_table,
			'internal_summary':internal_summary,
			'cases':case_paths,
			'comparison_long_table':self.path_nonad_reference_long_table,
			'comparison_table':self.path_nonad_reference_comparison_table,
			'comparison_summary':comparison_summary
		}


if __name__ == '__main__':
	tst_plt=all_plots()
	tst_plt.plot_asf_self_interp()

# also I think that the NUMAX is what is changing them 4.8.26, it is not at least I dont think so
# numax and nu_max are the same in li data
# sharma is both asf grid and apokasc (basically) 
