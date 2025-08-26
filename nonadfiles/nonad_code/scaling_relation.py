'''
this routine will take numax, and dnu to yield mass, radius, and g
calling sequence
./scaling_relation.py list.globalpars output
to return a file output with columns
file | M | Merr | g | gerr

otherwise, can be called like so:
from scaling_relation import scaling_relation
globalpars_file = 'list.globalpars'
output_file = 'asteroseismology_has_the_answers.out'
scaling_relation(globalpars_file, output_file)

the errors that are returned are from scaling relations
from hekker+2013, a grid-based approach only improves logg errors to 0.171% for 1% change in numax, compared to 0.4% change in g for 1% in numax 
from scaling relations

scaling relations used are
numax = g/sqrt(Teff)
dnu = sqrt(mean density)
scaled to solar values

the other scaling relations are from kallinger+2010

r = numax*(dnu)^(-2)*(teff)^1/2

m = r^3*dnu^2

solar parameters used are
dnu = 134.92 muhz (toutain&frohlich1992) pm 0.02muhz
numax = 3050 muhz
rsun = 6.958e10
msun = 1.988435e33
gsun = 2.7395e4
'''
from __future__ import print_function
import sys
import numpy as np
global teff_sun, logg_sun, numax_sun, dnu_sun, onepc, oneau, sigma, l_sun, dnu_sun_emp, numax_sun_emp
# Gravitational constant, in cgs
G = 6.67408e-8
# stefan-boltzmann constant, in cgs
sigma = 5.670367e-5
oneau = 1.495978707e13

onepc = 206265.*oneau
l_sun = 3.826e33
teff_sun = 5772. #(mamajek+ 2015)
numax_sun = 3050.

################################################################
# empirical values computed by marc (only actually dnu and numax, but i made teff and logg empirical variables just for completeness)
teff_sun_emp = 5777.
dnu_sun_emp = 135.146
numax_sun_emp = 3076.
logg_sun_emp = 2.7413e4
################################################################


dnu_sun = 134.92
R_sun = 6.958e10
M_sun = 1.988435e33
# logg_sun = 2.7395e4 # mamajek+2015 would be GM = 1.3271244e20 m^3s^{-2} / (R = 6.957e8 m)^2 = 2.74200e4
logg_sun = 2.74200e4
# below, the sharma value is 27413.240395773944 because that is 10.**(4.43796037457)
rho_sun = 1.4085514
L_sun = 3.828e33

def freq_dyn(M, R, sharma=False):
    '''
    Returns the dynamical frequency in muHz
    Inputs
    M : ndarray | float
     mass in solar masses
    R : ndarray | float
     radius in solar radii
    Outputs
    dyn_freq : ndarray | float
     \sqrt{\frac{R^3}/{GM}}, in muHz
    '''
    return np.sqrt(rho_sun*(M/R**3)*G)*1e6

def logg(numax, teff, sharma=False, emp=False, refs=[None, None]):
    '''
    logg ( = log10(g) in cgs)
    '''

    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4
    if refs[0] is not None:
        numax_sun = refs[0]
    if refs[1] is not None:
        dnu_sun = refs[1]

    return np.log10(numax) - np.log10(numax_sun) + 0.5*np.log10(teff) - 0.5*np.log10(teff_sun) + np.log10(logg_sun)

def numax(logg, teff, sharma=False, emp=False):
    '''
    return an expected numax given a log g and teff
    Inputs
    logg : float
     log10 surface gravity [cgs]
    teff : float
     effective temperature [K]
    [ emp : bool ]
    Outputs
    numax : float
     frequency of maximum oscillation [muhz]
     '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return 10.**(logg)/(logg_sun)*numax_sun*(teff/teff_sun)**(-0.5)


def numaxerr(logg, teff, logg_err, teff_err, sharma=False, emp=False):
    '''
    return an expected numax uncertainty given a log g and teff and their uncertainties
    Inputs
    logg : float
     log10 surface gravity [cgs]
    teff : float
     effective temperature [K]
    logg_err : float
    teff_err : float
    emp : bool
    Outputs
    numax : float
     frequency of maximum oscillation [muhz]
     '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return numax(logg, teff, sharma=sharma, emp=emp)*np.sqrt ( (teff_err/teff*0.5)**2 + (np.log(10.)*logg_err*logg)**2)

def dnu_MR(M, R, sharma=False):
    '''
    return an expected dnu given a mass and a radius
    Inputs
    M : float
     Mass [cgs]
    R : float
     Radius [cgs]
    Outputs
    numax : float
     large frequency separation [muhz]
     '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4
    rho = M/R**3
    return np.sqrt(rho/rho_sun)*dnu_sun

def loggerr(numax, teff, numax_err, teff_err, sharma=False, emp=False):
    '''
    error on logg ( = log10(g) in cgs)
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return np.sqrt((numax_err/numax)**2 + (0.5*teff_err/teff)**2)*np.log10(np.e)
# JCZ 290117
# added to convert galaxia log gs to a numax
def loggteff2numax(logg, teff, sharma=False):
    '''
    Inputs
    logg : float | ndarray
     log cgs surface gravity
    teff : float | ndarary
     effective temperature, K
    Outputs
     numax, muhz
     with reference value from sun of 3050muhz
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return 10.**(logg - np.log10(logg_sun) - (0.5*np.log10(teff) - 0.5*np.log10(teff_sun)))*numax_sun



def rho(dnu, sharma=False):
    '''
    mean density in linear cgs units
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return (dnu/dnu_sun)**2*rho_sun

def rhoerr(dnu, dnu_err, sharma=False):
    '''
    error of mean density in linear cgs units
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4

    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4
    return 2*dnu_err/dnu*rho(dnu, sharma=sharma)



def M(numax, teff, dnu, sharma=False, emp=False, refs=[None, None]):
    '''
    log10(M)
    Inputs
    [ refs : [numax_solar_ref, dnu_solar_ref], list of floats ]
     if you put in either numax_solar_ref or dnu_solar_ref, that will be used, no matter what sharma and emp are.
    [ emp : bool ]
     if True, use the empirical numax and dnu that marc calculated in APOKASC. Default False.
    
    '''
    
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4
    if refs[0] is not None:
        numax_sun = refs[0]
    if refs[1] is not None:
        dnu_sun = refs[1]
    # return 3.*np.log10(numax) - 3.*np.log10(numax_sun) - 4.*np.log10(dnu) + 4.*np.log10(dnu_sun) + 3./2.*np.log10(teff) - 3./2.*np.log10(teff_sun) + np.log10(M_sun)
    return 3.*(R(numax, teff, dnu, sharma=sharma, emp=emp)) - 3.*np.log10(R_sun) + 2.*np.log10(dnu) - 2.*np.log10(dnu_sun) + np.log10(M_sun)

def Merr(numax, teff, dnu, numax_err, teff_err, dnu_err, sharma=False, emp=False, refs=[None, None]):
    '''
    error of log10(M)
    Inputs
    [ refs : [numax_solar_ref, dnu_solar_ref], list of floats ]
     if you put in either numax_solar_ref or dnu_solar_ref, that will be used, no matter what sharma and emp are.    
    [ emp : bool ]
     if True, use the empirical numax and dnu that marc calculated in APOKASC. Default False.

    '''


    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return np.sqrt((((3.*Rerr(numax, teff, dnu, numax_err, teff_err, dnu_err, sharma=sharma, emp=emp, refs=refs)/np.log10(np.e))**2  + (2.*(dnu_err)/(dnu))**2)))*np.log10(np.e)
def M2(numax, dnu, teff, sharma=False, emp=False, refs=[None,None]):
    return 10.**(M(numax, teff, dnu, sharma=sharma, emp=emp, refs=refs))/M_sun
def Merr2(numax, dnu, teff, numax_err, dnu_err, teff_err, sharma=False, emp=False, refs=[None,None]):
    return np.log(10.0)*Merr(numax, teff, dnu, numax_err, teff_err, dnu_err, sharma=sharma, emp=emp, refs=refs)*M2(numax, dnu, teff, sharma=sharma, emp=emp, refs=refs)
def R(numax, teff, dnu, sharma=False, emp=False, refs=[None, None]):
    '''
    log10(R) in cgs
    Inputs
    [ refs : [numax_solar_ref, dnu_solar_ref], list of floats ]
     if you put in either numax_solar_ref or dnu_solar_ref, that will be used, no matter what sharma and emp are.
    [ emp : bool ]
     if True, use the empirical numax and dnu that marc calculated in APOKASC. Default False.

    '''
    import numpy as np
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4
    if refs[0] is not None:
        numax_sun = refs[0]
    if refs[1] is not None:
        dnu_sun = refs[1]
    return np.log10(numax) - np.log10(numax_sun) -2.*(np.log10(dnu/dnu_sun)) + 0.5*np.log10(teff/teff_sun) + np.log10(R_sun)

def Rerr(numax, teff, dnu, numax_err, teff_err, dnu_err, sharma=False, emp=False, refs=[None, None]):
    '''
    error on log10(R) in cgs
    Inputs
    [ refs : [numax_solar_ref, dnu_solar_ref], list of floats ]
     if you put in either numax_solar_ref or dnu_solar_ref, that will be used, no matter what sharma and emp are.
    [ emp : bool ]
     if True, use the empirical numax and dnu that marc calculated in APOKASC. Default False.

    '''


    
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4
    if refs[0] is not None:
        numax_sun = refs[0]
    if refs[1] is not None:
        dnu_sun = refs[1]

    return np.sqrt((0.5*teff_err/teff)**2 + (2.*dnu_err/dnu)**2 + (numax_err/numax)**2)*np.log10(np.e)

def L(numax, teff, dnu, d=1., sharma=False):
    '''
    Luminosity in cgs
    NB : NOT LOG10(L) BUT LINEAR
    ASSUMES DISTANCE OF 1pc
    Inputs
    [ d : float ]
    distance, in units of pc
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return 0.25 * R(numax, teff, dnu, sharma=sharma)**2 * sigma * teff**4 / (d*onepc)**2

def Lerr(numax, teff, dnu, numax_err, teff_err, dnu_err, d=1., derr=0., sharma=False):
    '''
    error in Luminosity, in cgs
    NB : assumes distance of 1pc and no error, by default
    Inputs
    [ d : float ]
    distance, in parsecs
    [ derr : float ]
    error on distance, in parsecs
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return np.sqrt((4*teff_err/teff)**2 + (2*derr/d)**2 + (2 * Rerr(numax, teff, dnu, numax_err, teff_err, dnu_err, sharma=sharma) / log10(np.e))**2)*L(numax, teff, dnu,d=d, sharma=sharma)



def derr_frac(numax, teff, dnu, l, numax_err, teff_err, dnu_err, l_err, sharma=False):
    '''
    per cent error in distance, given known luminosity, <l>, in cgs and error, <l_err>
    Outputs
    derr : float
     dimensionaless fractional error on distance, d_err/d
    Notes
    NOT ABSOLUTE ERROR
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return np.sqrt((l_err/l)**2 - (4.*teff_err/teff)**2 - (2.*Rerr(numax, teff, dnu, numax_err, teff_err, dnu_err, sharma=sharma) / log10(np.e))**2)/2.
    
def d(numax, teff, dnu, l, sharma=False):
    '''
    distance in pc
    Inputs
    Outputs
    d : float
     distance, in pc
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return np.sqrt(4.*l/R(numax, teff, dnu, sharma=sharma)**2/sigma/teff**4)/onepc

def derr(numax, teff, dnu, l, numax_err, teff_err, dnu_err, l_err, sharma=False):
    '''
    absolute error in distance, given known luminosity, <l>, in cgs and error, l_err>
    Outputs
    derr : float
     error on distance, derr, in pc
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return derr_frac(numax, teff, dnu, l, numax_err, teff_err, dnu_err, l_err, sharma=sharma)*d(numax, teff, dnu, l, sharma=sharma)


def d1(numax, teff, dnu, f=1., sharma=False, emp=False):
    '''
    distance in pc
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return 10.**(R(numax, teff, dnu, sharma=sharma)) * np.sqrt(sigma) * teff**2 / np.sqrt(f) /onepc

def d2(R, teff, f=1.0, sharma=False):
    '''
    distance in pc
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return 10.**(R) * np.sqrt(sigma) * teff**2 / np.sqrt(f) /onepc

def derr2(R, teff, Rerr, teff_err, f=1., ferr=0., sharma=False):
    '''
    error in flux, in cgs
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return np.sqrt((2*teff_err/teff)**2 + (0.5 * ferr/f)**2 + (Rerr / np.log10(np.e))**2)*d2(R, teff,f=f, sharma=sharma)


def R1(d, teff, f, sharma=False, emp=False):
    '''
    Given a distance (e.g., from Gaia), teff, and bolometric flux, returns LOG10(radius) of that star 
    Inputs
    d : float
     distance (in pc)
    teff : float
     effective temperature (K)
    f : float
     bolometric flux (cgs)
    Inputs
    [ emp : bool ]
     if True, use the empirical numax and dnu that marc calculated in APOKASC. Default False.
    Outputs
    R1 : float
     log10 radius, in cm
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return np.log10(d*onepc/np.sqrt(sigma)*np.sqrt(f)/teff**2)


def R1_err(d, teff, f, derr, tefferr, ferr, sharma=False, emp=False):
    '''
    return error on log10(R)
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return np.sqrt((0.5*ferr/f)**2 + (2.0*tefferr/teff) + (derr/d)**2)*np.log10(np.e)

def derr1(numax, teff, dnu, numax_err, teff_err, dnu_err, f=1., ferr=0., sharma=False, emp=False):
    '''
    error in flux, in cgs
    '''
    if sharma:
        teff_sun = 5777.
        dnu_sun = 135.1
        numax_sun = 3090.
        logg_sun = 2.7413e4
    elif emp:
        teff_sun = 5777.
        dnu_sun = 135.146
        numax_sun = 3076.
        logg_sun = 2.7413e4
    else:
        teff_sun = 5772. #(mamajek+ 2015)
        dnu_sun = 134.92
        numax_sun = 3050.
        logg_sun = 2.74200e4

    return np.sqrt((2*teff_err/teff)**2 + (0.5 * ferr/f)**2 + (Rerr(numax, teff, dnu, numax_err, teff_err, dnu_err, sharma=sharma) / np.log10(np.e))**2)*d1(numax, teff, dnu,f=f, sharma=sharma, emp=emp)
def R2(logg, mass):
    '''
    return log10(R), given logg and mass
    Inputs
    logg : float | ndarray
    mass : float | ndarray
    logg_err : float | ndarray
    mass_err : float | ndarray
    Outputs
    R1 : float
     log10 radius, in cm
    '''
    G = 6.67408e-8
    return np.log10(np.sqrt(G*mass*M_sun/10.**(logg)))

def R2_err(logg, mass, logg_err, mass_err):
    '''
    return error on log10(R), given logg and mass
    Inputs
    logg : float | ndarray
    mass : float | ndarray
    logg_err : float | ndarray
    mass_err : float | ndarray
    '''

    return np.sqrt((0.5*mass_err/mass)**2 + (logg_err/np.log10(np.e)*0.5)**2)*np.log10(np.e)

def logg_MR(M, R):
    '''
    given M and R in cgs, return log g
    Inputs
    M : ndarray | float
    R : ndarray | float
    Outputs
    log(10) g : ndarray | float
    '''
    return np.log10(G*M/R**2)
    

if __name__=='__main__':
    debug = False
    if debug:
        print ('found solar values of :')
        M_sun_ = M(numax_sun, teff_sun, dnu_sun)
        R_sun_ = R(numax_sun, teff_sun, dnu_sun)
        rho_sun_ = rho(dnu_sun)
        print( 'M : {}, R : {}, rho : {}'.format(M_sun_, R_sun_, rho_sun_))

    # read in numax and dnu
    data = pd.read_csv(globalpars_file, delimiter='\s+', engine='python', header=0)
    numax = data['numax']
    dnu = data['dnu']
    data['logg'] = logg(numax, teff)
    data['logg_e'] = loggerr(numax, teff, numax_err, teff_err)
    data['rho'] = rho(dnu)
    data['rho_e'] = rhoerr(dnu, dnu_err)
    data['M'] = M(numax, teff, dnu)
    data['M_e'] = Merr(numax, teff, dnu, numax_err, teff_err, dnu_err)
    data['R'] = R(numax, teff, dnu)
    data['R_e'] = Rerr(numax, teff, dnu, numax_err, teff_err, dnu_err)

    # np.savetxt((np.hstack(logg(
    gp_file = sys.argv[1]
