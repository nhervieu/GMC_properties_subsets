#import libraries
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from galaxies import Galaxy
import astropy.units as u
import powerlaw

#load file
mytable = Table.read('m100.co10.kkms_props_cprops_withSFR_subsets.fits')

##VARIABLE DEFINITIONS##

# Pull out the mass variable into a numpy array.
mass = mytable['MASS_EXTRAP'].data 

#Radius vs line width best fit line
R = np.arange(10,1000)
S = np.power(np.pi,1/2)*R/3.4
sigmav = np.sqrt(S)
one = np.arange(10000,10000000000,10000)

#Virial Mass best fit line
M_vir = 540*np.power(R,2)

#Luminous Mass best fit line
M_lum = 2000*np.power(np.power(M_vir/39,1.234567901)/130,0.8)

#Sigma_0, Mass Density, R_gal
sigma0 = mytable['VRMS_EXTRAP_DECONV']/np.sqrt(mytable['RADRMS_EXTRAP_DECONV'])  
M_den = mytable['MASS_EXTRAP']/(np.pi*np.power(mytable['RADRMS_EXTRAP_DECONV'],2)) 

#Star formation varialbes
Tdep = mytable['Tdep'].data
SFRSD = mytable['SFRSD'].data
mygalaxy = Galaxy("M100")
cpropstable = Table.read('m100.co10.kkms_props_cprops_withSFR_subsets.fits')
rgal=mygalaxy.radius(ra = cpropstable['XPOS'], dec = cpropstable['YPOS'])

#indices for subsets
centre_in = np.where(mytable['Nuclear'] == True)
arm_in = np.where(mytable['Arm'] == True)
interarm_in = np.where(mytable['Interarm'] == True)
outliers = np.where(mytable['Outliers'] == True)

#Calculate new depletion time
Tdep2 = M_den *1E6 / SFRSD

##PLOTS##
#Depletion time vs Surface Density
figure = plt.figure(figsize=(4.5,4))
plt.xlabel('$\Sigma\ ((M_{\odot})/pc^2)$')
plt.ylabel('$Depletion\ Time (yrs)$')
plt.loglog(M_den[centre_in],Tdep[centre_in],marker='v',c='b',linestyle='None')
plt.loglog(M_den[arm_in],Tdep[arm_in],marker='d',c='m',linestyle='None')
plt.loglog(M_den[interarm_in],Tdep[interarm_in],marker='o',c='g',linestyle='None')
plt.ylim(1E9, 1.5E10)
plt.xlim(10,1E4)
plt.tight_layout() 	
plt.savefig('Tdep_Mden_matplotlib.png')

#Depletion time vs Sigma_0
figure = plt.figure(figsize=(4.5,4))
plt.ylabel('$Depletion\ Time (yrs)$')
plt.xlabel('$\sigma_0$')
plt.loglog(sigma0[centre_in],Tdep[centre_in],marker='v',c='b',linestyle='None')
plt.loglog(sigma0[arm_in],Tdep[arm_in],marker='d',c='m',linestyle='None')
plt.loglog(sigma0[interarm_in],Tdep[interarm_in],marker='o',c='g',linestyle='None')
plt.ylim(1E9, 1.5E10)
plt.xlim(5E-2,4)
plt.tight_layout() 	
plt.savefig('Tdep_sigma0_matplotlib.png')

#Depletion time vs Surface Density
figure = plt.figure(figsize=(4.5,4))
plt.xlabel('$\Sigma\ ((M_{\odot})/pc^2)$')
plt.ylabel('$Depletion\ Time (actual) (yrs)$')
plt.loglog(M_den[centre_in],Tdep2[centre_in],marker='v',c='b',linestyle='None')
plt.loglog(M_den[arm_in],Tdep2[arm_in],marker='d',c='m',linestyle='None')
plt.loglog(M_den[interarm_in],Tdep2[interarm_in],marker='o',c='g',linestyle='None')
plt.ylim(4E8, 1E12)
plt.xlim(10,1E4)
plt.tight_layout() 	
plt.savefig('Tdep2_Mden_matplotlib.png')

#Depletion time vs Sigma_0
figure = plt.figure(figsize=(4.5,4))
plt.ylabel('$Depletion\ Time (actual) (yrs)$')
plt.xlabel('$\sigma_0$')
plt.loglog(sigma0[centre_in],Tdep2[centre_in],marker='v',c='b',linestyle='None')
plt.loglog(sigma0[arm_in],Tdep2[arm_in],marker='d',c='m',linestyle='None')
plt.loglog(sigma0[interarm_in],Tdep2[interarm_in],marker='o',c='g',linestyle='None')
plt.ylim(4E8, 1E12)
plt.xlim(5E-2,4)
plt.tight_layout() 	
plt.savefig('Tdep2_sigma0_matplotlib.png')




