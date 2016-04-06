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
H2SD = Tdep*SFRSD/1e6
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
plt.xlabel('$\Sigma\ (M_{\odot}\ \mathrm{pc}^{-2})$')
plt.ylabel(r'$\tau _\mathrm{H2}\ (\mathrm{yrs})$')
plt.loglog(M_den[centre_in],Tdep[centre_in],marker='v',c='b',linestyle='None')
plt.loglog(M_den[arm_in],Tdep[arm_in],marker='d',c='m',linestyle='None')
plt.loglog(M_den[interarm_in],Tdep[interarm_in],marker='o',c='g',linestyle='None')
plt.ylim(1E9, 1.5E10)
plt.xlim(10,1E4)
plt.tight_layout() 	
plt.savefig('Tdep_Mden_matplotlib.png')

#Depletion time vs Sigma_0
figure = plt.figure(figsize=(4.5,4))
plt.ylabel(r'$\tau _\mathrm{H2}\ (\mathrm{yrs})$')
plt.xlabel('$\sigma_0\ (\mathrm{km\ s^{-1}})$')
plt.loglog(sigma0[centre_in],Tdep[centre_in],marker='v',c='b',linestyle='None')
plt.loglog(sigma0[arm_in],Tdep[arm_in],marker='d',c='m',linestyle='None')
plt.loglog(sigma0[interarm_in],Tdep[interarm_in],marker='o',c='g',linestyle='None')
plt.ylim(1E9, 1.5E10)
plt.xlim(5E-2,4)
plt.tight_layout() 	
plt.savefig('Tdep_sigma0_matplotlib.png')

#Depletion time vs Surface Density
figure = plt.figure(figsize=(4.5,4))
plt.xlabel('$\Sigma\ (M_{\odot}\ \mathrm{pc}^{-2})$')
plt.ylabel(r'$\tau _\mathrm{GMA}\ (\mathrm{yrs})$')
plt.loglog(M_den[centre_in],Tdep2[centre_in],marker='v',c='b',linestyle='None')
plt.loglog(M_den[arm_in],Tdep2[arm_in],marker='d',c='m',linestyle='None')
plt.loglog(M_den[interarm_in],Tdep2[interarm_in],marker='o',c='g',linestyle='None')
plt.ylim(4E8, 1E12)
plt.xlim(10,1E4)
plt.tight_layout() 	
plt.savefig('Tdep2_Mden_matplotlib.png')

#Depletion time vs Sigma_0
figure = plt.figure(figsize=(4.5,4))
plt.ylabel(r'$\tau _\mathrm{GMA}\ (\mathrm{yrs})$')
plt.xlabel('$\sigma_0\ (\mathrm{km\ s^{-1}})$')
plt.loglog(sigma0[centre_in],Tdep2[centre_in],marker='v',c='b',linestyle='None')
plt.loglog(sigma0[arm_in],Tdep2[arm_in],marker='d',c='m',linestyle='None')
plt.loglog(sigma0[interarm_in],Tdep2[interarm_in],marker='o',c='g',linestyle='None')
plt.ylim(4E8, 1E12)
plt.xlim(5E-2,4)
plt.tight_layout() 	
plt.savefig('Tdep2_sigma0_matplotlib.png')

##SFE
figure = plt.figure(figsize=(4.5,4))
plt.plot(np.log10(np.power(mytable['VRMS_EXTRAP_DECONV'][centre_in],-0.32)*(np.power(mytable['MASS_EXTRAP'][centre_in]/np.power(mytable['RADRMS_EXTRAP_DECONV'][centre_in],3),0.5))),np.log10(1/Tdep[centre_in]),marker='v',c='b',linestyle='None')
plt.plot(np.log10(np.power(mytable['VRMS_EXTRAP_DECONV'][arm_in],-0.32)*(np.power(mytable['MASS_EXTRAP'][arm_in]/np.power(mytable['RADRMS_EXTRAP_DECONV'][arm_in],3),0.5))),np.log10(1/Tdep[arm_in]),marker='d',c='m',linestyle='None')
plt.plot(np.log10(np.power(mytable['VRMS_EXTRAP_DECONV'][interarm_in],-0.32)*(np.power(mytable['MASS_EXTRAP'][interarm_in]/np.power(mytable['RADRMS_EXTRAP_DECONV'][interarm_in],3),0.5))),np.log10(1/Tdep[interarm_in]),marker='o',c='g',linestyle='None')
plt.ylabel(r'$log_{10}(\tau _\mathrm{H2}^{-1})\ (\mathrm{yrs}^{-1})$')
plt.xlabel(r'$log_{10}(SFE)$')
plt.tight_layout() 	
plt.savefig('SFE_matplotlib.png')

# Krumholz/McKee SFE plot

import astropy.constants as con
# The typical sound speed in molecular gas is 0.2 km/s
mach = mytable['VRMS'].data/0.2
volume = 4*np.pi/3*(mytable['RADRMS_EXTRAP_DECONV'].data*u.pc)**3
rho = mytable['MASS_EXTRAP'].data*u.Msun/volume
tff = ((3*np.pi/(32*con.G*rho))**0.5).to(u.Myr)
virial_parameter = mytable['VIRMASS_EXTRAP_DECONV']/mytable['MASS_EXTRAP']
SFE = 0.014*(mach/100)**(-0.32)/tff*(virial_parameter)**(-0.68)
plt.clf()
figure = plt.figure(figsize=(4.5,4))
ax = figure.add_subplot(111)
plt.loglog(1e6/Tdep[centre_in],SFE[centre_in],marker='v',c='b',linestyle='None')
plt.plot(1e6/Tdep[arm_in],SFE[arm_in],marker='d',c='m',linestyle='None')
plt.plot(1e6/Tdep[interarm_in],SFE[interarm_in],marker='o',c='g',linestyle='None')
plt.xlabel(r'$\tau_{\mathrm{H2}}\ (\mathrm{Myr}^{-1})$')
plt.ylabel(r'$\tau_{\mathrm{KM}}\ (\mathrm{Myr}^{-1})$')
plt.xlim([1e-5,1e-2])
plt.ylim([1e-4,1e-1])
plt.plot([1e-5,1e-1],[1e-5,1e-1],lw=3,alpha=0.5)
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig('SFE_vs_KM.png')

tb2 = {'Clouds': ['All','Nuclear','Arm','Interarm'] 
	,'mean_T_dep': [
		'{:3.2f}'.format(np.nanmean(np.log10(Tdep))),
		'{:3.2f}'.format(np.nanmean(np.log10(Tdep[centre_in]))),
		'{:3.2f}'.format(np.nanmean(np.log10(Tdep[arm_in]))),
		'{:3.2f}'.format(np.nanmean(np.log10(Tdep[interarm_in])))]
	,'stde_T_dep': [
		'{:3.2f}'.format(np.nanstd(np.log10(Tdep),ddof=1)/np.sqrt(len(Tdep))),
		'{:3.2f}'.format(np.nanstd(np.log10(Tdep[centre_in]),ddof=1)/np.sqrt(len(Tdep[centre_in]))),
		'{:3.2f}'.format(np.nanstd(np.log10(Tdep[arm_in]),ddof=1)/np.sqrt(len(Tdep[arm_in]))),
		'{:3.2f}'.format(np.nanstd(np.log10(Tdep[interarm_in]),ddof=1)/np.sqrt(len(Tdep[interarm_in])))]
	,'mean_T_dep2': [
		'{:3.2f}'.format(np.nanmean(np.log10(Tdep2))),
		'{:3.2f}'.format(np.nanmean(np.log10(Tdep2[centre_in]))),
		'{:3.2f}'.format(np.nanmean(np.log10(Tdep2[arm_in]))),
		'{:3.2f}'.format(np.nanmean(np.log10(Tdep2[interarm_in])))]
	,'stde_T_dep2': [
		'{:3.2f}'.format(np.nanstd(np.log10(Tdep2),ddof=1)/np.sqrt(len(Tdep2))),
		'{:3.2f}'.format(np.nanstd(np.log10(Tdep2[centre_in]),ddof=1)/np.sqrt(len(Tdep2[centre_in]))),
		'{:3.2f}'.format(np.nanstd(np.log10(Tdep2[arm_in]),ddof=1)/np.sqrt(len(Tdep2[arm_in]))),
		'{:3.2f}'.format(np.nanstd(np.log10(Tdep2[interarm_in]),ddof=1)/np.sqrt(len(Tdep2[interarm_in])))]
	}


t2 = Table(tb2,names =('Clouds','mean_T_dep','stde_T_dep','mean_T_dep2','stde_T_dep2'))
print t2

t2.write('means_T.tex', format = 'latex')

















