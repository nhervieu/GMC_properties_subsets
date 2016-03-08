#import libraries
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from galaxies import Galaxy
import astropy.units as u
import powerlaw

#load file
mytable = Table.read('m100.co10.kkms_props_cprops_subsets_improved.fits')

##VARIABLE DEFINITIONS##

# Pull out the mass variable into a numpy array.
mass = mytable['MASS_EXTRAP'].data 
#remove two outliers    
mass = np.delete(mass,132)
mass = np.delete(mass,92)
 
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

mygalaxy = Galaxy("M100")
cpropstable = Table.read('m100.co10.kkms_props_cprops_subsets_improved.fits')
x, y = mygalaxy.radius(ra = mytable['XPOS'], dec = mytable['YPOS'], returnXY = True)
rgal=mygalaxy.radius(ra = cpropstable['XPOS'], dec = cpropstable['YPOS'])


#indexes for low R_gal group 
index = np.where(rgal.value < 1100)	
	
#centre
centre_in = np.where(mytable['Nuclear'] == True)
#arms
arm_in = np.where(mytable['Arms'] == True)
#print(arm_in)
#interarm	
interarm_in = np.where(mytable['Interarm'] == True)
#print(interarm_in)	
outliers = np.where(mytable['Outliers'] == True)
##PLOTS##

#Virial mass vs.luminous mass plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
plt.plot(one,one,linestyle = '-', c = 'k')
plt.loglog(mytable['MASS_EXTRAP'][arm_in],mytable['VIRMASS_EXTRAP_DECONV'][arm_in],marker='d',c='m',linestyle='None')
plt.loglog(mytable['MASS_EXTRAP'][interarm_in],mytable['VIRMASS_EXTRAP_DECONV'][interarm_in],marker='o',c='g',linestyle='None')
plt.loglog(mytable['MASS_EXTRAP'][outliers],mytable['VIRMASS_EXTRAP_DECONV'][outliers],marker='+',c='w',linestyle='None')
plt.loglog(mytable['MASS_EXTRAP'][centre_in],mytable['VIRMASS_EXTRAP_DECONV'][centre_in],marker='v',c='b',linestyle='None')

plt.xlabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$') 
plt.ylabel(r'$M_{\mathrm{vir}}\ (M_{\odot})$')
plt.ylim(4E5, 3E9)
plt.xlim(1E6,1E9)
#plt.axes().set_aspect('equal', 'datalim')
plt.tight_layout() 	
plt.savefig('MlumMvir_matplotlib.png')

#Radius vs. line width(velocity dispersion) plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
plt.plot(R,sigmav,linestyle = '-', c = 'k')
plt.ylabel(r'$\sigma\ (km\ s^{-1})$') 
plt.xlabel(r'$R\ (pc)$')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][arm_in],mytable['VRMS_EXTRAP_DECONV'][arm_in],marker='d',c='m',linestyle='None')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][interarm_in],mytable['VRMS_EXTRAP_DECONV'][interarm_in],marker='o',c='g',linestyle='None')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][centre_in],mytable['VRMS_EXTRAP_DECONV'][centre_in],marker='v',c='b',linestyle='None')
#plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][outliers],mytable['VRMS_EXTRAP_DECONV'][outliers],marker='+',c='m',linestyle='None')
plt.ylim(1, 1E2)
plt.xlim(35,1E3)
plt.tight_layout() 	
plt.savefig('LwRad_matplotlib.png')

#Luminous mass vs. radius plot
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
plt.plot(R,M_lum,linestyle = '-', c = 'k')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][arm_in],mytable['MASS_EXTRAP'][arm_in],marker='d',c='m',linestyle='None')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][interarm_in],mytable['MASS_EXTRAP'][interarm_in],marker='o',c='g',linestyle='None')
plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][centre_in],mytable['MASS_EXTRAP'][centre_in],marker='v',c='b',linestyle='None')
#plt.loglog(mytable['RADRMS_EXTRAP_DECONV'][outliers],mytable['MASS_EXTRAP'][outliers],marker='+',c='m',linestyle='None')
plt.xlabel(r'$R\ (pc)$') 
plt.ylabel(r'$M_{\mathrm{lum}}\ (M_{\odot})$')
plt.ylim(1E6, 1E9)
plt.xlim(40,1E3)
plt.tight_layout() 	
plt.savefig('MlumRad_matplotlib.png')

#Sigma_0 vs Mass Density
figure = plt.figure(figsize=(4.5,4))
plt.xlabel('$M/\pi R^2\ ((M_{\odot})/pc^2)$')
plt.ylabel('$\sigma_0$')
plt.loglog(M_den[arm_in],sigma0[arm_in],marker='d',c='m',linestyle='None')
plt.loglog(M_den[interarm_in],sigma0[interarm_in],marker='o',c='g',linestyle='None')
plt.loglog(M_den[centre_in],sigma0[centre_in],marker='v',c='b',linestyle='None')
plt.loglog(M_den[outliers],sigma0[outliers],marker='+',c='m',linestyle='None')
plt.ylim(3E-2, 6)
plt.xlim(1E1,1E4)
plt.tight_layout() 	
plt.savefig('Sigma0_Mden_matplotlib.png')

#Sigma_0 vs Galactocentric Radii
figure = plt.figure(figsize=(4.5,4))
plt.xlabel('$R_{gal} (pc)$')
plt.ylabel('$\sigma_0$')
plt.loglog(rgal.to(u.pc)[arm_in],sigma0[arm_in],marker='d',c='m',linestyle='None')
plt.loglog(rgal.to(u.pc)[interarm_in],sigma0[interarm_in],marker='o',c='g',linestyle='None')
plt.loglog(rgal.to(u.pc)[outliers],sigma0[outliers],marker='+',c='w',linestyle='None')
plt.loglog(rgal.to(u.pc)[centre_in],sigma0[centre_in],marker='v',c='b',linestyle='None')
plt.ylim(0.05, 5)
plt.xlim(1E2,5E4)
plt.axes().set_aspect('equal', 'datalim')
plt.tight_layout() 	
plt.savefig('sigma0_Rgal_matplotlib.png')

#X,Y position
figure = plt.figure(figsize=(4.5,4)) #figure size in inches
plt.xlabel('X position (kpc)') 
plt.ylabel('Y position (kpc)')
plt.plot(x.to(u.kpc)[centre_in],y.to(u.kpc)[centre_in],marker='v',c='b',linestyle='None')
plt.plot(x.to(u.kpc)[arm_in],y.to(u.kpc)[arm_in],marker='d',c='m',linestyle='None')
plt.plot(x.to(u.kpc)[interarm_in],y.to(u.kpc)[interarm_in],marker='o',c='g',linestyle='None')
plt.plot(x.to(u.kpc)[outliers],y.to(u.kpc)[outliers],marker='+',c='w',linestyle='None')
plt.tight_layout() 	
plt.axes().set_aspect('equal', 'datalim')
plt.savefig('xypos_matplotlib.png')



## MASS DISTRIBUTIONS ##

#Mass Distribution for All Clouds
figure = plt.figure(figsize=(4.5,4)) 
#Fit data 
myfit = powerlaw.Fit(mass,xmin = 1E7)
R, p = myfit.distribution_compare('power_law','truncated_power_law')
#Plot
myfit.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
myfit.power_law.plot_ccdf(label='Power Law')
myfit.plot_ccdf(drawstyle='steps',label='Data')
plt.legend(loc ='lower left',prop={'size':8})
plt.ylabel(r'$N$')
plt.xlabel(r'$Mass\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('powerlaw.png')



#Mass Distribution for Nuclear Clouds
figure = plt.figure(figsize=(4.5,4)) 
#Keep only clouds with rgal <1kpc
mass_nuc = mass[index]
myfit_nuc = powerlaw.Fit(mass_nuc, xmin = 1E7)
R_nuc, p_nuc = myfit_nuc.distribution_compare('power_law','truncated_power_law')
#Plot
myfit_nuc.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
myfit_nuc.power_law.plot_ccdf(label='Power Law')
myfit_nuc.plot_ccdf(drawstyle='steps',label='Data (Nuclear Clouds)')
plt.legend(loc ='lower left',prop={'size':8})
plt.ylabel(r'$N$')
plt.xlabel(r'$Mass\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('powerlaw_nuc.png')


#Mass Distribution for Disk Clouds
figure = plt.figure(figsize=(4.5,4)) 
#Ignore all clouds with rgal >1kpc
mass_disk = np.delete(mass,index)
#Fit Data
myfit_disk = powerlaw.Fit(mass_disk,xmin = 1E7)
R_disk, p_disk = myfit_disk.distribution_compare('power_law','truncated_power_law')
#Plot
myfit_disk.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
myfit_disk.power_law.plot_ccdf(label='Power Law')
myfit_disk.plot_ccdf(drawstyle='steps',label='Data (Disk Clouds)')
plt.legend(loc ='lower left',prop={'size':8})
plt.ylabel(r'$N$')
plt.xlabel(r'$Mass\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('powerlaw_disk.png')

#Mass Distribution for arm clouds
figure = plt.figure(figsize=(4.5,4)) 
mass_arm = mass[arm_in]
#Fit Data
myfit_arm = powerlaw.Fit(mass_arm,xmin = 1E7)
R_arm, p_arm = myfit_arm.distribution_compare('power_law','truncated_power_law')
#Plot
myfit_arm.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
myfit_arm.power_law.plot_ccdf(label='Power Law')
myfit_arm.plot_ccdf(drawstyle='steps',label='Data (Arm Clouds)')
plt.legend(loc ='lower left',prop={'size':8})
plt.ylabel(r'$N$')
plt.xlabel(r'$Mass\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('powerlaw_arm.png')

#Mass Distribution for inter arm clouds
figure = plt.figure(figsize=(4.5,4)) 
test = np.hstack([index,arm_in])
mass_interarm = np.delete(mass,test)
#mass_interarm = np.delete(mass_interarm,index)
#Fit Data
myfit_interarm = powerlaw.Fit(mass_interarm,xmin = 1E7)
R_interarm, p_interarm = myfit_interarm.distribution_compare('power_law','truncated_power_law')
#Plot
myfit_interarm.truncated_power_law.plot_ccdf(label='Trunc. Power Law')
myfit_interarm.power_law.plot_ccdf(label='Power Law')
myfit_interarm.plot_ccdf(drawstyle='steps',label='Data (Interarm Clouds)')
plt.legend(loc ='lower left',prop={'size':8})
plt.ylabel(r'$N$')
plt.xlabel(r'$Mass\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('powerlaw_interarm.png')

#plot on same window
figure = plt.figure(figsize=(4.5,4)) 
myfit_arm.truncated_power_law.plot_ccdf(label='Trunc. Power Law (Arm)')
myfit_arm.power_law.plot_ccdf(label='Power Law (Arm)')
myfit_arm.plot_ccdf(drawstyle='steps',label='Data (Arm)')
myfit_interarm.truncated_power_law.plot_ccdf(label='Trunc. Power Law (Interarm)')
myfit_interarm.power_law.plot_ccdf(label='Power Law (Interarm)')
myfit_interarm.plot_ccdf(drawstyle='steps',label='Data (Interarm)')

plt.legend(loc ='lower left',prop={'size':7})
plt.ylabel(r'$N$')
plt.xlabel(r'$Mass\ (M_{\odot})$')
plt.tight_layout() 	
plt.savefig('powerlaw_arm_interarm.png')


#print out table of alpha, R and p values for each mass distribution
tb = {'Clouds': ['All','Nuclear','Disk','Arm','Interarm'],'alpha': [myfit.alpha,myfit_nuc.alpha,myfit_disk.alpha,myfit_arm.alpha,myfit_interarm.alpha],'$\alpha$ (Truncated)': [myfit.truncated_power_law.alpha,myfit_nuc.truncated_power_law.alpha,myfit_disk.truncated_power_law.alpha,myfit_arm.truncated_power_law.alpha,myfit_interarm.truncated_power_law.alpha],'R':[R,R_nuc,R_disk,R_arm,R_interarm], 'p':[p,p_nuc,p_disk,p_arm,p_interarm], 'p':[p,p_nuc,p_disk,p_arm,p_interarm], 'p':[p,p_nuc,p_disk,p_arm,p_interarm], 'Cutoff Mass':[1/myfit.truncated_power_law.parameter2,1/myfit_nuc.truncated_power_law.parameter2,1/myfit_disk.truncated_power_law.parameter2,1/myfit_arm.truncated_power_law.parameter2,1/myfit_interarm.truncated_power_law.parameter2]}
t = Table(tb,names =('Clouds','alpha','alpha (Truncated)','R','p','Cutoff Mass'))
print t

t.write('table.tex', format = 'latex')
