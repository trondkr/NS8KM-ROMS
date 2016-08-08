import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *

import grd
#import seawater as sw
import eos80 as sw

import heatcap as hc
import glob
import os
from datetime import datetime, timedelta

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2016, 2, 8)
__modified__ = datetime(2016, 2, 8)
__version__  = "1.5"
__status__   = "Development"

# Calculate Heat content in the North Sea
#
# Script used for Havforsknigsrapporten 2016
#
# This script reads in data from the entire model domain but then only uses data values that are within
# the defined area for the North Sea. Everything else is masked and not used in the caculations.
# Masked area:
# if (0 < lat[j,i] < 59.5) and (0 < lon[j,i] <20) and (M[j,i]==1):

def contourHeatCap(lon,lat,hc):
	print "running contour map"
	mymap = Basemap(llcrnrlon=-18.0,
                      llcrnrlat=46.0,
                      urcrnrlon=25.5,
                      urcrnrlat=67.5,
                      resolution='i',projection='tmerc',lon_0=0,lat_0=50,area_thresh=50.)

	mymap.drawmeridians(np.arange(lon.min(),lon.max(),10),labels=[0,0,0,1])
	mymap.drawparallels(np.arange(lat.min(),lat.max(),4),labels=[1,0,0,0])
	levels = np.arange(np.min(hc)-0.1,np.max(hc)+0.1,(np.max(hc) - np.min(hc))/20.)

	x,y = mymap(lon,lat)

	CS = mymap.contourf(x, y, hc, levels,
                  cmap=plt.cm.RdYlBu_r,
                  origin='lower')

	#mymap.drawcoastlines()
	#mymap.drawcountries()
	mymap.fillcontinents(color='grey')
	plt.colorbar(CS,orientation='vertical',extend='both', shrink=0.5)
	CS.axis='tight'
	plotfile='figures/heatcap_NS8KM.pdf'
#    plt.savefig(plotfile, dpi='200')
	plt.show()
	

def extractTimeFromNetCDF(cdf,mytime):
  times = cdf.variables["ocean_time"][:]
  refDate=datetime(1948, 1, 1, 0, 0, 0)
  currentDate=refDate + timedelta(seconds=times[mytime])
  difference = currentDate - refDate
  #print "Days since refdate %s - %s"%(refDate,difference.days)
  if currentDate.month<10:
    month='0%s'%(currentDate.month)
  else:
    month='%s'%(currentDate.month)
  if currentDate.day<10:
    day='0%s'%(currentDate.day)
  else:
    day='%s'%(currentDate.day)
    
  return '%s-%s-%s'%(currentDate.year,month,day)

def findNorthSeaIndices(grdROMS):

	M = grdROMS.mask_rho
	lat = grdROMS.lat_rho
	lon = grdROMS.lon_rho
	Iw = np.zeros(np.shape(grdROMS.z_w))
	Ir = np.zeros(np.shape(grdROMS.z_r))
	for j in xrange(len(lat[:,0])):
		for i in xrange(len(lat[0,:])):
			if (0 < lat[j,i] < 59.5) and (0 < lon[j,i] <20) and (M[j,i]==1):
				Ir[:,j,i]=1
				Iw[:,j,i]=1
				

	print "Found %s valid indices of a total of %s (%s percent)"%(np.sum(Ir[0,:,:]),np.sum(M),float(sum(Ir[0,:,:])*1.0)/float(np.sum(M))*100.)
	return Iw, Ir

#romsgridpath = "/work/users/trondk/NS8km/FORCING/GRID/nordsjoen_8km_grid_hmax20m_v4.nc"
romsgridpath = "/Users/trondkr/Projects/is4dvar/Grid/nordsjoen_8km_grid_hmax20m_v3.nc"
fid=Dataset(romsgridpath)

grdROMS = grd.grdClass(romsgridpath,'NS8KM','NS8KM',False,'ocean','NS8KM')

Iw, Ir = findNorthSeaIndices(grdROMS)
Ir = np.asarray(Ir)
Iw = np.asarray(Iw)


M = grdROMS.mask_rho
Mns = np.ma.masked_where(Ir[0,:,:]==0,M)
pm = grdROMS.pm
pn = grdROMS.pn

pmns = np.ma.masked_where(Ir[0,:,:]==0,pm)
pnns = np.ma.masked_where(Ir[0,:,:]==0,pn)

# Horizontal area of the grid cells [m**2]
Area = 1./pmns * 1./pnns
# Vertical extent off the grid cells  [m]
z_w = grdROMS.z_w
z_wns = np.ma.masked_where(Iw==0,z_w)
	
dz = z_w[1:, :, :] - z_w[:-1, :, :]
# Volume  [m**3]
Vol = Area * dz

basedir='/imr/vol6/NS8KM/STORAGE-ASSIMILATION2010-2013/'

mypattern="ns8km_avg_*.nc"
argument="%s%s"%(basedir,mypattern)
allFiles = glob.glob(argument)
allFiles.sort()
    
outfile= 'Heatcontent_NorthSea_MyOcean.txt'
if os.path.exists(outfile): os.remove(outfile)
out=open(outfile,'a')

print "argument %s"%(argument)
print "Sorting %s files found in datadirectory"%(len(allFiles))

for myFile in allFiles:
	fidIn = Dataset(myFile)

	# If deep areas are included, should compute
	# in-situ temperature from the potential temperature
	# used by ROMS 

	# Mass  [kg]

	S = (fidIn.variables['salt'][0, :, :, :])
	T = (fidIn.variables['temp'][0, :, :, :])
	Tns = np.ma.masked_where(Ir==0,T)
	Sns = np.ma.masked_where(Ir==0,S)
	print "Temp original %s after: %s"%(np.sum(T),np.ma.mean(T))

	print "Temp masked %s after: %s"%(np.sum(Tns),np.ma.mean(Tns))

	z_r = grdROMS.z_r
	z_rns = np.ma.masked_where(Ir==0,z_r)
	
	Mass = Vol * sw.dens(Sns,Tns, z_rns)

	Mass = M*Mass   # Zero sea water mass on land

	# Heat capacity, absolute [J/K]
	abs_heat_cap = Mass * hc.heatcap(Sns, Tns, z_rns)

	#kg * (J/Kg*K) J/K
	# Total heat capacity, [J/kg]
	# Suppose north_sea_mask = 1 in North Sea, zero outside
	#north_sea_mask = np.ones(np.shape(M))
	total_heat = np.ma.sum(abs_heat_cap[:,:,:]*(Tns[:,:,:]))
	dataline = "%s\t %s\t %s\n"%(total_heat, extractTimeFromNetCDF(fidIn,0), np.ma.mean(Tns))
	out.writelines(dataline)
	print dataline
	fidIn.close()
	#contourHeatCap(grdROMS.lon_rho,grdROMS.lat_rho,abs_heat_cap[34,:,:])
out.close()

