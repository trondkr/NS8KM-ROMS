import os
import datetime as datetime
from netCDF4 import Dataset
import numpy as np
import glob
import string
from matplotlib.pyplot import cm 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *

__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime.datetime(2015, 6, 9)
__modified__ = datetime.datetime(2015, 10, 28)
__version__  = "1.0"
__status__   = "Production"

def docs():
	"""
	This script is run on Hexagon being called from script getstatus.py
	"""

def setupSystem():
	print "Initializing (func:setupSystem)"
	
	myremotedir="/work/shared/imr/NS8KM/FORWARD/Run/"
	mypattern = "ns8km_avg_000[0123]*.nc"
	mystartYear=2010
	mylocaldir="/Users/trondkr/Projects/is4dvar/RemoteStatus/"
	mygridfile="/work/shared/imr/NS8KM/FORCING/GRID/nordsjoen_8km_grid_hmax20m_v4.nc"
	myglorys2v3gridfile="/work/shared/imr/NS8KM/FORCING/GLORYS2V3/ftp.myocean.mercator-ocean.fr/Core/GLOBAL_REANALYSIS_PHYS_001_009/dataset-global-reanalysis-phys-001-009-ran-fr-glorys2v3-monthly-t/GLORYS2V3_ORCA025_19940115_R20130808_gridT.nc"
	myglorys2v3pathT="/work/shared/imr/NS8KM/FORCING/GLORYS2V3/ftp.myocean.mercator-ocean.fr/Core/GLOBAL_REANALYSIS_PHYS_001_009/dataset-global-reanalysis-phys-001-009-ran-fr-glorys2v3-monthly-t/"
	myglorys2v3pathS="/work/shared/imr/NS8KM/FORCING/GLORYS2V3/ftp.myocean.mercator-ocean.fr/Core/GLOBAL_REANALYSIS_PHYS_001_009/dataset-global-reanalysis-phys-001-009-ran-fr-glorys2v3-monthly-s/"
	
	return myremotedir, mypattern, mylocaldir,mygridfile,myglorys2v3gridfile,myglorys2v3pathT,myglorys2v3pathS,mystartYear

def getLocations(mygridfile):
	maxPositions=4
	
	xpos=np.zeros(shape=(maxPositions,1)); ypos=np.zeros(shape=(maxPositions,1)); xynames=np.zeros(shape=(maxPositions,1))
	lonpos=np.zeros(shape=(maxPositions,1)); latpos=np.zeros(shape=(maxPositions,1))
	xpos=[90,130,150,185]
	ypos=[30,25,90,35]
	xynames=["Soerlige Nordsjoen","Tyskebukten","Nordlige Nordsjoen","Skagerrak"]

	cdf = Dataset(mygridfile)
	lats=cdf.variables["lat_rho"][:]
	lons=cdf.variables["lon_rho"][:]
	lonpos = [lons[ypos[i],xpos[i]] for i in xrange(len(xpos[:]))]
	latpos = [lats[ypos[i],xpos[i]] for i in xrange(len(xpos[:]))]
	
	print "\nSystem setup with station locations:"
	for name,x,y,lo,la in zip(xynames,xpos,ypos,lonpos,latpos):
		print "%s => longitude=%3.2f  latitude=%3.2f (x=%s,y=%s)"%(name,lo,la,x,y)

	return xpos,ypos,xynames,lonpos,latpos

def getIndicesOfLocationsInGLORYS2V3Grid(lonpos,latpos,xynames,myglorys2v3gridfile):
	
	print "Creating indices for GLORYS2V3 grid (func:getIndicesOfLocationsInGLORYS2V3Grid)"
	
	cdf = Dataset(myglorys2v3gridfile)
	latsglorys=cdf.variables["nav_lat"][:]
	lonsglorys=cdf.variables["nav_lon"][:]
	numberOfPoints=1
	glorysIndicies=[]
	glorys=np.zeros(shape=(np.shape(lonpos)[0],2))

	for st_lon, st_lat in zip(lonpos,latpos):
		indices, dis = getStationIndices(st_lon, st_lat,lonsglorys,latsglorys,numberOfPoints)
		glorysIndicies.append(indices)

	for i in xrange(len(glorysIndicies)):
		glorys[i,0]=int(glorysIndicies[i][0])
		glorys[i,1]=int(glorysIndicies[i][1])
	
	return glorysIndicies 

def getStationIndices(st_lon,st_lat,longitudes,latitudes,numberOfPoints):
    
    distance = np.zeros((longitudes.shape),dtype=np.float64)
    listd=[]
    for eta in range(len(latitudes[:,0])):
        for xi in range(len(latitudes[0,:])):
            distance[eta,xi] = np.sqrt( (latitudes[eta,xi]-st_lat)**2.0 + (longitudes[eta, xi] - st_lon)**2.0 )
            listd.append(distance[eta,xi])

    listsIndexes=[]
    listd.sort()
    for i in range(numberOfPoints):
        value=listd[0]
        itemindex=np.where(distance==value)
        listsIndexes.append(itemindex)
       
        listd.pop(0)
  
    print '' 
    print '=====getStationIndices======'
    print 'Looking for longitude [%3.3f] and latitude [%3.3f]'%(st_lon,st_lat)
    print 'Result ===>'
    for i in range(numberOfPoints):
        print 'Found index pair in gridfile',listsIndexes[i]
        print 'Index corresponds to longitude [%3.3f] and latitude [%3.3f]'%(longitudes[listsIndexes[i][0],listsIndexes[i][1]],latitudes[listsIndexes[i][0],listsIndexes[i][1]])
    print '======================'
    print ''
      
    dis=[]
    
    for i in range(numberOfPoints):
        dis.append(np.sqrt( (latitudes[listsIndexes[i][0],listsIndexes[i][1]]-st_lat)**2.0 + (longitudes[listsIndexes[i][0],listsIndexes[i][1]] - st_lon)**2.0 ))
        
    return listsIndexes[0], dis


def getValuesFromGLORYS(glorysIndicies,myglorys2v3pathT,myglorys2v3pathS,xynames,startYear,debug):
	
	print "Extracting values from GLORYS2V3 (func:getValuesFromGLORYS)"
	
	mypaths=[myglorys2v3pathT,myglorys2v3pathS]
	variables=["votemper","vosaline"]
	first=True
	counter=0; varCounter=0

	for variable,datapath in zip(variables,mypaths):

		argument="%s*.nc"%(datapath)
		allFiles = glob.glob(argument)
		allFiles.sort()
		allFilesFiltered=[]

		# Remove files older than startYear (we only want 2009-2012)
		for afile in allFiles:
			(head, filename) = os.path.split(afile)
			l=string.split(filename,'_')
			date = str(l[2])
			if (startYear <= int(date[0:4])):
				allFilesFiltered.append(afile)
		allFiles = []
		allFiles = allFilesFiltered
		allFiles.sort()

		if first:
			allValues=np.zeros(shape=(len(glorysIndicies),len(allFiles),len(variables)))
			allDates=[]
			print "Sorting %s files found in GLORYS2V3 datadirectory"%(len(allFiles))
			first=False

		refdate=datetime.datetime(1948,1,1)
		
		for afile in allFiles:
			(head, filename) = os.path.split(afile)
			l=string.split(filename,'_')
			date = str(l[2])
			mydate=datetime.datetime(int(date[0:4]),int(date[4:6]),int(date[6:8]))
			if (varCounter==0):
				allDates.append(mydate)

			cdf = Dataset(afile)

			for station in xrange(len(glorysIndicies)):
				data = cdf.variables[variable][0,0,int(glorysIndicies[station][0]),int(glorysIndicies[station][1])]
				if debug:
					print "%s,%i : %s  - data for station (%s,%s) %s = > %s"%(counter,varCounter,mydate,int(glorysIndicies[station][0]),int(glorysIndicies[station][1]),variable,data)
				allValues[station,counter,varCounter]=data
			cdf.close()
			counter+=1
		counter=0
		varCounter+=1
	return allValues, allDates


def getValuesFromNS8KM(xpos,ypos,xynames,datapath,mypattern,debug):
	
	print "Extracting values from NS8KM (func:getValuesFromNS8KM)"
	
	first=True
	counter=0
	variables=["temp","salt"]
	argument="%s%s"%(datapath,mypattern)
	allFiles = glob.glob(argument)
	allFiles.sort()
		
	print "Sorting %s files found in NS8KM datadirectory"%(len(allFiles))

	refdate=datetime.datetime(1948,1,1)
		
	for afile in allFiles:
		
		cdf = Dataset(afile)
		times = cdf.variables["ocean_time"][:]
		if first:
			allValues=np.zeros(shape=(len(xpos),len(allFiles)*len(times),len(variables)))
			allDates=[]
			first = False
	
		for t,time in enumerate(times):
		

			mydate=(refdate + datetime.timedelta(seconds=int(time)))
			allDates.append(mydate)

			for station in xrange(len(xpos)):
				dataSST = cdf.variables["temp"][t,34,int(ypos[station]),int(xpos[station])]
				dataSSS = cdf.variables["salt"][t,34,int(ypos[station]),int(xpos[station])]
				if debug:
					print "%s : temp  - data for station (%s,%s) %s = > %s"%(counter,mydate,int(xpos[station]),int(ypos[station]),dataSST)
					print "%s : salt  - data for station (%s,%s) %s = > %s"%(counter,mydate,int(xpos[station]),int(ypos[station]),dataSSS)
				allValues[station,counter,0]=dataSST
				allValues[station,counter,1]=dataSSS
			counter+=1
		cdf.close()
		

	return allValues, allDates

def setupSubPlot(subplotIndex):
	from matplotlib.dates import  DateFormatter
	ax=plt.subplot(2, 1, subplotIndex)
	ax.tick_params(axis='both', which='major', labelsize=6)
	ax.tick_params(axis='both', which='minor', labelsize=6)

  	plt.xticks(rotation=45)
  	ax.xaxis.set_major_formatter( DateFormatter('%Y-%m-%d') )
  	return ax

def createTimeseriesPlot(stations,dataNS8KM,datesNS8KM,dataGLORYS2V3,datesGLORYS2V3,startDate,endDate,plotfileName):
	
	print "Creating timeseries plot (func:createTimeseriesPlot)"
	dates1 = matplotlib.dates.date2num(datesNS8KM)
  	dates2 = matplotlib.dates.date2num(datesGLORYS2V3)
  	
	ax = setupSubPlot(subplotIndex=1)
	fmtNS=["r-o","g-o","b-o","m-o"]
	fmtGL=["r-.","g-.","b-.","m-."]
	plots=[]
  	for i,station in enumerate(stations):
  		# SST
		myplt,=plt.plot_date(dates1, dataNS8KM[i,:,0],fmt=fmtNS[i],markersize=3,xdate=True,ydate=False)
		plots.append(myplt)

	for i,station in enumerate(stations):
  		# SST
  		plt.plot_date(dates2, dataGLORYS2V3[i,:,0],fmt=fmtGL[i],linewidth=2, xdate=True,ydate=False)
  	
  	legend(plots,stations,loc=3,prop={'size':6})

	ax.set_xlim(matplotlib.dates.date2num(startDate),matplotlib.dates.date2num(endDate))

	#ax.set_xlim(matplotlib.dates.date2num(startDate),matplotlib.dates.date2num(datetime.datetime(2013,6,23,0,0)))

	plt.ylabel('SST')

	ax2 = setupSubPlot(subplotIndex=2)
  	for i,station in enumerate(stations):
  		# SSS
		plt.plot_date(dates1, dataNS8KM[i,:,1],fmt=fmtNS[i],markersize=3,xdate=True,ydate=False)

	for i,station in enumerate(stations):
  		# SSS
		plt.plot_date(dates2, dataGLORYS2V3[i,:,1],fmt=fmtGL[i],linewidth=2, xdate=True,ydate=False)

	plt.ylabel('SSS')
	startLim=matplotlib.dates.date2num(startDate)
	endLim=matplotlib.dates.date2num(endDate)
	print startLim, endLim
	ax2.set_xlim(startLim,endLim)

	if os.path.exists(plotfileName): os.remove(plotfileName)
	plt.savefig(plotfileName,dpi=300)
	print 'Saved figure file %s\n'%(plotfileName)
   
def main():

	print "Starting program (func:main)"
	debug = False

	myremotedir, mypattern, mylocaldir, mygridfile, myglorys2v3gridfile,myglorys2v3pathT,myglorys2v3pathS,mystartYear = setupSystem()
	xpos,ypos,xynames,lonpos,latpos = getLocations(mygridfile)
	glorysIndicies = getIndicesOfLocationsInGLORYS2V3Grid(lonpos,latpos,xynames,myglorys2v3gridfile)
	allGLORYS2V3Values, allGLORYS2V3Dates = getValuesFromGLORYS(glorysIndicies,myglorys2v3pathT,myglorys2v3pathS,xynames,mystartYear,debug)
 	allNS8KMValues, allNS8KMDates = getValuesFromNS8KM(xpos,ypos,xynames,myremotedir,mypattern,debug)

 	# Create plot for date range of GLORYS2V3 defined dates
 	createTimeseriesPlot(xynames,allNS8KMValues,allNS8KMDates,allGLORYS2V3Values,allGLORYS2V3Dates,allGLORYS2V3Dates[0],allGLORYS2V3Dates[-1],"timeseries_NS8KM_vs_GLORYS2V3_humidityOK.png")

 	# Create zoomed plot for date range of NS8KM defined dates
 	createTimeseriesPlot(xynames,allNS8KMValues,allNS8KMDates,allGLORYS2V3Values,allGLORYS2V3Dates,allNS8KMDates[0],allNS8KMDates[-1],"timeseries_NS8KM_vs_GLORYS2V3_zoomed_humidityOK.png")



if __name__ == "__main__":
	main()


