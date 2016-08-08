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
import csv
import datetime

__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime.datetime(2016, 2, 9)
__modified__ = datetime.datetime(2016, 2, 9)
__version__  = "1.0"
__status__   = "Production"

def docs():
	"""
	This program plots the output of running calculateHeatContent.py
	"""

def readData(infileName,mytype):
	
	counter=0
	for row in csv.reader(open(infileName),delimiter='\t'):
		if (float(row[0])>0.6e21 and mytype=='Nowmaps'):
			counter+=1
		if (mytype in ['Forward_v35','Forward_v37']):
			counter+=1
	print "Reading in %s rows of data"%(counter)
	mydata=np.zeros((3,counter))
	mydates=[]
	counter=0
	for row in csv.reader(open(infileName),delimiter='\t'):
		
		if (float(row[0])>0.6e21 and mytype=='Nowmaps'):
			mydates.append(datetime.datetime.strptime(row[1], ' %Y-%m-%d')) # Parse date.
			mydata[1,counter] = float(row[0]) 
			mydata[2,counter] = float(row[2]) 
			counter+=1
		if (mytype in ['Forward_v35','Forward_v37']):
			mydates.append(datetime.datetime.strptime(row[1], ' %Y-%m-%d')) # Parse date.
			mydata[1,counter] = float(row[0]) 
			mydata[2,counter] = float(row[2]) 
			counter+=1
	return mydata,mydates

def setupSubPlot(subplotIndex):
	
	from matplotlib.dates import  DateFormatter
	ax=plt.subplot(2, 1, subplotIndex)
	
	ax.tick_params(axis='both', which='major', labelsize=9)
	ax.tick_params(axis='both', which='minor', labelsize=9)

  	plt.xticks(rotation=-30)
  	ax.xaxis.set_major_formatter( DateFormatter('%Y-%m-%d') )
  	return ax

def createTimeseriesPlot(dataCopernicus,datesCopernicus,dataForward,datesForward,dataForward35,datesForward35,dataForward37,datesForward37,
 		dataForward35KNI,datesForward35KNI,dataForward37HUM,datesForward37HUM,plotfileName,showRutgers):
	
	print "Creating timeseries plot (func:createTimeseriesPlot)"
	fig = plt.figure(figsize=(8, 6)) 
	ax = setupSubPlot(subplotIndex=1)
	fmtNS=["r-","m-","b-","c-","k-","y-"]
	plots=[]
	myWidth=2
  	myplt,=plt.plot_date(datesCopernicus, dataCopernicus[1,:],fmt=fmtNS[0],linewidth=myWidth,xdate=True,ydate=False)
	plots.append(myplt)

  	myplt,=plt.plot_date(datesForward, dataForward[1,:],fmt=fmtNS[1],linewidth=myWidth, xdate=True,ydate=False)  	
  	plots.append(myplt)

  	if showRutgers is True:
		myplt,=plt.plot_date(datesForward35, dataForward35[1,:],fmt=fmtNS[2],linewidth=myWidth, xdate=True,ydate=False)  	
  		plots.append(myplt)

	myplt,=plt.plot_date(datesForward35KNI, dataForward35KNI[1,:],fmt=fmtNS[3],linewidth=myWidth, xdate=True,ydate=False)  	
  	plots.append(myplt)

  	if showRutgers is True:
		myplt,=plt.plot_date(datesForward37, dataForward37[1,:],fmt=fmtNS[4],linewidth=myWidth, xdate=True,ydate=False)  	
  		plots.append(myplt)

  	myplt,=plt.plot_date(datesForward37HUM, dataForward37HUM[1,:],fmt=fmtNS[5],linewidth=myWidth, xdate=True,ydate=False)  	
  	plots.append(myplt)

  	mytypes=['Assimilasjon (v37)','Ismodell (Jon) v35','Rutgers (v35)', 'Ismodell uten is (v35)','Rutgers (v37)']
  	if showRutgers is True:
  		mytypes=['Assimilation (v37)','Icemodel v35','Rutgers (v35)', 'Icemodel no ice (v35)','Rutgers (v37)','Rutgers (v37) humidity']
  	else:
  		mytypes=['Assimilation (v37)','Icemodel v35', 'Icemodel no ice (v35)','Rutgers (v37) humidity']
  	
  	legend(plots,mytypes,loc=4,prop={'size':8})
  	ax.set_xlim(matplotlib.dates.date2num(datesCopernicus[0]),matplotlib.dates.date2num(datesCopernicus[-1]))
	ax.set_xticklabels([])
	plt.ylabel('Heat (Joule)')

	ax2 = setupSubPlot(subplotIndex=2)
  	plt.plot_date(datesCopernicus, dataCopernicus[2,:],fmt=fmtNS[0],linewidth=myWidth,xdate=True,ydate=False)
	plt.plot_date(datesForward, dataForward[2,:],fmt=fmtNS[1],linewidth=myWidth, xdate=True,ydate=False)
  	if showRutgers is True:
  		plt.plot_date(datesForward35, dataForward35[2,:],fmt=fmtNS[2],linewidth=myWidth, xdate=True,ydate=False)
  	plt.plot_date(datesForward35KNI, dataForward35KNI[2,:],fmt=fmtNS[3],linewidth=myWidth, xdate=True,ydate=False)
  	if showRutgers is True:
  		plt.plot_date(datesForward37, dataForward37[2,:],fmt=fmtNS[4],linewidth=myWidth, xdate=True,ydate=False)
  	
	plt.ylabel('Average temperature ($^\circ$C)')
	startLim=matplotlib.dates.date2num(datesCopernicus[0])
	endLim=matplotlib.dates.date2num(datesCopernicus[-1])
	
	ax2.set_xlim(startLim,endLim)

	if os.path.exists(plotfileName): os.remove(plotfileName)
	plt.savefig(plotfileName, dpi=300)
	print 'Saved figure file %s\n'%(plotfileName)
   	#plt.show()

def main():

	print "Starting program (func:main)"
	debug = True
	showRutgers=False
	if showRutgers is True:
		plotfilename="timeseries_Copernicus_heatcontent_wRutgers.jpeg"
	else:
		plotfilename="timeseries_Copernicus_heatcontent.jpeg"

	infileCopernicus='Heatcontent_NorthSea_MyOcean_2010_2013.txt'
	infileForward='Heatcontent_NorthSea_MyOcean_Jon.txt'
	infileForward35='Heatcontent_NorthSea_MyOcean_FORWARD35.txt'
	infileForward37='Heatcontent_NorthSea_MyOcean_FORWARD37.txt'
	infileForward37HUM='Heatcontent_NorthSea_MyOcean_FORWARD35_HUMIDITY.txt'
	infileForward35KNI='Heatcontent_NorthSea_MyOcean_FORWARD35-KATENOICE.txt'

	infileForwardKATENOICE='Heatcontent_NorthSea_MyOcean_FORWARD35-KATENOICE.txt'

	dataCopernicus,datesCopernicus = readData(infileCopernicus,'Nowmaps')
	dataForward,datesForward = readData(infileForward,'Forward_v37')
	dataForward35,datesForward35 = readData(infileForward35,'Forward_v35')
	dataForward37,datesForward37 = readData(infileForward37,'Forward_v37')
	dataForward35KNI,datesForward35KNI = readData(infileForward35KNI,'Forward_v35')
	dataForward37HUM,datesForward37HUM = readData(infileForward37HUM,'Forward_v37')

 	# Create plot for date range of GLORYS2V3 defined dates
 	createTimeseriesPlot(dataCopernicus,datesCopernicus,
 		dataForward,datesForward,
 		dataForward35,datesForward35,
 		dataForward37,datesForward37,
 		dataForward35KNI,datesForward35KNI,
 		dataForward37HUM,datesForward37HUM,
 		plotfilename,showRutgers)

 	


if __name__ == "__main__":
	main()


