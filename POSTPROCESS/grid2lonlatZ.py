# -*- coding: utf-8 -*-

# This program creates the final delivery files for NOWMAPS Copernicus reanalysis. The script takes files on the ROMS grid at Z_LEVELS
# as input and interpolates to a rectangulatr grid with the correct NetCDF attributes added according to MyOcean standard.
#
# Trond Kristiansen, 9.2.2016, 22.08.2015, 2014


import numpy as np
#from mathgrid import *
import os
from os import listdir
from os.path import isfile, join
from netCDF4 import Dataset, MFDataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
import string, copy
try:
    import ESMF
except ImportError:
    print "Could not find module ESMF"
    pass


def findWeights(inlevels, outlevels):
    weights={}
    debug=False

    for i,outlevel in enumerate(outlevels):
        if outlevel in inlevels:
            weights[outlevel]=[np.abs(np.subtract.outer(inlevels, outlevel)).argmin(),1.0]
            print "Found direct depth level %s %s"%(outlevel,weights[outlevel])
        else:
            inlevelscopy=inlevels[:]
            index1=np.abs(np.subtract.outer(inlevelscopy, outlevel)).argmin()
            index1=np.abs(np.subtract.outer(inlevels, inlevelscopy[index1])).argmin()
            if debug:
                print "Want depth: %s"%(outlevel)
                print "=> index1 %s value %s"%(index1,inlevels[index1])

            inlevelscopy.pop(index1)
            index2=np.abs(np.subtract.outer(inlevelscopy, outlevel)).argmin()
            index2=np.abs(np.subtract.outer(inlevels, inlevelscopy[index2])).argmin()

            if debug:
                print "=> index2 %s value %s"%(index2,inlevels[index2])
                print abs(outlevel-inlevels[index2]), float(abs(inlevels[index2]-inlevels[index1]))
                print abs(outlevel-inlevels[index1]), float(abs(inlevels[index2]-inlevels[index1]))

            weights[outlevel]=[index1,index2,abs(outlevel-inlevels[index2])/float(abs(inlevels[index2]-inlevels[index1])),
            abs(outlevel-inlevels[index1])/float(abs(inlevels[index2]-inlevels[index1]))]
            if debug:
                print "Found in-between depth levels %s %s"%(outlevel,weights[outlevel])
    return weights

def progressbar(message, percent, fileCounter):
    # http://stackoverflow.com/questions/3002085/python-to-print-out-status-bar-and-percentage 
    sys.stdout.write('\r')
    message = "%d%% (var: %s file number: %s)"%(percent,message,fileCounter)
    sys.stdout.write("  -> [%-100s] %s" % ('='*int(percent), message))
    sys.stdout.flush()


def getReferenceGrid():
    # Det rektangulære griddet NOWMAPS som data skal interpoleres til kommer
    # fra MetOffice (email from Inga Golbeck)
    # http://marine.copernicus.eu/web/69-interactive-catalogue.php?option=com_csw&view=details&product_id=NORTHWESTSHELF_REANALYSIS_PHYS_004_009

    gridfile="Reference_grid_edited.txt"
    fid=open(gridfile)

    lons=[]
    lats=[]
    lines=fid.readlines()
    fid.close()

    finishedLons=False
    for line in lines:
        line=line.rstrip()
        print line
        if (len(line)> 10 and finishedLons is False):

            l=string.split(line,',')

            for lon in l:
                try:
                    lon=lon.replace(" ", "")
                    lon=lon.replace(";", "")
                    print "lon",lon
                    lon=float(lon)
                    lons.append(lon)
                    print "=>",lon
                except:
                    print "Not a float=>", lon
        else:
            finishedLons=True
            l=string.split(line,',')
         
            for lat in l:
                try:
                    lat=lat.replace(" ", "")
                    lat=lat.replace(";", "")
                    
                    lat=float(lat)
                    lats.append(lat)
                    print "=>",lat
                except:
                    print "Not a float=>", lat

    return lons, lats
# ---------------------------------------------
#                   MAIN
# ---------------------------------------------

hexagon=False

if hexagon:
    mypath="/work/shared/imr/NS8KM/Z-LEVEL-ASSIMILATION2010-2013"
    outputdirname="/work/shared/imr/NS8KM/COPERNICUS_DELIVERY_30052016_REGRIDDED"
else:
    mypath="/Users/trondkr/Dropbox/delivery-COPERNICUS/salinity"
    outputdirname="/Users/trondkr/Projects/NOWMAPS/grid2lonlat/salt"

infiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
infiles.sort()

# -----------

UNDEF = -32768     # Short integer UNDEF
FUNDEF = -1.0E34   # Floating point UNDEF
FIRST = True

# For special cases we dont interpolate from ROMS grid but from the UK grid to UK grid. If 
# so the input and output grid are equal.

inputgridEqualToOutputgrid=True

# Dette kunne tas automatisk fra input fil
# På den annen side kan nå variable lett kuttes ut


# 2D variable, YX
varY = []  # konverterer som test på interpolasjonen

# 2D variable, YX
#varXY = ['lon', 'lat']  # konverterer som test på interpolasjonen
varXY=[]
# 3D variable, TYX
#varTYX =['Elev']
varTYX=[]

# 4D variable
varTZYX = ['vosaline']

# Levels
inlevels = [0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500, 600,
          700, 800, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]

outlevels = [0, 3, 10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600,
          750, 1000, 1500, 2000, 3000, 4000, 5000]

# Calculate the weights for vertical interpolation
weights=findWeights(inlevels,outlevels)

print "\nStarted converting netcdf files in polar stereographic projection to"
print "rectangular grid. Input files need to be on Z-level. Scrip must be run in"
print "directory where inputfiles are stored. Results will be stored in sub-diretory: RESULTS\n"

if not os.path.exists(outputdirname):
        os.makedirs(outputdirname)

for infile in infiles:
    print infile
    for myvar in varTZYX:
        # Read data
        outfilename=""
        infilename=""
        outfilename = outputdirname+"/"+infile[:-18]+str(myvar)+"-"+"20160621.nc"
        infilename = mypath+"/"+infile
        f0 = Dataset(infilename)
        print "Opened input file: %s"%(infilename)
        print "Results will be stored to: %s"%(outfilename)

        # Create lookup-table on first loop of infiles
        if (FIRST):
            temp0 = np.array(f0.variables[varTZYX[0]][0,0,:,:])
            lonin=f0.variables["longitude"][:]
            latin=f0.variables["latitude"][:]

            sea_mask = where(temp0 == UNDEF, 0, 1)

            # Lager her en look-up tabell fra lengde, bredde til grid-koord. Dette oppsettet
            # forutsetter at du kjenner grid koordinatene i lengde og breddegrad
            # 1. Definer output grid rektangulært
            lon,lat = getReferenceGrid()
            lon2=np.zeros(np.shape(lon))
            lat2=np.zeros(np.shape(lat))
            
            # Rydd opp slik at matrisene ser bra ut uten for mange desimaler
            for l in xrange(len(lon)):
                lon2[l]="%3.3f"%(lon[l])
            for l in xrange(len(lat)):
                lat2[l]="%3.3f"%(lat[l])

            lon=lon2; lat=lat2
            llon, llat = meshgrid(lon, lat)

            print "Turning on debugging for ESMF"
            if os.path.exists("grid2lonlat.ESMF_LogFile"): os.remove("grid2lonlat.ESMF_LogFile")
            ESMF.Manager(logkind=ESMF.LogKind.MULTI, debug=True)
            
            # Destination grid setup
            max_index = np.array([lat.size, lon.size])
            dstgrid = ESMF.Grid(max_index, coord_sys=ESMF.CoordSys.SPH_DEG,staggerloc=[ESMF.StaggerLoc.CENTER])

            gridCoordLon = dstgrid.get_coords(0)
            gridCoordLat = dstgrid.get_coords(1)
            gridCoordLon[...] = llon
            gridCoordLat[...] = llat  
        
            rectlons=dstgrid.coords[0][0][:]
            rectlats=dstgrid.coords[0][1][:]

            print "\nDESTINATION GRID properties"
            print "Created grid properties using ESMF"
            print "lons: ", np.min(rectlons), np.max(rectlons)
            print "lats", np.min(rectlats), np.max(rectlats)
            print "shapes", np.shape(rectlons), np.shape(rectlats)
            print "--------------------------------------------\n"

            # Source grid setup
            if inputgridEqualToOutputgrid:
                max_index = np.array([lat.size, lon.size])
                romsgrid = ESMF.Grid(max_index, coord_sys=ESMF.CoordSys.SPH_DEG,staggerloc=[ESMF.StaggerLoc.CENTER])
                gridCoordLon = romsgrid.get_coords(0)
                gridCoordLat = romsgrid.get_coords(1)
                gridCoordLon[...] = llon
                gridCoordLat[...] = llat  
            else:
                romsgrid = ESMF.Grid(filename=infilename, filetype=ESMF.FileFormat.GRIDSPEC,
                                          is_sphere=True, coord_names=['longitude','latitude'], add_mask=False)
        
            rectlons=romsgrid.coords[0][0][:]
            rectlats=romsgrid.coords[0][1][:]
            print "\nSOURCE GRID properties"
            print "Created grid properties using ROMS file"
            print "lons: ", np.min(rectlons), np.max(rectlons)
            print "lats", np.min(rectlats), np.max(rectlats)
            print "shapes", np.shape(rectlons), np.shape(rectlats)
            print "--------------------------------------------\n"

            print "  -> regridSrc2Dst - creating lookup table"
            fieldSrc = ESMF.Field(romsgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER,mask_values=sea_mask)
            fieldDst = ESMF.Field(dstgrid, "fieldDst", staggerloc=ESMF.StaggerLoc.CENTER)
            regridSrc2Dst = ESMF.Regrid(fieldSrc, fieldDst, regrid_method=ESMF.RegridMethod.BILINEAR,unmapped_action=ESMF.UnmappedAction.IGNORE)
            
            print "  -> Running interpolations"
            FIRST = False


        # --------------------------
        # Create output NetCDF file
        # --------------------------
        if os.path.exists(outfilename): os.remove(outfilename)
        f1 = Dataset(outfilename, mode='w', format='NETCDF4_CLASSIC')

        # Global attributes
        f1.Conventions = "CF-1.6"
        f1.title = "Monthly-mean (full water column) fields"
        f1.history = "Simulations were done using a 8x8 km polar stereographic grid projection, however the final " \
                     "data are presented using a reference grid. Conversion between grid " \
                     "projectionsby grid2lonlatZ.py"
        f1.source = "IMR, ROMSv3.7, IS4DVAR, NorthSea-8km reanalysis"
        f1.institution = "Institute of Marine Research, Norway"
        f1.references = "http://www.imr.no"
        f1.product_version = "1.0"
        f1.contact = "Trond.Kristiansen@imr.no"
        f1.netcdf_version_id = "netCDF-4 Classic"
        # Define dimensions
        f1.createDimension('time', None)
        f1.createDimension('depth', len(outlevels))
        f1.createDimension('longitude', len(lon))
        f1.createDimension('latitude', len(lat))

        v = f1.createVariable('time', 'd', ('time',))
        v0 = f0.variables['time']
        v.long_name = 'time'
        v.units = "Days since 1948-01-01 00:00:00"
        v.calendar = "Gregorian"

        ntimes = len(f0.dimensions['time'])
        v[:ntimes] = v0[:]

        v = f1.createVariable('depth', 'd', ('depth',))
        v.long_name = "depth"
        v.units = 'm'
        v.positive = 'down'
        v[:] = outlevels

        # Define coordinate variables
        v = f1.createVariable('longitude', 'd', ('longitude',))
        v.long_name = 'longitude'
        v.units = 'degrees_east'
        v[:] = lon

        v = f1.createVariable('latitude', 'd', ('latitude',))
        v.long_name = 'latitude'
        v.units = 'degrees_north'
        v[:] = lat

        # Data variable
        for var in varXY:
            print "Creating (x,y) variable", var
            v0 = f0.variables[var]
            v1 = f1.createVariable(var, 's', ('latitude', 'longitude'), fill_value=UNDEF)
            try:
                v1.long_name = v0.long_name
                v1.units = v0.units
            except:
                pass

            v1.scale_factor = v0.scale_factor
            v1.add_offset = v0.add_offset
            v1.missing_value = UNDEF
          #  v1._CoordinateAxes = "lat lon "

            # Modification
            if var == 'Topo':
                v1.standard_name = 'sea_floor_depth'


            F0 = np.array(v0[:,:])
            #undef = (F0==-32767)
            F1 = sample2D(F0, X, Y)
            # hadde trengt et valid-range atributt, U og V er ofte < null
            F1[F1<=-10.0] = FUNDEF              # ikke helt bra, temperature kan være < 0


            v1[:,:] = F1

        for var in varY:
            print var
            v0 = f0.variables[var]
            v1 = f1.createVariable(var, 'f', ('depth'),
                                   fill_value=UNDEF)
            try:
                v1.long_name = v0.long_name
                v1.units = v0.units
            except:
                pass

        for var in varTYX:
            print var
            v0 = f0.variables[var]
            v1 = f1.createVariable(var, 'i2', ('time', 'latitude', 'longitude'),
                                   fill_value=UNDEF)

            if var == 'ssflux':
                v1.standard_name = "surface_salinity_flux"
                v1.field="ssflux, scalar, series"
            elif var == 'shflux':
                v1.standard_name = "surface_downward_heat_flux_in_sea_water"
                v1.field="shflux, scalar, series"
            elif var=="zeta":
                 v1.standard_name = "sea_surface_height_above_geoid"
                 v1.field="zeta, scalar, series"
            else:
                v1.standard_name=v0.standard_name

            v1.long_name = v0.long_name
            v1.units = v0.units
            v1.scale_factor=v0.scale_factor
            v1.add_offset=v0.add_offset
           # v1.valid_min=v0.valid_min
           # v1.valid_max=v0.valid_max
            v1.missing_value=v0.missing_value
          #  v1._CoordinateAxes = "time lat lon "

            ntimes = len(f0.variables['time'][:])
            for i in xrange(ntimes):
                F0 = v0[i,:,:]
                # Extra trick if f0 has _FillValue
                F0 = np.where(F0.mask, FUNDEF, F0)
                F1 = sample2D(F0, X, Y)
                # hadde trengt et valid-range atributt, U og V er ofte < null
                #    F1[F1<=-10.0] = UNDEF
                #B1 = (F1 - v1.add_offset) / v1.scale_factor
                #B1[B1 < -1000] = UNDEF
                F1[F1 < -1000] = v1.add_offset + v1.scale_factor*UNDEF
                v1[i,:,:] = F1

        v0 = f0.variables[myvar]
      
        v1 = f1.createVariable(myvar, 'i2', ('time', 'depth', 'latitude', 'longitude'),fill_value=UNDEF)

        v1.standard_name=v0.standard_name
        v1.long_name = v0.long_name
        if (var=="vosaline"):
            v1.units = "psu"
        else:
            v1.units = v0.units
        v1.scale_factor=v0.scale_factor
        v1.add_offset=v0.add_offset
       # v1.valid_min=v0.valid_min
       # v1.valid_max=v0.valid_max
        v1.missing_value=v0.missing_value
      #  v1._CoordinateAxes = "time lev lat lon "


        # Modifications
        if myvar == 'vozocrtx':
             v1.standard_name = "eastward_sea_water_velocity"
             v1.field="eastward velocity, scalar, series"
        elif myvar == 'vomecrty':
             v1.standard_name = "northward_sea_water_velocity"
             v1.field="northward velocity, scalar, series"
        elif myvar == 'vosaline':
             v1.units = "1"
             v1.standard_name = "sea_water_salinity"
             v1.long_name="salinity"
             v1.field="salinity, scalar, series"
        elif myvar == 'votemper':
             v1.units = "degree_Kelvin"
             v1.standard_name = "sea_water_temperature"
             v1.long_name="temperature"
             v1.field="temperature, scalar, series"
        elif myvar== "shflux":
            v1.units=""
            v1.standard_name = "surface_downward_heat_flux_in_sea_water";

        elif myvar == 'zeta':
            # v1.units = "degree_Celsius"
             v1.standard_name = "sea_surface_height"
            # v1.long_name="temperature"
             v1.field="ssflux, scalar, series"

        ntimes = len(f0.variables['time'][:])
        for l in xrange(ntimes):
            Fz = np.squeeze(v0[l,:,:,:])
            Fz = np.where(Fz.mask, FUNDEF, Fz)
         
            for k in xrange(len(outlevels)):

                fieldDst[...] = v1.add_offset + v1.scale_factor*UNDEF

                weight=weights[outlevels[k]]
               
                ## Values are found at equal depth in infile as will be in target
                if len(weight)==2:
                    index=int(weight[0])

                    if inputgridEqualToOutputgrid:
                        fieldSrc[:,:] = Fz[index,:,:]
                        result=np.squeeze(Fz[index,:,:])
                    else:
                        fieldSrc[:,:] = np.fliplr(np.rot90(np.squeeze(Fz[index,:,:]),3))
                        F1 = regridSrc2Dst(fieldSrc, fieldDst, zero_region=ESMF.Region.SELECT)
                        result=F1.data
                    
                    print "Interpolating horisontally values for %s to find values at %s"%(inlevels[index],outlevels[k])
                   
                    result[result < -1000] =  v1.add_offset + v1.scale_factor*UNDEF
                    v1[l,k,:,:] = result

                else:
                    # Have to vertically interpolate to get the correct vertical depth
                    index1=int(weight[0])
                    index2=int(weight[1])
                    weight1=float(weight[2])
                    weight2=float(weight[3])
                    print "Interpolating vertically between %s and %s to find values at %s"%(inlevels[index1],inlevels[index2],outlevels[k])
                    if inputgridEqualToOutputgrid:
                        fieldSrc[:,:] = Fz[index1,:,:]
                        result1=np.squeeze(Fz[index1,:,:])
                    else:
                        fieldSrc[:,:] = np.fliplr(np.rot90(np.squeeze(Fz[index1,:,:]),3))
                        F1 = regridSrc2Dst(fieldSrc, fieldDst, zero_region=ESMF.Region.SELECT)
                        result1=F1.data
                    result1[result1 < -1000] =  v1.add_offset + v1.scale_factor*UNDEF
                    
                    if inputgridEqualToOutputgrid:
                        fieldSrc[:,:] = Fz[index1,:,:]
                        result2=np.squeeze(Fz[index2,:,:])
                    else:
                        fieldSrc[:,:] = np.fliplr(np.rot90(np.squeeze(Fz[index1,:,:]),3))
                        F2 = regridSrc2Dst(fieldSrc, fieldDst, zero_region=ESMF.Region.SELECT)
                        result2=F2.data
                    result2[result2 < -1000] =  v1.add_offset + v1.scale_factor*UNDEF
                    
                    v1[l,k,:,:] = result1*weight1 + result2*weight2
                    print  "Min %s and max %s weights %s and %s"%(np.min(v1[l,k,:,:]),np.max(v1[l,k,:,:]),weight1,weight2)
                # If the deepest depth in outlevels is deeper than what is found in inlevels,
                # the depth levels data array is set to undefined. We dont want to extrapolate values
                if outlevels[k] > np.max(inlevels):
                    print "Setting depth layer as undefined: deepest depth in outlevels is deeper than what is found in inlevels"
                    print "Inlevel max: %s outlevel: %s"%(inlevels[int(weight[0])],outlevels[k])
                    v1[l,k,:,:] = v1.add_offset + v1.scale_factor*UNDEF

        print "Results stored to: %s"%(outfilename)

        f0.close()
        f1.close()
    print "Finished\n"


