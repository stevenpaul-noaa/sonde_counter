#!/usr/bin/env python
import time
import sys
import os
from datetime import datetime
import re
'''
###########################

two optional arguments, mindate and maxdate in
yyyymmdd format.  

scans all drops in counter_all_drops_tail.txt 
(generated from sonde_counter_tial)

can read in names from counter_aliases.txt and print personal plots

2025-03-10 SCP: adjusted to read in tail numbers from all drops file, color coded dots based on tail number, added legend
2025-03-11 SCP: added selection of dates and tail numbers after running the script

needs:
numpy
matplotlib
cartopy
scipy

###########################
'''

# Assign unique colors for each tail number
tail_colors = {'N42RF': '#5CB200', 'N43RF': '#FF80A1', 'N49RF': '#305CDE'}

# Default color if an unexpected tail number appears
default_color = 'gray'

startscript=time.strftime('%Y%m%d_%H%M%S')
print("START")

mindate=19000101
maxdate = int(datetime.today().strftime('%Y%m%d'))  

if ( len(sys.argv) != 1) and ( len(sys.argv) != 3):
    print("USAGE: <python plot_tail.py> [yymmdd_startdate yymmdd_enddate ]")
    print("startdate defaults to 19010101, maxdate defaults to 99999999")
    exit()

if len(sys.argv) == 3:
    mindate=int(sys.argv[1])
    maxdate=int(sys.argv[2])
else:
    # Prompt the user if arguments are not provided
    mindate_input = input(f"Enter start date (YYYYMMDD) [default: {mindate}]: ").strip()
    if mindate_input:
        mindate = int(mindate_input)

    maxdate_input = input(f"Enter end date (YYYYMMDD) [default: {maxdate}]: ").strip()
    if maxdate_input:
        maxdate = int(maxdate_input)

print(f"Using date range: {mindate} to {maxdate}")

# Prompt user for aircraft selection
print("Select an option for aircraft filtering:")
print("1 - All aircraft (default)")
print("2 - All P-3 (N42RF and N43RF)")
print("3 - N42RF")
print("4 - N43RF")
print("5 - N49RF")
print("6 - N56RF")
print("7 - N677F (NCAR GV)")

selection = input("Enter your choice (1-7, default is 1): ").strip() or "1"

# Define aircraft filter based on selection
aircraft_filter = {
    "1": None,  # No filtering, include all
    "2": ["N42RF", "N43RF"],
    "3": ["N42RF"],
    "4": ["N43RF"],
    "5": ["N49RF"],
    "6": ["N56RF"],
    "7": ["N677F"]
}.get(selection, None)

# If the user enters an invalid selection, default to "All aircraft"
if aircraft_filter is None and selection not in {"1", "2", "3", "4", "5", "6", "7"}:
    print("Invalid selection. Defaulting to all aircraft.")
    aircraft_filter = None

f=open('counter_aliases.txt')
#names=['Steven Paul'] #can uncomment this line and comment six lines below to plot a single person
names=[]
names.append('all')
for line in f:
    if line[0] != '#' and line[0] != '\t':
        names.append(line.strip().strip(':'))
f.close()


import zipfile
ZipFile = zipfile.ZipFile("plot_"+startscript+".zip", "w" )
ZipFile.write(sys.argv[0],             compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write('counter_aliases.txt',   compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write('counter_all_drops_tail.txt', compress_type=zipfile.ZIP_DEFLATED)


for name in names:
    print('Starting User: '+ name)

    f=open('counter_all_drops_tail.txt','r',encoding="ascii", errors="surrogateescape")
    dropdata=[]
    for line in f:
        a=line.strip().split(',')
        tail = a[10].strip()

        # Filter based on selected aircraft
        if aircraft_filter is None or tail in aircraft_filter:
            if name == a[3] or name == 'all':
                #print(a[3],a[8],a[9])
                if float(a[8])!=0 and float(a[9])!=0:
                    if float(re.split(r"[_T]", a[0])[0]) >= mindate and float(re.split(r"[_T]", a[0])[0]) <= maxdate :
                        #dropdata is ( name, lat(float), lon(float) )
                        dropdata.append((a[3],float(a[8]),float(a[9]), a[10])) #appends name, lat, lon, tail number
    f.close()

    if not dropdata:
        print("No drops found for the selected aircraft/time range.")
        exit()

    lat_min=90
    lat_max=-90
    lon_min=999999
    lon_max=-999999
    
#    if len(dropdata) < 90:
#        continue

    #go through the drops
    for drop in dropdata:
        #id=drop[0]
        lat=drop[1]
        lon=(drop[2] + 360) %360
        if lon < 90:
            lon += 360
        #print('lat1: '+str(lat)+' ' + str(lat_min) + ' <-> '+ str(lat_max))
        #print('lon1: '+str(lon)+' ' + str(lon_min) + ' <-> '+ str(lon_max))
        if lat > lat_max:
            lat_max=lat
        if lat < lat_min:
            lat_min=lat
        if lon > lon_max:
            lon_max=lon
        if lon < lon_min:
            lon_min=lon
        #print('lat2: '+str(lat)+' ' + str(lat_min) + ' <-> '+ str(lat_max))
        #print('lon2: '+str(lon)+' ' + str(lon_min) + ' <-> '+ str(lon_max))
        tail = drop[3]


    lat_c=(lat_min+lat_max)/2
    lon_c=(lon_min+lon_max)/2

    lat_d=lat_c-lat_min
    lat_min=lat_min-0.1*lat_d
    lat_max=lat_max+0.1*lat_d

    lon_d=lon_c-lon_min
    lon_min=lon_min-0.01*lon_d
    if lon_min < -180:
        #lon_min=-180
        lon_min+=360
    lon_max=lon_max+0.01*lon_d


    #lat_min=-90
    #lat_max=90

    print('LAT_MIN: '+str(lat_min))
    print('LAT_CTR: '+str(lat_c))
    print('LAT_MAX: '+str(lat_max))
    print('LON_MIN: '+str(lon_min))
    print('LON_CTR: '+str(lon_c))
    print('LON_MAX: '+str(lon_max))

    import numpy as np
    import matplotlib.pyplot as plt
##    from mpl_toolkits.basemap import Basemap

    #fig = plt.figure(figsize=(15, 15))
    fig = plt.figure(dpi=200,figsize=(15, 15))
##    m = Basemap(projection='cyl', resolution='h', area_thresh=1000,
##                llcrnrlon=lon_min, llcrnrlat=lat_min,
##                urcrnrlon=lon_max, urcrnrlat=lat_max,
##                lat_0=lat_c, lon_0=lon_c)

                #lon_0=lon_c)

    import cartopy.crs as ccrs
    import cartopy.feature as cf

    #res='10m'
    res='50m'
    #res='110m'

    land=cf.NaturalEarthFeature('physical', 'land', res, edgecolor='face', facecolor='bisque')
    ocean=cf.NaturalEarthFeature('physical', 'ocean', res, edgecolor='face', facecolor='skyblue')
    coast=cf.NaturalEarthFeature('physical', 'coastline', res, edgecolor='black', facecolor='none')
    borders=cf.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', res, edgecolor='black', facecolor='none')

    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=lon_c))
    ax.set_extent([lon_min,lon_max,lat_min,lat_max])
    #ax.add_feature(cf.COASTLINE)
    #ax.add_feature(cf.BORDERS)
#    ax.add_feature(land)
#    ax.add_feature(ocean)
    #ax.add_feature(coast)
    #ax.add_feature(borders)

    #adds the image underneath, but crop goes away?
    #ax.stock_img()

    os.environ["CARTOPY_USER_BACKGROUNDS"] = "./background/"
    ax.background_img(name='BM', resolution='high') #resolution 'low' or 'high'

    xs=[]
    ys=[]
    
   
    #for drop in sorted(dropdata,reverse=True):
    for drop in sorted(dropdata):
        lat=drop[1]
        lon=(drop[2] + 360) %360
        if lon < 90:
            lon += 360
        tail = drop[3].strip()
        
        # Get the assigned color or use the default
        color = tail_colors.get(tail, default_color)
        plt.scatter(x=lon, y=lat, marker='o',s=1, color=color,transform=ccrs.PlateCarree())


    #plt.show()
    #exit()
    # Add a legend
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=tail_colors['N42RF'], markersize=5, label='N42RF', markeredgewidth=0),
                       Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=tail_colors['N43RF'], markersize=5, label='N43RF', markeredgewidth=0),
                       Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=tail_colors['N49RF'], markersize=5, label='N49RF', markeredgewidth=0),
                       Line2D([0], [0], marker='o', linestyle='None', markerfacecolor='gray', markersize=5, label='Other', markeredgewidth=0)]

    plt.legend(handles=legend_elements, loc='upper right')
        
    #fig.patch.set_facecolor('white')
    #fig.tight_layout(pad=0.5)
    #plt.gca().set_axis_off()
    #plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    #plt.margins(0,0)
    #plt.gca().xaxis.set_major_locator(plt.NullLocator())
    #plt.gca().yaxis.set_major_locator(plt.NullLocator())

    #if mindate != 19000101 or maxdate!=99999999:
    #    timeframe='_'+str(mindate)+'-'+str(maxdate)
    #else:
    #    timeframe=''

    # Determine the tail numbers to include in the filename
    if aircraft_filter and aircraft_filter != ["N42RF", "N43RF"]:  # If specific aircraft are selected
        tail_numbers_str = "_".join(aircraft_filter)
    else:
        tail_numbers_str = ""  # Don't add anything if "all aircraft" is selected

    # Modify the filename to include the tail numbers (only if not empty)
    if mindate != 19000101 or maxdate != int(datetime.today().strftime('%Y%m%d')):
        timeframe = f'_{mindate}-{maxdate}'
    else:
        timeframe = ''

    filename = f'plot_{startscript}_drops_{name.lower().replace(" ", "_").replace("?", "").replace("/", "")}{timeframe}'

    # Add tail numbers to the filename if not empty (only for selected aircraft)
    if tail_numbers_str:
        filename += f"_{tail_numbers_str}"

    filename += ".png"
    
    plt.savefig(filename, facecolor='red', bbox_inches='tight',pad_inches=0)
    plt.close()
    
    ZipFile.write(filename, compress_type=zipfile.ZIP_DEFLATED)
    
    #comment out below to get personalized plots
    if name=='all':
       exit()

exit()
