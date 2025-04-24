#!/usr/bin/env python



'''

2015-07-22 CFL: Updated to work on windows and python 3.4

2015-10-29 CFL:
updated to find dupes based on serial + operator name instead of date+operator
added output of bad and duplicate files
all_drops now outputs more data than before

2025-03-06 SCP: Updated to append tail number after lat and long



USAGE: python <this_file_name>.py "Start Directory"
Should have a valid 'aliases.txt' file

This program will scan through a tree of files looking 'R' AVAPS files 
and proceed to parse them for username and compile a total

drops with the same date/time and same operator are counted as duplicates, 
even if on different channels


'''
from __future__ import print_function

import os
import sys
import string
import operator
import time
import argparse
from netCDF4 import Dataset
from datetime import datetime
import json
import re



#read aliases from text file.  if a line begins with \t it is an alias for the previous name
def load_aliases(aliases):
    alias=""
    enddate = 0
    startdate = 0
    dropname = ""

    f=open('counter_aliases.txt','r')
    o=open(startscript+'_counter_aliases_out.txt','w')
    oo=open(startscript+'_counter_dropper_names.txt','w')
    ooo=open(startscript+'_counter_aliases_sorted.txt','w')
    line=f.readline().rstrip('\r\n')
    #o.write(line[:-1]+'\n')
    o.write(line+'\n')
    while line != "":
        if line[0] == "#":
            pass
        elif line[0] == "\t":
            splat = line[1:].split('\t',3)
            if len(splat) == 1:
                enddate = 99999999
                startdate = 0
                dropname=splat[0].upper()
            if len(splat) == 2:
                dropname = splat[0].upper()
                startdate = int(splat[1].replace('-',''))
                enddate = 99999999
            if len(splat) == 3:
                dropname = splat[0].upper()
                startdate = int(splat[1].replace('-',''))
                enddate = int(splat[2].replace('-',''))
            
            if dropname=='\'\'':
                dropname=''
                
            if dropname in aliases:
                entries = aliases[dropname][0]
                for entry in range(entries):
                    if (startdate >= aliases[dropname][1][entry][1]) and (startdate < aliases[dropname][1][entry][2]):
                        print("line: " + line)
                        print("OVERLAPPING DATE RANGES")
                        print("name: " + dropname)
                        print("startdate: " + str(startdate))
                        print("prevst: " + str(aliases[dropname][1][entry][1]))
                        print("prevend: " + str(aliases[dropname][1][entry][2]))
                        exit()
                    if (enddate > aliases[dropname][1][entry][1]) and (enddate <= aliases[dropname][1][entry][2]):
                        print("line: " + line)
                        print("OVERLAPPING DATE RANGES")
                        print("OVERLAPPING DATE RANGES")
                        print("name: " + dropname)
                        print("startdate: " + str(startdate))
                        print("prevst: " + str(aliases[dropname][1][entry][1]))
                        print("prevend: " + str(aliases[dropname][1][entry][2]))
                        exit()
                    if (startdate <= aliases[dropname][1][entry][1]) and (enddate >= aliases[dropname][1][entry][2]):
                        print("line: " + line)
                        print("OVERLAPPING DATE RANGES")
                        print("OVERLAPPING DATE RANGES")
                        print("name: " + dropname)
                        print("startdate: " + str(startdate))
                        print("prevst: " + str(aliases[dropname][1][entry][1]))
                        print("prevend: " + str(aliases[dropname][1][entry][2]))
                        exit() 
                aliases[dropname][0] += 1
                aliases[dropname][1].append((alias,startdate,enddate))
            else:
                aliases[dropname]=[1,[(alias,startdate,enddate)]]
                
        elif line[0] != "\n" and line[0] != '\r':
            alias=line.rstrip()
            alias=alias.rstrip(':')
        line=f.readline().rstrip()
        #o.write(line[:-1]+'\n')
        o.write(line+'\n')

    for name in sorted(aliases.keys()):
        for entry in range(aliases[name][0]):
            tmpstr="%30s  %30s %15d %15d\n" % (name, aliases[name][1][entry][0],aliases[name][1][entry][1],aliases[name][1][entry][2])
            oo.write(tmpstr)
    real_names={}
    for avaps_name in sorted(aliases.keys()):
        for entry in range(aliases[avaps_name][0]):
            real_name=aliases[avaps_name][1][entry][0]
            if not real_name in real_names:
                real_names[real_name]=[]
            
            start_date=aliases[avaps_name][1][entry][1]
            start_date_str=str(start_date)[0:4]+'-'+str(start_date)[4:6]+'-'+str(start_date)[6:8]
            end_date=aliases[avaps_name][1][entry][2]
            end_date_str=str(end_date)[0:4]+'-'+str(end_date)[4:6]+'-'+str(end_date)[6:8]
            if start_date==0 and end_date==99999999:
                real_names[real_name].append(avaps_name)
            elif end_date==99999999:
                real_names[real_name].append(avaps_name+'\t'+start_date_str)
            else:
                real_names[real_name].append(avaps_name+'\t'+start_date_str+'\t'+end_date_str)
            

    for real_name in sorted(real_names):
        ooo.write(real_name+':\n')
        for rna in sorted(real_names[real_name]):
            if rna=='':
                ooo.write('\t\'\'\n')
            else:
                ooo.write('\t'+rna+'\n')

    f.close()
    o.close()
    oo.close()
    ooo.close()
    #~ print('CREATED TXT FILE \''+startscript+'_aliases.txt\'')
    #~ print('CREATED TXT FILE \''+startscript+'_dropper_names.txt\'')



#read a D file and add to the totals in counts and add an entry to uid (if unique)
def process_file(root, filename, aliases, minmaxdate, dropdata, tots, badfiles, dupefiles, seen_ids,seen_uids):

    year='0000'
    if len(filename) == 16:
        filedate="19"+filename[1:14]
    else:
        filedate=filename[1:16]
    year=filedate[0:4]
    
    if not(filename[-1].isdigit() and filedate[0:8].isdigit() and filedate[9:13].isdigit()):
        tots['bad'] += 1
        badfiles.append(filename + ': BAD FILENAME: '+root+'\\'+filename)
        return


    #get a few files with bad date in filename
    if year < '1996':
        tots['bad'] += 1
        badfiles.append(filename + ': BAD FILENAME YEAR: '+root+'\\'+filename)
        return

    f=open(root+'\\'+filename,'r',encoding="ascii", errors="surrogateescape")
    text = f.read()
    f.close()
    

    #look for either of these strings to shorten the subsequent searches
    [_0,_,_1] = text.partition("Data Type:")
    if _1 != "":
        text=_1
        pretext=_0
    
    [_0,_,_1] = text.partition("Data Type/Data Channel:")
    if _1 != "":
        text=_1
        pretext=_0

    text=text.replace('\r','')

    #Get launch time fr0m file
    [_,_,_1] = text.partition("Launch Time (y,m,d,h,m,s):")
    if _1 == "":
        tots['bad'] += 1
        badfiles.append(filename + ': NO LAUNCH TIME: '+root+'\\'+filename)
        return
    else:
        [launch_date,_,_] = _1.partition("\n")
        launch_date=launch_date.strip(' ')
        launch_date=launch_date.strip('\r')
        launch_date=launch_date.replace('/','')
        launch_date=launch_date.replace('-','')
        launch_date=launch_date.replace(',','')
        launch_date=launch_date.replace(':','')
        launch_date=launch_date.replace(' ','_')
        [launch_date,_,_] = launch_date.partition('.')
        launch_date_int = int(launch_date[0:8])

    filedate_int = int(filedate[0:8])
    #check launch date and filedate to be the same
    #~ if launch_date_int != filedate_int:
        #~ tots['bad'] += 1
        #~ badfiles.append(filename + ': Launch Date and filedate mismatch: LD: '+launch_date+' FD: '+filedate+' '+root+'\\'+filename)
        #~ return

    #check for date range
    if launch_date_int < minmaxdate[0] or launch_date_int > minmaxdate[1]:
        tots['bad'] += 1
        badfiles.append(filename + ': DATE OUT OF RANGE: ' +root+'\\'+filename)
        return

    nd_text=text.upper()
    if 'NO DROP' in nd_text:
        #print('No Drop')
        tots['bad'] += 1
        badfiles.append(filename + ': "NO DROP" in text: '+root+'\\'+filename)
        return
    if 'NO SOUNDING' in nd_text:
        tots['bad'] += 1
        badfiles.append(filename + ': "NO SOUNDING" in text: '+root+'\\'+filename)
        return
    if 'DID NOT DROP' in nd_text:
        tots['bad'] += 1
        badfiles.append(filename + ': "DID NOT DROP" in text: '+root+'\\'+filename)
        return
    if 'GROUND TEST' in nd_text:
        tots['bad'] += 1
        badfiles.append(filename + ': "GROUND TEST" in text: '+root+'\\'+filename)
        return
    if 'NO LAUNCH' in nd_text and not "NO LAUNCH DETECT" in nd_text:
        tots['bad'] += 1
        badfiles.append(filename + ': "NO LAUNCH" in text: '+root+'\\'+filename)
        return
    if 'SOUNDING ABORTED' in nd_text:
        tots['bad'] += 1
        badfiles.append(filename + ': "SOUNDING ABORTED" in text: '+root+'\\'+filename)
        return
    if 'DID NOT LAUNCH' in nd_text:
        tots['bad'] += 1
        badfiles.append(filename + ': "DID NOT LAUNCH" in text: '+root+'\\'+filename)
        return
    [_,_,_1] = nd_text.partition("DROP NAME:")
    if _1 != '':
        if _1.partition("\n")[0].strip(' ') == 'GROUND':
            tots['bad'] += 1
            badfiles.append(filename + ': Drop Name is GROUND: '+root+'\\'+filename)
            return
    [_,_,_1] = nd_text.partition("SOUNDING NAME:")
    if _1 != '':
        if _1.partition("\n")[0].strip(' ') == 'GROUND':
            tots['bad'] += 1
            badfiles.append(filename + ': Sounding Name is GROUND: '+root+'\\'+filename)
            return
        
    #tail
    tail="-----"
    [_,_,_1] = text.partition("Aircraft Type/ID:")
    if _1 != '':
        [onc,_,_] = _1.partition("\n")
        tail=onc.split(",",1)[1]

    #name can be in one of 3 fields depending on format version
    name='---'
    [_,_,_1] = text.partition("Operator Name/Comments:")
    if _1 != '':
        [onc,_,_] = _1.partition("\n")
        name=onc.split(",",1)[0]

    [_,_,_1] = text.partition("System Operator:")
    if _1 != '':
        [name,_,_]=_1.partition("\n")

    [_,_,_1] = text.partition("Operator Name:")
    if _1 != '':
        [name,_,_]=_1.partition("\n")

    #no name found, skip file and pring 'X'
    if name=="---":
        tots['bad'] += 1
        badfiles.append(filename + ': NO NAME: '+root+'\\'+filename)
        return

    #uppercase name and look for in it aliases file, substitute name if found
    name=name.replace('\r','')
    name=name.replace('\n','')
    name=name.strip(' ')
    name=name.upper()
    origname=name
    name="'"+name+"'"
    if origname in aliases:
        for entry in range(aliases[origname][0]):
            if launch_date_int >= aliases[origname][1][entry][1] and launch_date_int <= aliases[origname][1][entry][2]:
                name = aliases[origname][1][entry][0]

    id='none'

    [_,_,_1] = text.partition('Sonde Type/ID/Tx Frequency:')
    if _1 != '':
        [id,_,_]=_1.partition("\n")
        id=id.split(',')[1]

    [_,_,_1] = text.partition('Sonde Type/ID/Ver/Tx Frequency:')
    if _1 != '':
        [id,_,_]=_1.partition("\n")
        id=id.split(',')[1]
    
    
    [_,_,_1] = text.partition('Sonde Type/ID/Ver/Battery/Tx Freq:')
    if _1 != '':
        [id,_,_]=_1.partition("\n")
        id=id.rsplit(',',4)[1]
    
    [_,_,_1] = text.partition('Sonde ID/ID$/Built/Ver/Type:')
    if _1 != '':
        [id,_,_]=_1.partition("\n")
        id=id.split(',')[0] 
    
    [_,_,_1] = text.partition('Sonde ID/ID$/Built/Firmware/Type:')
    if _1 != '':
        [id,_,_]=_1.partition("\n")
        id=id.split(',')[0] 

        
    [_,_,_1] = text.partition('Sonde ID/Built/Firmware/Type:')
    if _1 != '':
        [id,_,_]=_1.partition("\n")
        id=id.split(',')[0] 
    
    [_,_,_1] = text.partition('Sonde ID/Type/Rev/Built/Sensors:')
    if _1 != '':
        [id,_,_]=_1.partition("\n")
        id=id.split(',')[0] 
    
        
    id=id.strip()       
    
    if not ( 7 <= len(id) <= 10):
        tots['bad'] += 1
        badfiles.append(filename + ': BAD ID: '+id+' '+root+'\\'+filename)
        return
    elif id == 'none':
        tots['bad'] += 1
        badfiles.append(filename + ': BAD ID: '+id+' '+root+'\\'+filename)
        return
    
    if pretext.count('\n') < 10:
        #print("Data Len Good: ",filename,pretext.count('\n'))
        tots['bad'] += 1
        badfiles.append(filename + ': SHORT DATA PORTION: '+repr(pretext.count('\n'))+' '+root+'\\'+filename)
        return

    [_,_,_1] = nd_text.partition("PRE-LAUNCH OBS (LON,LAT,ALT):")
    if _1 == '':
        tots['bad'] += 1
        badfiles.append(filename + ': No Pre-Launch: '+root+'\\'+filename)
        return  
    else:
#        [_1,_,_] = _1.partition("\n")
#        _1=_1.strip(' ')
#        _1=_1.split(',')[2]
#        _1=_1.strip(' ')
#        _1=_1.strip(' M')
#        _1=_1.split('.')[0]
#        if int(_1) < 100 :
#            #print(int(_1))
#            tots['bad'] += 1
#            badfiles.append(filename + ': Launch Alt too low: '+_1+'meters '+root+'\\'+filename)
#            return  
        [_1,_,_] = _1.partition("\n")
        #_1=_1.strip(' ')
        _1=_1.split(',')
        if not len(_1) >= 3:
            badfiles.append(filename + ': BAD LAT_LON_ALT: '+root+'\\'+filename)
            return
        lon,lat,alt=_1[0:3]
        alt=alt.strip(' ')
        alt=alt.strip(' M')
        alt=alt.split('.')[0]
        if int(alt) < 100 :
            #print(int(_1))
            tots['bad'] += 1
            badfiles.append(filename + ': Launch Alt too low: '+str(alt)+' meters '+root+'\\'+filename)
            return
        lat=float(lat.split(' DEG')[0])
        lon=float(lon.split(' DEG')[0])
        if (lat > 90) or (lat < -90):
            lat=0
            lon=0
        elif (lon > 180) or (lon < -180):
            lat=0
            lon=0

    tots['good'] += 1
    if not id in seen_ids:
        seen_ids.append(id)

    UID=id+'_'+name
    if not UID in seen_uids:
        seen_uids.append(UID)
        dropdata.append((UID,filedate,name,origname,id,root+'\\'+filename,filename,lat,lon,tail))
    else:
        tots['dupe'] += 1
        dupefiles.append('DUPE: OF: '+dropdata[seen_uids.index(UID)][6]+' '+UID+' DF: '+ filename + ' NAME: '+name+' DUPE:'+root+'\\'+filename+ ' ORIG: '+dropdata[seen_uids.index(UID)][2]+' '+dropdata[seen_uids.index(UID)][5])
        return
    
    return

#read a netCDF file and add to the totals in counts and add an entry to uid (if unique)
def process_file_nc(root, filename, aliases, minmaxdate, dropdata, tots, badfiles, dupefiles, seen_ids,seen_uids):

    year='0000'
    
    # Regular expression to match the expected format
    pattern = r"^([A-Za-z0-9_-]+)-([A-Za-z0-9_-]+)-([0-9]+)-?(\d{8}T\d{6})?-([A-Za-z0-9_-]+)\.nc$"
    match = re.match(pattern, os.path.basename(filename))
    
    if match:
        project_name = match.group(1)
        mission_id = match.group(2)
        drop_number = match.group(3)
        filedate = match.group(4)  # YYYYMMDDTHHMMSS
        channel = match.group(5)
        if filedate == None:
            tots['bad'] += 1
            badfiles.append(filename + ': NO LAUNCH TIME: '+root+'\\'+filename)
            return
        else:
            launch_date_int = int(filedate[0:8])
            year=filedate[0:4] 

    else:
        tots['bad'] += 1
        badfiles.append(filename + ': BAD FILENAME: '+root+'\\'+filename)
        return
    
    #check for date range
    if launch_date_int < minmaxdate[0] or launch_date_int > minmaxdate[1]:
        tots['bad'] += 1
        badfiles.append(filename + ': DATE OUT OF RANGE: ' +root+'\\'+filename)
        return

    # import data from netcdf

    # set default values
    name = ''
    tail = 'No Tail Number'
    serial = 'No Serial Number'
    lat = 0
    lon = 0

    # read attributes
    try:    
        with Dataset(root+'\\'+filename, 'r') as nc:
            attrs={}
            for attr in nc.ncattrs():  # ncattrs() is the proper method for global attributes
                try:
                    attrs[attr] = nc.getncattr(attr)
                except AttributeError:
                    print(f"Warning: Global attribute '{attr}' is missing or unreadable in {filename}")
    except OSError as e:
        tots['bad'] += 1
        badfiles.append(filename + ': ERROR OPENING NETCDF FILE: ' +root+'\\'+filename)
        return

    name = attrs.get('DropOperator','')
    tail = attrs.get('DropTailNumber','No Tail Number')
    id = attrs.get('SerialNumber','No Serial Number')
    drop_launch_obs_str = attrs.get('DropLaunchObs', '{}')
    if drop_launch_obs_str.strip():
        try:
            drop_launch_obs = json.loads(drop_launch_obs_str)
        except json.JSONDecodeError:
            # If it's not JSON, assume it's a comma-separated string
            drop_launch_obs_list = drop_launch_obs_str.split(',')
        
            # Extract values based on known positions
            drop_launch_obs = {
                "latitude": float(drop_launch_obs_list[2]) if len(drop_launch_obs_list) > 2 and drop_launch_obs_list[2] else None,
                "longitude": float(drop_launch_obs_list[3]) if len(drop_launch_obs_list) > 3 and drop_launch_obs_list[3] else None,
            }
    else:
        drop_launch_obs = {}  # fallback to an empty dict
    lat = float(drop_launch_obs.get('latitude', 0))
    lon = float(drop_launch_obs.get('longitude', 0))

    #uppercase name and look for in it aliases file, substitute name if found
    name=name.replace('\r','')
    name=name.replace('\n','')
    name=name.strip(' ')
    name=name.upper()
    origname=name
    name="'"+name+"'"
    if origname in aliases:
        for entry in range(aliases[origname][0]):
            if launch_date_int >= aliases[origname][1][entry][1] and launch_date_int <= aliases[origname][1][entry][2]:
                name = aliases[origname][1][entry][0]

            

    id=id.strip()       
    
    if not ( 7 <= len(id) <= 10):
        tots['bad'] += 1
        badfiles.append(filename + ': BAD ID: '+id+' '+root+'\\'+filename)
        return
    elif id == 'No Serial Number':
        tots['bad'] += 1
        badfiles.append(filename + ': NO SERIAL NUMBER: '+id+' '+root+'\\'+filename)
        return

    if (lat > 90) or (lat < -90):
        lat=0
        lon=0
    elif (lon > 180) or (lon < -180):
        lat=0
        lon=0

    if lat == 0 and lon == 0:
        tots['bad'] += 1
        badfiles.append(filename + ': BAD LAT_LON_ALT: '+root+'\\'+filename)
        return

    tots['good'] += 1
    if not id in seen_ids:
        seen_ids.append(id)

    UID=id+'_'+name
    if not UID in seen_uids:
        seen_uids.append(UID)
        dropdata.append((UID,filedate,name,origname,id,root+'\\'+filename,filename,lat,lon,tail))
    else:
        tots['dupe'] += 1
        dupefiles.append('DUPE: OF: '+dropdata[seen_uids.index(UID)][6]+' DF: '+ filename + ' NAME: '+name+' DUPE:'+root+'\\'+filename+ ' ORIG: '+dropdata[seen_uids.index(UID)][2]+' '+dropdata[seen_uids.index(UID)][5])
        return
    
    return




###########################
#  START POINT
###########################

startscript=time.strftime('%Y%m%d_%H%M%S')
print("START")
if (len(sys.argv) < 2) or (len(sys.argv)==3) or (len(sys.argv) > 4):
    print ("Num Args: "+str(len(sys.argv)))
    for arg in sys.argv:
        print("arg: " + arg)
    print ("1) Count all D files in <startdir>:")
    print ("   <python> <scriptname> <startdir>")
    print ("2) Count only files in range min to max:")
    print ("   <python> <scriptname> <startdir> [min_yyymmdd max_yyymmdd]")
    exit()

from shutil import copyfile
copyfile(sys.argv[0],startscript + '_' + os.path.basename(sys.argv[0]))

start_dir = sys.argv[1]
#if start_dir[-1] != '/':
#   start_dir = start_dir + '/'
if start_dir[-1] != '\\':
    start_dir = start_dir + '\\'

print("STARTING DIR: " + start_dir)

#get min/max dates from cmd line
mindate_str = '00000000'
maxdate_str = '99999999'
if len(sys.argv) == 4:
    mindate_str = sys.argv[2]
    maxdate_str = sys.argv[3]
minmaxdate=(int(mindate_str.replace('-','')), int(maxdate_str.replace('-','')))


#load aliases from aliases.txt
aliases={}
load_aliases(aliases)

#list with messages about bad/duplicate files to be printed later
badfiles=[]
dupefiles=[]

#list of seen sonde id's for total sondes count
#1 sonde can have drops from it for different users
seen_ids=[]

#list of seen uids to determine if drop is a dupe or not
seen_uids=[]

#main list with drop information
dropdata=[]

#scan full directory tree and process matching files.
tots={'good':0,'dupe':0,'bad':0}
import os
time0=time.time()
time1=time0

os.system("")
filecount=0
print('SCANNING FOR FILECOUNT...',end='\r')
for root,subdirs,files in os.walk(start_dir):
    for filename in files:
        if filename.startswith("D"):
            if len(filename) in [16, 18]:
                if filename[len(filename)-2] == ".":
                    filecount+=1
        if filename.endswith(".nc"):
            if filename[-4].isdigit():
                filecount+=1

tot_filecount=filecount
print('\x1b[2KTOTAL FILES: {:,}'.format(filecount))

filecount=0

for root,subdirs,files in os.walk(start_dir):
    for filename in files:
        # process files that start with D, are 16-18 chars long, and second to last char is .
        if filename.startswith("D"):
            if len(filename) in [16, 18]:
                if filename[len(filename)-2] == ".":
                    filecount+=1
                    process_file(root, filename, aliases, minmaxdate, dropdata, tots, badfiles, dupefiles, seen_ids,seen_uids)
                    time2=time.time()
                    if (time2 - time1) > 0.25:
                        print('\x1b[2KPROCESSED: {:,}/{:,} (G: {:,} D: {:,} B: {:,}) {:s}'.format(filecount,tot_filecount, tots['good'],tots['dupe'], tots['bad'],root[len(start_dir):len(root)]+' '+filename), end='\r')
                        sys.stdout.flush()
                        time1 = time2
        # process files that end with .nc and the char before the .nc is a number (channel number) - ignores mission.nc files
        if filename.endswith(".nc"):
            if filename[-4].isdigit():
                filecount+=1
                process_file_nc(root, filename, aliases, minmaxdate, dropdata, tots, badfiles, dupefiles, seen_ids,seen_uids)
                time2=time.time()
                if (time2 - time1) > 0.25:
                    print('\x1b[2KPROCESSED: {:,}/{:,} (G: {:,} D: {:,} B: {:,}) {:s}'.format(filecount,tot_filecount, tots['good'],tots['dupe'], tots['bad'],root[len(start_dir):len(root)]+' '+filename), end='\r')
                    sys.stdout.flush()
                    time1 = time2

                    
print('\x1b[2KPROCESSED: {:,}/{:,} (G: {:,} D: {:,} B: {:,}) {:s}'.format(filecount,tot_filecount, tots['good'],tots['dupe'], tots['bad'],root[len(start_dir):len(root)]+' '+filename), end='\r')
print("");



#initialize vars
seen=list()
grand_total={}

#sort dropdata
dropdata_sorted = sorted(dropdata,key=lambda tup: tup[1])

#catch when no drops were found
if len(seen_uids) == 0:
    print('NO VALID DROPS FOUND')
    exit()
latest_date=dropdata_sorted[-1][1][0:4]+'-'+dropdata_sorted[-1][1][4:6]+'-'+dropdata_sorted[-1][1][6:8]

#first and last year/fy seen
earliest_year=int(dropdata_sorted[0][1][0:4])
latest_year=int(dropdata_sorted[-1][1][0:4])
num_years = latest_year - earliest_year + 2

earliest_fy=int(dropdata_sorted[0][1][0:4])
if int(dropdata_sorted[0][1][4:6]) >= 10:
    earliest_fy += 1
latest_fy=int(dropdata_sorted[-1][1][0:4])
if int(dropdata_sorted[-1][1][4:6]) >= 10:
    latest_fy += 1
num_fy = latest_fy - earliest_fy + 2

#get total number of sonde launches (not drops)
sondes_launched=len(seen_ids)

highest_total=0
highest_name=""
highest_name_len=0
highest_date=""
max_name_len=0
top_dropper=[]

#totals is a dict for holding
totals={}
totals['YEAR']=[0]*num_years
totals['YEAR'][num_years - 1]='TOTAL'
totals['----']=['----']*num_years
totals['----'][num_years - 1]='------'

#totals_fy is a dict for storing the fiscal year totals
totals_fy={}
totals_fy['FY']=[0]*num_fy
totals_fy['FY'][num_fy- 1]='TOTAL'
totals_fy['----']=['----']*num_fy
totals_fy['----'][num_fy - 1]='------'

outfile=open(startscript+'_counter_all_drops_tail.txt','w')
dropdata=sorted(dropdata)
dropdata_tot=len(dropdata)
dropdata_count=0
for drop in dropdata_sorted:
    print('\x1b[2KPROCESSED: {:,}/{:,} Drops'.format(dropdata_count,dropdata_tot), end='\r')
    dropdata_count=dropdata_count + 1;
    
    #write drop date 1st
    outfile.write(drop[1]+',')
    
    #loop through a drop entry and write the data
    for entry in range(len(drop)):
        outfile.write(str(drop[entry]))
        if entry < (len(drop)-1):
            outfile.write(',')
    outfile.write('\n')
    
    #get current year/fiscal year
    year=int(drop[1][0:4])
    fy=year
    if int(drop[1][4:6]) >= 10:
        fy+= 1
    
    #get aliased name and upate the field length if necessary
    name=drop[2]
    if len(name) > max_name_len:
        max_name_len=len(name)
    
    #if there is not a year/fy column header for current drop, add it
    if not year in totals['YEAR']:
        totals['YEAR'][year-earliest_year] = year
    if not fy in totals_fy['FY']:
        totals_fy['FY'][fy-earliest_fy] = fy

    #if there is not an entry for current drops' user, add it
    if not name in totals:
        totals[name]=['.']*(num_years)
    if not name in totals_fy:
        totals_fy[name]=['.']*(num_fy)
        
    #either convert current year total from . to 1 or add 1
    if totals[name][year-earliest_year] == '.':
        totals[name][year-earliest_year] = 1
    else:
        totals[name][year-earliest_year] += 1

    #either convert current fy total from . to 1 or add 1
    if totals_fy[name][fy-earliest_fy] == '.':
        totals_fy[name][fy-earliest_fy] = 1
    else:
        totals_fy[name][fy-earliest_fy] += 1

    #either convert current grand total from . to 1 or add 1
    if totals[name][num_years -1] == '.':
        totals[name][num_years -1] = 1
    else:   
        totals[name][num_years -1] += 1

    #either convert current grand fy total from . to 1 or add 1
    if totals_fy[name][num_fy -1] == '.':
        totals_fy[name][num_fy -1] = 1
    else:   
        totals_fy[name][num_fy -1] += 1

    #if the current drop's users' total is higher than the top dropper
    #update the top dropper field
    if totals[name][num_years -1] > highest_total:
        if len(name) > highest_name_len:
            highest_name_len = len(name)
        highest_total=totals[name][num_years -1]
        highest_date=drop[1][0:4]+'-'+drop[1][4:6]+'-'+drop[1][6:8]+' '+drop[1][9:11]+':'+drop[1][11:13]+':'+drop[1][13:15]
        
        #if the top dropper's name didnt change, dont update the list
        #but store to write at very end
        if highest_name != name:
            highest_name=name
            top_dropper.append((totals[highest_name][num_years -1], highest_name, highest_date ))

print('\x1b[2KPROCESSED: {:,}/{:,} Drops'.format(dropdata_count,dropdata_tot), end='\r')
print("")

#add the latest entry of the top dropper
top_dropper.append((totals[highest_name][num_years -1], highest_name, highest_date ))

#update the dashed line field to be as long as the longest name
totals['-'*max_name_len]=totals.pop('----')
totals_fy['-'*max_name_len]=totals_fy.pop('----')

#close out the all drops file
outfile.close()
#~ print('CREATED TXT FILE \''+startscript+'_all_drops.txt\'')

#sort totals by max col
#print("num years: " + str(num_years) + " --- " + str(type(num_years)))
#for ite in totals:
#    print(totals[ite][num_years - 1])

#remove two heading lines

#print(len(totals.keys()))

HDR1=totals.pop('YEAR')
HDR2=totals.pop('-'*max_name_len)
sorted_totals=sorted(totals.items(),key=lambda x:x[1][num_years - 1],reverse=True)
sorted_totals.insert(0,('-'*max_name_len,HDR2))
sorted_totals.insert(0,('YEAR',HDR1))
totals['YEAR']=HDR1
totals['-'*max_name_len]=HDR2

HDR1=totals_fy.pop('FY')
HDR2=totals_fy.pop('-'*max_name_len)
sorted_totals_fy=sorted(totals_fy.items(),key=lambda x:x[1][num_fy - 1],reverse=True)
sorted_totals_fy.insert(0,('-'*max_name_len,HDR2))
sorted_totals_fy.insert(0,('FY',HDR1))
totals_fy['FY']=HDR1
totals_fy['-'*max_name_len]=HDR2


#write the bad files data to a file
outfile=open(startscript+'_counter_bad_files.txt','w')
badfiles.sort()
for line in badfiles:
    outfile.write(line+'\n')
outfile.close()
#~ print('CREATED TXT FILE \''+startscript+'_bad_files.txt\'')


#write the duplicates data to a file
outfile=open(startscript+'_counter_duplicates.txt','w')
dupefiles.sort()
for line in dupefiles:
    outfile.write(line+'\n')
outfile.close()
#~ print('CREATED TXT FILE \''+startscript+'_duplicates.txt\'')


#write the gand totals to a file
outfile=open(startscript+'_counter_summary.txt','w')
outfile.write("\n")
outfile.write("-----------------------\n")
outfile.write("TOTALS BY CALENDAR YEAR\n")
outfile.write("-----------------------\n")
outfile.write("\n")
for line in sorted_totals:
    if line[0] == 'N/A':
        continue
    #print the name
    if line[0] == 'YEAR':
        tmpstr='%-*s' % (max_name_len,'NAME')
    else:
        tmpstr='%-*s' % (max_name_len,line[0])
    
    #loop through year entires for current user (skip last col)
    for year in range(len(line[1])-1):
        tmpstr='%s %4s' % (tmpstr,line[1][year])
    
    #for first two line, print the last col as string, others use comma
    # seperated numbers
    if line[0] == 'YEAR' or line[0][0] == '-':
        tmpstr='%s %6s' % (tmpstr,line[1][len(line[1])-1])  
    else:
        tmpstr='%s {:6,}'.format(line[1][len(line[1])-1]) % (tmpstr)    
    outfile.write(tmpstr+'\n')



#write the gand fy totals to a file
outfile.write("\n")
outfile.write("---------------------\n")
outfile.write("TOTALS BY FISCAL YEAR\n")
outfile.write("---------------------\n")
outfile.write("\n")
for line in sorted_totals_fy:
    if line[0] == 'N/A':
        continue

    #print the name
    if line[0] == 'FY':
        tmpstr='%-*s' % (max_name_len,'NAME')
    else:
        tmpstr='%-*s' % (max_name_len,line[0])
    
    #loop through year entires for current user (skip last col)
    for year in range(len(line[1])-1):
        #change year names on 1st line only
        if line[0] == 'FY':
            tmpstr='%s FY%2s' % (tmpstr,repr(line[1][year])[2:4])
        else:
            tmpstr='%s %4s' % (tmpstr,line[1][year])
    
    #for first two line, print the last col as string, others use comma
    # seperated numbers
    if line[0] == 'FY' or line[0][0] == '-':
        tmpstr='%s %6s' % (tmpstr,line[1][len(line[1])-1])  
    else:
        tmpstr='%s {:6,}'.format(line[1][len(line[1])-1]) % (tmpstr)    
    outfile.write(tmpstr+'\n')



#write the top dropper info to a file
outfile.write("\n")
outfile.write("-----------------\n")
outfile.write("TOP DROPPER TITLE\n")
outfile.write("-----------------\n")
outfile.write("\n")
outfile.write('%5s  %-*s  %s\n' % (' DROPS',highest_name_len,'NAME','DATE/TIME'))
outfile.write('%5s  %-*s  %s\n' % ('------',highest_name_len,'-'*highest_name_len,'-'*19))
for line in top_dropper:
    outfile.write('{:6,}  %-*s  %s\n'.format(line[0]) % (highest_name_len,line[1],line[2]))


#print some final summary information
outfile.write("\n")
outfile.write("Processing Time: %dm %5.2fs\n" % divmod(time.time() - time0, 60))
md='%08d' % minmaxdate[0]
md=md[0:4]+'-'+md[4:6]+'-'+md[6:8]
outfile.write("Minimum Date:    %s\n" % md)
md='%08d' % minmaxdate[1]
md=md[0:4]+'-'+md[4:6]+'-'+md[6:8]
outfile.write("Maximum Date:    %s\n" % md)
outfile.write("Latest Data:     %s\n" % latest_date)
outfile.write("Sondes Launched: {:,}\n".format(sondes_launched))

outfile.write("\n")
outfile.write("Good Files:  {:6,}\n".format(tots['good']))
outfile.write("Dupe Files:  {:6,}\n".format(tots['dupe']))
outfile.write("Bad Files:   {:6,}\n".format(tots['bad']))
outfile.write("TOTAL FILES: {:6,}\n".format(tots['good']+tots['bad']+tots['dupe']))
outfile.close()
#~ print('CREATED TXT FILE \''+startscript+'_summary.txt\'')





#print each FY and it's total
outfile=open(startscript+'_counter_totals_fy.txt','w')
outfile.write("\n")
outfile.write("----------------------------\n")
outfile.write("YEARLY TOTALS BY FISCAL YEAR\n")
outfile.write("----------------------------\n")
for year_idx in range(len(totals_fy['FY'])-1):
    outfile.write("\n")
    HDR1=totals_fy.pop('FY')
    HDR2=totals_fy.pop('-'*max_name_len)
    sorted_totals_fy=sorted(totals_fy.items(),key=lambda x: x[1][year_idx] if x[1][year_idx] != '.' else 0,reverse=True)
    sorted_totals_fy.insert(0,('-'*max_name_len,HDR2))
    sorted_totals_fy.insert(0,('FY',HDR1))
    totals_fy['FY']=HDR1
    totals_fy['-'*max_name_len]=HDR2
    outfile.write('%-*s  %4s\n' %(max_name_len,'NAME','FY'+repr(totals_fy['FY'][year_idx])[2:4]))
    outfile.write('%-*s  %4s\n' %(max_name_len,'-'*max_name_len,'----'))
    for entry in sorted_totals_fy:
        if entry[0]=='N/A':
            continue
        if entry[0][0:2]=='FY' or entry[0][0:2]=='--':
            continue
        name=entry[0]
        total=entry[1][year_idx]
        if total == '.':
            continue
        outfile.write('%-*s  %4s\n' %(max_name_len,name,total))
outfile.close()
#~ print('CREATED TXT FILE \''+startscript+'_totals_fy.txt\'')



#print each year and it's total
outfile=open(startscript+'_counter_totals_cy.txt','w')
outfile.write("\n")
outfile.write("------------------------------\n")
outfile.write("YEARLY TOTALS BY CALENDAR YEAR\n")
outfile.write("------------------------------\n")
for year_idx in range(len(totals['YEAR'])-1):
    outfile.write("\n")
    HDR1=totals.pop('YEAR')
    HDR2=totals.pop('-'*max_name_len)
    sorted_totals=sorted(totals.items(),key=lambda x: x[1][year_idx] if x[1][year_idx] != '.' else 0,reverse=True)
    sorted_totals.insert(0,('-'*max_name_len,HDR2))
    sorted_totals.insert(0,('YEAR',HDR1))
    totals['YEAR']=HDR1
    totals['-'*max_name_len]=HDR2
    outfile.write('%-*s  %4s\n' %(max_name_len,'NAME',repr(totals['YEAR'][year_idx])))
    outfile.write('%-*s  %4s\n' %(max_name_len,'-'*max_name_len,'----'))
    for entry in sorted_totals:
        if entry[0]=='N/A':
            continue
        if entry[0][0:4]=='YEAR' or entry[0][0:2]=='--':
            continue
        name=entry[0]
        total=entry[1][year_idx]
        if total == '.':
            continue
        outfile.write('%-*s  %4s\n' %(max_name_len,name,total))
outfile.close()
#~ print('CREATED TXT FILE \''+startscript+'_totals_cy.txt\'')


#~ print("")
#~ print('SUMMARY:')
#~ infile=open(startscript+'_summary.txt','r')
#~ for line in infile:
    #~ sys.stdout.write(line)
#~ infile.close()
#~ print("")



import zipfile
ZipFile = zipfile.ZipFile(startscript+".zip", "w" )
ZipFile.write(            'counter_aliases.txt',        compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write(startscript+'_counter_aliases_out.txt',   compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write(startscript+'_counter_aliases_sorted.txt',   compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write(startscript+'_counter_dropper_names.txt', compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write(startscript+'_counter_all_drops_tail.txt',     compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write(startscript+'_counter_bad_files.txt',     compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write(startscript+'_counter_duplicates.txt',    compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write(startscript+'_counter_summary.txt',       compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write(startscript+'_counter_totals_fy.txt',     compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write(startscript+'_counter_totals_cy.txt',     compress_type=zipfile.ZIP_DEFLATED)
ZipFile.write(startscript+'_'+sys.argv[0],      compress_type=zipfile.ZIP_DEFLATED)



copyfile(startscript+'_counter_dropper_names.txt', 'counter_dropper_names.txt' )
copyfile(startscript+'_counter_aliases_out.txt',   'counter_aliases_out.txt' )
copyfile(startscript+'_counter_aliases_sorted.txt','counter_aliases_sorted.txt' )
copyfile(startscript+'_counter_all_drops_tail.txt'    , 'counter_all_drops_tail.txt'     )
copyfile(startscript+'_counter_bad_files.txt'    , 'counter_bad_files.txt'     )
copyfile(startscript+'_counter_duplicates.txt'   , 'counter_duplicates.txt'    )
copyfile(startscript+'_counter_summary.txt'      , 'counter_summary.txt'       )
copyfile(startscript+'_counter_totals_fy.txt'    , 'counter_totals_fy.txt'     )
copyfile(startscript+'_counter_totals_cy.txt'    , 'counter_totals_cy.txt'     )


os.remove(startscript+'_counter_aliases_out.txt'   )
os.remove(startscript+'_counter_aliases_sorted.txt')
os.remove(startscript+'_counter_dropper_names.txt' )
os.remove(startscript+'_counter_all_drops_tail.txt'     )
os.remove(startscript+'_counter_bad_files.txt'     )
os.remove(startscript+'_counter_duplicates.txt'    )
os.remove(startscript+'_counter_summary.txt'       )
os.remove(startscript+'_counter_totals_fy.txt'     )
os.remove(startscript+'_counter_totals_cy.txt'     )
os.remove(startscript+'_'+sys.argv[0]  )

print("DONE")



exit()





























#loop through counts and get a list of years seen
for name in counts.keys():
    grand_total[name]=counts[name]['TOTAL']
    for year in counts[name].keys():
        if not year in seen:
            seen.append(year)


#sort seen list and grand_total list
seen=sorted(seen)
grand_total=sorted(grand_total.items(),key=operator.itemgetter(1), reverse=True)






    
    
print("")
print("")
print("####################")
print("# ALL YEARS TOTALS #")
print("####################")
print("")

csvout = open(startscript+'_totals.csv', 'w')
txtout = open(startscript+'_totals.txt', 'w')

sys.stdout.write("%30s" % ("NAME"))
txtout.write("%30s" % ("NAME"))
csvout.write("%s" % ("NAME"))

for year in seen:
    tmpstr="%6s" % (year)
    sys.stdout.write(tmpstr)
    txtout.write(tmpstr)
    csvout.write(",%s" % year)
    
sys.stdout.write("\n")
txtout.write("\n")
csvout.write("\n")

sys.stdout.write("%30s" % ("-----"))
txtout.write("%30s" % ("-----"))
for year in seen:
    sys.stdout.write("%6s" % ("-----"))
    txtout.write("%6s" % ("-----"))
sys.stdout.write("\n")
txtout.write("\n")

for name,tot in grand_total:
    sys.stdout.write("%30s" % (name))
    txtout.write("%30s" % (name))
    csvout.write("%s" % (name))
    for year in seen:
        if year in counts[name]:
            printval=counts[name][year]
            csvprintval=printval
        else:
            printval="."
            csvprintval=""
        sys.stdout.write("%6s" % (printval))
        txtout.write("%6s" % (printval))
        csvout.write(",%s" % (csvprintval))
    sys.stdout.write("\n")
    txtout.write("\n")
    csvout.write("\n")

txtout.write('\nLatest data: '+latest_date+'\n')
txtout.close()
csvout.close()
print("")
print('CREATED TXT FILE \''+startscript+'_totals.txt\'')
print('CREATED CSV FILE \''+startscript+'_totals.csv\'')



#print("")
#print("")
#print("#################")
#print("# YEARLY TOTALS #")
#print("#################")
#print("")
outfile=open(startscript+'_yearly_totals.txt','w')
for year in seen:
    yr_tot={}
    tmpstr="%5s  %s" % (year,"NAME")
    outfile.write(tmpstr+'\n')
    #print(tmpstr)
    tmpstr="%5s  %s" % ("-----","-----")
    outfile.write(tmpstr+'\n')
    #print(tmpstr)
    
    for name in counts.keys():
        if year in counts[name]:
            yr_tot[name]=counts[name][year]
    yr_tot=sorted(yr_tot.items(),key=operator.itemgetter(1), reverse=True)
    for name,tot in yr_tot:
        tmpstr="%5s  %s" % (tot,name)
        #print(tmpstr)
        outfile.write(tmpstr+'\n')
        
    #print('')
    outfile.write('\n')
outfile.write('\nLatest data: '+latest_date+'\n')
outfile.close()
print('CREATED TXT FILE \''+startscript+'_yearly_totals.txt\'')




#print("")
#print("")
#print("#############################")
#print("# TOP DROPPER TITLE HISTORY #")
#print("#############################")
#print("")
outfile=open(startscript+'_top_dropper.txt','w')
tmpstr="%5s  %-19s  %s" % ("TOTAL","TITLE EARNED", "NAME")
#print(tmpstr)
outfile.write(tmpstr+'\n')
tmpstr="%5s  %-19s  %s" % ("-----","------------", "----")
#print(tmpstr)
outfile.write(tmpstr+'\n')

#sort uids (date-order) and add totals for each person chronologically
top=""
uid_tots={}
uids=sorted(uids)
for uid in uids:
    name=uid[17:]
    if name in uid_tots:
        uid_tots[name]+=1
    else:
        uid_tots[name]=1
    if not top in uid_tots:
        top=name
        tmpstr="%5d  %s-%s-%s %s:%s:%s  %s" % (uid_tots[name], uid[1:5], uid[5:7], uid[7:9], uid[10:12],uid[12:14],uid[14:16], name)
        #print(tmpstr)
        outfile.write(tmpstr+'\n')
        
    if uid_tots[name] > uid_tots[top]:
        top=name
        tmpstr="%5d  %s-%s-%s %s:%s:%s  %s" % (uid_tots[name], uid[1:5], uid[5:7], uid[7:9], uid[10:12],uid[12:14],uid[14:16], name)
        #print(tmpstr)
        outfile.write(tmpstr+'\n')
    if name == top:
        last_date="%s-%s-%s %s:%s:%s" % (uid[1:5], uid[5:7], uid[7:9], uid[10:12],uid[12:14],uid[14:16])

tmpstr="%5d  %s  %s" % (uid_tots[top], last_date, top)
#print(tmpstr)
outfile.write(tmpstr+'\n')
outfile.write('\nLatest data: '+latest_date+'\n')
outfile.close()
print('CREATED TXT FILE \''+startscript+'_top_dropper.txt\'')





