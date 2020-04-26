#!/usr/bin/env python

import math
import re
import sys

if len(sys.argv)<2:
    print("""
Usage: %s [ -d ] \
          [ -h hostname ] [ -H header text ] \
          [ -s schemes_used ] [ -w slowdownmax ] 
          filename
""" % sys.argv[0])
    sys.exit(1)
# [ -t tag ] 

####
#### Option handling
####

debug = False
saveToPdf = True # no longer an option
slowdownmax = 6
hostname = "default"
header_text = None
schemes_used = "0,1,2,3,4,5,6,7"
tag = ""
sys.argv = sys.argv[1:]
#print(sys.argv); sys.exit(0)
while len(sys.argv)>1:
    if False:
        continue
    elif sys.argv[0] == "-d":
        debug = True
        sys.argv = sys.argv[1:]
        continue
    elif sys.argv[0]=="-h":
        sys.argv = sys.argv[1:]
        hostname = sys.argv[0]
        sys.argv = sys.argv[1:]
    elif sys.argv[0]=="-H":
        sys.argv = sys.argv[1:]
        header_text = re.sub("_"," ",sys.argv[0])
        sys.argv = sys.argv[1:]
    elif sys.argv[0]=="-s" :
        sys.argv = sys.argv[1:]
        schemes_used = sys.argv[0]
        sys.argv = sys.argv[1:]
    elif sys.argv[0]=="-t":
        sys.argv = sys.argv[1:]
        tag = sys.argv[0]
        sys.argv = sys.argv[1:]
    elif sys.argv[0]=="-w":
        sys.argv = sys.argv[1:]
        slowdownmax = int(sys.argv[0])
        sys.argv = sys.argv[1:]
    else:
        print("Strange option:",sys.argv[0])

import copy
import matplotlib
if saveToPdf:
    matplotlib.use('pdf')
import matplotlib.pyplot as plt

####
#### Data analysis
####

filename = sys.argv[0]
print(".. parsing file <<%s>>" % filename)

def reliable_time(line):
    splits = line.split(); times = []
    for t in splits:
        if re.match(r'[0-9]\.[0-9]+e[+-][0-9]+',t):
            times.append(float(t))
    #print(times)
    avg_time = sum(times)/len(times)
    deviate = sum( [ (t-avg_time)*(t-avg_time) for t in times] )
    sigma = math.sqrt( deviate/(len(times)-1) )
    good_times = [ t for t in times if abs(t-avg_time)<sigma ]
    if len(good_times)==0:
        good_time = -1.
    else:
        good_time = sum(good_times)/len(good_times)
    #print("average: %e becomes %e" % (avg_time,good_time))
    return good_time

schemes_used = [ int(s) for s in schemes_used.split(",") ]
print("Schemes used for plotting: %s" % str(schemes_used))
sizes = []; times = [ [] ]; bws = [ [ ] ]
with open(filename,"r") as timefile:
    for line in timefile:
        sending = re.match("^Sending ([0-9]+) words",line)
        if sending:
            words = sending.groups()[0]; words = int(words)
            sizes.append(8*words)
        scheme = re.search("Scheme ([0-9]+) ",line)
        if scheme:
            mode = scheme.groups()[0]; mode = int(mode)
            if mode>=len(times):
                times.append( [] ); bws.append( [] )
        measured_times = re.match("Measured",line);
        if measured_times:
            time = reliable_time(line)
            if mode in schemes_used:
                times[mode].append(time)
        data_rate = re.search("data rate=(.*)Gb/s",line)
        if data_rate:
            bwidth = data_rate.groups()[0]
            bwidth = float(bwidth)
            if mode in schemes_used:
                bws[mode].append(bwidth)

if debug:
    print("sizes: %d" % len(sizes))
    print("times: %d" % len(times))
    print("    %s" % str( [ len(t) for t in times ] ) )
    print("bws  : %d" % len(bws))

####
#### Plotting
####

scheme_names = [ "reference","copying","vector type","buffered","persistent",
            "subarray","onesided","packing(e)","packing(v)"]

fig,axes = plt.subplots(nrows=2,ncols=2)
timefig = axes[0][0]
bwfig = axes[0][1]
slowp = axes[1][0]
slowv = axes[1][1] # not actually used
slowv.axis('off')

#timefig.xlabel('words')
#timefig.ylabel('message time')
for ischeme in range(len(times)):
    if not ischeme in schemes_used : continue
    time = times[ischeme]
    print("Scheme {} has {} timings; #sizes={}".\
          format(ischeme,len(time),len(sizes)))
    if len(time)==0 or len(time) > len(sizes):
        continue
    print(".. plotting scheme %d" % ischeme)
    timesizes = sizes[:len(time)]
    timefig.loglog( timesizes,time, label=scheme_names[ischeme],color="C%d" % ischeme )
    timefig.text(0.1,0.5,"Time (sec)",transform=timefig.transAxes)

#bwfig.xlabel('words')
#bwfig.ylabel('bandwidth')
lines = []
for ischeme in range(len(times)):
    if not ischeme in schemes_used : continue
    bw = bws[ischeme]
    bwsizes = sizes[:len(bw)]
    if bwsizes==0: continue
    l, = bwfig.semilogx( bwsizes,bw, color="C%d" % ischeme )
    bwfig.text(0.1,0.5,"bwidth (Gb/s)",transform=bwfig.transAxes)
    lines.append(l)

slowp.set_ylim( [0,slowdownmax] )
for ischeme in range(len(times)):
    if not ischeme in schemes_used : continue
    bw = bws[ischeme]
    bwsizes = sizes[:len(bw)]
    if bwsizes==0 or len(bw)==0 or bw[0]==0: continue
    #print("Scheme",ischeme,"bw:",bw)
    try :
        bw = [ bws[0][i] / bw[i] for i in range(len(bw)) ]
    except ZeroDivisionError:
        bw = [ 0 for i in range(len(bw)) ]
    l, = slowp.semilogx( bwsizes,bw, color="C%d" % ischeme )
    slowp.text(0.1,0.5,"slowdown",transform=slowp.transAxes)

#plt.ylim( (0,3.5) )

plt.legend(lines,[ scheme_names[i] for i in schemes_used ],loc='upper left')
filename = hostname
if tag != "":
    filename += "-"+tag
if header_text is None:
    report_name = filename
else:
    report_name = filename + header_text
plt.suptitle('Packing on %s' % report_name)
if saveToPdf:
    pdffile = "%s.pdf" % (filename)
    plt.savefig(pdffile)
    print(".. writing <<%s>" % pdffile)
else:
    plt.show()
