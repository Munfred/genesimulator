# Simple gene simulator
# Does not support bistable promoter, and paused polymerase states 
# Please see `genesimulator_bistable.py` for that

import numpy as np
import matplotlib.pyplot as plt

class Polymerase(object):

    def __init__(self,id):
        self.id = id              # polymerase ID
        self.startpos = 0        # loci where it binded
        self.starttime = 0       # time when it bound
        self.endpos = 0          # where it terminated
        self.endtime = 0         # time of termination
        self.aborted = False     # if terminated though abortion, no RNA is producted
        self.times = []          # times, position and tags matrix [ time, pos, tags]

class Locus (object):
    def __init__(self,k=1.0,r=0.05,bindrate=0.0,terminationrate=0.0,abortionrate=0.0):
        self.occupied = 0                       # if there is a polymerase on the locus, occupied == Polymerase.ID
        self.k = k                              # locus elongation rate
        self.r = r                              # locus backtrack rate
        self.bindrate = bindrate                # locus polymerase binding rate (0 if not a promoter)
        self.terminationrate = terminationrate  # locus termination rate (0 if not a terminator)
        self.abortionrate = abortionrate        # locus abortion rate (0 if not a terminator)


# each gelement is declared like this:
# [element lenght, k,r, bindrate, termination rate,abortionrate]

elements = [
           [1, 1, 0, 0.1, 0, 0],      # 1 locus binding element
           [2, 1, 1, 0, 0, 0],        # 2 loci elongation element
           [1, 1, 0.2, 0, 0, 0],    # 1 loci pause
           [6, 1, 0.2, 0, 0, 0],      # 6 loci elongatio
           [1,  1, 0.2, 0, 0, 0.1],   # 1 loci abortion
           [2, 1, 1, 0, 0, 0],        # 2 loci elongation
           [1, 0, 0.2, 0, 1, 0],      # 1 loci termination
           ]

gene = []   #gene is a list that is initialized by creating Locus in it by reading the elements list
genematrix = [] #matrix containing the rates for each locus in the gene
for line in elements:
    for i in range(line[0]):
        gene.append(Locus(*line[1:6]))
        genematrix.append(line[1:6])

# necessary element so that I can do gene[i-1] for i =0 and gene [i+1] for gene = length. Must be "occupied" so that no pols try to go there
end = Locus(0,0,0,0,0)
end.occupied = True
gene.append(end)


totaltime = 100000    # total simulation time
pols = [0]          # list to store each of the Polymerase info (they're objects)
length = len(gene)  # gene length
z = 0               # normalizer weight number
t = 0               # current time
dt = 0              # timestep
x = 0               # random number used for timesteps
y = 0               # random number used for picking action
picker = 0          # picker = z*y for picking action
w = 0               # used for counting when picking
polcount = 0        # polymerase global counter
plotstep = totaltime/200       # time interval between plots
plotcounter = 0     # timer counter for the plots
plotmatrix =[]      # stores the values to be plotted
totalsteps = 0      #counts number of iterations of simulation
genestate = [None] * ((length)-1)   # gene state matrix for plotting

### ~~~ Main while loop for Gillespie ~~~ ###

while (t < totaltime):                  #while loop of gillespie algorithm
    z = 0                               #z is the normalizing number (but instead of normalizing i multiply the weights by z because its lets costly)
    x = np.random.random_sample()       #draws random number x for choosing time of event
    y = np.random.random_sample()       #draws random number y for choosing which event to happen
    w = 0                               #w is the accumulator used for picking the event to happen
    totalsteps += 1
    if t > plotcounter*plotstep:        #every plot timestep append the gene state to the genestate matix for plotting
        plotcounter = plotcounter +1
        for p in range(length - 1):
            genestate[p] = gene[p].occupied
        plotmatrix.append(genestate[:])
        print ('Time = ' + str(t))

    for i in range (length - 1):        #this loop calculates z by summing the rates of all events that can happen
        if gene[i].occupied != 0 and gene[(i+1)].occupied == 0:
            z += gene[i].k
        if gene[i].occupied != 0 and gene[i - 1].occupied == 0:
            z += gene[i].r
        if gene[i].bindrate and gene[i].occupied == 0:
            z += gene[i].bindrate
        if gene[i].terminationrate and gene[i].occupied !=0:
            z += gene[i].terminationrate
        if gene[i].abortionrate and gene[i].occupied != 0:
            z += gene[i].abortionrate

    dt = np.log(1/x)/z      #calculates the timestep size
    t = t + dt              #increases the time
    picker = z*y            #picker for picking the random number

    for i in range(length - 1): #this for loop keeps going and addings the weight of each event to w until w>picker, then it realizes the event and breaks off

        if gene[i].occupied != 0 and gene[i + 1].occupied == 0: #elongation event
            w += gene[i].k
            if w > picker:
                gene[i+1].occupied = gene[i].occupied   #moves to the location ahead
                gene[i].occupied = 0
                break

        if gene[i].occupied != 0 and gene[i - 1].occupied == 0: #backtrack event
            w += gene[i].r
            if w > picker:
                gene[i-1].occupied = gene[i].occupied   #moves to the location behind
                gene[i].occupied = 0
                break

        if gene[i].bindrate and gene[i].occupied == 0:          #binding event
            w += gene[i].bindrate
            if w > picker:
                polcount = polcount +1;
                pols.append(Polymerase(polcount))   #creates a new polymerase with id = polcount
                pols[polcount].starttime = t        #stores in the polymerase what time it bound
                pols[polcount].startpos = i         #stores in the polymerase where it bound
                gene[i].occupied = polcount         #puts in the binding location the id of the polymerase
                break

        if gene[i].terminationrate and gene[i].occupied != 0:   #termination event
            w += gene[i].terminationrate
            if w > picker:
                pols[gene[i].occupied].endtime = t  #stores in the polymerase the end time
                pols[gene[i].occupied].endpos = i   #stores in the polymerase the end position
                gene[i].occupied = 0                #frees the location
                break

        if gene[i].abortionrate and gene[i].occupied != 0:      #abortion event
            w += gene[i].abortionrate
            if w > picker:
                #print('loc ' + str(i) + ' aborted')
                gene[i].occupied = 0
                break

########### ~~~ Shenanigans for plotting ~~~ ##########

genematrix = list(map(list, zip(*genematrix)))
cmap = plt.cm.jet
cmap.set_under(color='white')
endtimes = []
for i in range(1,polcount):
    if pols[i].endtime:
        endtimes.append(pols[i].endtime)
endtimes = sorted(endtimes)
intervals = np.diff(endtimes)
poltimes = []
for i in range(1,polcount):
    if pols[i].endtime:
        poltimes.append(pols[i].endtime - pols[i].starttime)

plt.figure(0)
ax1 = plt.subplot2grid((3,3), (0,0),rowspan=3)
occplot = plt.imshow(plotmatrix,cmap=cmap,interpolation='nearest',vmin=0.5)
plt.title('Occupancy over time')
plt.ylabel('Time')
ax1.axes.get_xaxis().set_ticks([])

meantravel = np.mean(poltimes)
ax2 = plt.subplot2grid((3,3), (0,1), colspan=2)
plt.hist(poltimes,bins =50, color='#1eefff')
plt.title('Polymerase travel times histogram, mean = ' + str(meantravel))
plt.grid(True)

meaninterval = np.mean(intervals)
ax3 = plt.subplot2grid((3,3), (1, 1), colspan=2)
plt.hist(intervals,bins =50, color='#1ee000')
plt.grid(True)
plt.title('Termination Intervals, mean = ' + str(meaninterval))

ax4 = plt.subplot2grid((3,3), (2, 1), colspan=2)
polrange = range(len(endtimes))
plt.plot(endtimes,polrange, color = 'red')
plt.title('Terminated polymerases over time')
plt.grid(True)

plt.suptitle('Total time = ' + str(totaltime) + 's, Gene length = '+ str(length - 1) + 'L, ' + str(polcount) + ' pols, ' + str(totalsteps) + ' actions \n')
plt.show()

#plt.figure(1)
#genemap = plt.imshow(genematrix,cmap='Greens',interpolation='nearest')
#plt.axis('off')
plt.show()
