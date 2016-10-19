import csv
import numpy as np
import time
start_time = time.time()
import sys
import random
rng = random.SystemRandom()

class Polymerase(object):

    def __init__(self,id):
        self.id = id             # polymerase ID
        self.startpos = 0        # loci where it binded
        self.starttime = 0       # time when it bound
        self.endpos = 0          # where it terminated
        self.endtime = 0         # time of termination
        self.aborted = False     # if terminated though abortion, no RNA is producted
        self.times = []          # times, position and tags matrix [ time, pos, tags]
        self.paused = False      # keeps track is pol ever goes into paused state

class Locus (object):
    def __init__(self,k=80.0,r=0,bindrate=0.0,terminationrate=0.0,abortionrate=0.0):
        self.state = 'ON'                       # state of that loci (ON/OFF is for promoters)
        self.occupied = 0                       # if there is a polymerase on the locus, occupied == Polymerase.ID
        self.k = k                              # locus elongation rate
        self.r = r                              # locus backtrack rate
        self.bindrate = bindrate                # locus polymerase binding rate (0 if not a promoter)
        self.terminationrate = terminationrate  # locus termination rate (0 if not a terminator)
        self.abortionrate = abortionrate        # locus abortion rate (0 if not a terminator)


# each gelement is declared like this:
# [element lenght, k,r, bindrate, termination rate,abortionrate]
B_ON = 0.5 #binding rate when promoter is ON
B_OFF = 0 #binding rate when the promoter is OFF
footprint = 30 #polymerase footprint in base pairs (loci)
elements = [
           [1, 80, 0, B_ON, 0, 0],      # 1 locus binding element
           [5000, 80, 0, 0, 0, 0],      # 2000 loci elongation element
           [footprint + 1, 0, 0, 0, 80, 0],         # 30 loci termination element...to make sure the footprint clearance doesnt prevent pols form terminating
           ]
enter_pause_rate = 0.01 #rate at which polymerase goes into paused state from any locus
exit_pause_rate = 4 #rate at which polymerase exits the pause state

enter_ON_state = 0.5 #rate at which promoter off (polymerase can bind at rate B) if it was ion the OFF state
enter_OFF_state = 0.2 # rate at which the promoter goes into the OFF state

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


totaltime = 70000    # total simulation time
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

oldt = 0 # variable for keeping track of when to take snapshots
snapshotstep = 60 #interval between snapshots

### ~~~ Main while loop for Gillespie ~~~ ###
maxsteps = 20000
# open a file for writing.
randname = np.random.random_sample() # random string for filename


csv_out = open("snapshots_" + str(randname) + ".csv", 'wb')
mywriter = csv.writer(csv_out)

while (t< totaltime):                  #while loop of gillespie algorithm
    if totalsteps % 100 == 0:
        print (str(totalsteps) + ' Steps,' + ' Time  = ' + str(t))
    z = 0                               #z is the normalizing number (but instead of normalizing i multiply the weights by z because its lets costly)
    x = rng.random()       #draws random number x for choosing time of event
    y = rng.random()      #draws random number y for choosing which event to happen
    w = 0                               #w is the accumulator used for picking the event to happen
    totalsteps += 1
    if (t-oldt > snapshotstep and t > 500):        #every plot timestep append the gene state to the genestate matrix for plotting
        plotcounter = plotcounter +1
        for p in range(length - 1):
            genestate[p] = gene[p].occupied
        plotmatrix.append(genestate[:])

        mywriter.writerows([genestate]) # writes out the new row of genematrix to the csv file
        oldt = t #reserts our interval counter

        print "totalsteps =" + str(totalsteps)
        print "simulated time =" + str(t)
        print("--- %s seconds ---" % (time.time() - start_time))
        print "Total pols = " + str(polcount)

        endtimes = []
        for i in range(1, polcount):
            if pols[i].endtime:
                endtimes.append(pols[i].endtime)
        endtimes = sorted(endtimes)
        intervals = np.diff(endtimes)
        poltimes = []
        for i in range(1, polcount):
            if pols[i].endtime:
                poltimes.append(pols[i].endtime - pols[i].starttime)

        meantravel = np.mean(poltimes)
        meaninterval = np.mean(intervals)
        print "Polymerase mean travel time = " + str(meantravel)
        print 'Termination Intervals, mean = ' + str(meaninterval)

        text_file = open("Info_" + str(randname) + ".txt", "w")
        text_file.write("\n Gene length = " + str(length - footprint - 1))
        text_file.write("\n Totalsteps =" + str(totalsteps))
        text_file.write("\n Total polymerases =" + str(polcount))
        text_file.write("\n Simulated time =" + str(t))
        text_file.write('\n Mean Termination Intervals = ' + str(meaninterval))
        text_file.write("\n Polymerase mean travel time = " + str(meantravel))
        text_file.write("\n TERMINATION INTERVALS \n")
        text_file.write(str(intervals))
        text_file.write("\n TRAVEL TIMES \n")
        text_file.write(str(poltimes))
        text_file.write("\n TERMINATION TIMES \n")
        text_file.write(str(endtimes))
        text_file.close()


    if gene[0].state == 'ON':
        z += enter_OFF_state
    if gene[0].state == 'OFF':
        z += enter_ON_state

    for i in range (length - 1):        #this loop calculates z by summing the rates of all events that can happen

        if gene[i].occupied != 0 and (not pols[gene[i].occupied].paused) and (not any([gene[j].occupied for j in range(i+1, i + footprint+1)])):
        #if a loci is occupied and pol is not paused and the next 30 ones (footprint size)are free, consider elongation
            z += gene[i].k

        if gene[i].occupied != 0 and (not pols[gene[i].occupied].paused):
        # if a loci is occupied and pol is not paused, consider going into paused state
            z += enter_pause_rate

        if gene[i].occupied != 0 and (pols[gene[i].occupied].paused):
        # if a loci is occupied and pol is paused, consider exiting pause
            z += exit_pause_rate

        #BACKTRACK OPTION IS COMMENTED OUT BECAUSE WE DON'T USE IT HERE, SO MAKES THINGS A BIT FASTER
        #if gene[i].occupied != 0 and i>(footprint + 1) and (not pols[gene[i].occupied].paused) and (not all([gene[j].occupied for j in range(i - footprint, i+1)])):
        #if a loci is occupied and pol is not paused and the previous 30 ones are free, consider backtrack
        #    z += gene[i].r

        if gene[i].bindrate and (not any([gene[j].occupied for j in range(i, i + footprint)])):
        #if a loci  and the 30 ahead are is free, consider binding
            z += gene[i].bindrate

        if gene[i].terminationrate and gene[i].occupied !=0:
            #if a loci is occupied, consider termination
            z += gene[i].terminationrate

        # ABORTION OPTION IS COMMENTED OUT BECAUSE WE DON'T USE IT HERE, SO MAKES THINGS A BIT FASTER
        #if gene[i].abortionrate and gene[i].occupied != 0: #if a loci is occupied, consider abortion
        #    z += gene[i].abortionrate

            ''''            ~~~~~~~~ ROUTINE FOR PICKING WHICH EVENT ACTUALLY HAPPENS ~~~~~~~~~~~~          '''

    dt = np.log(1/x)/z      #calculates the timestep size
    t = t + dt              #increases the time
    picker = z*y            #picker for picking the random number


    # at this point w = 0 and since we only have one ON/OFF promoter, we can keep the switching routine outside the for loop

    if w < picker and gene[0].state == 'ON':
        w += enter_OFF_state
        if w > picker:
            gene[0].bindrate = B_OFF
            gene[0].state = 'OFF'
            #print 'turned ON time =' + str(t)


    if w < picker and  gene[0].state == 'OFF':
        w += enter_ON_state
        if w > picker:
            gene[0].bindrate = B_ON
            gene[0].state = 'ON'
            #print 'turned OFF time =' + str(t)


            ''''        ~~~~    for loop that executes whatever action       ~~~~        '''

    for i in range(length - 1): #this for loop keeps going and addings the weight of each event to w until w>picker, then it realizes the event and breaks off

        if w < picker and gene[i].occupied != 0 and (not pols[gene[i].occupied].paused) and (not any([gene[j].occupied for j in range(i+1, i + footprint+1)])):
        #elongation event: if a loci is occupied and pol is not paused and the next 30 ones (footprint size)are free
            w += gene[i].k
            if w > picker:
                gene[i+1].occupied = gene[i].occupied   #moves to the location ahead
                gene[i].occupied = 0
                break


        if w < picker and gene[i].occupied != 0 and (not pols[gene[i].occupied].paused):
        # goes into paused state - if a loci is occupied and pol is not paused
            #print 'picker = ' + str(picker) + '+ w =' + str(w) + 'poop'

            w += enter_pause_rate
            #print picker
            if w > picker:
                pols[gene[i].occupied].paused = True
                #print gene[i].occupied
                #print 'paused'
                break


        if w < picker and gene[i].occupied != 0 and (pols[gene[i].occupied].paused):
        # if a loci is occupied and pol is paused, consider exiting pause
            w += exit_pause_rate
            if w > picker:
                pols[gene[i].occupied].paused = False
                #print gene[i].occupied
                #print 'UNpaused'
                break

        # BACKTRACK OPTION IS COMMENTED OUT BECAUSE WE DON'T USE IT HERE, SO MAKES THINGS A BIT FASTER
        #if w < picker and gene[i].occupied != 0 and gene[i - 1].occupied == 0: #backtrack event
        #    w += gene[i].r
        #    if w > picker:
        #        gene[i-1].occupied = gene[i].occupied   #moves to the location behind
        #        gene[i].occupied = 0
        #        break


        if w < picker and gene[i].bindrate and (not any([gene[j].occupied for j in range(i, i + footprint)])):
        #if a loci  and the 30 ahead are is free, consider binding
            w += gene[i].bindrate
            if w > picker:
                polcount = polcount +1;
                pols.append(Polymerase(polcount))   #creates a new polymerase with id = polcount
                pols[polcount].starttime = t        #stores in the polymerase what time it bound
                pols[polcount].startpos = i         #stores in the polymerase where it bound
                gene[i].occupied = polcount         #puts in the binding location the id of the polymerase
                break

        if w < picker and gene[i].terminationrate and gene[i].occupied != 0:   #termination event
            w += gene[i].terminationrate
            if w > picker:
                pols[gene[i].occupied].endtime = t  #stores in the polymerase the end time
                pols[gene[i].occupied].endpos = i   #stores in the polymerase the end position
                gene[i].occupied = 0                #frees the location
                break

        # ABORTION OPTION IS COMMENTED OUT BECAUSE WE DON'T USE IT HERE, SO MAKES THINGS A BIT FASTER
        #if w < picker and gene[i].abortionrate and gene[i].occupied != 0:      #abortion event
        #    w += gene[i].abortionrate
        #    if w > picker:
        #        #print('loc ' + str(i) + ' aborted')
        #        gene[i].occupied = 0
        #        break


########### ~~~ Shenanigans for writing results ~~~ ##########





