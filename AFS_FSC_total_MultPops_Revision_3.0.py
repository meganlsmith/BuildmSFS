#!/usr/bin/env python

"""
This script will build a multidimensional allele frequency spectrum from an input 
matrix of SNPs. For use for as few as three and as many as ten populations. SNP matrix is 
filtered to all bi-allelic SNPs, that equal or surpass the population 
thresholds. The script will sample a single SNP per locus, and if a threshold that requires
subsampling is used, the user can replicate the observed AFS N times 
due to the subsampling of alleles per SNP. Output is an observed allele 
frequency spectrum for use in fastsimcoal2.

For use with pyRAD  output.
                                                     
python AFS_FSC_total_MultPops_Done.py traits.txt SNP_infile.txt 
       Threshold Monomorphics.txt/species.loci Replicate True

"""

import sys
import os
import shutil
import csv
import random
from scipy import ndarray
import numpy as np
from sys import argv
from sortedcontainers import SortedDict
from collections import Counter
from collections import OrderedDict
from operator import itemgetter
import heapq
import collections

script,Traits,file,Threshold,locus_file,nreps, popmonomorphic = argv #arguments the user must supply

Threshold=int(Threshold) # make sure threshold is an integer
nreps=int(nreps) # make sure nreps is an integer
#print Threshold
#print nreps

def pop_association(Traits):
    """Sets the individuals to their respective populations. 
    Also returns the sample counts per population."""
    print 'Processing traits file...'
    with open(Traits, 'r') as traits:
        Pops = OrderedDict() # initializes a dictionary
        Pop_counts = SortedDict() # initializes a dictionary
        next(traits) # iterates to next item in file, skipping the header
        for line in traits: # for each line
            line = line.strip().split() # strip the line based on the tab and split it into two parts (allele and population)
            Pops[line[0]] = line[1] # create a key in the Pops dictionary named for the allele and assign to it the name of the population to which the allele belongs
            if line[1] in Pop_counts: # check to see if the population is already present in the Pop_counts dictionary as a key
                Pop_counts[line[1]] += 1 # if it is, increase the count of alleles for that population
            else: # if it isn't
                Pop_counts[line[1]] = 1 # set the count equal to one
#        print Pops, Pop_counts
        return Pops, Pop_counts # return both of these dictionaries for use in other functions

#Pops, Pop_counts = pop_association(Traits)

def Get_Thresholds(Threshold,Pops,Pop_counts):
    """Get your thresholds organized"""
    print 'Storing thresholds...'
    Thresholds= SortedDict() # initialize a dictionary
    for i in range(0,len(Pop_counts)): # for each population
#        print i
        keys_thresh = '%s' % Pop_counts.iloc[i]# name a key after the population number
        Thresholds.update({keys_thresh:[]}) # add keys to empty dictionaries
        call= Pop_counts[keys_thresh] # call this index from Populations to figure out how many alleles are available for the population
        t = float(Threshold)/100 # convert threshold to proportion
        value = int(call * t) # calculate number of alleles needed for the population
        Thresholds[keys_thresh].append(value) # append value to proper dictionary key
#    print Thresholds
    return Thresholds
#Thresholds = Get_Thresholds(Threshold, Pops, Pop_counts)

#
def createpopcounts(Pops,Pop_counts):
    """This creates a dictionary with populations as keys for counting up to the thresholds."""
    CheckThresholds=Counter()
    for popid in range(0,len(Pop_counts)):
        Keyword = '%s_thresh' % Pop_counts.iloc[popid]
        CheckThresholds[Keyword] += 0
#    print CheckThresholds
    return CheckThresholds
#createpopcounts(Pops,Pop_counts)

#
def Biallelic_SNPs(file,Pops,Pop_counts):
#    """Filter SNP matrix and retain only bi-allelic SNPs that are
#       equal to or above both population thresholds."""
#   
    print "Filtering for Biallelic SNPs that meet thresholds..."    
    RawSNPs = csv.reader(open(file), delimiter = '\t') # open infile
    Data = zip(*RawSNPs) # make infile iterable
    indivs = Data[0] # assign values from column 1 to be names of individuals
    columns = Data[1:] # assign other column values to variable columns
#    print columns
##    #List of the unique polymorphic loci.
    PolyLoci = [] # initialize list to hold polymorphic loci
##
##    #Alleles allowed.
    allowed = ['A','C','G','T'] # specify what characters are allowed as alleles
##
##    #List of bi-allelic SNPs that meet both population thresholds.
    Bi_Thr = [] # initialize list to hold SNPs that meet population thresholds and are biallelic
    keep = 0
    dontkeep = 0
    for i in columns: # iterate through columns
        CheckThresholds = createpopcounts(Pops,Pop_counts)
        letscheck = 0
#        print '\n'
#        print 'Locus %r' % counttoprint        
        if i[0] not in PolyLoci: # If a column isn't in the list of polymorphic loci, add it
            PolyLoci.append(i[0])
        Alleles_set = [] # initialize dictionary to hold alleles that contain allowed character states
        for allele in range(1,len(i)): # for row in columns
            if i[allele] in allowed: # if value is allowed and column is not the first column
                if i[allele] not in Alleles_set: # if i[j] hasn't been added to Alleles_set already
                    Alleles_set.append(i[allele]) # add it
                popid=  Pops[indivs[allele]] # get population
                counterkey = '%s_thresh' % popid # get string to update count
                CheckThresholds.update([counterkey]) # update count for correct population 
        if len(Alleles_set) == 2: ## if the allele is biallelic
            for popid in range(0,len(Pop_counts)): # loop Through the list of populations and make keys with counts of zero for each population
                population= Pop_counts.iloc[popid]
                string = '%s_thresh' % population
#                print CheckThresholds[string]
#                print Thresholds[population][0]
                if CheckThresholds[string] >= Thresholds[population][0]:
#                    print 'in population %s, we have %s alleles, which is more than the threshold, %s' % (population, CheckThresholds[string], Thresholds[population][0])
                    letscheck +=1
#                    print letscheck
        if letscheck == len(Pop_counts):
            Bi_Thr.append(i)
            keep+=1
        elif letscheck < len(Pop_counts): 
            dontkeep+=1
        else:
            print 'something is not working.'
    length = len(columns)
    length_2 = len(PolyLoci)
    return Bi_Thr, indivs, length, length_2

#Bi_Thr, indivs, length, length_2 = Biallelic_SNPs(file,Pops,Pop_counts)

def subsample(Bi_Thr,indivs,length,length_2):
    """Randomly subsample a single linked SNP."""
    print "Subsampling a single linked SNP..."
    prevName = "" # we'll use this string to keep up with whether or not we're still at a certain locus
    Unlink = [] # we'll initialize a list to store our unlinked SNPs
    TempLink = [] # we need to temporarily store linked SNPs in a list, so that we can randomly subsample one
    Num = 1 # this counter will let us know that we're at the beginning of our list and need to assign something to a prevName before moving foreward

    for name in Bi_Thr: # let's look through all items in the list of Biallelic SNPs meeting our thresholds.
    
        GeneName = name[0].strip() # take the first entry in the line and call it the gene name. 
    
        if Num == 1: # if this is the first entry
            prevName = GeneName # assign prevName to the GeneName
            TempLink.append(name) # append this SNP to the TempLink
            Num += 1 # add one to our counter
        elif Num > 1 and Num < len(Bi_Thr): # if this isn't the first or last entry
        
#            print "Middle SNPs: %s\tNum: %d" % (GeneName, Num)
            if GeneName != prevName: # check to see if the GeneName is different than the previous name
            
                #Randomly sample single SNP from TempLink, 
                #and append to Unlink.
                Single = random.choice(TempLink)  # if it is, then we need to randomly choose one SNP from our TempLink
                Unlink.append(Single) # we then add this SNP to the Unlinked list
            
                #Empty list.
                TempLink = [] # we then empty our TempLink list to prepare for the next locus
                prevName = GeneName # we assign our new gene name as the prevName
            
                #Start new Temp list with this SNP.
                TempLink.append(name) # we append this SNP to our nearly cleared list Temp Link
                Num += 1 # and we'll add one to our Num counter to keep track of how many SNPs we've looked at
            
            else: # if the GeneName is the same as prevName
        
                TempLink.append(name) # add this SNP to the TempLink
                prevName = GeneName # set prevName equal to GeneName 
                Num += 1 # add to our counter
   
        elif Num == len(Bi_Thr): # if we're on the last entry of Bi_Thr
        
            if GeneName != prevName: # check to see if this is at the same locus as the previous SNP
                Single = random.choice(TempLink) # if it isn't, randomly choose a SNP from our Temp Link
                Unlink.append(Single) # add this randomly chosen SNP to the Unlink list
                Unlink.append(name) # also add the last entry to the Unlink list
            
            else:
                TempLink.append(name) # if it is the same, append the SNP to the Temp link
                Single = random.choice(TempLink) # randomly subsample one SNP from the Temp link
                Unlink.append(Single) # add this SNP to the Unlink list
    return Unlink

#Unlink = subsample(Bi_Thr,indivs,length,length_2)
#print Unlink
#
def DownSample(Unlink,Bi_Thr, indivs, Pop_counts,Thresholds):
    """Downsample SNPs for which we have more alleles than we need to meet the population thresholds."""
    print "Downsampling SNPs..."
    Downsampled = OrderedDict() # initialize an ordered dictionary to keep track of which individuals we want to sample at each locus.
    Allowed = ['A', 'C', 'G', 'T']
    for snp in Unlink:# loop through the snps
#        print snp
        key = snp[0] # the key is the name of the locus
        Downsampled.update({key:[]}) # add this key to the Downsampled dictionary
        Population_SNPs = Counter() # initialize counter to know how many alleles we have in each population
        for popid in range(0,len(Pop_counts)): # loop Through the list of populations and make keys with counts of zero for each population
            Population_SNPs[Pop_counts.iloc[popid]] += 0
        for population in range(0,len(Pop_counts)): # for each population
            for i in range (1,len(snp)):
#                if snp[i] not in Allowed:
#                    print 'this SNP isn\'t allowed'
                if snp[i] in Allowed:
                    if Pops[indivs[i]] == Pop_counts.iloc[population]: # check, and find individuals that are in that population
                        popper = "%s" % Pops[indivs[i]] 
                        Population_SNPs.update([popper]) # update the SNP count for that population
        for population in range(0,len(Pop_counts)): # for each population
            popper = "%s" % Pop_counts.iloc[population] # set the population name equal to the index
#            print 'In population %s, the number of alleles is %s, and we only need %s' % (popper, Population_SNPs[popper],Thresholds[popper][0])
            if Population_SNPs[popper] == Thresholds[popper][0]: # if this population snp count is equal to the threshold for the population
                for allele in range(1,len(snp)): # then loop through the alleles 
                    if snp[allele] in Allowed: # if the allele is in allowed
                        if Pops[indivs[allele]] == popper: # if the allele is from the current population
                            Downsampled[key].append(indivs[allele]) # add the allele to the downsampled list for that SNP.
            elif Population_SNPs[popper] > Thresholds[popper][0]: # if we exceed the threshold in a population
                counter = 0
                while counter < Thresholds[popper][0]: # While we haven't met the threshold
                    allele =  random.sample(xrange(len(snp)-1),1) # randomly choose a number that corresponds to an allele
                    allele = allele[0] + 1  # add one so we don't try to sample an allele from the locus line
                    if Pops[indivs[allele]] == popper: # if the allele belongs to the correct population
                        if snp[allele] in Allowed: # if the snp[allele] is in the list of allowed characters
                            if indivs[allele] not in Downsampled[key]: # if the individual is not already sampled
                                Downsampled[key].append(indivs[allele]) # add it to alleles to be sampled for the population
#                                print 'added it to %s' % (popper)
                                counter+=1 # use the counter to keep track of how many alleles you've sampled from a given population. 
#    for i in Downsampled:
#        print i, Downsampled[i], len(Downsampled[i])
#        print '\n\n\n'
    return Downsampled
#Downsampled = DownSample(Unlink,Bi_Thr, indivs, Pop_counts,Thresholds)
#print len(Downsampled)
#print Downsampled
#
def Empty_AFS(Pops,Pop_counts):
    print "Creating an empty dictionary for storing the AFS..."
    AFS_Empty = OrderedDict()
    npops = len(Pop_counts) # get count of number of populations
    i = 1 # initialize a counter
    for key in Pop_counts:
#         print i
         if i < npops+1:
             foo = "Max_%r" % i
#             print foo
             i+=1
             count = str(Pop_counts[key] *Threshold/100)
             exec(foo+ "= %s" % count)
    s = 0
    if npops==3:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    bin = '%r_%r_%r' % (q,r,s)
#                    print bin
                    AFS_Empty.update({bin:0})
#                   AFS_Empty[bin]+=1
    if npops==2:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                bin = '%r_%r' % (q,r)
#               print bin
                AFS_Empty.update({bin:0})
#               AFS_Empty[bin]+=1
    if npops==4:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        bin = '%r_%r_%r_%r' % (q,r,s,t)
#                        print bin
                        AFS_Empty.update({bin:0})
#        print AFS_Empty
#        print len(AFS_Empty)
    if npops==5:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            bin = '%r_%r_%r_%r_%r' % (q,r,s,t,u)
#                           print bin
                            AFS_Empty.update({bin:0})
#                           AFS_Empty[bin]+=1
    if npops==6:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            for v in range(Max_6+1):
                                bin = '%r_%r_%r_%r_%r_%r' % (q,r,s,t,u,v)
#                               print bin
                                AFS_Empty.update({bin:0})
#                               AFS_Empty[bin]+=1
    if npops==7:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            for v in range(Max_6+1):
                                for w in range(Max_7+1):
                                    bin = '%r_%r_%r_%r_%r_%r_%r' % (q,r,s,t,u,v,w)
#                                   print bin
                                    AFS_Empty.update({bin:0})
#                                   AFS_Empty[bin]+=1
    if npops==8:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            for v in range(Max_6+1):
                                for w in range(Max_7+1):
                                    for x in range(Max_8+1):
                                        bin = '%r_%r_%r_%r_%r_%r_%r_%r' % (q,r,s,t,u,v,w,x)
#                                       print bin
                                        AFS_Empty.update({bin:0})
#                                       AFS_Empty[bin]+=1
    if npops==9:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            for v in range(Max_6+1):
                                for w in range(Max_7+1):
                                    for x in range(Max_8+1):
                                        for y in range(Max_9+1):
                                            bin = '%r_%r_%r_%r_%r_%r_%r_%r_%r' % (q,r,s,t,u,v,w,x,y)
#                                           print bin
                                            AFS_Empty.update({bin:0})
#                                           AFS_Empty[bin]+=1
    if npops==10:
        for q in range(Max_1+1):
            for r in range(Max_2+1):
                for s in range(Max_3+1):
                    for t in range(Max_4+1):
                        for u in range(Max_5+1):
                            for v in range(Max_6+1):
                                for w in range(Max_7+1):
                                    for x in range(Max_8+1):
                                        for y in range(Max_9+1):
                                            for z in range(Max_10+1):
                                                bin = '%r_%r_%r_%r_%r_%r_%r_%r_%r_%r' % (q,r,s,t,u,v,w,x,y,z)
#                                               print bin
                                                AFS_Empty.update({bin:0})
#                                               AFS_Empty[bin]+=1

#    print len(AFS_Empty)
    return AFS_Empty
#
#AFS_Empty = Empty_AFS(Pops,Pop_counts)
#print AFS_Empty

def create_AFS(Bi_Thr, indivs, length, length_2, Unlink, AFS_Empty, Pop_counts, Thresholds,Downsampled):
    """Count the minor allele in each population 
       and populate the observed AFS."""
    print "Populating the AFS..."
    AFS_Full = AFS_Empty
#    print len(AFS_Full)
    Total_Bi_SNPs_used = 0 # initiate counter
    thehalfSNPs = 0
    Allowed = ['A', 'C', 'G', 'T']
    for snp in Unlink: # loop through snps (columns) in the unlinked group of SNPS
#        print '\n\n\n'
#        print 'results from %r:' % snp[0]
        Allele_Count = Counter() # initialize counter to determine major allele
        PopKeepCount = Counter() # initialize counter to count alleles in different bins
        PopKeepCount_2 = Counter() # initialize counter to count alleles in different bins
        PopKeepCount_3 = Counter() # initialize counter to count alleles in different bins
        for popid in range(0,len(Pop_counts)): # loop Through the list of populations and make keys with counts of zero for each population
            PopKeepCount[Pop_counts.iloc[popid]] == 0
            PopKeepCount_2[Pop_counts.iloc[popid]] == 0
            PopKeepCount_3[Pop_counts.iloc[popid]] == 0
        for allele in range(1,len(snp)): # then loop through the alleles 
            if indivs[allele] in Downsampled[snp[0]]:
                if snp[allele] in Allowed:
                    Allele_Count.update(snp[allele])  # update the SNP count'
                else:
                    print "that's not good"
        #print Allele_Count
        common_allele = Allele_Count.most_common(1)[0][1] # how many occurrences of the most common allele.
        value = Allele_Count.most_common(1)[0][0] # what is the most common allele
#        print(Allele_Count)
        total = sum(Allele_Count.values()) # how many alleles are there
        if int(common_allele) > 0.5 * total:  # is there a common allele
#            print "cool beans"
            for allele in range(1,len(snp)): # if this allele has a frequency > 1/2, loop through the alleles and
#                print '%s in %s' % (snp[allele],popper)
                if snp[allele] != value: # check that the allele is the minor allele
                    popper = Pops[indivs[allele]]
                    if snp[allele] in Allowed:
                        if indivs[allele] in Downsampled[snp[0]]:
                            PopKeepCount.update([popper])
#        else:
#            print "that's okay; it's half"
#        print PopKeepCount
#
#
            index = 1
            for key in Pop_counts:
#                print key
                foo = "store_%s" % str(index)
#                print foo
                x = '%s' % str(PopKeepCount[key])
#                print x
                exec(foo+"=%s" % x)
                index+=1
#                print store_1
#            print len(Pop_counts)
            if len(Pop_counts)==2:
#                print "hells yeah"
                topopulate= '%s_%s' % (store_1, store_2)
                if topopulate != '0_0=':
                     AFS_Full[topopulate]+=1
                     Total_Bi_SNPs_used+=1
            if len(Pop_counts)==3:
#                print "hells yeah"
                topopulate= '%s_%s_%s' % (store_1, store_2, store_3)
                if topopulate != '0_0_0':
                     AFS_Full[topopulate]+=1
                     Total_Bi_SNPs_used+=1
            if len(Pop_counts)==4:
                topopulate= '%s_%s_%s_%s' % (store_1, store_2, store_3,store_4)
                if topopulate != '0_0_0_0':
                     AFS_Full[topopulate]+=1
                     Total_Bi_SNPs_used+=1
            if len(Pop_counts)==5:
                topopulate= '%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5)
                if topopulate != '0_0_0_0_0':
                     AFS_Full[topopulate]+=1
                     Total_Bi_SNPs_used+=1
            if len(Pop_counts)==6:
                topopulate= '%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6)
                if topopulate != '0_0_0_0_0_0':
                     AFS_Full[topopulate]+=1
                     Total_Bi_SNPs_used+=1
            if len(Pop_counts)==7:
                topopulate= '%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7)
                if topopulate != '0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=1
                     Total_Bi_SNPs_used+=1
            if len(Pop_counts)==8:
                topopulate= '%s_%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7,store_8)
                if topopulate != '0_0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=1
                     Total_Bi_SNPs_used+=1
            if len(Pop_counts)==9:
                topopulate= '%s_%s_%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7,store_8,store_9,store_10)
                if topopulate != '0_0_0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=1
                     Total_Bi_SNPs_used+=1
            if len(Pop_counts)==10:
                topopulate= '%s_%s_%s_%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7,store_8,store_9,store_10)
                if topopulate != '0_0_0_0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=1
                     Total_Bi_SNPs_used+=1
        elif int(common_allele) == 0.5 * total:  # is there a common allele
            #print(int(common_allele))
            #print "we have a non minor allele"
#            store = list(Allele_Count.elements())
#            for item in store:
#                try:
#                    thefirstsnp
#                except NameError:
#                    thefirstsnp = item
##                    print thefirstsnp
#                else:
#                    if item == thefirstsnp:
#                        thefirstsnp=thefirstsnp
#                    else: 
#                        try:
#                            thesecondsnp
#                        except NameError:
#                            thesecondsnp = item
##                            print thesecondsnp
#                        else:
#                            thesecondsnp=thesecondsnp
            thefirstsnp = value
            thesecondsnp = Allele_Count.most_common(2)[1][0]
#            print common_allele
            for allele in range(1,len(snp)): # if this allele has a frequency > 1/2, loop through the alleles and
#                print '%s in %s' % (snp[allele],popper)
                #print(snp[allele])
                if snp[allele] == thefirstsnp: # check that the allele is the minor allele
                    popper = Pops[indivs[allele]]
                    if snp[allele] in Allowed:
                        if indivs[allele] in Downsampled[snp[0]]:
                            PopKeepCount_2.update([popper])
                if snp[allele] == thesecondsnp: # check that the allele is the minor allele
                    popper = Pops[indivs[allele]]
                    if snp[allele] in Allowed:
                        if indivs[allele] in Downsampled[snp[0]]:
                            PopKeepCount_3.update([popper])
            #print(PopKeepCount_3)
            #print(PopKeepCount_2)
            index = 1
            for key in Pop_counts:
                foo = "store_%s" % str(index)
#                print foo
                x = '%s' % str(PopKeepCount_2[key])
#                print x
                exec(foo+"=%s" % x)
                index+=1
#                print store_1
#            print len(Pop_counts)
            if len(Pop_counts)==2:
#                print "hells yeah"
                topopulate= '%s_%s' % (store_1, store_2)
#                if topopulate != '0_0_0':
                AFS_Full[topopulate]+=0.5
                thehalfSNPs+=0.5
            if len(Pop_counts)==3:
#                print "hells yeah"
                topopulate= '%s_%s_%s' % (store_1, store_2, store_3)
#                if topopulate != '0_0_0':
                AFS_Full[topopulate]+=0.5
                thehalfSNPs+=0.5
            if len(Pop_counts)==4:
                topopulate= '%s_%s_%s_%s' % (store_1, store_2, store_3,store_4)
                #print(topopulate)
                if topopulate != '0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==5:
                topopulate= '%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5)
                if topopulate != '0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==6:
                topopulate= '%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6)
                if topopulate != '0_0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==7:
                topopulate= '%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7)
                if topopulate != '0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==8:
                topopulate= '%s_%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7,store_8)
                if topopulate != '0_0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==9:
                topopulate= '%s_%s_%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7,store_8,store_9,store_10)
                if topopulate != '0_0_0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==10:
                topopulate= '%s_%s_%s_%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7,store_8,store_9,store_10)
                if topopulate != '0_0_0_0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            #print "The half count is %r" % thehalfSNPs
            index = 1
            for key in Pop_counts:
#                print key
                foo = "store_%s" % str(index)
#                print foo
                x = '%s' % str(PopKeepCount_3[key])
#                print x
                exec(foo+"=%s" % x)
                index+=1
#                print store_1
#            print len(Pop_counts)
            if len(Pop_counts)==2:
#                print "hells yeah"
                topopulate= '%s_%s' % (store_1, store_2)
#                if topopulate != '0_0_0':
                AFS_Full[topopulate]+=0.5
                thehalfSNPs+=0.5
            if len(Pop_counts)==3:
#                print "hells yeah"
                topopulate= '%s_%s_%s' % (store_1, store_2, store_3)
#                if topopulate != '0_0_0':
                AFS_Full[topopulate]+=0.5
                thehalfSNPs+=0.5
            if len(Pop_counts)==4:
                topopulate= '%s_%s_%s_%s' % (store_1, store_2, store_3,store_4)
                #print(topopulate)
                if topopulate != '0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==5:
                topopulate= '%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5)
                if topopulate != '0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==6:
                topopulate= '%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6)
                if topopulate != '0_0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==7:
                topopulate= '%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7)
                if topopulate != '0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==8:
                topopulate= '%s_%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7,store_8)
                if topopulate != '0_0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==9:
                topopulate= '%s_%s_%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7,store_8,store_9,store_10)
                if topopulate != '0_0_0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
            if len(Pop_counts)==10:
                topopulate= '%s_%s_%s_%s_%s_%s_%s_%s_%s_%s' % (store_1, store_2, store_3,store_4,store_5,store_6,store_7,store_8,store_9,store_10)
                if topopulate != '0_0_0_0_0_0_0_0_0_0':
                     AFS_Full[topopulate]+=0.5
                     thehalfSNPs+=0.5
#            print "The full count is %r" % thehalfSNPs

#    print thehalfSNPs
#    print Total_Bi_SNPs_used
    Total_Bi_SNPs_used = Total_Bi_SNPs_used + thehalfSNPs
#    print len(AFS_Full)
    return AFS_Full,Total_Bi_SNPs_used
#
#create_AFS(Bi_Thr, indivs, length, length_2, Unlink, AFS_Empty, Pop_counts, Thresholds,Downsampled)
def totalbp(locus_file):
    """Count total number of sequenced base pairs"""
    print "Counting sequenced basepairs..."
    
    with open(locus_file, 'r') as Infile:
        loci = ''
        Count = 0
        Loci_count = 0
        start = False
        for line in Infile:
            #count first locus
            if Loci_count == 0:
                line = line.strip().split()
                Count += len(line[1])+6
                Loci_count += 1
            
            #find line breaks
            if '//' in line:
                start = True
            else:
                if start == True and line.startswith(">"):
                    line = line.strip().split()
                    #6 takes into account RE site
                    Count += len(line[1])+6
                    Loci_count += 1
                    start = False                
    return Count, Loci_count

def get_mono_cell(locus_file,Count,Loci_count,TotalBi_SNPs_used, TotalSNPs):
    """Determine value to add to [0,0] cell"""
    print "Finding the value for the monomorphic cell..."
    TotalBP = Count
    if popmonomorphic == True:
        Monomorphics = int((TotalBi_SNPs_used * TotalBP) / TotalSNPs) - TotalBi_SNPs_used
    else:
        Monomorphics = 0
    return Monomorphics, \
           TotalBP, Loci_count

def Add_MONOmorphs(locus_file,Monomorphics, AFS_Full):
    print " Populating the monomorphic cell..."
    mono = int(Monomorphics)
#    print AFS_Full
    key = AFS_Full.keys()[0]
    AFS_Full[key] = mono
    return AFS_Full

def write_logfile(rep, Count, Loci_count, TotalBi_SNPs_used, TotalSNPs, Monomorphics):
    print "Writing the log file..."
    with open("Rep" + str(rep) + "_log.txt", 'w') as log: 
        log.write('Subsampling to populate the SFS')
        log.write('\n\n')
        log.write('Threshold: %r%s' % (int(Threshold),'%')) 
        log.write('\n\n\n')
        for popnum in range(0,len(Pop_counts)):
            log.write('Population %s: %s' % (popnum+1,Pop_counts.iloc[popnum]))
            log.write('\n')
            log.write('%r individuals' % Pop_counts[Pop_counts.iloc[popnum]])
            log.write('\n')
            counttocomma = 0
            for j in Pops:
                if Pops[j] == Pop_counts.iloc[popnum]:
                    if counttocomma == 0: 
                        log.write(j)
                        counttocomma+=1
                    elif counttocomma != 0: 
                        log.write(', ')
                        log.write(j)
            log.write('\n')
            log.write('Total Alleles: %s' % str(Pop_counts[Pop_counts.iloc[popnum]]))
            log.write('\t')
            log.write('Total Alleles Sampled: %s' % str(Pop_counts[Pop_counts.iloc[popnum]]* Threshold/100))
            log.write('\n\n\n')
        log.write('Total SNPs: %s' % TotalSNPs)
        log.write('\t')
        log.write('Total SNPs Used: %s' % TotalBi_SNPs_used)
        log.write('\n')
#        print TotalBi_SNPs_used
        log.write('Total BP: %s' % TotalBP)
        log.write('\t')
        AdjustedBP = int(Monomorphics) + int(TotalBi_SNPs_used)
        log.write('Adjusted BP: %s' % (AdjustedBP))

for rep in range(0,nreps):
    Pops, Pop_counts = pop_association(Traits)
    Thresholds = Get_Thresholds(Threshold,Pops,Pop_counts)
    Bi_Thr, indivs, TotalSNPs, length_2 = Biallelic_SNPs(file,Pops,Pop_counts)
    Unlink = subsample(Bi_Thr,indivs,TotalSNPs,length_2)
    Downsampled = DownSample(Unlink,Bi_Thr, indivs, Pop_counts,Thresholds)
#    print Pops
#    print Pop_counts
    AFS_Empty = Empty_AFS(Pops,Pop_counts)
#    print AFS_Empty
    AFS_Full,TotalBi_SNPs_used = create_AFS(Bi_Thr, indivs, TotalSNPs, length_2, Unlink, AFS_Empty,Pop_counts, Thresholds,Downsampled)
#    print AFS_Full
    Count, Loci_count = totalbp(locus_file)
    Monomorphics,TotalBP, Loci_count = get_mono_cell(locus_file,Count,Loci_count,TotalBi_SNPs_used, TotalSNPs)
    AFS_Full_Mono = Add_MONOmorphs(locus_file,Monomorphics, AFS_Full)
    with open("Rep" + str(rep) + "_MSFS.obs", 'w') as f:
#        print Pops
#        print Pop_counts
        write_logfile(rep, Count, Loci_count, TotalBi_SNPs_used, TotalSNPs, Monomorphics)
        with open ("Rep" + str(rep) + "_MSFS.obs", 'a') as f: 
            f.write('1 observations. No. of demes and sample sizes are on next line')
            f.write('\n')
            f.write(str(len(Pop_counts)))
            f.write('\t')
            for i in range(0,len(Pop_counts)):
                key = '%s' % Pop_counts.iloc[i]
                f.write(str(Pop_counts[key]))
                f.write('\t')
            f.write('\n')
        for bin in AFS_Full_Mono:
            with open("Rep" + str(rep) + "_MSFS.obs", 'a') as f:
                f.write(str(AFS_Full_Mono[bin]))
                f.write('\t')

