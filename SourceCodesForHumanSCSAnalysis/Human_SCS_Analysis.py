#!/usr/bin/env python
# coding: utf-8


import numpy as np
from matplotlib import pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

#aminoAcidCode = {'A':0, 'R':1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19, 'U': 20, 'X':21}
#aminoAcidCodeInverse = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U', 'X']
aminoAcidCode = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R':14, 'S': 15, 'T': 16, 'U': 17, 'V': 18, 'W': 19, 'X': 20, 'Y': 21}
aminoAcidCodeInverse = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y']
googleDrivePath = '/content/drive/My Drive/ncbi_dataset/'
googleColab = False

if googleColab == True:
    dataPath = googleDrivePath
else:
    dataPath = './ncbi_dataset/'

if googleColab == True:
    figPath = googleDrivePath + 'Figs/'
else:
    figPath = './ncbi_dataset/Figs/'

NumAminoAcids = 22
proteinData = []
proteinMetaData = []
aminoAcidFrequency = np.zeros(NumAminoAcids)
TotalAminoAcids = 0
TotalSCSs = [0, 0, 0, 0, 0, 0, 0, 0]
scsLength3 = 3
scsLength4 = 4
scsLength5 = 5
BigArraySize3 = NumAminoAcids**scsLength3
scsFrequency3 = np.zeros(BigArraySize3)
scsAvailability3 = np.zeros(BigArraySize3)
sortedScsFrequency3 = np.empty([0,0])
sortedScsAvailability3 = np.empty([0,0])
sortedScsAvailabilityUX3 = np.empty([0,0])

BigArraySize4 = NumAminoAcids**scsLength4
scsFrequency4 = np.zeros(BigArraySize4)
scsAvailability4 = np.zeros(BigArraySize4)
sortedScsFrequency4 = np.empty([0,0])
sortedScsAvailability4 = np.empty([0,0])
sortedScsAvailabilityUX4 = np.empty([0,0])
BigArraySize5 = NumAminoAcids**scsLength5
scsFrequency5 = np.zeros(BigArraySize5)
scsAvailability5 = np.zeros(BigArraySize5)
sortedScsFrequency5 = np.empty([0,0])
sortedScsAvailability5 = np.empty([0,0])
sortedScsAvailabilityUX5 = np.empty([0,0])

def scs2code(scs, scsLength):
  code = 0
  for i in range(len(scs)):
    code += aminoAcidCode[scs[i]] * NumAminoAcids**(scsLength - i - 1)
  return code

def code2scs(code, scsLength):
  value = code
  remainder = code
  stack = []
  scs = ''
  for i in range(scsLength):
    remainder = int(value % NumAminoAcids)
    value = value // NumAminoAcids
    stack.append(aminoAcidCodeInverse[remainder])
  while len(stack) != 0:
    scs += stack.pop()
  return scs




def showBasicInformation():
    print('The data was downloaded on Aug. 31, 2020.')
    print('The number of proteins included in the dataset:', '{:,}'.format(len(proteinData)))
    print('The total amino acids:', '{:,}'.format(TotalAminoAcids))
    print('The number of amino acids:')
    for i in range(NumAminoAcids):
        print('  ', aminoAcidCodeInverse[i], ':', '{:,}'.format(int(aminoAcidFrequency[i])))
    print('The number of SCSs of length 3:', '{:,}'.format(TotalSCSs[scsLength3]))
    print('The number of SCSs of length 4:', '{:,}'.format(TotalSCSs[scsLength4]))
    print('The number of SCSs of length 5:', '{:,}'.format(TotalSCSs[scsLength5]))



def readProteinDatafile():

    global dataPath
 
    #f = open('./ncbi_dataset/data/GCF_000001405.39/protein.faa')
    f = open(dataPath + 'protein.faa')
    #f = open('./proteinData/protein.faa')
    #readData = f.readlines()

    proteinMetaData.append(f.readline())
    aminoSeq = ''
    for line in f:
        if line[0] == '>':
            proteinData.append(aminoSeq)
            aminoSeq = ''
            proteinMetaData.append(line)       
        else:
            aminoSeq += line.rstrip()        
    proteinData.append(aminoSeq)

    
def objectDump():
    import joblib

    global TotalSCSs
    global dataPath

    joblib.dump(TotalSCSs, dataPath + "TotalSCSs.jb", compress=5)
    joblib.dump(aminoAcidFrequency, dataPath + "aminoAcidFrequency.jb", compress=5)
    joblib.dump(scsFrequency3, dataPath + "scsFrequency3.jb", compress=5)
    joblib.dump(scsFrequency4, dataPath + "scsFrequency4.jb", compress=5)
    joblib.dump(scsFrequency5, dataPath + "scsFrequency5.jb", compress=5)
    joblib.dump(scsAvailability3, dataPath + "scsAvailability3.jb", compress=5)
    joblib.dump(scsAvailability4, dataPath + "scsAvailability4.jb", compress=5)
    joblib.dump(scsAvailability5, dataPath + "scsAvailability5.jb", compress=5)
 
    
def objectLoad():
    import joblib
    global TotalAminoAcids
    global TotalSCSs
    global aminoAcidFrequency
    global scsFrequency3
    global scsFrequency4
    global scsFrequency5
    global scsAvailability3
    global scsAvailability4
    global scsAvailability5
    global dataPath

    TotalSCSs = joblib.load(dataPath + "TotalSCSs.jb")
    aminoAcidFrequency = joblib.load(dataPath + "aminoAcidFrequency.jb")
    scsFrequency3 = joblib.load(dataPath + "scsFrequency3.jb")
    scsFrequency4 = joblib.load(dataPath + "scsFrequency4.jb")
    scsFrequency5 = joblib.load(dataPath + "scsFrequency5.jb")
    scsAvailability3 = joblib.load(dataPath + "scsAvailability3.jb")
    scsAvailability4 = joblib.load(dataPath + "scsAvailability4.jb")
    scsAvailability5 = joblib.load(dataPath + "scsAvailability5.jb")
    
    makeSortedDataset()
    TotalAminoAcids = sum(aminoAcidFrequency)


def showProteinLengthDistribution(until = 0, save = 'n'):
    import statistics
    import math
    import matplotlib 
    get_ipython().run_line_magic('matplotlib', 'inline')
    import matplotlib.pyplot as plt
    from pylab import rcParams

    global figPath

    fig = plt.figure()
    rcParams['figure.figsize'] = 20,10
    matplotlib.rc('xtick', labelsize=12) 
    matplotlib.rc('ytick', labelsize=12)

    proteinLength = []
    for protein in proteinData:
        proteinLength.append(len(protein))
    #print('Sum:', '{:,}'.format(sum(proteinLength))))
    print('Max:', '{:,}'.format(max(proteinLength)))
    print('Min:', '{:,}'.format(min(proteinLength)))       
    print('Mean:', '{:5.1f}'.format(statistics.mean(proteinLength)))
    print('Median:', '{:5.1f}'.format(statistics.median(proteinLength)))
    print('Stdev:', '{:5.1f}'.format(statistics.stdev(proteinLength)))
    
    if until == 0:
        plt.hist(proteinLength, bins=1000)
    elif until < 10000:
        plt.hist(proteinLength, range = (0, until), bins=int(until/10))
    else:
        plt.hist(proteinLength, range = (0, until), bins=1000)

    plt.title('Distribution of Protein Length (Homo Sapiens)', fontsize = 12)
    plt.xlabel('Length of Protein', fontsize=12)
    plt.ylabel('Number of Proteins', fontsize=12)
    #plt.xticks(fontsize=14)
    #plt.yticks(fontsize=14)
    plt.grid(True)
    if save == 'y':
        import datetime
        fig.savefig(figPath + "protein_length_dist" + str(datetime.datetime.now().microsecond) + ".png", format="png", dpi=600)

    plt.show
    
    
    
def showScsFrequencyDistribution(scsLength, save = 'n'):
    import statistics
    import math
    import matplotlib 
    get_ipython().run_line_magic('matplotlib', 'inline')
    import matplotlib.pyplot as plt
    from pylab import rcParams
    rcParams['figure.figsize'] = 20,10
    matplotlib.rc('xtick', labelsize=20) 
    matplotlib.rc('ytick', labelsize=20)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    plt.subplots_adjust(wspace=0.4, hspace=0.6)



    if scsLength == scsLength3:
        scsData = scsFrequency3
    elif scsLength == scsLength4:
        scsData = scsFrequency4
    else:
        scsData = scsFrequency5
    
    print('Frequency of SCSs (length', '{:1d}'.format(scsLength), ')')
    print('Max:', '{:,}'.format(int(max(scsData))))
    print('Min:', '{:,}'.format(int(min(scsData))))       
    print('Mean:', '{:5.10f}'.format(statistics.mean(scsData)))
    print('Median:', '{:5.10f}'.format(statistics.median(scsData)))
    print('Stdev:', '{:5.10f}'.format(statistics.stdev(scsData)))

    ax1.hist(scsData, range=(0, 1000), bins=200)
    #ax1.gca().set_yscale("log") 
    ax1.set_yscale("log") 


    ax1.set_title('Distribution of SCS Frequency (0 to 1000)', fontsize = 12)
    ax1.set_xlabel('SCS Frequency', fontsize = 12)
    ax1.set_ylabel('Number of SCSs', fontsize = 12)
    ax1.grid(True)
    #plt.show
      
    ax2.hist(scsData, bins=np.logspace(np.log10(0.0000001), np.log10(max(scsData)), 200))
    #ax2.hist(scsData, bins = 10000)
    ax2.set_xscale("log") 
    ax2.set_yscale("log") 
    ax2.set_title('Distribution of SCS Frequency (0 to max)', fontsize = 12)
    ax2.set_xlabel('SCS Frequency', fontsize = 12)
    ax2.set_ylabel('Number of SCSs', fontsize = 12)
    ax2.grid(True)

    if save == 'y':
        import datetime
        fig.savefig(figPath + "freq_dist" + str(datetime.datetime.now().microsecond) + ".png", format="png", dpi=600)


    plt.show()




def showScsAvailabilityDistribution(scsLength, ux = 'n', save = 'n'):
    import statistics
    import math
    import matplotlib 
    get_ipython().run_line_magic('matplotlib', 'inline')
    import matplotlib.pyplot as plt
    from pylab import rcParams
    rcParams['figure.figsize'] = 20,10
    matplotlib.rc('xtick', labelsize=20) 
    matplotlib.rc('ytick', labelsize=20)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    plt.subplots_adjust(wspace=0.4, hspace=0.6)

    if ux == 'y':
        if scsLength == scsLength3:
            scsData = scsAvailabilityUX3
        elif scsLength == scsLength4:
            scsData = scsAvailabilityUX4
        else:
            scsData = scsAvailabilityUX5
    else:
        if scsLength == scsLength3:
            scsData = scsAvailability3
        elif scsLength == scsLength4:
            scsData = scsAvailability4
        else:
            scsData = scsAvailability5




    print('Availability of SCSs (length', '{:1d}'.format(scsLength), ')')
    print('Max:', '{:,}'.format(int(max(scsData))))
    print('Min:', '{:,}'.format(int(min(scsData))))       
    print('Mean:', '{:5.10f}'.format(statistics.mean(scsData)))
    print('Median:', '{:5.10f}'.format(statistics.median(scsData)))
    print('Stdev:', '{:5.10f}'.format(statistics.stdev(scsData)))

    #plt.hist(scsData, density=False, bins=[-1 + i*(10**k) if k == -1 else i*(10**k) for k in [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8] for i in range(1,10)])

    ax1.hist(scsData, range=(-1, 100), bins=200)
    #ax1.gca().set_yscale("log") 
    ax1.set_yscale("log") 


    ax1.set_title('Distribution of SCS Availability (-1 to 20)', fontsize = 12)
    ax1.set_xlabel('SCS Availability', fontsize = 12)
    ax1.set_ylabel('Number of SCSs', fontsize = 12)
    ax1.grid(True)
    #plt.show
    
    
    ax2.hist(scsData, bins=np.logspace(np.log10(0.0000001), np.log10(max(scsData)), 200))
    #ax2.hist(scsData, bins = 10000)
    ax2.set_xscale("log") 
    ax2.set_yscale("log") 
    ax2.set_title('Distribution of SCS Availability (0 to max)', fontsize = 12)
    ax2.set_xlabel('SCS Availability', fontsize = 12)
    ax2.set_ylabel('Number of SCSs', fontsize = 12)
    ax2.grid(True)


    if save == 'y':
        import datetime
        fig.savefig(figPath + "avail_dist" + str(datetime.datetime.now().microsecond) + ".png", format="png", dpi=600)


    plt.show()


    
    
def showAvailability(option, SCSs, aminoSeq, save = 'n'):
    
    get_ipython().run_line_magic('matplotlib', 'inline')
    import matplotlib 
    from matplotlib import pyplot as plt
    from pylab import rcParams
    rcParams['figure.figsize'] = int(len(aminoSeq)/2),10
    matplotlib.rc('xtick', labelsize=20) 
    matplotlib.rc('ytick', labelsize=20)
    
    availability = []
    index = 0
    for scs in SCSs:
        availability.append([])
        if scs == 3:
            scsAvailability = scsAvailability3
        elif scs == 4:
            scsAvailability = scsAvailability4
        elif scs == 5:
            scsAvailability = scsAvailability5
        for amino in range(len(aminoSeq) - scs + 1):
            availability[index].append(scsAvailability[scs2code(aminoSeq[amino:(amino + scs)], scs)])
        index += 1



    # Zero padding to the tail to equal the length of x_axis and y_axis.
    index = 0
    for scs in SCSs:
        for i in range(scs - 1):
            availability[index].append(0)
        index += 1

    x_axis = []
    for i in range(len(aminoSeq)):
        x_axis.append(aminoSeq[i])


    if option == 'vh':   #list horizontally the availability and the amino acid sequence.
        for i in range(len(SCSs)):
            print('')
            for position in range(len(aminoSeq)):
                print(x_axis[position], end = '')
            print('')
            for position in range(len(aminoSeq) - SCSs[i]):   
                print('{:.5g}'.format(availability[i][position]), end = ', ')
            print('{:.5g}'.format(availability[i][len(aminoSeq) - SCSs[i]]))

            

    elif option == 't': #list vartically the availability and the amino acid sequence.
        for i in range(len(SCSs)):
            print('')
            for position in range(len(aminoSeq) - SCSs[i] + 1):
                print(aminoSeq[position:position + SCSs[i]], end = ': ')
                print('{:.5g}'.format(availability[i][position]))
    else:                                #draw the graph
        x = range(len(aminoSeq))
        for i in range(len(SCSs)):
            plt.plot(x, availability[i], label = 'SCS Length = ' + str(SCSs[i]))
        plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=1, fontsize=18)    
        plt.xticks(x, x_axis)
        plt.grid(True)

        if save == 'y':
            import datetime
            fig.savefig(figPath + "avail" + str(datetime.datetime.now().microsecond) + ".png", format="png", dpi=600)

        plt.show()

 

    
def showFrequency(option, SCSs, aminoSeq, save = 'n'):
    
    get_ipython().run_line_magic('matplotlib', 'inline')
    import matplotlib 
    from matplotlib import pyplot as plt
    from pylab import rcParams
    rcParams['figure.figsize'] = int(len(aminoSeq)/2),10
    matplotlib.rc('xtick', labelsize=20) 
    matplotlib.rc('ytick', labelsize=20)
    
    frequency = []
    index = 0
    for scs in SCSs:
        frequency.append([])
        if scs == 3:
            scsFrequency = scsFrequency3
        elif scs == 4:
            scsFrequency = scsFrequency4
        elif scs == 5:
            scsFrequency = scsFrequency5
        for amino in range(len(aminoSeq) - scs + 1):
            frequency[index].append(scsFrequency[scs2code(aminoSeq[amino:(amino + scs)], scs)])
        index += 1

     
    index = 0
    for scs in SCSs:
        for i in range(scs - 1):
            frequency[index].append(0)
        index += 1
            
    x_axis = []
    for i in range(len(aminoSeq)):
        x_axis.append(aminoSeq[i])
 
    if option == 'vh':   #list horizontally the availability and the amino acid sequence.
        for i in range(len(SCSs)):
            print('')
            for position in range(len(aminoSeq)):
                print(x_axis[position], end = '')
            print('')
            for position in range(len(aminoSeq) - SCSs[i]):   
                print('{:.5g}'.format(frequency[i][position]), end = ', ')
            print('{:.5g}'.format(frequency[i][len(aminoSeq) - SCSs[i]]))     

    elif option == 't': #list vartically the availability and the amino acid sequence.
        for i in range(len(SCSs)):
            print('')
            for position in range(len(aminoSeq) - SCSs[i] + 1):
                print(aminoSeq[position:position + SCSs[i]], end = ': ')
                print('{:.5g}'.format(frequency[i][position]))
    else:                                #draw the graph
        x = range(len(aminoSeq))
        for i in range(len(SCSs)):
            plt.plot(x, frequency[i], label = 'SCS Length = ' + str(SCSs[i]))
        plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=1, fontsize=18)    
        plt.xticks(x, x_axis)
        plt.grid(True)

        if save == 'y':
            import datetime
            fig.savefig(figPath + "freq" + str(datetime.datetime.now().microsecond) + ".png", format="png", dpi=600)
        plt.show()
        

def bisect_left_reverse(a, x):
    '''
    Return idx at which x can be inserted for reversely ordered list a.
    Return the most left idx in case x exists in the list.
    '''
    if a[0] <= x:
        return 0
    if x < a[-1]:
        return len(a)
    # binary search
    ok = len(a) - 1
    ng = 0
    while (abs(ok - ng) > 1):
        mid = (ok + ng) // 2
        if a[mid] <= x:
            ok = mid
        else:
            ng = mid
    return ok


def findScsByFrequency(scsLength, frequency, bound = 5):
    
    global BigArraySize3
    global BigArraySize4
    global BigArraySize5
    
    if scsLength == 3:
        frequencyList = [d[1] for d in sortedScsFrequency3]
        rank = bisect_left_reverse(frequencyList, frequency)
        maxIndex = BigArraySize3
    elif scsLength == 4:
        frequencyList = [d[1] for d in sortedScsFrequency4]
        rank = bisect_left_reverse(frequencyList, frequency)
        maxIndex = BigArraySize4
    elif scsLength == 5:
        frequencyList = [d[1] for d in sortedScsFrequency5]
        rank = bisect_left_reverse(frequencyList, frequency)
        maxIndex = BigArraySize5
    else:
        print('No data for SCS', scsLength)

    lowerRank = rank - bound
    if lowerRank < 0:
        lowerRank = 0
    upperRank = rank + bound
    if upperRank >= maxIndex:
        upperRank = maxIndex - 1
        
 
    print("Frequency around", frequency, '->')
    if scsLength == 3:
        for i in range(lowerRank, upperRank + 1):
            print(code2scs(int(sortedScsFrequency3[i][0]), scsLength), ':', int(sortedScsFrequency3[i][1]))
    elif scsLength == 4:
        for i in range(lowerRank, upperRank + 1):
            print(code2scs(int(sortedScsFrequency4[i][0]), scsLength), ':', int(sortedScsFrequency4[i][1]))
    elif scsLength == 5:
        for i in range(lowerRank, upperRank + 1):
            print(code2scs(int(sortedScsFrequency5[i][0]), scsLength), ':', int(sortedScsFrequency5[i][1]))


def findScsByAvailability(scsLength, availability, bound = 5):
    
    global BigArraySize3
    global BigArraySize4
    global BigArraySize5
    
    if scsLength == 3:
        availabilityList = [d[1] for d in sortedScsAvailability3]
        rank = bisect_left_reverse(availabilityList, availability)
        maxIndex = BigArraySize3
    elif scsLength == 4:
        availabilityList = [d[1] for d in sortedScsAvailability4]
        rank = bisect_left_reverse(availabilityList, availability)
        maxIndex = BigArraySize4
    elif scsLength == 5:
        availabilityList = [d[1] for d in sortedScsAvailability5]
        rank = bisect_left_reverse(availabilityList, availability)
        maxIndex = BigArraySize5
    else:
        print('No data for SCS', scsLength)

    lowerRank = rank - bound
    if lowerRank < 0:
        lowerRank = 0
    upperRank = rank + bound
    if upperRank >= maxIndex:
        upperRank = maxIndex - 1
        
 
    print("Availability around", availability, '->')
    if scsLength == 3:
        for i in range(lowerRank, upperRank + 1):
            print(code2scs(int(sortedScsAvailability3[i][0]), scsLength), ':', int(sortedScsAvailability3[i][1]))
    elif scsLength == 4:
        for i in range(lowerRank, upperRank + 1):
            print(code2scs(int(sortedScsAvailability4[i][0]), scsLength), ':', int(sortedScsAvailability4[i][1]))
    elif scsLength == 5:
        for i in range(lowerRank, upperRank + 1):
            print(code2scs(int(sortedScsAvailability5[i][0]), scsLength), ':', int(sortedScsAvailability5[i][1]))




def findProteinsIncludingSCS(scs):
    import re
    #if option == 'count':
    proteinIndex = 0
    for aminoSeq in proteinData:
        count = len(re.findall(scs, aminoSeq))
        if count > 0:
            print(count, ' in ', proteinMetaData[proteinIndex])
        proteinIndex += 1
        

def listTopFrequency(until, scsLength):

    global sortedScsFrequency3
    global sortedScsFrequency4
    global sortedScsFrequency5


    if scsLength == 3:
        sortedScsFrequency = sortedScsFrequency3
    elif scsLength == 4:
        sortedScsFrequency = sortedScsFrequency4
    elif scsLength == 5:
        sortedScsFrequency = sortedScsFrequency5
    for i in range(until):
        print('{:7d}'.format(i+1), ': ', code2scs(sortedScsFrequency[i][0], scsLength), '->', '{:9.3f}'.format(sortedScsFrequency[i][1]))

        
def listTopAvailability(until, scsLength, ux = 'n'):

    global sortedScsAvailability3
    global sortedScsAvailability4
    global sortedScsAvailability5
    global sortedScsAvailabilityUX3
    global sortedScsAvailabilityUX4
    global sortedScsAvailabilityUX5


    if ux == 'y':
        if scsLength == 3:
            sortedScsAvailability = sortedScsAvailabilityUX3
        elif scsLength == 4:
            sortedScsAvailability = sortedScsAvailabilityUX4
        elif scsLength == 5:
            sortedScsAvailability = sortedScsAvailabilityUX5
    else:
        if scsLength == 3:
            sortedScsAvailability = sortedScsAvailability3
        elif scsLength == 4:
            sortedScsAvailability = sortedScsAvailability4
        elif scsLength == 5:
            sortedScsAvailability = sortedScsAvailability5



    for i in range(until):
        print('{:7d}'.format(i+1), ': ', code2scs(sortedScsAvailability[i][0], scsLength), '->', '{:9.3f}'.format(sortedScsAvailability[i][1]))

        
def listBottomFrequency(until, scsLength):
    
    global sortedScsFrequency3
    global sortedScsFrequency4
    global sortedScsFrequency5

    if scsLength == 3:
        sortedScsFrequency = sortedScsFrequency3
    elif scsLength == 4:
        sortedScsFrequency = sortedScsFrequency4
    elif scsLength == 5:
        sortedScsFrequency = sortedScsFrequency5
    for i in range(until):
        print('{:7d}'.format(i+1), ': ', code2scs(sortedScsFrequency[-1 - i][0], scsLength), '->', '{:9.3f}'.format(sortedScsFrequency[-1 - i][1]))



def listBottomAvailability(until, scsLength, ux = 'n'):
    
    global sortedScsAvailability3
    global sortedScsAvailability4
    global sortedScsAvailability5
    global sortedScsAvailabilityUX3
    global sortedScsAvailabilityUX4
    global sortedScsAvailabilityUX5

    if ux == 'y':
        if scsLength == 3:
            sortedScsAvailability = sortedScsAvailabilityUX3
        elif scsLength == 4:
            sortedScsAvailability = sortedScsAvailabilityUX4
        elif scsLength == 5:
            sortedScsAvailability = sortedScsAvailabilityUX5
    else:
        if scsLength == 3:
            sortedScsAvailability = sortedScsAvailability3
        elif scsLength == 4:
            sortedScsAvailability = sortedScsAvailability4
        elif scsLength == 5:
            sortedScsAvailability = sortedScsAvailability5

    for i in range(until):
        print('{:7d}'.format(i+1), ': ', code2scs(sortedScsAvailability[-1 - i][0], scsLength), '->', '{:9.3f}'.format(sortedScsAvailability[-1 - i][1]))

        
def listSortedFrequency(start, end, scsLength):
    
    global sortedScsFrequency3
    global sortedScsFrequency4
    global sortedScsFrequency5

    if scsLength == 3:
        sortedScsFrequency = sortedScsFrequency3
    elif scsLength == 4:
        sortedScsFrequency = sortedScsFrequency4
    elif scsLength == 5:
        sortedScsFrequency = sortedScsFrequency5
    for i in range(start-1, end):
        print('{:7d}'.format(i+1), ': ', code2scs(sortedScsFrequency[i][0], scsLength), '->', '{:9.3f}'.format(sortedScsFrequency[i][1]))



def listSortedAvailability(start, end, scsLength, ux = 'n'):
    
    global sortedScsAvailability3
    global sortedScsAvailability4
    global sortedScsAvailability5

    global sortedScsAvailabilityUX3
    global sortedScsAvailabilityUX4
    global sortedScsAvailabilityUX5

    if ux == 'y':
        if scsLength == 3:
            sortedScsAvailability = sortedScsAvailabilityUX3
        elif scsLength == 4:
            sortedScsAvailability = sortedScsAvailabilityUX4
        elif scsLength == 5:
            sortedScsAvailability = sortedScsAvailabilityUX5
    else:
        if scsLength == 3:
            sortedScsAvailability = sortedScsAvailability3
        elif scsLength == 4:
            sortedScsAvailability = sortedScsAvailability4
        elif scsLength == 5:
            sortedScsAvailability = sortedScsAvailability5

    for i in range(start-1, end):
        print('{:7d}'.format(i+1), ': ', code2scs(sortedScsAvailability[i][0], scsLength), '->', '{:9.3f}'.format(sortedScsAvailability[i][1]))



def showScsRank(scs, scsLength, ux = 'n'):
    from scipy.stats import rankdata
    UXAvailability = -2

    scsA3 = copy.deepcopy(scsAvailability3)
    if ux == 'n':
        for i in range(BigArraySize3):
            scs = code2scs(i, scsLength3)
            if 'U' in scs or 'X' in scs:
                scsA3[i] = UXAvailability

    scsA4 = copy.deepcopy(scsAvailability4)
    if ux == 'n':
        for i in range(BigArraySize4):
            scs = code2scs(i, scsLength4)
            if 'U' in scs or 'X' in scs:
                scsA4[i] = UXAvailability

    scsA5 = copy.deepcopy(scsAvailability5)
    if ux == 'n':
        for i in range(BigArraySize5):
            scs = code2scs(i, scsLength5)
            if 'U' in scs or 'X' in scs:
                scsA5[i] = UXAvailability



    if scsLength == 3:
        print('The rank of', scs, 'in Frequency: ', rankdata(-scsFrequency3, method = 'min')[scs2code(scs, scsLength)])
        print('The rank of', scs, 'in Availability: ', rankdata(-scsA3, method = 'min')[scs2code(scs, scsLength)])
    elif scsLength == 4:
        print('The rank of', scs, 'in Frequency: ', rankdata(-scsFrequency4, method = 'min')[scs2code(scs, scsLength)])
        print('The rank of', scs, 'in Availability: ', rankdata(-scsA4, method = 'min')[scs2code(scs, scsLength)])
    else:
        print('The rank of', scs, 'in Frequency: ', rankdata(-scsFrequency5, method = 'min')[scs2code(scs, scsLength)])
        print('The rank of', scs, 'in Availability: ', rankdata(-scsA5, method = 'min')[scs2code(scs, scsLength)])

              
        
def ZipfsLawAnalysis(save = 'n'):
    
    import copy    
    import matplotlib 
    get_ipython().run_line_magic('matplotlib', 'inline')
    import matplotlib.pyplot as plt
    from pylab import rcParams
    rcParams['figure.figsize'] = 10,10
    matplotlib.rc('xtick', labelsize=20) 
    matplotlib.rc('ytick', labelsize=20)
    
    fig = plt.figure()
    ax3 = fig.add_subplot(3, 1, 1)
    ax4 = fig.add_subplot(3, 1, 2)
    ax5 = fig.add_subplot(3, 1, 3)
    plt.subplots_adjust(wspace=0.4, hspace=0.6)


    scsF3 = sorted(copy.deepcopy(scsFrequency3), reverse = True)
    #scsA3.sort(reverse = True)
    scsF4 = sorted(copy.deepcopy(scsFrequency4), reverse = True)
    scsF5 = sorted(copy.deepcopy(scsFrequency5), reverse = True)
    #scsA5.sort(reverse = True)

    rank3 = [i for i in range(BigArraySize3)] 
    rank4 = [i for i in range(BigArraySize4)] 
    rank5 = [i for i in range(BigArraySize5)] 

    ax3.set_xscale("log") 
    ax3.set_yscale("log") 
    ax3.set_title('SCS Length 3', fontsize = 12)
    ax3.set_xlabel('log R', fontsize = 12)
    ax3.set_ylabel('log F', fontsize = 12)
    ax3.grid(True)
    ax3.plot(rank3, scsF3)

    ax4.set_xscale("log") 
    ax4.set_yscale("log") 
    ax4.set_title('SCS Length 4', fontsize = 12)
    ax4.set_xlabel('log R', fontsize = 12)
    ax4.set_ylabel('log F', fontsize = 12)
    ax4.grid(True)
    ax4.plot(rank4, scsF4)


    ax5.set_xscale("log") 
    ax5.set_yscale("log") 
    ax5.set_title('SCS Length 5', fontsize = 12)
    ax5.set_xlabel('log R', fontsize = 12)
    ax5.set_ylabel('log F', fontsize = 12)
    ax5.grid(True)
    ax5.plot(rank5, scsF5)

    #plt.title('log R vs log Frequency', fontsize = 12)

    if save == 'y':
        import datetime
        fig.savefig(figPath + "ZipfsLaw" + str(datetime.datetime.now().microsecond) + ".png", format="png", dpi=600)

    plt.show




def menu():
    print('initializeFromProteinDataset()')
    print('   requires the protein data file at ./ncbi_dataset/protein.aa')
    print('initializeFromScsDataset()')
    print('showBasicInformation()')
    print('showProteinLengthDistribution()')
    print('showScsFrequencyDistribution(scsLength), scsLength = 3, 4, or 5')
    print('   Example: scs.showScsFrequencyDistribution(3)')
    print('showScsAvailabilityDistribution(scsLength), scsLength = 3, 4, or 5')
    print('   Example: scs.showScsAvailabilityDistribution(5)')
    print('showScsFrequencyDistribution(scsLength), scsLength = 3, 4, or 5')
    print('   Example: scs.showScsFrequencyDistribution(5)')
    print('showScsAvailabilityDistribution(scsLength), scsLength = 3, 4, or 5')
    print('   Example: scs.showScsAvailabilityDistribution(4)')
    print("showAvailability(option, SCSs, aminoSeq)")
    print("   option = 't' (text) or g (graph), SCSs = [3,4,5], aminoSeq = 'AAACCC...'")
    print("   Example: scs.showAvailability('t', [5], 'ACCAACAVPDDCC...')")
    print("showFrequency(option, SCSs, aminoSeq)")
    print("   option = 't' (text) or 'g' (graph) SCSs = [3,4,5], aminoSeq = 'AAACCC...'")
    print("   Example: scs.showFrequency('g', [3, 4, 5], 'ACCAACAVPDDCC...')")
    print("findScsByFrequency(scsLength, frequency, around)")
    print("   scsLength = 3, 4, or 5, frequency = 1000 or any other, around = 10, 20, ...")
    print("   Example: scs.findScsByFrequency(5, 8000, 10)")
    print("findProteinsIncludingSCS(scs), scs = 'ACCAD'")
    print("   Example: scs.findProteinsIncludingSCS('ACCAD')")
    print("showScsRank(scs, scsLength), scs = 'ADMCC', scsLength = 3, 4, or 5")
    print("   Example: scs.showScsRank('ADMCC', 5)")
    print("listTopFrequency(until, scsLength), until = 100, scsLength = 3, 4, or 5")
    print("   Example: scs.listTopFrequency(100, 5)")
    print("listTopAvailability(until, scsLength), until = 100, scsLength = 3, 4, or 5")
    print("   Example: scs.listTopAvailability(80, 4)")
    print("listBottomFrequency(until, scsLength), until = 100, scsLength = 3, 4, or 5")
    print("   Example: scs.listBottomFrequency(200, 3)")
    print("listBottomAvailability(until, scsLength), until = 100, scsLength = 3, 4, or 5")
    print("   Example: scs.listBottomAvailability(300, 4)")
    print("listSortedFrequency(start, end, scsLength), start = 100, end = 200, scsLength = 3, 4, or 5")
    print("   Example: scs.listSortedFrequency(100, 200, 5)")
    print("listSortedAvailability(start, end, scsLength), start = 100, end = 200, scsLength = 3, 4, or 5")
    print("   Example: scs.listSortedAvailability(100, 300, 4)")
     
    print("ZipfsLawAnalysis()")

    

def makeScsFrequencyTable(scsLength, scsFrequency):
    for aminoSeq in proteinData:
        for amino in range(len(aminoSeq) - scsLength + 1):
            scsFrequency[scs2code(aminoSeq[amino:(amino + scsLength)], scsLength)] += 1
            TotalSCSs[scsLength] += 1
        
def makeAminoAcidFrequencyTable():
    global TotalAminoAcids
    for aminoSeq in proteinData:
        for amino in aminoSeq: 
            aminoAcidFrequency[aminoAcidCode[amino]] += 1
            TotalAminoAcids += 1

            
def numAcidsOnLeftSide(scs, index):
    count = 0
    for i in range(0, index - 1):
        if scs[index] == scs[i]:
            count += 1
    return count


def calcAvailability(scsAvailability, scsFrequency, scsLength):
    global TotalSCSs
    global TotalAminoAcids
        
    for scscode in range(len(scsFrequency)):
        scs = code2scs(scscode, scsLength)      
        E = TotalSCSs[scsLength]
        ratio = TotalSCSs[scsLength]/TotalAminoAcids
        for i in range(scsLength):
            E *= ( ratio * aminoAcidFrequency[aminoAcidCode[scs[i]]] - numAcidsOnLeftSide(scs, i))/(TotalSCSs[scsLength] - i)
        
        scsAvailability[scscode] = scsFrequency[scscode]/E - 1
        
        
def makeSortedDataset():

    import copy
    global sortedScsFrequency3
    global sortedScsAvailability3
    global sortedScsAvailabilityUX3
    global sortedScsFrequency4
    global sortedScsAvailability4
    global sortedScsAvailabilityUX4
    global sortedScsFrequency5
    global sortedScsAvailability5
    global sortedScsAvailabilityUX5
    global scsFrequency3
    global scsAvailability3
    global scsFrequency4
    global scsAvailability4
    global scsFrequency5
    global scsAvailability5

    UXAvailability = -2



    scscodeF3 = [i for i in range(BigArraySize3)]
    scsF3 = copy.deepcopy(scsFrequency3)
    scscodeA3 = [i for i in range(BigArraySize3)]
    scsA3 = copy.deepcopy(scsAvailability3)
    scsAUX3 = copy.deepcopy(scsAvailability3)
    for i in range(BigArraySize3):
        scs = code2scs(i, scsLength3)
        if 'U' in scs or 'X' in scs:
            scsA3[i] = UXAvailability
    sortedScsAvailability3 = sorted(np.column_stack((scscodeA3,scsA3)), key=lambda x: x[1], reverse = True)
    sortedScsAvailabilityUX3 = sorted(np.column_stack((scscodeA3,scsAUX3)), key=lambda x: x[1], reverse = True)
    sortedScsFrequency3 = sorted(np.column_stack((scscodeF3,scsF3)), key=lambda x: x[1], reverse = True)
     
    scscodeF4 = [i for i in range(BigArraySize4)]
    scsF4 = copy.deepcopy(scsFrequency4)
    scscodeA4 = [i for i in range(BigArraySize4)]
    scsA4 = copy.deepcopy(scsAvailability4)
    scsAUX4 = copy.deepcopy(scsAvailability4)
    for i in range(BigArraySize4):
        scs = code2scs(i, scsLength4)
        if 'U' in scs or 'X' in scs:
            scsA4[i] = UXAvailability
    sortedScsAvailability4 = sorted(np.column_stack((scscodeA4,scsA4)), key=lambda x: x[1], reverse = True)
    sortedScsAvailabilityUX4 = sorted(np.column_stack((scscodeA4,scsAUX4)), key=lambda x: x[1], reverse = True)   
    sortedScsFrequency4 = sorted(np.column_stack((scscodeF4,scsF4)), key=lambda x: x[1], reverse = True)

    scscodeF5 = [i for i in range(BigArraySize5)]
    scsF5 = copy.deepcopy(scsFrequency5)
    scscodeA5 = [i for i in range(BigArraySize5)]
    scsA5 = copy.deepcopy(scsAvailability5)
    scsAUX5 = copy.deepcopy(scsAvailability5)
    for i in range(BigArraySize5):
        scs = code2scs(i, scsLength5)
        if 'U' in scs or 'X' in scs:
            scsA5[i] = UXAvailability

    sortedScsAvailability5 = sorted(np.column_stack((scscodeA5,scsA5)), key=lambda x: x[1], reverse = True)
    sortedScsAvailabilityUX5 = sorted(np.column_stack((scscodeA5,scsAUX5)), key=lambda x: x[1], reverse = True)   
    sortedScsFrequency5 = sorted(np.column_stack((scscodeF5,scsF5)), key=lambda x: x[1], reverse = True)





def verifyFrequencyTables(scsLength, scsFrequency):
    global TotalSCSs
    global TotalAminoAcids
    global NumAminoAcids
 
    zero_count = 0
    for i in range(len(scsFrequency)):
        #TotalSCSs += int(scsFrequency[i])
        if scsFrequency[i] == 0:
            zero_count += 1

    print('# SCSs = ', TotalSCSs[scsLength], ' # of zero count SCSs = ', zero_count)

    TotalSCSsFromSeq = 0
    TotalAminoAcidsFromSeq = 0
    for aminoSeq in proteinData:
        TotalSCSsFromSeq += len(aminoSeq) - scsLength + 1
        TotalAminoAcidsFromSeq += len(aminoSeq)
    print('# SCSs from Sequences = ', TotalSCSsFromSeq)

    for i in range(NumAminoAcids):
        #TotalAminoAcids += int(aminoAcidFrequency[i])
        print('#', aminoAcidCodeInverse[i], ' = ', aminoAcidFrequency[i])
        
    print('# Amino Acids = ', TotalAminoAcids)
    print('# Amino Acids from Seq = ', TotalAminoAcidsFromSeq)
        


def getAvailability(scsLength, scsAvailability, aminoSeq, save = 'n'):
    
    get_ipython().run_line_magic('matplotlib', 'inline')
    import matplotlib 
    from matplotlib import pyplot as plt
    from pylab import rcParams
    rcParams['figure.figsize'] = int(len(aminoSeq)/2),10
    matplotlib.rc('xtick', labelsize=20) 
    matplotlib.rc('ytick', labelsize=20)
    
    availability = []
    for amino in range(len(aminoSeq) - scsLength + 1):
        availability.append(scsAvailability[scs2code(aminoSeq[amino:(amino + scsLength)], scsLength)])

    for i in range(scsLength - 1):
        availability.append(0)
        
    x_axis = []
    for i in range(len(aminoSeq)):
        x_axis.append(aminoSeq[i])
            
    x = range(len(aminoSeq))
    plt.plot(x, availability)
    plt.xticks(x, x_axis)

    if save == 'y':
        import datetime
        fig.savefig(figPath + "avail" + str(datetime.datetime.now().microsecond) + ".png", format="png", dpi=600)

    plt.show()



def initializeFromProteinDataset():

    global scsLength3
    global scsLength4
    global scsLength5
    global TotalAminoAcids
    
    readProteinDatafile()

    makeScsFrequencyTable(scsLength5, scsFrequency5)
    makeScsFrequencyTable(scsLength4, scsFrequency4)
    makeScsFrequencyTable(scsLength3, scsFrequency3)
    makeAminoAcidFrequencyTable()
    calcAvailability(scsAvailability3, scsFrequency3, scsLength3)
    calcAvailability(scsAvailability4, scsFrequency4, scsLength4)
    calcAvailability(scsAvailability5, scsFrequency5, scsLength5)
    makeSortedDataset()
    TotalAminoAcids = sum(aminoAcidFrequency)

def initializeFromScsDataset():
    
    global scsLength3
    global scsLength4
    global scsLength5

    readProteinDatafile()
    objectLoad()
    TotalAminoAcids = sum(aminoAcidFrequency)
     

    



