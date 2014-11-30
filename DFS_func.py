'''File storing the functions useful for DFS_TOOL.py'''

import re
import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as spp
# import iminuit


def parser(FILE):
    '''Function doing the bulk of the parsing, with regular expressions'''
    MAX_INT = []
    MAX_MASSES = []
    INT = []
    MASSES = []
    SCANS = 0
    for line in FILE:
        line = line.rstrip()
        # Find line with number of counts in one position
        if 'PeakIntensity' in line:
            index = 0
            SCANS += 1
            INT.append([0] * 3000)
            MASSES.append([0] * 3000)
            data = re.findall('[0-9]{0,20}\.[0-9]{0,10}', line)
            MAX_INT.append(float(data[0]))
            MAX_MASSES.append(float(data[1]))
        # Find start of a new packet
        if 'Packet #' in line:
            data = re.findall('[0-9]{0,20}\.[0-9]{0,10}', line)
            INT[SCANS - 1][index] = float(data[0])
            MASSES[SCANS - 1][index] = float(data[1])
            index += 1
    # Output of number of scans, of packets
    print 'Number of data points per scan', index
    for i in range(len(MASSES)):
        MASSES[i] = MASSES[i][0:index - 1]
        INT[i] = INT[i][0:index - 1]
    print 'Number of scans :', len(INT)

    return MAX_INT, INT, MASSES, index


def refmassfinder(MASSES):
    '''Function to find a reference list of masses'''
    # Look for extremes in initial and final masses
    # Note : observation that initial and final masses
    # are often messed up by Xcalibur.
    # MIN_M = min([MASSES[i][0] for i in range(len(MASSES))])
    # MAX_M = max([MASSES[i][-1] for i in range(len(MASSES))])
    # REF_M = []
    # Look at all mass ranges, find one which has both extremes
    # for i in range(len(MASSES)):
    #    if min(MASSES[i]) == MIN_M and max(MASSES[i]) == MAX_M:
    #        print 'Scan number', i, 'has been taken as reference '
    #        REF_M = MASSES[i]
    #        break
    # If non found, take first one which has minimum and whose
    # last value is not zero
    # if REF_M == []:
    #    for i in range(len(MASSES)):
    #        if min(MASSES[i]) == MIN_M and MASSES[-1] != 0:
    #            print 'Scan number', i, 'has been taken as reference '
    #            REF_M = MASSES[i]
    #            break
    REF_M = []
    for i in MASSES:
        for j in i:
            if j not in REF_M:
                REF_M.append(j)
    REF_M.sort()
    return REF_M


def datacorrecter(INT, MASSES, REF_M):
    '''Function to correct recorded data for repeated and omitted values,
    not perfect but solves most data reporting errors'''
    for i in range(len(MASSES)):
        if MASSES[i] != REF_M:
            for j in range(len(MASSES[i])):
                # Omitted values, at start, assign 0
                if MASSES[i][j] > REF_M[j]:
                    MASSES[i].insert(j, REF_M[j])
                    MASSES[i].pop(-1)
                    INT[i].insert(j, 0)
                    INT[i].pop(-1)
                # Repeated values, data removed
                if MASSES[i][j] < REF_M[j] and MASSES[i][j] != 0:
                    MASSES[i].pop(j)
                    MASSES[i].append(0)
                    INT[i].pop(j)
                    INT[i].append(0)
                # End of mass range errors with 0 reported
                if MASSES[i][j] == 0:
                    MASSES[i][j] = REF_M[j]
                    INT[i][j] = 0
    return INT, MASSES


def exporter(dtt, filename):
    '''Export the data to a csv file'''
    export = open(filename, 'wb')
    wrt = csv.writer(export, dialect='excel')
    for item in dtt:
        wrt.writerow(item)
    export.close()
    return


def selecter(MAX_INT, INCLUDED, STRIPPED_BEFORE, KEPT_WIDTH):
    '''Function selecting which data to keep'''
    # Assume start in stable conditions with gas flowing
    FLATTOP = 1
    # If there are drops in maximum intensity
    # Assume they are due to changing from one
    # Bellow to the other, cut according to variables
    for i in range(len(MAX_INT)):
        if (MAX_INT[i] < max(MAX_INT) / 2 and i > STRIPPED_BEFORE
                + KEPT_WIDTH and FLATTOP == 1):
            INCLUDED += [
                a for a in range(
                    i - STRIPPED_BEFORE - KEPT_WIDTH, i - STRIPPED_BEFORE)]
            FLATTOP = 0
        if MAX_INT[i] > 0.7 * max(MAX_INT) and FLATTOP == 0:
            FLATTOP = 1
    # If no drops in intensity, keep everything
    if INCLUDED == []:
        INCLUDED = range(len(MAX_INT))
    return INCLUDED

def lowSelecter(MAX_INT, INCLUDED, INT, STRIPPED_BEFORE, KEPT_WIDTH):
    '''Function with an low noise method for finding the bellows switches,
    tuned for low signal settings'''
    #Fuction works by detecting large drops in total ion intensity
    #Averaged over a window of three scans
    #First, makes a list of the signal sum for every scan
    SUM_INT=[]
    for i in range(len(MAX_INT)):
        SUM_INT.append(sum(INT[i]))
    LAST_SIGNAL=SUM_INT[0]
    FLATTOP=1
    #Now, observes the mean change in total signal, with a three-cycle window
    #Normalized to the cumulative intensity signal
    #Compared to the most recent accepted signal, to account for pressure bleed-out
    for i in range(3,len(MAX_INT)-3):
        drop=(np.mean(SUM_INT[i:i+3])-np.mean(SUM_INT[i-3:i]))/LAST_SIGNAL
        if (drop < -0.3 and i > STRIPPED_BEFORE + KEPT_WIDTH and FLATTOP == 1):
            INCLUDED += [
                a for a in range(i - STRIPPED_BEFORE - KEPT_WIDTH, i - STRIPPED_BEFORE)]
            FLATTOP = 0
            LAST_SIGNAL=SUM_INT[i-3]
        if SUM_INT[i] > 0.7 * LAST_SIGNAL and FLATTOP == 0:
            FLATTOP = 1
    #If signal is stable at the end of the scan for long enough, include the last cycle
    #if FLATTOP == 1 and INCLUDED[-1] + KEPT_WIDTH + STRIPPED_BEFORE < len(SUM_INT):
    #    INCLUDED += [
    #        a for a in range(i - STRIPPED_BEFORE - KEPT_WIDTH, i - STRIPPED_BEFORE)]

    # If no drops in intensity, keep everything
    if INCLUDED ==[]:
        INCLUDED = range(len(MAX_INT))
    return INCLUDED



def excluder(MAX_INT, INT, INCLUDED, REF_M, MASSES):
    '''Function using the included variable to trim the data'''
    CLN_INT = list(MAX_INT)
    CYC_NUM = 0
    # Create an empty result variable
    RSLT = np.zeros(shape=(100, len(REF_M)))
    counts = np.zeros(shape=(100, len(REF_M)))
    # Fill in, with one row for each cycle
    for i in range(len(MAX_INT)):
        if i not in INCLUDED:
            CLN_INT[i] = 0
        if i in INCLUDED and i - 1 not in INCLUDED:
            CYC_NUM += 1
        if i in INCLUDED:
            for j in range(len(MASSES[i])):
                if MASSES[i][j] != 0:
                    index = REF_M.index(MASSES[i][j])
                    counts[CYC_NUM - 1][index] += 1
                    RSLT[CYC_NUM - 1][index] += np.array(INT[i][j])
    for i in range(len(RSLT)):
        for j in range(len(RSLT[i])):
            if counts[i][j] != 0:
                RSLT[i][j] = RSLT[i][j]/counts[i][j] * np.mean(counts[i])
    return RSLT, CLN_INT


def bg_correction(RSLT):
    print('Background correction :')
    # Create an object to store clicks on canvas, see definition
    xvales = getxvals(RSLT.T)
    # Pause
    a = raw_input(
        'Press a key after you are done with your selection of points')
    plt.close()
    for i in range(len(RSLT)):
        if max(RSLT[i]) == 0:
            break
    bgs = []
    for peak_shape in RSLT[0:i]:
        bgs += [np.mean(peak_shape[int(xvales.x[0]): int(xvales.x[1])])]
    print('Background values: ', bgs)
    print('Average: ', np.mean(bgs))
    print('Observed standard variation: ', np.std(bgs))
    while True:
        a = raw_input('Use (A) for all, (I)ndvidual values'
                      ' or (M)anual input for all?')
        if a == 'A':
            for peak_shape in RSLT[0:i]:
                peak_shape -= np.mean(bgs)
            break
        if a == 'I':
            j = 0
            for peak_shape in RSLT[0:i]:
                peak_shape -= bgs[j]
                j += 1
            break
        if a == 'M':
            nbg = float(raw_input('Input manual background value : '))
            for peak_shape in RSLT[0:i]:
                peak_shape -= nbg
            break
    return RSLT


class getxvals:

    '''Used for manual selection of windows'''

    def __init__(self, MAX_INT):
        '''Creator function'''
        figWH = (8, 5)  # inches
        self.fig = plt.figure(figsize=figWH)
        # Plot max intensities as a canvas
        plt.plot(MAX_INT, '--')
        self.ax = self.fig.get_axes()[0]
        self.x = []  # will contain "x" values
        self.y = []  # will contain "y" values for plotting purposes
        a = raw_input('Press the (Z) key if you need to zoom ')
        print(a)
        if a == 'Z':
            print('Use the lens to zoom on the feature which'
                  ' will be used to measure the background.')
            while True:
                a = raw_input('Press the (Z) key when done')
                if a == 'Z':
                    break
        print('You can now click to select an adequate window.')
        print('Left click to select a point, right click to close the window.')
        # Calls to class functions
        self.connect = self.ax.figure.canvas.mpl_connect
        self.disconnect = self.ax.figure.canvas.mpl_disconnect
        self.clickCid = self.connect("button_press_event", self.onClick)

    def onClick(self, event):
        '''Click event function'''
        if event.inaxes:
            # Left click : add coordinates
            if event.button == 1:
                self.x.append(event.xdata)
                self.y.append(event.ydata)
                # Plot click position
                plt.plot(event.xdata, event.ydata, 'gs')
                # If even number of click, plot segment
                if len(self.x) % 2 == 0:
                    plt.plot([int(self.x[-2]), int(self.x[-1])], [
                             int(self.y[-2]), int(self.y[-1])],
                             linewidth=4, color='g')
            # Right click : disconnect
            if event.button == 3:
                self.cleanup()

    def cleanup(self):
        '''Cleaner/disconnecter'''
        self.disconnect(self.clickCid)
        plt.close()


def masker(REF_M, RSLT):
    '''To exclude part of the mass spectrum '''
    mask = [0, 0]
    mask[0] = float(raw_input('Low-mass end of the window to exclude : '))
    mask[1] = float(raw_input('High-mass end of the window to exclude : '))
    tobemasked = []
    for i in range(len(REF_M)):
        if REF_M[i] > mask[0] and REF_M[i] < mask[1]:
            tobemasked += [i]
    print tobemasked
    REF_M = np.delete(REF_M, tobemasked)
    a = []
    for i in range(len(RSLT)):
        a.append(np.delete(RSLT[i], tobemasked))
    RSLT = np.array(a)
    return REF_M, RSLT


def nearest_point(m, masses):
    '''Simple tool for masking'''
    for i in range(len(masses) - 1):
        if masses[i] <= m and masses[i + 1] > m:
            break
    return i


def peakshape(mass, center, sigma, cupwidth):
    '''Model of the recorded peak shapes, assume Gaussian ion beam'''
    return ((spp.erf((mass - center + cupwidth) / sigma) -
             spp.erf((mass - center - cupwidth) / sigma)))


def permildiff(solutions, i, j):
    calc = permilcalc(solutions, 'amp'+str(i), 'amp'+str(j))
    print('Differences in the ratio of peaks ' + str(i) + ' and ' + str(j))
    print('List in permil ')
    for rsl in calc:
        print rsl
    print('Average : ' + str(np.mean([c for c in calc])))
    print('Standard deviation : ' + str(np.std([c for c in calc])))
    print('Standard error : ' + str(
        np.std([c for c in calc])/np.sqrt(len(calc)-1)))
    return calc


def permilcalc(solutions, var1, var2):
    permils = np.zeros((len(solutions)-1)/2)
    j = 0
    for i in range(0, len(solutions)-2, 2):
        ref1 = solutions[i][var1]/solutions[i][var2]
        smp = solutions[i+1][var1]/solutions[i+1][var2]
        ref2 = solutions[i+2][var1]/solutions[i+2][var2]
        permils[j] = 1000*(2*smp/(ref1+ref2)-1)
        j += 1
    return permils


def iscan(m, var, mini, maxi, nb):
    list = np.arange(mini, maxi, (maxi - mini) / nb)
    var_mini = 0
    print list
    return var_mini
