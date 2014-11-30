#!/usr/bin/python2.7

'''Tool to extract and fit the raw data
from the text files created by Xcalibur software'''

# import re
# import csv
import sys
import matplotlib.pyplot as plt
import numpy as np
# import scipy.special as spp
import iminuit
import probfit
# import minuit
import DFS_func
from itertools import combinations

# Dynamic plottig in the whole program
plt.ion()

# For drag-and-drop capacity on windows systems
# Ask for file name otherwise
if (len(sys.argv) < 2):
    FILENAME = raw_input('Name of the file to parse ? ')
else:
    FILENAME = sys.argv[1]

FILE = open(FILENAME, 'r')
FILENAME = FILENAME.rstrip('.txt')

# EXTRACTING DATA FROM THE FILE
MAX_INT, INT, MASSES, index = DFS_func.parser(FILE)

# CORRECTION OF DATA REPORTING ERRORS
# STEP 1 : FINDING A REFERENCE LIST OF MASSES
REF_M = DFS_func.refmassfinder(MASSES)
if REF_M[0] == 0:
    REF_M = REF_M[1:]

# STEP 2 : CORRECTING MASSES AND INTENSITIES WITH REFERENCE LIST
# INCORRECT AS INSTRUMENT DOES NOT ALWAYS START/STOP AT SAME PLACES
# INT, MASSES = DFS_func.datacorrecter(INT, MASSES, REF_M)

# STEP 3 : OPTIONAL EXPORT
while True:
    choice = raw_input('Export the data ? Y/N ')
    if choice == 'Y':
        print('Masses and intensities exported as .csv files')
        DFS_func.exporter(MASSES, FILENAME + '_m.csv')
        DFS_func.exporter(INT, FILENAME + '_i.csv')
    else:
        print('No export')
        break

# DATA SELECTION
# Beginning of the data processing part
# Aims :
# 1) remove all records where the intensity is not stable yet
# 2) sort reference and sample scans,
# to bin them or to represent them on their own
# 3) export the relevant information
INCLUDED = []
STRIPPED_BEFORE = 5
KEPT_WIDTH = 5
# MASKING DATA NEAR CHANGES OF BELLOW
INCLUDED = DFS_func.selecter(MAX_INT, INCLUDED, STRIPPED_BEFORE, KEPT_WIDTH)

# CONSTRUCTING THE MASKED DATA
# Adjusting the width of the selection window
# Or choosing arbitrary time windows
while True:
    # Showing last selection done, or default
    RSLT, CLN_INT = DFS_func.excluder(MAX_INT, INT, INCLUDED, REF_M, MASSES)
    f1, axarr1 = plt.subplots(2)
    axarr1[0].plot(MAX_INT, '--')
    axarr1[0].plot(CLN_INT)
    axarr1[0].set_title('Kept values in green')
    axarr1[1].plot(REF_M, RSLT.T)
    axarr1[1].set_title('Averaged mass scans for each cycle')
    print 'Result from automatic window selection in green.'
    # Offering choice to change type of selection
    choice = raw_input(
        '(O)kay, change (W)idth or (M)anual selection of the limits ? ')
    INCLUDED = []
    plt.close()
    # Simple width change
    if choice == 'W':
        KEPT_WIDTH = int(
            raw_input('Enter new width (current is %d)' % KEPT_WIDTH))
        INCLUDED = DFS_func.selecter(
            MAX_INT, INCLUDED, STRIPPED_BEFORE, KEPT_WIDTH)
        print 'Result with the new width'
    # Manual selection of bounds
    if choice == 'M':
        plt.close()
        # Create an object to store clicks on canvas, see definition
        xvales = DFS_func.getxvals(MAX_INT)
        # Pause
        a = raw_input(
            'Press a key after you are done with your selection of points')
        plt.close()
        # Calculate new INCLUDED from click positions
        for i in range(len(xvales.x)):
            xvales.x[i] = int(xvales.x[i])
            if i % 2 != 0:
                INCLUDED += range(int(xvales.x[i - 1]), int(xvales.x[i]))
        # If odd number of click, include everything till last point
        if len(xvales.x) % 2 != 0:
            INCLUDED += [range(xvales.x[-1], len(REF_M))]
    # Keep last selection
    if choice == 'O':
        print('Mask selected')
        f1, axarr1 = plt.subplots(2)
        axarr1[0].plot(MAX_INT, '--')
        axarr1[0].plot(CLN_INT)
        axarr1[0].set_title('Kept values in green')
        axarr1[1].plot(REF_M, RSLT.T)
        axarr1[1].set_title('Averaged mass scans for each cycle')
        break
    if choice not in ['O', 'W', 'M']:
        print('Wrong entry, try again')
        INCLUDED = DFS_func.selecter(
            MAX_INT, INCLUDED, STRIPPED_BEFORE, KEPT_WIDTH)

# REMOVNG THE BACKGROUND
plt.close()
answer = raw_input('Do you want to correct for background? Y/N ')
if answer == 'Y':
    RSLT = DFS_func.bg_correction(RSLT)
    plt.close()

f1, axarr1 = plt.subplots(2)
axarr1[0].plot(MAX_INT, '--')
axarr1[0].plot(CLN_INT)
axarr1[0].set_title('Kept values in green')
axarr1[1].plot(REF_M, RSLT.T)
axarr1[1].set_title('Averaged mass scans for each cycle')


# REMOVING PART OF THE MASS SPECTRUM
print 'Do you want to keep the entire mass spectrum?'
while True:
    # Allowing to mask part of the data to fit only relevant peaks
    choice = raw_input('(O)kay, (E)xlcude part of the data ? ')
    if choice == 'E':
        REF_M, RSLT = DFS_func.masker(REF_M, RSLT)
        plt.close()
        f, axarr = plt.subplots(2)
        axarr[0].plot(MAX_INT)
        axarr[0].plot(CLN_INT, '--')
        axarr[1].plot(REF_M, RSLT.T)
    else:
        plt.close()
        f, axarr = plt.subplots(2)
        axarr[0].plot(MAX_INT)
        axarr[0].plot(CLN_INT, '--')
        axarr[1].plot(REF_M, RSLT.T)
        break

# ERRORS ESTIMATION AND NORMALIZATION
# Rationale for errors :
# Add average value to avoid giving too much weight
# to low intensity parts of the mass spectrum
# Assume reported count on 0.1 seconds scan and
# that reported is recalculated to 1 sec, not actual counts
# Removing negative values
# RSLT = abs(RSLT)
errs = 10 * index * np.sqrt(abs(RSLT) / (10 * index)) + 1
# Normalization
for i in range(len(RSLT)):
    if max(RSLT[i]) > 0:
        errs[i] = (errs[i] + np.mean(errs[i]))/2 / max(RSLT[i])
        RSLT[i] = RSLT[i] / max(RSLT[i])
    if max(RSLT[i]) == 0:
        break
# CUT EMPTY ROWS
RSLT = RSLT[0:i]
errs = errs[0:i]
REF_M = np.array(REF_M)

# EXPORT OF CLEANED DATA AND ERRORS
while True:
    choice = raw_input('Export the cleaned-up data ? Y/N ')
    if choice == 'Y':
        print('Masses and intensities exported as .csv files')
        DFS_func.exporter(MASSES, FILENAME + '_RSLT.csv')
        DFS_func.exporter(INT, FILENAME + '_ERRS.csv')
    else:
        print('No export')
        break

# MODELLING
# NB OF PEAKS
PEAKS = int(raw_input('Number of peaks to model ? '))
if type(PEAKS) == int and PEAKS > 1 and PEAKS < 5:
    print 'Trying to fit with', PEAKS, 'peaks'
else:
    PEAKS = 2
    print('Invalid number, defaulting to 2 peaks')

# Building initial guess
mass_guess = np.zeros(PEAKS)
# Asking for rough estimate of center of peaks
for i in range(PEAKS):
    string = 'Mass for peak number ' + str(i + 1) + ' ? '
    mass_guess[i] = raw_input(string)
# Looking at data to get estimates of peak heigth
# Not great if overlapping by more than 50% of width
abund_guess = 0.5 * \
    np.array([RSLT[0][DFS_func.nearest_point(i, REF_M)] for i in mass_guess])
# Arbitrary values for resolution and cup width, easy to tweak with
# a rough scan after the fit objects are created
res_guess = np.array([float(30000)])
sigma_guess = ((min(REF_M) + max(REF_M)) / 2) / \
    res_guess / 2.0 / 1.645 / np.sqrt(2.0)
cup_guess = np.array([float(0.00015 * REF_M[0] / 18)])

# CONSTRUCTING MODEL
# Building one function with arguments for each peak, multiple concurrent fits
peak0 = peak1 = peak2 = peak3 = []
for i in range(PEAKS):
    name = 'peak' + str(i)
    vars()[name] = probfit.Extended(
        probfit.rename(DFS_func.peakshape,
                       ['mass', 'center' + str(i), 'sigma', 'cupwidth']),
        extname='amp' + str(i))

# Sum those with the special iminut functions
if PEAKS == 1:
    MODEL = peak0
if PEAKS == 2:
    MODEL = probfit.AddPdf(peak0, peak1)
if PEAKS == 3:
    MODEL = probfit.AddPdf(peak0, peak1, peak2)
if PEAKS == 4:
    MODEL = probfit.AddPdf(peak0, peak1, peak2, peak3)

# STEP ONE : fit the first row (or only row)
print('Fitting the first peak shape')
# Setting up a chi2 function
chi2 = probfit.Chi2Regression(MODEL, REF_M, RSLT[0], errs[0])
# Setting up dictionnary of initial values
pars = dict()
pars['cupwidth'] = cup_guess
pars['sigma'] = sigma_guess
for i in range(PEAKS):
    pars['amp' + str(i)] = abund_guess[i]
    pars['center' + str(i)] = mass_guess[i]
    pars['error_center' + str(i)] = 0.0001
    pars['error_amp' + str(i)] = 0.001
    pars['error_cupwidth'] = 0.001
    pars['error_sigma'] = 0.001

print(abund_guess)

m = iminuit.Minuit(chi2, **pars)
plt.close()
f, axarr = plt.subplots(3)
axarr[0].plot(MAX_INT)
axarr[0].plot(CLN_INT, '--')
axarr[1].plot(REF_M, RSLT.T)
chi2.draw(m, parts=False)
print('Initial guess in red, data points and error bars in blue')
# a = raw_input('Press a key to begin playing with input parameters')
# print cup_guess
# best_try = probfit.try_chi2(
#     MODEL, RSLT[0], amp0=abund_guess[0], amp1=abund_guess[1],
#    sigma=sigma_guess, center0=mass_guess[0], center1=mass_guess[1],
#    cupwidth=[cup_guess[0]*0.1, cup_guess[0]*10])
# print best_try
# for i in np.arange(-10,10,1):
#   .....:   pars['cupwidth']= cup_guess + i * pars['error_cupwidth']
#   .....:   m = iminuit.Minuit(chi2, **pars)
#   .....:   a = chi2.draw(m ,parts = False)
#   .....:   list += [(i, np.sum((a[0][1]-a[2][1])**2))]

a = raw_input('Press a key to begin fitting the data')
m.migrad()
a = raw_input('Press a key to display the fit result')
plt.close()
chi2.draw(m, parts=False)
a = raw_input('Press a key to begin to fit all peak shapes')
plt.close()

# STEP TWO : fit all rows one by one
# Variables to store all fit parameters
variables = []
errors = []
m = []
chi2 = []
for i in range(len(RSLT)):
    chi2 += [probfit.Chi2Regression(MODEL, REF_M, RSLT[i], errs[i])]
    pars = dict()
    pars['cupwidth'] = cup_guess
    pars['sigma'] = sigma_guess
    for j in range(PEAKS):
        pars['amp' + str(j)] = abund_guess[j]
        pars['center' + str(j)] = mass_guess[j]
        pars['error_center' + str(j)] = 0.0001
        pars['error_amp' + str(j)] = 0.001
        pars['error_cupwidth'] = 0.001
        pars['error_sigma'] = 0.001
    print i
    m += [iminuit.Minuit(chi2[i], pedantic=False, print_level=0, **pars)]
    m[i].migrad()
    variables += [m[i].values]
    errors += [m[i].errors]


for j in range(min([2, 1 + len(variables) / 5])):
    for i in range(5):
        if i + j * 5 < len(variables):
            plt.subplot(2, 5, 1 + i + j * 5)
            chi2[i + j * 5].draw(m[i + j * 5])

# STEP THREE : OUTPUTS

if len(variables) > 1:
    print('Calculations of the differences between sample and reference')
    print('Done as Isodat would do')
    if PEAKS == 1:
        print('Only one peak')
        print('Average intensity : ' +
              str(np.mean([i['amp0'] for i in variables])))
        print('Precision : ' + str(np.std([i['amp0'] for i in variables])))
    if PEAKS > 1:
        A = combinations(range(PEAKS), 2)
        for i in A:
            sol = DFS_func.permildiff(variables, *i)
if len(variables) == 1:
    print('Only one continuous record of peak shape')
    for i in range(PEAKS):
        print('Peak number ' + str(i + 1))
        print('Intensity : ' +
              str(variables[0]['amp' + str(i)]) + ' +/- ' + str(
                  errors[0]['amp' + str(i)]))
