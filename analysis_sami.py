import numpy as np
import common

mlow = 0
mupp = 1
dm = 0.05
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

catID, lambdaR_SAMI_nocorr = np.loadtxt('jvds_stelkin_cat_v012_mge_private.dat', unpack=True, usecols=[0, 41])
ids_table = np.arange(0,len(lambdaR_SAMI_nocorr))
mstar_SAMI, lambdaR_SAMI, ellip_SAMI, sigma_SAMI = np.loadtxt('jvds_stelkin_cat_v012_mge_seecorr_kh20_v190620_private.dat', unpack=True, usecols=[6, 41, 15, 35])
gn, sgn, mstar_eag, lambdaR_eag_eo, lambdaR_eag, ellip_eag_eo, ellip_eag, sigma_eag = np.loadtxt('LambdaGalaxiesz0_RefL100N1504.dat', unpack=True, usecols=[0,1,2,4,5,11,12,14])
galidl = np.loadtxt('GalaxyIDsLambdaGalaxiesz0_RefL100N1504.dat', unpack=True, usecols=[0])


massbins = [10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9]
ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <=11))
n_gals = len(mstar_SAMI[ind])

props_selected_eagle = np.zeros(shape = (5,n_gals))

g = 0
for i in range(0,len(massbins)-1):
    ind = np.where((mstar_SAMI >= massbins[i]) & (mstar_SAMI < massbins[i+1]) & (lambdaR_SAMI >= 0))
    len_sami = len(mstar_SAMI[ind]) 
    ind = np.where((mstar_eag >= 10.0**massbins[i]) & (mstar_eag < 10.0**massbins[i+1]))
    print("number of SAMI galaxies in bin", len_sami, " number of EAGLE galaxies in bin", len(mstar_eag[ind]))
    mstarin = mstar_eag[ind]
    lambdaR_eagin = lambdaR_eag[ind]
    ellip_eagin = ellip_eag[ind]
    sigma_eagin = sigma_eag[ind]
    galid_in = galidl[ind]
    ids = np.arange(len(mstarin))
    selected = np.random.choice(ids, size=len_sami)
    for j in range(0,len(selected)):
        props_selected_eagle[0,g+j] = np.log10(mstarin[selected[j]])
        props_selected_eagle[1,g+j] = lambdaR_eagin[selected[j]]
        props_selected_eagle[2,g+j] = ellip_eagin[selected[j]]
        props_selected_eagle[3,g+j] = sigma_eagin[selected[j]]
        props_selected_eagle[4,g+j] = galid_in[selected[j]]

    g = g + len(selected)

ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <= 10.5) & (lambdaR_SAMI >= 0))
H, _ = np.histogram(lambdaR_SAMI[ind],bins=np.append(mbins, mupp))
ylambdasamiL = H
H, _ = np.histogram(ellip_SAMI[ind],bins=np.append(mbins, mupp))
yellipsamiL = H

ind = np.where((props_selected_eagle[0,:] >= 10) & (props_selected_eagle[0,:] <= 10.5))
H, _ = np.histogram(props_selected_eagle[1,ind],bins=np.append(mbins, mupp))
ylambdaeagL = H
H, _ = np.histogram(props_selected_eagle[2,ind],bins=np.append(mbins, mupp))
yellipeagL = H

ind = np.where((mstar_SAMI >= 10.5) & (mstar_SAMI <= 10.9) & (lambdaR_SAMI >= 0))
H, _ = np.histogram(lambdaR_SAMI[ind],bins=np.append(mbins, mupp))
ylambdasamiH = H
H, _ = np.histogram(ellip_SAMI[ind],bins=np.append(mbins, mupp))
yellipsamiH = H

ind = np.where((props_selected_eagle[0,:] >= 10.5) & (props_selected_eagle[0,:] <= 10.9))
H, _ = np.histogram(props_selected_eagle[1,ind],bins=np.append(mbins, mupp))
ylambdaeagH = H
H, _ = np.histogram(props_selected_eagle[2,ind],bins=np.append(mbins, mupp))
yellipeagH = H


plt = common.load_matplotlib()

xtit="$\\lambda_{\\rm r_{50}}$"
ytit="$N$"

fig = plt.figure(figsize=(8.5,8.5))
plt.subplots_adjust(left=0.18, bottom=0.17)

xmin, xmax, ymin, ymax = 0, 1, 0, 100

ax = fig.add_subplot(221)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.5, 0.5, 20, 20))
ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <= 10.5) & (lambdaR_SAMI >= 0))
ax.hist(lambdaR_SAMI[ind], bins=20, range=(0,1), facecolor='blue', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k', label='SAMI')
ax.hist(lambdaR_SAMI_nocorr[ind], bins=20, range=(0,1), edgecolor='blue', histtype='step', linewidth=2, alpha=0.85, label='SAMI no corr')

ind = np.where((props_selected_eagle[0,:] >= 10) & (props_selected_eagle[0,:] <= 10.5))
lin = props_selected_eagle[1,ind]
ax.hist(lin[0], bins=20, range=(0,1), facecolor='red', histtype='step', linewidth=2, alpha=0.25, fill=True, edgecolor='k', label='EAGLE')


ax = fig.add_subplot(222)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ' ', locators=(0.5, 0.5, 20, 20))
ind = np.where((mstar_SAMI >= 10.5) & (mstar_SAMI <= 10.9) & (lambdaR_SAMI >= 0))
ax.hist(lambdaR_SAMI[ind], bins=20, range=(0,1), facecolor='blue', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k', label='SAMI')
ax.hist(lambdaR_SAMI_nocorr[ind], bins=20, range=(0,1), edgecolor='blue', histtype='step', linewidth=2, alpha=0.85, label='SAMI no corr')

ind = np.where((props_selected_eagle[0,:] >= 10.5) & (props_selected_eagle[0,:] <= 10.9))
lin = props_selected_eagle[1,ind]
ax.hist(lin[0], bins=20, range=(0,1), facecolor='red', histtype='step', linewidth=2, alpha=0.25, fill=True, edgecolor='k', label='EAGLE')
common.prepare_legend(ax, ['b','b','r'], loc='upper left')

xtit="$\\epsilon_{\\rm r_{50}}$"
ax = fig.add_subplot(223)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.5, 0.5, 20, 20))
ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <= 10.5) & (lambdaR_SAMI >= 0))
ax.hist(ellip_SAMI[ind], bins=20, range=(0,1), facecolor='blue', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')
ind = np.where((props_selected_eagle[0,:] >= 10) & (props_selected_eagle[0,:] <= 10.5))
lin = props_selected_eagle[2,ind]
ax.hist(lin[0], bins=20, range=(0,1), facecolor='red', histtype='step', linewidth=2, alpha=0.25, fill=True, edgecolor='k')

ax.text(0.32,90,'$\\rm log_{10}(M_{\\rm star}/M_{\\odot})=[10,10.5)$')

ax = fig.add_subplot(224)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ' ', locators=(0.5, 0.5, 20, 20))
ind = np.where((mstar_SAMI >= 10.5) & (mstar_SAMI <= 10.9) & (lambdaR_SAMI >= 0))
ax.hist(ellip_SAMI[ind], bins=20, range=(0,1), facecolor='blue', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')
ind = np.where((props_selected_eagle[0,:] >= 10.5) & (props_selected_eagle[0,:] <= 10.9))
lin = props_selected_eagle[2,ind]
ax.hist(lin[0], bins=20, range=(0,1), facecolor='red', histtype='step', linewidth=2, alpha=0.25, fill=True, edgecolor='k')

ax.text(0.27,90,'$\\rm log_{10}(M_{\\rm star}/M{\\odot})=[10.5,10.9]$')


obsdir = 'Plots/'
common.savefig(obsdir, fig, "lambdaR_dist.pdf")

xtit="$\\epsilon_{\\rm r_{50}}$"
ytit="PDF"

fig = plt.figure(figsize=(5,7.5))
plt.subplots_adjust(left=0.18, bottom=0.17)
xmin, xmax, ymin, ymax = 0, 0.4, 0, 7.5
ax = fig.add_subplot(211)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 0.1, 1, 1))

ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <= 10.9) & (lambdaR_SAMI >= 0) & (lambdaR_SAMI <= 0.08 + ellip_SAMI/4.0) & (ellip_SAMI <=0.4))
ax.hist(ellip_SAMI[ind], bins=20, range=(0,0.5), density=True, facecolor='blue', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k',label='SAMI')
ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <= 10.9) & (lambdaR_SAMI >= 0) & (lambdaR_SAMI_nocorr <= 0.08 + ellip_SAMI/4.0) & (ellip_SAMI <=0.4))
ax.hist(ellip_SAMI[ind], bins=20, range=(0,0.5), density=True, edgecolor='blue', histtype='step', linewidth=2, alpha=0.85, label='SAMI no corr')

ind = np.where((props_selected_eagle[0,:] >= 10) & (props_selected_eagle[0,:] <= 10.9) & (props_selected_eagle[1,:] <= 0.08 + props_selected_eagle[2,:]/4.0) & (props_selected_eagle[2,:] <= 0.4))
ein = props_selected_eagle[2,ind]
ax.hist(ein[0], bins=20, range=(0,0.5), density=True, facecolor='red', histtype='step', linewidth=2, alpha=0.25, fill=True, edgecolor='k', label='EAGLE')
common.prepare_legend(ax, ['b','b','r'], loc='upper left')

ax = fig.add_subplot(212)
xtit="$\\sigma_{\\rm r_{50}}$"
xmin, xmax, ymin, ymax = 50, 240, 0, 0.035
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(25, 25, 0.01, 0.01))

ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <= 10.9) & (lambdaR_SAMI >= 0) & (lambdaR_SAMI <= 0.08 + ellip_SAMI/4.0) & (ellip_SAMI <=0.4) & (sigma_SAMI > 0))
ax.hist(sigma_SAMI[ind], bins=20, range=(50,240), density=True, facecolor='blue', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k',label='SAMI')
ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <= 10.9) & (lambdaR_SAMI >= 0) & (lambdaR_SAMI_nocorr <= 0.08 + ellip_SAMI/4.0) & (ellip_SAMI <=0.4) & (sigma_SAMI > 0))
ax.hist(sigma_SAMI[ind], bins=20, range=(50,240), density=True, edgecolor='blue', histtype='step', linewidth=2, alpha=0.85, label='SAMI no corr')
ind = np.where((props_selected_eagle[0,:] >= 10) & (props_selected_eagle[0,:] <= 10.9) & (props_selected_eagle[1,:] <= 0.08 + props_selected_eagle[2,:]/4.0) & (props_selected_eagle[2,:] <= 0.4))
ein = props_selected_eagle[3,ind]
galidsmatched = props_selected_eagle[4,ind]
ax.hist(ein[0], bins=20, range=(60,180), density=True, facecolor='red', histtype='step', linewidth=2, alpha=0.25, fill=True, edgecolor='k', label='EAGLE')

common.savefig(obsdir, fig, "lambdaR_Eps_plane.pdf")



#################################################### analysis of classification ###########################################

#match GalaxyId with Group and SubGroupNumber
def match_galid(gnl, sgnl):
    galid, gn, sgn = np.loadtxt('/home/clagos/Trabajo/EAGLE/SelectionGalaxiesz0-L100N1504.dat', unpack=True, usecols=[0,1,2])
    galid_to_write = np.zeros(shape = len(gnl))
    g=0
    for a,b in zip(gnl, sgnl):
        match = np.where((gn == a) & (sgn == b))
        if(len(gn[match]) > 0):
           galid_to_write[g] = galid[match]
           g = g + 1
    return galid_to_write

def obtain_ellipticity_and_lambda(galid):
    props = np.zeros(shape = (2,len(galid)))
    for i in range(0,len(galid)):
        match = np.where(galidl == galid[i])
        if(len(galidl[match]) > 0):
           props[0,i] = lambdaR_eag_eo[match]
           props[1,i] = ellip_eag_eo[match]
    return props

#read table with classification
galid, class1, class2, class3 = np.loadtxt('Classification/Classification.txt', unpack=True, usecols=[0,1,3,5])
commonids = np.in1d(galid, galidsmatched)
matched = np.where(commonids == True)
#print(galid[matched])

all_classified_srs = np.in1d(galidl, galid)
#print(all_classified_srs)
#for a,b,c in zip(galid, ellip_eag_eo[all_classified_srs], ellip_eag[all_classified_srs]):
#    print(a,b,c)


classification = np.zeros(shape = (len(galid), 3))
classification[:,0] = class1
classification[:,1] = class2
classification[:,2] = class3

props = obtain_ellipticity_and_lambda(galid)
#determine the number of galaxies that have three matches, two matches, or all different answers, use flags=2, 1, 0 to determine this
successmatching = np.zeros(shape = (len(galid)))
successmatching[:] = -1

for i in range(0,len(successmatching)):
    if((class1[i] == class2[i]) & (class2[i] == class3[i])):
        successmatching[i] = 2
    elif((class1[i] == class2[i]) | (class2[i] == class3[i]) | (class1[i] == class3[i])):
        successmatching[i] = 1
    elif((class1[i] != class2[i]) & (class2[i] != class3[i]) & (class1[i] != class3[i])):
        successmatching[i] = 0

###############################################################################3
#plot success rate of matching procedure

xtit="Matching success"
ytit="N"
x_values = ['$0\%$', '$66\%$', '$100\%$']
fig = plt.figure(figsize=(3,5))
plt.subplots_adjust(left=0.25, bottom=0.17)
xmin, xmax, ymin, ymax = -0.5, 2.5, 0, 150
ax = fig.add_subplot(111)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 50, 50))
ax.hist(successmatching, bins=3, range=(-0.5,2.5), facecolor='blue', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')

plt.xticks([0,1,2], x_values)
obsdir = 'Classification/Plots/'
common.savefig(obsdir, fig, "MatchesInClass.pdf")


#################################################################
#plot distribution of kinematic classes for the three match success %

xtit="Kinematic class"
ytit="PDF"
x_values = ['FSR', 'RSR', '2$\sigma$', 'Prol','Uncl','R']
fig = plt.figure(figsize=(3,5))
plt.subplots_adjust(left=0.2, bottom=0.15)
xmin, xmax, ymin, ymax = -0.5, 2.5, 0, 1
def plot_kin_class(ax, match=0):
    ind = np.where(successmatching[:] == match)
    ax.hist(class1[ind], bins=6, range=(-0.5,5.5), density=True, facecolor='orange', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')
    ax.hist(class2[ind], bins=6, range=(-0.5,5.5), density=True, facecolor='red', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')
    ax.hist(class3[ind], bins=6, range=(-0.5,5.5), density=True, facecolor='salmon', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')

subplots = (311, 312, 313)
matches = (0, 1, 2)
text_insert = ['$0\%$', '$66\%$', '$100\%$']
xpos, ypos = xmax - 0.75 * (xmax - xmin)/2.0, ymax - 0.15 * (ymax - ymin)
for i in range(0,len(matches)):
    ax = fig.add_subplot(subplots[i])
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 50, 50))
    ax.text(xpos, ypos, 'Match succ = %s' % text_insert[i])
    plot_kin_class(ax, match=matches[i])
    if(i <= 1):
       plt.xticks([0,1,2,3,4,5],[" ", " "," "," "," "," "])
    else:
       plt.xticks([0,1,2,3,4,5], x_values)

    plt.yticks([])
    
common.savefig(obsdir, fig, "KinClassVsMatchingSuccess.pdf")

###################################################################
#plot ellipticity distribution for slow rotators in classes 0 and 1
xtit="$\\epsilon_{\\rm r_{50},edge-on}$"
ytit="N"
fig = plt.figure(figsize=(3,5))
plt.subplots_adjust(left=0.25, bottom=0.17)
xmin, xmax, ymin, ymax = 0, 0.6, 0, 30
subplots = (211, 212)
matches = (0, 1)
labels = ('FSR', 'RSR')
xpos, ypos = xmin + 0.1 * (xmax - xmin)/2.0, ymax - 0.15 * (ymax - ymin)
def plot_elip_class(ax, classification, match=0, color='red'):
    ind = np.where(classification[:] == match)
    ein = props[1,ind]
    ax.hist(ein[0], bins=10, range=(0,0.6), facecolor=color, histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')

for i in range(0,len(matches)):
    ax = fig.add_subplot(subplots[i])
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 0.1, 5, 5))
    ax.text(xpos, ypos, 'Kin class = %s' % labels[i])
    plot_elip_class(ax, class1, match=matches[i], color='orange')
    plot_elip_class(ax, class2, match=matches[i], color='red')
    plot_elip_class(ax, class3, match=matches[i], color='salmon')

obsdir = 'Classification/Plots/'
common.savefig(obsdir, fig, "EllipticityInClass.pdf")


###################################################################################
###################################################################################
###################################################################################
# Analysis of merger history of different kinematic classes

mergerdata = np.loadtxt('MergerHistoryCumulative.dat')
def find_classification(selective = True):
    if(selective):
       matchedobjs = np.where(successmatching >= 0)
       classificationin = classification[matchedobjs,:]
       classificationin = classificationin[0]
       finalclass = np.zeros(shape = len(galid[matchedobjs]))
       finalgalid = galid[matchedobjs]
       for i in range(0, len(galid[matchedobjs])):
           finalclass[i] = np.median(classificationin[i,:])
       matched = np.in1d(galidl, finalgalid)
       match = np.where(matched == True)
       finalgn = gn[match]
       finalsgn = sgn[match]
       return(finalclass, finalgn, finalsgn, finalgalid)
    else:
       matched = np.in1d(galidl, galid)
       match = np.where(matched == True)
       finalgn = gn[match]
       finalsgn = sgn[match]
       return(class3, finalgn, finalsgn, galid)


finalclass, finalgn, finalsgn, finalgalid = find_classification(selective = True)


mergerhist_srs = np.zeros(shape = (12,len(finalgn)))
for i in range(0,len(finalgn)):
    match = np.where((mergerdata[:,0] == finalgn[i]) & (mergerdata[:,1] == finalsgn[i]))
    mergerhist_srs[0:8,i] = mergerdata[match,2:10]
    mergerhist_srs[8,i]   = np.sum(mergerdata[match,3:6]) #all minor mergers
    mergerhist_srs[9,i]   = np.sum(mergerdata[match,6:9]) #all major mergers
    mergerhist_srs[10,i]  = mergerdata[match,3] + mergerdata[match,6] #dry mergers
    mergerhist_srs[11,i]  = mergerdata[match,4] + mergerdata[match,5] + mergerdata[match,7] + mergerdata[match,8] #wet mergers

############################
# plot distribution of kinematic classes depending on merger history
subplots = (611, 612, 613, 614, 615, 616)
mergers_explore = (9, 8, 0, 7, 10, 11)
labels = ('Major mergers', 'Minor mergers', 'Very minor mergers', 'No mergers', 'Dry mergers', 'Wet mergers')
xtit="Kinematic class"
ytit="N"
fig = plt.figure(figsize=(3.5,8))
plt.subplots_adjust(left=0.25, bottom=0.17)
x_values = ['FSR', 'RSR', '2$\sigma$', 'Prol','Uncl','R']
colors = ['k', 'darkgreen', 'blue', 'red', 'Navy', 'Yellow']
xmin, xmax, ymin, ymax = -0.5, 5.5, 0, 60

xpos, ypos = xmax - 0.63 * (xmax - xmin), ymax - 0.15 * (ymax - ymin)

def plot_kin_class_mergers(ax, classi, color='red'):
    ax.hist(classi, bins=6, range=(-0.5,5.5), facecolor=color, histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')

for  i in range(0,len(mergers_explore)):
    ax = fig.add_subplot(subplots[i])
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 15, 15))
    ax.text(xpos, ypos, labels[i])
    if i == 0:
       ind = np.where(mergerhist_srs[mergers_explore[i],:] > 0)
    elif i == 1:
       ind = np.where((mergerhist_srs[mergers_explore[i-1],:]== 0) & (mergerhist_srs[mergers_explore[i],:] > 0))
    elif i == 2:
       ind = np.where((mergerhist_srs[mergers_explore[i-2],:]== 0) & (mergerhist_srs[mergers_explore[i-1],:]== 0) & (mergerhist_srs[mergers_explore[i],:] > 0))
    elif i == 3:
       ind = np.where(mergerhist_srs[mergers_explore[i],:] == 0)
    elif i == 4:
       ind = np.where(mergerhist_srs[mergers_explore[i],:] > 0)
    elif i == 5:
       ind = np.where((mergerhist_srs[mergers_explore[i-1],:] == 0) & (mergerhist_srs[mergers_explore[i],:] > 0))

    plot_kin_class_mergers(ax, finalclass[ind], color=colors[i])
    if(i <= 4):
       plt.xticks([0,1,2,3,4,5],[" ", " "," "," "," "," "])
    else:
       plt.xticks([0,1,2,3,4,5], x_values)

common.savefig(obsdir, fig, "MergerHistoryDistributionKinClasses.pdf")


################
#find distribution of gas fractions and mass ratios of the last merger of prolate galaxies

def find_props_last_merger(gn, sgn, merger_thresh = 0):
    groupnumbers = np.loadtxt('GroupNumberHistory-Ref-L100.dat')
    subgroupnumbers = np.loadtxt('SubGroupNumberHistory-Ref-L100.dat')
    gasfractions = np.loadtxt('GasFractionMergerHistory-Ref-L100.dat')
    massfractions = np.loadtxt('MergerRatioHistory-Ref-L100.dat')
    snap, z, lbt = np.loadtxt('SnapshotsEAGLE.dat', unpack = True, usecols = (0,1,2))

    gasfracs_last_merger = np.zeros(shape = len(gn)) 
    massrat_last_merger = np.zeros(shape = len(gn))
    lbt_last_merger = np.zeros(shape = len(gn))
    z0gn = groupnumbers[:,0]
    z0sgn = subgroupnumbers[:,0]

    #loop for all input galaxies and find their last merger
    for i in range(len(gn)):
        match = np.where((z0gn == gn[i]) & (z0sgn == sgn[i]))
        massfracthisgal = massfractions[match[0],:]
        gasfracthisgal = gasfractions[match[0],:]
        lastmer = 0
        s = 0
        while ((lastmer <= merger_thresh) & (s < 20)):
               lastmer = massfracthisgal[0,s] 
               s = s + 1
        gasfracs_last_merger[i] = gasfracthisgal[0,s-1]
        massrat_last_merger[i] = massfracthisgal[0,s-1]
        lbt_last_merger[i] = lbt[s-1]

    #apply maximum gas fraction
    ind = np.where(gasfracs_last_merger > 1.5)
    gasfracs_last_merger[ind] = 1.5
    return(massrat_last_merger, gasfracs_last_merger, lbt_last_merger)

def info_last_merger_kinclass(kinclass, gn, sgn, selection=0, merger_thresh = 0):
    prolates = np.where(kinclass == selection)
    gnpro = gn[prolates]
    sgnpro = sgn[prolates]
    props = np.zeros(shape = (3, len(gnpro)))
    (massfrac_pro, gasfrac_pro, lbt_pro) = find_props_last_merger(gnpro, sgnpro, merger_thresh = merger_thresh)
    props[0,:] = massfrac_pro[:]
    props[1,:] = gasfrac_pro[:]
    props[2,:] = lbt_pro[:]

    return props

############################
# plot distribution of merger properties of kinematic classes
subplots = (311, 312, 313)
fig = plt.figure(figsize=(3,7))
plt.subplots_adjust(left=0.25, bottom=0.17)
kinclasses = [3,2,1,0] #,1,2,3]
colors = ['Chocolate', 'YellowGreen', 'SteelBlue', 'DarkMagenta']
labels = ['Prol', '2$\sigma$', 'RSR', 'FSR'] #, '2$\sigma$', 'Prol']

props=[0,1,3]
ymin = (0.17, 0.18, 0)
xmin = (0, 0, 1)
xmax = (1, 1.5, 12.)
ymax = (1.05, 1.05, 1.05)
ticks = (0.2, 0.2, 2)
xtit=["$M_{\\star,\\rm sec}/M_{\\star,\\rm prim}$", "$M_{\\rm SF,total}/M_{\\star,\\rm total}$", "lookback time [Gyr]"]
ytit="Cumulative dist."

for p in range(0,3):
    ax = fig.add_subplot(subplots[p])
    common.prepare_ax(ax, xmin[p], xmax[p], ymin[p], ymax[p], xtit[p], ytit, locators=(ticks[p], ticks[p], 100, 100))
    for i in range(0,len(kinclasses)):
        props = info_last_merger_kinclass(finalclass, finalgn, finalsgn, selection = kinclasses[i])
        true_mergers = np.where(props[0,:] > 0)
        xin = props[p,true_mergers]
        ax.hist(xin[0], bins=15, range=(xmin[p], xmax[p]), density = True, cumulative = True, facecolor=colors[i], histtype='step', linewidth=2, alpha = 0.85, fill=False, edgecolor=colors[i], label=labels[i])
        ax.axes.yaxis.set_ticks([])
    if p == 0:
       xin = [0.3, 0.3]
       yin = [0, 1.05]
       ax.plot(xin, yin, linestyle='dotted', color='k')
       ax.arrow(0.3, 0.95, 0.13, 0)
       ax.text(0.3,0.965, 'major', color='k')
       xin = [0.1, 0.1]
       ax.plot(xin, yin, linestyle='dotted', color='k')
       ax.arrow(0.1, 0.95, 0.07, 0)
       ax.arrow(0.3, 0.95, -0.07, 0)
       ax.text(0.11,0.965, 'minor', color='k')
    if p == 1:
       xin = [0.1, 0.1]
       yin = [0, 1.05]
       ax.plot(xin, yin, linestyle='dotted', color='k')
       ax.arrow(0.1, 0.95, 0.13, 0)
       ax.text(0.1,0.965, 'wet', color='k')

       common.prepare_legend(ax, colors, loc='lower right')
plt.tight_layout()
common.savefig(obsdir, fig, "CumDistributionMergerParametersKinClasses.pdf")

fig = plt.figure(figsize=(3,7))
plt.subplots_adjust(left=0.25, bottom=0.1)
ymin =  (0.1, 0.1, 0.1) #0.29, 0.3, 0.04)
ymax = (6, 5, 0.5)
ytit="PDF"

for p in range(0,3):
    ax = fig.add_subplot(subplots[p])
    common.prepare_ax(ax, xmin[p], xmax[p], ymin[p], ymax[p], xtit[p], ytit, locators=(ticks[p], ticks[p], 15, 15))
    for i in range(0,len(kinclasses)):
        props = info_last_merger_kinclass(finalclass, finalgn, finalsgn, selection = kinclasses[i])
        true_mergers = np.where(props[0,:] > 0)
        xin = props[p,true_mergers]
        ax.hist(xin[0], bins=15, range=(xmin[p], xmax[p]), density = True, facecolor=colors[i], histtype='step', linewidth=2, alpha = 0.85, fill=False, edgecolor=colors[i], label=labels[i])
    if p == 0:
       xin = [0.3, 0.3]
       yin = [0, 6]
       ax.plot(xin, yin, linestyle='dotted', color='k')
       xin = [0.1, 0.1]
       ax.plot(xin, yin, linestyle='dotted', color='k')
       ax.arrow(0.3, 5.1, 0.13, 0, head_width=0.12, head_length=0.025)
       ax.text(0.3,5.3, 'major', color='k')
       ax.arrow(0.1, 5.1, 0.07, 0, head_width=0.12, head_length=0.025)
       ax.arrow(0.3, 5.1, -0.07, 0, head_width=0.12, head_length=0.025)
       ax.text(0.11,5.3, 'minor', color='k')
    if p == 1:
       xin = [0.1, 0.1]
       yin = [0, 5]
       ax.plot(xin, yin, linestyle='dotted', color='k')
       ax.arrow(0.1, 4.2, 0.18, 0, head_width=0.12, head_length=0.025)
       ax.text(0.1,4.4, 'wet', color='k')
       common.prepare_legend(ax, colors, loc='upper right')
plt.tight_layout()
common.savefig(obsdir, fig, "DistributionMergerParametersKinClasses.pdf")

############### properties of last major merger #########################

thresholds = [0.099, 0.299]
titles = ['Minor', 'Major']
kinclasses = [3,1,0] #,1,2,3]
colors = ['Chocolate', 'SteelBlue', 'DarkMagenta']
labels = ['Prol', 'RSR', 'FSR'] #, '2$\sigma$', 'Prol']


for t in range(0,len(thresholds)):
    fig = plt.figure(figsize=(3,7))
    plt.subplots_adjust(left=0.25, bottom=0.1)
    ymin =  (0.1, 0.1, 0.1) #0.29, 0.3, 0.04)
    ymax = (6, 4, 0.4)
    ytit="PDF"
   
    for p in range(0,3):
        ax = fig.add_subplot(subplots[p])
        common.prepare_ax(ax, xmin[p], xmax[p], ymin[p], ymax[p], xtit[p], ytit, locators=(ticks[p], ticks[p], 15, 15))
        for i in range(0,len(kinclasses)):
            if(t == 0):
               ind = np.where(mergerhist_srs[9,:] == 0)
               props = info_last_merger_kinclass(finalclass[ind], finalgn[ind], finalsgn[ind], selection = kinclasses[i], merger_thresh = thresholds[t])
            else:
               props = info_last_merger_kinclass(finalclass, finalgn, finalsgn, selection = kinclasses[i], merger_thresh = thresholds[t])
            true_mergers = np.where(props[0,:] > 0)
            xin = props[p,true_mergers]
            ax.hist(xin[0], bins=15, range=(xmin[p], xmax[p]), density = True, facecolor=colors[i], histtype='step', linewidth=2, alpha = 0.85, fill=False, edgecolor=colors[i], label=labels[i])
        if p == 0:
           xin = [0.3, 0.3]
           yin = [0, 6]
           ax.plot(xin, yin, linestyle='dotted', color='k')
           xin = [0.1, 0.1]
           ax.plot(xin, yin, linestyle='dotted', color='k')
           ax.arrow(0.3, 5.1, 0.13, 0, head_width=0.12, head_length=0.025)
           ax.text(0.3,5.3, 'major', color='k')
           ax.arrow(0.1, 5.1, 0.07, 0, head_width=0.12, head_length=0.025)
           ax.arrow(0.3, 5.1, -0.07, 0, head_width=0.12, head_length=0.025)
           ax.text(0.11,5.3, 'minor', color='k')
        if p == 1:
           xin = [0.1, 0.1]
           yin = [0, 5]
           ax.plot(xin, yin, linestyle='dotted', color='k')
           ax.arrow(0.1, 3.2, 0.18, 0, head_width=0.12, head_length=0.025)
           ax.text(0.1,3.4, 'wet', color='k')
           common.prepare_legend(ax, colors, loc='upper right')
           ax.text(0.55,ymax[p]+0.05, titles[t], color='k')

    plt.tight_layout()
    common.savefig(obsdir, fig, "Distribution" + titles[t] + "MergerParametersKinClasses.pdf")

