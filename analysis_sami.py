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
gn, sgn, mstar_eag, lambdaR_eag, ellip_eag, sigma_eag = np.loadtxt('LambdaGalaxiesz0_RefL100N1504.dat', unpack=True, usecols=[0,1,2,5,12,14])


massbins = [10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9]
ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <=11))
n_gals = len(mstar_SAMI[ind])

props_selected_eagle = np.zeros(shape = (4,n_gals))

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
    ids = np.arange(len(mstarin))
    selected = np.random.choice(ids, size=len_sami)
    for j in range(0,len(selected)):
        props_selected_eagle[0,g+j] = np.log10(mstarin[selected[j]])
        props_selected_eagle[1,g+j] = lambdaR_eagin[selected[j]]
        props_selected_eagle[2,g+j] = ellip_eagin[selected[j]]
        props_selected_eagle[3,g+j] = sigma_eagin[selected[j]]

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


obsdir = '/opt/claudia/ObsData/SAMI'
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
print(max(sigma_SAMI[ind]))
ind = np.where((props_selected_eagle[0,:] >= 10) & (props_selected_eagle[0,:] <= 10.9) & (props_selected_eagle[1,:] <= 0.08 + props_selected_eagle[2,:]/4.0) & (props_selected_eagle[2,:] <= 0.4))
ein = props_selected_eagle[3,ind]
ax.hist(ein[0], bins=20, range=(60,180), density=True, facecolor='red', histtype='step', linewidth=2, alpha=0.25, fill=True, edgecolor='k', label='EAGLE')

common.savefig(obsdir, fig, "lambdaR_Eps_plane.pdf")



#################################################### analysis of classification ###########################################

#read table with classification
galid, class1, class2, class3 = np.loadtxt('Classification/Classification.txt', unpack=True, usecols=[0,1,3,5])
classification = np.zeros(shape = (len(galid), 3))
classification[:,0] = class1
classification[:,1] = class2
classification[:,2] = class3

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
obsdir = '/opt/claudia/ObsData/SAMI/Classification/Plots/'
common.savefig(obsdir, fig, "MatchesInClass.pdf")

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
       plt.xticks([0,1,2,3,4,5]," ")
    else:
       plt.xticks([0,1,2,3,4,5], x_values)

    plt.yticks([])
    
obsdir = '/opt/claudia/ObsData/SAMI/Classification/Plots/'
common.savefig(obsdir, fig, "KinClassVsMatchingSuccess.pdf")


print(successmatching)





