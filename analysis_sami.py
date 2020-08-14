import numpy as np
import common

AdoptProcedureSimSpin = False

mlow = 0
mupp = 1
dm = 0.05
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

catID, lambdaR_SAMI_nocorr = np.loadtxt('jvds_stelkin_cat_v012_mge_private.dat', unpack=True, usecols=[0, 41])
ids_table = np.arange(0,len(lambdaR_SAMI_nocorr))
mstar_SAMI, lambdaR_SAMI, ellip_SAMI, sigma_SAMI = np.loadtxt('jvds_stelkin_cat_v012_mge_seecorr_kh20_v190620_private.dat', unpack=True, usecols=[6, 41, 15, 35])
gn, sgn, mstar_eag, lambdaR_eag_eo, lambdaR_eag, ellip_eag_eo, ellip_eag, sigma_eag, r50_eag = np.loadtxt('LambdaGalaxiesz0_RefL100N1504.dat', unpack=True, usecols=[0,1,2,4,5,11,12,14,6])
#gn, sgn, mstar_eag, lambdaR_eag_eo, lambdaR_eag, ellip_eag_eo, ellip_eag, sigma_eag = np.loadtxt('EAGLE_L100N1504_ForLuca.dat', unpack=True, usecols=[0,1,2,12,13,8,9,6])


def calculate_galid(gng,sgng):
    gid, gn, sgn = np.loadtxt('/opt/claudia/EAGLE/Selection_Galaxies_Forensics_L0100N1504z0.dat', unpack=True, usecols=[0,2,3])
    galid = np.zeros(shape = len(gng))
    for i in range(0,len(gng)):
        match = np.where((gn == gng[i]) & (sgn == sgng[i]))
        galid[i] = gid[match]

    return galid



galidl = calculate_galid(gn,sgn)


if(AdoptProcedureSimSpin == True):
   ind = np.where(mstar_eag >= 1e10)
   gn = gn[ind]
   sgn = sgn[ind]
   mstar_eag = mstar_eag[ind]
   galidl = galidl[ind]
   lambdaR_eag_eo = lambdaR_eag_eo[ind]
   ellip_eag_eo = ellip_eag_eo[ind]
   
   def take_properties_simspin(galid, inclination=60):
       gid, inc, r50, e50,  l50 = np.loadtxt('Classification/keh_EAGLE_kinematics_ASGR_method.txt', unpack=True, usecols=[0,1,2,3,6])
       r50g = np.zeros(shape = len(galid))
       er50g = np.zeros(shape = len(galid))
       l50g = np.zeros(shape = len(galid))
   
       for i in range(0,len(galidl)):
           match = np.where((gid == galid[i]) & (inc == inclination))
           if(len(gid[match]) > 0):
              r50g[i] = r50[match]
              er50g[i] = e50[match]
              l50g[i] = l50[match]
           
       return (r50g, er50g, l50g)
  
   (r50_eag, ellip_eag, lambdaR_eag) = take_properties_simspin(galidl, inclination=60) 

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
xmin, xmax, ymin, ymax = 0, 0.6, 0, 7.5
ax = fig.add_subplot(211)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 0.1, 1, 1))

def classify_srs(l,e):
    classification = np.zeros(shape = len(l))
    sr = np.where((l < 0.08 + e/4.0) & (e <= 0.4))
    classification[sr] = 1
    sr = np.where((l < 0.18) & (e > 0.4))
    classification[sr] = 1
    return classification

sami_class = classify_srs(lambdaR_SAMI, ellip_SAMI)
sami_class_nocorr = classify_srs(lambdaR_SAMI_nocorr, ellip_SAMI)
ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <= 10.9) & (lambdaR_SAMI >= 0) & (sami_class > 0))
ax.hist(ellip_SAMI[ind], bins=20, range=(xmin,xmax), density=True, facecolor='blue', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k',label='SAMI')
ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <= 10.9) & (lambdaR_SAMI >= 0) & (sami_class_nocorr > 0))
ax.hist(ellip_SAMI[ind], bins=20, range=(xmin,xmax), density=True, edgecolor='blue', histtype='step', linewidth=2, alpha=0.85, label='SAMI no corr')

eagle_class = classify_srs(props_selected_eagle[1,:], props_selected_eagle[2,:])
ind = np.where((props_selected_eagle[0,:] >= 10) & (props_selected_eagle[0,:] <= 10.9) & (eagle_class > 0))
ein = props_selected_eagle[2,ind]
ax.hist(ein[0], bins=20, range=(xmin,xmax), density=True, facecolor='red', histtype='step', linewidth=2, alpha=0.25, fill=True, edgecolor='k', label='EAGLE')
common.prepare_legend(ax, ['b','b','r'], loc='upper left')

ax = fig.add_subplot(212)
xtit="$\\sigma_{\\rm r_{50}}$"
xmin, xmax, ymin, ymax = 50, 240, 0, 0.035
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(25, 25, 0.01, 0.01))

ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <= 10.9) & (lambdaR_SAMI >= 0) & (sami_class > 0))
ax.hist(sigma_SAMI[ind], bins=20, range=(50,240), density=True, facecolor='blue', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k',label='SAMI')
ind = np.where((mstar_SAMI >= 10) & (mstar_SAMI <= 10.9) & (lambdaR_SAMI >= 0) & (sami_class_nocorr > 0))
ax.hist(sigma_SAMI[ind], bins=20, range=(50,240), density=True, edgecolor='blue', histtype='step', linewidth=2, alpha=0.85, label='SAMI no corr')
ind = np.where((props_selected_eagle[0,:] >= 10) & (props_selected_eagle[0,:] <= 10.9) & (eagle_class > 0))
ein = props_selected_eagle[3,ind]
galidsmatched = props_selected_eagle[4,ind]
ax.hist(ein[0], bins=20, range=(50,240), density=True, facecolor='red', histtype='step', linewidth=2, alpha=0.25, fill=True, edgecolor='k', label='EAGLE')

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

def obtain_ellipticity_and_lambda(galid, inclination=60):
    props = np.zeros(shape = (3,len(galid)))
    props2 = np.zeros(shape = (5,len(galid)))

    gid, inc, r50, e50,  l50 = np.loadtxt('Classification/keh_EAGLE_kinematics_ASGR_method.txt', unpack=True, usecols=[0,1,2,3,6])
    for i in range(0,len(galid)):
        match = np.where((gid == galid[i]) & (inc == inclination))
        if(len(gid[match]) > 0):
           props[0,i] = l50[match]
           props[1,i] = e50[match]
           props[2,i] = r50[match]
    for i in range(0,len(galid)):
        match = np.where(galidl == galid[i])
        if(len(galidl[match]) > 0):
           props2[0,i] = lambdaR_eag_eo[match]
           props2[1,i] = ellip_eag_eo[match]
           props2[2,i] = r50_eag[match]
           props2[3,i] = gn[match]
           props2[4,i] = sgn[match]

    return (props, props2)

#read table with classification
galid, class1, class2, class3, class4, class5 = np.loadtxt('Classification/Classification.txt', unpack=True, usecols=[0,1,3,5,7,9])
commonids = np.in1d(galid, galidsmatched)
matched = np.where(commonids == True)

gnsr, sgnsr = np.loadtxt('SlowRotatorsEAGLE.txt', unpack=True, usecols=[0,1])
galidsr = calculate_galid(gnsr, sgnsr)
#find missing IDs
missing = 0
found = 0
ids_missing = galidsr
for i in range(0,len(galidsr)):
    match = np.where(galid == galidsr[i])
    if len(galid[match]) < 1:
       ids_missing[i] = galidsr[i]
       missing = missing + 1
    else:
       found = found + 1
       ids_missing[i] = 0


Nclassifiers = 5
all_classified_srs = np.in1d(galidl, galid)
#print(all_classified_srs)
#for a,b,c in zip(galid, ellip_eag_eo[all_classified_srs], ellip_eag[all_classified_srs]):
#    print(a,b,c)


classification = np.zeros(shape = (len(galid), Nclassifiers))
classification[:,0] = class1
classification[:,1] = class2
classification[:,2] = class3
classification[:,3] = class4
classification[:,4] = class5
classification = classification.astype(int)

(props, props2) = obtain_ellipticity_and_lambda(galid, inclination=60)
(props3, props2) = obtain_ellipticity_and_lambda(galid, inclination=75)

#select outliers:
ind = np.where((props[1,:] > 0.5) & (props2[1,:] < 0.2))
b =  props2[3,ind]
c = props2[4,ind]
d = props2[2,ind]
for a,b,c,d in zip(galid[ind], b[0], c[0], d[0]):
    print(a,b,c,d)

print(np.median(props2[2,:]))
ind = np.where(props2[1,:]/props[1,:] < 0.5)
print(np.median(props2[2,ind]))
ind = np.where(props2[1,:]/props[1,:] < 0.35)
print(np.median(props2[2,ind]))
ind = np.where(props2[1,:]/props[1,:] < 0.2)
print(np.median(props2[2,ind]))
ind = np.where(props2[1,:]/props[1,:] < 0.1)
print(np.median(props2[2,ind]))

ind = np.where(props2[1,:]/props[1,:] > 0.5)
print(np.median(props2[2,ind]))
ind = np.where(props2[1,:]/props[1,:] > 1)
print(np.median(props2[2,ind]))

#plot ellipticities and lambdaR for the two catalogues
xtit="$\\lambda_{\\rm r_{50}}$ edge-on claudia"
ytit="$\\lambda_{\\rm r_{50}}$ 75deg kate"
fig = plt.figure(figsize=(4,7))
plt.subplots_adjust(left=0.25, bottom=0.17)
xmin, xmax, ymin, ymax = 0,0.2,0,0.2
ax = fig.add_subplot(311)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 0.1, 0.1, 0.1))
ax.plot(props2[0,:], props3[0,:], 'ko', alpha=0.5)
#ax.hexbin(props2[0,:], props3[0,:], props2[2,:], xscale='linear', yscale='linear', gridsize=(30,30), cmap='Spectral', mincnt=0)
ax.plot([0,0.2],[0,0.2],linestyle='solid')
ax = fig.add_subplot(312)
xtit="$\\epsilon_{\\rm r_{50}}$ edge-on claudia"
ytit="$\\epsilon_{\\rm r_{50}}$ 75deg kate"
xmin, xmax, ymin, ymax = 0,1,0,1
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.2, 0.2, 0.2, 0.2))
ax.plot(props2[1,:], props3[1,:], 'ko', alpha=0.5)
#im = ax.hexbin(props2[1,:], props3[1,:], props2[2,:], xscale='linear', yscale='linear', gridsize=(20,20), cmap='Spectral', mincnt=0)
ax.plot([0,1],[0,1],linestyle='solid')

ax = fig.add_subplot(313)
xtit="${\\rm r_{50}}/\\rm kpc$ claudia"
ytit="${\\rm r_{50}}$ kate"
xmin, xmax, ymin, ymax = 0,10,0,10
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))
ax.plot(props2[2,:], props3[2,:], 'ko', alpha=0.5)
ax.plot([0,10],[0,10],linestyle='solid')
#cbar_ax = fig.add_axes([0.86, 0.4, 0.025, 0.25])
#cbar = fig.colorbar(im, cax=cbar_ax)
#cbar.ax.set_ylabel('$\\rm r_{\\rm 50}/kpc$')

obsdir = 'Classification/Plots/'
plt.tight_layout()
common.savefig(obsdir, fig, "TestEllipticitiesAndLambdaR.pdf")


#determine the number of galaxies that have three matches, two matches, or all different answers, use flags=2, 1, 0 to determine this
successmatching = np.zeros(shape = (len(galid)))
successmatching[:] = -1

for i in range(0,len(successmatching)):
    successmatching[i] = Nclassifiers - len(np.unique(classification[i,:]))

###############################################################################3
#plot success rate of matching procedure

xtit="Matching success"
ytit="N"
x_values = ['$0$', '$0.4$', '$0.6$', '$0.8$','$1$']
fig = plt.figure(figsize=(3,5))
plt.subplots_adjust(left=0.25, bottom=0.17)
xmin, xmax, ymin, ymax = -0.5, 4.5, 0, 105
ax = fig.add_subplot(111)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 50, 50))
ax.hist(successmatching, bins=5, range=(-0.5,4.5), facecolor='blue', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')

plt.xticks([0,1,2,3,4], x_values)
common.savefig(obsdir, fig, "MatchesInClass.pdf")


#################################################################
#plot distribution of kinematic classes for the three match success %

xtit="Kinematic class"
ytit="PDF"
x_values = ['FSR', 'RSR', '2$\sigma$', 'Prol','Uncl','R']
fig = plt.figure(figsize=(3,5))
plt.subplots_adjust(left=0.2, bottom=0.15)
xmin, xmax, ymin, ymax = -0.5, 2.5, 0, 1
colors = ['orange', 'red', 'salmon', 'Crimson','DarkRed']
def plot_kin_class(ax, match=0):
    ind = np.where(successmatching[:] == match)
    for i in range(0,Nclassifiers):
        y = classification[ind, i]
        ax.hist(y[0], bins=6, range=(-0.5,5.5), density=True, facecolor=colors[i], histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')
    #ax.hist(class2[ind], bins=6, range=(-0.5,5.5), density=True, facecolor='red', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')
    #ax.hist(class3[ind], bins=6, range=(-0.5,5.5), density=True, facecolor='salmon', histtype='step', linewidth=2, alpha=0.5, fill=True, edgecolor='k')

subplots = (411, 412, 413, 414)
matches = (1, 2, 3, 4)
text_insert = ['$40\%$', '$60\%$', '$80\%$', '$100\%$']
xpos, ypos = xmax - 0.75 * (xmax - xmin)/2.0, ymax - 0.15 * (ymax - ymin)
for i in range(0,len(matches)):
    ax = fig.add_subplot(subplots[i])
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 50, 50))
    ax.text(xpos, ypos, 'Match succ = %s' % text_insert[i])
    plot_kin_class(ax, match=matches[i])
    if(i < len(matches)-1):
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
    plot_elip_class(ax, class4, match=matches[i], color='crimson')
    plot_elip_class(ax, class5, match=matches[i], color='DarkRed')

obsdir = 'Classification/Plots/'
common.savefig(obsdir, fig, "EllipticityInClass.pdf")


###################################################################################
###################################################################################
###################################################################################
# Analysis of merger history of different kinematic classes

mergerdata = np.loadtxt('MergerHistoryCumulative.dat')
def find_classification(selective = True):
    if(selective):
       matchedobjs = np.where(successmatching > 1)
       classificationin = classification[matchedobjs,:]
       classificationin = classificationin[0]
       finalclass = np.zeros(shape = len(galid[matchedobjs]))
       finalgalid = galid[matchedobjs]
       for i in range(0, len(galid[matchedobjs])):
           counts = np.bincount(classificationin[i,:])
           finalclass[i] = np.argmax(counts)
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

