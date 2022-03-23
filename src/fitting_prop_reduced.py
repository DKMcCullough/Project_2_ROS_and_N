from numpy import *
from scipy import *
from pylab import *
from read_dthief import *
from scipy.integrate import *
import matplotlib as plt 
from plotnine import *


def get_Qn():
        return 9.4e-15*14.0*1e+6 # mumol cell-1

def f(u,t,mumax,ks,nrk,ndk):
        Qn = get_Qn()
        N,P = u[0],u[1]
        if (N > 1e-4):
                df = nrk
        else:
                df = ndk
        dPdt = mumax*N/(N+ks)*P - df*P
        dNdt = -mumax*N/(N+ks)*P*Qn
        return concatenate([r_[[dNdt]],r_[[dPdt]]])

def integrate(N,P,ntimes,ptimes,pars,delt=1.0/24.0,forshow=False):
        days = amax(ptimes)
        times = linspace(0,days,int(days/delt))
        u = odeint(f,r_[[N,P]],times,args=(exp(pars[0]),exp(pars[1]),exp(pars[2]),exp(pars[3]))).T
        if forshow==False:
                ninds = r_[[where(abs(a-times)==min(abs(a-times)))[0][0] for a in ntimes]] # get time indices where host abundances were measured
                pinds = r_[[where(abs(a-times)==min(abs(a-times)))[0][0] for a in ptimes]] # same for viruses
                nnt = u[0][ninds] 
                pnt = u[1][pinds] 
        else:
                nnt = u[0]
                pnt = u[1]
        return nnt,pnt

def get_chi(nts,pts,ncs,pcs,pstds,pars):
        chisum = 0
        mumax,ks,nrks,ndks = pars[0],pars[1],pars[2:7],pars[7:12]
        for (nt,pt,nc,pc,pstd,nrk,ndk) in zip(nts,pts,ncs,pcs,pstds,nrks,ndks):
                parsloc = r_[[mumax,ks,nrk,ndk]]
                nnt,pnt = integrate(nc[0],pc[0],nt,pt,parsloc)	
                chi = sum((pnt - pc) ** 2 / (pstd ** 2))
                chisum = chisum + chi
        return chisum

ddir = './morris/data/'




plt.rc('xtick', color='k')
plt.rc('lines', lw=2.0)
plt.rc('grid', lw=4)
rc('font', weight='bold',size=16)





pt1,pc1 = read_dthief(ddir+'UH18301_005HOOH.txt')
pt2,pc2 = read_dthief(ddir+'UH18301_02HOOH.txt')
pt3,pc3 = read_dthief(ddir+'UH18301_04HOOH.txt')
pt4,pc4 = read_dthief(ddir+'UH18301_08HOOH.txt')
pt5,pc5 = read_dthief(ddir+'UH18301_10HOOH.txt')
pc1,pc2,pc3,pc4,pc5 = 10**pc1*1.01,10**pc2*0.99,10**pc3,10**pc4*0.999,10**pc5

nc1 = r_[[(pc1[-2]-pc1[0])*get_Qn()]]
nc2 = r_[[(pc1[-2]-pc1[0])*get_Qn()]]
nc3 = r_[[(pc1[-2]-pc1[0])*get_Qn()]]
nc4 = r_[[(pc1[-2]-pc1[0])*get_Qn()]]
nc5 = r_[[(pc1[-2]-pc1[0])*get_Qn()]]

nt1 = r_[[pt1[0]]]
nt2 = r_[[pt2[0]]]
nt3 = r_[[pt3[0]]]
nt4 = r_[[pt4[0]]]
nt5 = r_[[pt5[0]]]

pts,pcs = r_[[pt1,pt2,pt3,pt4,pt5]]/24.0,r_[[pc1,pc2,pc3,pc4,pc5]]
nts,ncs = r_[[nt1,nt2,nt3,nt4,nt5]],r_[[nc1,nc2,nc3,nc4,nc5]]
pstds = r_[[0.1*mean(pc) for pc in pcs]]
fa,ax = subplots(1,2,figsize=[12,6])
f1,ax1 = subplots(figsize=[9,6])

ax[1].errorbar(pt1,pc1/1e+5,label=r'0.05 $\mu$M HOOH', fmt='o',yerr=pstds[0]/1e+5)
ax[1].errorbar(pt2,pc2/1e+5,label=r'0.2 $\mu$M HOOH',fmt='o',yerr=pstds[1]/1e+5)
ax[1].errorbar(pt3,pc3/1e+5,label=r'0.4 $\mu$M HOOH',fmt='o',yerr=pstds[2]/1e+5)
#ax[1].errorbar(pt4,pc4/1e+5,label=r'0.8 $\mu$M HOOH',fmt='o',yerr=pstds[3]/1e+5)
#ax[1].errorbar(pt5,pc5/1e+5,label=r'10 $\mu$M HOOH',fmt='o',yerr=pstds[4]/1e+5)  

l = ax[1].legend(loc='upper left',fontsize=12,ncol=2)
l.draw_frame(False)

HOOHs = r_[[0.05,0.2,0.4,0.8,10.0]]

mumax = 0.3368234145366085
ks = 4.8350570572782885e-09
kdam1 = 0.06736321395954317
kdam2 = 0.19474866062938384
kdam3 = 0.6925560296486394
kdam4 = 1.4097180432026337
kdam5 = 3.947694942792268
kddam1 = 0.030527209915541417
kddam2 = 0.5029222415746297
kddam3 = 4.408103330893407
kddam4 = 0.7929504739796394
kddam5 = 1.3663649550894377

pars = log(r_[[mumax,ks,kdam1,kdam2,kdam3,kdam4,kdam5,kddam1,kddam2,kddam3,kddam4,kddam5]])
npars = pars.shape[0]

chi = get_chi(nts,pts,ncs,pcs,pstds,pars)

stds = zeros(npars) + .05 # this controls how large the random steps are in the parameter search (see below)
opt = r_[[1,1,1,1,1,1,1,1,1,1,1,1]] # this allows you to control whether each parameter is imposed or fitted
names = ['mumax','ks','kdam1','kdam2','kdam3','kdam4','kdam5','kddam1','kddam2','kddam3','kddam4','kddam5'] # name each parameter array - useful for printing the output

nits = 500 # number of iterations
pits = 100  # frequency with which to print results to the user about the progress of the algorithm
burnin = 200 

# distribution arrays and acceptance ratios - these are containers to be added to
ar = 0.0
ars = r_[[]]
mumaxs,kss,kdam1s,kdam2s,kdam3s,kdam4s,kdam5s,kddam1s,kddam2s,kddam3s,kddam4s,kddam5s = r_[[]],r_[[]],r_[[]],r_[[]],r_[[]],r_[[]],r_[[]],r_[[]],r_[[]],r_[[]],r_[[]],r_[[]]
pall = r_[[mumaxs,kss,kdam1s,kdam2s,kdam3s,kdam4s,kdam5s,kddam1s,kddam2s,kddam3s,kddam4s,kddam5s]]

# now actually do the fitting
for it in arange(1,nits,1):
        parsnew = pars + opt*normal(0,stds,npars) # this is where we randomly change the parameter values 
        lmax = amax(pars[3:7])
        parsnew = r_[[min(p,log(4.0)) for p in parsnew]]
        chinew = get_chi(nts,pts,ncs,pcs,pstds,parsnew) 
        if exp(chi-chinew) > rand(): # KEY STEP
                chi = chinew
                pars = parsnew #  new parameters can be a little bit 'wrong'
                if it > burnin: # only store the parameters if you've gone through the burnin period
                        pall = append(pall,pars[:,None],1)
                ar = ar + 1.0 # acceptance ratio
        if (it % pits == 0):
                #print (it,chi,ar/pits)
                ars = append(ars,ar/pits)
                ar = 0.0
'''
# print output to screen
print ('Optimal parameters for high temp conditions')
pars = r_[[ mean(p) for p in pall]]
for (p,l) in zip(pars,names):
        print (l,'=',exp(p))

print (' ')
print ('Standard deviations')
for (p,l) in zip(pall,names):
        print (l+'std','=',std(exp(p)))
print (' ')
'''

# redefine times for nicer looking plots
# run again just for nicer looking plots (more even timesteps)
mumax,ks,nrks = log(0.4498586),log(5.07530e-09),log(r_[[0.053584,0.333407,0.92427,1.482875,3.93573]])
ndks = nrks
ax[1].set_prop_cycle(None)
count=0
for (nt,pt,nc,pc,pstd,nr,nd) in zip(nts,pts,ncs,pcs,pstds,nrks,ndks):
        parsloc = r_[[mumax,ks,nr,nd]]
        delt = 1.0/24.0
        nnt,pnt = integrate(nc[0],pc[0],nt,pt,parsloc,forshow=True)
        mtimes = linspace(0,amax(pt),(amax(pt) /delt ))
        ax[0].plot(mtimes*24,nnt, markersize=9)
        ax[1].plot(mtimes*24,pnt/1e+5, markersize=9)
        count=count+1
        #print(count) 
        if count == 3:
                break

ax[0].set_prop_cycle(None)
ax[1].set_prop_cycle(None)
count=0
mumax,ks,nrks,ndks = pars[0],pars[1],pars[2:7],pars[7:12]
for (nt,pt,nc,pc,pstd,nr,nd) in zip(nts,pts,ncs,pcs,pstds,nrks,ndks):
        parsloc = r_[[mumax,ks,nr,nd]]
        delt = 1.0/24.0
        nnt,pnt = integrate(nc[0],pc[0],nt,pt,parsloc,forshow=True)
        mtimes = linspace(0,amax(pt),(amax(pt) /delt ))
        ax[0].plot(mtimes*24,nnt,ls='--')
        ax[1].plot(mtimes*24,pnt/1e+5,ls='--')
        count=count+1
        #print(count) 
        if count == 3:
                break
#HOOHs = r_[[0.05,0.2,0.4,0.8,10.0]]
#nrks = 0.053584,0.333407,0.92427,1.482875,3.93573


ax1.errorbar(0.05,0.053584,c='b',marker='o',fmt='o',markersize=10,label='Damage rate, $\kappa_{dam}$')
ax1.errorbar(0.2,0.333407,c='orange',marker='o',fmt='o',markersize=10,label='Damage rate, $\kappa_{dam}$')
ax1.errorbar(0.4,0.92427,c='g',marker='o',fmt='o',markersize=10,label='Damage rate, $\kappa_{dam}$')
#ax1.errorbar(0.8,1.482875,c='r',marker='o',fmt='o',markersize=10,label='Damage rate, $\kappa_{dam}$')
#ax1.errorbar(10.0,3.93573,c='purple',marker='o',fmt='o',markersize=10,label='Damage rate, $\kappa_{dam}$')


ax1.errorbar(0.05,0.030527209915541417,c='b',marker='o',fmt='o',markersize=10,label='Deplete N Damage rate, $\kappa_{dam,d}$')
ax1.errorbar(0.05,0.5029222415746297,c='orange',marker='o',fmt='o',markersize=10,label='Deplete N Damage rate, $\kappa_{dam,d}$')
ax1.errorbar(0.05,4.408103330893407,c='g',marker='o',fmt='o',markersize=10,label='Deplete N Damage rate, $\kappa_{dam,d}$')





#ax1.errorbar(HOOHs,exp(nrks),yerr=std(exp(pall)[2:7],1),c='k',marker='o',fmt='o',markersize=10,label='High N damage rate, $\kappa_{dam,r}$')  #High N damage rate, $\kappa_{dam,r}
#ax1.errorbar(HOOHs[:3],exp(ndks[:3]),yerr=std(exp(pall)[2:7],1)[:3],c='k',marker='*',fmt='o',markersize=10,label='Low N damage rate $\kappa_{dam,d}$')
mh = linspace(0,10,1000)
mhd = linspace(0,1,1000)
dmax,alp = 4.6,3.0
md = dmax*mh/(mh+(dmax/alp))
mdd = dmax*mhd/(mhd+(dmax/(alp*4)))
ax1.plot(mh,md,c='k')
ax1.plot(mhd,mdd,c='k',ls='--')
ax[1].set_ylim(ymax=2.5)
l = ax1.legend(loc='lower right', fontsize=16)
l.draw_frame(False)


ax[0].set_xlabel('Time (hours)', fontweight='bold', fontsize=20)
ax[1].set_xlabel('Time (hours)', fontweight='bold', fontsize=20)
ax[0].set_ylabel('Nitrogen Concentration ($\mu$mol ml$^{-1}$)',fontweight='bold',  fontsize=20)
ax[1].set_ylabel('Cells (x10$^{5}$ ml$^{-1}$)',fontweight='bold', fontsize=20)
ax1.set_xlabel(r'HOOH ($\mu$M)', fontweight='bold', fontsize=20)
ax1.set_ylabel(r'Damage rate (day$^{-1}$)', fontweight='bold',fontsize=20)

fa.subplots_adjust(wspace=0.3)

#fa.savefig('figures/tseries')
#f1.savefig('figures/damage')

show()
