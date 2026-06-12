from scipy.stats import norm
from statistics import mode
import seaborn as sns
import  pickle
import matplotlib.mlab as mlab
import illustris_python as il
import numpy as np
from scipy.stats import binned_statistic_2d
import matplotlib.pyplot as plt
import math
from collections import namedtuple
#from scipy.spatial.transform import Rotation as R
import matplotlib.mlab as mlab
from scipy.interpolate import interp1d



def Distance_1D(X, X_POS, BoxSize):
    '''This function takes as input a 1D vector containing the positions of particles, the X_POS(float) that is the position respect to where we will compute the Distance and the BoxSize, this function consider a periodical Box. The output is a 1D vector with the same size than X containing the distance from X to X_POS'''
    D=X-X_POS
    D=np.where(D>BoxSize/2, D-BoxSize, D)
    D=np.where(D<-BoxSize/2, D+BoxSize, D)
    return D  


def CM_1D(X, M, X_POS, BoxSize, BoxType=0):
    D=Distance_1D(X,X_POS, BoxSize)
    CM= np.inner(D,M)/np.sum(M)
    CM=CM+X_POS
    if BoxType==1: 
        if CM > BoxSize/2: CM=CM-BoxSize 
        if CM <-BoxSize/2: CM=CM+BoxSize
    elif BoxType==0:
        if CM > BoxSize: CM=CM-BoxSize 
        if CM <0: CM=CM+BoxSize
    else:
        print("El BoxType= 0 para caja de 0 a L y BoxType=1 para caja de --L/2 a L/2")
    return CM 


def CM_2D(R,M,R_POS,BoxSize, BoxType=0):
    try: return np.array([CM_1D(R[:,0],M, R_POS[0],BoxSize, BoxType),CM_1D(R[:,1],M, R_POS[1],BoxSize, BoxType) ])
    except: print("El vector debe contener las posiciones de al menos dos partículas. ")
    
    
def Distance_2D(R,R_POS,BoxSize, BoxType=0):
    '''This function takes as input a 2D vector containing the positions of particles, the X_POS(2D point) that is the position respect to where we will compute the Distance and the BoxSize, this function consider a periodical Box. The output is a 2D vector with the same size than X containing the euclidean distance from X to X_POS'''
    try: return np.transpose(np.array([Distance_1D(R[:,0], R_POS[0],BoxSize),Distance_1D(R[:,1], R_POS[1],BoxSize) ]))
    except: return np.transpose(np.array([Distance_1D(R[0], R_POS[0],BoxSize),Distance_1D(R[1], R_POS[1],BoxSize) ]))


def CM_3D(R,M,R_POS,BoxSize, BoxType=0):
    try: return np.array([CM_1D(R[:,0],M, R_POS[0],BoxSize, BoxType),CM_1D(R[:,1],M, R_POS[1],BoxSize, BoxType), CM_1D(R[:,2],M, R_POS[2],BoxSize, BoxType) ])
    except: print("El vector debe contener las posiciones de al menos dos partículas. ")
    
    
def Distance_3D(R,R_POS,BoxSize):
    '''This function takes as input a 3D vector containing the positions of particles, the X_POS(3D point) that is the position respect to where we will compute the Distance and the BoxSize, this function consider a periodical Box. The output is a 3D vector with the same size than X containing the euclidean distance from X to X_POS'''
    try: return np.transpose(np.array([Distance_1D(R[:,0], R_POS[0],BoxSize),Distance_1D(R[:,1], R_POS[1],BoxSize), Distance_1D(R[:,2], R_POS[2],BoxSize) ]))
    except: return np.transpose(np.array([Distance_1D(R[0], R_POS[0],BoxSize),Distance_1D(R[1], R_POS[1],BoxSize), Distance_1D(R[2], R_POS[2],BoxSize) ]))


def Distance_3D_vectors(vector1, vector2, BoxSize):
    ''' This function takes as input two 3D vectors of the same size and computes the 3D distance between those two vectors component by component. The output of this function is a 3D array containing the distance between those points '''
    if np.size(vector1)==np.size(vector2):
        for i in range(int(np.size(vector1)/3)): 
            try: aux_offset= np.append(aux_offset, [Distance_3D(vector1[i], vector2[i], BoxSize)], axis=0)
            except: aux_offset=Distance_3D(vector1[i], vector2[i], BoxSize).reshape(1,3)
    return aux_offset


def Comp_CM_PartType(basePath, N, Header, Halos_cat, selected_halos,partType='gas'):
    global aux_Halo_cat
    i=0
    for NHalo in selected_halos:
        print(f'Computing the center of mass of {partType} from the halo number: {NHalo}')
        halo= il.snapshot.loadHalo(basePath, N, NHalo, partType=partType,fields=['Coordinates','Masses'])
        
        SZ=CM_3D(halo['Coordinates'],halo['Masses'],Halos_cat['GroupPos'][i,:],Header['BoxSize'])
        try: aux_Halo_cat['GroupSZPos']=np.append(aux_Halo_cat['GroupSZPos'],[SZ],axis=0)
        except: aux_Halo_cat['GroupSZPos']=SZ.reshape(1,3)
        print("El tamaño del array del SZ es ", np.size(aux_Halo_cat['GroupSZPos'])/3)
        i+=1
    print("We finished the omputing of CM")


def plot_functiont(PartPos, Weights, header):
    x=PartPos[:,0]
    y=PartPos[:,1]
    a=np.max([np.abs(np.min(x)), np.abs(np.max(x)), np.abs(np.min(y)), np.abs(np.max(y)) ])
    nPixels =[900, 900]
    minMax =[-a,a]
    grid, _, _, _ = binned_statistic_2d(x, y, Weights, 'sum', bins=nPixels, range=[minMax,minMax])
    pxSize = 2*a / nPixels[0] # code units
    pxSize_kpc = pxSize * header['Time'] / header['HubbleParam']
    pxArea = pxSize_kpc**2
    grid_log_msun_per_kpc2 = np.log10(grid.T * 1e10 / header['HubbleParam']) #/ pxArea)
    fig = plt.figure(figsize=(16,10))
    ax = fig.add_subplot(111)
    extent = 0.5*np.array([-a, a, -a, a])
    #ax.imshow(grid_log_msun_per_kpc2.T, extent=extent, cmap='rainbow', aspect=nPixels[1]/nPixels[0])#,vmin=5,vmax=8)
    ax.contourf(grid_log_msun_per_kpc2, cmap='rainbow')#,vmin=5,vmax=8)
    ax.contour(grid_log_msun_per_kpc2,cmap='rainbow')#,vmin=5,vmax=8)
    #plt.colorbar()
    #sns.kdeplot(grid_log_msun_per_kpc2, shade=True,ax=ax)
    plt.show()
#######################################################################################################################################################################################################################################################################################################
def plot_function(x,y, weights,header,x_cm,y_cm,r200,x_BCG,y_BCG,x_SZ,y_SZ):
    #a=math.ceil(np.max([np.abs(np.min(x)), np.abs(np.max(x)), np.abs(np.min(y)), np.abs(np.max(y)) ])/1000)*1000
    a=np.max([np.abs(np.min(x)), np.abs(np.max(x)), np.abs(np.min(y)), np.abs(np.max(y)) ])
    nPixels = [600,600]
    minMax = [-a, a]
    #header = il.groupcat.loadHeader(basePath, N)
    print(np.shape(x), np.shape(y), np.shape(weights))
    grid, _, _, _ = binned_statistic_2d(x, y, weights, 'sum', bins=nPixels, range=[minMax,minMax])
    pxSize = 2*a / nPixels[0] # code units
    pxSize_kpc = pxSize * header['Time'] / header['HubbleParam']
    pxArea = pxSize_kpc**2
    grid_log_msun_per_kpc2 = np.log10(grid * 1e10 / header['HubbleParam']) #/ pxArea)
    fig = plt.figure(figsize=(16,10))
    ax = fig.add_subplot(111)
    extent = [-a, a, -a, a]
    plt.imshow(grid_log_msun_per_kpc2, extent=extent, cmap='rainbow', aspect=nPixels[1]/nPixels[0])#,vmin=5,vmax=8)
    print("Corre bien")
   # ax.autoscale(False)
   # ax.set_xlabel('y [ckpc/h]')
   # ax.set_ylabel('z [ckpc/h]')
   # plt.colorbar(label=r'$masstype Density [log M$_\odot$ kpc$^{-3}$]');
   # plt.plot(y_cm,-x_cm,"rX", markersize='15',label="Gas CM")
   # plt.plot(y_SZ,-x_SZ,"g+", markersize='15',label="Density peak")
   # xc = np.linspace(-r200,r200,1000)
   # yc = np.sqrt(-xc**2+r200**2)
   # plt.plot(xc, yc,'k')
    #plt.plot(xc,-yc,'k')
    #plt.plot(y_BCG,-x_BCG,"b*", markersize='15',label="BCG position")
    plt.legend() 
    plt.show()
    return fig
########################################################################################################################################################################################################################################################################
def dictionary_select_items(x, ind):
    '''This function needs as input a dictionary x and an array of indices ind
    and returns a new dictionary with all fields, but just with the selected data'''
    keys_to_update=[i for i in x.keys() -{'count'}]
    z=x
    for i in keys_to_update:
        z[i]=x[i][ind]
    z['count']= np.size(np.where(ind))
    return z


def join_dictionaries(x,y):
    ''' This function joints two dictionaries x and y and joints the data of y after the data of x dictionaries, to joint those two dictionaries is it mandatory that x and y have the same keywords '''
    if x.keys()==y.keys():
        z=x
        keys_to_update=[i for i in x.keys() -{'count'}]
        for i in keys_to_update:
            z[i]=np.concatenate((x[i], y[i]), axis=0)
        z['count']= np.int(x['count'])+np.int(y['count'])
    else:
        print('To joint two dictionaries is mandatory they have the same keys.')
    return z


def TNG2PROPnum(Halos_cats, Subhalos_cats,snap):
    ''' THis function takes the numeration of halos in TNG Illustris simulations and update then of stack information of different redshifts adding 1e10 times the snapshot number to the halo r subhalo number '''
    Halos_cats['GroupNumber']= np.array([snap*1e10+i for i in range(int(Halos_cats['count']))]).astype(int)
    Subhalos_cats['SubhaloNumber']= np.array([snap*1e10+i for i in range(int(Subhalos_cats['count']))]).astype(int)
    Subhalos_cats['SubhaloGrNr']= (Subhalos_cats['SubhaloGrNr']+snap*1e10).astype(int)
    Halos_cats['GroupFirstSub']= (Halos_cats['GroupFirstSub']+snap*1e10).astype(int)
    return Halos_cats, Subhalos_cats


def select_halos(Halos_cat, Halos_limits):
    for i in Halos_limits.keys():
        aux_key=i.split('/')
        if np.size(aux_key)==1:
            if 'selected_halos' in locals():
                selected_halos=np.logical_and(np.logical_and(selected_halos,Halos_cat[i]>=Halos_limits[i][0]), Halos_cat[i]<=Halos_limits[i][1])
            else:
                selected_halos=np.logical_and(Halos_cat[i]>=Halos_limits[i][0], Halos_cat[i]<=Halos_limits[i][1])
        else:
            if 'selected_halos' in locals():
                selected_halos=np.logical_and(np.logical_and(selected_halos,Halos_cat[aux_key[0]][:,int(aux_key[1])]>=Halos_limits[i][0]), Halos_cat[aux_key[0]][:,int(aux_key[1])]<=Halos_limits[i][1])
            else:
                selected_halos=np.logical_and(Halos_cat[aux_key[0]][:,int(aux_key[1])]>=Halos_limits[i][0], Halos_cat[aux_key[0]][:,int(aux_key[1])]<=Halos_limits[i][1])
    return selected_halos


def select_subhalos(Subhalos_cat, Subhalos_limits, selected_halos):
    for i in Subhalos_limits.keys():
        aux_key=i.split('/')
        if np.size(aux_key)==1:
            if 'selected_subhalos' in locals():
                selected_subhalos= np.logical_and(np.logical_and(selected_subhalos, Subhalos_cat[i]>=Subhalos_limits[i][0]), Subhalos_cat[i]<=Subhalos_limits[i][1])
            else:
                selected_subhalos=np.logical_and(np.logical_and(np.logical_and(Subhalos_cat['SubhaloFlag'],np.in1d(Subhalos_cat['SubhaloGrNr'], np.where(selected_halos)[0])),Subhalos_cat[i]>=Subhalos_limits[i][0]), Subhalos_cat[i]<=Subhalos_limits[i][1])
        else:
            if 'selected_subhalos' in locals():
                selected_subhalos= np.logical_and(np.logical_and(selected_subhalos, Subhalos_cat[aux_key[0]][:,int(aux_key[1])]>=Subhalos_limits[i][0]), Subhalos_cat[aux_key[0]][:,int(aux_key[1])]<=Subhalos_limits[i][1])
            else:
                selected_subhalos=np.logical_and(np.logical_and(np.logical_and(Subhalos_cat['SubhaloFlag'],np.in1d(Subhalos_cat['SubhaloGrNr'], np.where(selected_halos)[0])),Subhalos_cat[aux_key[0]][:,int(aux_key[1])]>=Subhalos_limits[i][0]), Subhalos_cat[aux_key[0]][:,int(aux_key[1])]<=Subhalos_limits[i][1])
    return selected_subhalos



def add_center_mass(aux_Halo_cat, basePath, N, Header, Halos_cat, selected_halos, partType):
    i=0
    print(f'Computing the center of mass of the gas in the snapshot {N}')
    for NHalo in np.where(selected_halos)[0]:
        #print(f'Computing the center of mass of gas from the halo number:   {NHalo}')
        halo_gas= il.snapshot.loadHalo(basePath, N, NHalo, partType=partType,fields=['Coordinates','Masses'])
        SZ=CM_3D(halo_gas['Coordinates'],halo_gas['Masses'],Halos_cat['GroupPos'][i,:],Header['BoxSize'])
        if 'GroupSZPos'  in aux_Halo_cat.keys():
            #print("Se crea la keyword GroupSZPos")
            aux_Halo_cat['GroupSZPos']=np.append(aux_Halo_cat['GroupSZPos'],[SZ],axis=0)
        else:
            #print("Se crea la keyword GroupSZ")
            aux_Halo_cat['GroupSZPos']=SZ.reshape(1,3)
        #print("El tamaño del array del SZ es ", np.size(aux_Halo_cat['GroupSZPos'])/3)
        i+=1
    #print("This Snapshot have been finished")
    return aux_Halo_cat


def add_center_mass_DM(aux_Halo_cat, basePath, N, Header, Halos_cat, selected_halos, partType='dm'):
    i=0
    print(f'Computing the center of mass of the dark matter in the snapshot {N}')
    for NHalo in np.where(selected_halos)[0]:
        #print(f'Computing the center of mass of gas from the halo number:   {NHalo}')
        halo_dm= il.snapshot.loadHalo(basePath, N, NHalo, partType=partType,fields=['Coordinates'])
        masses_dm=np.ones(int(np.shape(halo_dm)[0]))#As all DM particles have the same mass we put as 1
        DM_CM=CM_3D(halo_dm,masses_dm,Halos_cat['GroupPos'][i,:],Header['BoxSize'])
        #print(DM_CM)
        if 'Group_DM_CM'  in aux_Halo_cat.keys():
            #print("Se crea la keyword GroupSZPos")
            aux_Halo_cat['Group_DM_CM']=np.append(aux_Halo_cat['Group_DM_CM'],[DM_CM],axis=0)
        else:
            #print("Se crea la keyword GroupSZ")
            aux_Halo_cat['Group_DM_CM']=DM_CM.reshape(1,3)
        #print("El tamaño del array del SZ es ", np.size(aux_Halo_cat['GroupSZPos'])/3)
        i+=1
    #print("This Snapshot have been finished")
    return aux_Halo_cat


def add_SubhalosDist(Subhalos_cat, Halos_cat, Header):
    Subhalos_cat['SubhaloRelPos']=Subhalos_cat['SubhaloPos'].copy()
    Subhalos_cat['SubhaloDist']=Subhalos_cat['SubhaloSFR'].copy()
    Subhalos_cat['SubhaloDistR200']=Subhalos_cat['SubhaloSFR'].copy()
    for hal in Halos_cat['GroupNumber']:
        substructures=np.where(Subhalos_cat['SubhaloGrNr']==hal)[0]
        Subhalos_cat['SubhaloRelPos'][substructures]=Distance_3D(Subhalos_cat['SubhaloPos'][substructures,:], Halos_cat['GroupPos'][Halos_cat['GroupNumber']==hal].flatten(), Header['BoxSize'])
        Subhalos_cat['SubhaloDistR200'][substructures]=np.linalg.norm(Subhalos_cat['SubhaloRelPos'][substructures],axis=1)/Halos_cat['Group_R_Crit200'][Halos_cat['GroupNumber']==hal].flatten()
    Subhalos_cat['SubhaloDist']=np.linalg.norm(Subhalos_cat['SubhaloRelPos'],axis=1)
    print(hal)
    return Subhalos_cat


def beta_mod(z):
    betas=[-8.15,-7.72,-7.64]
    redshifts=[0,0.75,1]
    alpha_quench_mod=interp1d(redshifts, betas)
    return(alpha_quench_mod(z))


def quench_value(LogMStellar, z):
    alpha=0.8
    return((alpha-1)*LogMStellar+(beta_mod(z)-1))


def get_redshift(GrNr):
    aux=int(str([GrNr])[1:3])
    if aux==99: z=0
    elif aux==91: z=0.1
    elif aux==84: z=0.2
    elif aux==78: z=0.3
    elif aux==72: z=0.4
    elif aux==67: z=0.5
    elif aux==59: z=0.7
    elif aux==50: z=1
    return z


def add_SubhalosDist2(Subhalos_cat, Halos_cat, Header):
    Subhalos_cat['SubhaloRelPos']=np.empty(np.shape(Subhalos_cat['SubhaloPos'].copy()))
    Subhalos_cat['SubhaloDist']=np.empty(np.shape(Subhalos_cat['SubhaloSFR'].copy()))
    Subhalos_cat['SubhaloDistR200']=np.empty(np.shape(Subhalos_cat['SubhaloSFR'].copy()))
    Halos_cat['GroupSigmaR200']=np.empty(np.shape(Halos_cat['GroupFirstSub'].copy()))
    Subhalos_cat['SubhaloVel3DSigma']=np.empty(np.shape(Subhalos_cat['SubhaloSFR'].copy()))
    Subhalos_cat['SubhaloVel3D']=np.empty(np.shape(Subhalos_cat['SubhaloPos'].copy()))
    for hal in Halos_cat['GroupNumber']:
        #print(hal)
        substructures=Subhalos_cat['SubhaloGrNr']==hal
        Subhalos_cat['SubhaloRelPos'][substructures]=Distance_3D(Subhalos_cat['SubhaloPos'][substructures,:], Halos_cat['GroupPos'][Halos_cat['GroupNumber']==hal].flatten(), Header['BoxSize'])
        Subhalos_cat['SubhaloDistR200'][substructures]=np.linalg.norm(Subhalos_cat['SubhaloRelPos'][substructures],axis=1)/Halos_cat['Group_R_Crit200'][Halos_cat['GroupNumber']==hal].flatten()
        substructuresR200=np.logical_and(Subhalos_cat['SubhaloGrNr']==hal, Subhalos_cat['SubhaloDistR200']<=1 )
        #Subhalos_cat['SubhaloVel3D'][substructures]=np.linalg.norm(Subhalos_cat['SubhaloVel'][substructures,:]-np.mean(Subhalos_cat['SubhaloVel'][substructuresR200],axis=0), axis=1)
        Subhalos_cat['SubhaloVel3D'][substructures,:]=Subhalos_cat['SubhaloVel'][substructures,:]-np.mean(Subhalos_cat['SubhaloVel'][substructures,:],axis=0).flatten()
        #(mu_s, sigma_s) = norm.fit(Subhalos_cat['SubhaloVel3D'][substructures])
        sigma_s=np.sqrt(np.sum(np.square(Subhalos_cat['SubhaloVel3D'][substructures,:]))/np.size(Subhalos_cat['SubhaloVel3D'][substructures,:]))
        Halos_cat['GroupSigmaR200'][Halos_cat['GroupNumber']==hal]=sigma_s
        Subhalos_cat['SubhaloVel3DSigma'][substructures]=np.linalg.norm(Subhalos_cat['SubhaloVel3D'][substructures, :], axis=1)/sigma_s
    Subhalos_cat['SubhaloDist']=np.linalg.norm(Subhalos_cat['SubhaloRelPos'],axis=1)
    return Subhalos_cat, Halos_cat


def comp_Inertia_tensor(particles,SubhaloPos,Header, DistanceLIM):
    Dist= Distance_3D(particles['Coordinates'],SubhaloPos, Header['BoxSize'])
    #print(DistanceLIM)
    IDSP=np.linalg.norm(Dist,axis=1)<=float(DistanceLIM)
    #print(IDSP)
    Ixx=np.sum(particles['Masses'][IDSP]*(np.square(Dist[:,1][IDSP])+ np.square(Dist[:,2][IDSP])))
    Iyy=np.sum(particles['Masses'][IDSP]*(np.square(Dist[:,0][IDSP])+ np.square(Dist[:,2][IDSP])))
    Izz=np.sum(particles['Masses'][IDSP]*(np.square(Dist[:,0][IDSP])+ np.square(Dist[:,1][IDSP])))
    Ixy=-np.sum(particles['Masses'][IDSP]*np.multiply(Dist[:,0][IDSP],Dist[:,1][IDSP]))
    Ixz=-np.sum(particles['Masses'][IDSP]*np.multiply(Dist[:,0][IDSP],Dist[:,2][IDSP]))
    Iyz=-np.sum(particles['Masses'][IDSP]*np.multiply(Dist[:,1][IDSP],Dist[:,2][IDSP]))
    I=np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
    eig_vals=np.sqrt(np.linalg.eigvals(I))
    return np.sort(eig_vals)[::-1]


def add_InertiaTensor(basePath, Subhalos_cat, selected_subhalos,  N,Header, FieldLIM, ValLIM):
    eigen=np.empty((np.size(np.where(selected_subhalos)[0]),3))
    CoA=np.empty(np.size(np.where(selected_subhalos)[0]))
    NumStarPart=np.empty(np.size(np.where(selected_subhalos)[0]))
    
    print("Computing the inertia tensor ")
    kl=0
    for selH in np.where(selected_subhalos)[0]:
        #print(f'Loading data for Subhalo {selH}')
        IT_Fields=['Masses', 'Coordinates']
        stars_particles=il.snapshot.loadSubhalo(basePath,N,selH,'stars',fields=IT_Fields)
        NumStarPart[kl]=np.size(stars_particles['Masses'])
        aux_eig=comp_Inertia_tensor(stars_particles,Subhalos_cat['SubhaloPos'][selH,:], Header, ValLIM*Subhalos_cat[FieldLIM][:,4][selH]) 
        #print(aux_eig)
        eigen[kl,:]=aux_eig
        CoA[kl]=aux_eig[-1]/aux_eig[0]
        kl+=1
    return CoA, eigen, NumStarPart


def add_ReducedInertiaTensor(basePath, Subhalos_cat, selected_subhalos, N, Header, FieldLIM, ValLIM):
    eigen=np.empty((np.size(np.where(selected_subhalos)[0]),3))
    qs=np.empty((np.size(np.where(selected_subhalos)[0]),2))
    NumStarPart=np.empty(np.size(np.where(selected_subhalos)[0]))
    
    print("Computing the inertia tensor ")
    kl=0
    for selH in np.where(selected_subhalos)[0]:
        #print(f'Loading data for Subhalo {selH}')
        IT_Fields=['Masses', 'Coordinates']
        stars_particles=il.snapshot.loadSubhalo(basePath,N,selH,'stars',fields=IT_Fields)
        NumStarPart[kl]=np.size(stars_particles['Masses'])
        #print(np.shape(Subhalos_cat['SubhaloCM'][selH,:]))
        aux_eig, aux_qs=comp_Reduced_Inertia_tensor(stars_particles,Subhalos_cat['SubhaloPos'][selH,:], Header, ValLIM*Subhalos_cat[FieldLIM][:,4][selH]) 
        #print(aux_eig)
        eigen[kl,:]=aux_eig
        qs[kl,:]=aux_qs
        kl+=1
    return qs, eigen, NumStarPart


def comp_Reduced_Inertia_tensor(particles, SubhaloPos, Header, DistanceLIM):
    Dist= Distance_3D(particles['Coordinates'],SubhaloPos, Header['BoxSize'])
    #print(DistanceLIM)
    DISTMOD=np.linalg.norm(Dist,axis=1)
    IDSP=np.logical_and(DISTMOD<=float(DistanceLIM), DISTMOD>0)
    Wrn=1/np.square(np.linalg.norm(Dist,axis=1))
    MTOT=np.sum(particles['Masses'][IDSP])
    Ixx=np.sum(particles['Masses'][IDSP]*np.square(Dist[:,0][IDSP])*Wrn[IDSP])/MTOT
    Iyy=np.sum(particles['Masses'][IDSP]*np.square(Dist[:,1][IDSP])*Wrn[IDSP])/MTOT
    Izz=np.sum(particles['Masses'][IDSP]*np.square(Dist[:,2][IDSP])*Wrn[IDSP])/MTOT
    Ixy=np.sum(particles['Masses'][IDSP]*np.multiply(Dist[:,0][IDSP],Dist[:,1][IDSP])*Wrn[IDSP])/MTOT
    Ixz=np.sum(particles['Masses'][IDSP]*np.multiply(Dist[:,0][IDSP],Dist[:,2][IDSP])*Wrn[IDSP])/MTOT
    Iyz=np.sum(particles['Masses'][IDSP]*np.multiply(Dist[:,1][IDSP],Dist[:,2][IDSP])*Wrn[IDSP])/MTOT
    I=np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
    try: a, b, c=np.sqrt(np.sort(np.linalg.eigvals(I))[::-1])
    except:
        print("Problem detected with the inertia tensor: ", Wrn)
        a,b,c=np.array([np.nan, np.nan, np.nan])
    q=b/a
    s=c/a
    return np.array([a,b,c]), np.array([q,s])
#def comp_Reduced_Inertia_tensor(stars_particles,Subhalos_cat['SubhaloCM'][selH,:], Header):
#def comp_Reduced_Inertia_tensor(stars_particles,Subhalos_cat, Header):
#    l=0
#    Dist= Distance_3D(particles['Coordinates'],SubhaloCM, Header['BoxSize'])
#    if l=0:
#        IDSP=np.linalg.norm(Dist,axis=1)<=np.PINF
#    else:
#        "compute the IDS"
#        pass
#    while np.logical_or(l=0, q_error> 0.01, s_error> 0.01):
#        qa=q
#        sa=s
#        Wrn=1/np.square(np.linalg.norm(Dist,axis=1))
#        MTOT=np.sum(particles['Masses'][IDSP])
#        Ixx=np.sum(particles['Masses'][IDSP]*np.square(Dist[:,0][IDSP])*Wrn[IDSP])/MTOT
#        Iyy=np.sum(particles['Masses'][IDSP]*np.square(Dist[:,1][IDSP])*Wrn[IDSP])/MTOT
#        Izz=np.sum(particles['Masses'][IDSP]*np.square(Dist[:,2][IDSP])*Wrn[IDSP])/MTOT
#        Ixy=np.sum(particles['Masses'][IDSP]*np.multiply(Dist[:,0][IDSP],Dist[:,1][IDSP])*Wrn[IDSP])/MTOT
#        Ixz=np.sum(particles['Masses'][IDSP]*np.multiply(Dist[:,0][IDSP],Dist[:,2][IDSP])*Wrn[IDSP])/MTOT
#        Iyz=np.sum(particles['Masses'][IDSP]*np.multiply(Dist[:,1][IDSP],Dist[:,2][IDSP])*Wrn[IDSP])/MTOT
#        I=np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
#        abc=np.sqrt(np.sort(np.linalg.eigvals(I))[::-1])
#        a=abc[0]
#        b=abc[1]
#        c=abc[2]
#        q=b/a
#        s=c/a
#        qerror=abs(q-qa)/q
#        serror=abs(s-sa)/s
#    
#    return np.array([a,b,c])
    
########################################
##### Compute the Gas Temperature Fractions

def temperature(u, mu, UEnergyOUMass=1.0E10, kB=1.38E-16, gamma=5/3):
    return (gamma-1)*u/kB*UEnergyOUMass*mu


def mu(xe, xH=0.76, mp=1.6726E-24):
    return 4/(1+3*xH+4*xH*xe)*mp


def compute_mass_gas_temp_limits(masses, temp, temp_min=0, temp_max=np.PINF):
    ids_fits_temp=np.logical_and(temp>=temp_min, temp<temp_max)
    return np.sum(masses[ids_fits_temp])


def compute_mass_gas_temp(masses, temperature):
    ### 0 Hot gas T>=1E6
    HOTGAS=compute_mass_gas_temp_limits(masses, temperature, temp_min=1E6)
    ### 1 Warm gas 1E5<=T<1E6
    WARMGAS=compute_mass_gas_temp_limits(masses, temperature, temp_min=1E5, temp_max=1E6)
    ### 2 Cool gas 1E5<=T<1E4
    COOLGAS=compute_mass_gas_temp_limits(masses, temperature, temp_min=1E4, temp_max=1E5)
    ### 3 Cold gas T<=1E4
    COLDGAS=compute_mass_gas_temp_limits(masses, temperature, temp_max=1E4)
    return np.array([HOTGAS, WARMGAS, COOLGAS, COLDGAS, HOTGAS+WARMGAS+COOLGAS+COLDGAS])

#######Based on Torrey2017
def compute_mass_gas_temp_limitsT2017(masses, temp, density, temp_min=0, temp_max=np.PINF, density_min=np.NINF, density_max=np.PINF):
    ids_fits_temp=np.logical_and(np.logical_and(np.logical_and(temp>=temp_min, temp<temp_max), density>=density_min), density<density_max)
    return np.sum(masses[ids_fits_temp])


def compute_mass_gas_tempT2017(masses, temperature, density):
    """This function takes as input the masses of the particles, tenperature in K, and log of density in nH/cm**-3"""
    ### 0 Hot gas T>=1E7
    HOTGAS=compute_mass_gas_temp_limitsT2017(masses, temperature, density,  temp_min=1E7)
    ### 1 WHIM gas 1E5<=T<1E7
    WARMGAS=compute_mass_gas_temp_limitsT2017(masses, temperature, density, temp_min=1E5, temp_max=1E7)
    ### 2 Cold and diffuse gas T<=1E5 and log(density) < -3.6
    CDIF=compute_mass_gas_temp_limitsT2017(masses, temperature, density,  temp_max=1E5, density_max=2.5E-4)
    ### 3 Condensed gas T<=1E5 and log(density) > -3.6
    COND=compute_mass_gas_temp_limitsT2017(masses, temperature, density, temp_max=1E5, density_min=2.5E-4)
    return np.array([HOTGAS, WARMGAS, CDIF, COND, HOTGAS+WARMGAS+CDIF+COND])
                

def add_SubhaloMassTemp(Subhalos_cat, selected_subhalos, snapNum, basePath):
    num_selec_subhalos=np.where(selected_subhalos)[0]
    aux_subhalos={'GasContentTemp': np.empty((np.shape(num_selec_subhalos)[0],5))}
    cont=0
    for i in num_selec_subhalos:
        part=il.snapshot.loadSubhalo(basePath, snapNum, id=i, partType='gas',fields= ['Masses', 'InternalEnergy', 'ElectronAbundance'])
        if 'Masses' in part.keys():
            mu_v=mu(part['ElectronAbundance'])
            temp=temperature(part['InternalEnergy'], mu_v)
            aux_subhalos['GasContentTemp'][cont,:]=compute_mass_gas_temp(part['Masses'], temp)
        else:
            aux_subhalos['GasContentTemp'][cont,:]=np.array([0.0,0.0,0.0,0.0,0.0])
        cont=cont+1
    print(aux_subhalos['GasContentTemp'])
    return aux_subhalos['GasContentTemp']
###############################################
############ Compute the RAM PRESSURE AND VELOCITIES

def add_RAM_P(Halos_cat, Subhalos_cat, selected_halos, selected_subhalos, snapNum, Header, basePath, aux_Subhalos_cat):
    aux_Subhalos_cat['SubhaloVMedium']=np.empty(np.shape(aux_Subhalos_cat['SubhaloVel']))
    aux_Subhalos_cat['SubhaloRhoMedium']=np.empty(np.shape(aux_Subhalos_cat['SubhaloGrNr']))
    aux_Subhalos_cat['SubhaloRAMPress']=np.empty(np.shape(aux_Subhalos_cat['SubhaloGrNr']))
    cont=0
    count_halos=0
    for i in np.where(selected_halos)[0]:
        print(f"Loading the particles of the halo number {i} ")
        Halo_gas=il.snapshot.loadHalo(basePath, snapNum, id=i, partType='gas', fields=['Masses', 'ParticleIDs', 'Coordinates', 'Density', 'Velocities'])
        print(f"Gas particles of the halo number {i} have been loaded")
        IDs_Sat_particle_list=[]
        AUX_HLI=np.where(np.logical_and(Subhalos_cat['SubhaloMassType'][:,0]>0, Subhalos_cat['SubhaloGrNr']==Halos_cat['GroupNumber'][count_halos]))[0]
        AUX_HLI=np.delete(AUX_HLI, np.in1d(AUX_HLI, np.where(np.in1d(Subhalos_cat['SubhaloNumber'],Halos_cat['GroupFirstSub']))[0]))
        HBLACK_LIST=[]
        for k in AUX_HLI:
            PLHalo=il.snapshot.loadSubhalo(basePath, snapNum, id=k, partType='gas', fields=['ParticleIDs'])
            if np.logical_not(type(PLHalo) is dict):
                IDs_Sat_particle_list=np.append(IDs_Sat_particle_list, PLHalo).astype(int)
            else:
                HBLACK_LIST.append(k)
        #print(np.shape(IDs_Sat_particle_list))
        #print("We finished to load the IDS of particles beolonging to all satellites")
        IDs_not_in_Sat=np.logical_not(np.in1d(Halo_gas['ParticleIDs'], IDs_Sat_particle_list))
        #print("We finished to compare the halo and satellite particles")

        #print( np.size(np.where(np.logical_and(selected_subhalos, Subhalos_cat['SubhaloGrNr']==Halos_cat['GroupNumber'][count_halos]))[0]), "Galaxies to be compared")
        #print("There are:", np.size(Halo_gas['ParticleIDs']), "particles in the Halo")
        Halo_gas_WOUT_sat=dictionary_select_items(Halo_gas, IDs_not_in_Sat)
        aux_Subhalos_cat['SubhaloVMedium'][cont,:]=np.array([np.nan, np.nan, np.nan])
        aux_Subhalos_cat['SubhaloRhoMedium'][cont]=np.nan
        aux_Subhalos_cat['SubhaloRAMPress'][cont]=np.nan
        cont=cont+1
        IDS_halos=np.where(np.logical_and(selected_subhalos, Subhalos_cat['SubhaloGrNr']==Halos_cat['GroupNumber'][count_halos]))[0]
        IDS_satelites=np.delete(IDS_halos, np.in1d(IDS_halos, np.where(np.in1d(Subhalos_cat['SubhaloNumber'],Halos_cat['GroupFirstSub']))[0]))
        for k in IDS_satelites:
            if k not in HBLACK_LIST:
                Halo_gas_WOUT_sat['TempDist']= Distance_3D(Halo_gas_WOUT_sat['Coordinates'], Subhalos_cat['SubhaloPos'][k,:],Header['BoxSize'] )
                IDs_selected_particles=Halo_gas_WOUT_sat['ParticleIDs'][np.linalg.norm(Halo_gas_WOUT_sat['TempDist'],axis=1)<20*Subhalos_cat['SubhaloHalfmassRadType'][k,4]]
                if IDs_selected_particles.size!=0:
                    aux_Subhalos_cat['SubhaloVMedium'][cont,:]=np.array([np.median(Halo_gas_WOUT_sat['Velocities'][np.in1d(Halo_gas_WOUT_sat['ParticleIDs'], IDs_selected_particles),0]),np.median(Halo_gas_WOUT_sat['Velocities'][np.in1d(Halo_gas_WOUT_sat['ParticleIDs'], IDs_selected_particles),1]), np.median(Halo_gas_WOUT_sat['Velocities'][np.in1d(Halo_gas_WOUT_sat['ParticleIDs'], IDs_selected_particles),2])])
                    #print(Halo_gas_WOUT_sat['Density'][np.in1d(Halo_gas_WOUT_sat['ParticleIDs'], IDs_selected_particles)])
                    aux_Subhalos_cat['SubhaloRhoMedium'][cont]=mode(Halo_gas_WOUT_sat['Density'][np.in1d(Halo_gas_WOUT_sat['ParticleIDs'], IDs_selected_particles)])
                else:
                    aux_Subhalos_cat['SubhaloVMedium'][cont,:]=np.array([0,0,0])
                    aux_Subhalos_cat['SubhaloRhoMedium'][cont]=0
                if Subhalos_cat['SubhaloMassType'][k,0]>0:
                    Subhalo_Gas_Particles=il.snapshot.loadSubhalo(basePath, snapNum, id=k, partType='gas', fields=['Coordinates', 'Velocities'])
                    Vrel=mode(np.linalg.norm(Subhalo_Gas_Particles['Velocities']-aux_Subhalos_cat['SubhaloVMedium'][cont,:], axis=1))
                else: 
                    Vrel=np.nan
            else:
                aux_Subhalos_cat['SubhaloVMedium'][cont,:]=np.array([np.nan, np.nan, np.nan])
                aux_Subhalos_cat['SubhaloRhoMedium'][cont]=np.nan

            aux_Subhalos_cat['SubhaloRAMPress'][cont]=aux_Subhalos_cat['SubhaloRhoMedium'][cont]*Vrel**2
            cont=cont+1
        #print("We finished to load the IDS of particles beolonging to all satellites")
        count_halos=count_halos+1
    print("Finished the process")
    aux_Subhalos_cat['SubhaloVMediumPC']=aux_Subhalos_cat['SubhaloVMedium'].copy()*np.sqrt(Header['Time'])
    aux_Subhalos_cat['SubhaloRhoMediumPC']=aux_Subhalos_cat['SubhaloRhoMedium'].copy()*6.7677E-12*(Header['HubbleParam'])**2/(Header['Time'])**3
    aux_Subhalos_cat['SubhaloRAMPressPC']=aux_Subhalos_cat['SubhaloRAMPress'].copy()*6.7677E-12*(Header['HubbleParam'])**2/(Header['Time'])**2
    return aux_Subhalos_cat