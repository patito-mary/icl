import numpy as np

def sim_params(SIMTYPE):
    if SIMTYPE=="TNG100-1":
        L_Box=110.7 ###Mpc
        N_DM=1820^3 
        M_DM=7.5E6 ##Solar masses
        M_BP=1.4E6 ##Solar masses
        N_snap=100
        N_Subfind=4371211 #At (z=0)
        h=0.6774
    if SIMTYPE=="TNG300-1":
        L_Box=302.6 ###Mpc
        N_DM=2500^3 
        M_DM=5.9E7
        M_BP=1.1E7
        N_snap=100
        N_Subfind=14485709 ##At(z=0)
        h=0.6774
    return L_Box,N_DM,M_DM,M_BP,N_snap,N_Subfind,h
#################################################################################################################################################################################File to input the parameters to Create the Catalogue to be used in the Analysis. 
#######################################################################################
#Path where the simulations files are storaged. Use TNG100-1 for the 100Mpc simuations.
basePath = '/home/tnguser/sims.TNG/TNG100-1/output/'
#List of snapshots to be used for the creation of the Catalogue
N_list=[99]
#Parameters to select the clusters to be used....
[_,_,M_DM,M_BP,_,_,h]=sim_params("TNG100-1")
MP_DM=1000
MP_BP=715
Halos_limits={"Group_M_Crit200":np.array([1.1E13, np.PINF])*h/1E10}
Subhalos_limits={"SubhaloMassType/1":np.array([MP_DM*M_DM, np.PINF])*h/1E10,
    "SubhaloMassType/4":np.array([MP_BP*M_BP, np.PINF])*h/1E10}
#List of fields to storage in the new Halos and Subhalos Catalogue (Each of those fields must exist in the original Catalogue.
Halos_Fields=['GroupFirstSub','GroupNsubs','Group_M_Crit200','GroupCM','Group_R_Crit200','GroupPos', 'GroupGasMetalFractions','GroupGasMetallicity','GroupMass', 'GroupMassType', 'GroupNsubs','GroupSFR','GroupStarMetallicity','GroupWindMass','Group_M_Crit200','Group_R_Crit200']
#Halos_Fields=['GroupFirstSub','GroupNsubs','Group_M_Crit200','GroupCM','Group_R_Crit200','GroupPos']
Subhalos_Fields=['SubhaloFlag', 'SubhaloCM','SubhaloPos','SubhaloGrNr','SubhaloMassType','SubhaloSFR','SubhaloMass']
#Halos and Subhalos catalogue files names to save the outputs. 
OHalos="Haloshift.pkl"
OSubhalos="SubHaloshift.pkl"