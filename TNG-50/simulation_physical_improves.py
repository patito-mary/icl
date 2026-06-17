def Distance_1D(X, X_POS, BoxSize):
    """This function calculate the physical distance between centers/particles of galaxies
    avoiding the periodic box shifting 

    Args:
        X (float): position of the particle/center of the galaxy [kpc/h]
        X_POS (float): position of the box respecto to the particle/center [kpc/h]
        BoxSize (float): large of the box in kpc/h
    """
    
    D = X-X_POS
    D = np.where(D > BoxSize/2, D-BoxSize, D)
    D = np.where(D < -BoxSize/2, D+BoxSize, D)
    return D

def classify_fossil(group, mags_all, central_mag, 
                    normpos, gals_id, m200_val, r200_val, F_list, NF_list, 
                    r_frac=0.5, gap_threshold=2.0):
    """CLassify halos as fossil and non fossil, by the magnitude gap and R_cut:
    Args:
        group (int): group number
        mags_all (array): array with the magnitudes of all the subhalos
        central_mag (float): magnitude of the central galaxy
        normpos (array): normalized position of the satellites respect to R200
        gals_id (array): array with the IDs of the galaxies
        m200_val (float): M200 value of the halo
        r200_val (float): R200 value of the halo
        F_list (list): list to append the fossil halos
        NF_list (list): list to append the non fossil halos
        r_frac (float, optional): fraction of R200 to consider for the classification. Defaults to 0.5.
        gap_threshold (float, optional): magnitude gap threshold for classification. Defaults to 2.0.
    """
    sats_mask = np.where((normpos <= r_frac) & normpos != 0)[0]
    
    if len(sats_mask) >= 1:
        next_satellite = mags_all[sats_mask].min()
        diff = abs(central_mag - next_satellite)
        
        is_fossil = diff >= gap_threshold
        target = F_list if is_fossil else NF_list
        target[0].append(group)
        target[1].append(diff)
        target[2].append(normpos[sats_mask])
        target[3].append(m200_val)
        target[4].append(r200_val if r200_val > 0 else np.nan)
        target[5].append(gals_id[sats_mask])
        
        label: f'{r_frac:g}'
        print(f'Group {group} classified as {"fossil" if is_fossil else "non-fossil"} with gap {diff:.2f} and R/R200 {normpos[sats_mask]}')
        return True
    return False