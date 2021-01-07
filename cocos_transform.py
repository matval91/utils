import numpy as np

def cocos_transform(eqd, COCOSin, COCOSout, sigma_ip_out=0, sigma_b0_out=0):
    """ cocos transformation
    This function converts the magnetic input from their starting cocos to cocos 3 (needed by ascot5)

    https://www.sciencedirect.com/science/article/pii/S0010465512002962

    Table 4 is not updated. EFIT has usually cocos=7, but if people try to keep q>0, then cocos changes.

    Parameters:
        COCOSin (int): input cocos
        COCOSout (int): output cocos
        sigma_ip_out (int): output sign of ip requested. If 0, it will keep the ones in the eqdsk
        sigma_b0_out (int): output sign of b0 requested. If 0, it will keep the ones in the eqdsk
        
    Attributes:
        None
    """
    print("COCOS tranformation from "+str(COCOSin)+" to "+str(COCOSout))
    eqdout = eqd

    cocosin = fill_cocosdict(COCOSin)
    cocosin['sigma_ip'] = np.sign(eqd.Ip)
    cocosin['sigma_b0'] = np.sign(eqd.B0EXP)

    #These cocos are for COCOS 3
    cocosout = fill_cocosdict(COCOSout)
    #Checking the signs of the current and field desired as output
    cocosout['sigma_ip'] = np.sign(eqd.Ip)    if sigma_ip_out == 0 else sigma_ip_out
    cocosout['sigma_b0'] = np.sign(eqd.B0EXP) if sigma_b0_out == 0 else sigma_b0_out

    # Define effective variables: sigma_Ip_eff, sigma_B0_eff, sigma_Bp_eff, exp_Bp_eff as in Appendix C
    sigma_Bp_eff = cocosin['sigma_Bp'] * cocosout['sigma_Bp']
    exp_Bp_eff   = cocosout['exp_Bp']  - cocosin['exp_Bp']
    sigma_rhothetaphi_eff = cocosin['sigma_rhothetaphi'] * cocosout['sigma_rhothetaphi']
    if 'sigma_ip' in cocosout.keys() and 'sigma_b0' in cocosout.keys():
        print(f"Requested sign(Ip)= {cocosout['sigma_ip']}, sign(b0)= {cocosout['sigma_b0']}")
        sigma_Ip_eff = cocosin['sigma_ip']*cocosout['sigma_ip']
        sigma_B0_eff = cocosin['sigma_b0']*cocosout['sigma_b0']
    else:
        print('No sign of Ip nor B0 requested')
        sigma_Ip_eff = cocosin['sigma_RphiZ']*cocosout['sigma_RphiZ']
        sigma_B0_eff = cocosin['sigma_RphiZ']*cocosout['sigma_RphiZ']

    # Define input
    F_in       = eqd.T
    FFprime_in = eqd.TTprime
    pprime_in  = eqd.pprime

    psirz_in   = eqd.psi
    psiaxis_in = eqd.psiaxis
    psiedge_in = eqd.psiedge
    psigrid_in = eqd.psi_grid

    q_in  = eqd.q
    b0_in = eqd.B0EXP
    ip_in = eqd.Ip
    
    # Transform
    F         = F_in       * sigma_B0_eff
    FFprime   = FFprime_in * sigma_Ip_eff * sigma_Bp_eff / (2*np.pi)**exp_Bp_eff
    pprime    = pprime_in  * sigma_Ip_eff * sigma_Bp_eff / (2*np.pi)**exp_Bp_eff
    
    _fact_psi = sigma_Ip_eff * sigma_Bp_eff * (2*np.pi)**exp_Bp_eff
    psirz     = psirz_in   * _fact_psi 
    psi_grid  = psigrid_in * _fact_psi
    psiaxis   = psiaxis_in * _fact_psi
    psiedge   = psiedge_in * _fact_psi
    
    q  = q_in  * sigma_Ip_eff * sigma_B0_eff * sigma_rhothetaphi_eff
    b0 = b0_in * sigma_B0_eff
    ip = ip_in * sigma_Ip_eff

    # Define output
    eqdout.T       = F
    eqdout.TTprime = FFprime
    eqdout.pprime  = pprime
    
    eqdout.psi      = psirz
    eqdout.psi_grid = psi_grid
    eqdout.psiaxis  = psiaxis
    eqdout.psiedge  = psiedge
    
    eqdout.q     = q
    eqdout.B0EXP = b0
    eqdout.Ip    = ip
    return eqdout

def fill_cocosdict(COCOS):
    """
    Function to fill the dictionary with the COCOS variables

    Parameters:
        COCOS (int): input cocos
    Args:
        cocosdict (dict): dictionary with cocos variables
    """
    cocos_keys = ['sigma_Bp', 'sigma_RphiZ', 'sigma_rhothetaphi',\
          'sign_q_pos', 'sign_pprime_pos', 'exp_Bp']
    cocosdict = dict.fromkeys(cocos_keys)

    cocosdict['exp_Bp'] = 1 if COCOS > 10 else 0 # if COCOS>=10, this should be 1

    if COCOS==1 or COCOS==11:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==2 or COCOS==12:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==3 or COCOS==13:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = +1
    elif COCOS==4 or COCOS==14:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = +1
    elif COCOS==5 or COCOS==15:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==6 or COCOS==16:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==7 or COCOS==17:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = +1
    elif COCOS==8 or COCOS==18:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = +1
    else:
        raise ValueError(f"COCOS {COCOS} does not exists \n")

    return cocosdict
