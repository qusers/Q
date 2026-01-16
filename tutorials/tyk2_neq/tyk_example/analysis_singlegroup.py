import os
import time
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

#this is a script that carries out BAR analysis for a given replicate number
#it iterates over the folders in the water and protein group, and reads in the work values produced by the NEQ script
#then the csv files can be examined individually, or collected using other scripts.


#calculation time

from scipy.optimize import brentq
from scipy.stats import gaussian_kde

rep_trial = int(sys.argv[1])

#parameters for sampling

def get_works_from_logs(folder):
    Wf = []
    Wr = []
    init_fwd = 0
    init_rev = 0
    forward_logs = glob.glob(os.path.join(folder, "neq_1_*.log"))
    reverse_logs = glob.glob(os.path.join(folder, "neq_0_*.log"))
    file_default_array = []
    # Process forward logs
    for fname in forward_logs:
        second_flag = False
        with open(fname, "r") as f:
            lines = f.readlines()
            if not lines:
                continue
            if lines[-1].strip() == "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" or lines[-1].strip() != "================================================================================":
                continue
            for line in lines:
                    if "work" in line:
                        try:
                            if second_flag:
                                init_fwd = (float(line.split()[6]))
                                break
                            else:
                                second_flag = True
                        except Exception as e:
                            print(f"Failed to parse work in {fname}, {e}")
            
            for line in lines[-5:]:
                if "work" in line:
                    try:
                        Wf.append(float(line.split()[6]))
                    except Exception:
                        print(f"Failed to parse work in {fname}")

    # Process reverse logs
    for fname in reverse_logs:
            second_flag = False
            with open(fname, "r") as f:
                lines = f.readlines()
                if not lines:
                    continue
                if lines[-1].strip() == "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" or lines[-1].strip() != "================================================================================":
                    continue
                for line in lines:
                    if "work" in line:
                        try:
                            if second_flag:
                                init_rev = (float(line.split()[6]))
                                break
                            else:
                                second_flag = True
                        except Exception:
                            print(f"Failed to parse work in {fname}")
                for line in lines[-5:]:
                    if "work" in line:
                        try:
                            Wr.append(float(line.split()[6]))
                            #print(float(line.split()[6]))
                        except Exception:
                            print(f"Failed to parse work in {fname}")

    print(f'work values in {folder}: Forward={len(forward_logs)}, Reverse={len(reverse_logs)}')
    print(len(Wf),len(Wr))
    #print(Wf)
    return Wf, Wr

def do_calc(Wf, Wr,image_prefix=''):

    
  
  # conver Wf and Wr to np.array
    Wf, Wr = np.array(Wf), np.array(Wr)
    Nf, Nr = len(Wf), len(Wr)
  
    if np.min(Wf) > np.max(-Wr) or np.min(-Wr) > np.max(Wf):
        print("Warning: minimal overlap — BAR unreliable.")
  #set units here
    T = 300.0          # temperature in K
    units = "kBT"     # choose: "kBT", "kcal", or "kJ". Note that these are units in which WORK is provided, not the units of the output
    
    # Boltzmann constants
    kB_kcal = 0.0019872041   # kcal/(mol*K)
    kB_kJ   = 0.0083144621   # kJ/(mol*K)


    # -------------------------
    # UNIT HANDLING
    # -------------------------
    if units == "kBT":
        beta = 1.0
        # works already in kBT
    elif units == "kcal":
        beta = 1.0 / (kB_kcal * T)
        # assume Wf/Wr provided in kcal/mol
    elif units == "kJ":
        beta = 1.0 / (kB_kJ * T)
        # assume Wf/Wr provided in kJ/mol
    else:
        raise ValueError("units must be 'kBT', 'kcal', or 'kJ'")



    # BAR function with correct sign convention
    def bar_func(dF, Wf, Wr, beta):
        a = 1.0 / (1.0 + np.exp(beta * (Wf - dF)))
        b = 1.0 / (1.0 + np.exp(beta * (Wr + dF)))
        return np.sum(a) - np.sum(b)

    # Solve for ΔF
    xs = np.linspace(-100, 100, 10001)
    vals = np.array([bar_func(x, Wf, Wr, beta) for x in xs])
    idxs = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    dF_root = brentq(bar_func, xs[idxs[0]], xs[idxs[0]+1], args=(Wf, Wr, beta))
    
    
    # Bootstrap uncertainty
    rng = np.random.default_rng(12345)
    dF_boot = []
    for _ in range(1000):
        Wf_bs = Wf[rng.integers(0, Nf, Nf)]
        Wr_bs = Wr[rng.integers(0, Nr, Nr)]
        try:
            root_bs = brentq(bar_func, dF_root-10, dF_root+10, args=(Wf_bs, Wr_bs, beta))
            dF_boot.append(root_bs)
        except:
            pass

    dF_boot = np.array(dF_boot)
    mean_boot = np.mean(dF_boot)
    std_boot = np.std(dF_boot, ddof=1)
    ci_low, ci_high = np.percentile(dF_boot, [2.5, 97.5])

    # -------------------------
    # Unit conversions if input was kBT
    # -------------------------
    if units == "kBT":
        dF_kBT = dF_root
        dF_kcal = dF_root * (kB_kcal * T)
        dF_kJ = dF_root * (kB_kJ * T)

        mean_kcal = mean_boot * (kB_kcal * T)
        mean_kJ   = mean_boot * (kB_kJ * T)

        ci_low_kcal = ci_low * (kB_kcal * T)
        ci_high_kcal = ci_high * (kB_kcal * T)

        ci_low_kJ = ci_low * (kB_kJ * T)
        ci_high_kJ = ci_high * (kB_kJ * T)

    # -------------------------

    # KDEs for overlap visualization
    WnegR = -Wr
    kde_f = gaussian_kde(Wf)
    kde_rneg = gaussian_kde(WnegR)
    xgrid = np.linspace(min(Wf.min(), WnegR.min())-1, max(Wf.max(), WnegR.max())+1, 2000)
    overlap = np.trapz(np.minimum(kde_f(xgrid), kde_rneg(xgrid)), xgrid)
    print(f'the overlap is approx. {overlap}')

    # Plot 1: Histograms + KDEs + ΔF line
    plt.figure(figsize=(8,6))
    plt.hist(Wf, bins=100, alpha=0.5, density=True, label='Forward Wf')
    plt.hist(WnegR, bins=100, alpha=0.5, density=True, label='Negated reverse -Wr')
    plt.plot(xgrid, kde_f(xgrid), label='KDE Wf')
    plt.plot(xgrid, kde_rneg(xgrid), label='KDE -Wr')
    plt.axvline(dF_root, linestyle='--', color='red', label=f'BAR ΔF = {dF_root:.3f}')
    plt.xlabel('Work')
    plt.ylabel('Density')
    #plt.xlim(-100,100)
    plt.legend()
    plt.title('Forward & Negated-Reverse Work Distributions')
    plt.tight_layout()
    plt.savefig(image_prefix+'bar_hist.png', dpi=200)

    # Plot 2: Bootstrap ΔF distribution
    plt.figure(figsize=(8,5))
    plt.hist(dF_boot, bins=40, density=True, alpha=0.7)
    plt.axvline(mean_boot, linestyle='--', color='red', label=f'mean = {mean_boot:.3f}')
    plt.axvline(ci_low, linestyle=':', color='black', label=f'95% CI = [{ci_low:.3f}, {ci_high:.3f}]')
    plt.axvline(ci_high, linestyle=':')
    plt.xlabel('ΔF (bootstrap samples)')
    plt.ylabel('Density')
    plt.legend()
    plt.title('Bootstrap Distribution of BAR ΔF')
    plt.tight_layout()
    plt.savefig(image_prefix+'bar_bootstrap.png', dpi=200)
    plt.close()

    # -------------------------
    # Results
    # -------------------------
    print(f"BAR ΔF = {dF_root:.6f} {units}")
    print(f"Bootstrap mean = {mean_boot:.6f}, std = {std_boot:.6f} {units}")
    print(f"95% CI = [{ci_low:.6f}, {ci_high:.6f}] {units}")

    if units == "kBT":
        print("\nConverted results:")
        print(f"ΔF = {dF_kBT:.3f} kBT")
        print(f"ΔF = {dF_kcal:.3f} kcal/mol")
        print(f"ΔF = {dF_kJ:.3f} kJ/mol")
        print(f"95% CI kcal/mol = [{ci_low_kcal:.3f}, {ci_high_kcal:.3f}]")
        print(f"95% CI kJ/mol   = [{ci_low_kJ:.3f}, {ci_high_kJ:.3f}]")
        print("\n")
        print("\n")
        print("\n")
    return(dF_kBT, dF_kBT-ci_low, overlap)
    
for rep_to_try in range(rep_trial,rep_trial+1):      
    columns = ["folder", "ddF", 'dF_protein', 'dF_water', "ddF_error", 'overlap protein', 'overlap water']
    df = pd.DataFrame(columns=columns)  
    counter=0
    list_of_folders = os.popen(f'ls -d neq_water/*rep{rep_to_try}').read().split()
    
    for folder in list_of_folders:
        try:
            print(f"Processing folder: {folder}")
            folder = folder +'/'
            print(folder)
            folder_prot = folder.replace("neq_water", "neq_protein")
            print(folder_prot)
            folder_water = folder
            Wf_p, Wr_p = get_works_from_logs(folder_prot)
            if len(Wf_p) == 0 or len(Wr_p) == 0:
                print(f"No work values found for group {folder_prot}, skipping...")
                df.loc[counter] = [folder_water, 'NAN','NAN','protein problem','NAN', 'NAN', 'NAN']
                counter+=1
                continue
            dF_p, dF_error_p, overlap_p = do_calc(Wf_p, Wr_p, folder_prot)
            Wf_p = [] 
            Wr_p = []                        
            Wf_w, Wr_w = get_works_from_logs(folder_water)
            if len(Wf_w) == 0 or len(Wr_w) == 0:
                print(f"No work values found for group {folder_water}, skipping...")
                df.loc[counter] = [folder_water, 'NAN','NAN','NAN','water problem', 'NAN', 'NAN']
                counter +=1

                continue
            dF_w, dF_error_w, overlap_w = do_calc(Wf_w, Wr_w, folder_water)
            Wf_w = [] 
            Wr_w = []  
            #calculated ddF
            ddF = dF_p - dF_w
            ddF_error = np.sqrt(dF_error_p**2 + dF_error_w**2)
    
            print(folder_water.split('/')[1])
            print(f"ddF = {ddF} +/- {ddF_error}")
            
            print('\n')
            #time to append the array
            df.loc[counter] = [folder_water, ddF, dF_p, dF_w, ddF_error, overlap_p, overlap_w]
            plt.close('all')
            counter +=1
        except Exception as e:
            print(f"Error processing folder {folder}: {e}")
            df.loc[counter] = [folder_water, 'NAN','NAN','general problem','NAN', 'NAN', 'NAN']
            counter +=1
          
    
    
    
    
    df.to_csv(f"rep{rep_to_try}.csv", index_label="simulation_id")
      
    


# Now do the whole group
