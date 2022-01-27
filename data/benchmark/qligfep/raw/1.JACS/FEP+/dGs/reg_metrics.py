import numpy as np
import glob
from scipy import stats
import matplotlib.pyplot as plt 
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import matthews_corrcoef

ddG = []
categorised = {}

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def pearsonr_ci(x,y,alpha=0.05):
    ''' calculate Pearson correlation along with the confidence interval using scipy and numpy
    Parameters
    ----------
    x, y : iterable object such as a list or np.array
      Input for correlation calculation
    alpha : float
      Significance level. 0.05 by default
    Returns
    -------
    r : float
      Pearson's correlation coefficient
    pval : float
      The corresponding p value
    lo, hi : float
      The lower and upper bound of confidence intervals
    '''

    r, p = stats.pearsonr(x,y)
    r_z = np.arctanh(r)
    se = 1/np.sqrt(x.size-3)
    z = stats.norm.ppf(1-alpha/2)
    lo_z, hi_z = r_z-z*se, r_z+z*se
    lo, hi = np.tanh((lo_z, hi_z))
    return r, p, lo, hi

def read_csv(input_file):
    exp = []
    calc = []
    with open(input_file, 'r', encoding='utf-8-sig') as infile:
        for line in infile:
            line = line.split(',')
            if line[0] == 'Exp. dG':
              continue
            exp.append(float(line[0])) 
            calc.append(float(line[1]))
    infile.close
    return exp, calc
            
    
def analyse(experiment, calculated):
    x = np.array([float(i) for i in experiment])
    y = np.array([float(i) for i in calculated])
    x2 = np.sign(x)
    y2 = np.sign(y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    MAE = mean_absolute_error(x, y)
    RMSE = rmse(x,y)
    #MCC = matthews_corrcoef(x2,y2)
    tau, p_value = stats.kendalltau(x, y)
    rho, pval = stats.spearmanr(x,y)
    r, p, lo, hi = pearsonr_ci(x,y)
    lo = lo**2
    hi = hi**2
    print("r2", round(r**2, 3))
    print("r2 CI", round(lo, 3), "-", round(hi, 3))
    print("MAE", round(MAE, 3))
    print("RMSE", round(RMSE, 3))
    #print("MCC", round(MCC, 3))
    print("slope", round(slope, 3))
    print("intercept", round(intercept, 3))
    print("p_value", round(p_value, 3))
    print("std_error", round(std_err, 3))
    print("kendall_tau", round(tau, 3))
    print("spearman_rho", round(rho, 3))
    print("n", len(x))

csvfiles = glob.glob('*.csv')
for csv in csvfiles:
    name = csv[:-4]
    print('Regression metrics for dG for {} target'.format(name))
    exp, calc = read_csv(csv)
    analyse(exp, calc)
  
#print(csvfiles)

