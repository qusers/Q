import matplotlib.pyplot as plt

def linplot(storedir='.',include='All', d='ddG',x=None,y=None,error=None):
    targs = []
    preds = []
    sems = []
    k = 0
    markers = ['s', 'o']
    fig, ax = plt.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [15, 1]}, figsize=(8,8), facecolor='white')
    ax[1].axis('off')
    if d == 'ddG':
        type = 'relat'
        ax[0].set_xlim(-8, 8)
        ax[0].set_ylim(-8, 8)
        ax[0].plot((-8, 8), (-8, 8), color='black')
        ax[0].plot((-8, 7), (-7, 8), color='black', linestyle='--', alpha=0.7)
        ax[0].plot((-7, 8), (-8, 7), color='black', linestyle='--', alpha=0.7)
        ax[0].plot((-8, 6), (-6, 8), color='black', linestyle='-.', alpha=0.7)
        ax[0].plot((-6, 8), (-8, 6), color='black', linestyle='-.', alpha=0.7)
    elif d == 'dG':
        type = 'abs'
        ax[0].set_xlim(-16, -2)
        ax[0].set_ylim(-16, -2)
        ax[0].plot((-18, 0), (-18, 0), color='black')
        ax[0].plot((-18, -2), (-17, -1), color='black', linestyle='--', alpha=0.7)
        ax[0].plot((-17, -1), (-18, -2), color='black', linestyle='--', alpha=0.7)
        ax[0].plot((-18, -3), (-16, -1), color='black', linestyle='-.', alpha=0.7)
        ax[0].plot((-16, -1), (-18, -3), color='black', linestyle='-.', alpha=0.7)

    ax[0].set_xlabel('{} exp [$kcal mol^{}$]'.format(d, -1), fontsize=18)
    ax[0].set_ylabel('{} pred [$kcal mol^{}$]'.format(d, -1), fontsize=18)
    ax[0].tick_params(axis='both', which='major', labelsize=16)

    
    targs =  [i for target in targs for i in target]
    preds =  [i for target in preds for i in target]
    sems =  [i for target in sems for i in target]
    
    ax[0].scatter(x, y)
   
    #targs = [j for i, j in enumerate(targs) if i not in ban_idx]
    #preds = [j for i, j in enumerate(preds) if i not in ban_idx]
    #sems = [j for i, j in enumerate(sems) if i not in ban_idx]
    #m = analysis(targs, preds, sems)
    #reg = [m['slope'] * x + m['intercept'] for x in np.arange(-18, 9, 1)]
    #ax[0].plot(np.arange(-18, 9, 1), reg, color='orange', linewidth=3)
    #ax[0].legend(loc=4)
    ax[0].grid(alpha=0.3)
    plt.savefig('{}/plot_{}.png'.format(storedir, d), dpi=300)