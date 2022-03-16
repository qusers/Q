import matplotlib.pyplot as plt
import math

def linplot(storedir='.', d='dG', x=None, y=None, error=None):
    unit = {
        'ddG' : '\u0394\u0394G',
        'dG'  : '\u0394G'    
    }

    # create the dimensions for the figure
    dim = math.ceil(abs(max(x+y, key=abs))+abs(max(error, key=abs)))

    ticks = range(-dim, dim + 1)

    fig, ax = plt.subplots(figsize=(dim+2,dim+2))
    ax.set_xlim(-dim, dim)
    ax.set_ylim(-dim, dim)    

    #plot reference lines
    ax.plot((-dim, dim), (-dim, dim), color='black')
    ax.plot((-dim, dim-1), (-dim+1, dim), color='black', linestyle='--', alpha=0.7)
    ax.plot((-dim+1, dim), (-dim, dim-1), color='black', linestyle='--', alpha=0.7)
    ax.plot((-dim, dim-2), (-dim+2, dim), color='black', linestyle='-.', alpha=0.7)
    ax.plot((-dim+2, dim), (-dim, dim-2), color='black', linestyle='-.', alpha=0.7)

    # ticks, labels and grid
    ax.set_xlabel('{} exp [$kcal mol^{}$]'.format(unit[d], -1), fontsize=18)
    ax.set_ylabel('{} pred [$kcal mol^{}$]'.format(unit[d], -1), fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xticks (ticks)
    ax.set_yticks (ticks)
    ax.errorbar(x, y, error, fmt='o')
    ax.grid(alpha=0.3)

    # this is to plot the fit line
    # will move to metrics generation part
    #targs =  [i for target in targs for i in target]
    #preds =  [i for target in preds for i in target]
    #sems =  [i for target in sems for i in target]   
    #targs = [j for i, j in enumerate(targs) if i not in ban_idx]
    #preds = [j for i, j in enumerate(preds) if i not in ban_idx]
    #sems = [j for i, j in enumerate(sems) if i not in ban_idx]
    #m = analysis(targs, preds, sems)
    #reg = [m['slope'] * x + m['intercept'] for x in np.arange(-18, 9, 1)]
    #ax[0].plot(np.arange(-18, 9, 1), reg, color='orange', linewidth=3)
    #ax[0].legend(loc=4)
    plt.savefig('{}/plot_{}.png'.format(storedir, d), dpi=300, bbox_inches='tight')