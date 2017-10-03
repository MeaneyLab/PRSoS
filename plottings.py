"""
This file holds the plotting components for the PRS_on_Spark
@author: zhuyuecai
@note: the main function is only for testing purpose
@requires: matplotlib, for testing, you need nupmy
"""
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl


import math
import os


"""
@note: normalize the pvalue to make a proper color

"""
def __norm_p_value(p_value):
    SIG = 0.05
    re = 0
    if p_value > SIG :
        re = ((p_value - SIG)/(1-SIG))*0.4+0.6
    else:
        re = (p_value/SIG)*0.4
    return re

def __norm_p_list(p_values):
    return [__norm_p_value(p) for p in p_values]


"""
@author: zhuyuecai
@param pheno_rs: a dictionary with the pheno type name as keys, and a list of ordered R square as values
@param p_values: a list of ordered p values range from 0 to 1
@param width: number of subplot in one row, default 2
"""

def r_square_plots(pheno,rs,p_for_rs, p_values,outputName,width = 2,bar_width = 0.01):
    
    n_plot = len(pheno)
    row = 0
    if n_plot%width == 0 :
        row = n_plot/width
    else:
        row = n_plot/width+1
    p_w = 0.9
    f_h = 0.8
    h_sp = 0.1
    sp = 0.05
    l = 0.05
    b = 0.05
    s_p_h = (f_h - (row-1)*h_sp)/row
    s_p_w = (p_w -(width-1)*sp)/width

    fig = plt.figure(figsize=(20,10))
    axs = []


    n_p = 0
    for i in range(row-1,-1,-1):
        for j in range(width):
            if n_p < n_plot:
                axs.append(fig.add_axes([l+j*(s_p_w+sp), b+i*s_p_h+i*h_sp, s_p_w, s_p_h, ]))
                n_p+=1
    #axs.append(fig.add_axes([l+s_p_w+sp, b+s_p_h+h_sp, s_p_w, s_p_h, ]))
    #axs.append(fig.add_axes([l, b, s_p_w, s_p_h, ]))
   # axs.append(fig.add_axes([l+s_p_w+sp, b, s_p_w, s_p_h, ]))
    axs.append(fig.add_axes([l, 0.95, p_w, 0.05, ]))

    #for i in range(n_plot):
    #    axs.append(fig.add_subplot(math.ceil(float(n_plot)/width+1),width, (i%width)+1))
    #axs.append(fig.add_subplot(n_plot/width+1,1,n_plot%width))
    counter = 0
    anno = '%1.4f'
    highest_correlate = {}
    for i in range(n_plot):
        ms = _get_max_positions_(rs[i])
        
        #axs[counter].bar(p_values,rs[i],(max(p_values)-p_values[0])/len(p_values),color =  cm.get_cmap('cool')(p_for_rs[i]))
        axs[counter].bar(p_values,rs[i],bar_width,color =  cm.get_cmap('seismic')(__norm_p_list(p_for_rs[i])))
        ms_ps= []
        for m in ms:
            #axs[counter].scatter(p_values[m],rs[i][m]+bar_width,color = 'red',marker = "d")
            axs[counter].annotate(anno%(p_values[m]),[p_values[m],rs[i][m]])
            ms_ps.append(p_values[m]) 

        highest_correlate[pheno[i]]=ms_ps
        
        axs[counter].set_title(pheno[i])
        axs[counter].set_ylim([0,max(rs[i])*1.1])
        axs[counter].set_xlabel('p value threshold')
        axs[counter].set_ylabel('R square explained')
        

        counter+=1
    #axs[n_plot]

    cmap = mpl.cm.seismic
    #norm = mpl.colors.Normalize(vmin=0, vmax=1.0)
    
    bounds = [0.0,0.0001,0.001, 0.05,0.1,0.5, 1]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    
    cb1 = mpl.colorbar.ColorbarBase(axs[counter], cmap=cmap,
                                norm=norm,orientation='horizontal')
    axs[counter].set_xlabel('correlation p-value')
    #fig.text(0.5, 0.0, 'p value threshold', ha='center')
    #fig.text(0.0, 0.5, 'R Square from regression', va='center', rotation='vertical')
    #plt.colorbar()
    plt.savefig('{}.png'.format(outputName), bbox_inches='tight')
    print("R-squared plot saved to {}".format(os.path.basename('{}.png'.format(outputName))))
    return highest_correlate



"""
@author: zhuyuecai
@note: helper function, used in other functions
@param rs: a list of numerical values
@return: a list of the positions of the maximums in rs
"""


def _get_max_positions_(rs):
    m = max(rs)
    return [i for i, j in enumerate(rs) if j == m]



if __name__ == '__main__':
    import numpy as np
    pheno = ['t1','t2','t3','t1','t2','t3']
    prs = [[1,2,5,4,1,8,2],[6,7,5,4,1,8,2],[6,7,5,4,1,8,2],[1,2,5,4,1,8,2],[6,7,5,4,1,8,2],[6,7,5,4,1,8,2]]
    p_rs = [np.linspace(0,1,7),np.linspace(0,1,7),np.linspace(0,1,7),np.linspace(0,1,7),np.linspace(0,1,7),np.linspace(0,1,7)]
    ps = np.linspace(0,1,7)
    r_square_plots(pheno,prs,p_rs,ps,"Testing",width=3,bar_width=0.1)
