import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib import patches as mpatches
from matplotlib import colors as mcolors
#from matplotlib import colors
from matplotlib.ticker import PercentFormatter

def cmap_from_txt_file(file):
    return pd.read_csv(file, header=None, skiprows=1,names=['value','red','green','blue','alpha','label'])

def ctable_to_LScolormap(ctable, skipFirstLine=True):
    # list of positions (values) where colours are assigned to
    if skipFirstLine:
        rgb_list=(list(ctable['red'][1:]/255.),list(ctable['green'][1:]/255.),list(ctable['blue'][1:]/255.))
        values=ctable['value'][1:] # skip the no-data value
    else:
        rgb_list=(list(ctable['red']/255.),list(ctable['green']/255.),list(ctable['blue']/255.))
        values=ctable['value']
    
    vmin=values.astype(np.float64).min()
    vmax=values.astype(np.float64).max()
     
    vrange = vmax-vmin
 
    # pos_list should position the colours in the range from 0 to 1
    #   for PPI (0...30K), with no offset, we divide by vrange (30K)
    #   for NDVI (-0.08...0.92), with offset, we first subtract vmin value (-0.08) and then divide by vrange (1)

    voffset = values.astype(np.float64) - vmin
    pos_list = round(voffset / vrange,5).tolist()  
    
    nodatacolour=(rgb_list[0][0],rgb_list[1][0],rgb_list[2][0])
    
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[pos_list[i], rgb_list[num][i], rgb_list[num][i]] for i in range(len(pos_list))]
        cdict[col] = col_list
       
    mycmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=30000)
    mycmp.set_bad(color=(220/255.,220/255.,220/255.))
    
    return mycmp