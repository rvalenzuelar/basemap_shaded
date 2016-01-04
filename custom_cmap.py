
'''
Chris Slocum
http://schubert.atmos.colostate.edu/~cslocum/custom_cmap.html

    
make_cmap takes a list of tuples which contain RGB values. The RGB
values may either be in 8-bit [0 to 255] (in which bit must be set to
True when called) or arithmetic [0 to 1] (default). make_cmap returns
a cmap with equally spaced colors.
Arrange your tuples so that the first color is the lowest value for the
colorbar and the last is the highest.
position contains values from 0 to 1 to dictate the location of each color. 
'''

import matplotlib as mpl
import numpy as np

def make_cmap(colors=None, position=None, bit=False):

    if isinstance(colors,str):
        colors=tint(name=colors)
        sea=np.asarray([0])
        land=np.linspace(0.01,1,len(colors)-1)
        position=np.append(sea,land)
        bit=True
   
    if position is None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")

    bit_rgb = np.linspace(0,1,256)
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap

def tint(name=None):

    if name == 'arid':
        colors=[(245,245,245),(235,235,237),(220,220,220),(212,207,204),
                (212,193,179),(212,184,163),(212,201,180),(202,190,174),
                (180,170,158),(170,160,150),(160,152,141),(146,136,129),]
        colors.reverse()

    elif name == 'warm_humid':
        colors=[(245,245,245),(235,235,237),(220,220,220),(212,207,204),
                (212,193,179),(212,184,163),(212,201,180),(169,192,166),
                (134,184,159),(120,172,149),(114,164,141),(106,153,135),
                (152,221,250)]
        colors.reverse()

    elif name == 'cold_humid':
        colors=[(245,245,245),(235,235,237),(220,220,220),(212,207,204),
                (212,193,179),(212,184,163),(212,201,180),(180,192,180),
                (145,177,171),(130,165,159),(120,159,152),(112,147,141)]
        colors.reverse()


    return colors

