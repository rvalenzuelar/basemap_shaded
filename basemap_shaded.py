
import matplotlib.pyplot as plt
import numpy as np
# import rof_denoise as denoise
import Circlem

from matplotlib.colors import LightSource
from scipy.interpolate import interp1d, Rbf
from scipy import ndimage
from custom_cmap import make_cmap
from mpl_toolkits.basemap import Basemap
from geographiclib.geodesic import Geodesic
from rv_utilities import add_colorbar


locations = {'BBY': (38.32, -123.07),
             'FRS': (38.52, -123.25),
             'CZD': (38.61, -123.22),
             'Petaluma': (38.232, -122.636),
             'SCK': (37.93, -121.22)
             }

class elevation:

    def __init__(self, file_name=None, domain_num=None,
                 domain=None, sigma=None, source=None):
        
        self.file_name = file_name
        self.domain_num = domain_num
        self.domain = domain
        self.sigma = sigma
        self.source = source
        self.dtm = None
        self.latn = None
        self.latx = None
        self.lonn = None
        self.lonx = None
        self.profile_line = None
        self.profile_elevation = None
        self.line_x = None
        self.line_y = None
        self.lats = None
        self.lons = None

    def get_elevation(self, line=None):

        import gdal

        ''' store dtm in data '''
        datafile = gdal.Open(self.file_name)
        geotransform = datafile.GetGeoTransform()
        self.geotransform = geotransform
        cols = datafile.RasterXSize
        rows = datafile.RasterYSize

        ''' geographic axes '''
        originX = geotransform[0]
        originY = geotransform[3]
        pixelW = geotransform[1]
        pixelH = geotransform[5]

        endingX = originX+cols*pixelW
        endingY = originY+rows*pixelH

        xg = np.linspace(originX, endingX, cols)
        yg = np.linspace(originY, endingY, rows)
        
        fx = interp1d(xg, range(len(xg)))
        fy = interp1d(yg, range(len(yg)))

        ' set of domains to plot '
        if self.domain_num is not None:
            param = []
            param.append([-124.0, -122.0, 39.5, 38.0, 10])
            param.append([-124.0, -122.9, 39.1, 38.1, 10])
            param.append([-123.3, -122.9, 38.8, 38.3, 8])
            param.append([-123.3, -123.0, 38.7, 38.4, 5])
            param.append([-123.8, -122.55, 39.1, 38.2, 10])
            lonn, lonx, latx, latn, sigma = param[self.domain_num]
        else:
            lonn, lonx, latx, latn, sigma = self.domain

        xini = int(fx(lonn))
        xend = int(fx(lonx))
        yini = int(fy(latx))
        yend = int(fy(latn))

        band = datafile.GetRasterBand(1)
        if self.dtm is None:
            source = band.ReadAsArray(xini, yini, (xend-xini), (yend-yini))
            dtm = gaussian_filter(source, sigma=sigma)
            nlats,nlons = dtm.shape
            self.dtm = np.flipud(dtm)
            self.latn = latn
            self.latx = latx
            self.lonn = lonn
            self.lonx = lonx
            self.lats = np.linspace(latn,latx,nlats)
            self.lons = np.linspace(lonn,lonx,nlons)


        if line is not None:
            profile = []
            for p in line:
                profile.append(getDtmElevation(p[1], p[0], band, geotransform))
            self.profile = np.array(profile)

    def get_dem_matfile(self,resample=None):

        import scipy.io as sio
        from scipy.interpolate import RectBivariateSpline as spline
        
        mat = sio.loadmat(self.file_name)
        dtm = np.flipud(mat['Z1'])
        rows, cols = dtm.shape

        latn, latx = mat['latlim'][0]
        lonn, lonx = mat['lonlim'][0]
        xlat = np.linspace(latn, latx, rows)
        interpLat = interp1d(xlat, range(rows))
        xlon = np.linspace(lonn, lonx, cols)
        interpLon = interp1d(xlon, range(cols))

        inX0 = int(interpLon(self.domain[0]))
        inX1 = int(interpLon(self.domain[1]))
        inY0 = int(interpLat(self.domain[2]))
        inY1 = int(interpLat(self.domain[3]))
        dtm_out = dtm[inY0:inY1, inX0:inX1]
        
        if resample is None:
            self.dtm = dtm_out
        else:
            rows,cols = dtm_out.shape
            spl = spline(range(rows),range(cols),dtm_out)            
            ri = np.linspace(0,rows,resample[0])
            ci = np.linspace(0,cols,resample[1])
            Xi,Yi = np.meshgrid(ci,ri)
            xif,yif = Xi.flatten(), Yi.flatten()
            zi = np.zeros((resample[0],resample[1])).flatten()
            for n,x,y in zip(range(zi.size),yif,xif):
                zi[n] = spl.ev(x,y)
            self.dtm = np.reshape(zi,(resample[0],resample[1]))

        self.latn, self.latx = self.domain[2:]
        self.lonn, self.lonx = self.domain[:2]
            

    def plot_elevation_map(self, ax=None, shaded=True,
                           cmap=None, add_loc=None,
                           colorbar=True, blend_mode='overlay',
                           grid=True, gridcolor='w',
                           altitude_range=[-10,600],
                           contour_lines=None,
                           latdelta=None, londelta=None,
                           homed=None, figsize=None,
                           addrivers=False,
                           add_geoline=None):

        import os
        
        if homed is None:
            homed = os.path.expanduser('~')


        ' get elevation model '
        if self.source is not None:
            fname = self.source
            dtmfile = homed + '/' + fname
            self.file_name = dtmfile
            if self.dtm is None:
                self.get_dem_matfile()
        else:
            fname = 'merged_dem_38-39_123-124_extended.tif'
            dtmfile = homed + '/Github/RadarQC/' + fname
            self.file_name = dtmfile
            if self.dtm is None:
                self.get_elevation()

        if ax is None:
            if figsize is None:
                fig,ax = plt.subplots(figsize=(10,10))
            else:
                fig,ax = plt.subplots(figsize=figsize)

        ' make map axis '
        m = Basemap(projection='merc',
                    llcrnrlat=self.latn,
                    urcrnrlat=self.latx,
                    llcrnrlon=self.lonn,
                    urcrnrlon=self.lonx,
                    resolution='h',
                    ax=ax)

        vmin, vmax = altitude_range

        if shaded:
            ' make hill shaded image '
            ls = LightSource(azdeg=15, altdeg=60)
            rgb = ls.shade(self.dtm,
                           cmap=getattr(plt.cm,cmap),
                           vmin=vmin,
                           vmax=vmax,
                           blend_mode='soft',
                           fraction=0.7)
            m.imshow(rgb)
        else:
            m.imshow(self.dtm, cmap=cmap, vmin=vmin, vmax=vmax)

        ' add parallels and meridians '
        parallels = np.arange(-90, 90, latdelta)
        meridians = np.arange(-180, 180, londelta)
        lw = 0
        if grid:
            lw = 1.
        parallels = m.drawparallels(parallels,
                        labels=[1, 0, 0, 0], fontsize=10,
                        labelstyle='+/-', fmt='%2.1f',linewidth=lw,
                        color=gridcolor)
        m.drawmeridians(meridians,
                        labels=[0, 0, 0, 1], fontsize=10,
                        labelstyle='+/-', fmt='%2.1f', linewidth=lw,
                        color=gridcolor)

        for p in parallels:
            try:
                parallels[p][1][0].set_rotation(90)
            except:
                pass

        ' add locations '
        if add_loc:
            fsize = 15
            psize = 50
            ec = (1.0,0,0,1)
            fc = (0.5,0,0,1)
            if isinstance(add_loc,dict):
                for loc, coord in locations.iteritems():
                    x, y = m(*coord[::-1])
                    m.scatter(x, y, psize, color='r')
                    ax.text(x, y, loc, ha='right', va='bottom',
                            color='r',fontsize=fsize,weight='bold')
            elif isinstance(add_loc,list):
                for loc in add_loc:
                    coord = locations[loc]
                    x, y = m(*coord[::-1])
                    m.scatter(x, y, psize, facecolor=fc,edgecolor=ec)
                    ax.text(x, y, loc, ha='right', va='bottom',
                            color='r',fontsize=fsize,weight='bold')                

        ' add section line '
        if add_geoline:
            if isinstance(add_geoline,dict):
                geolines = list(add_geoline)
            elif isinstance(add_geoline,list):
                geolines = add_geoline
                
            for line in geolines:                
                x, y = get_line_ini_end(line)
                self.line_x = x
                self.line_y = y
                x, y = m(*[x, y])
    
                if 'color' in line:
                    color = line['color']
                else:
                    color='k'
                    
                m.plot(x, y, color=color, linewidth=2)
    
                if 'ndiv' in line:
                    ndiv = line['ndiv']
                    xp = np.linspace(x[0], x[1], ndiv)
                    yp = np.linspace(y[0], y[1], ndiv)
                    m.plot(xp, yp, marker='s',
                           markersize=5, color=color)
                
                if 'label' in line:
                    if 'center' in line:
                        x,y = m(*line['center'][::-1])
                        ax.text(x,y,line['label']+r'$^\circ$',
                                color=color,
                                weight='bold',
                                fontsize=15,
                                rotation=90-line['az'])

        ' add rivers '
        if addrivers:
            rivers = get_rivers(mmap=m)
            ax.add_collection(rivers)

        ' add colorbar '
        if colorbar:
            im = m.imshow(self.dtm, cmap=cmap, vmin=vmin, vmax=vmax)
            im.remove()
            add_colorbar(ax,im,label='Meters')

        ' add contour line(s) '
        if contour_lines is not None:
            nlats, nlons = self.dtm.shape
            lons, lats = np.meshgrid(
                np.linspace(self.lonn, self.lonx, nlons),
                np.linspace(self.latn, self.latx, nlats))
            x, y = m(lons, lats)
            m.contour(x, y, self.dtm, contour_lines, colors='k')

        ' add coastline'
        m.drawcoastlines(color='w')


    def plot_elevation_profile(self, ax=None, npoints=500):

        x = self.line_x
        y = self.line_y
        c0 = [y[0], x[0]]
        c1 = [y[1], x[1]]
        line = interpolateLine(c0, c1, npoints)
        self.get_elevation(line)
        y = self.profile.astype(float)
        y[y == 0] = np.nan
        dist = self.linesec['dist']
        x = np.linspace(0, dist, len(y))
        ax.plot(x, y, color='k', lw=2)
        ax.set_xlim([0, dist])
        ax.set_xlabel('Distance [km]')
        ax.set_ylabel('Altitude [m]')


def plot_terrain_profile():

    from matplotlib import gridspec
    warm = make_cmap(colors='warm_humid')

    linesec = {'origin': (38.29, -123.59),
               'az': 50, 'dist': 110, 'ndiv': 0}
               
    elev = elevation(domain_num=4)

    plt.figure(figsize=(8.5, 11))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    elev.plot_elevation_map(ax=ax0, cmap=warm, shaded=False,
                            add_loc=['BBY','FRS','CZD'],
                            colorbar=True, grid=False,
                            altitude_range=[0, 800],
                            contour_lines=[800, 1000],
                            latdelta=0.1, londelta=0.2,
                            add_geoline=linesec,
                            )
    elev.plot_elevation_profile(ax=ax1)


def plot_obs_domain(ax=None, cmap=None):

    if cmap is None:
        warm = make_cmap(colors='warm_humid')

    linesec = {'origin': (38.29, -123.59),
               'az': 50, 'dist': 110, 'ndiv': 0}

    if ax is None:
        fig, ax = plt.subplots()

    elev = elevation(linesec=linesec,
                     domain=[-124.0, -122.0, 37.8, 39.38],
                     source='BNCaliforniaDEM.mat')

    m = elev.plot_elevation_map(ax=ax, cmap=cmap, shaded=False,
                                add_loc=['FRS','BBY','CZD','Petaluma'],
                                colorbar=True,
                                add_geoline=linesec,
                                grid=False,
                                altitude_range=[-5, 800],
                                contour_lines=[800],
                                latdelta=0.2, londelta=0.2,
                                addrivers=False,
                                homed='/home/raul/Dropbox/NOCAL_DEM')

    add_rings(ax, space_km=10, color='k', mapping=[m, 38.51, -123.25])


def plot_petaluma_gap(ax=None,cmap=None):

    if cmap is None:
        warm = make_cmap(colors='warm_humid')

    if ax is None:
        fig, ax = plt.subplots()

    elev = elevation(domain=[-123.5, -121, 37.7, 38.7],
                     source='NCalDEMforGapFlow.mat')
    elev.plot_elevation_map(ax=ax, cmap=cmap, shaded=False,
                            locations=locs, colorbar=False, grid=False,
                            altitude_range=[-9, 800],
                            contour_lines=[800],
                            latdelta=0.2, londelta=0.2,
                            addrivers=False,
                            homed='/home/raul/Dropbox/NOCAL_DEM')


def interpolateLine(start_point, finish_point, number_points):
    gd = Geodesic.WGS84.Inverse(start_point[0], start_point[1],
                                finish_point[0], finish_point[1])
    line = Geodesic.WGS84.Line(gd['lat1'], gd['lon1'], gd['azi1'])
    line_points = []
    for i in range(number_points):
        point = line.Position(gd['s12'] / number_points * i)
        line_points.append((point['lat2'], point['lon2']))

    return line_points


def getDtmElevation(x, y, band, gt):
    col = []
    px = int((x - gt[0]) / gt[1])
    py = int((y - gt[3]) / gt[5])
    win_xsize, win_ysize = [1, 1]
    data = band.ReadAsArray(px, py, win_xsize, win_ysize)
    col.append(data[0][0])
    col.append(0)
    return col[0]


def get_line_ini_end(geoline):

    if isinstance(geoline, dict):
        if 'origin' in geoline:
            x1, y1 = geoline['origin'][::-1]
            azi = geoline['az']
            dist = geoline['dist']*1000  # [m]            
            x2, y2 = get_arrival_point(y1,x1,azi,dist)[::-1]
        elif 'center' in geoline:
            x0, y0 = geoline['center'][::-1]
            azi = geoline['az']
            dist = geoline['dist']*1000  # [m]            
            x1, y1 = get_arrival_point(y0,x0,azi,dist)[::-1]            
            x2, y2 = get_arrival_point(y0,x0,azi+180,dist)[::-1]            
            
    elif isinstance(geoline, list):
        x1, y1 = geoline[0][::-1]
        x2, y2 = geoline[1][::-1]
        
    return [[x1, x2], [y1, y2]]


def get_arrival_point(lat,lon,azi,dist):

    gd = Geodesic.WGS84.Direct(lat, lon, azi, dist)

    return gd['lat2'], gd['lon2']


def get_rivers(mmap=None):

    import shapefile
#    from matplotlib.patches import Polygon
    from matplotlib.collections import LineCollection

    shf = '/home/raul/Github/basemap_shaded/sonoma_rivers/sonoma_rivers'
    s = shapefile.Reader(shf)

    shapes = s.shapes()
    Nshp = len(shapes)

    segs = []
    for n in range(Nshp):
        pts = shapes[n].points
        lons, lats = zip(*pts)
        x, y = mmap(lons, lats)
        segs.append(zip(x, y))
    lines = LineCollection(segs)
    lines.set_edgecolors('b')
    lines.set_linewidth(1)

    return lines


def plot_map(ax=None, domain=2, cmap='gist_earth', blend_mode='overlay'):


    if ax is None:
        fig,ax=plt.subplots(figsize=(10,10))

#    dtmfile = '/home/raul/Github/RadarQC/merged_dem_38-39_123-124_extended.tif'
    dtmfile = '/home/raul/Dropbox/NOCAL_DEM/merged_dem_38-39_123-124_extended.tif'

    elev = elevation(file_name = dtmfile,
                     domain    = [-123.4, -122.8, 38.8, 38.1, 8])

    elev.get_elevation()

    elev_lims=[0,700]

    ls = LightSource(azdeg=15, altdeg=60)
    rgb = ls.shade(elev.dtm,
                   cmap=getattr(plt.cm,cmap),
                   vmin=elev_lims[0],
                   vmax=elev_lims[1],
                   blend_mode='soft',
                   fraction=0.7)
    
    ax.imshow(rgb)
    
    'Use a proxy artist for the colorbar'
    im = ax.imshow(elev.dtm,
                   cmap=cmap,
                   vmin=elev_lims[0],
                   vmax=elev_lims[1],
                   origin='lower',
                   )
    
    im.remove()
    add_colorbar(ax,im,label='Meters')


def rgb2gray(rgb):

    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray

def interp_rbf(array, x, y, res=10):

    xm, ym = np.meshgrid(x, y)

    xx = xm[::res, ::res]
    yy = ym[::res, ::res]
    zz = array[::res, ::res]

    f = Rbf(xx, yy, zz)

    zi = f(xm, ym)
    return zi


def gaussian_filter(array, sigma=10):

    filtered = ndimage.filters.gaussian_filter(array, sigma)
    return filtered


def hp_filter(array):

    # A very simple and very narrow highpass filter
    kernel = np.array([[-1, -1, -1],
                       [-1,  8, -1],
                       [-1, -1, -1]])
    hp = ndimage.convolve(array, kernel)
    filtered = array-hp
    return filtered


def add_rings(ax, space_km=10, color='k', mapping=False):

    # textdirection=225
    textdirection = -5

    for r in range(0, 60 + space_km, space_km):

        ring = add_ring(ax=ax, radius=r, mapping=mapping, color=color)
        vert = ring.get_path().vertices
        x, y = vert[textdirection]
        # ax.text(x * r, y * r, str(r), ha='center', va='center',
        #         bbox=dict(fc='none', ec='none', pad=2.),
        #         clip_on=True)


def add_ring(ax=None, radius=None, mapping=False, color=None, lw=1):

    from shapely.geometry import Polygon
    from descartes import PolygonPatch

    if mapping:
        m = mapping[0]
        olat = mapping[1]
        olon = mapping[2]
        c = Circlem.circle(m, olat, olon, radius * 1000.)
        circraw = Polygon(c)
        circ = PolygonPatch(circraw, fc='none', ec=color)
    else:
        circ = plt.Circle((0, 0), radius, fill=False,
                          color=color, linewidth=lw)

    ax.add_patch(circ)

    return circ
