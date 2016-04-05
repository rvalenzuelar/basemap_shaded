
import matplotlib.pyplot as plt
import numpy as np
import rof_denoise as denoise

from matplotlib.colors import LightSource
from scipy.interpolate import interp1d, Rbf
from scipy import ndimage
from custom_cmap import make_cmap
from mpl_toolkits.basemap import Basemap, cm
from geographiclib.geodesic import Geodesic


class elevation:

    def __init__(self, file_name=None, domain_num=None,
                 linesec=None):
        self.file_name = file_name
        self.domain_num = domain_num
        self.dtm = None
        self.latn = None
        self.latx = None
        self.lonn = None
        self.lonx = None
        self.profile_line = None
        self.profile_elevation = None
        self.linesec = linesec
        self.line_x = None
        self.line_y = None

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
        param = []
        param.append([-124.0, -122.0, 39.5, 38.0, 10])
        param.append([-124.0, -122.9, 39.1, 38.1, 10])
        param.append([-123.3, -122.9, 38.8, 38.3, 8])
        param.append([-123.3, -123.0, 38.7, 38.4, 5])
        param.append([-123.8, -122.55, 39.1, 38.2, 10])

        lonn, lonx, latx, latn, sigma = param[self.domain_num]

        xini = int(fx(lonn))
        xend = int(fx(lonx))
        yini = int(fy(latx))
        yend = int(fy(latn))

        band = datafile.GetRasterBand(1)
        if self.dtm is None:
            source = band.ReadAsArray(xini, yini, (xend-xini), (yend-yini))
            dtm = gaussian_filter(source, sigma=sigma)
            self.dtm = np.flipud(dtm)
            self.latn = latn
            self.latx = latx
            self.lonn = lonn
            self.lonx = lonx

        if line is not None:
            profile = []
            for p in line:
                profile.append(getDtmElevation(p[1], p[0], band, geotransform))
            self.profile = np.array(profile)

    def plot_elevation_map(self, ax=None, shaded=True,
                           cmap=None, locations=None, colorbar=False,
                           blend_mode='overlay', grid=True,
                           altitude_range=None, contour_lines=None,
                           latdelta=None, londelta=None):

        import os
        homed = os.path.expanduser('~')

        ' get elevation model '
        fname = 'merged_dem_38-39_123-124_extended.tif'
        dtmfile = homed+'/Github/RadarQC/'+fname

        self.file_name = dtmfile
        if self.dtm is None:
            self.get_elevation()

        ' make map axis '
        m = Basemap(projection='merc',
                    llcrnrlat=self.latn,
                    urcrnrlat=self.latx,
                    llcrnrlon=self.lonn,
                    urcrnrlon=self.lonx,
                    resolution='c',
                    ax=ax)

        vmin, vmax = altitude_range

        if shaded:
            ' make hill shaded image '
            ls = LightSource(azdeg=15, altdeg=60)
            rgb = ls.shade(self.dtm, cmap=cmap, vmin=vmin, vmax=vmax,
                           blend_mode='soft', fraction=0.7)
            m.imshow(rgb)
        else:
            m.imshow(self.dtm, cmap=cmap, vmin=vmin, vmax=vmax)

        ' add parallels and meridians '
        parallels = np.arange(-90, 90, latdelta)
        meridians = np.arange(-180, 180, londelta)
        lw = 0
        if grid:
            lw = 1.
        m.drawparallels(parallels,
                        labels=[1, 0, 0, 0], fontsize=10,
                        labelstyle='+/-', fmt='%2.1f', linewidth=lw)
        m.drawmeridians(meridians,
                        labels=[0, 0, 0, 1], fontsize=10,
                        labelstyle='+/-', fmt='%2.1f', linewidth=lw)

        ' add locations '
        if locations:
            for loc, coord in locations.iteritems():
                x, y = m(*coord[::-1])
                m.scatter(x, y, 80, color='red')
                ax.text(x, y, loc, ha='right', va='bottom')

        ' add section line '
        if self.linesec is not None:
            x, y = get_line_ini_end(self)
            self.line_x = x
            self.line_y = y
            x, y = m(*[x, y])
            m.plot(x, y, color='k', linewidth=2)
            ndiv = self.linesec['ndiv']
            if ndiv > 0:
                xp = np.linspace(x[0], x[1], ndiv)
                yp = np.linspace(y[0], y[1], ndiv)
                # m.plot(x1, y1, marker='s', markersize=15, color='g')
                # m.plot(x2, y2, marker='s', markersize=15, color='r')
                m.plot(xp, yp, marker='s', markersize=5, color='k')

        ' add rivers '
        rivers = get_rivers(mmap=m)
        ax.add_collection(rivers)

        ' add colorbar '
        if colorbar:
            im = m.imshow(self.dtm, cmap=cmap, vmin=vmin, vmax=vmax)
            im.remove()
            cb = m.colorbar(im)

        ' add contour line(s) '
        if contour_lines is not None:
            nlats, nlons = self.dtm.shape
            lons, lats = np.meshgrid(
                np.linspace(self.lonn, self.lonx, nlons),
                np.linspace(self.latn, self.latx, nlats))
            x, y = m(lons, lats)
            m.contour(x, y, self.dtm, [800, 1000], colors='k')

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


def main(plot=False):
    from matplotlib import gridspec

    if plot:
        # arid = make_cmap(colors='arid', bit=True)
        # cold = make_cmap(colors='cold_humid', bit=True)
        warm = make_cmap(colors='warm_humid')

        loc1 = {'BBY': (38.32, -123.07), 'CZD': (38.61, -123.22)}
        loc2 = {'BBY': (38.32, -123.07), 'FRS': (38.51, -123.25)}
        linesec = {'origin': (38.29, -123.59),
                   'az': 50, 'dist': 110, 'ndiv': 0}
        elev = elevation(linesec=linesec, domain_num=4)

        fig = plt.figure(figsize=(8.5, 11))
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
        elev.plot_elevation_map(ax=ax0, cmap=warm, shaded=False,
                                locations=loc2, colorbar=True, grid=False,
                                altitude_range=[0, 800],
                                contour_lines=[800, 1000],
                                latdelta=0.1, londelta=0.2
                                )
        elev.plot_elevation_profile(ax=ax1)
        plt.show(block=False)


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


def get_line_ini_end(elev):

    linesec = elev.linesec
    if isinstance(linesec, dict):
        x1, y1 = linesec['origin'][::-1]
        x2, y2 = get_arrival_point(linesec)[::-1]
    elif isinstance(linesec, list):
        x1, y1 = linesec[0][::-1]
        x2, y2 = linesec[1][::-1]
    return [[x1, x2], [y1, y2]]


def get_arrival_point(linesec):

    lat, lon = linesec['origin']
    azi = linesec['az']
    dist = linesec['dist']*1000  # [m]
    gd = Geodesic.WGS84.Direct(lat, lon, azi, dist)

    return gd['lat2'], gd['lon2']


def get_rivers(mmap=None):

    import shapefile
    from matplotlib.patches import Polygon
    from matplotlib.collections import LineCollection

    shf = '/home/rvalenzuela/Github/basemap_shaded/sonoma_rivers/sonoma_rivers'
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


def plot_map(ax=None, domain=0, cmap=plt.cm.gist_earth, blend_mode='overlay'):

    dtmfile = '/home/raul/Github/RadarQC/merged_dem_38-39_123-124_extended.tif'

    dtm, _ = get_elevation(dtmfile=dtmfile, domain=domain)

    ls = LightSource(azdeg=15, altdeg=60)
    rgb = ls.shade(dtm, cmap=cmap, vmin=0, vmax=1000,
                   blend_mode='soft', fraction=0.7)
    ax.imshow(rgb)

    'Use a proxy artist for the colorbar'
    im = ax.imshow(dtm, cmap=cmap, vmin=0, vmax=1000,)
    im.remove()
    plt.colorbar(im)


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
