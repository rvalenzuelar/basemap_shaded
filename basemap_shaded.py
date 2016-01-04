
import gdal
import matplotlib.pyplot as plt
import numpy as np
import rof_denoise as denoise

from matplotlib.colors import LightSource
from scipy.interpolate import interp1d, Rbf
from scipy import ndimage
from custom_cmap import make_cmap
from mpl_toolkits.basemap import Basemap,cm


def main():

	# arid = make_cmap(colors='arid', bit=True)	
	# cold = make_cmap(colors='cold_humid', bit=True)	
	warm = make_cmap(colors='warm_humid')	

	# fig,ax=plt.subplots(1,3,sharey=True)
	# plot_map(ax=ax[0],domain=2,cmap=arid)
	# plot_map(ax=ax[1],domain=2,cmap=cold)
	# plot_map(ax=ax[2],domain=2,cmap=warm)
	# plt.show()


	fig,ax=plt.subplots()
	plot_geomap(ax=ax,domain=0,cmap=warm)
	plt.show()


def plot_geomap(ax=None,domain=0,cmap=None,blend_mode='overlay'):


	dtmfile='/home/raul/Github/RadarQC/merged_dem_38-39_123-124_extended.tif'

	dtm,geobound = get_elevation(dtmfile=dtmfile,domain=domain)
	dtm=np.flipud(dtm)

	m=Basemap(projection='merc',
					llcrnrlat=geobound[0],
					urcrnrlat=geobound[1],
					llcrnrlon=geobound[2],
					urcrnrlon=geobound[3],
					resolution='c',
					ax=ax)
	
	ls = LightSource(azdeg=15,altdeg=60)
	rgb=ls.shade(dtm,cmap=cmap,vmin=0,vmax=1000,
				blend_mode='soft',fraction=0.7)
	m.imshow(rgb)
	
	rivers=get_rivers(mmap=m)
	# m.plot(rivers)
	ax.add_collection(rivers)


	'Use a proxy artist for the colorbar'
	im = m.imshow(dtm, cmap=cmap, vmin=0,vmax=1000,)
	im.remove()
	plt.colorbar(im)


def get_rivers(mmap=None):

	import shapefile
	from matplotlib.patches import Polygon
	# from matplotlib.collections import PatchCollection
	from matplotlib.collections import LineCollection
	
	s=shapefile.Reader('./sonoma_rivers/sonoma_rivers')

	shapes=s.shapes()
	Nshp=len(shapes)

	# lines=[]
	# for n in range(Nshp):
	# 	ptchs   = []
	# 	pts     = np.asarray(shapes[n].points)
	# 	prt     = shapes[n].parts
	# 	par     = list(prt) + [pts.shape[0]]
	# 	for pij in xrange(len(prt)):
	# 		ptchs.append(Polygon(pts[par[pij]:par[pij+1]]))
	# 	lines.append(PatchCollection(ptchs,facecolor='blue',edgecolor='blue', linewidths=.1))
	segs=[]
	for n in [0]:
		pts = shapes[n].points
		# print pts
		lons,lats=zip(*pts)
		# print lons
		x, y = mmap(lons,lats)
		segs.append(zip(x,y))
	lines=LineCollection(segs)
	lines.set_facecolors('b')
	lines.set_edgecolors('b')
	lines.set_linewidth(0.3)

	return lines


def plot_map(ax=None,domain=0,cmap=plt.cm.gist_earth,blend_mode='overlay'):


	dtmfile='/home/raul/Github/RadarQC/merged_dem_38-39_123-124_extended.tif'

	dtm,_ = get_elevation(dtmfile=dtmfile,domain=domain)

	
	ls = LightSource(azdeg=15,altdeg=60)
	rgb=ls.shade(dtm,cmap=cmap,vmin=0,vmax=1000,
				blend_mode='soft',fraction=0.7)
	ax.imshow(rgb)
	
	'Use a proxy artist for the colorbar'
	im = ax.imshow(dtm, cmap=cmap, vmin=0,vmax=1000,)
	im.remove()
	plt.colorbar(im)


def get_elevation(dtmfile=None,domain=None):

	''' store dtm in data '''
	datafile = gdal.Open(dtmfile)
	geotransform=datafile.GetGeoTransform()
	cols=datafile.RasterXSize
	rows=datafile.RasterYSize

	''' geographic axes '''
	originX=geotransform[0]
	originY=geotransform[3]
	pixelW=geotransform[1]
	pixelH=geotransform[5]

	endingX=originX+cols*pixelW
	endingY=originY+rows*pixelH

	xg=np.linspace(originX,endingX,cols)
	yg=np.linspace(originY,endingY,rows)

	fx=interp1d(xg,range(len(xg)))
	fy=interp1d(yg,range(len(yg)))

	if domain == 0:
		lonn=-124.0
		lonx=-122.9
		latx=39.1
		latn=38.1
		sigma=10
	elif domain == 1:
		lonn=-123.3
		lonx=-122.9
		latx=38.8
		latn=38.3
		sigma=8
	elif domain == 2:
		lonn=-123.3
		lonx=-123.0
		latx=38.7
		latn=38.4
		sigma=5

	xini=int(fx(lonn))
	xend=int(fx(lonx))
	yini=int(fy(latx))
	yend=int(fy(latn))

	band=datafile.GetRasterBand(1)		
	source=band.ReadAsArray(xini,yini,(xend-xini),(yend-yini))

	dtm = gaussian_filter(source,sigma=sigma)

	return dtm,[latn,latx,lonn,lonx]

def interp_rbf(array,x,y,res=10):
	
	xm,ym = np.meshgrid(x,y)

	xx=xm[::res,::res]
	yy=ym[::res,::res]
	zz=array[::res,::res]

	f=Rbf(xx,yy,zz)

	zi=f(xm,ym)
	return zi

def gaussian_filter(array,sigma=10):

	filtered = ndimage.filters.gaussian_filter(array,sigma)

	return filtered

def hp_filter(array):

	# A very simple and very narrow highpass filter
	kernel = np.array([[-1, -1, -1],
	                   [-1,  8, -1],
	                   [-1, -1, -1]])
	hp = ndimage.convolve(array, kernel)
	filtered=array-hp
	return filtered


main()