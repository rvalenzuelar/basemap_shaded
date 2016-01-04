
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

	loc1={'BBY':(38.32,-123.07), 'CZD':(38.61,-123.22)}
	loc2={'BBY':(38.32,-123.07), 'CZD':(38.61,-123.22), 'FRS':(38.51,-123.25)}

	# fig,ax=plt.subplots(figsize=(10,12))
	# plot_geomap(ax=ax,domain=1,cmap=warm, locations=loc2,colorbar=True)
	# plt.show(block=False)

	fig,ax=plt.subplots(1,2,figsize=(12,10))
	plot_geomap(ax=ax[0],domain=0,cmap=warm, locations=loc2)
	plot_geomap(ax=ax[1],domain=1,cmap=warm, locations=loc2,colorbar=True)
	plt.subplots_adjust(wspace=0.05)
	plt.show(block=False)


def plot_geomap(ax=None,domain=0,cmap=None,locations=None, colorbar=False,blend_mode='overlay'):


	# dtmfile='/home/raul/Github/RadarQC/merged_dem_38-39_123-124_extended.tif'
	dtmfile='/home/rvalenzuela/Github/RadarQC/merged_dem_38-39_123-124_extended.tif'

	dtm,geobound = get_elevation(dtmfile=dtmfile,domain=domain)
	dtm=np.flipud(dtm)

	m=Basemap(projection='merc',
					llcrnrlat=geobound[0],
					urcrnrlat=geobound[1],
					llcrnrlon=geobound[2],
					urcrnrlon=geobound[3],
					resolution='i',
					ax=ax)
	
	ls = LightSource(azdeg=15,altdeg=60)
	rgb=ls.shade(dtm,cmap=cmap,vmin=0,vmax=800,
				blend_mode='soft',fraction=0.7)
	m.imshow(rgb)
	# m.drawcoastlines()

	if domain == 0:
		deg_delta=0.2
	else:
		deg_delta=0.1
	parallels=np.arange(-90,90,deg_delta)
	meridians=np.arange(-180,180,deg_delta)
	m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10,labelstyle='+/-', fmt='%2.1f')
	m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10,labelstyle='+/-', fmt='%2.1f')

	if locations:
		for loc,coord in locations.iteritems():
			x,y = m(*coord[::-1])
			m.scatter(x,y, 80,color='red')
			plt.text(x,y,loc,ha='right',va='top')

	rivers=get_rivers(mmap=m)
	ax.add_collection(rivers)

	if colorbar:
		im = m.imshow(dtm, cmap=cmap, vmin=0,vmax=800)
		im.remove()
		cb=m.colorbar(im)

def get_rivers(mmap=None):

	import shapefile
	from matplotlib.patches import Polygon
	# from matplotlib.collections import PatchCollection
	from matplotlib.collections import LineCollection
	
	s=shapefile.Reader('./sonoma_rivers/sonoma_rivers')

	shapes=s.shapes()
	Nshp=len(shapes)

	segs=[]
	for n in range(Nshp):
		pts = shapes[n].points
		lons,lats=zip(*pts)
		x, y = mmap(lons,lats)
		segs.append(zip(x,y))
	lines=LineCollection(segs)
	lines.set_edgecolors('b')
	lines.set_linewidth(1)

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