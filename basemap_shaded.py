
import gdal
import matplotlib.pyplot as plt
import numpy as np
import rof_denoise as denoise

from matplotlib.colors import LightSource
from scipy.interpolate import interp1d, Rbf
from scipy import ndimage
from custom_cmap import make_cmap
from mpl_toolkits.basemap import Basemap,cm


def main(plot=False):

	if plot:
		# arid = make_cmap(colors='arid', bit=True)	
		# cold = make_cmap(colors='cold_humid', bit=True)	
		warm = make_cmap(colors='warm_humid')	


		loc1={'BBY':(38.32,-123.07), 'CZD':(38.61,-123.22)}
		loc2={'BBY':(38.32,-123.07), 'CZD':(38.61,-123.22), 'FRS':(38.51,-123.25)}

		linesec={'origin':(38.1,-122.95), 'az':338.0, 'dist':160}

		fig,ax=plt.subplots(figsize=(10,12))
		plot_geomap(ax=ax,domain=0,cmap=warm, locations=loc2,colorbar=True,linesec=linesec)
		plt.show(block=False)

		# fig,ax=plt.subplots(1,2,figsize=(12,10))
		# plot_geomap(ax=ax[0],domain=0,cmap=warm, locations=loc2)
		# # plot_geomap(ax=ax[1],domain=1,cmap=warm, locations=loc2,colorbar=True)
		# plt.subplots_adjust(wspace=0.05)
		# plt.show(block=False)


def plot_geomap(ax=None,shaded=True,domain=3,cmap=None,locations=None, 
					colorbar=False,blend_mode='overlay',linesec=None):

	' get elevation model '
	# dtmfile='/home/raul/Github/RadarQC/merged_dem_38-39_123-124_extended.tif'
	dtmfile='/home/rvalenzuela/Github/RadarQC/merged_dem_38-39_123-124_extended.tif'
	dtm,geobound = get_elevation(dtmfile=dtmfile,domain=domain)
	dtm=np.flipud(dtm)

	' make map axis '
	m=Basemap(projection='merc',
					llcrnrlat=geobound[0],
					urcrnrlat=geobound[1],
					llcrnrlon=geobound[2],
					urcrnrlon=geobound[3],
					resolution='c',
					ax=ax)
	
	if shaded:
		' make hill shaded image '
		ls = LightSource(azdeg=15,altdeg=60)
		rgb=ls.shade(dtm,cmap=cmap,vmin=0,vmax=800,
					blend_mode='soft',fraction=0.7)
		m.imshow(rgb)
	else:
		m.imshow(dtm,cmap=cmap,vmin=0,vmax=800)

	' add parallels and meridians '
	if domain == 0:
		deg_delta=0.3
	elif domain == 1:
		deg_delta=0.2
	else:
		deg_delta=0.1
	parallels=np.arange(-90,90,deg_delta)
	meridians=np.arange(-180,180,deg_delta)
	m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10,labelstyle='+/-', fmt='%2.1f')
	m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10,labelstyle='+/-', fmt='%2.1f')

	' add locations '
	if locations:
		for loc,coord in locations.iteritems():
			x,y = m(*coord[::-1])
			m.scatter(x,y, 80,color='red')
			plt.text(x,y,loc,ha='right',va='top')

	' add section line '
	if isinstance(linesec,dict):
		pass
	elif isinstance(linesec,list):
		x1,y1=m(*linesec[0][::-1])
		x2,y2=m(*linesec[1][::-1])
		x,y=[[x1,x2], [y1,y2]]
		xp=np.linspace(x1,x2,9)
		yp=np.linspace(y1,y2,9)
		m.plot(x1,y1,marker='s',markersize=15,color='g')
		m.plot(x2,y2,marker='s',markersize=15,color='r')
		m.plot(xp,yp,marker='s',markersize=5,color='k')
		m.plot(x,y,color='k',linewidth=2)

	' add rivers '
	rivers=get_rivers(mmap=m)
	ax.add_collection(rivers)

	' add colorbar '
	if colorbar:
		im = m.imshow(dtm, cmap=cmap, vmin=0,vmax=800)
		im.remove()
		cb=m.colorbar(im)




def get_rivers(mmap=None):

	import shapefile
	from matplotlib.patches import Polygon
	from matplotlib.collections import LineCollection
	
	shf='/home/rvalenzuela/Github/basemap_shaded/sonoma_rivers/sonoma_rivers'
	s=shapefile.Reader(shf)

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

	param=[]
	param.append([-124.0, -122.0, 39.5, 38.0, 10])	
	param.append([-124.0, -122.9, 39.1, 38.1, 10])	
	param.append([-123.3, -122.9, 38.8, 38.3, 8])
	param.append([-123.3, -123.0, 38.7, 38.4, 5])
	
	lonn,lonx,latx,latn,sigma=param[domain]

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