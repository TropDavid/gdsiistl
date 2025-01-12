# PPLN ring test layot v01
# 10/11/2018 - Boris Desiatov - Harvard
##
#
#
import os
import numpy as np
import gdspy
import math
from scipy import special
import uuid
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath



exportname = "spiral_test_02"

chipW=6000
chipH=6000

def render_text(text, size=None, position=(0, 0), font_prop=None, tolerance=0.1):
    path = TextPath(position, text, size=size, prop=font_prop)
    polys = []
    xmax = position[0]
    for points, code in path.iter_segments():
        if code == path.MOVETO:
            c = gdspy.Curve(*points, tolerance=tolerance)
        elif code == path.LINETO:
            c.L(*points)
        elif code == path.CURVE3:
            c.Q(*points)
        elif code == path.CURVE4:
            c.C(*points)
        elif code == path.CLOSEPOLY:
            poly = c.get_points()
            if poly.size > 0:
                if poly[:, 0].min() < xmax:
                    i = len(polys) - 1
                    while i >= 0:
                        if gdspy.inside(poly[:1], [polys[i]], precision=0.1 * tolerance)[0]:
                            p = polys.pop(i)
                            poly = gdspy.boolean([p], [poly], 'xor', precision=0.1 * tolerance,
                                                 max_points=0).polygons[0]
                            break
                        elif gdspy.inside(polys[i][:1], [poly], precision=0.1 * tolerance)[0]:
                            p = polys.pop(i)
                            poly = gdspy.boolean([p], [poly], 'xor', precision=0.1 * tolerance,
                                                 max_points=0).polygons[0]
                        i -= 1
                xmax = max(xmax, poly[:, 0].max())
                polys.append(poly)
    return polys


def spiralWG(x0=0,y0=0,cellName='',wgW=0.8,bendR=50,turns=12,wgSpace=0.15,layer=0):
	block_cell = gdspy.Cell(cellName+str(uuid.uuid4()) )
	x=x0
	y=y0+(bendR*(turns+1)*wgSpace)

	wg1 = gdspy.Path(wgW, (x,y))

	wg1.turn(bendR, angle=np.pi,layer=layer,number_of_points=500)	# taper
	wg1.turn(bendR, angle=-np.pi,layer=layer,number_of_points=500)	# taper

	for ii in range(1,turns+1):
		wg1.turn(bendR*(2.+wgSpace*(ii)), angle=-np.pi,layer=layer,number_of_points=500)	# taper


	wg2 = gdspy.Path(wgW, (x,y))
	wg2 = wg2.segment(0,direction='-x',layer=0)
	wg2.turn(bendR*(2.0+wgSpace/2), angle=-np.pi,layer=layer,number_of_points=500)	# taper
	for ii in range(0,turns-1):
		wg2.turn(bendR*((2.+wgSpace*2)+wgSpace*ii), angle=-np.pi,layer=layer,number_of_points=500)	# taper

	block_cell.add(wg1)
	block_cell.add(wg2)
	
	spiralLength=wg1.length+wg2.length
	print("Spiral Length : "+str(spiralLength ))

	return [block_cell, [wg1.x,wg1.y],spiralLength]

def spiralDeviceWG(x0=0,y0=0,cellName='',wgW=0.8,bendR=50,turns=12,wgSpace=0.15,layer=0,taperW=0.8,taperL=100):
	block_cell = gdspy.Cell(cellName+str(uuid.uuid4()))
	path1 = gdspy.Path(taperW, (x0, y0))
	# block_cell.add(gdspy.Text(num, 20, (path1.x+200,path1.y+10), layer=layer))
	path1.segment(taperL,direction='+x', final_width=wgW,layer=layer) # taper
	path1.segment(chipW/2.5,direction='+x', final_width=wgW,layer=layer) # taper
	[wgSpiralCell,[x,y],spiralLength] = spiralWG(x0=path1.x,y0=path1.y,cellName="cell_01",wgW=wgW,bendR=bendR,turns=turns,wgSpace=wgSpace,layer=layer )
	
	path2 = gdspy.Path(wgW, (x, y))
	path2.segment(chipW-taperL-x,direction='+x', final_width=wgW,layer=layer) # taper
	path2.segment(taperL,direction='+x', final_width=taperW,layer=layer) # taper
	
	block_cell.add([path1,path2])
	top_cell.add(gdspy.CellReference(wgSpiralCell, (0, 0)))
	return block_cell



def spiralDeviceRing(x0=0,y0=0,cellName='',wgW=0.8,bendR=50,ringPos=1000
	,shiftPos=4000,shiftLen=400,turns=12,wgSpace=0.15,layer=0
	,taperW=0.8,gap=0.8,taperL=100):
	block_cell = gdspy.Cell(cellName+str(uuid.uuid4()))
	path1 = gdspy.Path(taperW, (x0, y0))
	# block_cell.add(gdspy.Text(num, 20, (path1.x+200,path1.y+10), layer=layer))
	path1.segment(taperL,direction='+x', final_width=wgW,layer=layer) # taper
	path1.segment(shiftPos,direction='+x', final_width=wgW,layer=layer) # taper

	# xx=path1.x-1000
	xx=ringPos
	yy=path1.y+gap+wgW

	[wgSpiralCell,[x,y],spiralLength] = spiralWG(x0=xx,y0=yy,cellName="cell_01"
		,wgW=wgW,bendR=bendR,turns=turns,wgSpace=wgSpace,layer=layer )
	path3 = gdspy.Path(wgW, (xx, yy	))
	path3.segment(0,direction='-x',layer=0)

	path3.turn( bendR*(2.+wgSpace*(turns+3))   , angle=-np.pi,layer=layer,number_of_points=500)	# taper
	path3.turn(bendR, angle=-np.pi,layer=layer,number_of_points=500)	# taper



	path1.turn(bendR,angle=np.pi/2,layer=layer,number_of_points=100)		# bend
	path1.segment(shiftLen,direction='+y',layer=layer) 	#wg after bendR
	path1.turn(bendR,angle=-np.pi/2,layer=layer,number_of_points=100)		# bend
	path1.segment(chipW-taperL-path1.x,direction='+x', final_width=wgW,layer=layer) # taper
	path1.segment(taperL,direction='+x', final_width=taperW,layer=layer) # taper
	

	block_cell.add([path1,path3])
	block_cell.add(gdspy.CellReference(wgSpiralCell, (0, 0)))
	textStr = 'L=' +str(np.round(spiralLength)) + ' ,w=' + str(round(wgW,2)) + ' ,g=' + str(round(gap,2))	
	# textStr=str(123213)+'4324'
	text = gdspy.PolygonSet(render_text(textStr, 20,position=(path1.x-600,path1.y+15), font_prop=fp), layer=layer)
	top_cell.add(text)  
	return block_cell





top_cell = gdspy.Cell("topcell")
fp = FontProperties(family='serif', style='normal')

for ii in range (0,4):
	wgSpir = spiralDeviceWG(0,ii*750,"cell_01",wgSpace=0.4)
	top_cell.add(gdspy.CellReference(wgSpir, (0, 0)))


# for ii in range (0,5):
for ii in range(0,5):
	for kk in range(0,5):
		gap=0.7+0.05*kk
		wgSpir = spiralDeviceRing(0,5000+ii*1200+kk*40
			,ringPos=chipW-(1200+kk*1000),shiftPos=chipW-(800+kk*1000)
			,shiftLen=800,cellName="cell_01",wgSpace=0.4, turns=8+ii,gap=gap)
		top_cell.add(gdspy.CellReference(wgSpir, (0, 0)))







## ------------------------------------------------------------------ ##
##  OUTPUT                                                            ##
## ------------------------------------------------------------------ ##

## Output the layout to a GDSII file (default to all created cells).
gdspy.write_gds(exportname+'.gds', unit=1.0e-6, precision=1.0e-9)


faltten_cell=top_cell.flatten()
gdspy.write_gds(exportname+'_flattened.gds', cells=[faltten_cell], unit=1.0e-6, precision=1.0e-9)
# print('Using gdspy module version ' + gdspy.__version__)


#Cell area
textarea = 'Area = %d um^2' % (top_cell.area()) 
exptime=(top_cell.area()*0.00008961209);
exptime1=(top_cell.area()*0.000256);
print(textarea + ' ZEP exposure time :'+ str(round(exptime)) +' min, HSQ exposure time:'+str(round(exptime1))+' min')



