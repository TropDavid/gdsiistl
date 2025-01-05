# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 12:52:15 2024

@author: user
"""
import gdspy
import numpy as np
import uuid


ld_X1 = {"layer": 2,    "datatype": 0} # 800 Nitride
ld_X3 = {"layer":70 , "datatype":0} # 350 Nitride
ld_LNAP = {"layer": 301,    "datatype": 0} # LN Path
ld_LNAX = {"layer":302 , "datatype": 4} # LN Not 
ld_LNASIZE = {"layer": 303 , "datatype": 0} # mark of where LN is
ld_LNARFVIA = {"layer": 320 , "datatype": 0} # Conection between LN and Electrode
ld_LNARFP = {"layer": 330 , "datatype": 0} # Electrode layer
ld_LNARFPAD = {"layer": 331 , "datatype": 0} # Opening to the Electrode


lib = gdspy.GdsLibrary()
cell =lib.new_cell('TOP')

lib.read_gds("‏‏My_Cells_2_use.GDS")

# the notation is from left to right and the (0,0) is at the (left, middle) of the BB (So 5 overlap in required) , length of 400
Taper_end=lib.cells['ligentecInvertedTaper_w1.0BB'] 

# the notation is from left to right and the (0,0) is at the (left, middle) of the BB (So 5 overlap in required)
# width 30 , length 220 , the 1 is at center and the 2 are at +- 12.27
MMI_1X2=lib.cells['ligentecMMI1x2BB']


# the notation is from left to right and the (0,0) is at the (left, middle) of the BB (So 5 overlap in required)
# the Top and Bottom Metal are 50 width, and middle is 30 , the 2 WG are at +-19 , The Top start at 23 and the Bottom at -23
# the overall width of the BB is 156 and length is 10,752 
GSGM=lib.cells['ligentecLNA15ModulatorPushPullCbandLongBB'] 

############################################################################################################################################
#  50 micron from CHS layer if less then 100 from CSL and 10 if its more
#  X1 width 0.2 , dis 0.3
#  X3 width 0.2 , dis 0.3
#  LNAP width 2 , dis 5
#  LNAP must be enclosed by LNASIZE
#  LNARFP width 1.5, dis 3
#  LNARFP must be enclosed by LNA by 5 , lLNA = (LNASIZE NOT LNAX) OR LNAP
#  LNARFVIA width 3, dis 3
#  LNARFVIA  must be enclosed by LNARFP by 5
#  LNARFPAD width 10, dis 10
#  LNARFPAD  must be enclosed by LNARFP by 10
#  Dis between LNARFPAD and LNARFVIA is 10
############################################################################################################################################

Cell_Length = 15810
Cell_Width = 4850

WG_Width = 1
G_RF_Width = 50
S_RF_Width = 30
Taper_Length = 400
Overlap_Length = 5
MMI_Length = 220
radius_bend = 100


def sbendPath(wgsbend,L=100,H=50,info = ld_X1):
# the formula for cosine-shaped s-bend is: y(x) = H/2 * [1- cos(xpi/L)]
# the formula for sine-shaped s-bend is: y(x) = xH/L - H/(2pi) * sin(x2*pi/L)
    def sbend(t):
        y = H/2 * (1- np.cos(t*np.pi))
        x =L*t
        
        return (x,y)
    
    def dtsbend(t):
        dy_dt = H/2*np.pi*np.sin(t*np.pi)
        dx_dt = L

        return (dx_dt,dy_dt)

    wgsbend.parametric(sbend ,dtsbend , number_of_evaluations=100,**info)  
    return wgsbend   
 

def sbendPathM(wgsbend,L=100,H=50,info = ld_X1):

    def sbend(t):
        y = H/2 * (np.cos(t*np.pi))
        x = L*t
        
        return (x,y)
    
    def dtsbend(t):
        dy_dt =  -H/2*np.pi*np.sin(t*np.pi)
        dx_dt = L

        return (dx_dt,dy_dt )

    wgsbend.parametric(sbend ,dtsbend , number_of_evaluations=100,**info)  
    return wgsbend      


def sbendPathMBetter(wgsbend,L=100,H=50,info = ld_X1):

    def sbend(t):
        y = -(H/2 * (1- np.cos(t*np.pi)))
        x = L*t
        
        return (x,y)
    
    def dtsbend(t):
        dy_dt =  -H/2*np.pi*np.sin(t*np.pi)
        dx_dt = L

        return (dx_dt,dy_dt )

    wgsbend.parametric(sbend ,dtsbend , number_of_evaluations=100,**info)  
    return wgsbend      

def a2r(ang):  # angle to radian
    return np.pi/180*ang



cell.add(gdspy.CellReference(Taper_end, (400,0) , rotation = 180))
path1 = gdspy.Path(width = WG_Width , initial_point = (Taper_Length - Overlap_Length ,0))
path1.segment(length = 1000 , direction = "+x" , **ld_X1)

cell.add(gdspy.CellReference(MMI_1X2, (path1.x,path1.y)))

pathTop = gdspy.Path(width = WG_Width , initial_point = (path1.x + MMI_Length - 10, 12.27))
pathBottom = gdspy.Path(width = WG_Width , initial_point = (path1.x + MMI_Length - 10, -12.27))

pathTop.segment(length = 100 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 100 , direction = "+x" , **ld_X1)

pathTop = sbendPath(pathTop , L = 200 , H = 19 - 12.27)
pathBottom = sbendPathMBetter(pathBottom , L = 200 , H = 19 - 12.27)

pathTop.segment(length = 100 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 100 , direction = "+x" , **ld_X1)

##############################################################################################################
Left_Electrode_x = pathTop.x 
Right_Electrode_x = pathTop.x + 10752 
G_Top_y = 23 + 25
G_Bottom_y = -23 - 25
S_y = 0
##############################################################################################################

cell.add(gdspy.CellReference(GSGM, (pathTop.x-Overlap_Length ,path1.y)))


pathTop.x = pathTop.x + 10752 - Overlap_Length 
pathBottom.x = pathBottom.x + 10752 - Overlap_Length 

pathTop.segment(length = 100 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 100 , direction = "+x" , **ld_X1)

pathTop = sbendPathMBetter(pathTop , L = 200 , H = 19 - 12.27)
pathBottom = sbendPath(pathBottom , L = 200 , H = 19 - 12.27)

pathTop.turn(radius = radius_bend , angle = 'l' , **ld_X1)
pathTop.turn(radius = radius_bend , angle = -np.pi , **ld_X1)
pathTop.turn(radius = radius_bend , angle = 'l' , **ld_X1)

pathBottom.turn(radius = radius_bend , angle = 'r' , **ld_X1)
pathBottom.segment(length = 100 , direction = "-y" , **ld_X1)
pathBottom.turn(radius = radius_bend , angle = np.pi , **ld_X1)
pathBottom.segment(length = 100 , direction = "+y" , **ld_X1)
pathBottom.turn(radius = radius_bend , angle = 'r' , **ld_X1)

pathTop.segment(length = 100 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 100 , direction = "+x" , **ld_X1)

cell.add(gdspy.CellReference(MMI_1X2, (pathTop.x + 220  - 10 ,path1.y) , rotation = 180))

path1.x = pathTop.x + MMI_Length  - 2* Overlap_Length

path1.segment(length = Cell_Length - path1.x - 395 , direction = "+x" , **ld_X1)

cell.add(gdspy.CellReference(Taper_end, (path1.x - Overlap_Length ,path1.y)))


cell.add(path1)
cell.add(pathTop)
cell.add(pathBottom)

#########################################################################################################################################################################################
FlexPath = gdspy.FlexPath(points = [(Left_Electrode_x + Overlap_Length ,S_y),(Left_Electrode_x - Overlap_Length*4 , S_y)], width = [50, 30 , 50] , offset = 48 , **ld_LNARFP)
FlexPath.arc(radius = 100, initial_angle = a2r(-90), final_angle = a2r(-180))
FlexPath.segment(end_point = (0,500) ,width = (80,80,80) , offset = (-150,0,150) , relative = True)
cell.add(gdspy.Rectangle((Left_Electrode_x - Overlap_Length*4 - 100  - 50, S_y + 100 + 480), (Left_Electrode_x - Overlap_Length*4 - 100  + 50, S_y + 100 + 580),**ld_LNARFP))
cell.add(gdspy.Rectangle((Left_Electrode_x - Overlap_Length*4 - 100 - 150  - 50, S_y + 100 + 480), (Left_Electrode_x - Overlap_Length*4 - 100  - 150 + 50, S_y + 100 + 580),**ld_LNARFP))
cell.add(gdspy.Rectangle((Left_Electrode_x - Overlap_Length*4 - 100 +150 - 50, S_y + 100 + 480), (Left_Electrode_x - Overlap_Length*4 - 100 + 150 + 50, S_y + 100 + 580),**ld_LNARFP))

cell.add(gdspy.Rectangle((Left_Electrode_x - Overlap_Length*4 - 100  - 50, S_y + 100 + 480), (Left_Electrode_x - Overlap_Length*4 - 100  + 50, S_y + 100 + 580),**ld_LNARFPAD))
cell.add(gdspy.Rectangle((Left_Electrode_x - Overlap_Length*4 - 100 - 150  - 50, S_y + 100 + 480), (Left_Electrode_x - Overlap_Length*4 - 100  - 150 + 50, S_y + 100 + 580),**ld_LNARFPAD))
cell.add(gdspy.Rectangle((Left_Electrode_x - Overlap_Length*4 - 100 +150 - 50, S_y + 100 + 480), (Left_Electrode_x - Overlap_Length*4 - 100 + 150 + 50, S_y + 100 + 580),**ld_LNARFPAD))




FlexPath1 = gdspy.FlexPath(points = [(Right_Electrode_x - Overlap_Length ,S_y),(Right_Electrode_x + Overlap_Length*4 , S_y)], width = [50, 30 , 50] , offset = 48 , **ld_LNARFP)
FlexPath1.arc(radius = 100, initial_angle = a2r(90), final_angle = a2r(0))
FlexPath1.segment(end_point = (0,-500) ,width = (80,80,80) , offset = (-150,0,150) , relative = True)
cell.add(gdspy.Rectangle((Right_Electrode_x + Overlap_Length*4 + 100  - 50, S_y - 100 - 480), (Right_Electrode_x + Overlap_Length*4 + 100  + 50, S_y - 100 - 580),**ld_LNARFP))
cell.add(gdspy.Rectangle((Right_Electrode_x + Overlap_Length*4 + 100 - 150  - 50, S_y - 100 - 480), (Right_Electrode_x + Overlap_Length*4 + 100  - 150 + 50, S_y - 100 - 580),**ld_LNARFP))
cell.add(gdspy.Rectangle((Right_Electrode_x + Overlap_Length*4 + 100 + 150 - 50, S_y - 100 - 480), (Right_Electrode_x + Overlap_Length*4 + 100 + 150 + 50, S_y - 100 - 580),**ld_LNARFP))

cell.add(gdspy.Rectangle((Right_Electrode_x + Overlap_Length*4 + 100  - 50, S_y - 100 - 480), (Right_Electrode_x + Overlap_Length*4 + 100  + 50, S_y - 100 - 580),**ld_LNARFPAD))
cell.add(gdspy.Rectangle((Right_Electrode_x + Overlap_Length*4 + 100 - 150  - 50, S_y - 100 - 480), (Right_Electrode_x + Overlap_Length*4 + 100  - 150 + 50, S_y - 100 - 580),**ld_LNARFPAD))
cell.add(gdspy.Rectangle((Right_Electrode_x + Overlap_Length*4 + 100 + 150 - 50, S_y - 100 - 480), (Right_Electrode_x + Overlap_Length*4 + 100 + 150 + 50, S_y - 100 - 580),**ld_LNARFPAD))




cell.add(FlexPath)
cell.add(FlexPath1)
###########################################################################################################################################################################################


lib.write_gds('Modulator.gds')   





