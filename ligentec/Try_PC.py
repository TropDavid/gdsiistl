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

# the notation is from left to right and the (-5,0) is at the
# width 30 , length 220 , the 1 is at center and the 2 are at +- 12.27
MMI_1X2=lib.cells['ligentecMMI1x2BB']


# the notation is from left to right and the (0,0) is at the (left, middle) of the BB (So 5 overlap in required)
# the Top and Bottom Metal are 50 width, and middle is 30 , the 2 WG are at +-19 , The Top start at 23 and the Bottom at -23
# the overall width of the BB is 156 and length is 10,752 
GSGM=lib.cells['ligentecLNA15ModulatorPushPullCbandLongBB'] 


# WG at (-5,0) width 1, Contact are of width 90 at (-5,+-49)
# The width is 210 and length 3552
PM = lib.cells['ligentecLNA15PhaseShifterCbandShortBB']


# WG input (-5,0) width 1 , wg Output at (515,0) and (515,-54.6)
# length 520 and width 70
PBS = lib.cells['ligentecPBSBB']

# DC input at (-5,10) , (-5,0) Otput (220+Length , 10) , (220+Length, 0)
# width 60 and length 220 + Length
DC_L124 = lib.cells['ligentecDC_L124.0BB']

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
C_RF_Width = 90
Taper_Length = 400
Overlap_Length = 5
PBS_Length = 520
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


cell.add(gdspy.CellReference(Taper_end, (400,27.3) , rotation = 180))
path1 = gdspy.Path(width = WG_Width , initial_point = (Taper_Length - Overlap_Length ,27.3))
path1.segment(length = 50 , direction = "+x" , **ld_X1)

###########################################PBS###################################################################
cell.add(gdspy.CellReference(PBS, (path1.x,path1.y)))

pathTop = gdspy.Path(width = WG_Width , initial_point = (path1.x + PBS_Length - 10, path1.y))
pathBottom = gdspy.Path(width = WG_Width , initial_point = (path1.x + PBS_Length - 10, path1.y - 54.6))

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

pathTop  = sbendPath(pathTop , L = 250 , H = (430 - 54.6)/2)
pathBottom = sbendPathMBetter(pathBottom , L = 250 , H = (430 - 54.6)/2)

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

####################################################################PM1##############################################
cell.add(gdspy.CellReference(PM , (pathTop.x , pathTop.y)))
cell.add(gdspy.CellReference(PM , (pathBottom.x , pathBottom.y)))

pathTop.x = pathTop.x + 3552 - 10
pathBottom.x = pathBottom.x + 3552 - 10

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

pathTop = sbendPathMBetter(pathTop,L = 250 , H = pathTop.y - 10 - 25)
pathBottom = sbendPath(pathBottom, L = 250 , H = abs(pathBottom.y) + 25)

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

################################################################DC1##################################################
cell.add(gdspy.CellReference(DC_L124,((pathBottom.x,pathBottom.y))))

pathTop.x = pathTop.x + 220+ 124 - 10
pathBottom.x = pathBottom.x + 220 + 124 - 10

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

pathTop  = sbendPath(pathTop , L = 250 , H = (430 - 10)/2)
pathBottom = sbendPathMBetter(pathBottom , L = 250 , H = (430 - 10)/2)

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

#############################################################PM2######################################################
cell.add(gdspy.CellReference(PM , (pathTop.x , pathTop.y)))
cell.add(gdspy.CellReference(PM , (pathBottom.x , pathBottom.y)))

#######################################################################################################################
C_x = pathTop.x + 3552/2
C_y = pathTop.y + 210/2 + 150
#######################################################################################################################

pathTop.x = pathTop.x + 3552 - 10
pathBottom.x = pathBottom.x + 3552 - 10



pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

pathTop = sbendPathMBetter(pathTop,L = 250 , H = pathTop.y - 10 - 25)
pathBottom = sbendPath(pathBottom, L = 250 , H = abs(pathBottom.y) + 25)

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

####################################################################DC2###############################################
cell.add(gdspy.CellReference(DC_L124,((pathBottom.x,pathBottom.y))))

pathTop.x = pathTop.x + 220+ 124 - 10
pathBottom.x = pathBottom.x + 220 + 124 - 10

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

pathTop  = sbendPath(pathTop , L = 250 , H = (430 - 10)/2)
pathBottom = sbendPathMBetter(pathBottom , L = 250 , H = (430 - 10)/2)

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

#############################################################PM3######################################################
cell.add(gdspy.CellReference(PM , (pathTop.x , pathTop.y)))
cell.add(gdspy.CellReference(PM , (pathBottom.x , pathBottom.y)))

pathTop.x = pathTop.x + 3552 - 10
pathBottom.x = pathBottom.x + 3552 - 10

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

pathTop = sbendPathMBetter(pathTop,L = 250 , H = pathTop.y - 27.3)
pathBottom = sbendPath(pathBottom, L = 250 , H = abs(pathBottom.y) - 27.3)

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

cell.add(gdspy.CellReference(PBS, (pathBottom.x + PBS_Length - Overlap_Length*2 ,pathBottom.y),rotation = 180 ))

pathBottom.x = pathBottom.x + PBS_Length - 10
pathBottom.segment(length = Cell_Length - pathBottom.x - 400 , direction = "+x" , **ld_X1)
cell.add(gdspy.CellReference(Taper_end, (pathBottom.x,pathBottom.y)))

cell.add(path1)
cell.add(pathTop)
cell.add(pathBottom)
############################################################################Contact#####################################
C_x = C_x - 4*150
for i in range(9):
     cell.add(gdspy.Rectangle((C_x+i*150 - 50,C_y - 50), (C_x+i*150 + 50 ,C_y + 50) , **ld_LNARFP))
     cell.add(gdspy.Rectangle((C_x+i*150 - 50,C_y - 50), (C_x+i*150 + 50 ,C_y + 50) , **ld_LNARFPAD))







lib.write_gds('Try_PC.gds') 











