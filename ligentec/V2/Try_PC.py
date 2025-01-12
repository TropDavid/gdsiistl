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

# MMI 2X2 width 1 at (-5,23.74) , (-5,0) OUT at (365, 23.74) , (365,0)
# length 370 width 60
MMI2X2 = lib.cells['ligentecMMI2x2BB']

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
MMI_Length = 365


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

pathTop  = sbendPath(pathTop , L = 250 , H = (230 - 54.6)/2)
pathBottom = sbendPathMBetter(pathBottom , L = 250 , H = (230 - 54.6)/2)

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

####################################################################PM1##############################################
cell.add(gdspy.CellReference(PM , (pathTop.x , pathTop.y)))
cell.add(gdspy.CellReference(PM , (pathBottom.x , pathBottom.y)))

pathTop.x = pathTop.x + 3552 - 10
pathBottom.x = pathBottom.x + 3552 - 10

x_C1 = pathTop.x 
y_C1_T_Upper = pathTop.y + 4 + 45
y_C1_T_Lowwer = pathTop.y - 4 - 45
y_C1_B_Upper = pathBottom.y + 4 + 45
y_C1_B_Lowwer = pathBottom.y - 4 - 45


pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

pathTop = sbendPathMBetter(pathTop,L = 250 , H = pathTop.y - 23.74)
pathBottom = sbendPath(pathBottom, L = 250 , H = abs(pathBottom.y))

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

################################################################MMI1##################################################
cell.add(gdspy.CellReference(MMI2X2,((pathBottom.x,pathBottom.y))))

pathTop.x = pathTop.x + MMI_Length - 10
pathBottom.x = pathBottom.x + MMI_Length - 10

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

pathTop  = sbendPath(pathTop , L = 250 , H = (230 - 10)/2)
pathBottom = sbendPathMBetter(pathBottom , L = 250 , H = (230 - 10)/2)

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

#############################################################PM2######################################################
cell.add(gdspy.CellReference(PM , (pathTop.x , pathTop.y)))
cell.add(gdspy.CellReference(PM , (pathBottom.x , pathBottom.y)))

x_C2 = pathTop.x 
y_C2_T_Upper = pathTop.y + 4 + 45
y_C2_T_Lowwer = pathTop.y - 4 - 45
y_C2_B_Upper = pathBottom.y + 4 + 45
y_C2_B_Lowwer = pathBottom.y - 4 - 45

#######################################################################################################################
C_x = pathTop.x + 3552/2
C_y = pathTop.y + 210/2 + 200
#######################################################################################################################

pathTop.x = pathTop.x + 3552 - 10
pathBottom.x = pathBottom.x + 3552 - 10


pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

pathTop = sbendPathMBetter(pathTop,L = 250 , H = pathTop.y - 23.74)
pathBottom = sbendPath(pathBottom, L = 250 , H = abs(pathBottom.y))

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

####################################################################MMI2###############################################
cell.add(gdspy.CellReference(MMI2X2,((pathBottom.x,pathBottom.y))))

pathTop.x = pathTop.x + MMI_Length - 10
pathBottom.x = pathBottom.x + MMI_Length - 10

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

pathTop  = sbendPath(pathTop , L = 250 , H = (230 - 10)/2)
pathBottom = sbendPathMBetter(pathBottom , L = 250 , H = (230 - 10)/2)

pathTop.segment(length = 50 , direction = "+x" , **ld_X1)
pathBottom.segment(length = 50 , direction = "+x" , **ld_X1)

#############################################################PM3######################################################
cell.add(gdspy.CellReference(PM , (pathTop.x , pathTop.y)))
cell.add(gdspy.CellReference(PM , (pathBottom.x , pathBottom.y)))

x_C3 = pathTop.x 
y_C3_T_Upper = pathTop.y + 4 + 45
y_C3_T_Lowwer = pathTop.y - 4 - 45
y_C3_B_Upper = pathBottom.y + 4 + 45
y_C3_B_Lowwer = pathBottom.y - 4 - 45


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
############################################################################Contact##########################################################################
C_x = C_x - 4*150
for i in range(9):
     cell.add(gdspy.Rectangle((C_x+i*150 - 50,C_y - 50), (C_x+i*150 + 50 ,C_y + 50) , **ld_LNARFP))
     cell.add(gdspy.Rectangle((C_x+i*150 - 50,C_y - 50), (C_x+i*150 + 50 ,C_y + 50) , **ld_LNARFPAD))
     
############################################################################C1##############################################################################     
pathC1_T_Upper = gdspy.Path(width = 90 , initial_point = (x_C1 - 5 , y_C1_T_Upper))
pathC1_T_Upper.segment(length = 10 , direction = "+x" , **ld_LNARFP)
pathC1_T_Upper.segment(length = 200 , direction = "+x" ,final_width = 5 , **ld_LNARFP)
pathC1_T_Upper = sbendPath(pathC1_T_Upper , L = 300, H = C_y - pathC1_T_Upper.y , info = ld_LNARFP)
pathC1_T_Upper.segment(length = C_x - pathC1_T_Upper.x , direction = '+x' , **ld_LNARFP)
cell.add(pathC1_T_Upper)


pathC1_T_Lowwer = gdspy.Path(width = 90 , initial_point = (x_C1 -5, y_C1_T_Lowwer))
pathC1_T_Lowwer.segment(length = 100 , direction = "+x" , **ld_LNARFP)

pathC1_B_Upper = gdspy.Path(width = 90 , initial_point = (x_C1 -5, y_C1_B_Upper))
pathC1_B_Upper.segment(length = 100 , direction = "+x" , **ld_LNARFP)

cell.add([pathC1_T_Lowwer, pathC1_B_Upper])
cell.add(gdspy.Rectangle((pathC1_B_Upper.x - 10 , pathC1_B_Upper.y - 90/2),(pathC1_T_Lowwer.x + 100 , pathC1_T_Lowwer.y + 90/2),**ld_LNARFP))
pathC1_Ground = gdspy.Path(width = 5 , initial_point = (pathC1_B_Upper.x + 95 , (pathC1_B_Upper.y + pathC1_T_Lowwer.y)/2))
pathC1_Ground = sbendPath(pathC1_Ground , L = 500 - 95 - 100 , H = C_y - 50 - 10 - pathC1_Ground.y , info = ld_LNARFP)
pathC1_Ground.segment(length = C_x + 1*150 - pathC1_Ground.x  - 10, direction = '+x' , **ld_LNARFP)
pathC1_Ground.turn(radius = 10 ,angle =  'l' ,tolerance = 0.00001 , **ld_LNARFP )
cell.add(pathC1_Ground)


pathC1_B_Lowwer = gdspy.Path(width = 90 , initial_point = (x_C1 - 5 , y_C1_B_Lowwer))
pathC1_B_Lowwer.segment(length = 10 , direction = "+x" , **ld_LNARFP)
pathC1_B_Lowwer.segment(length = 200 , direction = "+x" ,final_width = 5 , **ld_LNARFP)
pathC1_B_Lowwer = sbendPath(pathC1_B_Lowwer , L = 500, H = C_y - 50 - 20 - pathC1_B_Lowwer.y , info = ld_LNARFP)
pathC1_B_Lowwer.segment(length = C_x + 2*150 - pathC1_B_Lowwer.x  - 20, direction = '+x' , **ld_LNARFP)
pathC1_B_Lowwer.turn(radius = 20 ,angle =  'l' ,tolerance = 0.00001 , **ld_LNARFP )
cell.add(pathC1_B_Lowwer)

#############################################################################C2#####################################################################################
pathC2_B_Lowwer = gdspy.Path(width = 90 , initial_point = (x_C2 , y_C2_B_Lowwer))
pathC2_B_Lowwer.segment(length = 10 , direction = "-x" , **ld_LNARFP)
pathC2_B_Lowwer.segment(length = 200 , direction = "-x" ,final_width = 5 , **ld_LNARFP)
pathC2_B_Lowwer.arc(radius = (C_y - 50 - 30 - pathC2_B_Lowwer.y)/2 , initial_angle = a2r(-90) , final_angle = a2r(-270) , tolerance = 0.0001 , **ld_LNARFP)
pathC2_B_Lowwer.segment(length = C_x + 3*150 - pathC2_B_Lowwer.x  - 30, direction = '+x' , **ld_LNARFP)
pathC2_B_Lowwer.turn(radius = 30 ,angle =  'l' ,tolerance = 0.00001 , **ld_LNARFP )
cell.add(pathC2_B_Lowwer)


pathC2_T_Lowwer = gdspy.Path(width = 90 , initial_point = (x_C2 , y_C2_T_Lowwer))
pathC2_T_Lowwer.segment(length = 100 , direction = "-x" , **ld_LNARFP)

pathC2_B_Upper = gdspy.Path(width = 90 , initial_point = (x_C2 , y_C2_B_Upper))
pathC2_B_Upper.segment(length = 100 , direction = "-x" , **ld_LNARFP)

cell.add([pathC2_T_Lowwer, pathC2_B_Upper])
cell.add(gdspy.Rectangle((pathC2_B_Upper.x + 10 , pathC2_B_Upper.y - 90/2),(pathC2_T_Lowwer.x - 100 , pathC2_T_Lowwer.y + 90/2),**ld_LNARFP))
pathC2_Ground = gdspy.Path(width = 5 , initial_point = (pathC2_B_Upper.x - 95 , (pathC2_B_Upper.y + pathC2_T_Lowwer.y)/2))
pathC2_Ground.arc(radius = (C_y - 50 - 40 - pathC2_Ground.y)/2 , initial_angle = a2r(-90) , final_angle = a2r(-270) , tolerance = 0.0001 , **ld_LNARFP)
pathC2_Ground.segment(length = C_x + 4*150 - pathC2_Ground.x  - 40, direction = '+x' , **ld_LNARFP)
pathC2_Ground.turn(radius = 40 ,angle =  'l' ,tolerance = 0.00001 , **ld_LNARFP )
cell.add(pathC2_Ground)


pathC2_T_Upper = gdspy.Path(width = 90 , initial_point = (x_C2 , y_C2_T_Upper))
pathC2_T_Upper.segment(length = 10 , direction = "-x" , **ld_LNARFP)
pathC2_T_Upper.segment(length = 200 , direction = "-x" ,final_width = 5 , **ld_LNARFP)
pathC2_T_Upper.arc(radius = (C_y - 50 - 50 - pathC2_T_Upper.y)/2 , initial_angle = a2r(-90) , final_angle = a2r(-270) , tolerance = 0.0001 , **ld_LNARFP)
pathC2_T_Upper.segment(length = C_x + 5*150 - pathC2_T_Upper.x  - 50, direction = '+x' , **ld_LNARFP)
pathC2_T_Upper.turn(radius = 50 ,angle =  'l' ,tolerance = 0.00001 , **ld_LNARFP )
cell.add(pathC2_T_Upper)

#############################################################################C3#####################################################################################
pathC3_T_Upper = gdspy.Path(width = 90 , initial_point = (x_C3 , y_C3_T_Upper))
pathC3_T_Upper.segment(length = 10 , direction = "-x" , **ld_LNARFP)
pathC3_T_Upper.segment(length = 200 , direction = "-x" ,final_width = 5 , **ld_LNARFP)

pathC3_T_Upper_mid = gdspy.Path(width = 5 , initial_point = (pathC3_T_Upper.x , pathC3_T_Upper.y))
pathC3_T_Upper_mid = sbendPath(pathC3_T_Upper_mid , L = 300, H = C_y - pathC3_T_Upper_mid.y , info = ld_LNARFP)
pathC3_T_Upper_mid.mirror((pathC3_T_Upper.x,pathC3_T_Upper.y+5),(pathC3_T_Upper.x,pathC3_T_Upper.y-5))
pathC3_T_Upper_mid.segment(length = pathC3_T_Upper_mid.x -C_x - 8*150 , direction = '-x' , **ld_LNARFP)
cell.add(pathC3_T_Upper)
cell.add(pathC3_T_Upper_mid)


pathC3_T_Lowwer = gdspy.Path(width = 90 , initial_point = (x_C3, y_C3_T_Lowwer))
pathC3_T_Lowwer.segment(length = 100 , direction = "-x" , **ld_LNARFP)

pathC3_B_Upper = gdspy.Path(width = 90 , initial_point = (x_C3, y_C3_B_Upper))
pathC3_B_Upper.segment(length = 100 , direction = "-x" , **ld_LNARFP)

cell.add([pathC3_T_Lowwer, pathC3_B_Upper])

cell.add(gdspy.Rectangle((pathC3_B_Upper.x + 10 , pathC3_B_Upper.y - 90/2),(pathC3_T_Lowwer.x - 100 , pathC3_T_Lowwer.y + 90/2),**ld_LNARFP))

pathC3_Ground = gdspy.Path(width = 5 , initial_point = (pathC3_T_Lowwer.x - 100 , (pathC3_B_Upper.y + pathC3_T_Lowwer.y)/2))
pathC3_Ground = sbendPath(pathC3_Ground , L = 500 - 95 - 100 , H = C_y - 50 - 10 - pathC3_Ground.y , info = ld_LNARFP)
pathC3_Ground.mirror((pathC3_B_Upper.x - 100,pathC3_B_Upper.y+5),(pathC3_B_Upper.x - 100,pathC3_B_Upper.y-5))
pathC3_Ground.segment(length =  pathC3_Ground.x  - 10 - C_x - 7*150  , direction = '-x' , **ld_LNARFP)
pathC3_Ground.arc(radius = 10 , initial_angle = a2r(-90) , final_angle = a2r(-180) , tolerance = 0.0001 , **ld_LNARFP)
cell.add(pathC3_Ground)



pathC3_B_Lowwer = gdspy.Path(width = 90 , initial_point = (x_C3 , y_C3_B_Lowwer))
pathC3_B_Lowwer.segment(length = 10 , direction = "-x" , **ld_LNARFP)
pathC3_B_Lowwer.segment(length = 200 , direction = "-x" ,final_width = 5 , **ld_LNARFP)

pathC3_B_Lowwer_mid = gdspy.Path(width = 5 , initial_point = (pathC3_B_Lowwer.x , pathC3_B_Lowwer.y))
pathC3_B_Lowwer_mid = sbendPath(pathC3_B_Lowwer_mid , L = 500, H = C_y - pathC3_B_Lowwer_mid.y - 50 - 20 , info = ld_LNARFP)
pathC3_B_Lowwer_mid.mirror((pathC3_B_Lowwer.x,pathC3_B_Lowwer.y+5),(pathC3_B_Lowwer.x,pathC3_B_Lowwer.y-5))
pathC3_B_Lowwer_mid.segment(length = pathC3_B_Lowwer_mid.x -C_x - 6*150 - 20 , direction = '-x' , **ld_LNARFP)
pathC3_B_Lowwer_mid.arc(radius = 20 , initial_angle = a2r(-90) , final_angle = a2r(-180) , tolerance = 0.0001 , **ld_LNARFP)
cell.add(pathC3_B_Lowwer)
cell.add(pathC3_B_Lowwer_mid)







lib.write_gds('Try_PC.gds') 











