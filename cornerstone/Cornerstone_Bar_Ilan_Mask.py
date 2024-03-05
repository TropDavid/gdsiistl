# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 12:54:31 2024

@author: user
"""
import gdspy
import numpy as np
import uuid


ld_WG = {"layer": 3 , "datatype": 0}
ld_GC = {"layer":6 , "datatype":0}

chip_L = 6000
chip_W = 3000

WG = 0.45
WG_Broad = 1.5
Broad_Length = 20
period = 0.63
FF = 0.5
MMI_L = 32.7
MMI_W = 6
MMI_diss_From_center = 1.57
Radius = 25
Gap =np.array([0.2,0.225,0.25])
diss_in_bundel = 127
Deg = 50

lib = gdspy.GdsLibrary()
cell = lib.new_cell('TOP')



def grating(period, number_of_teeth, fill_frac, width, position, direction, lda=1, sin_theta=0,
            focus_distance=-1, focus_width=-1, tolerance=0.001, layer=0, datatype=0):
    '''
    Straight or focusing grating.

    period          : grating period
    number_of_teeth : number of teeth in the grating
    fill_frac       : filling fraction of the teeth (w.r.t. the period)
    width           : width of the grating
    position        : grating position (feed point)
    direction       : one of {'+x', '-x', '+y', '-y'}
    lda             : free-space wavelength
    sin_theta       : sine of incidence angle
    focus_distance  : focus distance (negative for straight grating)
    focus_width     : if non-negative, the focusing area is included in
                      the result (usually for negative resists) and this
                      is the width of the waveguide connecting to the
                      grating
    tolerance       : same as in `path.parametric`
    layer           : GDSII layer number
    datatype        : GDSII datatype number

    Return `PolygonSet`
    '''
    if focus_distance < 0:
        p = gdspy.L1Path((position[0] - 0.5 * width,
                          position[1] + 0.5 * (number_of_teeth - 1 + fill_frac) * period),
                         '+x', period * fill_frac, [width], [], number_of_teeth, period,
                         **ld_WG)
    else:
        neff = lda / float(period) + sin_theta
        qmin = int(focus_distance / float(period) + 0.5)
        p = gdspy.Path(period * fill_frac, position)
        c3 = neff**2 - sin_theta**2
        w = 0.5 * width
        for q in range(qmin, qmin + number_of_teeth):
            c1 = q * lda * sin_theta
            c2 = (q * lda)**2
            p.parametric(lambda t: (width * t - w,
                                    (c1 + neff * np.sqrt(c2 - c3 * (width * t - w)**2)) / c3),
                         tolerance=tolerance, max_points=0, **ld_WG )
            p.x = position[0]
            p.y = position[1]
        sz = p.polygons[0].shape[0] // 2
        if focus_width == 0:
            p.polygons[0] = np.vstack((p.polygons[0][:sz, :], [position]))
        elif focus_width > 0:
            p.polygons[0] = np.vstack((p.polygons[0][:sz, :],
                                          [(position[0] + 0.5 * focus_width, position[1]),
                                           (position[0] - 0.5 * focus_width, position[1])]))
        p.fracture()
        #circle = gdspy.Round((position[0],position[1]+25), 15, tolerance=0.01,**ld_Silox)
        #bound = gdspy.Round((position[0],position[1]+25), 30, tolerance=0.01,**ld_GC )
       
    if direction == '-x':
        return p.rotate(0.5 * np.pi, position)
    elif direction == '+x':
        return p.rotate(-0.5 * np.pi, position)
    elif direction == '-y':
        return p.rotate(np.pi, position)
    else:
        return p


def a2r(ang):  # angle to radian
    return np.pi/180*ang
    

def U_Shape_Rings(Straight_WG_Ring = 100,diss_in_bundel = 127 , position = (0,0) , number_of_rings = 1 , gap = 0.2): 
    x = position[0]
    y = position[1]
    
    G = gap 
    text = gdspy.Text("G = "+str(np.round(G,4)),20,(position[0] - 60, position[1] - 0.5*diss_in_bundel))
    cell.add(text)
    
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=position, direction='-x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(x, y - diss_in_bundel), direction='-x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    path1 = gdspy.Path(width = WG , initial_point = (x,y))
    path_d = gdspy.Path(width = WG , initial_point = (x,y - diss_in_bundel))
    
    path1.segment(length = Straight_WG_Ring , direction = '+x' , **ld_WG)
    path_d.segment(length = Straight_WG_Ring - (1+number_of_rings)*gap - (1+number_of_rings)*WG - 2*Radius*number_of_rings , direction = '+x' , **ld_WG)
    
    path1.arc(radius = Radius , initial_angle = a2r(90) , final_angle = a2r(0) , **ld_WG)
    path1.segment(length = diss_in_bundel - Radius  , direction = '-y' , **ld_WG)
    path1.segment(length = diss_in_bundel/2  , direction = '-y' , **ld_WG)
    
    
    x = path1.x    
    for i in range(0,number_of_rings):
        x_r_b = x - WG - gap - Radius*2*i
        path2 = gdspy.Path(width = WG , initial_point = (x_r_b,path1.y))
        path2.arc(radius = Radius , initial_angle = a2r(0) , final_angle = a2r(360) , tolerance = 0.00001 ,**ld_WG)
        cell.add(path2)
        x = x_r_b  
    
    path1.segment(length = diss_in_bundel/2  , direction = '-y' , **ld_WG)
    path1.segment(length = diss_in_bundel - Radius  , direction = '-y' , **ld_WG)
    path1.arc(radius = Radius , initial_angle = a2r(360) , final_angle = a2r(270) , **ld_WG)
    path1.segment(length = Straight_WG_Ring , direction = '-x' , **ld_WG)       
   
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path1.x,path1.y), direction='-x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    path_d.arc(radius = Radius , initial_angle = a2r(90) , final_angle = a2r(0) , **ld_WG)
    path_d.segment(length = diss_in_bundel - Radius*2 , direction = '-y' , **ld_WG)
    path_d.arc(radius = Radius , initial_angle = a2r(360) , final_angle = a2r(270) , **ld_WG)
    path_d.segment(length = Straight_WG_Ring - (2*gap + 2*WG + 2*Radius)*number_of_rings , direction = '-x' , **ld_WG)
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path_d.x,path_d.y), direction='-x', lda=1.55,
                         sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))  
       
        
    cell.add(path1)
    cell.add(path2)
    cell.add(path_d)

def U_Shape_Disk(Straight_WG_Ring = 100,diss_in_bundel = 127 , position = (0,0) , number_of_rings = 1 , gap = 0.2): 
    x = position[0]
    y = position[1]
    
    G = gap 
    text = gdspy.Text("G = "+str(np.round(G,4)),20,(position[0] - 60, position[1] - 0.5*diss_in_bundel))
    cell.add(text)
    
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=position, direction='-x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(x, y - diss_in_bundel), direction='-x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    path1 = gdspy.Path(width = WG , initial_point = (x,y))
    path_d = gdspy.Path(width = WG , initial_point = (x,y - diss_in_bundel))
    
    path1.segment(length = Straight_WG_Ring , direction = '+x' , **ld_WG)
    path_d.segment(length = Straight_WG_Ring - (1+number_of_rings)*gap - (1+number_of_rings)*WG - 2*Radius*number_of_rings , direction = '+x' , **ld_WG)
    
    path1.arc(radius = Radius , initial_angle = a2r(90) , final_angle = a2r(0) , **ld_WG)
    path1.segment(length = diss_in_bundel - Radius  , direction = '-y' , **ld_WG)
    path1.segment(length = diss_in_bundel/2  , direction = '-y' , **ld_WG)
    
    
    x = path1.x    
    for i in range(0,number_of_rings):
        x_r_b = x - WG - gap - Radius -Radius*i
        cell.add(gdspy.Round(center = (x_r_b,path1.y),radius = Radius ,**ld_WG))
        x = x_r_b  
    
    path1.segment(length = diss_in_bundel/2  , direction = '-y' , **ld_WG)
    path1.segment(length = diss_in_bundel - Radius  , direction = '-y' , **ld_WG)
    path1.arc(radius = Radius , initial_angle = a2r(360) , final_angle = a2r(270) , **ld_WG)
    path1.segment(length = Straight_WG_Ring , direction = '-x' , **ld_WG)       
   
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path1.x,path1.y), direction='-x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    path_d.arc(radius = Radius , initial_angle = a2r(90) , final_angle = a2r(0) , **ld_WG)
    path_d.segment(length = diss_in_bundel - Radius*2 , direction = '-y' , **ld_WG)
    path_d.arc(radius = Radius , initial_angle = a2r(360) , final_angle = a2r(270) , **ld_WG)
    path_d.segment(length = Straight_WG_Ring - (2*gap + 2*WG + 2*Radius)*number_of_rings , direction = '-x' , **ld_WG)
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path_d.x,path_d.y), direction='-x', lda=1.55,
                         sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))  
       
        
    cell.add(path1)
    cell.add(path_d)

def Straight_Rings (WG = WG,Straight_WG = 100 , position = (0,0) , number_of_rings_ontop = 1 , number_of_ring_side = 1 , gap = 0.2):
    
    G = gap 
    text = gdspy.Text("G = "+str(np.round(G,4)),10,(position[0] - 60, position[1] + Radius))
    cell.add(text)
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=position, direction='-x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    path1 = gdspy.Path(width = WG , initial_point = position)
    
    for j in range(0,number_of_ring_side):
        path1.segment(length = Straight_WG/(1+number_of_ring_side) , direction = '+x' , **ld_WG)
        x = path1.x
        y = path1.y
        for i in range(0,number_of_rings_ontop):
            y_r_b = y + WG + gap + 2*Radius*i
            path2 = gdspy.Path(width = WG , initial_point = (x,y_r_b))
            path2.arc(radius = Radius , initial_angle = a2r(270) , final_angle = a2r(630) , tolerance = 0.00001 ,**ld_WG)
            cell.add(path2)
            y = y_r_b 
    path1.segment(length = Straight_WG/(1+number_of_ring_side) , direction = '+x' , **ld_WG)
    y_drop = y_r_b + 2*Radius + gap + WG
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path1.x , y_drop), direction='+x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))    
    
    path_d = gdspy.Path(width = WG , initial_point = (path1.x , y_drop))
    path_d.segment(length = Straight_WG, direction = '-x' , **ld_WG)
    
    cell.add(path1)
    cell.add(path_d)
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path1.x,path1.y), direction='+x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path_d.x,path_d.y), direction='-x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))

def MZI (Straight_in_bend = 20 , MMI_L = 32.7 , MMI_W = 6 , MMI_diss_From_center = 1.57 , Deg = 45 , position = (0,0) , number_of_bend = 1):
    L = int(4*Radius*(np.pi-1)*number_of_bend + Straight_in_bend)
    text = gdspy.Text("L = "+str(np.round(L,0)),20,(position[0] - 60, position[1] + 30))
    cell.add(text)
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=position, direction='-x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))    
    
    path1 = gdspy.Path(width = WG , initial_point = position)
    path1.segment(length = Broad_Length , direction = '+x' , final_width = WG_Broad , **ld_WG)
    cell.add(path1)
    cell.add(gdspy.Rectangle((path1.x , path1.y - MMI_W/2), (path1.x + MMI_L , path1.y + MMI_W/2) , **ld_WG))
    x = path1.x + MMI_L
    y_T = path1.y + MMI_diss_From_center
    y_B = path1.y - MMI_diss_From_center
    
    path_T = gdspy.Path(width = WG_Broad , initial_point = (x,y_T))
    path_B = gdspy.Path(width = WG_Broad , initial_point = (x,y_B))
    
    path_T.segment(length = Broad_Length , direction = '+x' , final_width = WG , **ld_WG)
    path_B.segment(length = Broad_Length , direction = '+x' , final_width = WG , **ld_WG)
    
    path_T.arc(Radius  , a2r(270) , a2r(360) , **ld_WG)
    path_T.arc(Radius   , a2r(180) , a2r(90) , **ld_WG)
    path_T.segment(length = 4*Radius*number_of_bend + 50, **ld_WG)
    path_T.arc(Radius , a2r(90) , a2r(0) , **ld_WG)
    path_T.arc(Radius, a2r(180) , a2r(270) , **ld_WG)
    
    path_B.arc(Radius, a2r(90) , a2r(0) , **ld_WG)
    path_B.arc(Radius , a2r(180) , a2r(270) , **ld_WG)
    for i in range(0,number_of_bend):
        path_B.segment(length = 50/(number_of_bend*2),direction = '+x' , **ld_WG)
        path_B.arc(radius = Radius , initial_angle = a2r(270) , final_angle = a2r(360) , **ld_WG)
        path_B.segment(length = Straight_in_bend/2 ,direction = '+y' , **ld_WG)
        path_B.arc(radius = Radius , initial_angle = a2r(180) , final_angle = a2r(90), **ld_WG)
        path_B.arc(radius = Radius , initial_angle = a2r(90) , final_angle = a2r(0), **ld_WG)
        path_B.segment(length = Straight_in_bend/2 ,direction = '-y' , **ld_WG)
        path_B.arc(radius = Radius , initial_angle = a2r(180) , final_angle = a2r(270), **ld_WG)
        path_B.segment(length = 50/(number_of_bend*2),direction = '+x' , **ld_WG)
    
    path_B.arc(Radius  , a2r(270) , a2r(360) , **ld_WG)
    path_B.arc(Radius  , a2r(180) , a2r(90) , **ld_WG)  
    
    path_B.segment(length = Broad_Length , direction = '+x' ,final_width = WG_Broad , **ld_WG)
    path_T.segment(length = Broad_Length , direction = '+x' ,final_width = WG_Broad , **ld_WG)
    
    x = path_T.x + MMI_L
    y = path_T.y - MMI_diss_From_center
    cell.add(gdspy.Rectangle((x , y - MMI_W/2), (x - MMI_L , y + MMI_W/2) , **ld_WG))
    path1 = gdspy.Path(width = WG_Broad , initial_point = (x,y))
    path1.segment(length = Broad_Length , direction = '+x' , final_width = WG , **ld_WG)
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path1.x , path1.y), direction='+x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    cell.add(path1)
    cell.add(path_T)
    cell.add(path_B)


def MMI1X4 (MMI_W = 10.5 , MMI_L = 45.75 , dis_from_edge = 0.5 , position = (0,0) , diss_bet_WG = 1):
    
    W = MMI_W
    L = MMI_L
    text = gdspy.Text("L = "+str(np.round(L,2)),20,(position[0] - 60, position[1] + 50))
    cell.add(text)
    text = gdspy.Text("W = "+str(np.round(W,2)),20,(position[0] - 60, position[1] + 30))
    cell.add(text)
    
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=position, direction='-x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    path1 = gdspy.Path(width = WG , initial_point = position)
    path1.segment(length = 20 , direction = "+x" , **ld_WG)
    path1.segment(length = Broad_Length , direction = "+x" , final_width = WG_Broad , **ld_WG)
    cell.add(path1)
    rec = gdspy.Rectangle((path1.x ,path1.y - MMI_W/2), (path1.x + MMI_L, path1.y + MMI_W/2) , **ld_WG)
    cell.add(rec)
    
    y_1 = path1.y + MMI_W/2 - dis_from_edge - WG_Broad/2
    y_4 = path1.y - MMI_W/2 + dis_from_edge + WG_Broad/2
    diff = (y_1 - y_4 - WG_Broad - WG_Broad*2)/3
    y_2 = y_1 - diff - WG_Broad
    y_3 = y_2 - diff - WG_Broad
    
    path_1 = gdspy.Path(width = WG_Broad , initial_point = (path1.x + MMI_L , y_1))
    path_1.segment(length = Broad_Length , direction = "+x" , final_width = WG , **ld_WG)
    path_1.arc(Radius  , a2r(270) , a2r(270 + Deg*2) , **ld_WG)
    path_1.arc(Radius  , a2r(90 + Deg*2) , a2r(90) , **ld_WG)
    path_1.segment(length = 100 , direction = "+x", **ld_WG)
    cell.add(path_1)
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path_1.x,path_1.y), direction='+x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    
    
    path_2 = gdspy.Path(width = WG_Broad , initial_point = (path1.x +  MMI_L , y_2))
    path_2.segment(length = Broad_Length , direction = "+x" , final_width = WG , **ld_WG)
    path_2.segment(length = 100  , direction = "+x", **ld_WG)
    path_2.arc(Radius  , a2r(270) , a2r(270 + Deg) , **ld_WG)
    path_2.arc(Radius  , a2r(90 + Deg) , a2r(90) , **ld_WG)
    path_2.segment(length = path_1.x - path_2.x , direction = "+x", **ld_WG)
    cell.add(path_2)
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path_2.x,path_2.y), direction='+x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    
    
    path_3 = gdspy.Path(width = WG_Broad , initial_point = (path1.x +  MMI_L , y_3))
    path_3.segment(length = Broad_Length , direction = "+x" , final_width = WG , **ld_WG)
    path_3.segment(length = 100 , direction = "+x", **ld_WG)
    path_3.arc(Radius , a2r(90) , a2r(90 - Deg) , **ld_WG)
    path_3.arc(Radius , a2r(270 - Deg) , a2r(270) , **ld_WG)
    path_3.segment(length = path_1.x - path_3.x , direction = "+x", **ld_WG)
    cell.add(path_3)
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path_3.x,path_3.y), direction='+x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    
    
    path_4 = gdspy.Path(width = WG_Broad , initial_point = (path1.x +  MMI_L , y_4))
    path_4.segment(length = Broad_Length , direction = "+x" , final_width = WG , **ld_WG)
    path_4.arc(Radius , a2r(90) , a2r(90 - Deg*2) , **ld_WG)
    path_4.arc(Radius , a2r(270 - Deg*2) , a2r(270) , **ld_WG)
    path_4.segment(length = 100 , direction = "+x", **ld_WG)
    cell.add(path_4)
    cell.add(grating(period=period, number_of_teeth=20, fill_frac=FF, width=21.5, position=(path_4.x,path_4.y), direction='+x', lda=1.55,
                             sin_theta=np.sin(np.pi * -8 / 180), focus_distance=21.5,focus_width= WG,tolerance=0.001,**ld_WG))
    
    


for i in range (0,len(Gap)):
    U_Shape_Rings(Straight_WG_Ring = 150,diss_in_bundel = 127 ,position = (300*i , 0) ,number_of_rings = 1 , gap = Gap[i] )
    U_Shape_Rings(Straight_WG_Ring = 150,diss_in_bundel = 127 ,position = (300*i , -450) ,number_of_rings = 2 , gap = Gap[i])

for i in range (0,len(Gap)):
    U_Shape_Disk(Straight_WG_Ring = 150,diss_in_bundel = 127 ,position = (900 + 300*i , 0) ,number_of_rings = 1 , gap = Gap[i] )
    U_Shape_Disk(Straight_WG_Ring = 150,diss_in_bundel = 127 ,position = (900 + 300*i , -450) ,number_of_rings = 2 , gap = Gap[i])


x = 850
y = 900

text = gdspy.Text("WG = "+str(0.45),20,(780,-950))
cell.add(text)
text = gdspy.Text("WG = "+str(1),20,(800,-1250))
cell.add(text)

for i in range (0,len(Gap)):
    Straight_Rings(Straight_WG = 100 ,position = (300*i , -(y + 1*Radius)) ,number_of_rings_ontop = 1 , gap = Gap[i] )
    Straight_Rings(Straight_WG = 100 ,position = (300*i , -(y + 6*Radius + 20)) ,number_of_rings_ontop = 2 , gap = Gap[i])



for i in range (0,len(Gap)):
    Straight_Rings(Straight_WG = 200 ,position = (1000 +300*i , -(y + 1*Radius)) ,number_of_rings_ontop = 1 , number_of_ring_side = 2 ,  gap = Gap[i] )
    Straight_Rings(Straight_WG = 200 ,position = (1000 +300*i , -(y + 6*Radius + 20)) ,number_of_rings_ontop = 2, number_of_ring_side = 2  , gap = Gap[i])

y = y + 300

for i in range (0,len(Gap)):
    Straight_Rings(WG = 1,Straight_WG = 100 ,position = (300*i , -(y + 1*Radius)) ,number_of_rings_ontop = 1 , gap = Gap[i] )
    Straight_Rings(WG = 1,Straight_WG = 100 ,position = (300*i , -(y + 6*Radius + 20)) ,number_of_rings_ontop = 2 , gap = Gap[i])



for i in range (0,len(Gap)):
    Straight_Rings(WG = 1 ,Straight_WG = 200 ,position = (1000 +300*i , -(y + 1*Radius)) ,number_of_rings_ontop = 1 , number_of_ring_side = 2 ,  gap = Gap[i] )
    Straight_Rings(WG = 1 ,Straight_WG = 200 ,position = (1000 +300*i , -(y + 6*Radius + 20)) ,number_of_rings_ontop = 2, number_of_ring_side = 2  , gap = Gap[i])


MZI(position = (1800 , -50))
MZI(Straight_in_bend = 30 ,position = (1800 , -200))

MZI(position = (1800 , -350) ,number_of_bend= 2)
MZI(Straight_in_bend = 30 , position = (1800 , -500) ,number_of_bend= 2)

MZI(position = (1800 , -650) , number_of_bend = 3)
MZI(Straight_in_bend = 30 , position = (1800 , -800) ,number_of_bend= 3)


for i in range(0,5):
    MMI1X4(MMI_W = 10.5 , MMI_L = 43.75 + i*1 , dis_from_edge = 0.5 ,position = (0 + 400*i , -1500) , diss_bet_WG = 1)

lib.write_gds("Cornerstone_Bar_Ilan_BD_Mask.gds")

