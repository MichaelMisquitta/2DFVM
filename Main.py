import triangle as tr
import numpy as np
import matplotlib.pyplot as plt
from typing import Literal
import numpy.typing as npt
import math

#vertices = np.array(((0,0),(1,0),(1,1),(0,1)))
#A = dict(vertices)
#B = tr.triangulate(A)


def main():

    pts=np.array(((0, 0), (1, 0), (1, 1), (0, 1)))
    seg = [[0,1],[1,2],[2,3],[3,0]]
    markers = np.array([1,2,3,4])
    A = dict(vertices=pts, segments=seg,segment_markers=markers)
    
    #pass p to input string to triangulate a planar straight line graph which will tallow segment_markers
    B = tr.triangulate(A,'pqa0.1e')  #qa0.01    0.01 is the max area of the triangle
    #tr.compare(plt, A, B)  #broken function in triangle, don't use
    vertices = B["vertices"]
    triangles= B["triangles"]
    boundary_segments = np.sort(B["segments"],axis=1) #sort to get same segments as those in sides
    
    boundary_markers = B["segment_markers"]
    #testz = [[1 ,1, 1], [2, 2, 2]]
    
    plt.triplot(B["vertices"][:,0],B["vertices"][:,1],triangles)
    plt.show()
    
    #Triangles list is list of each cell
    #Make a list of all unique sides from the triangles list
    #sides will be stored as 2 pts in increasing index order
    #loop through all triangles and form sides to add to list
    #Each time side is added eg 2 pts are added, in 
    #array for cells per side, add cell index then after all cells have 
    #been run through, sort the cell indices from low to high
    
    #2 data structures for sides,
        #first is pts that make side, sorted
        #second is cells that border side, sorted
    sides = []
    side_cells = []
    cell_centers = []
    side_normals = []
    side_area = []
    
    calc_sides_connectivity(sides,side_cells,triangles)
    calc_cell_centers(cell_centers,vertices,triangles)
    
    cell_centers = np.array(cell_centers)
    
    #maybe change to not take in np array
    calc_side_normals(side_normals,side_area,sides,side_cells,cell_centers,vertices)
    
    
    sides = np.array(sides)
    side_cells = np.array(side_cells)
    #set boundary condition array for boundary sides
    
    num_sides = sides.shape[0]
    boundary_index = np.zeros([num_sides,1])  #0 is interior, 1 is left, 2 is right, 3 is bottom, 4 is top
    set_boundary_index(boundary_index,sides,boundary_segments)
    #I should make boundary index list based off of index in the side array and 
    #the value from segments that Triangle gives. I should also move boundary sides 
    #to separate array to iterate over so I dont need to go through all of them and constantly
    #check if its a boundary condition or not
    
    #convert to np arrays for efficiency, or maybe start with np arrays instead of lists
    sides = np.array(sides)

    
    num_cells = triangles.shape[0]
    print("Done Building Connectivity")
    
    #Data structures for flow variables
    u = np.array([num_cells,1])
    v = np.array([num_cells,1])
    p = np.array([num_cells,1])      
    
def calc_sides_connectivity(sides,side_cells,triangles):
    for ind,cell in enumerate(triangles):
        tri_sides = calc_sides(cell)
        for side in tri_sides:
            try:
                found_index = sides.index(side)
            except ValueError:
                found_index = -1
            if found_index == -1:
                sides.append(side)
                side_cells.append([ind])
            else:
                side_cells[found_index].append(ind)
    #Dont need to sort because code iterates over cells in order but have it anyways
    #side_cells = [sorted(row) for row in side_cells]
    #fill 2nd position of border cells with -1 so list can be converted to numpy array
    for ind,row in enumerate(side_cells):
        
        if len(row) == 1:
            side_cells[ind].append(-1)
    
    return
                
def calc_sides(triangle):
    #pts are sorted from low to high
    triangle = sorted(triangle)
    side_a = [triangle[0],triangle[1]]
    side_b = [triangle[0],triangle[2]]
    side_c = [triangle[1],triangle[2]]
    sides = [side_a,side_b,side_c]
    return(sides)

def calc_cell_centers(cell_centers,vertices,triangles):
    for cell_ind,cell in enumerate(triangles):
            tri_pts = []
            for pt_ind,point in enumerate(cell):
                tri_pts.append(vertices[point,:].tolist())
            center = calc_center(np.array(tri_pts))
            cell_centers.append(center)

def calc_side_normals(side_normals,side_area,sides,side_cells,cell_centers: npt.NDArray[np.dtype],vertices):
    for side_ind,side in enumerate(sides):
        pt_a = vertices[side[0],:]
        pt_b = vertices[side[1],:]
        dx = pt_a[0]-pt_b[0]
        dy = pt_a[1]-pt_b[1]
        magnitude = np.sqrt(dx**2+dy**2)
        #store side area
        area = magnitude
        
        #works if we go through sides in order from 0 -> n
        side_area.append(area)
        
        guess_normal = np.array([dy/magnitude,-dx/magnitude])
        
        #get cell_center value for first cell in side_cells at specific index
        first_cell_center = cell_centers[side_cells[side_ind][0]]
        check_line = first_cell_center-pt_a
   
        angle = calc_angle(guess_normal,check_line)
        if abs(angle) < 90:
            #normal direction is correct for first cell
            side_normals.append([guess_normal,-guess_normal])
        else:
            side_normals.append([-guess_normal,guess_normal])
            
def set_boundary_index(boundary_index,sides,boundary_segments,boundary_markers):
    boundary_segments_list = boundary_segments.tolist()
    #compare given segments to sides because fewer comparisons and I know each one is in there, rest can be set to 0 as default
    for ind,side in enumerate(boundary_segments_list):
        idx = boundary_segments_list.index(side)
        boundary_value = boundary_markers[idx]
        boundary_index
        
            
    
def calc_angle(line_a,line_b):
    mag_a = np.linalg.norm(line_a)
    mag_b = np.linalg.norm(line_b)
    dot = np.dot(line_a,line_b)
    angle_rad = np.arccos(np.clip(dot/(mag_a*mag_b),-1.0,1.0))
    return np.degrees(angle_rad)
    
def calc_center(tri_pts):
    x_center = np.mean(tri_pts[:,0])
    y_center = np.mean(tri_pts[:,1])
    center = [x_center,y_center]
    return center

def threshold(pt,val,thres):
    diff = abs(pt-val)
    if diff < thres:
        return True
    
    return False
    
main()


















