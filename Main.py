import triangle as tr
import numpy as np
import matplotlib.pyplot as plt
from typing import Literal
import numpy.typing as npt
import math
from scipy.sparse.linalg import spsolve

#vertices = np.array(((0,0),(1,0),(1,1),(0,1)))
#A = dict(vertices)
#B = tr.triangulate(A)


def main():

    pts=np.array(((0, 0), (1, 0), (1, 1), (0, 1)))
    seg = [[0,1],[1,2],[2,3],[3,0]]
    markers = np.array([0,1,2,3])  # bottom, right, top, left
    A = dict(vertices=pts, segments=seg,segment_markers=markers)
    
    #pass p to input string to triangulate a planar straight line graph which will tallow segment_markers
    B = tr.triangulate(A,'pqa0.05e')  #qa0.01    0.01 is the max area of the triangle
    #tr.compare(plt, A, B)  #broken function in triangle, don't use
    vertices = B["vertices"]
    triangles= B["triangles"]
    boundary_segments = np.sort(B["segments"],axis=1) #sort to get same segments as those in faces
    
    boundary_markers = B["segment_markers"]
    #testz = [[1 ,1, 1], [2, 2, 2]]
    
    plt.triplot(B["vertices"][:,0],B["vertices"][:,1],triangles)
    
    
    #Triangles list is list of each cell
    #Make a list of all unique faces from the triangles list
    #faces will be stored as 2 pts in increasing index order
    #loop through all triangles and form faces to add to list
    #Each time face is added eg 2 pts are added, in 
    #array for cells per face, add cell index then after all cells have 
    #been run through, sort the cell indices from low to high
    
    #2 data structures for faces,
        #first is pts that make face, sorted
        #second is cells that border face, sorted
    faces = []
    boundary_faces = boundary_segments.tolist()
    face_cells = []
    face_distances = []
    cell_centers = []
    face_centers = []
    face_normals = []
    face_area = []
    
    calc_faces_connectivity(faces,face_cells,triangles)
    calc_cell_centers(cell_centers,vertices,triangles)
    
    cell_centers = np.array(cell_centers)
 
    #fig,ax = plt.subplots()
    for cell_ind,coords in enumerate(cell_centers):
        plt.text(coords[0],coords[1],str(cell_ind),ha="center",va="center",color="red",fontsize=10)
    plt.show()
    
    #maybe change to not take in np array
    calc_face_normals_area(face_normals,face_area,faces,face_cells,cell_centers,vertices)
    calc_face_centers(faces,vertices,face_centers)
    calc_face_distances(face_cells,face_centers,cell_centers,face_distances)
    
    
    face_cells = np.array(face_cells)
    #set boundary condition array for boundary faces
    
    num_faces = len(faces)
    boundary_index = np.zeros([num_faces,1])  #0 is interior, 1 is left, 2 is right, 3 is bottom, 4 is top
    set_boundary_index(boundary_index,faces,boundary_segments,boundary_markers)
    #I should make boundary index list based off of index in the face array and 
    #the value from segments that Triangle gives. I should also move boundary faces 
    #to separate array to iterate over so I dont need to go through all of them and constantly
    #check if its a boundary condition or not
    
    #convert to np arrays for efficiency, or maybe start with np arrays instead of lists
    faces = np.array(faces)

    
    num_cells = triangles.shape[0]
    print("Done Building Connectivity")
    
    #boundary conditions
    bc = [0,1,0,1] # 1000 pa exit and 1m/s entrance
    
    
    u = np.ones([num_cells])
    v = np.zeros([num_cells])
    p = np.zeros([num_cells])  
    
    #initialize starting values
    p += bc[1] #in place modification of p
    u -= 0.5

    
    #SIMPLE
    #velocity guess
    u_new,a_matrix = solve_momentum(face_cells,boundary_index,cell_centers,face_distances,face_area,face_normals,u,v,p,num_cells,bc,0)
    v_new,a_matrix = solve_momentum(face_cells,boundary_index,cell_centers,face_distances,face_area,face_normals,u,v,p,num_cells,bc,1)
    
    #solve pressure correction equation
    p_new = pressure_correction(a_matrix,face_cells)
    #adjust pressures and velocity using calculated pressure adjustment
    print("done")
    
    

def pressure_correction(a_matrix,num_cells,face_cells):
    a_diag = np.diagonal(a_matrix)
    b = np.zeros([num_cells])
    pressure_matrix = np.zeros([num_cells,num_cells])
    for face_ind,cells in enumerate(face_cells):
        if boundary_index[face_ind] == 0:
            for cell in cells:
                cell_A = cells[0]
                cell_B = cells[1]
            
                area = face_area[face_ind]
                
                a_A = a_diag[cell_A]
                a_B = a_diag[cell_B]
                
                pressure_matrix[cell_A,cell_A] += area/(a_A)*area
                pressure_matrix[cell_A,cell_B] += -area/(a_B)*area
                
                pressure_matrix[cell_B,cell_B] += area/(a_B)*area
                pressure_matrix[cell_B,cell_A] += -area/(a_A)*area
                
            
        else:
            
        
        print()
    

#this function should work for both u and v momentum since everything is same except unknown velocity
#currently set up for vetical inlet and outlet, need to change contribution to source for angled inlet/outlet
def solve_momentum(face_cells,boundary_index,cell_centers,face_distances,face_area,face_normals,u_old,v_old,p_old,num_cells,bc,norm_ind):
    
    momentum = np.zeros([num_cells,num_cells])
    source = np.zeros([num_cells])
    #print(face_cells)
    for face_ind,cells in enumerate(face_cells):
        #print(cells)
        if boundary_index[face_ind] == 0:
            for cell in cells:
                cell_A = cells[0]
                cell_B = cells[1]
                #calculate fluxes
                #approximate velocity at face using distance weighed average from the cell centers
                dist_A = face_distances[face_ind][0]
                dist_B = face_distances[face_ind][1]
                A_normal = face_normals[face_ind][0][norm_ind] #component of normal that is being used
                B_normal = face_normals[face_ind][1][norm_ind]
                #get component that goes in cell A and cell B
                area = face_area[face_ind]
                #avg_velocity = (velocity[cell_A]*dist_A+velocity[cell_B]*dist_B)/(dist_A+dist_B)
                
                #rho*u*u term
                A_comp_u = area*(u_old[cell_A]*dist_A/(dist_A+dist_B))
                B_comp_u = area*(u_old[cell_B]*dist_B/(dist_A+dist_B))
                
                #rho*v*u term
                A_comp_v = area*(v_old[cell_A]*dist_A/(dist_A+dist_B))
                B_comp_v = area*(v_old[cell_B]*dist_B/(dist_A+dist_B))
                
                #can be cleaned up to combine lines
                
                #pressure terms added to source
                #pressure_term = -(p_old[cell_A]*dist_A + p_old[cell_B]*dist_B)/(dist_A+dist_B)*area*normal
                
                source[cell_A] += -p_old[cell_A]*dist_A/(dist_A+dist_B)*area*A_normal
                source[cell_B] += -p_old[cell_B]*dist_B/(dist_A+dist_B)*area*B_normal
                
                momentum[cell_A,cell_A] += (A_comp_u+A_comp_v)*A_normal
                momentum[cell_A,cell_B] += (B_comp_u+B_comp_v)*A_normal
                momentum[cell_B,cell_B] += (B_comp_u+B_comp_v)*B_normal
                momentum[cell_B,cell_A] += (A_comp_u+A_comp_v)*B_normal
        else:  #boundary cell
            cell_A = cells[0]
            A_normal = face_normals[face_ind][0][0]
            match boundary_index[face_ind]:
                case 0: #wall
                    pass #0 flux so nothing to add to matrix
                case 1: #outlet
                    #normal gradient of velocity is 0 so v at wall is same as center, ish 
                    #flux = area*velocity[cell_A]
                    momentum[cell_A,cell_A] += area*u_old[cell_A]*A_normal
                    
                case 2: #wall
                    pass #0 flux
                    
                case 3: #inlet
                    #flux added to source vector
                    source[cell_A] = -bc[3]*A_normal  #goes in negative because on opposite side of = sign

    #calculate u_new
    u_new = spsolve(momentum,source)
    return u_new,momentum
                    
            
            
    return
            
            
def calc_face_centers(faces,vertices,face_centers):
    for face in faces:
        pt_a = face[0]
        pt_b = face[1]
        a_coords = vertices[pt_a]
        b_coords = vertices[pt_b]  #indexing with one value slices row as 1D array
        face_x = (a_coords[0]+b_coords[0])/2
        face_y = (a_coords[1]+b_coords[1])/2
        face_centers.append([face_x,face_y])
    return

def calc_face_distances(face_cells,face_centers,cell_centers,face_distances):
    for face_ind,cells in enumerate(face_cells):
        face_distance_row = []
        for cell in cells:
            face_pt = face_centers[face_ind]
            cell_pt = cell_centers[cell]
            distance = math.sqrt((face_pt[0]-cell_pt[0] )**2 + (face_pt[1]-cell_pt[1] )**2)
            face_distance_row.append(distance)
        face_distances.append(face_distance_row)
        
    return

def calc_faces_connectivity(faces,face_cells,triangles):
    for ind,cell in enumerate(triangles):
        tri_faces = calc_faces(cell)
        for face in tri_faces:
            try:
                found_index = faces.index(face)
            except ValueError:
                found_index = -1
            if found_index == -1:
                faces.append(face)
                face_cells.append([ind])
            else:
                face_cells[found_index].append(ind)
    #Dont need to sort because code iterates over cells in order but have it anyways
    #face_cells = [sorted(row) for row in face_cells]
    #fill 2nd position of border cells with -1 so list can be converted to numpy array
    for ind,row in enumerate(face_cells):
        
        if len(row) == 1:
            face_cells[ind].append(-1)
    
    return
                
def calc_faces(triangle):
    #pts are sorted from low to high
    triangle = sorted(triangle)
    face_a = [triangle[0],triangle[1]]
    face_b = [triangle[0],triangle[2]]
    face_c = [triangle[1],triangle[2]]
    faces = [face_a,face_b,face_c]
    return(faces)

def calc_cell_centers(cell_centers,vertices,triangles):
    for cell_ind,cell in enumerate(triangles):
            tri_pts = []
            for pt_ind,point in enumerate(cell):
                tri_pts.append(vertices[point,:].tolist())
            center = calc_center(np.array(tri_pts))
            cell_centers.append(center)

def calc_face_normals_area(face_normals,face_area,faces,face_cells,cell_centers: npt.NDArray[np.dtype],vertices):
    for face_ind,face in enumerate(faces):
        pt_a = vertices[face[0],:]
        pt_b = vertices[face[1],:]
        dx = pt_a[0]-pt_b[0]
        dy = pt_a[1]-pt_b[1]
        magnitude = np.sqrt(dx**2+dy**2)
        #store face area
        area = magnitude
        
        #works if we go through faces in order from 0 -> n
        face_area.append(area)
        
        guess_normal = np.array([dy/magnitude,-dx/magnitude])
        
        #get cell_center value for first cell in face_cells at specific index
        first_cell_center = cell_centers[face_cells[face_ind][0]]
        check_line = first_cell_center-pt_a
   
        angle = calc_angle(guess_normal,check_line)
        if abs(angle) < 90:
            #normal direction is correct for first cell
            face_normals.append([guess_normal,-guess_normal])
        else:
            face_normals.append([-guess_normal,guess_normal])
            
def set_boundary_index(boundary_index,faces,boundary_segments,boundary_markers):
    boundary_segments_list = boundary_segments.tolist()
    #compare given segments to faces because fewer comparisons and I know each one is in there, rest can be set to 0 as default
    for ind,face in enumerate(boundary_segments_list):
        idx = faces.index(face)
        boundary_value = boundary_markers[ind]
        boundary_index[idx] = boundary_value
    
    return
            
    
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


















