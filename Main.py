import triangle as tr
import numpy as np
import matplotlib.pyplot as plt

#vertices = np.array(((0,0),(1,0),(1,1),(0,1)))
#A = dict(vertices)
#B = tr.triangulate(A)



A = dict(vertices=np.array(((0, 0), (1, 0), (1, 1), (0, 1))))
B = tr.triangulate(A,'qa0.01')  #qa0.01    0.01 is the max area of the triangle
#tr.compare(plt, A, B)  #broken function in triangle, don't use
triangles= B["triangles"][:,2]
plt.triplot(B["vertices"][:,0],B["vertices"][:,1],B["triangles"])
plt.show()