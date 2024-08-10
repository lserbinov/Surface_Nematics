# grab self geometry index
node = hou.pwd()
geo = node.geometry()
import scipy.io

# Load the .mat file
path = '/home/cuncheng/Dev/active_nematics_shell/matlab/data/geo{:d}.mat'.format(int(hou.frame()))
print(path)
mat = scipy.io.loadmat(path)

# Extract the matrices
M = mat['M']
P = mat['P']
velocity = mat['velocity']
pressure = mat['pressure']

# define a funciton to create triangels
def addTriangle(p,q,r):
    # create polygon of 3 corners
    f = geo.createPolygon()
    f.addVertex(p)
    f.addVertex(q)
    f.addVertex(r)
 

# Loop through the face matrix M and vertex matrix P
points = []
for i in range(len(P)):
    # Create points
    point = geo.createPoint()
    point.setPosition(P[i])
    points.append(point)
print(M)
print(points[0])
print(points[M[0][0] - 1])
    # Create triangle faces
for j in range(len(M)):
    addTriangle(points[M[j][2]-1], points[M[j][1]-1], points[M[j][0]-1])
# addTriangle(points[M[0][0]-1], points[M[0][1]-1], points[M[0][2]-1])
