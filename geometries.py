import numpy as np

def calc_distance(coord1,coord2):
    dist = ((coord1[0]-coord2[0])**2 + 
            (coord1[1]-coord2[1])**2 + 
            (coord1[2]-coord2[2])**2)
    return(dist)

def calc_angle(coord1,coord2,coord3):
    "From the interwebz, needs fixing"
    a = [coord1]
    b = [coord2]
    c = [coord3]
    
    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return(np.degrees(angle))

def calc_torsion(coord1,coord2,coord3,coord4):
    "From the interwebz, needs fixing"
    p0 = coord1
    p1 = coord2
    p2 = coord3
    p3 = coord4

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)
    
    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    
    print(np.dot(b0, b1))
    print(np.dot(b2, b1))
    
    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    
    #print(np.arctan2(y, x))
    #print(np.degrees(np.arctan2(y, x)))
    return np.degrees(np.arctan2(y, x))