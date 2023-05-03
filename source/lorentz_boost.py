import numpy as np

c = 299792458 # speed of light (m/s) (NIST)

def boost(gamma, four_vector, x, y, z): # take the gamma factor of the frame, the 4-vector that you want to boost, and the unit vector of the frame that you are boosting into
    # Define the angles
    theta = np.arccos(z)
    if   x>0:
        phi = np.arctan(y/x)
    elif x<0  and y>=0:
        phi = np.arctan(y/x) + np.pi
    elif x<0  and y<0:
        phi = np.arctan(y/x) - np.pi
    elif x==0 and y>0:
        phi = np.pi/2
    elif x==0 and y<0:
        phi = -np.pi/2
    elif x==0 and y==0:
        phi = 0

    # Calculate velocity info from the provided gamma
    beta = np.sqrt(1-gamma**(-2))
    v = c * beta

    # Calculate velocity components for boost matrix
    vx = v * np.sin(theta) * np.cos(phi)
    vy = v * np.sin(theta) * np.sin(phi)
    vz = v * np.cos(theta)

    # Define the boost matrix
    boost_matrix = np.array([[ gamma     ,      -gamma*vx/c      ,      -gamma*vy/c      ,      -gamma*vz/c      ],
                             [-gamma*vx/c, 1+(gamma-1)*vx**2/v**2,   (gamma-1)*vx*vy/v**2,   (gamma-1)*vx*vz/v**2],
                             [-gamma*vy/c,   (gamma-1)*vy*vx/v**2, 1+(gamma-1)*vy**2/v**2,   (gamma-1)*vy*vz/v**2],
                             [-gamma*vz/c,   (gamma-1)*vz*vx/v**2,   (gamma-1)*vz*vy/v**2, 1+(gamma-1)*vz**2/v**2]])

    # Calculate the boost
    boosted_vector = np.dot(boost_matrix, four_vector)
    return boosted_vector

# Testing nonsense below
###vec = np.array([52*931.49+0.5*45**2/(52*931.49),0,0,45])
###v_new = boost(1.196999163270353, vec, 0, 0, -1)
###v_new
###
###52*931.45
###(v_new[3]**2)*0.5/(52*931.45)
###np.sqrt((v_new[0]**2)-(v_new[3]**2))
