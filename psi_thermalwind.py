import numpy as np
def psi_thermalwind(z, pb, pn, f):
    """
    Desription: 
    
    Computes a vertical profile of the overturning streamfunction in the North Atlantic 
    based on Nikurashin and Vallis, (2012) and Jansen et al., (2018) with boundary 
    conditions \psi = 0 at the surface and bottom of the ocean.
    
    Input:  
        z  = depth                                          [ m ]
        pn = vertical profile of density in northern region [ kg m^{-3} ]
        pb = vertical profile of density in basin           [ kg m^{-3} ]
        f  = coriolis parameter                             [ s^{-1} ]
    
    Output: 
        psi = vertical profile of overturning               [ Sv ]
    """
    
    # Remove bathymetric issues
    isn  = np.isnan(pn)
    zt   = z[~isn]
    pnt  = pn[~isn]
    pbt  = pb[~isn]
    
    # Constants 
    g = 9.81           # Gravitational acceleration
    rho0 = 1025        # Reference density
    cons = g/(f*rho0)  # Constant in thermal wind expression
    H    = zt[-1]      # Depth of ocean
    
    # Velocity grid
    zW   = zt

    # Number of vertical grid points, etc
    nz   = len(zt)
    zs   = np.zeros(nz+1)
    zs[1:nz] = 0.5 * (zt[1:]+zt[:-1])
    dz  = zs[1:] - zs[:-1]
    dz[-1] = dz[-2]

    # Compute first integral
    int_1 = np.zeros(nz)
    for i in range(1, nz):
        int_1[i] = int_1[i-1] + dz[i]*cons*(pbt[i-1]-pnt[i-1])
        
    # Compute second integral
    int_2  = np.zeros(nz)
    for i in range(1, nz):
        int_2[i] = int_2[i-1] + dz[i]*0.5*(int_1[i] + int_1[i-1])

    # Integrations terms
    C = -int_2[-1]/H 
    
    D = 0
    
    # Vertial profile of overturning
    psi = int_2 + C*zW + D

    return psi/1e6    
