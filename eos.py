def eos(T, S, T0, S0, rho0, alpha, beta):
    """
    Description:
    
    Computes potential density using a linear equation of state.
        
    Input: 
        T    = temperature                            [ °C ]
        S    = salinity                               [ g/kg ]
        T0   = reference temperature                  [ °C ]
        S0   = reference salinity                     [ g/kg ]
        rho0 = reference density                      [ kg/m^3 ]
        
    Output:
        P  = potential density                        [ kg/m^3 ]
        PT = thermal component of potential density   [ kg/m^3 ]
        PH = haline component of potential density    [ kg/m^3 ]
    """
    
    # Thermal component 
    PT = rho0*(alpha*(T-T0)
    
    # Haline component 
    PH = rho0*(beta*(S-S0))
               
    # Total
    P  = rho0 + PT + PH
    
    return P, PT, PH
