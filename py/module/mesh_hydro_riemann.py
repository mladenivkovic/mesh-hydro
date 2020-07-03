#!/usr/bin/env python3

# An exact Riemann solver and solution sampler


import numpy as np


gamma = 5./3.

GP1 = gamma + 1
GM1 = gamma - 1
GP1OGM1 = GP1/GM1
GM1OGP1 = 1/GP1OGM1
GM1HALF = 0.5 * GM1
alpha = GM1HALF/gamma

epsilon = 1e-6




def riemann_solver(rho, u, p, t):
    """
    Solve the riemann problem for ideal gases
    and sample the solution at time t.
    Assumes box size is = 1, 
    and that the left and right state are separated
    at index nx/2, where nx = rho.shape[0]

    rho, u, p are assumed to be numpy arrays of floats
    containing density, velocity, and pressure of the gas

    returns arrays of equal shapes with the sampled solution
    for density, velocity and pressure at time t
    """

    if t == 0:
        return rho, u, p

    nx = rho.shape[0]
    dx = 1./nx
    i = nx // 2 - 1 # in hydro_io.py/read_twostate_ic: rho[:nxhalf] = rhoL, nxhalf = nx // 2

    rhoL = rho[i]
    rhoR = rho[i+1]
    uL = u[i]
    uR = u[i+1]
    pL = p[i]
    pR = p[i+1]

    print("Solving Riemann problem with")
    print("rhoL = {0:12.6f}, rhoR = {1:12.6f}".format(rhoL, rhoR))
    print("  uL = {0:12.6f},   uR = {1:12.6f}".format(uL, uR))
    print("  pL = {0:12.6f},   pR = {1:12.6f}".format(pL, pR))
    print("  nx = ", nx)
    print("   t = ", t)

    # check if we have vacuum
    is_vacuum = False
    if rhoL == 0:
        is_vacuum = True
    elif rhoR == 0:
        is_vacuum = True
    else:
        # don't compute square roots of zero, so compute
        # soundspeeds only now
        aL = soundspeed(pL, rhoL)
        aR = soundspeed(pR, rhoR)
        if uR - uL >= 2/GM1 * (aL + aR):
            is_vacuum = True

    if not is_vacuum:
        pstar, ustar = find_star_state(rhoL, uL, pL, rhoR, uR, pR)

    rho_sol = np.empty(rho.shape, dtype = np.float)
    u_sol = np.empty(u.shape, dtype = np.float)
    p_sol = np.empty(p.shape, dtype = np.float)

    center = (nx//2)*dx # position of the center
    for i in range(nx):
        x = (i+0.5)*dx - center
        xt = x / t
        if not is_vacuum: 
            rho_sol[i], u_sol[i], p_sol[i] = sample_solution(rhoL, rhoR, uL, uR, ustar, pL, pR, pstar, xt)
        else:
            rho_sol[i], u_sol[i], p_sol[i] = sample_vacuum_solution(rhoL, rhoR, uL, uR, pL, pR, xt)

        

    return rho_sol, u_sol, p_sol   





def sample_solution(rhoL, rhoR, uL, uR, ustar, pL, pR, pstar, xt):
    """
    Sample the solution at the place xt = x/t
    """

    aL = soundspeed(pL, rhoL)
    aR = soundspeed(pR, rhoR)

    if xt < ustar:
        # we are in the left region
        psopL = pstar / pL

        if pstar > pL:
            # left shock
            SL = uL - aL * np.sqrt(0.5*GP1/gamma *psopL + alpha)
            if xt < SL:
                # outside left shock
                rho = rhoL
                u = uL
                p = pL
            else:
                # inside left shock
                rho = (psopL + GM1OGP1) / (GM1OGP1 * psopL + 1) * rhoL
                u = ustar
                p = pstar

        else:
            # left rarefaction
            SHL = uL - aL
            if xt < SHL:
                # outside the fan
                rho = rhoL
                u = uL
                p = pL
            else:
                astarL = aL * psopL ** (alpha)
                STL = ustar - astarL
                if xt > STL:
                    # in central region outside fan
                    rho = rhoL * psopL ** (1./gamma)
                    u = ustar
                    p = pstar
                else:
                    # inside the fan
                    fact = ( 2 / GP1 + GM1OGP1 / aL *(uL - xt) )**(2./GM1)
                    rho = rhoL * fact
                    u =  2 / GP1 * (GM1HALF * uL + aL + xt)
                    p = pL * fact**gamma


    else:
        # we are in the right region
        psopR = pstar / pR

        if pstar > pR:
            # right shock
            SR = uR + aR * np.sqrt(0.5*GP1/gamma *psopR + alpha)
            if xt > SR:
                # outside right shock
                rho = rhoR
                u = uR
                p = pR
            else:
                # inside right shock
                rho = (psopR + GM1OGP1) / (GM1OGP1 * psopR + 1) * rhoR
                u = ustar
                p = pstar

        else:
            # right rarefaction
            SHR = uR + aR
            if xt > SHR:
                # outside the fan
                rho = rhoR
                u = uR
                p = pR
            else:
                astarR = aR * psopR**alpha
                STR = ustar + astarR
                if xt < STR:
                    # in central region outside fan
                    rho = rhoR * psopR ** (1./gamma)
                    u = ustar
                    p = pstar
                else:
                    # inside the fan
                    fact = ( 2 / GP1 - GM1OGP1 / aR *(uR - xt) )**(2/GM1)
                    rho = rhoR * fact
                    u =  2 / GP1 * (GM1HALF * uR - aR + xt)
                    p = pR * fact**gamma


    return rho, u, p






def find_star_state(rhoL, uL, pL, rhoR, uR, pR):
    """
    Find the star state pressure and velocities following Toro 1999
    returns: pstar, ustar: stare state pressure and velocity
    """

    aL = soundspeed(pL, rhoL)
    aR = soundspeed(pR, rhoR)
    AL = A_K(rhoL)
    AR = A_K(rhoR)
    BL = B_K(pL)
    BR = B_K(pR)

    # find initial guess for pstar.
    # use Two Rarefaction Approximation
    ppv = 0.5*(pL + pR) - 0.125 * (uR - uL)*(rhoL + rhoR) * (aL + aR)
    pstar = max(epsilon, ppv)
    # use two-shock guess (needs two rarefaction solution as well)
    #  gL = np.sqrt(AL / (pstar + BL))
    #  gR = np.sqrt(AR / (pstar + BR))
    #  pstar = (gL*pL + gR*pR + uL - uR)/(gL + gR)
    #  pstar = max(epsilon, pstar)

    diff = 1
    i = 0
    while diff > epsilon:

        f = f_K(pstar, pL, AL, BL, aL) + f_K(pstar, pR, AR, BR, aR) + uR - uL
        dfdp = df_Kdp(pstar, rhoL, pL, AL, BL, aL) + df_Kdp(pstar, rhoR, pR, AR, BR, aR)
        pstar_new = pstar - f / dfdp
        diff = 2 * abs(pstar_new - pstar) / abs(pstar_new + pstar)
        pstar = pstar_new

        if i > 1000:
            print("Got 1k iterations for exact riemann solver. Exiting now")
            break
        i+=1

        # don't allow negative pressure
        if pstar < epsilon:
            pstar = epsilon

    print("Found star state pressure after", i, "iterations")


    ustar = uL - f_K(pstar, pL, AL, BL, aL)
    print("Got pstar = {0:12.6f}, ustar = {1:12.6f}".format(pstar, ustar))

    return pstar, ustar

    




def f_K(pstar, pK, AK, BK, aK):
    """
    Compute f_{K=L, R} (See section 3.2 in equations_and_implementation_details.pdf)
    pstar:  p in star region
    pK: pLeft or pRight
    AK : 2/(gamma + 1) / rhoK
    BK : (gamma - 1)/(gamma + 1) * pK
    aK : sound speed in region K
    """
    if pstar > pK:
        # shock relation
        return (pstar - pK) * np.sqrt(AK/(pstar + BK))
    else:
        # rarefaction relation
        return 2 * aK / GM1 * ((pstar/pK)**alpha - 1)






def df_Kdp(pstar, rhoK, pK, AK, BK, aK):
    """
    Compute  del f_{K=L, R}/dp (See section 3.2 in equations_and_implementation_details.pdf)
    pstar:  p in star region
    pK: pLeft or pRight
    AK : 2/(gamma + 1) / rhoK
    BK : (gamma - 1)/(gamma + 1) * pK
    aK : sound speed in region K
    """
    if pstar > pK:
        # shock relation
        return ( 1. - 0.5 * (pstar - pK)/(pstar + BK) ) * np.sqrt(AK/(pstar + BK))
    else:
        # rarefaction relation
        return 1./(aK * rhoK) * (pstar/pK)**(-0.5*GP1/gamma)





def A_K(rhoK):
    """
    Compute A_{L, R}
    """
    return 2 / (GP1 * rhoK)





def B_K(pK):
    """
    Compute B_{L,R}
    """
    return pK / GP1OGM1





def soundspeed(p, rho):
    """
    Compute the sound speed of the gas
    """
    return np.sqrt(p * gamma / rho)





def sample_vacuum_solution(rhoL, rhoR, uL, uR, pL, pR, xt):
    """
    Sample the solution in the presence of vacuum
    """

    if rhoL == 0 and rhoR == 0:
        return 0., 0., 0.

    if rhoL == 0:
        # left vacuum state
        aR = soundspeed(pR, rhoR)
        SR = uR - 2 * aR / GM1
        SHR = uR + aR

        if xt <= SR:
            # left vacuum
            rho = 0
            u = SR
            p = 0
        elif xt < SHR:
            # inside right rarefaction
            fact = ( 2 / GP1 - GM1OGP1 / aR * (uR - xt) )**(2/GM1)
            rho = rhoR * fact
            u =  2 / GP1 * (GM1HALF * uR - aR + xt)
            p = pR * fact**gamma
        else:
            # in right state
            rho = rhoR
            u = uR
            p = pR

    elif rhoR == 0:
        # Right vacuum state
        aL = soundspeed(pL, rhoL)
        SL = uL + 2 * aL / GM1
        SHL = uL - aL

        if xt >= SL:
            rho = 0
            u = SL
            p = 0
        elif xt > SHL:
            fact = ( 2 / GP1 + GM1OGP1 / aL *(uL - xt) )**(2/GM1)
            rho = rhoL * fact
            u =  2 / GP1 * (GM1HALF * uL + aL + xt)
            p = pL * fact**gamma
        else:
            rho = rhoL
            u = uL
            p = pL
    else:
        # Vacuum generating state
        aL = soundspeed(pL, rhoL);
        aR = soundspeed(pR, rhoR);
        SL = uL + 2*aL/GM1
        SR = uR - 2*aR/GM1
        SHL = uL - aL;
        SHR = uR + aR;

        if xt <= SHL:
            # outside left rarefaction, original state
            rho = rhoL
            u = uL
            p = pL
        elif xt < SL:
            # inside left rarefaction fan
            fact = ( 2 / GP1 + GM1OGP1 / aL * (uL - xt) )**(2/GM1)
            rho = rhoL * fact
            u =  2 / GP1 * (GM1HALF * uL + aL + xt)
            p = pL * fact**gamma
        elif xt < SR:
            # vacuum region
            rho = 0
            u = 0.5 * (SL + SR)
            p = 0
        elif xt < SHR:
            # inside right rarefaction fan
            fact = ( 2 / GP1 - GM1OGP1 / aR * (uR - xt) )**(2/GM1)
            rho = rhoR * fact
            u =  2 / GP1 * (GM1HALF * uR - aR + xt)
            p = pR * fact**gamma
        else:
            # right original state
            rho = rhoR
            u = uR
            p = pR

    return rho, u, p
