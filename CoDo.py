# CoDo has been developed by Daina Romeo and is available as open source model.

import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time
from datetime import date
import pyautogui
import sys

# PLEASE NOT THAT THE FILEPATH USES "/" AND NOT "\"
# -------- INDICATE THE FULL PATH OF THE MPPD DATASET CSV FILE -------------------------

MPPD_dataset = 'C:/Users/roda/Documents/Python Scripts/MPPD_df_auto.csv'

# ------ INDICATE THE LOCATION OF THE MPPD SOFTWARE FOLDER - USUALLY FOUND IN DOCUMENTS
MPPD_results = 'C:/Users/roda/Documents/MPPD'

# ------- INDICATE THE FULL PATH OF THE TEMPLATE FILE --------------------------------------
template_file = 'C:/Users/roda/Documents/Python Scripts/CoDo/CoDo template Veno.xlsx'

# ------- INDICATE THE FULL PATH OF THE MPPD ICON IMAGE --------------------------------

MPPD_icon = 'C:/Users/roda/Documents/Python Scripts/MPPD_icon.png'

# - INDICATE THE FIRST AND LAST LINE TO COMPUTE (IF ONLY ONE LINE, START_VALUE=END_VALUE)
start_value = 4
end_value = 43


# ------------------------------------------------------------------------------------
#            END OF USER'S INPUTS - DO NOT MODIFY AFTER THIS LINE
# ------------------------------------------------------------------------------------

####Added: check time at which deposited concentration becomes constant and returns it

def in_vitro_simulation(ID=None, Output_file='missing name', Substance=None,
                        media_viscosity=0.00089, media_density=1.0104, media_temperature=37,
                        pp_density=None, agg_diameter=None,
                        agg_effective_density=None, pp_diameter=None, column_height=None,
                        initial_concentration=None,
                        simulation_time=None, subcompartment_height=0.005, simulation_time_interval=0.5,
                        output_time_interval=30, eq_change_rate=0.002, output_compartment_height=0.01,
                        to_plot=False, save=False, sedimentation_cdependence=0.0,
                        diffusion_cdependence=0.0,
                        initial_dissolution=0, dissolution_rate_type=0, dissolution_rate=0,
                        dissolution_times=0, dissolution_fractions=0, sticky=False,
                        adsorption_dissociation_constant=1e-9, Eq=False, **kwargs):
    Today = date.today()
    Date_to_report = Today.strftime("%d-%m-%Y")
    if Output_file == 'missing name':
        Output_file = f"{Substance}_{ID}_{Date_to_report}"

    # The in vitro simulation is based on the DG_NANOTRANSPORT_SIMULATE - DISTORTED GRID NANOTRANSPORT SIMULATOR
    # from   Glen DeLoid 07/13/15   (Copyright 2015, Harvard T. H. Chan School of Public Health)
    # Which has been converted to Python programming language and modified by Daina Romeo, 12/10/2020

    # Input data - constants

    # define constants

    # Boltzman constant
    kB = 1.3806503e-23  # m^2 kg s^-2 K^-1
    # absolute zero (degrees C below 0C)
    absZ = 273.15
    # gravitational acceleration (m/s^2)
    Ga = 9.8
    # Na Avigadro's number
    Na = 6.0221413e23

    # Media and Particle/agglomerate variables are stored in a dictionary
    Media = dict()
    P = dict()

    # Media properties
    # absolute temperature
    Media['absT'] = absZ + media_temperature
    # media density ( in kg/m^3)
    Media['dens'] = media_density * 1.0e3

    # Particle/Agglomerate-related properties

    # agglomerate radius in m
    P['rad'] = 0.5e-9 * agg_diameter
    # agglomerate cross-sectional areas
    P['CrossArea'] = P['rad'] ** 2 * math.pi

    # density of material (convert to mg/cm3 == kg/m3)
    P['pp_density'] = pp_density * 1.0e3
    # agglomerate effective density converted in kg/m^3 e.g. g/cm3 x 1000
    P['agg_effective_density'] = agg_effective_density * 1.0e3
    # mass of the agglomerates (kg)
    P['mass'] =  P['agg_effective_density'] * (4 / 3) * math.pi * P['rad'] ** 3
    # agglomerate molar concentration per mass concentration
    # n/m^3 = (kg/m^3)/kg = C/P.mass
    # n/L = 1e-3 * n/m^3 = 1e-3*C/P.mass
    # moles/L = (n/L)/Na = (1/Na)*1e-3*C/P.mass
    P['molpermass'] = 1e-3 / (P['mass'] * Na)
    # agglomerate cross-sectional area per unit mass
    P['CrossAreaMass'] = P['CrossArea'] / P['mass']
    # Coefficient for sedimentation concentration dependence
    P['SedCod'] = sedimentation_cdependence
    # Coefficient for diffusion concentration dependence (Dc = Do/(1 + Kd*c)), usually 0.1, or 0;
    P['DiffCod'] = diffusion_cdependence
    # type of dissolution rate
    # 0 = no further dissolution after initial dissolution (stable)
    # 1 = fraction of original per hour (constant)
    # 2 = specified times and fractions (specified curve, interpolated linearly)
    P['DissRateType'] = int(dissolution_rate_type)
    # initial dissolution fraction
    P['initial_dissolution'] = initial_dissolution
    # rate of dissolution (after initial)
    # only relevant for rate_type = 0 or 1
    P['DissRate'] = dissolution_rate
    # times for dissolution fraction data (h)
    # only relevant for rate_type = 2;
    P['DissT'] = dissolution_times
    # dissolution fractions corresponding to specified times
    P['TDissFrx'] = dissolution_fractions

    # Simulation/experimental parameters

    # column height (m)
    colH = column_height * 1.0e-3
    # compartment height(m)
    subH = subcompartment_height * 1.0e-3
    # output data compartment height
    BottomH = output_compartment_height * 1.0e-3
    # number of compartments - use floor to round
    # top of top compartment will be at or <subH below top of column
    # requires assumption of initial homogeneous mixture
    ncomp = round(colH / subH)
    # number of compartments to include in "bottom fraction"
    BottN = round(BottomH / subH)
    # number of BottomH sized compartments (combined/averaged compartments)
    ncompBott = ncomp / BottN
    # simulation time (sec)
    Tend = simulation_time * 3600
    # simulation time interval(sec)
    dt = simulation_time_interval
    # initial total concentration (mg/ml) (same as kg/m^3) of primary material
    C0p = initial_concentration
    # output time interval (min)
    OutputIntervalMin = output_time_interval
    # output time interval (h)
    OutputIntervalHrs = OutputIntervalMin / 60.0
    # number of output points
    NOutPoints = math.ceil(simulation_time * 60 / OutputIntervalMin) + 1
    # N x g
    Ng = 1

    # material/media/agglomerate fractional and multiplier values

    # Fraction of agglomerate that is particle
    P['PFracAgg'] = (P['agg_effective_density'] - Media['dens']) / (P['pp_density'] - Media['dens'])
    # fraction in agglomerate that is media
    P['MediaFracAgg'] = 1 - P['PFracAgg']
    # mass of media per unit mass (e.g. kg) of material (in agglomerate)
    P['mm'] = (P['MediaFracAgg'] * Media['dens']) / (P['PFracAgg'] * P['pp_density'])
    # total mass of agglomerate per mass of material (multiplier for nm
    # concentration for calculating agglomerate mass/concentration for a given
    # raw material mass/concentration). Equal to the mass of media per unit mass
    # of material (P.mm) plus the unit mass (1.0) of material
    P['Magg'] = P['mm'] + 1.0

    # initial transport coefficients for all agglomerate size species

    # agglomerate sed coeff (at C=0)
    P['Scoeff'] = 2.0 * (P['rad'] ** 2) * (P['agg_effective_density'] - Media['dens']) / (9.0 * media_viscosity)  # sec
    # agglomerate diff coeff (at C=0)
    P['Dcoeff'] = (kB * Media['absT']) / (6.0 * math.pi * media_viscosity * P['rad'])
    # P.dcoeff = (kB * Sol.tabs * P.DiffCorrGD)./(6.0 * pi * Sol.visc * P.rad);

    # initial concentration of agglomerate species and dissolved material

    # initial agglomerates concentrations without dissolution
    C0a = C0p * P['Magg']
    # Initial agglomerates concentration when accounting for initially dissolved fraction
    C0a = C0a * (1.0 - P['initial_dissolution'])
    # initial concentration of dissolved fraction
    DissC = C0p * P['initial_dissolution']
    # Starting dissolved and undissolved fractions (above and beyond initial
    # dissolution fraction)
    DissFrx = P['initial_dissolution']
    UndissFrx = 1.0 - P['initial_dissolution']

    # initial bound concentrations of agglomerates at bottom
    Cabound = np.array([0])
    # initial fraction of bottom area occupied
    FrxBottomOccupied = 0.0
    # kill back diffusion from bottom
    BotDoff = 0

    # Surface adsorption (Langmuir) approach to stickiness

    # Surface molar dissociation constant (mol m
    Kd = adsorption_dissociation_constant

    # starting compartment agglomerate concentrations - equal concentration in every compartment
    Ca = np.tile(C0a, [ncomp, 1])

    # total mass of material per area
    MPAp = C0p * ncomp

    # initial concentration-adjusted sedimentation coefficients
    # (at centers of compartments, always)
    S = P['Scoeff'] / (1 + P['SedCod'] * Ca)

    # initial concentration-adjusted diffusion coefficients
    # at boundaries
    D = np.zeros((ncomp + 1, 1))

    for k in range(1, ncomp):
        D[k, :] = P['Dcoeff'] / (1 + P['DiffCod'] * np.mean(Ca[(k - 1):(k + 1), ], axis=0))

    if BotDoff != 0:
        D[ncomp - 1, ] *= 0.5

    # Z positions of subcompartment centers for output data (mm)
    # gets middle position, than adds the subcompartment height to go to next middle position

    Z = np.arange(subH / 2, (ncomp * subH), subH) * 1e3

    # Calculate maximum fractional compartment distance sedimented in dt
    MaxSdx = abs(P['Scoeff']) * (Ng * Ga * dt / subH)
    # Calculate max fractional concentration change from diffusion in dt
    MaxDdcdx = P['Dcoeff'] * (dt / subH ** 2)
    # if either max is above 0.5 (.49 to be safe), throw error and recommend
    if MaxSdx >= 0.49 and MaxDdcdx >= 0.49:
        # calculate maximum time for sedimentation
        maxdts = 0.49 * subH / (Ng * Ga * abs(P['Scoeff']))
        # calculate maximum time for diffusion
        maxdtd = (0.49 * subH ** 2) / P['Dcoeff']
        # maximum time for a step (maxdt) is min of max diffusion and sedimentation
        maxdt = min(maxdts, maxdtd)
        raise ValueError(f'the simulation time interval is too long, please use a value < {maxdt} seconds')


    # diffusion concentration change calculation constant
    Diffint = dt / subH ** 2
    # sedimentation distance sedimented calculation constant
    Sedint = Ng * Ga * dt / subH
    # SedCoeffConst (s = SConst * rad^2)
    SConst = 2.0 * (P['agg_effective_density'] - Media['dens']) / (9.0 * media_viscosity)
    # DiffCoeffConst (D = DConst/rad)
    DConst = (kB * Media['absT']) / (6.0 * math.pi * media_viscosity)
    # number of time intervals (counts) before output
    outcount = int(round(OutputIntervalMin * 60 / dt))

    # output time counter
    outcounter = 0

    # Initialize output data time point counter
    I = 0
    # Initial time in seconds
    T = 0.0

    ##################################################################################################################
    # Simulation at T0
    # At T=0, adjust for dissolution if needed and output first data point

    ###Create empty result lists of lists to set size before loop

    ConcPerBottHcomp = []
    MassPerBottHcomp = []
    DepMassBott = []
    DepMassDissBottcomp = []
    NBottcomp = []
    NDepBottcomp = []
    OutSATpdzbot = []
    SADepBottcomp = []
    OutFrxMassdzbot = []
    OutT = []
    OutCa = []
    OutSumCa = []
    OutCp = []
    OutTotalCp = []
    OutDissC = []
    OutFrxOcc = []
    OutMBound = []
    OutTotalCpd = []
    OutFrxMass = []
    DepMasscomp= []
    DepMassDisscomp = []
    OutNp = []
    OutNTp = []
    OutNTpDep = []
    NDepBottcomp = []
    OutSAp = []
    OutSATp = []
    OutSATpDep = []
    OutCpdzbot = []
    OutMpdzbot = []
    OutNpdzbot = []
    OutNpDepdzbot = []
    OutSApdzbot = []
    OutSADepdzbot = []

    # Perform Dissolution (if necessary) for T=0
    if P['DissRateType'] > 0:
        diss_result = perform_dissolution(P, T, C0p, UndissFrx, Ca, SConst, DConst)
        DissC = diss_result[0]
        Ca = diss_result[1]
        P['rad'] = diss_result[2]
        P['Scoeff'] = diss_result[3]
        P['Dcoeff'] = diss_result[4]
        DissFrx = diss_result[5]
        UndissFrx = diss_result[6]

    # Calculate output point
    result = calculate_output_point(P, T, Ca, DissC, FrxBottomOccupied, Cabound,
                                    subH, ncompBott, BottN, MPAp, ncomp)

    OutT.append(result[0])
    OutCa.append(result[1])
    OutSumCa.append(result[2])
    OutCp.append(result[3])
    OutTotalCp.append(result[4])
    OutDissC.append(result[5])
    OutFrxOcc.append(result[6])
    OutMBound.append(result[7])
    ConcPerBottHcomp.append(result[8])
    OutTotalCpd.append(result[9])
    MassPerBottHcomp.append(result[10])
    OutFrxMass.append(result[11])
    OutFrxMassdzbot.append(result[12])
    DepMasscomp.append(result[13])
    DepMassBott.append(result[14])
    DepMassDisscomp.append(result[15])
    DepMassDissBottcomp.append(result[16])
    OutNp.append(result[17])
    OutNTp.append(result[18])
    NBottcomp.append(result[19])
    OutNTpDep.append(result[20])
    NDepBottcomp.append(result[21])
    OutSAp.append(result[22])
    OutSATp.append(result[23])
    OutSATpdzbot.append(result[24])
    OutSATpDep.append(result[25])
    SADepBottcomp.append(result[26])
    OutCpdzbot.append(result[27])
    OutMpdzbot.append(result[28])
    OutNpdzbot.append(result[29])
    OutNpDepdzbot.append(result[30])
    OutSApdzbot.append(result[31])
    OutSADepdzbot.append(result[32])

    # Increment output data time point counter
    I = I + 1

    ## variables to check for reaching equiibrium in deposition over time
    eq_count = []

    print(f'Loop starts - particle ID {ID}')

    #########START LOOP#################################################################################################
    # Loop through time in seconds

    while T <= Tend:
        # Perform round of SEDIMENTATION

        # calculate sedimentation distances (in compartment heights (dx's)) of compartment
        # center= (note, compartment center!) 'pseudo boundaries.'
        Sd = S * Sedint

        CaSd = Ca * Sd
        Ca1mSd = Ca * (1 - Sd)
        firstCapr = Ca1mSd[0, np.newaxis]
        middleCapr = CaSd[0:ncomp - 2] + Ca1mSd[1:ncomp - 1]
        lastCapr = (CaSd[ncomp - 2] + Ca[ncomp - 1])[np.newaxis, :]
        Capr = np.concatenate((firstCapr, middleCapr, lastCapr), axis=0)

        # to allow diffusion to act on same concentrations used in
        # sedimentation
        # concentration deltas for sedimentation
        dCaS = Capr - Ca

        # Perform round of DIFFUSION

        # concentration differences across boundaries (i+1 - i)
        CaDiff = Ca[1:] - Ca[0:ncomp - 1]

        if sticky:
            CaDiff[-1] = CaDiff[-1] - Cabound

        # diffusion for first compartment (movement only from 2nd compartment
        # back to first (this) compartment)
        Capr[0] = Ca[0] + (Diffint * D[1] * CaDiff[0])

        # diffusion for last compartment (movement only from this compartment
        # to compartment above) with CaDiff adjusted for Cabound

        Capr[-1] = Ca[-1] - (Diffint * D[ncomp - 1] * CaDiff[-1])

        # diffusion for all other compartments (from next compartment (+) into
        # current compartment and from current compartment (-) into compartment
        # above
        Capr[1:-1] = Ca[1:-1] + Diffint * (
                (D[2:ncomp] * CaDiff[1:(ncomp - 1)]) - (D[1:(ncomp - 1)] * CaDiff[0:(ncomp - 2)]))

        # to allow diffusion to act on same concentrations used in
        # sedimentation, set Ca to diffused Capr and add sed deltas
        Ca = Capr + dCaS
        for conc, name in zip([Ca, Capr, dCaS, Sd], ['Ca', 'Capr', 'dCaS', 'Sd']):
            if np.isnan(conc).any():
                raise ValueError(name)
        # Adjust for dissolution if needed
        if P['DissRateType'] > 0:
            diss_result = perform_dissolution(P, T, C0p, UndissFrx, Ca, SConst, DConst)
            DissC = diss_result[0]
            Ca = diss_result[1]
            P['rad'] = diss_result[2]
            P['Scoeff'] = diss_result[3]
            P['Dcoeff'] = diss_result[4]
            DissFrx = diss_result[5]
            UndissFrx = diss_result[6]

            # Update bound at bottom
        if sticky:
            bottom_result = update_bottom_bound_by_adsorption(P, Ca, Kd, subH)
            MaxFrxBottomOccupied = bottom_result[0]
            Frxavail = bottom_result[1]
            Frxbound = bottom_result[2]
            FrxBottomOccupied = bottom_result[3]
            Cabound = bottom_result[4]

        # Adjust sedimentation and diffusion coefficients for concentration

        # new concentration-adjusted sedimentation coefficients
        # at centers of compartments
        if P['SedCod'] != 0.0:
            S = P['Scoeff'] / (1 + P['SedCod'] * Ca)

        # new concentration-adjusted diffusion coefficients
        # at boundaries
        if P['DiffCod'] != 0.0:
            for k in range(1, ncomp):
                D[k] = P['Dcoeff'] / (1 + P['DiffCod'] * np.mean(np.array([Ca[k - 1], Ca[k]]), axis=0))

        # other concentration or time-dependent changes can be made here

        # Increment time and counter
        T += dt
        outcounter += 1
        # Calculate and plot output data if time
        # reset counter each time it reaches the reporting time, and calculate and save output
        if outcounter == outcount:
            outcounter = 0

            #     if mod(T,60*OutputIntervalMin)==0
            #         # if exactly at output interval mark

            # Calculate output point
            result = calculate_output_point(P, T, Ca, DissC, FrxBottomOccupied,
                                            Cabound, subH, ncompBott, BottN, MPAp, ncomp)
            # raise error if variables become nan - indicates non-convergence
            for conc, name in zip([T, Ca, DissC, FrxBottomOccupied,
                                   Cabound, subH, ncompBott, BottN, MPAp, ncomp], ['T', 'Ca', 'DissC', 'FrxBottomOccupied',
                                                                              'Cabound', 'subH', 'ncompBott', 'BottN',
                                                                              'MPAp', 'ncomp']):
                if np.isnan(conc).any():
                    raise ValueError(name)

            OutT.append(result[0])
            OutCa.append(result[1])
            OutSumCa.append(result[2])
            OutCp.append(result[3])
            OutTotalCp.append(result[4])
            OutDissC.append(result[5])
            OutFrxOcc.append(result[6])
            OutMBound.append(result[7])
            ConcPerBottHcomp.append(result[8])
            OutTotalCpd.append(result[9])
            MassPerBottHcomp.append(result[10])
            OutFrxMass.append(result[11])
            OutFrxMassdzbot.append(result[12])
            DepMasscomp.append(result[13])
            DepMassBott.append(result[14])
            DepMassDisscomp.append(result[15])
            DepMassDissBottcomp.append(result[16])
            OutNp.append(result[17])
            OutNTp.append(result[18])
            NBottcomp.append(result[19])
            OutNTpDep.append(result[20])
            NDepBottcomp.append(result[21])
            OutSAp.append(result[22])
            OutSATp.append(result[23])
            OutSATpdzbot.append(result[24])
            OutSATpDep.append(result[25])
            SADepBottcomp.append(result[26])
            OutCpdzbot.append(result[27])
            OutMpdzbot.append(result[28])
            OutNpdzbot.append(result[29])
            OutNpDepdzbot.append(result[30])
            OutSApdzbot.append(result[31])
            OutSADepdzbot.append(result[32])

            # Check if the rate of deposition (as deposited fraction tx+1 - deposited fraction tx) is below the threshold
            # defined by eq_change_rate

            if (OutFrxMassdzbot[-1][-1, -1] - OutFrxMassdzbot[-2][-1, -1]) / OutputIntervalHrs <= eq_change_rate:
                # print(OutFrxMassdzbot[-1][-1, -1] - OutFrxMassdzbot[-2][-1, -1])
                eq_count.append(I)

            # Increment output data time point counter
            I = I + 1

        print(f'ID {ID}, end loop {T},calculating point {I}')

        # increase simulation time if deposition rate is not under the threshold
        if T >= Tend:
            if not eq_count and Eq:
                Tend += OutputIntervalMin * 60
            else:
                if not Eq:
                    eq_time = 'not defined'
                elif len(eq_count) == len(OutT[eq_count[0]:]):
                    begin_eq = eq_count[0]
                    eq_time = OutT[begin_eq]
                else:
                    print(
                        f'Deposited fraction decreases over time after reaching threshold value. Please check it out. ID {ID}')
                    eq_time = 'not defined - decresing deposited fraction'
    if to_plot:
        plot_output(ConcPerBottHcomp, OutputIntervalHrs, DepMassBott, OutFrxMassdzbot,
                    Output_file, SADepBottcomp)

    print(f'ID {ID}: simulation complete')

    # Save output data
    if save:
        write_output_data(OutT, sticky, OutFrxMassdzbot, DepMassBott, DepMassDissBottcomp, NDepBottcomp,
                          SADepBottcomp, OutDissC, OutFrxOcc, OutMBound, Output_file, eq_time, eq_change_rate)

    if not Eq:
        eq_frac, eq_mass, eq_sa, eq_num = np.nan, np.nan, np.nan, np.nan
    else:
        eq_frac, eq_mass, eq_sa, eq_num = float(OutFrxMassdzbot[begin_eq][-1]), float(DepMassBott[begin_eq][-1]), \
                                          float(SADepBottcomp[begin_eq][-1]), float(NDepBottcomp[begin_eq][-1])
    one_line_result = [ID, Substance, pp_diameter, float(agg_diameter), initial_concentration,
                       simulation_time, eq_time, float(OutFrxMassdzbot[NOutPoints - 1][-1]),
                       float(DepMassBott[NOutPoints - 1][-1]), float(DepMassDissBottcomp[NOutPoints - 1][-1]),
                       float(SADepBottcomp[NOutPoints - 1][-1]), float(NDepBottcomp[NOutPoints - 1][-1]),
                       eq_frac, eq_mass, eq_sa, eq_num]

    print(f'ID {ID}: In vitro simulation complete')
    return one_line_result


# ----------------------END OF MAIN FUNCTION-----------------------#


# -----------------------------------------------------------------#
#
# OUTPUT DATA POINT FUNCTION
#
# -----------------------------------------------------------------#
def calculate_output_point(P, T, Ca, DissC, FrxBottomOccupied, Cabound, subH, ncompBott, BottN, MPAp, ncomp):
    # save the time, in hours
    OutT = T / 3600

    # mass concentrations of each agglomerate species (mg/ml)
    OutCa = Ca

    # total mass agglomerage concentrations
    OutSumCa = np.sum(Ca, axis=1, keepdims=True)

    # mass raw material/particle concentrations for individual agglomerates
    OutCp = Ca / P['Magg']

    # total undissolved mass material/particle concentrations
    OutTotalCp = np.sum(Ca, axis=1, keepdims=True) / P['Magg']

    # dissolved concentration
    OutDissC = DissC

    # fraction of floor occupied
    OutFrxOcc = FrxBottomOccupied

    # mass of material bound per area = c bound * height of compartment
    # (mg/cm^2). multiply concentration bound mg/cm^3 by subH in cm (subH * 100)
    OutMBound = (np.sum(Cabound) / P['Magg']) * subH * 1.0e2

    #
    # MASS AND FRACTION OF MASS
    #

    # CONCENTRATION

    # undissolved particle mass at BottomH per diameter size (shape num BottomH intervals * number of particle sizes)
    OutCpdzbot = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        OutCpdzbot[int(k / BottN)] = np.sum(OutCp[k:k + BottN], axis=0) / BottN

    # average total undissolved mass material/particle concentrations at BottomH
    # intervals (e.g. 10 microns)
    ConcPerBottHcomp = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        ConcPerBottHcomp[int(k / BottN), :] = np.sum(OutTotalCp[k:k + BottN, ]) / BottN

    # total mass material/particle concentrations including dissolved (mg/cm^3 == kg/m^3)
    OutTotalCpd = OutTotalCp + DissC

    # average total mass material/particle concentrations at BottomH
    # intervals (e.g. 10 microns)
    MassPerBottHcomp = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        MassPerBottHcomp[int(k / BottN), :] = np.sum(OutTotalCpd[k:k + BottN, ]) / BottN

    # FRACTION PER COMPARTMENT AND "DEPOSITED"

    # Fraction of material/particle mass in each compartment
    # fx = mass in compartment/total mass in all compartments initially
    # fx = (c * pi*r^2*h)/(c0 * ncompartments * pi*r^2*h)
    OutFrxMass = OutTotalCp / MPAp

    # total Fraction mass material/particle at BottomH
    # intervals (e.g. 10 microns)
    OutFrxMassdzbot = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        OutFrxMassdzbot[int(k / BottN), :] = np.sum(OutFrxMass[k:k + BottN, ])

    # MASS DEPOSITED PER FLOOR AREA

    # total mass of material "deposited" per unit floor area in each
    # compartment (mg/cm^2)
    DepMasscomp = OutTotalCp * subH * 1.0e2

    # total mass "deposited" per unit floor at BottomH
    # intervals (e.g. 10 microns)
    DepMassBott = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        DepMassBott[int(k / BottN), :] = np.sum(DepMasscomp[k:k + BottN, ])

    ## mass deposited per unit floor at BottomH intervals, per particle diameter

    OutMpdzbot = OutCpdzbot * subH * 1.0e2 * BottN

    # total "deposited" mass including dissolved per unit floor in each
    # compartment (mg/cm^2) *mult by 1e2 to convert from kg/m^3 to mg/cm^2
    DepMassDisscomp = OutTotalCpd * subH * 1.0e2

    # total mass "deposited" per unit floor including dissolved at BottomH
    # intervals (e.g. 10 microns)
    DepMassDissBottcomp = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        DepMassDissBottcomp[int(k / BottN), :] = np.sum(DepMassDisscomp[k:k + BottN, ])

        #
    # PARTICLE NUMBER
    #

    # CONCENTRATION

    # number concentration for each agglomerate
    # converted from m^-3 to cm^-3 by dividing by 1e6
    OutNp = OutCp * 3 / (math.pi * 4 * (P['rad'] ** 3) * 1.0e6 * P['pp_density'] * P['PFracAgg'])

    # total number concentration (a bit silly, combining all different sizes, but for consistency
    OutNTp = np.sum(OutNp, axis=1, keepdims=True)

    ## number concentration at BottomH per particle diameter

    OutNpdzbot = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        OutNpdzbot[int(k / BottN)] = np.sum(OutNp[k:k + BottN], axis=0) / BottN

    # average number concentration at BottomH
    # intervals (e.g. 10 microns)
    NBottcomp = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        NBottcomp[int(k / BottN)] = np.sum(OutNTp[k:k + BottN, ]) / BottN

    # DEPOSITED

    # total number concentration of material "deposited" per unit floor area
    # (cm^-2). multiply concentration cm^-3 by subH in cm (subH * 100)
    OutNTpDep = OutNTp * subH * 1.0e2

    # total number concentration "deposited" per floor area at BottomH
    # intervals (e.g. 10 microns)
    NDepBottcomp = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        NDepBottcomp[int(k / BottN), :] = np.sum(OutNTpDep[k:k + BottN, ])

    ## number of particle of different diameter deposited per cm2
    OutNpDepdzbot = OutNpdzbot * subH * 1.0e2 * BottN
    #
    # SURFACE AREA
    #

    # CONCENTRATION

    # surface area concentration of material/particle for each agglomerate
    # converted from (m^2/m^3)to (cm^2/cm^3) by dividing by 10e2
    OutSAp = OutCp * 3 / (P['rad'] * 1.0e2 * P['pp_density'])

    # total surface area concentration of material (cm^2/cm^3)
    OutSATp = np.sum(OutSAp, axis=1, keepdims=True)

    # average surface area concentration at BottomH
    # intervals (e.g. 10 microns)
    OutSATpdzbot = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        OutSATpdzbot[int(k / BottN), :] = np.sum(OutSATp[k:k + BottN, ]) / BottN

    ## surface area deposited at BottomH per particle size
    OutSApdzbot = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        OutSApdzbot[int(k / BottN)] = np.sum(OutSAp[k:k + BottN], axis=0) / BottN

    # DEPOSITED

    # total surface area concentration of material "deposited" per unit floor area
    # (cm^2/cm^2). multiply concentration cm^2/cm^3 by subH in cm (subH * 100)
    OutSATpDep = OutSATp * subH * 1.0e2

    # total surface area concentration "deposited" per floor area at BottomH
    # intervals (e.g. 10 microns)
    SADepBottcomp = np.zeros((int(ncompBott), 1))
    for k in range(0, ncomp, BottN):
        SADepBottcomp[int(k / BottN), :] = np.sum(OutSATpDep[k:k + BottN, ])

    ## surface area deposited per floor area at BottomH per particle size
    OutSADepdzbot = OutSApdzbot * subH * 1.0e2 * BottN

    return (OutT, OutCa, OutSumCa, OutCp, OutTotalCp, OutDissC, OutFrxOcc,
            OutMBound, ConcPerBottHcomp, OutTotalCpd, MassPerBottHcomp,
            OutFrxMass, OutFrxMassdzbot, DepMasscomp, DepMassBott,
            DepMassDisscomp, DepMassDissBottcomp, OutNp, OutNTp, NBottcomp,
            OutNTpDep, NDepBottcomp, OutSAp, OutSATp, OutSATpdzbot,
            OutSATpDep, SADepBottcomp, OutCpdzbot, OutMpdzbot, OutNpdzbot,
            OutNpDepdzbot, OutSApdzbot, OutSADepdzbot)


# -----------------------------------------------------------------#
#
# PLOT TO DATA POINT FUNCTION
#
# -----------------------------------------------------------------#

def plot_output(ConcPerBottHcomp, OutputIntervalHrs, DepMassBott, OutFrxMassdzbot, Output_file,
                SADepBottcomp):
    # get number of timepoints from data
    ntime = len(ConcPerBottHcomp)

    figdepositfinal, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))

    x1 = list(np.arange(float(ntime)) * OutputIntervalHrs)
    y1 = np.array(DepMassBott)[:, -1].flatten()
    # plot bottom mass concentration vs time
    ax1.plot(x1, y1)
    ax1.set_title('Bottom mass per area over time')
    ax1.set_xlabel('Time [h]')
    ax1.set_ylabel('Total particle mass per bottom area [mg/cm2]')

    # plot frx deposited vs time
    y2 = np.array(OutFrxMassdzbot)[:, -1].flatten()
    ax2.plot(x1, y2)
    ax2.set_title('Fraction deposited over time')
    ax2.set_xlabel('Time [h]')
    ax2.set_ylabel('Fraction deposited')

    y3 = np.array(SADepBottcomp)[:, -1].flatten()
    # plot bottom mass concentration vs time
    ax3.plot(x1, y3)
    ax3.set_title('Bottom particle surface area per area over time')
    ax3.set_xlabel('Time [h]')
    ax3.set_ylabel('Total particle surface area per bottom area [cm2/cm2]')
    plt.tight_layout()
    figdepositfinal.savefig('Deposited_' + Output_file + '.pdf')
    plt.close(figdepositfinal)


# -----------------------------------------------------------------#
#
# WRITE OUTPUT DATA FUNCTION
#
# -----------------------------------------------------------------#
def write_output_data(OutT, sticky, OutFrxMassdzbot, DepMassBott,
                      DepMassDissBottcomp, NDepBottcomp,
                      SADepBottcomp, OutDissC, OutFrxOcc, OutMBound, Output_file, eq_time, eq_change_rate):
    print('Writing output data, please wait.')

    # BOTTOM SUMMARY
    name_row = [f'Time at which the marginal increase in deposited fraction per h is < {eq_change_rate} (h)', 'Deposited fraction (of total mass)',
                'Mass deposited per area (mg cm^-2)', 'Mass deposited and dissolved fraction per area (mg cm^-2)', 'Number of particles per area (cm^-2)',
                'Particles surface area per area (cm^2 cm^-2)', 'Concentration of dissolved fraction (mg cm^-3)']
    name_col = [('Time: ' + str(OutT[x]) + ' hours') for x in range(len(OutT))]
    time_row = [0 for x in range(len(OutT))]
    time_row[0] = eq_time

    if sticky:
        name_row += ['% floor occ. (%)', 'mass bound area^-1 (mg cm^-2)']
        matrix = np.concatenate((np.array(time_row)[np.newaxis, :],
                                 np.array(OutFrxMassdzbot)[:, -1].swapaxes(1, 0),
                                 np.array(DepMassBott)[:, -1].swapaxes(1, 0),
                                 np.array(DepMassDissBottcomp)[:, -1].swapaxes(1, 0),
                                 np.array(NDepBottcomp)[:, -1].swapaxes(1, 0),
                                 np.array(SADepBottcomp)[:, -1].swapaxes(1, 0), np.array(OutDissC)[np.newaxis, :],
                                 np.array(OutFrxOcc)[np.newaxis, :], np.array(OutMBound)[np.newaxis, :]), axis=0)
    else:
        matrix = np.concatenate((np.array(time_row)[np.newaxis, :],
                                 np.array(OutFrxMassdzbot)[:, -1].swapaxes(1, 0),
                                 np.array(DepMassBott)[:, -1].swapaxes(1, 0),
                                 np.array(DepMassDissBottcomp)[:, -1].swapaxes(1, 0),
                                 np.array(NDepBottcomp)[:, -1].swapaxes(1, 0),
                                 np.array(SADepBottcomp)[:, -1].swapaxes(1, 0), np.array(OutDissC)[np.newaxis, :]),
                                axis=0)

    file_to_save = pd.DataFrame(matrix, columns=name_col, index=name_row)
    filename = Output_file + '.csv'
    file_to_save.to_csv(filename)
    print('File saved')


# -----------------------------------------------------------------#
#
# DISSOLUTION FUNCTION
#
# -----------------------------------------------------------------#
def perform_dissolution(P, T, C0p, UndissFrx, Ca, SConst, DConst):
    if P['DissRateType'] > 0:
        # if dissolution occurs beyond initial dissolution that occured
        # prior to start of simulation

        # calculate current dissolution fraction over and above
        # P['dissFrx0, the dissolution that has occured since start of
        # simulation
        if P['DissRateType'] == 1:
            # if constant dissolution rate per hour
            # never allow to reach 1, but very close - negligible undissolved)
            dissfrac = P['initial_dissolution'] + (P['DissRate'] * T / 3600.0)
            NewDissFrx = min(0.9999, dissfrac)
        else:
            # dissolution fractions over time (curve) provided by user, allow
            # more complex/realistic dissolution models
            interp = np.interp(xp=P['DissT'], fp=P['TDissFrx'], x=T / 3600.0)

            NewDissFrx = min(0.9999, (P['initial_dissolution'] + interp))

        # New value of undissolved fraction (of what was undissolved after
        # initial dissolution start point)
        NewUndissFrx = 1.0 - NewDissFrx
        # new total dissolved material concentration
        DissC = C0p * NewDissFrx

        # fraction of undissolved at previous iteration remaining
        # undissolved at present iteration or time point = c2/c1
        FrxRemUndiss = NewUndissFrx / UndissFrx
        # new compartment agglomerate concentrations
        Ca *= FrxRemUndiss

        #  concentration change shrinks agglomerates

        # calculate new agglomerate radii
        # v = 4/3*pi*r^3
        # v2/v1 = c2/c1 = r2^3/r1^3
        # r2/r1 = (c2/c1)^(1/3)
        P['rad'] *= FrxRemUndiss ** (1 / 3)

        # calculate new agglomerate species diffusion and
        # sedimentation coefficients

        # agglomerate sed coeff (at C=0)
        P['Scoeff'] = SConst * P['rad'] ** 2  # sec
        # agglomerate diff coeff (at C=0)
        P['Dcoeff'] = DConst / P['rad']

        # set current dissolved and undissolved fraction values to new values
        DissFrx = NewDissFrx
        UndissFrx = NewUndissFrx

    return DissC, Ca, P['rad'], P['Scoeff'], P['Dcoeff'], DissFrx, UndissFrx


# -----------------------------------------------------------------

def update_bottom_bound_by_adsorption(P, Ca, Kd, subH):
    # molar concentration of all agglomerates in bottom compartment
    L = Ca[-1, :] * P['molpermass']
    # sum of molar concentrations
    Lsum = sum(L)
    # fraction occupied theta
    MaxFrxBottomOccupied = Lsum / (Kd + Lsum)

    Frxavail = sum(P['CrossAreaMass'] * Ca[-1, :] * subH)
    # fraction of available material bound
    Frxbound = min(1.0, (MaxFrxBottomOccupied / Frxavail))
    # actual fraction bottom occupied
    FrxBottomOccupied = min(MaxFrxBottomOccupied, Frxavail)
    # Concentrations bound
    Cabound = Ca[-1, :] * Frxbound

    return MaxFrxBottomOccupied, Frxavail, Frxbound, FrxBottomOccupied, Cabound


# ----------------------------------------------------------------

# calculate specific surface area from diameter and density of primary particle
def ssa(diameter, density):
    # radius in cm: nm/0.0000001
    radius = diameter * 0.0000001 / 2
    volume = 4 / 3 * np.pi * radius ** 3
    surface_area = 4 * np.pi * radius ** 2
    # mass in grams
    mass = volume * density
    # ssa in m2/g: cm2/g = 0.0001 m2/g
    ssa = surface_area * 0.0001 / mass
    return ssa


# ---------------------------------------------------------
# Density calculation functions
# calculate agglomerate effective density if missing
def sterling_density(aggl_diam, pp_diam, media_dens, pp_dens):
    if np.isnan(pp_diam):
        raise ValueError('The diameter of the primary particle is missing. '
                         'Impossible to calculate the agglomerate effective density.')
    DF = 2.1
    agglomerate_porosity = 1 - (aggl_diam / pp_diam) ** (DF - 3)
    agg_dens = (1 - agglomerate_porosity) * pp_dens + agglomerate_porosity * media_dens
    return agg_dens


# calculate density of agglomerates in air
def air_density(agg_effective_density, media_density, pp_density):
    return (agg_effective_density - media_density) / (pp_density - media_density) * pp_density


# ---------------------------------------------------------------
# Find location of sample in MPPD database
def find_in_MPPD(MPPD_database, Input_parameters):
    if Input_parameters['air_type'] == 'pp':
        air_diameter = float(Input_parameters['pp_diameter'])
    else:
        air_diameter = float(Input_parameters['agg_diameter'])
    MPPD_line = MPPD_database[(MPPD_database.Substance == Input_parameters['Substance']) & (
            MPPD_database.agg_diameter_nm == air_diameter) & (
                                          MPPD_database.pp_diameter_nm == float(Input_parameters['pp_diameter'])) &
                              (MPPD_database.Exposure_time_h == Input_parameters['simulation_time']) & (
                                      MPPD_database.sex == Input_parameters['sex'])
                              & (abs(
        MPPD_database['density_g/cm2'] - Input_parameters['air_agg_density']) < 0.00000001)]
    MPPD_line_long = MPPD_database[(MPPD_database.Substance == Input_parameters['Substance']) & (
            MPPD_database.agg_diameter_nm == air_diameter) & (
                                               MPPD_database.pp_diameter_nm == float(Input_parameters['pp_diameter'])) &
                                   (MPPD_database.Exposure_time_w == 1820) & (
                                               MPPD_database.sex == Input_parameters['sex']) &
                                   (abs(MPPD_database['density_g/cm2'] - Input_parameters[
                                       'air_agg_density']) < 0.00000001)]
    return MPPD_line, MPPD_line_long


# ------------------------------------------------
# functions for running the MPPD model
def coordinates(x, y):
    global BASE_LEFT, BASE_TOP
    return {
        'x': x / 1.703225806 + BASE_LEFT,
        'y': y / 1.703225806 + BASE_TOP
    }


def clear_input(x, y, val):
    pyautogui.click(x=x, y=y)
    pyautogui.hotkey('ctrl', 'a')
    pyautogui.write(str(val))


######################################################################
###EXECUTING CODE
#######################################################################

start = time.time()

Input_df = pd.read_excel(template_file, sheet_name='Data',
                         skiprows=[0, 1], dtype={'ID': str, 'Output_file': str, 'Substance': str,
                                                 'media_viscosity': float, 'media_density': float,
                                                 'media_temperature': float,
                                                 'pp_density': float, 'agg_diameter': float,
                                                 'agg_effective_density': float, 'column_height': float,
                                                 'initial_concentration': float,
                                                 'simulation_time': int, 'eq_change_rate': float,
                                                 'subcompartment_height': float, 'simulation_time_interval': float,
                                                 'output_time_interval': float,
                                                 'output_compartment_height': float,
                                                 'sedimentation_cdependence': float,
                                                 'diffusion_cdependence': float,
                                                 'initial_dissolution': float, 'dissolution_rate_type': float,
                                                 'dissolution_rate': float,
                                                 'dissolution_times': str, 'dissolution_fractions': str,
                                                 'air_type': str,
                                                 'adsorption_dissociation_constant': float})

lines_to_compute = [x - 4 for x in range(start_value, end_value + 1)]
to_simulate_df = Input_df.loc[lines_to_compute, :].copy()
to_simulate_df.reset_index(drop=True, inplace=True)
required_df = to_simulate_df[['ID', 'Substance', 'sex', 'air_type', 'media_viscosity', 'media_density', 'pp_density',
                              'agg_diameter', 'column_height', 'initial_concentration', 'simulation_time']]
if required_df.isnull().values.any():
    raise ValueError('One or more required parameters are missing. Please check the template input file.')
if (required_df['agg_diameter'] < 26).any():
    raise ValueError('At least one agglomerate diameter is < 26nm. Please remove/modify entry.')
to_simulate_df['air_type'] = to_simulate_df['air_type'].map(lambda x: x.lower())
to_simulate_df['air_type'] = to_simulate_df['air_type'].str.strip()
if not ((to_simulate_df['air_type'] == 'agg') | (to_simulate_df['air_type'] == 'pp')).all():
    raise ValueError('Air type parameter is different from pp or agg. Please check input file.')
# calculate agglomerate density if missing via Sterling equation
to_simulate_df['agg_effective_density'] = to_simulate_df.apply(
    lambda x: sterling_density(x['agg_diameter'], x['pp_diameter'], x['media_density'], x['pp_density']) if np.isnan(
        x['agg_effective_density']) else x['agg_effective_density'], axis=1)
to_simulate_df['air_agg_density'] = to_simulate_df.apply(
    lambda x: air_density(x['agg_effective_density'], x['media_density'], x['pp_density']) if x['air_type'] == 'agg' else x['pp_density'], axis=1)
to_simulate_df['sex'] = to_simulate_df['sex'].map(lambda x: x.lower())
to_simulate_df['sex'] = to_simulate_df['sex'].str.strip()
if not ((to_simulate_df['sex'] == 'male') | (to_simulate_df['sex'] == 'female')).all():
    raise ValueError('Sex parameter is different from male or female. Please check input file.')

###### UNIQUE OUTPUT FILE
final_df_columns = ['ID', 'Substance', 'Primary_particle_diameter_nm', 'Agglomerate_diameter_nm',
                    'Initial concentration_mg/cm3', 'End_time_h', 'Stable_time_h', 'Fraction_deposited_end',
                    'Mass_per_Area_end_mg/cm2', 'Mass+Dissolved_per_Area_end_mg/cm2', 'SA_per_Area_end_cm2/cm2',
                    'N_per_Area_end_n/cm2', 'Fraction_deposited_eq', 'Mass_per_Area_eq_mg/cm2',
                    'SA_per_Area_eq_cm2/cm2',
                    'N_per_Area_eq_n/cm2', 'Agglomerate_density_air', 'Sticky_bottom',
                    'SA_per_area_end_pp_cm2/cm2', 'sex', 'air_particle_type', 'exposure_concentration_air_end_mg/m3',
                    'exposure_concentration_air_eq_mg/m3', 'exposure_concentration_air_year_mg/m3',
                    'exposure_concentration_air_worklife_mg/m3']

final_file = []

Input_dict = to_simulate_df.to_dict(orient='records')
# Check if the cases evaluated are already present in the MPPD database, otherwise calculates them and saves the file
MPPD = pd.read_csv(MPPD_dataset)

to_compute_MPPD = []
for line in range(len(lines_to_compute)):
    Input_parameters = Input_dict[line]
    MPPD_line, MPPD_line_long = find_in_MPPD(MPPD, Input_parameters)
    if MPPD_line.empty or MPPD_line_long.empty:
        to_compute_MPPD.append(line)

if to_compute_MPPD:
    MPPD_list = []
    MPPD_data = []
    columns_MPPD = ['Substance', 'pp_diameter_nm', 'agg_diameter_nm', 'density_g/cm2', 'Lung_Model', 'sex',
                    'Exposure', 'Clearance', 'Exposure_time_h', 'Exposure_time_d', 'Exposure_time_w',
                    'Alveolar_retention_h_mg',
                    'Alveolar_retention_y_mg', 'Alveolar_retention_w_mg', 'Alveolar_surface_cm2', 'Deposited_per_area_h_mg/cm2',
                    'Deposited_per_area_y_mg/cm2',
                    'Deposited_per_area_w_mg/cm2']
    list_filenames = []
    # pulmonary parameters woman from Brown 2013
    URT = 40
    TV = 464
    BF = 14
    FRC = 2680
    # exposure parameters
    days = [1, 5]
    weeks = [1, 1820]
    scenario_label = ['s', 'w']

    # size parameters
    diameter_label = ['agg', 'pp']

    print(
        'Not all particles are available in the MPPD database. They will be calculated now. \nPlease open the MPPD program,'
        'and make sure the resolution of your main screen is 2560 x 1440 (100% scaling). Keep the MPPD program on your main screen. \n'
        'To run the analysis using the MPPD software, you will need to have the MPPD software window on screen (not covered by any other program/window).'
        'You should not use or move the mouse, as the MPPD will be run by automatically moving the mouse to run the program multiple times. \n'
        'When you are ready, press ENTER; you will have 2 seconds to have the MPPD software window open and visible on screen (if not already), afterwards please do not move the mouse anymore.')
    input('Please, press ENTER to start running the calculation of MPPD results')

    time.sleep(2)
    icon = pyautogui.locateOnScreen(MPPD_icon)
    if not icon:
         print('MPPD window not located on screen, please make sure that the resolution is 2560 x 1440 (100% scaling),')
         input('If all points mentioned before have been fixed, press ENTER to retry. You will have 2 seconds to have '
             'the MPPD software window on screen and to stop moving the mouse')
         time.sleep(2)
         icon = pyautogui.locateOnScreen(MPPD_icon)
         if not icon:
             print('Impossible to find the MPPD window. The model will stop running.')
             sys.exit()

    BASE_LEFT, BASE_TOP = icon.left, icon.top

    for line in to_compute_MPPD:
        density = to_simulate_df.loc[line, 'air_agg_density']
        diameter = [to_simulate_df.loc[line, 'agg_diameter'] / 1000, to_simulate_df.loc[line, 'pp_diameter'] / 1000]
        hours = [to_simulate_df.loc[line, 'simulation_time'], 8]
        ID = int(to_simulate_df.loc[line, 'ID'])
        Substance = to_simulate_df.loc[line, 'Substance']

        # Begins autogui
        pyautogui.click(**coordinates(22, 75))
        pyautogui.click()

        # selecting airway morphometry
        pyautogui.click(**coordinates(107, 76))
        time.sleep(0.5)
        pyautogui.moveTo(**coordinates(107, 114))
        pyautogui.click()

        # insert values deposition
        if to_simulate_df.loc[line, 'sex'] == 'female':
            time.sleep(0.5)
            clear_input(**coordinates(372, 496), val=FRC)
            clear_input(**coordinates(372, 596), val=URT)
        elif to_simulate_df.loc[line, 'sex'] == 'male':
            time.sleep(0.5)
            clear_input(**coordinates(372, 496), val=3300)
            clear_input(**coordinates(372, 596), val=50)
        time.sleep(0.5)
        pyautogui.click(**coordinates(726, 701))
        #
        if to_simulate_df.loc[line, 'air_type'] == 'agg':
            size = 0
        elif to_simulate_df.loc[line, 'air_type'] == 'pp':
            size = 1
        else:
            raise ValueError('The type of particle in air should be either pp or agg. Please check the template file.')

        # selecting inhalant properties- aerosol
        pyautogui.click(**coordinates(107, 76))
        pyautogui.moveTo(**coordinates(107, 156))
        pyautogui.click()
        pyautogui.moveTo(**coordinates(422, 156), duration=0.15)
        pyautogui.click()
        # insert values in aerosol window
        clear_input(**coordinates(422, 256), val=density)
        clear_input(**coordinates(422, 406), val=diameter[size])
        pyautogui.click(**coordinates(792, 1066))
        # selecting exposure conditions
        pyautogui.click(**coordinates(107, 76))
        pyautogui.moveTo(**coordinates(107, 196))
        pyautogui.moveTo(**coordinates(422, 196), duration=0.15)
        pyautogui.click()
        # insert values exposure conditions
        if to_simulate_df.loc[line, 'sex'] == 'female':
            clear_input(**coordinates(552, 666), val=BF)
            clear_input(**coordinates(552, 726), val=TV)
        elif to_simulate_df.loc[line, 'sex'] == 'male':
            clear_input(**coordinates(552, 666), val=12)
            clear_input(**coordinates(552, 726), val=625)
        pyautogui.click(**coordinates(792, 1006))

        for exposure in [0, 1]:
            scenario = scenario_label[exposure]
            diameter_used = diameter_label[size]
            if to_simulate_df.loc[line, 'sex'] == 'female':
                sex_label = 'f'
            else:
                sex_label = 'm'
            filename = f'{ID}{Substance}_{diameter[size]}um_{diameter_used}_{scenario}_{sex_label}'
            print(filename)
            list_filenames.append(filename)

            # selecting clearance
            pyautogui.click(**coordinates(107, 76))
            pyautogui.moveTo(**coordinates(107, 236))
            pyautogui.moveTo(**coordinates(422, 236), duration=0.15)
            pyautogui.moveTo(**coordinates(422, 276), duration=0.15)
            pyautogui.click()
            # insert values deposition
            clear_input(**coordinates(622, 676), val=hours[exposure])
            clear_input(**coordinates(622, 716), val=days[exposure])
            clear_input(**coordinates(622, 756), val=weeks[exposure])
            pyautogui.click(**coordinates(782, 841))

            # accept settings and run
            pyautogui.click(**coordinates(236, 76))
            pyautogui.click(**coordinates(236, 104))
            pyautogui.click(**coordinates(236, 76))
            pyautogui.click(**coordinates(236, 146))
            time.sleep(4)
            pyautogui.click(**coordinates(467, 526))

            # create report file
            pyautogui.click(**coordinates(392, 76))
            pyautogui.click(**coordinates(392, 104))
            time.sleep(0.5)
            clear_input(**coordinates(822, 586), val=filename)
            time.sleep(0.5)
            pyautogui.click(x=1400, y=873)
            time.sleep(0.5)

            MPPD_list.append([Substance, float(diameter[1] * 1000), float(diameter[size] * 1000), float(density),
                              'Yeh/Schum Symmetric', Input_df.loc[line, 'sex'], 'constant', 'Yes', int(hours[exposure]),
                              int(days[exposure]),
                              int(weeks[exposure])])
    print('MPPD results obtained: you can now close the MPPD program and use the mouse again.')
    for count in range(len(list_filenames)):
        parameters = MPPD_list[count]
        file = list_filenames[count]
        result = pd.read_csv(f'{MPPD_results}//{file}.rpt', skiprows=[0,1], encoding='ANSI', sep='\t', skipinitialspace=True,
                             names=['hour', 'day', 'week', 'minute', 'TB_cleared', 'TB_retention', 'TB_retention2',
                                    'Alveolar_cleared', 'empty', 'Alveolar_retention', 'Alveolar_retention_2'])
        # find the line at which the clearance results start
        start_line = result.index[result.hour == 'Hour'].tolist()[0]
        # results are split on two columns, fill one with all of them
        mask = np.isnan(result['Alveolar_retention'])
        result['Alveolar_retention'] = np.where(mask, result['Alveolar_retention_2'], result['Alveolar_retention'])
        if parameters[5] == 'male':
            alv_sa = 792000
        else:
            alv_sa = 559000
        # find results of alveolar retention for either exposure time or year and worklife time
        if parameters[10] == 1820:
            Alv_retention_y = result['Alveolar_retention'][
                (result.index > start_line) & (result.hour == str(parameters[8]))
                & (result.day == str(parameters[9])) & (
                        result.week == '52')]
            Alv_retention_w = result['Alveolar_retention'][
                (result.index > start_line) & (result.hour == str(parameters[8]))
                & (result.day == str(parameters[9])) & (result.week == str(parameters[10]))]
            alv_dep_sa_y = float(Alv_retention_y) / alv_sa
            alv_dep_sa_w = float(Alv_retention_w) / alv_sa
            parameters.extend(
                [np.nan, float(Alv_retention_y), float(Alv_retention_w), alv_sa, np.nan, alv_dep_sa_y, alv_dep_sa_w])
        else:
            Alv_retention = result['Alveolar_retention'][
                (result.index > start_line) & (result.hour == str(parameters[8]))
                & (result.day == str(parameters[9])) & (
                        result.week == str(parameters[10]))]
            alv_dep_sa = float(Alv_retention) / alv_sa
            parameters.extend([float(Alv_retention), np.nan, np.nan, alv_sa, alv_dep_sa, np.nan, np.nan])
        empty_items = 18 - len(parameters)
        if empty_items < 0:
            raise ValueError('more parameters than column names')
        MPPD_data.append(parameters)
    MPPD_df = pd.DataFrame(data=MPPD_data, columns=columns_MPPD)
    MPPD_df.drop_duplicates(inplace=True)
    MPPD_updated = pd.concat([MPPD, MPPD_df])
    MPPD_updated.to_csv(MPPD_dataset, index=False)

too_small = []
# calculate the deposited amounts in vitro and the corresponding exposure levels in vivo
for line in range(len(lines_to_compute)):
    Input_parameters = Input_dict[line]
    if Input_parameters['dissolution_rate_type'] == 2:
        Input_parameters['dissolution_fractions'] = np.asarray(
            [float(x) for x in Input_parameters['dissolution_fractions'].split(',')])
        Input_parameters['dissolution_times'] = np.asarray(
            [float(x) for x in Input_parameters['dissolution_times'].split(',')])
    # load MPPD database that has been updated with missing samples
    MPPD = pd.read_csv(MPPD_dataset)
    # identify if MPPD database includes the data for the specific substance
    if Input_parameters['air_type'] == 'pp':
        air_diameter = float(Input_parameters['pp_diameter'])
    else:
        air_diameter = float(Input_parameters['agg_diameter'])
    MPPD_line, MPPD_line_long = find_in_MPPD(MPPD, Input_parameters)
    MPPD_line = MPPD_line.to_dict(orient='records')[0]
    MPPD_line_long = MPPD_line_long.to_dict(orient='records')[0]
    if pd.isnull(Input_parameters['sticky']):
        Input_parameters['sticky'] = 'False'
    Input_parameters = {k: v for k, v in Input_parameters.items() if type(v) != float or (type(v) == float and v == v)}

    # Run in vitro simulation
    line_for_df = in_vitro_simulation(**Input_parameters)

    # if there is a result for "equilibrium" time
    if isinstance(line_for_df[6], float):
        # consider the increase per hour of deposited dose given by the difference between the yearly deposited and
        # the deposited at simulation time; this is an approximation assuming that the increase in deposition is linear
        hourly_add_dep = (MPPD_line_long['Deposited_per_area_y_mg/cm2'] - MPPD_line['Deposited_per_area_h_mg/cm2']) / (
                8760 - Input_parameters['simulation_time'])
        dep_area_eq = MPPD_line['Deposited_per_area_h_mg/cm2'] + (
                line_for_df[6] - Input_parameters['simulation_time']) * hourly_add_dep
        exposure_eq = float(line_for_df[13] / dep_area_eq)
    else:
        exposure_eq = 'not defined'
    exposure_end = float(line_for_df[8] / MPPD_line['Deposited_per_area_h_mg/cm2'])
    exposure_year = float(line_for_df[8] / MPPD_line_long['Deposited_per_area_y_mg/cm2'])
    exposure_worklife = float(line_for_df[8] / MPPD_line_long['Deposited_per_area_w_mg/cm2'])
    if 'ssa' not in Input_parameters.keys():
        Input_parameters['ssa'] = ssa(Input_parameters['pp_diameter'], Input_parameters['pp_density'])
    SA_dep_ssa = Input_parameters['ssa'] * line_for_df[8]
    line_for_df.extend([Input_parameters['air_agg_density'], str(Input_parameters['sticky']), SA_dep_ssa, Input_parameters['sex'],
                        Input_parameters['air_type'], exposure_end, exposure_eq, exposure_year, exposure_worklife])

    final_file.append(line_for_df)

final_df = pd.DataFrame(data=final_file, columns=final_df_columns)
now = time.strftime("%Y-%m-%d_%H_%M")
name = f'CoDo_results_{now}.csv'
final_df.to_csv(name, index=False)
end = time.time()
print(f'Task completed. Time for simulations: {round(end - start)} seconds')
