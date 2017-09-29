#
#   GEOCARB
#
# This code is a python translation of David Archer's version of Robert Berner's GEOCARB III model
# Original code is available at http://forecast.uchicago.edu/Projects/geocarb.doc.html and
# http://forecast.uchicago.edu/Projects/geocarb_web.F
#
# I have not added anything original. Merely translated it from FORTRAN to Python,
# cleaned it up a little, and started adding detailed comments to explain how it works.
#
# I also noted what seem to be typos in the temperature dependence of the weathering
# feedback in lines 106-115 of Archer's code (calculating f_bt_co2).
#
# Jonathan Gilligan <jonthan.gilligan@vanderbilt.edu>
#
#
import math, numpy
import os, sys
import string, time
import pandas

trace = False
trace_interval = 100000

# Definition of GEOG (Berner 2004, pp. 35-36; Berner & Kothavala 2001, pp. 184, 193)
#
# GEOG = temperature difference from today due to meridional continental drift.
#
# For formulation, see Otto-Bleisner 1995, p. 11,544, "For the early Phanerozoic,
# a 1.5 degree decrease in mean latitude is equivalent to a 1 degree C increase
# in continental temperature."
#
# Today, mean continental latitude is 30N, so zero occurs at 30N or 30S
#
# It's not clear where the numerical factors 0.5 and 10 come from.
# Reading Otto-Bleisner and fitting her data to a linear model
# T = a cos(lat) + b, I get a =  18.6 and b = -1.9 when I fit all the data
# or a = 57.1 and b = -31.8 when I omit two outliers
# (T=0 is radically different from earlier times).
# This formula seems to be using a = 10.0 and b = -8.0, which don't match
# Otto-Bleisner's data at all.
#
#
def calc_geog(mean_latitude):
    return (math.cos(mean_latitude / 180.0 * math.pi) - 0.5) * 10.0 - 3.0

#
# Function NEWT
#
# Calculate partitioning of carbon into CO2, dissolved CO2/carbonic acid, carbonate, and bicarbonate.
#
# This follows the method in the appendix of Peng 1987
#
# alk = total alkalinity, in micromoles per kg
# tco2 = total carbon in the ocean, in micromoles per kg
#       Peng Eq. A1: tco2 = solubility * pco2 + [hco3-] + [co3-2]
#       solubility * pco2 = saturated dissolved co2.
#
# sal = salinity in ppt.
# temp = ocean temperature in celsius
# k1 and k2 are the dissociation constants of carbonic acid and bicarbonate.
#
def newt(alk, tco2, sal, temp, k1, k2):
    # coefficient for dissociation of boric acid
    #
    # Peng's version cites Lyman 1958, giving
    # log10(kb) = -9.26 + 0.00886 * sal + (tkt - 2.73)
    kb=1.E-9
    # total silicon
    tsi = 0.
    tp = 0.
    # total boric acid,
    # Peng et al. take this from Culkin 1965.
    # See also, Millero 1995, eq. 53
    tbor = 4.106E-4 * sal / 35.
    # Kelvin temp
    tkt = temp + 273.
    # Kelvin temp / 100
    tk = tkt / 100.

    # coefficient of dissociation of silicic acid
    #
    # Peng gives ksi = 4E-10 and assumes it's constant wrt salinity and temp.
    # "Since the alkalinity due to silicate dissociation in seawater is commonly
    #  less than 0.2% of the total alkalinity ... the error caused by this
    #  assumption is very small." (p. 454).
    ksi = 1.E-10
    # coefficients of dissociation of phosphoric acid
    #
    # kp1 ~ 2E-2, so it can be neglected.
    # expressions for kp2 and kp3 are fitted by Peng et al. to data from
    # Kester 1967
    kp2 = math.exp( -9.039 - 1450./tkt)
    kp3 = math.exp(  4.466 - 7276./tkt)
    # dissociation constant for water.
    # Peng et al. take this from Culberson  1973
    # See also, Millero 1995, Eq. 63
    kw  = math.exp( 148.9802 - 13847.26/tkt -23.6521 * math.log(tkt) - 0.019813 * sal \
          + math.sqrt(sal) * (-79.2447 + 3298.72/tkt + 12.0408 * math.log(tkt)))
    # Total activity coefficientfor hydrogen
    # Peng et al. take this from Culberson 1973 and Takahashi 1982a
    #
    # [H+] = ah/fh
    fh = 1.29 -0.00204 * tkt + 4.6E-4 * sal**2 - 1.48 * 1.E-6 * sal**2 * tkt

    c1 = k1/2.0
    c2 = 1.0 - 4.0 * k2/k1
    c4 = tbor * kb

    # Total activity for hydrogen.
    # Initial guess, which will be refined in the loop below
    #
    # aht = (1/2)(k1/ac)((tco2 - ac)
    #        + sqrt( (tco2 - ac)^2 + 4 (ac * k2/k1) * (2 * tco2 - ac) ),
    # where ac = carbonate alkalinity.
    # Eq. (A10)
    #
    # Total alkalinity is known (fixed), and ac depends on aht, so we iterate
    # solving for ac as a function of aht, and then solving for aht as a function
    # of the new ac, and repeating.
    #
    # aht = test value for ah
    aht= 1.E-8

    for icnt in range(100):
        #
        # bm = borate alkalinity = [H2BO3-]
        # Eq. (A5)
        bm = tbor * kb / (aht + kb)
        #
        # sim = silicate alkalinity = [H3SiO4-]
        # Eq. (A6)
        sim = tsi * 4E-10 / (aht + 4E-10)
        # pm = phosphate alkalinity = [H2PO4-] + 2 [HPO4-2] + 3 [PO4-3]
        # Eq. (A7)
        pm = tp * (1 / (1 + kp2/aht + kp2*kp3 / aht**2) \
                 + 2 / (1 + aht/kp2 + kp3/aht) \
                 + 3 / (1 + aht/kp3 + aht**2 / (kp3*kp2)) )
        #
        # wm = water alkalinity = [OH-] - [H+]
        # Eq. (A8)
        wm = ( (kw*fh/aht) - (aht/fh) )

        # a = carbonate alkalinity
        #
        # Total alkalinity is the sum of the alkalinity due to each species.
        # carbonate, borate, silicate, phosphate, and water,
        # so carbonate alkalinity = total alk. - alk. due to other species.
        a = alk - bm - sim - pm - wm

        x = a/tco2
        #
        # ah1 = new value of ah from Eq. (A10)
        #
        ah1 = c1/x * (1.0 - x + math.sqrt(1.0 + c2 * x * (x - 2.0)))
        # did it converge?
        if math.fabs(1.0 -aht/ah1) > 5E-5:
            # if it hasn't converged, the new test value is ah1
            aht=ah1
        else:
            # pass
            break
            #when we're done, ah1 is our best estimate for ah.

    # Once we know ah and ac we can partition inorganic carbon:
    #
    # a = hco3 + 2 co3
    # tco2 = co2 + hco3 + co3
    # (a-tco2) = co3 - co2
    #
    # k1 = ah hco3/co2 (note: there's a typo on p. 454 of Peng)
    # k2 = ah co3/hco3
    #
    # k1*k2 = ah^2 * co3/co2
    # (a-tco2) / (1 - ah^2/k1k2) = (co3 - co2) / (1 - co2/co3)
    #       = co3 * (1 - co2/co3) /(1 - co2/co3) = co3
    co3 = (a-tco2)/(1.0 - (ah1*ah1)/(k1*k2))
    #
    # ah/k1 = co2/hco3
    #
    # tco2 / (1 + ah/k1 + k2/ah) = tco2 / (1 + co2/hco3 + co3/hco3)
    #   = hco3 * (co2 + hco3 + co3) / (hco3 + co2 + co3) = hco3
    #
    hco3 = tco2/(1.0 + ah1/k1 + k2/ah1)
    #
    # tco2 = co2 * (1+k1/ah + k1*k2/ah^2) (Eq. A9)
    #
    # invert the equation to get saturated pco2:
    #
    # co2 = tco2 / (1 + k1/ah + k1*k2/ah^2)
    co2 = tco2/(1.0 + k1/ah1 + k1*k2/(ah1*ah1))

    # Units: micromoles per kilogram seawater.
    #
    return co2, hco3, co3

#
# Calculate pk1, pk2: log of rate constants. Numbers from
# C. Mehrbach et al., "Measurement of the Apparent Dissociation
# Constants of Carbonic Acid in Seawater at Atmospheric Pressure"
# Limnology & Oceanography 18, 897 (1973),
# as cited by Peng 1979
#
def logk1(tk,sal):  # Mehrbach Table 1, eq. 10
    lk1 = 13.7201 - 0.031334 * tk - 3235.76 / tk - 1.3E-5 * sal * tk + 0.1032 * math.sqrt(sal)
    return lk1

def logk2(tk, sal): # Mehrbach Table 1, eq. 11
    lk2 = -5371.9645 -1.671221 * tk + 128375.28 / tk + 2194.3055 * math.log10( tk ) - 0.22913 * sal \
          - 18.3802 * math.log10( sal ) + 8.0944E-4 * sal * tk  + 5617.11 * math.log10( sal ) / tk \
          - 2.136 * sal / tk
    return lk2

#
# Calculate partitioning of inorganic carbon in the oceans.
# Follow the methods (and numbers) of Peng 1979.
#
#
# conc is a 2-element vector.
#
# conc[0] is total dissolved inorganic carbon
# conc[1] is ocean alkalinity
# Both are in units of moles per cubic meter.
#
# temp is the ocean temperature, in celsius.
#
def calc_co2(conc, temp):
    # salinity in ppt.
    sal = 35.0
    tk = temp + 273.15

    k1 = 10**logk1(tk,sal)
    k2 = 10**logk2(tk,sal)

    # Convert conc from moles per cubic meter to micromoles per kg
    co2, hco3, co3 = newt(conc[0]/1000., conc[1]/1000., sal, temp, k1, k2)

    # Henry's law coefficient.
    # Weiss, 1974 Eq. 12, with coefficients in moles/kg*atm from Table 1,
    # as cited by Peng 1979.
    # See also, Millero 1995, Eq. 26.
    #
    # Note that Peng calls khco2 "alpha_s".
    #
    # Units mole/kg/atm.
    khco2 = math.exp( -60.2409 + 9345.17 / tk + 23.3585 * math.log(tk/100.)
                + sal * (0.023517 - 2.3656e-4 * tk + 4.7036e-7 * tk * tk) )

    # Convert from Peng's units (micromole per kg and microatmosphere)
    # to mole per cubic meter and atmosphere.
    pco2 = co2 / khco2 * 1.e6
    co3 = co3 * 1.e6

    return pco2, co3

#
# This is the main function.
#
# GEOL_YEAR = year before present (used to set solar constant
# MEAN_LATITUDE = mean latitude of continents. Used to determine runoff parameters from mean solar zenith
# DT2X = climate sensitivity (celsius warming for doubled CO2)
# B_PLANTS = presence vascular plants
# F_LAND = fraction of land (relative to today)
# TIMESTEP = years per timestep.
#
# Basic code is from Berner and Kothavala, GEOCARB III. Archer added code to determine
# partitioning of inorganic carbon in the oceans.
#
#
def geocarb(co2_slug = (0, 1000, 0), co2_degas = 7.5,
            timestep = (50.,1., 50.),
            duration = (5E6, 100, 2E6),
            dt2x  = 3.0,
            geol_year = 0.0, mean_latitude = 3.0,
            b_plants = True, f_land = 1,
            lead_in = 0,
            co2_ts = None):
    global trace
    global trace_interval

    np = len(duration) # number of periods. Default = 2: a 5-million year spinup and then a 2 million year experimental run.

    if type(co2_slug) in (int, float):
        co2_slug = [co2_slug] * np
    if type(co2_degas) in (int, float):
        co2_degas = [co2_degas] * np
    if type(timestep) in (int, float):
        timestep = [timestep] * np
    if type(duration) in (int, float):
        duration = [duration] * np
    if type(geol_year) in (int, float):
        geol_year = [geol_year] * np
    if type(mean_latitude) in (int, float):
        mean_latitude = [mean_latitude] * np
    if type(b_plants) in (bool, int, float):
        b_plants = [bool(b_plants)] * np
    if type(f_land) in (int, float):
        f_land = [f_land] * np

    # Slug of CO2 released at the end of a period. For the default 2-period run, co2_slug = (x,0),
    # where x is the CO2 release in gigatons.
    #
    co2_slug = list(map(lambda x : x * 1E+3/12., co2_slug)) # Convert GTon C to 10^12 moles

    n_steps = list(map(lambda i : int(duration[i] / timestep[i]), range(np)))

    if co2_ts is None:
        co2_ts = list(map(lambda n : [0.] * n, n_steps))
    co2_ts = list(co2_ts)
    for i in range(np):
        if type(co2_ts[i]) in (int, float):
            co2_ts[i] = [co2_ts[i]] * n_steps[i]
    co2_ts = list(map(lambda x : numpy.array(x) * 1E+3/12., co2_ts)) # Convert GTon C to 10^12 moles


    # Carbonate weathering constant, from Berner 1991 (GEOCARB I), p. 344
    weat_carb_k = 0.00261 * 3000
    # I'm not clear where this constant comes from.
    weat_sil_k = 7.0
    #
    # dt2x = 3.0
    # time_step = 50.0    # years
    #
    MOLpPPM = 2.0E+3 / 12.0 # conversion between 10^12 moles carbon and ppm in atmosphere
    ocn_vol = 3.0E+2 * 4000 # ocean volume (units of 10^12 cubic meters)
    MOLpMM3 = ocn_vol   # conversion between ocean inventory in 10^12 moles and moles per cubic meter
    GEpPPM = 1.0/300.0 * MOLpPPM # gas exchange in 10^12 moles to ppm in atmosphere

    #
    # solar_response = solar response factor Ws from Berner 2004 (p. 31)
    # and Berner & Kothavala 2001 (p. 184, Eq. 3 and pp. 191-92)
    #
    solar_response = -7.4
    #
    # river runoff parameter: transport of dissolved weathering products to ocean
    # Berner + Kothavala 2001, p. 184, Eq. 2 and pp. 191-92.
    #
    run = 0.045 # Bermer & Kothavala 2001, Fig 3, cold periods
    #
    # tau  = gamma from Berner 1994 Eq. 32, Berner & Kothavala 2001, (Eq. 3)
    #
    # z = ACT in Berner & Kothavala 2001, Eq. 2.
    #
    # See also Eq. 29, 31, 34 in Berner 1994.
    # z = Delta E / (R T T0), where
    #       Delta E = activation energy for silicate dissolution
    #               = 15,000 cal/mole
    #       R = gas constant = 1.99 calorie/mole-kelvin
    #       T, T0 = current and "original" temperature in Kelvin.
    #       We approximate T = T0 = constant = 288K
    #
    # So z = ACT = 0.09
    #
    # original: tau = dt2x, but I think it should be this:
    tau = dt2x / math.log(2)
    z = 0.09
    #
    # CO2 fertilization effect on vascular plants: B&K, Eq. 5,
    # and the discussion on pp. 188-90. See also, Berner 1994
    # Eq. 34-37
    #
    fert = 0.4

    # run_steps = int(run_time / timestep)
    # leadin_index = spinup_steps - int(1000 / timestep)
    leadin_index = int(lead_in / timestep[0])
    data_length = sum(n_steps) - leadin_index

    i_alk = 0
    i_tc = 1

    # initial conditions
    pco2 = 200.0    # ppm
    atm_temp_0 = 15.0 # celsius
    ocn_temp_0 = 4.0
    ocn_temp = ocn_temp_0

    # ocn_conc = [ total alkalinity, total carbon ], where total alkalinity includes all species
    # (water dissociation, phosphate, silicate and borate as well as carbonate.
    #
    # alk = [OH-] - [H+] + [H2BO3-] + [H3SiO4-] + [H2PO4-] + 2 [HPO4-2] + 3 [PO4-3]
    #        + [HCO3-] + 2 [CO3-2]
    #
    # tco2 = co2 + hco3- + co3-2
    #
    # See Peng 79, eq (A5)-(A8)
    #
    ocn_conc = [2.3,0.0]   # mol/m^3
    ocn_tco2_top = 5.0     # two-box model of ocean: deep and mixed layer.
    ocn_tco2_bot = 1.0     # tco2 in moles per cubic meter.

    ocn_inv = [0.0,0.0]

    for it in range(10):
        #
        # Mixing deep and shallow ocean
        #
        ocn_conc[i_tc] = (ocn_tco2_top + ocn_tco2_bot)/2
        pco2_test, co3 = calc_co2(ocn_conc, ocn_temp_0)
        if pco2_test > pco2:
            ocn_tco2_top = ocn_conc[i_tc]
        else:
            ocn_tco2_bot = ocn_conc[i_tc]

    year = timestep[0] - duration[0]
    data_names = ("year", "tco2", "pCO2", "alk", "CO3", "WeatC", "WeatS",
                  "TotW", "BurC", "Degas", "Emiss", "Tatm", "Tocn")
    data_cols = len(data_names)
    data = numpy.zeros((data_length, data_cols),numpy.double)
    data_index = 0

    for period in range(np):
        time_step = timestep[period]
        nsteps = int(duration[period] / time_step)
        co2_series = co2_ts[period]
        degas = co2_degas[period]
        gyear = geol_year[period]
        geog = calc_geog(mean_latitude[period])
        has_plants = b_plants[period]
        frac_land = f_land[period]
        last_year = year
        delta_time = time_step
        dt = last_year % time_step
        if dt > 0:
          delta_time -= dt
          dt = delta_time
          nsteps += 1
          co2_series = numpy.insert(co2_series, 0, co2_series[0])
          data = numpy.append(data, numpy.zeros((1,data_cols)), 0)
        print("nsteps = ", nsteps, ", len(series) = ", len(co2_series),
              ", last_year = ", last_year, ", step = ", time_step, ", dt =", dt)
        for itime in range(nsteps):
            year = last_year + dt + itime * time_step
            co2_emissions = co2_series[itime]
            # r_co2 = ratio of atmos. CO2 to late preindustrial holocene
            r_co2 = pco2 / 280.
            co2_doublings = math.log(r_co2,2)
            # air temperature adjusts instantly.
            #
            # I wonder what would happen if we put in a delay time on this
            # to take thermal inertia into account.
            atm_temp = atm_temp_0 + dt2x * co2_doublings
            ocn_temp_target = atm_temp - 10.
            # ocean temperaeture has 1000-year time constant
            ocn_temp = ocn_temp + (ocn_temp_target - ocn_temp) / 1000.0 * delta_time
            #
            # Basic strategy is laid out in Berner & Kothavala 2001 p. 184:
            # Eq(1): weathering feedback f_bt_co2(T,CO2) = f_t(t) f_co2(co2)
            # Eq(2): f_t(t) = exp(z Delta_t) (1 + run * Delta_t)^0.65
            # Eq(3): Delta_t = tau * co2_doublings - solar_response * (geol_year/570) + geog
            # Eq(4): f_co2(co2) = sqrt(r_co2) for prevasculare plant weathering
            # Eq(5): f_co2(co2) = (2 r_co2 / (1 + r_co2))^fert for weathering affected
            #                                                   by vascular plants
            #
            delta_t = co2_doublings * dt2x + geog + solar_response * (gyear / 570.)
            #
            # Original version:
            # delta_t = r_co2**tau + solar_response(geol_year / 579) + geog
            #
            #
            # Basic weathering: CO2 thermostat. Berner 1994 Eq. 33,
            # Substitute Eq. (3) into the first part of Eq(2)
            #
            # This part is f_t( T(CO2, geol_year, geog) )
            #
            f_bt_co2 = math.exp(delta_t *z) * (1 + run  * delta_t )**0.65
            #
            #
            # Here we calculate f_co2(co2).
            # Archer's version is:
            #    with plants, f_co2(co2) = 1
            #    without plants, f_co2(co2) = 0.8
            # Vascular plants accelerate weathering
            #
            # Berner and Kothavala have carbon dioxide fertilization effect for vascular plants
            # See B&K 2001 Eq. 4,5 and the discussion on pp. 188-90.
            # See also, Berner 1994 Eq. 34-37
            #
            if True:
                # Archer's version.
                if not has_plants:
                  f_bt_co2 *= 0.8
            else:
                # Berner and Kothavala version
                if not has_plants:
                    f_bt_co2 *= math.sqrt(r_co2)
                else:
                    f_bt_co2 *= (2 * r_co2 / (1. + r_co2))**fert
            #
            # Carbonate and silicate weathering
            #
            weat_carb = weat_carb_k * f_bt_co2 * frac_land # 10^12 mol/yr
            weat_sil = weat_sil_k * f_bt_co2 * frac_land
#
#   GAS EXCHANGE
#
            pco2_sat, co3 = calc_co2(ocn_conc, atm_temp)
            #
            # Gas exchange: atmosphere saturates with
            # CO2 relative to dissolved vapor pressure.
            #
            ge = GEpPPM * (pco2_sat - pco2) # + up
#
#   CaCO3 precipitation/burial
#
            bc = 3E-3 * (co3-150.) * 1E3/12. # Archer 1997 and Archer et al., 2009: 150 is offset between ocn mean and deep. Units: 10^12 mole/year.

            # Direct air-temperature dependence of carbonate burial
            f_caco3_t = 0.0
            if f_caco3_t > 0:
                bc *= math.exp((atm_temp - atm_temp_0) * f_caco3_t)

#
#   Convert ocean and atmospheric concentrations to inventories.
#
            for i_s in range(2):
                ocn_inv[i_s] = ocn_conc[i_s] * ocn_vol # moles
            atm_inv = pco2 * MOLpPPM # atmospheric inventory (10^12 moles)

#
#   Adjust inventories by burial, weathering, and gas-exchange fluxes.
#   Atmosphere-ocean exchange is constrained by conservation of matter,
#   but exhanges with solid earth are not because there's no solid carbon inventory.
#
            ocn_inv[i_alk] += ( weat_carb + weat_sil - bc) * 2.0 * delta_time # ocean inventory of alkalinity (10^12 moles)
            ocn_inv[i_tc] +=  ( weat_carb + weat_sil - bc - ge ) * delta_time # ocean inventory of total carbon (10^12 moles)
            atm_inv += (ge + degas + co2_emissions - weat_sil) * delta_time       # atmospheric inventory of CO2 (10^12 moles)

#
#   Now that we've adjusted the inventories, we need to update the concentrations to match.
#

            for i_s in range(2):
                ocn_conc[i_s] = ocn_inv[i_s] / ocn_vol # ocean concentrations of alkalinity and total carbon in moles per cubic meter
            pco2 = atm_inv / MOLpPPM #atmospheric concentration of CO2, ppm.
#
#   Finally. We're done with the time step and we can record the data.
#
            if data_index >= leadin_index:
                # "year", "tco2", "pCO2", "alk", "CO3", "WeatC", "WeatS", "TotW", "BurC", "Degas", "Emiss", "Tatm", "Tocn")

                data[data_index - leadin_index,:] = (year, ocn_conc[i_tc],
                    pco2, ocn_conc[i_alk], co3, weat_carb, weat_sil, weat_carb + weat_sil,
                    bc, degas, co2_emissions * 12/1.E+3, atm_temp, ocn_temp)
            if trace:
                if ((period > 0) or (itime > 0)) and (int(year) % trace_interval == 0):
                    print ("Period %2d, Step %5d, Year %8d" % (period, itime, year))

            data_index += 1
            delta_time = time_step
        #
        # We've reached the end of a period.
        # Do we dump a big slug of CO2 into the atmosphere at the transition to the next period?
        #
        atm_inv += co2_slug[period]
        pco2 = atm_inv/MOLpPPM
    df = pandas.DataFrame(data, columns = data_names)
    return df

def save(df, filename):
    df.to_csv(filename, index = False)

#
# BIBLIOGRAPHY:
#
#
# Archer, D., Kheshgi, H., and Maier-Reimer, E., 1997, "Multiple timescales for neutralization of fossil fuel CO2" GRL 24, 405-08.
# Archer, D., et al., 2009, "Atmospheric Lifetime of Fossil Fuel Carbon Dioxide" Annu. Rev. Earth Planet. Sci. 37, 117-34.
# Berner, R.A., 1991, "A Model for amospheric CO2 over Phanerozoic time" Am. J. Sci. 291, 339-76.
# Berner, R.A., 1994, "GEOCARB II: A revised model of atmospheric CO2 over Phanerozoic time" Am. J. Sci. 294, 56-91.
# Berner, R.A., 2001, "GEOCARB III: A revised model of atmospheric CO2 over Phanerozoic time" Am. J. Sci. 301, 182-204.
# Berner, R.A., Lenton, T., Tester, M., and Beerling, D.J., 1998, "The carbon cycle and CO2 over Phanerozoic time: The role of land plants" Phil. Trans.: Biol. Sci. 353, 75-82.
# Culberson, C.H., and Pytkowicz, R.M., 1973, "Ionization of water in seawater," Marine Chem., 1, 309-16.
# Culkin, F., 1965, "The Major Constituents of Seawater" in J.P. Riley & G. Skirrow, eds., Chemical Oceanography, Vol. 1, Ch. 4, 1st ed. (Academic Press), 121-61.
# Kester, D.R., and Pytkowicz, R.M., "Determination of the apparent dissociation constants of phosphoric acid in seawater," Limnol & Oceanogr. 12, 243-52 (1967).
# Mehrbach, C., Culberson, C.H., Hawley, J.E., and Pytkowicz, R.M., 1973, "Measurement of the apparent dissociation constnts of carbonic acid in seawater at atmospheric pressure," Limnol. & Oceanogr. 18, 897-907.
# Peng, T-H, Takahashi, T., Broecker, W.S., and Olaffson, J., 1987, "Seasonal variability of carbon dioxide, nutrients, and oxxygen in the northern North Atlantic surface water: Observations and a model," Tellus 39B, 439-58.
# Otto-Bleisner, B.L., 1995 "Continental drift, runoff, and weathering feedbacks: Implications from climate model experiments," JGR 100, 11,527-48
# Royer, D.L., Berner, R.A., and Beerling, D.J., "Phanerozoic atmospheric CO2 change: evaluating geochemical and paleobiological approaches" Earth-Science Rev. 54, 349-92.
# Takahashi, T., Williams, R.T., and Boss, D.L., 1982a, "Carbonate Chemistry," in GEOSECS Pacific Expedition, vol. 3. Hydrographic Data, 1973-1974. W.S. Broecker, D.W. Spencer, and H. Craig, eds, (US Gov't Printing Office), 78-82.
# Weiss, R.F., 1974 "Carbon DIoxide in Water and Seawater: The Solubility of a Non-Ideal Gas" Marine Chem. 2, 203-15
#

def main():
    global g_data
    g_data = geocarb((0.,1000.,0.), (7.5,7.5,7.5),
                       (50,1,50),(5E6,100,2E6),
                       3.0, (0,0,0), (30,30,30),
                       (True,True,True), (1.0,1.0,1.0),
                       3E6,
                       (0, numpy.arange(100) * 10, 0))

if __name__ == "__main__":
    main()
