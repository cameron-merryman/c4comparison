
# calc_c4_photosynthesis.R
## Description: This script defines a function to simulate leaf-level photosynthesis of C4 plants, 
## using user-specified options.
## Originator: Nick Smith

### function variable descriptions
#### calculated
# vpmax: the maximum rate of PEP carboxylation (µmol m-2 s-1)
#### prescribed
# func: the photosynthesis calculation function to use, options:
  # 1. 'collatz': Collatz et al. (1992)
# tleaf: leaf temperature (°C)
# ca: atmospheric co2 (µmol mol-1)
# par: photosynthetically active radiation at leaf surface (µmol m-2 s-1)
# q_yield: quantum yield (mol mol-1)
# abs: absoptance (unitless)
# z: elevation (m.a.s.l.)
# ci_frac: ratio of intercellular to atmospheric co2 (unitless)
# func: ci response function to use, options:
  # 1. 'ca_frac': ci is calculated as a fraction of ca
# tref_vcmax: reference temperature for vcmax temperature correction (°C)
# vcmax_ref: vcmax at reference temperature (µmol m-2 s-1)
# func_vcmax_tresp: the vcmax temperature response function to use, options:
  # 1. 'collatz': Collatz et al. (1992)
# tref_vpmax: reference temperature for vpmax temperature correction (°C)
# vpmax_ref: vpmax at reference temperature (µmol m-2 s-1)
# func_vcmax_tresp: the vpmax temperature response function to use, options:
  # 1. 'collatz': Collatz et al. (1992)

## packages
library(R.utils)

## source functions
sourceDirectory('Functions')

calc_c4_photosynthesis = function(func = 'collatz',
                                  tleaf = 25, ca = 400, par = 500,
                                  q_yield = 0.066, abs = 0.8,
                                  z = 0,
                                  ca_frac = 0.8, func_ci = 'ca_frac',
                                  tref_vcmax = 25, vcmax_ref = 39, func_vcmax_tresp = 'collatz',
                                  tref_vpmax = 25, vpmax_ref = 0.78, func_vpmax_tresp = 'collatz'){
  
  if(func == 'collatz'){
    
    vcmax = calc_vcmax_tresp(tleaf = tleaf, tref = tref_vcmax, vcmax_ref = vcmax_ref, func = func_vcmax_tresp)
    vpmax = calc_vpmax_tresp(tleaf = tleaf, tref = tref_vpmax, vpmax_ref = vpmax_ref, func = func_vpmax_tresp)
    ci = ci_calc(ca = ca, ca_frac = ca_frac, func = func_ci)
    
    ac = vcmax
    aj = q_yield * abs * par
    ap = vpmax * ci
    
    a = pmin(ac, aj, ap)
    
    return(data.frame(tleaf, ca, par,
                      q_yield, abs,
                      vcmax, vpmax, ci, 
                      ac, aj, ap,
                      a))
    
  } 
  
  else if(func == 'von caemmerer') {
    
    vcmax = calc_vcmax_tresp(tleaf = tleaf, tref = tref_vcmax, vcmax_ref = vcmax_ref, func = func_vcmax_tresp)
    vpmax = calc_vpmax_tresp(tleaf = tleaf, tref = tref_vpmax, vpmax_ref = vpmax_ref, func = func_vpmax_tresp)
    ci = ci_calc(ca = ca, ca_frac = ca_frac, func = func_ci)
    patm = calc_patm(z)
    cm <- ci # check that this is how von caemmerer does this
    kp <- calc_kp_temp_pa(tleaf, z)
    Vp = (cm * vpmax)/(cm + kp)
    ac = ((cm * vpmax)/(cm + kp)) - (0.5 * (0.01 * vcmax)) + (gbs * cm) # Eqn 26 in Von Caemerer
    # leakage <- leakiness * Al # need to fix
    leakage = gbs * ((cm + ((Vp - ac)/gbs)) - cm) #Pulled from Eqn. 22 in Scott paper and Eqn. 4 in Von Caemmerer
    cbs <- calc_cbs(cm, leakage) # Eqn. 2.41?
    chi_bs <- cbs / ca
    oa <- 209460 * 1e-6 * patm
    oi <- oa * (cm/ca) #Pulled from Helens model
    obs <- oi # om + os in Von Caemmerer..? Partial pressure of O2 in the leaf. Gave up and pulled from helens model
    # ac = vcmax * ((cbs - gammastar) / (kr * (1 + obs/ko) + cbs)) #Vpr in von caemmerer?
    # ac = ((cm * vpmax)/(cm + kp))-0.5(0.01 * vcmax) + (gbs * cm)
    # 
    # aj = q_yield * abs * par # need to fix 
    z = 1.25 #??? When Fcyc is equal to 0.5 z = 1.25, I dont really know what that means
    x = 0.4 # Partitioning factor of electron transport rate
    Vc = (cbs * vcmax)/(cbs + 1.210 * ((1 + obs)/ 292)) # Eqn 7 Von Caemmerer
    J <- ((2/z) * Vp) + ((3/z) * ((1 + 7 * (0.5/1310) * obs)/ 3 * cbs) * Vc) # Eqn 29 Von Caemmerer
    aj = (z/2) * x * J - gbs * (cbs - cm) - (0.5 * (0.01 * vcmax)) # I know this is wrong...
    
    a = pmin(ac, aj)
    
    return(data.frame(tleaf, ca, par,
                      q_yield, abs,
                      vcmax, vpmax, cm, 
                      ac, aj,
                      a))
    
  }
  
  else{
    
    return('invalid photosynthesis function specified')
    
  }
  
}
calc_c4_photosynthesis('von caemmerer')
calc_c4_photosynthesis('collatz')
