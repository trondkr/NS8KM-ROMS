/*
** Svn $Id: basin.h 8 2007-02-06 19:00:29Z arango $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for NS8KM
*/

#define RST_SINGLE         /*NB anne change define if single precision restart fields */
#undef  PERFECT_RESTART    /* use to include perfect restart variables */
#define CURVGRID           /* define if using  curvilinear coordinate grid*/

#define UV_ADV             /* turn ON or OFF advection terms */
#define UV_COR             /* turn ON or OFF Coriolis term */
#undef  UV_VIS2            /* turn ON or OFF Laplacian horizontal mixing, Frode off */
#undef  UV_VIS4            /* turn ON or OFF biharmonic horizontal mixing */
#undef  UV_SADVECTION      /* turn ON or OFF splines vertical advection */
#define TS_C4VADVECTION    /* define if 4th-order centered vertical advection */
#define UV_QDRAG           /* turn ON or OFF quadratic bottom friction */

#undef  VISC_GRID          /* viscosity coefficient scaled by grid size */
#define NONLIN_EOS         /* define if using nonlinear equation of state */
#undef  MIX_GEO_UV         /* mixing on geopotential (constant Z) surfaces */
#undef  WJ_GRADP           /* Weighted density Jacobian (Song, 1998) */
#define DJ_GRADPS          /* Splines density  Jacobian (Shchepetkin, 2000) */
#undef  DIFF_GRID          /* diffusion coefficient scaled by grid size */

#undef  TS_DIF2            /* turn ON or OFF Laplacian horizontal mixing */
#undef  TS_DIF4            /* turn ON or OFF biharmonic horizontal mixing */
#define TS_U3HADVECTION    /* define if 3rd-order upstream horiz. advection */
#undef  TS_A4HADVECTION    /* define if 4th-order Akima horiz. advection */
#undef  TS_C4HADVECTION    /* define if 4th-order centered horizontal advection */

#undef  TS_MPDATA          /* define if recursive MPDATA 3D advection */

#undef  TS_A4VADVECTION    /* define if 4th-order Akima vertical advection */
#undef  TS_SVADVECTION     /* define if splines vertical advection */
#define TS_C4VADVECTION    /* define if 4th-order centered vertical advection */
#undef  TS_SMAGORINSKY     /* define if Smagorinsky-like diffusion */

#undef  MIX_GEO_TS         /* mixing on geopotential (constant Z) surfaces, Frode undef */
#undef  MIX_S_UV           /* mixing along constant S-surfaces,   Frode -*/

#define SALINITY           /* define if using salinity */
#define SOLVE3D            /* define if solving 3D primitive equations */
#undef  BODYFORCE          /* define if applying stresses as bodyforces */
#undef SPLINES            /* turn ON or OFF parabolic splines reconstruction */
#define MASKING            /* define if there is land in the domain */
#define AVERAGES           /* define if writing out time-averaged data */
#define AVERAGES_AKT       /*undef Frode, define if writing out time-averaged AKt*/
#define AVERAGES_AKS       /*undef Frode, define if writing out time-averaged AKs*/
#undef  AVERAGES_AKV       /* define if writing out time-averaged AKv */
#undef STATIONS           /* define if writing out station data */
#undef STATIONS_CGRID     /* define if extracting data at native C-grid */

#undef  BVF_MIXING         /* define if Brunt_Vaisala frequency mixing */
#undef  LMD_MIXING         /* define if Large et al. (1994) interior closure */
#define  MY25_MIXING        /* define if Mellor/Yamada level-2.5 mixing */
#undef GLS_MIXING         /* Activate Generic Length-Scale mixing */

#ifdef GLS_MIXING
# define N2S2_HORAVG       /* Activate horizontal smoothing of buoyancy/shear */
#endif
#ifdef LMD_MIXING
# undef  LMD_BKPP          /* use if bottom boundary layer KPP mixing */
# undef  LMD_CONVEC        /* use to add convective mixing due to shear instability */
# undef  LMD_DDMIX         /* use to add double-diffusive mixing */
# undef  LMD_NONLOCAL      /* use if nonlocal transport */
# undef  LMD_RIMIX         /* use to add diffusivity due to shear instability */
# undef  LMD_SHAPIRO       /* use if Shapiro filtering boundary layer depth */
# undef  LMD_SKPP          /* use if surface boundary layer KPP mixing */
#endif

#define ANA_BSFLUX         /* analytical bottom salinity flux */
#define ANA_BTFLUX         /* analytical bottom temperature flux */
#undef  ANA_GRID           /* analytical grid set-up */
#undef  ANA_INITIAL        /* analytical initial conditions */
#undef  ANA_MASK           /* analytical land/sea masking */
#undef  ANA_MEANRHO
#undef  ANA_SSFLUX         /* analytical surface salinity flux */
#undef  ANA_STFLUX         /* analytical surface temperature flux */
#undef  ANA_SMFLUX         /* analytical surface momentum stress */

#undef  WESTERN_WALL
#undef  NORTHERN_WALL
#undef  SOUTHERN_WALL
#undef  EASTERN_WALL

#define EAST_FSCHAPMAN     /* free-surface Chapman condition */
#define EAST_M2FLATHER     /* 2D momentum Flather condition */
#define EAST_M3NUDGING     /* 3D momentum passive/active nudging term */
#define EAST_M3RADIATION   /* 3D momentum radiation condition */
#define EAST_TNUDGING      /* tracers passive/active nudging term */
#define EAST_TRADIATION    /* tracers radiation condition */

#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_M3NUDGING
#define WEST_M3RADIATION
#define WEST_TNUDGING
#define WEST_TRADIATION

#define NORTH_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_M3NUDGING
#define NORTH_M3RADIATION
#define NORTH_TNUDGING
#define NORTH_TRADIATION

/* CLIMATOLOGY */
#define M2CLIMATOLOGY      /* define 2D momentum climatology */
#define M3CLIMATOLOGY      /* define 3D momentum climatology */
#define TCLIMATOLOGY       /* define tracers climatology */
#define M3CLM_NUDGING      /* nudging 3D momentum climatology */
#define TCLM_NUDGING       /* nudging tracers climatology */

#undef  WRF_COUPLING       /* coupling to WRF atmospheric model */

/* SEA ICE */
#define ICE_MODEL          /* turn ON or OFF sea ice module */

/* ATMOSPHERIC FORCING */
#define BULK_FLUXES        /* turn ON or OFF bulk fluxes computation */
#ifdef BULK_FLUXES
# ifdef ICE_MODEL
#  define ICE_BULK_FLUXES
# endif
# undef  ANA_RAIN          /* analytical rain fall rate */
# undef  ANA_PAIR          /* analytical surface air pressure */
# undef  ANA_HUMIDITY      /* analytical surface air humidity */
# undef  ANA_CLOUD         /* analytical cloud fraction */
# undef  ANA_TAIR          /* analytical surface air temperature */
# undef  ANA_WINDS         /* analytical surface winds */
# define EMINUSP           /* turn ON internal calculation of E-P */
# define ANA_SRFLUX        /* analytical surface shortwave radiation flux */
# define ALBEDO            /* use albedo equation for shortwave radiation */
# undef  LONGWAVE_OUT      /* compute outgoing longwave radiation */
# define LONGWAVE          /* Compute net longwave radiation internally */
# define COOL_SKIN         /* turn ON or OFF cool skin correction *//* Ikke def hos Frode*/
#endif

#define ATM_PRESS          /* use to impose atmospheric pressure onto sea surface */
#define SOLAR_SOURCE       /* define solar radiation source term */
#define SPECIFIC_HUMIDITY  /* if input is specific humidity in kg/kg */

/* TIDES */
#define SSH_TIDES          /* turn on computation of tidal elevation */
#define UV_TIDES           /* turn on computation of tidal currents */
#define ADD_FSOBC          /* Add tidal elevation to processed OBC data */
#define ADD_M2OBC          /* Add tidal currents  to processed OBC data */
#undef  RAMP_TIDES         /* Spin up tidal forcing */


/* ELVER */
#undef UV_PSOURCE         /* turn ON or OFF point Sources/Sinks */
#undef TS_PSOURCE         /* turn ON or OFF point Sources/Sinks */

/* SEA ICE */
#ifdef ICE_MODEL
# define ICE_THERMO
# define ICE_MK
# undef  ICE_ALB_EC92
# define ICE_MOMENTUM
# undef  ICE_MOM_BULK
# define ICE_EVP
# define ICE_ADVECT
# define ICE_SMOLAR
# define ICE_UPWIND

# define WEST_AICLAMPED
# define WEST_HICLAMPED
# define WEST_HSNCLAMPED
# define WEST_TICLAMPED
# define WEST_SFWATCLAMPED
# define WEST_SIG11CLAMPED
# define WEST_SIG22CLAMPED
# define WEST_SIG12CLAMPED
# define WEST_MIGRADIENT

# define NORTH_AICLAMPED
# define NORTH_HICLAMPED
# define NORTH_HSNCLAMPED
# define NORTH_TICLAMPED
# define NORTH_SFWATCLAMPED
# define NORTH_SIG11CLAMPED
# define NORTH_SIG22CLAMPED
# define NORTH_SIG12CLAMPED
# define NORTH_MIGRADIENT

# define EAST_AICLAMPED
# define EAST_HICLAMPED
# define EAST_HSNCLAMPED
# define EAST_TICLAMPED
# define EAST_SFWATCLAMPED
# define EAST_SIG11CLAMPED
# define EAST_SIG22CLAMPED
# define EAST_SIG12CLAMPED
# define EAST_MIGRADIENT

# define SOUTH_AICLAMPED
# define SOUTH_HICLAMPED
# define SOUTH_HSNCLAMPED
# define SOUTH_TICLAMPED
# define SOUTH_SFWATCLAMPED
# define SOUTH_SIG11CLAMPED
# define SOUTH_SIG22CLAMPED
# define SOUTH_SIG12CLAMPED
# define SOUTH_MIGRADIENT

#endif

/* --------------------------- */

