/*
** svn $Id: double_gyre.h 172 2007-04-11 01:45:48Z arango $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for NS8KM
**
** Application flag:   NS8KM
** Input script:       ocean_ns8km.in
*/
#define NLM_DRIVER              /* Nonlinear Basic State trajectory */

/*
**-----------------------------------------------------------------------------
**  Nonlinear basic state tracjectory.
**-----------------------------------------------------------------------------
*/
#if defined NLM_DRIVER
#define UV_ADV
#define TS_C4VADVECTION
#define DJ_GRADPS
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define TS_U3HADVECTION
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define SPLINES
#define POWER_LAW
#define MASKING
#define AVERAGES
#define SOLAR_SOURCE
#define NUDGING_COFF                  /* use ana_nudgcoef.h */
#define TCLIMATOLOGY
#define TCLM_NUDGING
#define SRELAXATION
#define MY25_MIXING
#define NETCDF4
#define DEFLATE
#define CHARNOK /*Charnok Surface Roughness From Wind Stress */

#ifdef GLS_MIXING
# define KANTHA_CLAYSON
# undef  CANUTO_A
# define N2S2_HORAVG /*Horizontal Smoothing of Buoyancy/Shea */
#endif

#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
#endif

#define SPONGE
#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define FORWARD_MIXING
#define FORWARD_WRITE
#define OUT_DOUBLE
#endif

/* ATMOSPHERIC FORCING */
#define BULK_FLUXES        /* turn ON or OFF bulk fluxes computation */
#ifdef BULK_FLUXES
# undef  ANA_RAIN          /* analytical rain fall rate */
# undef  ANA_PAIR          /* analytical surface air pressure */
# undef  ANA_HUMIDITY      /* analytical surface air humidity */
# undef  ANA_CLOUD         /* analytical cloud fraction */
# undef  ANA_TAIR          /* analytical surface air temperature */
# undef  ANA_WINDS         /* analytical surface winds */
# define EMINUSP           /* turn ON internal calculation of E-P */
# define ANA_SRFLUX        /* analytical surface shortwave radiation flux */
# define ALBEDO            /* use albedo equation for shortwave radiation */
# define CLOUDS
# undef  LONGWAVE_OUT      /* compute outgoing longwave radiation */
# define LONGWAVE          /* Compute net longwave radiation internally */
# define COOL_SKIN         /* turn ON or OFF cool skin correction *//* Ikke def hos Frode*/
# define SHORTWAVE
#endif

#define ATM_PRESS          /* use to impose atmospheric pressure onto sea surface */
#define SOLAR_SOURCE       /* define solar radiation source term */
#define SPECIFIC_HUMIDITY  /* if input is specific humidity in kg/kg */

/* TIDES */
#define SSH_TIDES          /* turn on computation of tidal elevation */
#define UV_TIDES           /* turn on computation of tidal currents */
#define ADD_FSOBC          /* Add tidal elevation to processed OBC data */
#define ADD_M2OBC          /* Add tidal currents  to processed OBC data */
#undef RAMP_TIDES         /* Spin up tidal forcing */
#undef WET_DRY

/* ELVER */
#define UV_PSOURCE         /* turn ON or OFF point Sources/Sinks */
#define TS_PSOURCE         /* turn ON or OFF point Sources/Sinks */
