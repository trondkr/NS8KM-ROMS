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

#define IS4DVAR                 /* Incremental, strong constraint 4DVAR */

#define DJ_GRADPS
#define UV_COR
#define UV_ADV
#define UV_QDRAG
#define UV_VIS2
#define POWER_LAW
#define HDF5
#define MIX_S_UV
#define TS_U3HADVECTION
#define UV_U3HADVECTION
#undef TS_C4VADVECTION
#undef UV_C4ADVECTION
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define SPLINES
#define MASKING
#define DJ_GRADPS
#define  MY25_MIXING
#define SRELAXATION

/* Vertical subgridscale turbulence closure */
#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
#endif
#ifdef GLS_MIXING
# define KANTHA_CLAYSON
# undef  CANUTO_A
# define N2S2_HORAVG
#endif

#define SOLAR_SOURCE
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
#define IMPLICIT_VCONV
#define FORWARD_READ
#define FORWARD_MIXING
#define FORWARD_WRITE
#define OUT_DOUBLE


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

#define NUDGING_COFF                  /* use ana_nudgcoef.h */
#define TCLIMATOLOGY
#define TCLM_NUDGING

