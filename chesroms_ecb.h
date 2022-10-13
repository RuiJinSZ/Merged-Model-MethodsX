/*
** Include file "chesroms_ecb.h"
*******************************************************************************
** Copyright (c) 2007  ChesROMS Group, version 1.0                           **
********************************************************** Wen Long *********** 
**                                                                           **
**  This is the C-preprocessing options by using the command #define and     **
**  #undef to activate and deactivate ChesROMS options                       **
**                                                                           **
*******************************************************************************
*/

/*psl20180805: After discussion with M.A.M.Friedrichs, the variables and blocks
  of code that used to be associated with the CPP flags
    CARBON,OXYGEN,USECOS_DOM,BIO_SEDIMENT,DENITRIFICATION
  are now permanently activated in ChesROMS-ECB. We still #define these flags
  below as they are recognized by Rutger's trunk (i.e. they show up in the run's
  metadata) and because CARBON is necessary to unlock the variable OCEAN%pH.*/

/*CAREFUL, only use this option in conjunction with the forcing files
    frc_ches_era5_with_downwelling_longwave*.nc. It indicates to ROMS that you
    are feeding the downwelling longwave (rather than the net longwave).
  Moreover, do not forget to set the scaling factor of swrad to "1.0d0" in
    file varinfo.dat if you are prescribing shortwave from ERA5.*/
#define LONGWAVE_OUT

/* Options for 1900s, psl20180508 */
#undef  ATM_CO2_1900
#undef  WATER_TEMP_1900

/* Basic physics options */
#define UV_ADV
#define UV_COR 
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS

/* Basic numerics options */
#define TS_MPDATA
#define DJ_GRADPS
#undef  SPLINES
#define CURVGRID
#define MASKING

/* Outputs, Diagnostics output set below */
#define AVERAGES
#define AVERAGES_FLUXES /* write out surface and bottom average fluxes */
#undef  AVERAGES_AKV
#undef  AVERAGES_AKT
#define STATIONS
#undef  FLOATS

/* Surface and bottom boundary conditions */
#define WIND_MINUS_CURRENT  /*psl20180225, take into account v_water*/
#define LIMIT_STFLX_COOLING /*psl20171215, suppress surf. cooling if freezing*/
#define UV_QDRAG
#define BULK_FLUXES
#define SOLAR_SOURCE
#define ANA_BSFLUX /*Bottom Salinity    Flux*/
#define ANA_BTFLUX /*Bottom Temperature Flux*/
#define EMINUSP

/* Vertical subgridscale turbulence closure */
#define GLS_MIXING
#ifdef  GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
#endif

/* Open boundary condition settings */
#define SSH_TIDES
#ifdef  SSH_TIDES
# define RAMP_TIDES
# define ADD_FSOBC
#endif
#define UV_TIDES
#ifdef  UV_TIDES
# define ADD_M2OBC
#endif

/*restart options*/
#undef  PERFECT_RESTART

#define BIO_FENNEL
#ifdef  BIO_FENNEL
# define CARBON /*This flag is necessary to define OCEAN%pH*/
# define OXYGEN
# define BIO_SEDIMENT
# define DENITRIFICATION
# define pCO2_RZ

# define H_SULF

# ifdef  pCO2_RZ
#  define pCO2_RZ_MILLERO_2010
# endif
# define TALK_NONCONSERV

/*usecos*/
# define USECOS_BURIAL
# define USECOS_SLOPPY

# undef  FIXED_Vp0
/*chesroms-bgc iss and term and light attenuation option*/
# define CHESROMS_BGC_BO2LIM
# define CHESROMS_BGC_ATT
# define WC_NITRIFICATION
# define WC_DENITRIFICATION
# define ISS_2_SIZE_CLASSES

/*Atmospheric Nitrogen Deposition  ABever March 2014*/
/*Use BIO_..._SFLUX to read the values from a netcdf file, otherwise use BIO_ANA_...
to use analytical options. (Sept2015, values in analytical files need modified by the user, currently testing values)*/

# define BIO_NH4_SFLUX
# define BIO_NO3_SFLUX
# define BIO_DON_SFLUX

/*
# define ANA_SPFLUX
# define BIO_ANA_NH4_SFLUX
# define BIO_ANA_NO3_SFLUX
# define BIO_ANA_DON_SFLUX
*/

/*define diagnostic term*/
# define DIAGNOSTICS
# ifdef  DIAGNOSTICS
#  define DIAGNOSTICS_BIO
#  define DIAGNOSTICS_TS
#  define DIAGNOSTICS_UV
#  if defined DIAGNOSTICS_TS && defined DIAGNOSTICS_BIO
/*  PSL 20170424: Use vertically-integrated diagnostics*/
#   define VERTICALLY_INTEGRATED_DIAGNOSTICS
#  endif
# endif

#endif /*BIO_FENNEL*/

/* Hypoxia options: For Malcolm Scully 1-Term DO Passive Tracer Model
   T_PASSIVE and HYPOXIA are both necessary for 1-Term DO, HYPOXIA_MAXVAL_CHECK is optional. */
#define T_PASSIVE
#define HYPOXIA
#undef  HYPOXIA_MAXVAL_CHECK

#if defined T_PASSIVE || defined BIO_FENNEL
# define ANA_SPFLUX /*Surface Passive tracer Flux*/
# define ANA_BPFLUX /*Bottom  Passive tracer Flux*/
#endif

/*NetCDF-4 is NOT enabled on the NetCDF-Fortran build of XSEDE (!)*/
/*  To convince yourself: nf-config --all*/
/*  --has-nc4   -> no*/
#define HDF5
#define DEFLATE

/*define from kalev's file*/





#define ANA_SSS
#define ANA_SST

#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX



#define ANA_BMFLUX
#define ANA_HUMID

#define ANA_LRFLUX