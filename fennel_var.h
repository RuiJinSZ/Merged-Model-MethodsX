/*
** svn $Id: redoxh_var.h 795 2016-05-11 01:42:43Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the Fennel et al. (2006) ecosystem   **
**  model variables that are used in input and output NetCDF files.   **
**  The metadata information is read from "varinfo.dat".              **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iPhyt)')
                idTvar(iPhyt)=varid
              CASE ('idTvar(iZoop)')
                idTvar(iZoop)=varid
              CASE ('idTvar(iLDeN)')
                idTvar(iLDeN)=varid
              CASE ('idTvar(iSDeN)')
                idTvar(iSDeN)=varid
              CASE ('idTvar(iChlo)')
                idTvar(iChlo)=varid
              CASE ('idTvar(iN2__)')
                idTvar(iN2__)=varid
              CASE ('idTvar(iDON_)')
                idTvar(iDON_)=varid
              CASE ('idTvar(irfDON_)')
                idTvar(irfDON_)=varid  
              CASE ('idTvar(iDOP_)')
                idTvar(iDOP_)=varid     
# ifdef CARBON
              CASE ('idTvar(iTIC_)')
                idTvar(iTIC_)=varid
              CASE ('idTvar(iTAlk)')
                idTvar(iTAlk)=varid
              CASE ('idTvar(iLDeC)')
                idTvar(iLDeC)=varid
              CASE ('idTvar(iSDeC)')
                idTvar(iSDeC)=varid
# endif
# ifdef OXYGEN
              CASE ('idTvar(iOxyg)')
                idTvar(iOxyg)=varid
# endif
# ifdef H_SULF
              CASE ('idTvar(iH2S_)')
                idTvar(iH2S_)=varid
              CASE ('idTvar(iSO4_)')
                idTvar(iSO4_)=varid
              CASE ('idTvar(iNO2_)')
                idTvar(iNO2_)=varid
              CASE ('idTvar(iPO4_)')
                idTvar(iPO4_)=varid
              CASE ('idTvar(iDiaz)')
                idTvar(iDiaz)=varid
              CASE ('idTvar(iLDeP)')
                idTvar(iLDeP)=varid
              CASE ('idTvar(iSDeP)')
                idTvar(iSDeP)=varid
# endif
/*
**  Adjoint sensitivity state biological tracers.
*/

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
              CASE ('idTads(iNO3_)')
                idTads(iNO3_)=varid
              CASE ('idTads(iNH4_)')
                idTads(iNH4_)=varid
              CASE ('idTads(iPhyt)')
                idTads(iPhyt)=varid
              CASE ('idTads(iZoop)')
                idTads(iZoop)=varid
              CASE ('idTads(iLDeN)')
                idTads(iLDeN)=varid
              CASE ('idTads(iSDeN)')
                idTads(iSDeN)=varid
              CASE ('idTads(iChlo)')
                idTads(iChlo)=varid
              CASE ('idTads(iN2__)')
                idTads(iN2__)=varid
              CASE ('idTads(iDON_)')
                idTads(iDON_)=varid 
              CASE ('idTads(irfDON_)')
                idTads(irfDON_)=varid
              CASE ('idTads(iDOP_)')
                idTads(iDOP_)=varid     
# ifdef CARBON
              CASE ('idTads(iTIC_)')
                idTads(iTIC_)=varid
              CASE ('idTads(iTAlk)')
                idTads(iTAlk)=varid
              CASE ('idTads(iLDeC)')
                idTads(iLDeC)=varid
              CASE ('idTads(iSDeC)')
                idTads(iSDeC)=varid
# endif
# ifdef OXYGEN
              CASE ('idTads(iOxyg)')
                idTads(iOxyg)=varid
# endif
#endif

/*
**  Biological tracers open boundary conditions.
*/

              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid

              CASE ('idTbry(iwest,iNH4_)')
                idTbry(iwest,iNH4_)=varid
              CASE ('idTbry(ieast,iNH4_)')
                idTbry(ieast,iNH4_)=varid
              CASE ('idTbry(isouth,iNH4_)')
                idTbry(isouth,iNH4_)=varid
              CASE ('idTbry(inorth,iNH4_)')
                idTbry(inorth,iNH4_)=varid

              CASE ('idTbry(iwest,iPhyt)')
                idTbry(iwest,iPhyt)=varid
              CASE ('idTbry(ieast,iPhyt)')
                idTbry(ieast,iPhyt)=varid
              CASE ('idTbry(isouth,iPhyt)')
                idTbry(isouth,iPhyt)=varid
              CASE ('idTbry(inorth,iPhyt)')
                idTbry(inorth,iPhyt)=varid

              CASE ('idTbry(iwest,iZoop)')
                idTbry(iwest,iZoop)=varid
              CASE ('idTbry(ieast,iZoop)')
                idTbry(ieast,iZoop)=varid
              CASE ('idTbry(isouth,iZoop)')
                idTbry(isouth,iZoop)=varid
              CASE ('idTbry(inorth,iZoop)')
                idTbry(inorth,iZoop)=varid

              CASE ('idTbry(iwest,iLDeN)')
                idTbry(iwest,iLDeN)=varid
              CASE ('idTbry(ieast,iLDeN)')
                idTbry(ieast,iLDeN)=varid
              CASE ('idTbry(isouth,iLDeN)')
                idTbry(isouth,iLDeN)=varid
              CASE ('idTbry(inorth,iLDeN)')
                idTbry(inorth,iLDeN)=varid

              CASE ('idTbry(iwest,iSDeN)')
                idTbry(iwest,iSDeN)=varid
              CASE ('idTbry(ieast,iSDeN)')
                idTbry(ieast,iSDeN)=varid
              CASE ('idTbry(isouth,iSDeN)')
                idTbry(isouth,iSDeN)=varid
              CASE ('idTbry(inorth,iSDeN)')
                idTbry(inorth,iSDeN)=varid

              CASE ('idTbry(iwest,iChlo)')
                idTbry(iwest,iChlo)=varid
              CASE ('idTbry(ieast,iChlo)')
                idTbry(ieast,iChlo)=varid
              CASE ('idTbry(isouth,iChlo)')
                idTbry(isouth,iChlo)=varid
              CASE ('idTbry(inorth,iChlo)')
                idTbry(inorth,iChlo)=varid

              CASE ('idTbry(iwest,iN2__)')
                idTbry(iwest,iN2__)=varid
              CASE ('idTbry(ieast,iN2__)')
                idTbry(ieast,iN2__)=varid
              CASE ('idTbry(isouth,iN2__)')
                idTbry(isouth,iN2__)=varid
              CASE ('idTbry(inorth,iN2__)')
                idTbry(inorth,iN2__)=varid

              CASE ('idTbry(iwest,iDON_)')
                idTbry(iwest,iDON_)=varid
              CASE ('idTbry(ieast,iDON_)')
                idTbry(ieast,iDON_)=varid
              CASE ('idTbry(isouth,iDON_)')
                idTbry(isouth,iDON_)=varid
              CASE ('idTbry(inorth,iDON_)')
                idTbry(inorth,iDON_)=varid

              CASE ('idTbry(iwest,irfDON_)')
                idTbry(iwest,irfDON_)=varid
              CASE ('idTbry(ieast,irfDON_)')
                idTbry(ieast,irfDON_)=varid
              CASE ('idTbry(isouth,irfDON_)')
                idTbry(isouth,irfDON_)=varid
              CASE ('idTbry(inorth,irfDON_)')
                idTbry(inorth,irfDON_)=varid  

              CASE ('idTbry(iwest,iDOP_)')
                idTbry(iwest,iDOP_)=varid
              CASE ('idTbry(ieast,iDOP_)')
                idTbry(ieast,iDOP_)=varid
              CASE ('idTbry(isouth,iDOP_)')
                idTbry(isouth,iDOP_)=varid
              CASE ('idTbry(inorth,iDOP_)')
                idTbry(inorth,iDOP_)=varid  

#ifdef CARBON
              CASE ('idTbry(iwest,iTIC_)')
                idTbry(iwest,iTIC_)=varid
              CASE ('idTbry(ieast,iTIC_)')
                idTbry(ieast,iTIC_)=varid
              CASE ('idTbry(isouth,iTIC_)')
                idTbry(isouth,iTIC_)=varid
              CASE ('idTbry(inorth,iTIC_)')
                idTbry(inorth,iTIC_)=varid

              CASE ('idTbry(iwest,iTAlk)')
                idTbry(iwest,iTAlk)=varid
              CASE ('idTbry(ieast,iTAlk)')
                idTbry(ieast,iTAlk)=varid
              CASE ('idTbry(isouth,iTAlk)')
                idTbry(isouth,iTAlk)=varid
              CASE ('idTbry(inorth,iTAlk)')
                idTbry(inorth,iTAlk)=varid

              CASE ('idTbry(iwest,iLDeC)')
                idTbry(iwest,iLDeC)=varid
              CASE ('idTbry(ieast,iLDeC)')
                idTbry(ieast,iLDeC)=varid
              CASE ('idTbry(isouth,iLDeC)')
                idTbry(isouth,iLDeC)=varid
              CASE ('idTbry(inorth,iLDeC)')
                idTbry(inorth,iLDeC)=varid

              CASE ('idTbry(iwest,iSDeC)')
                idTbry(iwest,iSDeC)=varid
              CASE ('idTbry(ieast,iSDeC)')
                idTbry(ieast,iSDeC)=varid
              CASE ('idTbry(isouth,iSDeC)')
                idTbry(isouth,iSDeC)=varid
              CASE ('idTbry(inorth,iSDeC)')
                idTbry(inorth,iSDeC)=varid
#endif
#ifdef OXYGEN
              CASE ('idTbry(iwest,iOxyg)')
                idTbry(iwest,iOxyg)=varid
              CASE ('idTbry(ieast,iOxyg)')
                idTbry(ieast,iOxyg)=varid
              CASE ('idTbry(isouth,iOxyg)')
                idTbry(isouth,iOxyg)=varid
              CASE ('idTbry(inorth,iOxyg)')
                idTbry(inorth,iOxyg)=varid
#endif
#ifdef H_SULF
              CASE ('idTbry(iwest,iH2S_)')
                idTbry(iwest,iH2S_)=varid
              CASE ('idTbry(ieast,iH2S_)')
                idTbry(ieast,iH2S_)=varid
              CASE ('idTbry(isouth,iH2S_)')
                idTbry(isouth,iH2S_)=varid
              CASE ('idTbry(inorth,iH2S_)')
                idTbry(inorth,iH2S_)=varid
              
              CASE ('idTbry(iwest,iSO4_)')
                idTbry(iwest,iSO4_)=varid
              CASE ('idTbry(ieast,iSO4_)')
                idTbry(ieast,iSO4_)=varid
              CASE ('idTbry(isouth,iSO4_)')
                idTbry(isouth,iSO4_)=varid
              CASE ('idTbry(inorth,iSO4_)')
                idTbry(inorth,iSO4_)=varid

              CASE ('idTbry(iwest,iNO2_)')
                idTbry(iwest,iNO2_)=varid
              CASE ('idTbry(ieast,iNO2_)')
                idTbry(ieast,iNO2_)=varid
              CASE ('idTbry(isouth,iNO2_)')
                idTbry(isouth,iNO2_)=varid
              CASE ('idTbry(inorth,iNO2_)')
                idTbry(inorth,iNO2_)=varid
              
              CASE ('idTbry(iwest,iPO4_)')
                idTbry(iwest,iPO4_)=varid
              CASE ('idTbry(ieast,iPO4_)')
                idTbry(ieast,iPO4_)=varid
              CASE ('idTbry(isouth,iPO4_)')
                idTbry(isouth,iPO4_)=varid
              CASE ('idTbry(inorth,iPO4_)')
                idTbry(inorth,iPO4_)=varid
              
              CASE ('idTbry(iwest,iDiaz)')
                idTbry(iwest,iDiaz)=varid
              CASE ('idTbry(ieast,iDiaz)')
                idTbry(ieast,iDiaz)=varid
              CASE ('idTbry(isouth,iDiaz)')
                idTbry(isouth,iDiaz)=varid
              CASE ('idTbry(inorth,iDiaz)')
                idTbry(inorth,iDiaz)=varid
              
              CASE ('idTbry(iwest,iLDeP)')
                idTbry(iwest,iLDeP)=varid
              CASE ('idTbry(ieast,iLDeP)')
                idTbry(ieast,iLDeP)=varid
              CASE ('idTbry(isouth,iLDeP)')
                idTbry(isouth,iLDeP)=varid
              CASE ('idTbry(inorth,iLDeP)')
                idTbry(inorth,iLDeP)=varid

              CASE ('idTbry(iwest,iSDeP)')
                idTbry(iwest,iSDeP)=varid
              CASE ('idTbry(ieast,iSDeP)')
                idTbry(ieast,iSDeP)=varid
              CASE ('idTbry(isouth,iSDeP)')
                idTbry(isouth,iSDeP)=varid
              CASE ('idTbry(inorth,iSDeP)')
                idTbry(inorth,iSDeP)=varid
#endif

/*
**  Biological tracers point Source/Sinks (river runoff).
*/

              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iNH4_)')
                idRtrc(iNH4_)=varid
              CASE ('idRtrc(iPhyt)')
                idRtrc(iPhyt)=varid
              CASE ('idRtrc(iZoop)')
                idRtrc(iZoop)=varid
              CASE ('idRtrc(iLDeN)')
                idRtrc(iLDeN)=varid
              CASE ('idRtrc(iSDeN)')
                idRtrc(iSDeN)=varid
              CASE ('idRtrc(iChlo)')
                idRtrc(iChlo)=varid
              CASE ('idRtrc(iN2__)')
                idRtrc(iN2__)=varid
              CASE ('idRtrc(iDON_)')
                idRtrc(iDON_)=varid
              CASE ('idRtrc(irfDON_)')
                idRtrc(irfDON_)=varid 
              CASE ('idRtrc(iDOP_)')
                idRtrc(iDOP_)=varid      
#ifdef CARBON
              CASE ('idRtrc(iTIC_)')
                idRtrc(iTIC_)=varid
              CASE ('idRtrc(iTAlk)')
                idRtrc(iTAlk)=varid
              CASE ('idRtrc(iLDeC)')
                idRtrc(iLDeC)=varid
              CASE ('idRtrc(iSDeC)')
                idRtrc(iSDeC)=varid
#endif
#ifdef OXYGEN
              CASE ('idRtrc(iOxyg)')
                idRtrc(iOxyg)=varid
#endif
#ifdef H_SULF
              CASE ('idRtrc(iH2S_)')
                idRtrc(iH2S_)=varid
              CASE ('idRtrc(iSO4_)')
                idRtrc(iSO4_)=varid
              CASE ('idRtrc(iNO2_)')
                idRtrc(iNO2_)=varid
              CASE ('idRtrc(iPO4_)')
                idRtrc(iPO4_)=varid
              CASE ('idRtrc(iDiaz)')
                idRtrc(iDiaz)=varid
              CASE ('idRtrc(iLDeP)')
                idRtrc(iLDeP)=varid
              CASE ('idRtrc(iSDeP)')
                idRtrc(iSDeP)=varid
#endif

#ifdef DIAGNOSTICS_BIO

/*
**  Biological tracers term diagnostics.
*/
# ifdef DENITRIFICATION
              CASE ('iDbio2(iDNIT)')
                iDbio2(iDNIT)=varid
# endif
# ifdef CARBON
              CASE ('iDbio2(iCOfx)')
                iDbio2(iCOfx)=varid
              CASE ('iDbio2(ipCO2)')
                iDbio2(ipCO2)=varid
# endif
# ifdef OXYGEN
              CASE ('iDbio2(iO2fx)')
                iDbio2(iO2fx)=varid
# endif
# ifdef H_SULF
              CASE ('iDbio2(iSRR_)')
                iDbio2(iSRR_)=varid
              CASE ('iDbio2(iBFlx)')
                iDbio2(iBFlx)=varid
              CASE ('iDbio2(iSRRb)')
                iDbio2(iSRRb)=varid

              CASE ('iDbio3(iSOx_)')
                iDbio3(iSOx_)=varid
              CASE ('iDbio3(iOxic)')
                iDbio3(iOxic)=varid
              CASE ('iDbio3(iSRRa)')
                iDbio3(iSRRa)=varid
              CASE ('iDbio3(iDno2)')
                iDbio3(iDno2)=varid
              CASE ('iDbio3(iAnx_)')
                iDbio3(iAnx_)=varid
              CASE ('iDbio3(iNit1)')
                iDbio3(iNit1)=varid
              CASE ('iDbio3(iNit2)')
                iDbio3(iNit2)=varid
              CASE ('iDbio3(iDno3)')
                iDbio3(iDno3)=varid
              CASE ('iDbio3(iSNx_)')
                iDbio3(iSNx_)=varid
              CASE ('iDbio3(iSNy_)')
                iDbio3(iSNy_)=varid     
# endif

# ifdef USECOS_BURIAL
              CASE ('iDbio2(iNbur)')
                iDbio2(iNbur)=varid
              CASE ('iDbio2(iNbot)')
                iDbio2(iNbot)=varid  
              CASE ('iDbio2(iCbot)')
                iDbio2(iCbot)=varid
              CASE ('iDbio2(iCbur)')
                iDbio2(iCbur)=varid 
              CASE ('iDbio2(iPbot)')
                iDbio2(iPbot)=varid
              CASE ('iDbio2(iPbur)')
                iDbio2(iPbur)=varid     
# endif

              CASE ('iDbio3(iPPro)')
                iDbio3(iPPro)=varid
              CASE ('iDbio3(iNO3u)')
                iDbio3(iNO3u)=varid
              CASE ('iDbio3(iNFix)')
                iDbio3(iNFix)=varid 
              CASE ('iDbio3(iPO4u)')
                iDbio3(iPO4u)=varid  
#endif  


