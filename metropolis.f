c **************************************************************
c
c This file contains the subroutines:  metropolis
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku Hu
c
c **************************************************************

      subroutine metropolis(eol,currtem,acepta)

c===============================================================
c== SUBROUTINE FOR METROPOLIS UPDATE OF CONFIGURATIONS 	  ==
c==									 	  ==
c== CALLS: energy,addang,(rand),					  ==
c== dummy (function provided as argument)				  ==
c===============================================================


      include 'INCL.H'

	     real*8 e
        dimension aux_vlvr(mxvr)
        jv = idvr(1+int(nvr*rand()))  ! select var.
        
c================================
c== Get Proposal configuration ==
c================================
         aux_vlvr = vlvr  ! save old
        
         e = rand()
         dv = axvr(jv)*e 
         vlvr(jv) = addang(vrol,dv) 
         enw = energy()
         delta = enw - eol 

c================================
c== check acceptance criteria ===
c================================
        if (delta.LE.0.0d0) then
          eol = enw
          acepta = acepta + 1.0d0
        elseif (exp(-delta/currtem).GT.rand()) then
          eol = enw
          acepta = acepta + 1.0d0
        else
c          vlvr(jv) = vrol
          vlvr = aux_vlvr
        endif
      return
      end


