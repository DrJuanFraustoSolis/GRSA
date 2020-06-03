c Golden Ratio Simulated annealing (GRSA)
c Copyright (C) 2020  Dr. Juan Paulo Sánchez Hernández, and Dr. Juan Frausto Solis
c Copyright (C) 2005 Frank Eisenmenger, U.H.E. Hansmann, Shura Hayryan, Chin-Ku Hu

c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or (at
c your option) any later version.
c 
c This program is distributed in the hope that it will be useful, but
c WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c General Public License for more details.
c 
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
c USA.


c***********************GRSA************************************
c
c -This file contains a hybridization of simulated annealing called
c  Golden Ratio Simulated annealing (GRSA)
c -The GRSA has three fundamental parts: 1) The cooling scheme, 
c  2) The overheating strategies, and 3) The convergence criteria
c 
c -The firts version was developed by Dr. Juan Paulo Sánchez Hernández, 
c and Dr. Juan Frausto Solis
c  2015-2020
c***************************************************************


      subroutine  grsa(energia)

C --------------------------------------------------------------
C PURPOSE: SIMULATED ANNEALING SEARCH OF LOWEST-POTENTIAL-ENERGY
C          CONFORMATIONS OF PROTEINS
C
C CALLS: addang,energy,metropolis,outvar,outpdb,rgyr,setvar,zimmer
C
C ---------------------------------------------------------------
C ---------------------------------------------------------------

C------------------------------------------
C---------------------GRSA-----------------
C------------------------------------------

      include 'INCL.H'

      character(8) x2
      logical lrand
      parameter(lrand=.true.)
	logical bandera
	logical bandera1
	logical bandera2
	logical bandera3
	logical bandera4

        dimension vlvrm(mxvr)
	dimension idvrm(mxvr)
	dimension axvrm(mxvr)

        dimension vlvr_aux(mxvr)
	dimension idvr_aux(mxvr)
	dimension axvr_aux(mxvr)

	real*8 energia,e,bangl
	real inicio,final,total
	real inicioT,namefile
	integer seed1,p,THREADS,ID,MaxThreads
 	
c	  integer time2
c      character*10 a(3)
c      integer datetime(8)
c      character*50 name2      
      
        bandera = .true.
 	bandera1 = .true.
 	bandera2 = .true.
 	bandera3 = .true.
 	bandera4 = .true.
	val_pen = 90.0d0

c	seed=86456
	paro_recta = 0.000000000001d0
	recta = 0
	C = 0
	D = 0
	cnt = 1
	cnt2 = 0

	inicioT = secnds(0.0)
       nresi=irsml2(1)-irsml1(1) + 1     
c ==================
c == random start ==
c ==================

       if(lrand) then
        	do i=1,nvr
         	iv=idvr(i)  ! provides index of non-fixed variable
        	 e = rand(int(inicioT))
        	 dv=axvr(iv)*e
	         vr=addang(pi,dv)
          	vlvr(iv)=vr
        	enddo

  	  end if		
        eol = energy()
	ymin = eol
 	epsilon = eol

c =======================================================
c == Write start configuration in pdb-format into file ==
c =======================================================
c      do i=1,ntlml
c        call outpdb(i,11)
c      end do

c==========================
c==	Simulated Annealing  ==
c==========================
      write(*,*) 'SIMULATED ANNEALING'
	call cpu_time(start)

c ====================================================
c == Initial parameters for Simulated Annealing ==
c ====================================================
	
c ====================================================
c ========== Parameters of Met-Enkephaline ==============
c ====================================================	          

	currtem = 226624331780759000000000000000000000.0d0 !Initial Temperature Met-enkephaline
	paro = 0.00000014573301692707600 !Final Temperature Met-enkephaline
      

	temp_aureo = currtem*0.618
	temp_aureo1 = temp_aureo*0.618
	temp_aureo2 = temp_aureo1*0.618
	temp_aureo3 = temp_aureo2*0.618
	temp_aureo4 = temp_aureo3*0.618

        alfa = 0.70
	blmax = 360.0d0
	bbeta = 1.00067911076425 !

	temp1=1.9E-4
  	ban1=1
	ban2=1
	ban3=1
	ban4=1
	band5=1
	band6=1	
	band7=1
c =================================
c == Start the Simulated Annealing == 
c =================================
      do while(currtem.GE.paro)
        propon = 0.0d0
        acepta = 0.0d0
        ycurr = 0.0d0
c =====================================
c == Start Metropolis ==
c == Crescent Markov Chain ==
c =====================================
        do while(propon.LE.NINT(blmax))
          propon = propon + 1.0d0
          ycurr = eol

          call metropolis(eol,currtem,acepta)
c =================================================
c == Store the local lowest-energy conformation  ==
c =================================================
          if (eol.LT.ycurr) then
            ycurr = eol
          end if
c =================================================
c == Store the global lowest-energy conformation ==
c =================================================
          if (eol.LT.ymin) then
            ymin = eol
            write(*,*) 'MINIMA:', currtem, ymin

c ==================================
c == backup of Minimum point   ==
c == idvrm, axvrm, vlvrm  ==
c ==================================
       	do i=1,nvr
		  idvrm(i)  = idvr(i)
		  iv        = idvr(i)
		  axvrm(iv) = axvr(iv)
		  vlvrm(iv) = vlvr(iv)
       	enddo
          end if
c==================================================
c   The convergence criteria for dynamic equilibrium
c==================================================

			if ((currtem.LT.val_pen).AND. (cnt.LT.3)) then
				D = D + ymin
				C = C + cnt * ymin
				cnt = cnt + 1
			end if

			if (cnt.EQ.3) then
		        recta = (((12 * C)-(6*(cnt-1)*D))/(cnt**3 - cnt))	
			cnt = 1	

			 if (recta.LT.paro_recta) then
c----------------------------Overheating strategy-------------------------------------------
				currtem = currtem*0.618
				bangl= bangl+1
		         endif		
			
			  if (bangl.EQ.2) then
				currtem = paro
			  endif

			endif


        enddo
c ========================
c == End of Metropolis ==
c ========================

c ==============================
c == Lower the Temperature ==
c ==============================
        currtem = alfa * currtem
	write(*,*) currtem,ymin
c ==============================
c Cooling Scheme by golden ratio
c ==============================

	if ((currtem.LT.temp_aureo).AND. (bandera)) then
		alfa = 0.75d0
		bandera = .false.
	else if ((currtem.LT.temp_aureo1).AND. (bandera1)) then
		alfa = 0.80d0
		bandera1 = .false.
	else if ((currtem.LT.temp_aureo2).AND. (bandera2)) then
		alfa = 0.85d0
		bandera2 = .false.
	else if ((currtem.LT.temp_aureo3).AND. (bandera3)) then
		alfa = 0.90d0
		bandera3 = .false.
	else if ((currtem.LT.temp_aureo4).AND. (bandera4)) then
		alfa = 0.95d0
		bandera4 = .false.
c--------------------Overheating strategy-------------------------------------------------
		currtem = currtem+currtem
	endif

c ==============================
c == Crescent Markov Chain ==
c ==============================
        blmax = bbeta * blmax

      enddo
	energia = ymin
c ====================================
c == End of the Simulated Annealing ==
c ====================================
	  namefile=energia*10000
	  write(x2,'(I8)') int(namefile) 
	  open(17, file='./RESULTADOS/deltas'//x2//'.txt',status='new')
	  write(17,*) 'Delta min y max',deltamin,deltamax	
          open(16, file='./RESULTADOS/Estruc'//x2//'.pdb',status='new')	  
	  close(17)
c =======================================================
c == Write start configuration in pdb-format into file ==
c =======================================================
       do i=1,ntlml
         call outpdb(i,16)
       end do

       close(16)

      end
