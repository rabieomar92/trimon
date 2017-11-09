!     MODULE CriticalityCalculation
!
!     THIS MODULE CONTAINS THE NECESSARY METHODS THAT ENABLE THE CODE USER AND
!     PROGRAMMERS TO DO CRITICALITY CALCULATION OF A TRIGA REACTOR CORE.
!
!     AUTHOR: M. R. OMAR, ALL COPYRIGHTS RESERVED.
!     (c) UNIVERSITI SAINS MALAYSIA
!     (c) AGENSI NUKLEAR MALAYSIA

      module CriticalityCalculation
      
         use RandomNumberGenerator
         use NeutronBankManager
         use NeutronTrack
         use CellGeometry
         use TXSReader
         
         implicit none 
         
         ! KEFF_GUESS STORES THE INITIAL GUESS OF K-EFF.
         real :: KEFF_GUESS = 1.0
         
         ! KEFF_WGT IS THE PARTICLE WEIGHT THAT HELPS TO COMPENSATE THE PROBLEM
         ! IN WHICH THE FISSION NEUTRON  SITE COUNT  IS  NOT ENOUGH TO INITIATE
         ! NEXT NEUTRON GENERATIN CYCLE.
         real :: KEFF_WGT   = 1.0
         
      contains
         
!        ----------------------------------------------------------------------
!        THIS IS THE MAIN SUBROUTINE THAT CALCULATES THE K-EFF  OF  THE REACTOR
!        CORE.  IT SHOULD  BE  NOTED THAT  THE  USER HAS TO PROVIDE THE INITIAL
!        GUESS OF K-EFF.  THE INITIAL GUESS PROVIDED BY THE USER WILL BE STORED 
!        IN  KEFF_GUESS.  THE  USER  HAS  TO  PROVIDE  THE  NUMBER  OF PARTICLE
!        HISTORIES PER CYCLE  (WHICH  IS  ALSO KNOWN AS THE NUMBER  OF  NEUTRON 
!        SOURCE), THE NUMBER OF CYCLES AND ALSO THE SKIP COUNT.  THE SKIP COUNT
!        ALLOWS  THE CODE  TO  SKIP THE MOVING AVERAGE CALCULATION  AS  WELL AS
!        TALLYING PROCESSES TO AVOID ERRORS BEFORE THE CONVERGENCE OF K-EFF AND
!        THE FISSION SITE DISTRIBUTION ψ(n).
!        ----------------------------------------------------------------------
         subroutine CalculateKeff(pnHistories, pnCycles, pnSkip,
     &                            prKeffGuess)
         
            ! ROUTINE INPUT PARAMETERS
            integer, intent(in) :: pnHistories    ! NO. OF HISTORIES PER CYCLE.
            integer, intent(in) :: pnCycles       ! NO. OF  CALCULATION CYCLES.
            integer, intent(in) :: pnSkip         ! CALCULATION SKIP COUNT.
            real   , intent(in) :: prKeffGuess    ! INITIAL GUESS OF K-EFF
            
            ! PRIVATE VARIABLES
            integer :: iNID
            integer :: iCycle
            real    :: rKeff                ! VALUE OF K-EFF FOR CURRENT CYCLE.
            real    :: rKeffAvg             ! A MOVING AVERAGE VALUE  OF K-EFF.
            real    :: rKeffSqAvg           ! A MOVING AVERAGE VALUE OF K-EFF^2
            real    :: rKeffError           ! STANDARD ERROR OF K-EFF.
            integer :: i                    ! MULTI-PURPOSE INTEGER VAR.
            
            ! SETTING  THE  AVERAGE VALUES  TO ZERO BEFORE BEGIN MOVING AVERAGE 
            ! CALCULATION.
            rKeffAvg   = 0.0
            rKeffSqAvg = 0.0
            
            ! SETTING THE INITIAL GUESS OF K-EFF.
            KEFF_GUESS = prKeffGuess
            ! HERE WE SEGREGATE THE INITIAL GUESS OF FISSION  SRC. DISTRIBUTION
            ! BY SETTING A MONO-ENERGETIC POINT SRC. AT THE CENTRE OF EACH UNIT
            ! CELL OF THE TRIGA CORE.
            call SegregateNeutronSource(pnHistories)
            
            ! PRINTING K-EFF TABLE HEADER.
            print*
            write(*,'(A,45A1)') '  ', ('-',i=1,45)
            write(*,'(A10,2X,A15,2X,A15)') 'CYCLE',
     &         'K-EFF', ' AVERAGE K-EFF'
            write(*,'(A,45A1)') '  ', ('-',i=1,45)
            ! PRINT THE INITIAL GUESS OF K-EFF (CYCLE-0)
            write(*,'(I10,2X,F15.5,2X,A15)') 
     &             0, KEFF_GUESS, 'SKIP'
     
            ! BEGIN THE CYCLE LOOP, THE NUMBER OF CYCLES IS pnCycles.
            do iCycle=1, pnCycles, 1     
            
            ! WE CLEAN-UP THE FISSION BANK.  THIS WILL ALSO  RESETS  THE NO. OF
            ! FISSION NEUTRONS IN FISSION BANK (N_KEFFSRC) TO ZERO.
            call ClearFissionBank()
            
            ! HERE WE LOOP  TO  SIMULATE ALL NEUTRON SOURCES.  WE TRACK NEUTRON
            ! HISTORIES  TO  OBTAIN  ψ(n+1).  iNID   IS THE  NEUTRON ID THAT IS 
            ! CURRENTLY BEING TRACKED. WE SET THE INITIAL  NEUTRON ID TO 1.  
            ! THE RULE OF THUMB: 0 < NID <= N_KEFFHIS.
            iNID = 0
            
            ! WHILE NOT ALL NEUTRONS IN THE MAIN BANK ARE TRACKED.
            do while (iNID .lt. N_KEFFHIS)

               iNID = iNID + 1
             
               ! BEGIN TRACK NEUTRON ID = iNID.  THIS WILL SET THE NEUTRON BANK 
               ! ARRAYS CURSOR TO THE DESIRED NEUTRON.
               call BeginTracking(iNID)
               
               ! BEGIN NEUTRON  HISTORY  SIMULATION UNTIL IT IS DEAD OR ESCAPED.
               ! IF THE NEUTRON IS STILL ALIVE AND INSIDE THE REACTOR, ITS ISTAT
               ! VALUE IS ZERO (0) (N_WAIT). CELL ID EQUALS TO -2 INDICATES THAT
               ! THE NEUTRON IS OUTSIDE THE REACTOR CORE.
               do while((C_CEL(iNID,1) .ne. -2) .and. 
     &                  (ISTAT(iNID,1) .eq. N_WAIT))
     
                  ! HERE WE USED A SPECIAL TRANSPORT CODE DESIGNED FOR KEFF CALCU
                  ! LATION. PROCEED TO KeffTransport() SUBROUTINE FOR MORE INFO.
                  call KeffTransport()
                  
               enddo
               
               ! IF THE NEUTRON HAS ESCAPED FROM THE REACTOR, WE SET ITS STATUS
               ! TO N_LOST (=2).
               if(ISTAT(iNID,1) .eq. N_WAIT) ISTAT(iNID,1) = N_LOST
              
            ! END HISTORY ITERATIONS FOR A NEUTRON SOURCE.   
            enddo
            
            ! WE ESTIMATE CURRENT K-EFF USING THE FORMULA
            !
            !              NO. OF NEUTRONS IN GENERATION (N+1)     N_KEFFSRC
            !    K-EFF = -------------------------------------- =  ---------
            !               NO. OF NEUTRONS IN GENERATION (N)      N_KEFFHIS
            !
            rKeff = real(N_KEFFSRC) / real(N_KEFFHIS)
            
            ! WE SET THE GUESS K-EFF TO CURRENT K-EFF VALUE.
            KEFF_GUESS = rKeff
            
            ! CALCULATE MOVING AVERAGES
            if(iCycle .gt. pnSkip) then
               rKeffAvg = (1.0 - 1.0 / real(iCycle-pnSkip)) * rKeffAvg +
     &                    rKeff / real(iCycle-pnSkip)
               rKeffSqAvg = (1.0 - 1.0 / real(iCycle-pnSkip)) 
     &                     * rKeffSqAvg +
     &                    rKeff**2 / real(iCycle-pnSkip)
            endif
            
            ! IF CUrRENT CYCLE IS LESS THAN SKIP COUNT,  WE IGNORE PRINTING THE 
            ! MOVING AVERAGE OF K-EFF.
            if(iCycle .le. pnSkip) then
               write(*,'(I10,2X,F15.5,2X,A15,I10)') 
     &             iCycle, rKeff, 'SKIP', N_KEFFSRC  
            else
               write(*,'(I10,2X,F15.5,2X,F15.5,I10)') 
     &             iCycle, rKeff, rKeffAvg, N_KEFFSRC
            endif
            ! NOW WE SET FISSION NEUTRONS AT THE COLLISION SITES AS THE NEUTRON
            ! SOURCE OF THE NEXT CRITICALITY CALCULATION CYCLE.
            KEFF_WGT = real(N_KEFFHIS) / real(N_KEFFSRC)
            call SetNextGenerationSource()
            
            ! END OF CYCLE LOOP
            enddo
            
            ! WE PRINT THE AVERAGE VALUE OF K-EFF WITH ITS ERROR.
            rKeffError = sqrt((rKeffSqAvg - rKeffAvg**2) / 
     &               real(iCycle-pnSkip-1))
            write(*,'(A,F7.5,A,F7.5,A)') 'Estimated k-eff : (', 
     &         rKeffAvg, ' +/- ', rKeffError, ')'
     
            return
            
         end subroutine
         
!        ----------------------------------------------------------------------
!        THIS IS THE  TRANSPORT SUBROUTINE  MODIFIED FOR K-EFF CALCULATION. THE
!        DIFFERENCE  IS THAT WHEN NEUTRON UNDERGOES FISSION, THE CHILD NEUTRONS
!        ARE BANKED FOR LATER USE OF FISSION SOURCE OF  NEXT CYCLE.  WHEREAS IN 
!        THE STANDARD  FIXED  NEUTRON  SOURCE,  THE CHILD NEUTRONS HISTORY  ARE
!        IMMEDIATELY TRACKED AFTER BORN.
!        ----------------------------------------------------------------------
         subroutine KeffTransport()
            
            ! PRIVATE VARIABLES
            real :: rDistCol        ! DISTANCE TO THE NEXT NEUTRON COLLISION.
            real :: rDistBoundary   ! NEAREST DISTANCE TO THE CELL BOUNDARY.
            real :: rX, rY, rZ      ! NEUTRON POSITION (x,y,z).
            real :: rU, rV, rW      ! NEUTRON DIRECTION Ω = (u,v,w).
            
            ! AT THE MOMENT iStatus IS JUST A DUMMY VARIABLE STORING THE STATUS
            ! OF NEUTRON REACTION OPERATIONS.
            integer :: iStatus
            
            ! WE SAMPLE THE DISTANCE TO THE NEXT NEUTRON COLLISION.
            rDistCol = DistanceToNextCollision()
            
            ! WE OBTAIN THE DISTANCE TO THE NEAREST CELL BOUNDARY.
            rDistBoundary = DistanceToNearestBoundary()
            
            ! WE OBTAIN THE CURRENT NEUTRON POSITION IN THE REACTOR.
            rX = C_POS(TRACK_NID, 1)
            rY = C_POS(TRACK_NID, 2)
            rZ = C_POS(TRACK_NID, 3)
            rU = C_DIR(TRACK_NID, 1)
            rV = C_DIR(TRACK_NID, 2)
            rW = C_DIR(TRACK_NID, 3)
            
            
            if(rDistCol .gt. rDistBoundary) then 
               ! THE DISTANCE TO THE NEAREST BOUNDARY IS LESS THAN THE DISTANCE 
               ! TO THE  NEXT COLLISION,  WE TRANSPORT  THE NEUTRON TO THE CELL 
               ! BOUNDARY. r(k+1) = r(k) + D Ω
               rX = rX + rDistBoundary * rU
               rY = rY + rDistBoundary * rV
               rZ = rZ + rDistBoundary * rW
            else
               ! IF COLLISION IS POSSIBLE, WE TRANSPORT THE NEUTRON TOWARDS THE
               ! COLLISION SITE. r(k+1) = r(k) + D_col Ω
               rX = rX + rDistCol * rU
               rY = rY + rDistCol * rV
               rZ = rZ + rDistCol * rW
               
               ! IF COLLISION IS POSSIBLE WE FIRST  SELECT THE NEUTRON REACTION
               ! WE SEE WHETHER THE NEUTRONS WILL UNDERGO SCATTERING OR ABSORP-
               ! TION.
               if(WillScatter() .eqv. .true.) then
                  call ScatterNeutron()
               else
                  ! IF ABSORPTION IS POSSIBLE, WE CHECK WHETHER THE NEUTRON IS 
                  ! ABSORBED IN FUEL OR BEING CAPTURED. 
                  if(WillFission() .eqv. .true.) then
                     call KeffFission()
                  else 
                     call AbsorbNeutron()
                  endif
               endif
            endif
            
            ! FINALLY WE REGISTER  THE LOCATION  OF THE NEUTRON TO THE  NEUTRON
            ! NEUTRON BANK RECORD.
            call SetNeutronPosition(TRACK_NID, rX, rY, rZ, iStatus)
            
            return
            
         end subroutine
         
!        ----------------------------------------------------------------------         
!        THIS SUBROUTINE DEFINES  THE FISSION PROCESS FOR K-EFF CALCULATION. IN 
!        K-EFF CALCULATION, ON A FISSION REACTION IS HAPPENING AT THE COLLISION
!        SITE, THE FISSION NEUTRONS PRODUCED ARE BANKED  IN THE FISSION NEUTRON
!        BANK. THE FISSION SITE  IS  ALSO SAVED TO  THE NEUTRON BANK SO THAT IT 
!        CAN BE USED  TO SAMPLE THE NEUTRON FISSION SOURCE DISTRIBUTION FOR THE
!        NEXT CYCLE.
!        ----------------------------------------------------------------------
         subroutine KeffFission()
         
            ! rP IS A TABLE THAT CONSISTS OF THE CUMMULATIVE DISTRIBUTION FUNC.
            ! OF THE FISSION NEUTRON ENERGY GROUP. WE SELECT THE OUTGOING ENER-
            ! GY GROUP OF FISSION NEUTRON BY INVERSING THE CUMMULATIVE DIST.
            real :: rP(TXS_GRP+1) 
            
            ! CURRENT ENERGY GROUP OF THE INCIDENT NEUTRON (CURRENT NEUTRON).
            integer :: iGroup
            
            ! A REAL VARIABLE THAT STORES A RANDOM NUMBER.
            real    :: rEpsilon
            
            ! AN INTEGER VARIABLE THAT STORES THE OUTGOING ENERGY GROUP OF
            ! FISSION NEUTRONS.
            integer :: iFissionEnergyGroup
            
            ! A DUMMY VARIABLE THAT STORES THE STATUS OF A NEUTRON REACTION.
            integer :: iStatus
            
            ! A  RANDOM ANGLE [0,2π)  THAT IS USED  WHEN SAMPLING  THE OUTGOING
            ! DIRECTION Ω(u,v,w) OF THE FISSION NEUTRONS.
            real    :: rRandomPhi
            
            ! A SCATTERING COSINE, μ, USED TO SAMPLE  THE OUTGOING ANGLE OF THE
            ! FISSION NEUTRON.
            real    :: rMu
            
            ! INTEGERS THAT DEFINE THE CELL INDICES.
            integer :: iCell, iLayer
            
            ! MULTI-PURPOSE INTEGERS USED FOR ITERATIONS.
            integer :: j, k
            
            ! REAL VARIABLES THAT STORE CURRENT INCIDENT NEUTRON DIRECTION.
            real    :: rU, rV, rW
            
            ! REAL VARIABLES THAT STORE  THE OUTGOING  NEUTRON  DIRECTION FOR A
            ! NEWBORN FISSION NEUTRON.
            real    :: rUNew, rVNew, rWNew
            
            real, parameter :: PI = 3.141592654
            
            ! VARS. USED FOR SAMPLING  THE NUMBER OF OUTGOING FISSION NEUTRONS.
            real    :: rINu, rNu, iNu

            real    :: rYield 
            ! HERE WE OBTAIN ALL OF THE CURRENT NEUTRON IDENTITY.
            iGroup = C_GRP(TRACK_NID, 1)
            iCell  = C_CEL(TRACK_NID, 1)
            iLayer = C_CEL(TRACK_NID, 2)
            rU     = C_DIR(TRACK_NID, 1)
            rV     = C_DIR(TRACK_NID, 2)
            rW     = C_DIR(TRACK_NID, 3)
            rYield = TXS_TABLE(iLayer, iCell, iGroup, 3)
            ! PREPARE THE FISSION NEUTRON COUNT. WE TOSS A RANDOM NUMBER.
            rEpsilon = Rnd(int(TRACK_NID, 8))
            
            ! CALC. THE NUMBER  OF  FISSION NEUTRONS  ACCORDING TO THE INCOMING
            ! NEUTRON ENERGY.  iNu IS THE INTEGER VERSION OF THE rNu. WE OBTAIN
            ! THE AVG. NUMBER OF NEUTRONS (rNu)  FROM THE NUBAR-ENERGY CORRELA-
            ! TION. THIS CAN BE DONE BY USING THE GetNuBar FUNCTION.  SEE FUNC.
            ! DEFINITION IN NeutronInteractions MODULE.
            rNu  = GetNuBar()
            
            ! NEXT WE SAMPLE THE INTEGER VERSION OF AVG. NUMBER OF FISSION NEU-
            ! TRONS RELEASED PER FISSION.
            rINu = floor(rNu)
            if(rEpsilon .le. (rNu - rINu)) then
               iNu = int(rINu) + 1
            elseif(rEpsilon .gt. (rNu - rINu)) then
               iNu = int(rINu)
            endif
            
            ! THIS PART IS IMPORTANT TO PREVENT THE SHORTAGE OF NEW FISSION SRC
            ! SITES FOR THE USE OF NEXT CYCLE. THIS IS THE WEIGHTED SAMPLING OF
            ! THE NUMBER OF NEUTRONS TO BE RELEASE DURING  FISSION  REACTION AT 
            ! THE COLLISION SITE.
            !                
            !               _    Σf   1
            !         ν = W v  ----- ---- + ε
            !                   Σtot k(n) 
            !
            ! SEE LOS ALAMOS MCNP MANUAL la-ur-03-1987 (PUBLIC), PG 2-168.
            
            rNu = KEFF_WGT * (rYield / GetSigTot()) 
     &            * (1.0 / KEFF_GUESS) +
     &             Rnd(int(TRACK_NID,8))
c            rNu = KEFF_WGT * real(iNu) / KEFF_GUESS
            iNu = ceiling(rNu)


            ! WE EXIT THIS SUBROUTINE OF THERE IS NO FISSION NEUTRON IS PRODUCED.
            if(iNu .eq. 0) return
            
            ! FOR EACH FISSION NEUTRONS PRODUCED, WE SAMPLE ITS OUTGOING ENERGY
            ! AND ITS OUTGOING DIRECTION. THEN WE STORE IT IN THE FISSION BANK.
            do k=1, iNu, 1
            
            ! **PART I: SAMPLING THE OUTGOING ENERGY OF THE FISSION NEUTRON.
            
            ! FIRST WE TOSS FOR A RANDOM NUMBER.
            rEpsilon = Rnd(int(TRACK_NID, 8))
            
            ! NEXT WE PREPARE THE CUMMULATIVE DISTRIBUTION FUNC. OF THE OUTGOING
            ! NEUTRON ENERGY GROUP.
            rP(1) = 0.0
            do j=1, TXS_GRP, 1
               rP(j+1) = TXS_TABLE(iLayer, iCell, j, 4) + rP(j)
            enddo
            
            ! THEN USING THE PREPARED RANDOM NUMBER,  rEpsilon,  WE USE THE INVERSE 
            ! METHOD TO SAMPLE THE OUTGOING ENERGY GROUP FROM THE CUMMULATIVE DIST.
            do j=1, TXS_GRP, 1
               if((rEpsilon .ge. rP(j)) .and.
     &            (rEpsilon .lt. rP(j+1))) then
                  iFissionEnergyGroup = j
               endif
            enddo
            
            ! **PART II: SAMPLING THE OUTGOING DIRECTION OF THE FISSION NEUTRON.
            ! FIRST WE SAMPLE THE SCATTERING COSINE, μ ∈ (-1, 1).
            rMu = 2.0 * Rnd(int(TRACK_NID,8))  - 1.0
            
            ! NEXT WE SAMPLE A RANDOM NUMBER, ϕ ∈ [0,2π).
            rRandomPhi = 2.0 * PI * Rnd(int(TRACK_NID,8))      
            
            ! NEXT WE CALCULATE THE OUTGOING NEUTRON DIRECTION USING (μ, ϕ).
            rUNew = rMu * rU + (sqrt(1.0 - rMu**2) * 
     &              (rU * rW * cos(rRandomPhi) - rV * sin(rRandomPhi))) 
     &              / sqrt(1.0 - rW**2)

            rVNew = rMu * rV + (sqrt(1.0 - rMu**2) * 
     &              (rV * rW * cos(rRandomPhi) + rU * sin(rRandomPhi))) 
     &              / sqrt(1. - rW**2)     

            rWNew = rMu * rW - sqrt(1.0 - rMu**2) *
     &              sqrt(1. - rW**2) * cos(rRandomPhi)
     
            ! THEN WE ADD THE FISSION NEUTRON INTO  THE FISSION BANK SO THAT WE
            ! CAN USE IT LATER AS A SOURCE FOR NEXT CYCLE.
            call AddFissionNeutron(C_POS(TRACK_NID, 1),
     &                               C_POS(TRACK_NID, 2),
     &                               C_POS(TRACK_NID, 3),
     &                               rUNew, rVNew, rWNew,
     &                               iFissionEnergyGroup, 1.0)
            
            enddo
            ! END LOOP FOR EACH NEWBORN FISSION NEUTRON

            call KillNeutron(TXS_GRP, N_FISS, iStatus)
            return
         end subroutine         
         
!        ----------------------------------------------------------------------        
!        THIS SUBROUTINE RANDOMLY SELECTS FISSION SOURCES FROM THE FISSION NEU-
!        TRON BANK.
!        ----------------------------------------------------------------------
         subroutine SetNextGenerationSource()
            integer :: i, j, n
            real    :: rEpsilon
            integer :: iRandomNID
            
            ! HERE WE SELECT A NEUTRON FOR N_KEFFHIS TIMES,SO THAT THE NEXT GE-
            ! NERATION HAS THE SAME AMOUNT OF NEUTRON SOURCE. THE SELECTED NEU-
            ! TRON IS THEN TRANSFERED TO THE MAIN NEUTRON BANK.
            do n=1, N_KEFFHIS, 1
            
               ! HERE WE RANDOMLY SHUFFLE THE NEUTRON SELECTIONS.
1              i = int(ceiling(N_KEFFSRC * Rnd(int(N_KEFFHIS+1,8))))

               ! IN CASE IF THE COMPUTER GONE CRAZY  AND  ACCIDENTALLY SELECTS 
               ! THE WRONG NEUTRON ID, WE RE-SAMPLE THE NEUTRON ID AGAIN. LOL.
               if(i .gt. N_KEFFHIS) goto 1
               
               ! HERE WE  TRANSFER ALL OF THE FISSION NEUTRON INFORMATION INTO
               ! THE MAIN NEUTRON BANK.
               C_DIR(n, 1) = FC_DIR(i, 1)
               C_DIR(n, 2) = FC_DIR(i, 2)
               C_DIR(n, 3) = FC_DIR(i, 3)
               L_DIR(n, 1) = FL_DIR(i, 1)
               L_DIR(n, 2) = FL_DIR(i, 2)
               L_DIR(n, 3) = FL_DIR(i, 3)
               C_POS(n, 1) = FC_POS(i, 1)
               C_POS(n, 2) = FC_POS(i, 2)
               C_POS(n, 3) = FC_POS(i, 3)
               L_POS(n, 1) = FL_POS(i, 1)
               L_POS(n, 2) = FL_POS(i, 2)
               L_POS(n, 3) = FL_POS(i, 3)
               B_POS(n, 1) = FB_POS(i, 1)
               B_POS(n, 2) = FB_POS(i, 2)
               B_POS(n, 3) = FB_POS(i, 3)
               C_CEL(n, 1) = FC_CEL(i, 1)
               C_CEL(n, 2) = FC_CEL(i, 2)
               L_CEL(n, 1) = FL_CEL(i, 1)
               L_CEL(n, 2) = FL_CEL(i, 2)
               B_CEL(n, 1) = FB_CEL(i, 1)
               B_CEL(n, 2) = FB_CEL(i, 2)
               C_GRP(n, 1) = FC_GRP(i, 1)
               NCHLD(n, 1) = FNCHLD(i, 1)
               PARNT(n, 1) = FPARNT(i, 1)
               ISTAT(n, 1) = N_WAIT   
               C_WGT(n, 1) = 1.0
               
              
            enddo

            return
            
         end subroutine

!        ----------------------------------------------------------------------         
!        THIS SUBROUTINE PREPARES THE CELL CENTER POINT FOR INITIAL SOURCE DIST
!        PREDICTION.
!        ----------------------------------------------------------------------
         subroutine CellCentrePoint(piCellID, piLayerID, prX, prY, prZ)
            integer, intent(in)  :: piCellID
            integer, intent(in)  :: piLayerID
            real   , intent(out) :: prX, prY, prZ
            
            ! REACTOR CORE CELL
            integer :: iI, iK, iL
            ! RING OUTER RADIUS, R(i) AND R(i-1), OF THE CELL.
            real    :: rRi, rRi1
            ! ANGLE(THETA) OF THE CELL, θ(k) AND θ(k-1)
            real    :: rThk, rThk1
            
            ! WE PREPARE THE CELL INDICES (i,k,l) = (iI,iK,iL).
            call CellIndex(piCellID, iI, iK)
            iL    = piLayerID
            
            ! WE PREPARE THE VALUES OF RADIUS AND  THETA REQUIRED FOR  THE CALC
            ! OF THE CELL CENTRE POINT.
            rRi   = RingRadius(iI)                 ! R(i)
            rRi1  = RingRadius(iI-1)               ! R(i-1)
            rThk  = ThetaFromIndex(iI, iK)         ! θ(k)
            rThk1 = ThetaFromIndex(iI, iK-1)       ! θ(k-1)
            
            ! X = -1/2[ R(i) + R(i-1) ] cos[ 1/2[ θ(k) + θ(k-1) ]
            prX = -0.5 * (rRi + rRi1) *
     &             cos( 0.5 * (rThk + rThk1) )
            ! Y =  1/2[ R(i) + R(i-1) ] sin[ 1/2[ θ(k) + θ(k-1) ]
            prY =  0.5 * (rRi + rRi1) *
     &             sin( 0.5 * (rThk + rThk1) )
            ! HERE WE OBTAIN THE Z POINT INTUITIVELY.
            if(iL .eq. 1) then
               prZ = 0.5 * LayerHeight(iL)
            else
               prZ = 0.5 * (LayerHeight(iL) + LayerHeight(iL-1))
            endif
            
            return
            
         end subroutine

!        ----------------------------------------------------------------------
!        THIS SUBROUTINE  ENABLES  US TO SAMPLE ISOTROPIC NEUTRON DIRECTION. WE 
!        WILL  BE  USING THIS SUB TO GENERATE INITIAL SOURCE NEUTRON DIRECTION. 
!        THE CALCULATION IS BASED  ON  MONTE CARLO  METHODS  FOR RAD. TRANSPORT
!        FUNDAMENTALS AND ADVANCED TOPICS, OLEG N. VASSILIEV, PAGE 40.
!        ----------------------------------------------------------------------
         subroutine SampleIsotropicDirection(prU, prV, prW)
            real, intent(out) :: prU, prV, prW
            
            real :: rEpsilon
            real :: rMu, rPhi
            real, parameter :: PI = 3.141592654
            
            rEpsilon = Rnd(int(TRACK_NID,8))
            rMu      = 2.0 * rEpsilon - 1.0
            rEpsilon = Rnd(int(TRACK_NID,8))
            rPhi     = 2.0 * PI * rEpsilon
            
            prU = sqrt(1.0 - rMu**2) * cos(rPhi)
            prV = sqrt(1.0 - rMu**2) * sin(rPhi)
            prW = rMu
            
            return
         end subroutine

!        ----------------------------------------------------------------------
!        THIS SUBROUTINE SEGREGATES  THE  INITIAL SOURCE DISTRIBUTION Ψ(0) PRO-
!        VIDED BEFORE STARTING THE FIRST  CYCLE OF K-EFF CALCULATION. THE DIST.
!        IS ASSUMED TO BE SEGREGATED  AT  THE CENTRE POINT OF EACH REACTOR CORE
!        CELLS, AND THE INITIAL NEUTRON SOURCE IS ISOTROPIC AND HAS A MONOGROUP
!        ENERGY. FORTUNATELY, REGARDLES OF  THE  INITIAL ENERGY AND DIRECTIN OF
!        THE NEUTRON SOURCE, THE NEUTRON SOURCE DISTRIBUTION  WILL CONVERGES TO 
!        THE ACTUAL FISSION DISTRIBUTION  CHARACTERISED BY THE REACTOR GEOMETRY
!        AND MATERIALS.   
!        ----------------------------------------------------------------------      
         subroutine SegregateNeutronSource(pnHistories)
            ! THE NUMBER OF NEUTRON HISTORIES PER CYCLE.
            integer, intent(in) :: pnHistories
            
            ! THE MAXIMUM VALUE OF CELL ID.
            integer :: nMaxCellID
            
            ! THE MAMXIMUM VALUE OF LAYER ID.
            integer :: nMaxLayerID
            
            ! INTEGER INDICES STORING THE CURRENT CELL ID AND LAYER ID.
            integer :: iCell, iLayer
            
            ! THE POSITION OF THE  CELL  CENTRE  POINT (ALSO THE INITIAL SOURCE
            ! POINT).
            real    :: rX, rY, rZ
            
            ! THE INITIAL DIRECTION OF THE NEUTRON SRC.  THIS  IS ASSUMED TO BE
            ! ISOTROPIC.
            real    :: rU, rV, rW
            
            ! THE NUMBER OF NEUTRONS TO BE EJECTED PER UNIT CELL.
            integer :: nNeutronPerCell
            
            integer :: i
            
            ! FIRST WE OBTAIN THE MAXIMUM NUMBER OF CELL ID. IT IS A FIXED VALUE.
            if(RingCount() .eq. 6) nMaxCellID = 91
            if(RingCount() .eq. 7) nMaxCellID = 127
            
            ! THEN WE OBTAIN THE MAXIMUM LAYER ID VALUE FROM THE GEOMETRY MODULE.
            ! WE OBTAIN THIS FROM THE GEOMETRY MODULE BECAUSE THE USER CAN SET 
            ! THE MAXIMUM NUMBER OF CORE LAYERS TO ANY  INTEGER VALUE > 0.
            nMaxLayerID = LayerCount()
            
            ! HERE WE CALCULATE THE NUMBER OF NEUTRONS TO BE EJECTED PER CELL SO
            ! THAT WE CAN ADJUST THE NUMBER OF NEUTRON SOURCE TO BE THE SAME NUM
            ! BER OF NEUTRON HISTORIES REQUESTED BY THE CODE USER.
            nNeutronPerCell = int(pnHistories /
     &                       (nMaxCellID * nMaxLayerID))
     
            ! WE SET THE HISTORY COUNT.
            call SetKeffHistoryCount(nMaxCellID * nMaxLayerID * 
     &            nNeutronPerCell)
            
            ! WE BEGIN SEGREGATING THE INITIAL NEUTRON SOURCE AND STORE THEM INT 
            ! THE PRIMARY NEUTRON BANK.
            do iLayer=1, nMaxLayerID, 1
               do iCell=1, nMaxCellID, 1
                  call CellCentrePoint(iCell, iLayer, rX, rY, rZ)
                  do i=1, nNeutronPerCell, 1
                     call SampleIsotropicDirection(rU, rV, rW)
                     call AddPrimaryNeutron(rX, rY, rZ, 
     &                                      rU, rV, rW,
     &                                      5, 1.0)
                  enddo         
               enddo
            enddo
            
            return
            
         end subroutine
         
         
      end module
