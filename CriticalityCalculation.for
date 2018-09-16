!                                                                       
!     Eigenvalue Calculation Module, Revision 180915-1.
!     Author M. R. Omar, October 2017. All copyrights reserved, 2017.
!     (c) Universiti Sains Malaysia
!     (c) Malaysian Nuclear Agency
!
!     NOTICE:  All information contained herein is, and remains the pro-
!     perty of the copyright owners  and  their suppliers,  if any.  The 
!     intellectual and technical concepts contained herein are  proprie-
!     tary to the copyright owners and their suppliers and are protected
!     by  trade  secret  or  copyright  law.    Dissemination  of   this 
!     source  code  or  reproduction  of  this  source  code is strictly
!     forbidden  unless  prior written permission  is  obtained from the 
!     owners.
!
!     This source code was created on 1/11/2017 11:15 PM by M. R. Omar.
!     Last revision date 15/9/2018.
!
!     THIS MODULE CONTAINS THE NECESSARY METHODS THAT ENABLE THE CODE USER AND
!     PROGRAMMERS TO DO CRITICALITY CALCULATION OF A TRIGA REACTOR CORE.
!

      module CriticalityCalculation
      
         use RandomNumberGenerator
         use NeutronBankManager
         use NeutronTrack
         use CellGeometry
         use TXSReader
         use Tally
         
         implicit none 
         
         
         ! KEFF_GUESS STORES THE INITIAL GUESS OF K-EFF.
         real :: KEFF_GUESS = 1.0
         real :: NOMINAL_POWER = 1000.0
         real :: OUTLIER_CONTROL = 0.020
         ! KEFF_WGT IS THE PARTICLE WEIGHT THAT HELPS TO COMPENSATE THE PROBLEM
         ! IN WHICH THE FISSION NEUTRON  SITE COUNT  IS  NOT ENOUGH TO INITIATE
         ! NEXT NEUTRON GENERATION CYCLE.
         real    :: KEFF_WGT     = 1.0
         ! KEFF_CYCLE INDICATES THE CURRENT CALCULATION CYCLE
         integer :: KEFF_CYCLE   = 1
         ! TOTAL_CYCLE STORES THE TOTAL NUMBER OF CYCLE SET BY THE USER. THIS
         ! SHARED VARIABLE WILL BE SET OUTSIDE THIS MODULE.
         integer :: TOTAL_CYCLE  = 1
         ! INITIAL_HSRC IS THE INITIAL VALUE OF SOURCE SHANNON ENTROPY FOR CAL-
         ! CULATION.
         real    :: INITIAL_HSRC = 0.0
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
     &                            prKeffGuess,pcRunID)
         
            ! ROUTINE INPUT PARAMETERS
            integer, intent(in) :: pnHistories    ! NO. OF HISTORIES PER CYCLE.
            integer, intent(in) :: pnCycles       ! NO. OF  CALCULATION CYCLES.
            integer, intent(in) :: pnSkip         ! CALCULATION SKIP COUNT.
            real   , intent(in) :: prKeffGuess    ! INITIAL GUESS OF K-EFF
            character(len=30), intent(in) :: pcRunID            
            ! PRIVATE VARIABLES
            integer :: iNID
            integer :: iCycle
            real    :: rKeff                ! VALUE OF K-EFF FOR CURRENT CYCLE.
            real    :: rKeffAvg             ! A MOVING AVERAGE VALUE  OF K-EFF.
            real    :: rKeffSqAvg           ! A MOVING AVERAGE VALUE OF K-EFF^2
            real    :: rKeffError           ! STANDARD ERROR OF K-EFF.
            integer :: i, j, k, nTotCycles          ! MULTI-PURPOSE INTEGER VAR.
            integer :: NFO = 24
            integer :: NFT = 26
            character(len=40) :: cNumText
            real    :: rCPUStart, rCPUFinish
            integer :: iCPUMinutes, iCPUSeconds, iCPUMinTot, iCPUSecTot
            real    :: rX, rX2, rxBar, rx2Bar, rSBar, rRErr, rFOM
            real    :: rCurrentError
            character(len=80) :: caKeffLine(1000000)

            ! SETTING  THE  AVERAGE VALUES  TO ZERO BEFORE BEGIN KEFF AVERAGING. 
            ! CALCULATION.
            rKeffAvg   = 0.0
            rKeffSqAvg = 0.0
            TOTAL_CYCLE = pnCycles
            iCPUMinTot = 0
            iCPUSecTot = 0
            ! OPEN FILE FOR WRITING THE K-EFF VALUE FOR ALL CYCLES.
            open(unit=NFO,file='KEFF.OUT',status='unknown')
            write(NFO,'(2A)') '  EFFECTIVE MULTIPLICATION ',
     &                       'FACTOR CALCULATION RESULT'
            write(NFO,'(2A)') '  =========================',
     &                       '========================='
            write(NFO, '(6A)')  '  Run ID                        = ',
     &          pcRunID
            write(NFO,'(A,I0)') '  Total number of cycles        = ',
     &          pnCycles
            write(NFO,'(A,I0)') '  Total number of active cycles = ',
     &          pnCycles-pnSkip
            write(NFO,'(A,I0)') '  Total number of histories     = ',
     &          pnHistories
            write(NFO,'(A,I0)') '  Total number of energy groups = ',
     &          TXS_GRP
            write(NFO,'(A,I0)') '  Number of cells per layer     = ',
     &          TXS_CEL
            write(NFO,'(A,I0)') '  Number of layers              = ',
     &          TXS_LAY
     
            write(NFO,'(A)') ''
            write(NFO,'(A,58A1)') '  ', ('-',i=1,58)
            write(NFO,'(A10,2X,A10,2X,2A10,A13)') 'CYCLE',
     &         'K-EFF', ' AVG K-EFF', 'H-SRC', 'CPU TIME'
            write(NFO,'(A,58A1)') '  ', ('-',i=1,58)
            ! OPEN FILE FOR WRITING NEUTRON TRACKS OF THE FIRST 1000 HISTORIES.
            open(unit=NFT,file='TRACK.OUT',status='unknown')
            ! SETTING THE INITIAL GUESS OF K-EFF.
            KEFF_GUESS = prKeffGuess
            rKeff = KEFF_GUESS
            nTotCycles = pnCycles
            ! HERE WE SEGREGATE THE INITIAL GUESS OF FISSION  SRC. DISTRIBUTION
            ! BY SETTING A MONO-ENERGETIC POINT SRC. AT THE CENTRE OF EACH UNIT
            ! CELL OF THE TRIGA CORE.
            call SegregateNeutronSource(pnHistories)
            
            ! PRINTING K-EFF TABLE HEADER.
            print*
            write(*,'(A,67A1)') '  ', ('-',i=1,67)
            write(*,'(A10,2X,A10,2X,3A10,A13)') 'CYCLE',
     &         'K-EFF', ' AVG K-EFF', 'TALLY', 'H-SRC', 'CPU TIME'
            write(*,'(A,67A1)') '  ', ('-',i=1,67)
            ! PRINT THE INITIAL GUESS OF K-EFF (CYCLE-0)
            write(*,'(I10,2X,F10.5,2X,A10,A10,F10.3)') 
     &            0, KEFF_GUESS, 'SKIP', 'GUESS', -1.0
     &            * log(1.0/real(TXS_LAY*TXS_CEL))/log(2.0)
            write(NFO,'(I10,2X,F10.5,2X,F10.5,2X,F10.5)') 0, KEFF_GUESS,
     &            0.0, -1.0
     &            * log(1.0/real(TXS_LAY*TXS_CEL))/log(2.0)
            ! BEGIN THE CYCLE LOOP, THE NUMBER OF CYCLES IS pnCycles.
            do iCycle=1, nTotCycles, 1     
            
!-----------HERE WE START THE CPU TIME COUNTER
            call cpu_time(rCPUStart)
            
            ! WE CLEAN-UP THE FISSION BANK.  THIS WILL ALSO  RESETS  THE NO. OF
            ! FISSION NEUTRONS IN FISSION BANK (N_KEFFSRC) TO ZERO.
            call ClearFissionBank()
312         call ResetTally()
            
            ! HERE WE LOOP  TO  SIMULATE ALL NEUTRON SOURCES.  WE TRACK NEUTRON
            ! HISTORIES  TO  OBTAIN  ψ(n+1).  iNID   IS THE  NEUTRON ID THAT IS 
            ! CURRENTLY BEING TRACKED. WE SET THE INITIAL  NEUTRON ID TO 1.  
            ! THE RULE OF THUMB: 0 < NID <= N_KEFFHIS.
            iNID = 0
            
            
            ! WHILE NOT ALL NEUTRONS IN THE MAIN BANK ARE TRACKED.
            !$OMP PARALLEL DO
            do while (iNID .lt. N_KEFFHIS)

               iNID = iNID + 1
             
               ! BEGIN TRACK NEUTRON ID = iNID.  THIS WILL SET THE NEUTRON BANK 
               ! ARRAYS CURSOR TO THE DESIRED NEUTRON.
               call BeginTracking(iNID)
               
               ! SHOW USERS CURRENT (NEUTRON-ID)-OF-(NUMBER OF HISTORIES).
               if(mod(iNID,100) .eq. 0) then
                  write(*,'(A,I0,A,I0,A,I0,A)',advance='no') 
     &            '  Simulating neutron ', iNID, ' of ', 
     &             N_KEFFHIS, '. [', 100*iNID/N_KEFFHIS, '%]'//char(13)
               endif
               ! BEGIN NEUTRON  HISTORY  SIMULATION UNTIL IT IS DEAD OR ESCAPED.
               ! IF THE NEUTRON IS STILL ALIVE AND INSIDE THE REACTOR, ITS ISTAT
               ! VALUE IS ZERO (0) (N_WAIT). CELL ID EQUALS TO -2 INDICATES THAT
               ! THE NEUTRON IS OUTSIDE THE REACTOR CORE.
               !$OMP PARALLEL DO
               do while((C_CEL(iNID,1) .ne. -2) .and.
     &                  (C_CEL(iNID,2) .ne. -2) .and.
     &                  (ISTAT(iNID,1) .eq. N_WAIT))
     
                  ! HERE WE USED A SPECIAL TRANSPORT CODE DESIGNED FOR KEFF CALCU
                  ! LATION. PROCEED TO KeffTransport() SUBROUTINE FOR MORE INFO.
                  call KeffTransport()
                  if(iCycle .eq. pnCycles) then
                     if(iNID .le. 1000) then
                        write(NFT, '(I10,3F20.3)') iNID, C_POS(iNID,1),
     &                     C_POS(iNID,2), C_POS(iNID,3)
                     endif
                  endif
                  
               enddo
               !$OMP END PARALLEL DO
               ! IF THE NEUTRON HAS ESCAPED FROM THE REACTOR, WE SET ITS STATUS
               ! TO N_LOST (=2).
               if(ISTAT(iNID,1) .eq. N_WAIT) ISTAT(iNID,1) = N_LOST
              
            ! END HISTORY ITERATIONS FOR A NEUTRON SOURCE.   
            enddo
            !$OMP END PARALLEL DO
            
            ! WE ESTIMATE CURRENT K-EFF USING THE FORMULA
            !
            !              NO. OF NEUTRONS IN GENERATION (N+1)     N_KEFFSRC
            !    K-EFF = -------------------------------------- =  ---------
            !               NO. OF NEUTRONS IN GENERATION (N)      N_KEFFHIS
            !


            if (iCycle .le. pnSkip) then
               rKeff = real(N_KEFFSRC) / real(N_KEFFHIS)
            endif
     
            if (abs(rKeff - (real(N_KEFFSRC) / real(N_KEFFHIS))) <
     &          (OUTLIER_CONTROL)) then     
               rKeff = real(N_KEFFSRC) / real(N_KEFFHIS)
            else
               nTotCycles = nTotCycles + 1
            endif
            ! WE SET THE GUESS K-EFF OF NEXT CYCLE TO CURRENT K-EFF VALUE.
            KEFF_GUESS = rKeff

            ! CALCULATE AVERAGES
            if(iCycle .gt. pnSkip) then
                 
               rKeffAvg = (1.0 - 1.0 / real(iCycle-pnSkip)) * rKeffAvg +
     &                    rKeff / real(iCycle-pnSkip)

               rKeffSqAvg = (1.0 - 1.0 / real(iCycle-pnSkip)) 
     &                     * rKeffSqAvg +
     &                    rKeff**2 / real(iCycle-pnSkip)
            endif

            !-----------HERE WE END THE CPU TIME COUNTER
            call cpu_time(rCPUFinish)
            ! NOW CALCULATE THE CPU LEAD TIME FOR CURRENT FISSION CYCLE
            iCPUMinutes = int(floor((rCPUFinish-rCPUStart)/60.0))
            iCPUSeconds = int(rCPUFinish-rCPUStart)-60*iCPUMinutes
            iCPUMinTot = iCPUMinTot + iCPUMinutes
            iCPUSecTot = iCPUSecTot + iCPUSeconds

            
            ! IF CURRENT CYCLE IS LESS THAN SKIP COUNT,  WE IGNORE PRINTING THE 
            ! MOVING AVERAGE OF K-EFF.
            if(iCycle .le. pnSkip) then
               write(*,'(I10,2X,F10.5,2X,2A10,F10.3,2X,I6,A2,I2,A1)') 
     &             iCycle, rKeff, 'SKIP', 'INACTIVE',
     &             ShannonEntropy(), iCPUMinutes, 'm ', iCPUSeconds, 's'
               write(NFO,
     &         '(I10,2X,F10.5,2X,F10.5,2X,F10.5,I6,A2,I2,A1)') 
     &             iCycle, rKeff, 0.0,
     &             ShannonEntropy(), iCPUMinutes, 'm ', iCPUSeconds, 's'
            else
            write(*,'(I10,2X,F10.5,2X,F10.5,A10,F10.3,2X,I6,A2,I2,A1)')
     &             iCycle, rKeff, rKeffAvg , 'ACTIVE',
     &             ShannonEntropy(), iCPUMinutes, 'm ', iCPUSeconds, 's'
               write(NFO,
     &         '(I10,2X,F10.5,2X,F10.5,2X,F10.5,I6,A2,I2,A1)') 
     &            iCycle, rKeff, rKeffAvg, ShannonEntropy(), 
     &            iCPUMinutes, 'm ', iCPUSeconds, 's'
            endif
            ! NOW WE SET FISSION NEUTRONS AT THE COLLISION SITES AS THE NEUTRON
            ! SOURCE OF THE NEXT CRITICALITY CALCULATION CYCLE.
            KEFF_WGT = real(N_KEFFHIS) / real(N_KEFFSRC)
            
            !call SShannonEntropy()
            call SetNextGenerationSource()
            
            KEFF_CYCLE = KEFF_CYCLE + 1
            if(iCycle .gt. pnSkip) then
               call UpdateTallyAverages()
            endif

            
            ! END OF CYCLE LOOP
            enddo
            

            
            write(*,'(A,67A1)') '  ', ('-',i=1,67)
            print*
            ! WE PRINT THE AVERAGE VALUE OF K-EFF WITH ITS ERROR.
            rKeffError = sqrt((rKeffSqAvg - rKeffAvg**2) / 
     &               real(iCycle-pnSkip-1))
            write(*,'(A,F7.5,A,F7.5,A)') '  Estimated k-eff : (', 
     &         rKeffAvg, ' +/- ', rKeffError, ')'

            write(NFO,'(A,58A1)') '  ', ('-',i=1,58)
            
            rewind NFO
            
            i = 0
            do 
               i = i + 1
               read(NFO,'(A80)',end=298) caKeffLine(i)
            enddo
298         continue

            rewind NFO
            
            do j=1, i-1, 1
               write(NFO,'(A80)') caKeffLine(j)
               if(j .eq. 10) then
                  write(NFO,'(A,F7.5,A,F7.5,A)') 
     &               '  Estimated k-eff : (', 
     &               rKeffAvg, ' +/- ', rKeffError, ')'
                  if(iCPUSecTot .gt. 60) then
                     iCPUMinTot = iCPUMinTot + 
     &                 int(floor(real(iCPUSecTot)/60.0))
                     iCPUSecTot = iCPUSecTot - 
     &                 int(floor(real(iCPUSecTot)/60.0))
                     write(NFO,'(A,I0,A,I0,A)')     
     &                   '  Elapse time     : ',
     &                   iCPUMinTot, 'm ', iCPUSecTot, 's'     
                  else
                     write(NFO,'(A,I0,A,I0,A)')     
     &                   '  Elapse time     : ',
     &                   iCPUMinTot, 'm ', iCPUSecTot, 's'     
                  endif
                  write(NFO,'(A,G0.3,A)') 
     &                   '  Processing rate : ',
     &              real(pnHistories)/(iCPUMinTot*60.0 +
     &              iCPUSecTot*1.0), ' neutrons/s'
               endif
            enddo
            

               
 
            ! CLOSE NEUTRON TRACK FILE AND KEFF FILE.
            close(unit=NFT)
            close(unit=NFO)
            
            call WriteFluxToFile(pcRunID)
            call WritePowerToFile(NOMINAL_POWER,pcRunID)
            call WritePowerDistToFile(NOMINAL_POWER,pcRunID)
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
               call SetNeutronPosition(TRACK_NID, rX, rY, rZ, iStatus)
            else
               ! IF COLLISION IS POSSIBLE, WE TRANSPORT THE NEUTRON TOWARDS THE
               ! COLLISION SITE. r(k+1) = r(k) + D_col Ω
               rX = rX + rDistCol * rU
               rY = rY + rDistCol * rV
               rZ = rZ + rDistCol * rW
               call SetNeutronPosition(TRACK_NID, rX, rY, rZ, iStatus)
               
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

            if(KEFF_CYCLE .le. TOTAL_CYCLE) then
               call ScoreTrackLengthTally(TRACK_NID, 
     &               KEFF_WGT*real(N_KEFFHIS))
               return
            endif
            
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
            real    :: rINu, rNu, rNuRes
            integer :: iNu

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
            !rEpsilon = Rnd(int(TRACK_NID, 8))
            
            ! CALC. THE NUMBER  OF  FISSION NEUTRONS  ACCORDING TO THE INCOMING
            ! NEUTRON ENERGY.  iNu IS THE INTEGER VERSION OF THE rNu. WE OBTAIN
            ! THE AVG. NUMBER OF NEUTRONS (rNu)  FROM THE NUBAR-ENERGY CORRELA-
            ! TION. THIS CAN BE DONE BY USING THE GetNuBar FUNCTION.  SEE FUNC.
            ! DEFINITION IN NeutronInteractions MODULE.
            rNu  = GetNuBar()
         
            ! NEXT WE SAMPLE THE INTEGER VERSION OF AVG. NUMBER OF FISSION NEU-
            ! TRONS RELEASED PER FISSION.
c            rINu = floor(rNu)
c            if(rEpsilon .le. (rNu - rINu)) then
c               iNu = int(rINu) + 1
c            elseif(rEpsilon .gt. (rNu - rINu)) then
c               iNu = int(rINu)
c            endif
            
            ! THIS PART IS IMPORTANT TO PREVENT THE SHORTAGE OF NEW FISSION SRC
            ! SITES FOR THE USE OF NEXT CYCLE. THIS IS THE WEIGHTED SAMPLING OF
            ! THE NUMBER OF NEUTRONS TO BE RELEASE DURING  FISSION  REACTION AT 
            ! THE COLLISION SITE.
            !                
            !               _    Σf   1
            !         ν = W v  ----- ----
            !                   Σtot k(n) 
            !
            ! REFER ROMANO & FORGET (OPENMC)
            

               rNu = KEFF_WGT * (rYield / GetSigTot())
     &            * (1.0 / KEFF_GUESS) 
c               rEpsilon = Rnd(int(TRACK_NID,8))
c               if(rEpsilon .lt. (rNu - floor(rNu))) then
c                  rNu = floor(rNu)
c               elseif(rEpsilon .ge. (rNu - floor(rNu))) then
c                  rNu = floor(rNu) + 1.0
c               endif
c               iNu = int(rNu)
           
               iNu = ceiling(rNu)
      

            ! HERE WE UPDATE THE NUMBER OF DAUGHTER NEUTRONS PRODUCED
             if(iNu .le. 0) then
                NCHLD(TRACK_NID, 1) = 0
             else
                NCHLD(TRACK_NID, 1) = iNu
             endif
             
            ! WE EXIT THIS SUBROUTINE OF THERE IS NO FISSION NEUTRON IS PRODUCED.
            if(iNu .le. 0) return
            
            
            ! FOR EACH FISSION NEUTRONS PRODUCED, WE SAMPLE ITS OUTGOING ENERGY
            ! AND ITS OUTGOING DIRECTION. THEN WE STORE IT IN THE FISSION BANK.
            !$OMP PARALLEL DO
            do k=1, iNu, 1
            
            ! **PART I: SAMPLING THE OUTGOING ENERGY OF THE FISSION NEUTRON.
            
            ! FIRST WE TOSS FOR A RANDOM NUMBER.
            rEpsilon = Rnd(int(TRACK_NID, 8))
            
            ! NEXT WE PREPARE THE CUMMULATIVE DISTRIBUTION FUNC. OF THE OUTGOING
            ! NEUTRON ENERGY GROUP.
            rP(1) = 0.0
            !$OMP PARALLEL DO
            do j=1, TXS_GRP, 1
               rP(j+1) = TXS_TABLE(iLayer, iCell, j, 4) + rP(j)
            enddo
            !$OMP END PARALLEL DO
            ! THEN USING THE PREPARED RANDOM NUMBER,  rEpsilon,  WE USE THE INVERSE 
            ! METHOD TO SAMPLE THE OUTGOING ENERGY GROUP FROM THE CUMMULATIVE DIST.
            !$OMP PARALLEL DO
            do j=1, TXS_GRP, 1
               if((rEpsilon .ge. rP(j)) .and.
     &            (rEpsilon .lt. rP(j+1))) then
                  iFissionEnergyGroup = j
               endif
            enddo
            !$OMP END PARALLEL DO
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
     &                               iFissionEnergyGroup, KEFF_WGT)

            enddo
            !$OMP END PARALLEL DO
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
            !$OMP PARALLEL DO
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
               C_WGT(n, 1) = KEFF_WGT
               
              
            enddo
            !$OMP END PARALLEL DO
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
            !$OMP PARALLEL DO
            do iLayer=1, nMaxLayerID, 1
               do iCell=1, nMaxCellID, 1
                  call CellCentrePoint(iCell, iLayer, rX, rY, rZ)
                  do i=1, nNeutronPerCell, 1
                     call SampleIsotropicDirection(rU, rV, rW)
                     call AddPrimaryNeutron(rX, rY, rZ, 
     &                                      rU, rV, rW,
     &                                      1, KEFF_WGT)
                  enddo         
               enddo
            enddo          
            !$OMP END PARALLEL DO
            return
            
         end subroutine
         
         real function ShannonEntropy()
            integer :: i, j, k
            integer :: nMaxLayerID
            integer :: nMaxCellID
            real    :: rTotEntropy, rEntropy
            integer :: nTotSitesPerCell
            nMaxLayerID = LayerCount()
            rTotEntropy = 0.0
            rEntropy = 0.0
            if(RingCount() .eq. 6) nMaxCellID = 91
            if(RingCount() .eq. 7) nMaxCellID = 127
            !$OMP PARALLEL DO
            do i=1, nMaxLayerID, 1
               do j=1, nMaxCellID, 1
                  nTotSitesPerCell = 0
                  do k=1, N_KEFFSRC, 1
                     if((FC_CEL(k,1) .eq. j) .and.
     &                  (FC_CEL(k,2) .eq. i)) then
                        nTotSitesPerCell = nTotSitesPerCell + 1
                     endif
                  enddo
                  rEntropy = real(nTotSitesPerCell)/real(N_KEFFSRC)
                
                  if(nTotSitesPerCell .gt. 0) then
                     rTotEntropy = rTotEntropy - rEntropy * 
     &                              log(rEntropy) / log(2.0)                  
                  endif
                  
               enddo
            enddo
            !$OMP END PARALLEL DO
            ShannonEntropy = rTotEntropy
            return
         end function
         
         subroutine SShannonEntropy()
            integer :: i, j, k
            integer :: nMaxLayerID
            integer :: nMaxCellID
            real    :: rTotEntropy, rEntropy
            integer :: nTotSitesPerCell
            nMaxLayerID = LayerCount()
            rTotEntropy = 0.0
            rEntropy = 0.0
            if(RingCount() .eq. 6) nMaxCellID = 91
            if(RingCount() .eq. 7) nMaxCellID = 127
            
            do i=1, nMaxLayerID, 1
               do j=1, nMaxCellID, 1
                  nTotSitesPerCell = 0
                  do k=1, N_KEFFSRC, 1
                     if((FC_CEL(k,1) .eq. j) .and.
     &                  (FC_CEL(k,2) .eq. i)) then
                        nTotSitesPerCell = nTotSitesPerCell + 1
                     endif
                  enddo
                  
                  rEntropy = real(nTotSitesPerCell)/real(N_KEFFSRC)
                  print*, 'cell=',j,' layer=',i, ' ', rEntropy * 
     &                           log(rEntropy) / log(2.0)
                  if(nTotSitesPerCell .gt. 0) then
                     rTotEntropy = rTotEntropy - rEntropy * 
     &                              log(rEntropy) / log(2.0)                  
                  endif

                  
               enddo
            enddo

            return
         end subroutine        
         
      end module
