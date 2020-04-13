!                                                                       
!     Eigenvalue Calculation Module, Revision 200413-1.
!     Author M. R. Omar, October 2020. All copyrights reserved, 2020.
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
!     Last revision date 13/4/2020. Adapted survival biasing technique.
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
         real    :: KEFF_GUESS = 1.0
         real    :: NOMINAL_POWER = 1000.0
         ! KEFF_WGT IS THE PARTICLE WEIGHT THAT HELPS TO COMPENSATE THE PROBLEM
         ! IN WHICH THE FISSION NEUTRON  SITE COUNT  IS  NOT ENOUGH TO INITIATE
         ! NEXT NEUTRON GENERATION CYCLE.
         real    :: KEFF_WGT     = 1.0
         real    :: KEFF_WCUT    = 0.00001
         real    :: KEFF_WSUR    = 1.0
         
         ! DEFINING K-EFF ESTIMATORS.
         ! TRACKLENGTH KEFF ESTIMATOR
         real :: KEFF_TRLEST = 0.0
         ! ABSORPTION KEFF ESTIMATOR
         real :: KEFF_ABSEST = 0.0
         ! COLLISION KEFF ESTIMATOR
         real :: KEFF_COLEST = 0.0 
         
         ! KEFF_CYCLE INDICATES THE CURRENT CALCULATION CYCLE
         integer :: KEFF_CYCLE   = 1
         integer :: KEFF_SKIP    = 1
         
         ! TOTAL_CYCLE STORES THE TOTAL NUMBER OF CYCLE SET BY THE USER. THIS
         ! SHARED VARIABLE WILL BE SET OUTSIDE THIS MODULE.
         integer :: TOTAL_CYCLE  = 1
         ! INITIAL_HSRC IS THE INITIAL VALUE OF SOURCE SHANNON ENTROPY FOR CAL-
         ! CULATION.
         real    :: INITIAL_HSRC = 0.0
         real    :: SRC_H(50,127)
         
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
            integer :: nTotal
            real    :: rKeff                ! VALUE OF K-EFF FOR CURRENT CYCLE.
            real    :: rKeffAvg             ! A MOVING AVERAGE VALUE  OF K-EFF.
            real    :: rKeffSqAvg           ! A MOVING AVERAGE VALUE OF K-EFF^2
            real    :: rKeffError           ! STANDARD ERROR OF K-EFF.
            real    :: rHSrcAvg             ! A MOVING AVERAGE VALUE  OF H-SRC.
            real    :: rHSrcSqAvg           ! A MOVING AVERAGE VALUE OF H-SRC^2
            real    :: rHSrcError           ! STANDARD ERROR OF H-SRC.
            integer :: i, j, k, n, nTotCycles          ! MULTI-PURPOSE INTEGER VAR.
            integer :: NFO = 24
            integer :: NFT = 26
            character(len=40) :: cNumText
            real    :: rCPUStart, rCPUFinish
            integer :: iCPUMinutes, iCPUSeconds, iCPUMinTot, iCPUSecTot
            real    :: rX, rX2, rxBar, rx2Bar, rSBar, rRErr, rFOM
            real    :: rCurrentError
            real    :: rTotalSeconds = 0.0            ! TOTAL CPU TIME TAKEN WHEN RUNNING THE WHOLE BATCHES
            real    :: rTotalActiveSeconds = 0.0      ! TOTAL CPU TIME TAKEN WHEN RUNNING ACTIVE BATCHES
            character(len=80) :: caKeffLine(1000000)
            real :: rShannonH
            character(len=4) cDummy
            real    :: raHSrc(pnCycles)
            
            
            ! SETTING  THE  AVERAGE VALUES  TO ZERO BEFORE BEGIN KEFF AVERAGING. 
            ! CALCULATION.
            rKeffAvg   = 0.0
            rKeffSqAvg = 0.0
            rHSrcAvg   = 0.0
            rHSrcSqAvg = 0.0
            TOTAL_CYCLE = pnCycles
            iCPUMinTot = 0
            iCPUSecTot = 0
            KEFF_SKIP = pnSkip
            raHSrc(:) = 0.0d0
            
            ! OPEN FILE FOR WRITING THE K-EFF VALUE FOR ALL CYCLES.
            open(unit=NFO, file='KEFF.OUT', status='unknown')
            open(unit=NFT, file='TRACK.OUT', status='unknown')
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
            write(NFO,'(A,68A1)') '  ', ('-',i=1,68)
            write(NFO,'(A10,6(2X,A10))') 'CYCLE',
     &         'K-EFF', ' AVG K-EFF', 'STDEV', 'H-SRC', 'CPU TIME'
            write(NFO,'(A,68A1)') '  ', ('-',i=1,68)

            ! SETTING THE INITIAL GUESS OF K-EFF.
            KEFF_GUESS = prKeffGuess
            rKeff = KEFF_GUESS
            nTotCycles = pnCycles
            N_KEFFHIS = pnHistories

            ! HERE WE SEGREGATE THE INITIAL GUESS OF FISSION  SRC. DISTRIBUTION
            ! BY SETTING A MONO-ENERGETIC POINT SRC. AT THE CENTRE OF EACH UNIT
            ! CELL OF THE TRIGA CORE.
            call SegregateNeutronSource(pnHistories)
            
            ! PRINTING K-EFF TABLE HEADER IN CONSOLE OUTPUT
            print*
            write(*,'(A,75A1)') '  ', ('-',i=1,75)
            write(*,'(A10,2X,A10,2X,A10,2X,2A10,A11,A10)') 'CYCLE',
     &         'K-EFF', ' AVG K-EFF', 'STDEV', 'H-SRC', '   T/CYCLE',
     &         'CTM'
            write(*,'(A,75A1)') '  ', ('-',i=1,75)
            ! PRINT THE INITIAL GUESS OF K-EFF (CYCLE-0)
            call cpu_time(rCPUStart)
            write(*,
     &      '(I10,2X,F10.5,2X,A10,2X,A10,F10.3,2X,F8.1,A,F9.3)')
     &            0, KEFF_GUESS, 'GUESS', 'INACTIVE', -1.0
     &            * log(1.0/real(1*1))/log(2.0),
     &             0.0, 's ', rCPUStart/60.0
            write(NFO,'(I10,2X,F10.5,2X,A10,2X,A10,2X,F10.5,2X,A10)')
     &             0, KEFF_GUESS,
     &            '-', '-', -1.0
     &            * log(1.0/real(1*1))/log(2.0), '-'
     
            ! BEGIN THE CYCLE LOOP, THE NUMBER OF CYCLES IS pnCycles.
            do iCycle=1, nTotCycles, 1     
            
            ! RESET ALL K-EFF ESTIMATORS
            KEFF_TRLEST = 0.0
            KEFF_ABSEST = 0.0
            KEFF_COLEST = 0.0
            
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
            do while (iNID .lt. N_KEFFHIS)

               iNID = iNID + 1
               ! BEGIN TRACK NEUTRON ID = iNID.  THIS WILL SET THE NEUTRON BANK 
               ! ARRAYS CURSOR TO THE DESIRED NEUTRON.
               call BeginTracking(iNID)
               
               ! SHOW USERS CURRENT (NEUTRON-ID)-OF-(NUMBER OF HISTORIES).
c               if(mod(iNID,10) .eq. 0) then
c                  call cpu_time(rCPUFinish)
c                     write(*,'(A,I0,A,F0.2,A,I0,A,A1)',advance='no') 
c     &            '  > Progress: ',
c     &             100*iNID/N_KEFFHIS,
c     &            '%   Rate: ',
c     &            real(iNID)/(real(rCPUFinish-rCPUStart)/3600.0)
c     &            /1000000, 'M n/hr  NID: ', TRACK_NID, 
c     &            '                         ', char(13)                   
c
c               else
c                     write(*,'(A,I0,A,F0.2,A,I0,A1)',advance='no') 
c     &            '  > Progress: ', 100*iNID/N_KEFFHIS,
c     &            '%   Rate: ',
c     &            real(iNID)/(real(rCPUFinish-rCPUStart)/3600.0)
c     &            /1000000, 'M n/hr  NID: ', TRACK_NID, char(13)
c
c               endif
               
               ! BEGIN NEUTRON  HISTORY  SIMULATION UNTIL IT IS DEAD OR ESCAPED.
               ! IF THE NEUTRON IS STILL ALIVE AND INSIDE THE REACTOR, ITS ISTAT
               ! VALUE IS ZERO (0) (N_WAIT). CELL ID EQUALS TO -2 INDICATES THAT
               ! THE NEUTRON IS OUTSIDE THE REACTOR CORE.
               do while((ISTAT(iNID, 1) .eq. 0) .and.
     &                  (C_CEL(iNID, 1) .gt. -2) .and.
     &                  (C_CEL(iNID, 2) .gt. -2) .and.
     &                  (C_CEL(iNID, 1) .le. TXS_CEL) .and.
     &                  (C_CEL(iNID, 2) .le. TXS_LAY) .and.
     &                  (C_GRP(iNID, 1) .gt. 0) .and.
     &                  (C_GRP(iNID, 1) .le. TXS_GRP))

                  ! HERE WE USED A SPECIAL TRANSPORT CODE DESIGNED FOR KEFF CALCU
                  ! LATION. PROCEED TO Transport() SUBROUTINE FOR MORE INFO.
                  call Transport()
                  
                  ! HERE WE RECORD THE POSITION INTO A FILE FOR NEUTRON TRACK VI-
                  ! SUALIZATION.
                  if(iCycle .eq. pnCycles) then
                     if(iNID .le. 1000) then
                        write(NFT, '(I10,3F20.3)') iNID, L_POS(iNID,1),
     &                     L_POS(iNID,2), L_POS(iNID,3)
                     endif
                  endif
                  
               ! END TRANSPORT TRACKING OF A NEUTRON.
               enddo
               ISTAT(iNID,1) = N_LOST
            ! END HISTORY ITERATIONS FOR A GENERATION CYCLE.   
            enddo

            ! WE SET THE GUESS K-EFF OF NEXT CYCLE TO CURRENT K-EFF VALUE.

            rKeff = KEFF_COLEST / real(N_KEFFHIS)  ! real(N_KEFFSRC) / (real(N_KEFFHIS))
            KEFF_GUESS = rKeff
            
            ! CALCULATE AVERAGES AND STANDARD ERROR
            if(iCycle .gt. pnSkip) then                 
               rKeffAvg = (1.0 - 1.0 / real(iCycle-pnSkip)) * rKeffAvg +
     &                    rKeff / real(iCycle-pnSkip)

               rKeffSqAvg = (1.0 - 1.0 / real(iCycle-pnSkip)) 
     &                     * rKeffSqAvg +
     &                    rKeff**2 / real(iCycle-pnSkip)
               rKeffError = sqrt((rKeffSqAvg - rKeffAvg**2) / 
     &               real(iCycle-pnSkip-1))
               rKeffError = rKeffError / sqrt(real(iCycle-pnSkip-1))
            endif

            ! CALCULATE THE AVERAGE OF H-SRC FOR DECIDING CONVERGENCE

            
            ! HERE WE END THE CPU TIME COUNTER
            call cpu_time(rCPUFinish)
            rTotalSeconds = rTotalSeconds + rCPUFinish-rCPUStart
            rShannonH = ShannonEntropy()
            
            ! IF CURRENT CYCLE IS LESS THAN SKIP COUNT,  WE IGNORE PRINTING THE 
            ! MOVING AVERAGE OF K-EFF.
            if(iCycle .le. pnSkip) then
              
               write(*,
     &         '(I10,2X,F10.5,2X,A10,2X,A10,F10.6,2X,F8.1,A1,F10.3)')
     &             iCycle, rKeff, 'SKIP', 'INACTIVE',
     &             rShannonH,rCPUFinish-rCPUStart, 's',
     &            real(rCPUFinish)/60.0
               write(NFO,
     &         '(I10,2X,F10.5,2X,A10,2X,A10,2X,F10.5,2X,F9.1,A1)') 
     &             iCycle, rKeff, '-', '-',
     &             rShannonH, rCPUFinish-rCPUStart, 's'               

            else
               rTotalActiveSeconds = rTotalActiveSeconds + 
     &            rCPUFinish-rCPUStart
               write(*,
     &      '(I10,2X,F10.5,2X,F10.5,2X,F10.5,F10.6,2X,F8.1,A1,F10.3)')
     &             iCycle, rKeff, rKeffAvg , rKeffError,
     &             rShannonH, rCPUFinish-rCPUStart, 's',
     &            real(rCPUFinish)/60.0
c     &            int(1.0/((rKeffError/rKeffAvg)**2 * 
c     &            rTotalActiveSeconds/60.0))
               write(NFO,
     &         '(I10,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5,2X,F9.1,A1)') 
     &            iCycle, rKeff, rKeffAvg, rKeffError, rShannonH, 
     &            rCPUFinish-rCPUStart, 's'
            endif
            
            raHSrc(iCycle) = rShannonH
            ! NOW WE SET FISSION NEUTRONS AT THE COLLISION SITES AS THE NEUTRON
            ! SOURCE OF THE NEXT CRITICALITY CALCULATION CYCLE.
            KEFF_WGT = 1.0 ! real(N_KEFFHIS) / real(N_KEFFSRC)
          
            ! HERE WE SET THE NEXT GENERATION SOURCE.
            call SetNextGenerationSource()
            
            KEFF_CYCLE = KEFF_CYCLE + 1
            if(iCycle .gt. pnSkip) then
               call UpdateTallyAverages()
            endif

            ! END OF CYCLE LOOP
            enddo
            

            ! PRINT THE LINE
            write(*,'(A,75A1)') '  ', ('-',i=1,75)
            print*
            n = (pnCycles-pnSkip)/2
            ! COMPUTE THWE LAST HALF OF THE HSRC AVERAGE VALUE
            do i=1, n, 1
               rHSrcAvg = (1.0 - 1.0 / real(i)) * rHSrcAvg +
     &                    raHSrc(i+n+pnSkip-1) / real(i)

               rHSrcSqAvg = (1.0 - 1.0 / real(i)) 
     &                     * rHSrcSqAvg +
     &                    raHSrc(i+n+pnSkip-1)**2 / real(i)
               rHSrcError = sqrt((rHSrcSqAvg - rHSrcAvg**2) / 
     &               real(i-1))             
            enddo
            
            do i=1, pnCycles, 1
               if(abs(raHSrc(i)-rHSrcAvg) .lt. rHsrcError) then
                  write(*,'(A,I0,2A)') '  HGMC: Cycle ', i,
     &          ' is the first cycle that falls within one s.d.',
     &          ' of the source '
                  write(*,'(2A)') '        entropy population.',
     &          ' The number of skip cycles should be at least '
                  write(*,'(A)') '        more than this cycle.'
                goto 379
               endif
            enddo
379         continue         
            write(*,'(A,F0.5,A,F0.5,A)') 
     &      '  HGMC: H-SRC converged to an average value of ',
     &         rHSrcAvg, ' +/- ', rHSrcError / sqrt(real(n)), '.'

            ! WE PRINT THE AVERAGE VALUE OF K-EFF WITH ITS ERROR.
            rKeffError = sqrt((rKeffSqAvg - rKeffAvg**2) / 
     &               real(iCycle-pnSkip-1)) / 
     &               sqrt(real(iCycle-pnSkip-1))
            write(*,'(A,F0.5,A,F0.5,A)') '  HGMC: Estimated k-eff : ',
     &         rKeffAvg, ' +/- ', rKeffError, '.'
            write(NFO,'(A,68A1)') '  ', ('-',i=1,68)
            
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
                     write(NFO,'(A,G0.3,A)')     
     &               '  Elapse time     : ',   rTotalSeconds, 's'   
                  write(NFO,'(A,G0.3,A)') 
     &                   '  Processing rate : ',
     &              real(pnHistories)/(rTotalSeconds/pnCycles),
     &                ' neutrons/s'
               endif
            enddo
            
            call WriteFluxToFile(pcRunID)
            call WritePowerToFile(NOMINAL_POWER,pcRunID)
            call WritePowerDistToFile(NOMINAL_POWER,pcRunID)
               
 
            ! CLOSE NEUTRON TRACK FILE AND KEFF FILE.
            close(unit=NFT)
            close(unit=NFO)
            

            return
           
         end subroutine
         
!        ----------------------------------------------------------------------
!        THIS IS THE  TRANSPORT SUBROUTINE  MODIFIED FOR K-EFF CALCULATION. THE
!        DIFFERENCE  IS THAT WHEN NEUTRON UNDERGOES FISSION, THE CHILD NEUTRONS
!        ARE BANKED FOR LATER USE OF FISSION SOURCE OF  NEXT CYCLE.  WHEREAS IN 
!        THE STANDARD  FIXED  NEUTRON  SOURCE,  THE CHILD NEUTRONS HISTORY  ARE
!        IMMEDIATELY TRACKED AFTER BORN.
!        ----------------------------------------------------------------------
         subroutine Transport()
             
            ! PRIVATE VARIABLES
            real :: rDistCol        ! DISTANCE TO THE NEXT NEUTRON COLLISION.
            real :: rDistBoundary   ! NEAREST DISTANCE TO THE CELL BOUNDARY.
            real :: rX, rY, rZ      ! NEUTRON POSITION (x,y,z).
            real :: rU, rV, rW      ! NEUTRON DIRECTION Ω = (u,v,w).
            integer :: iCell, iLayer
            ! AT THE MOMENT iStatus IS JUST A DUMMY VARIABLE STORING THE STATUS
            ! OF NEUTRON REACTION OPERATIONS.
            integer :: iStatus
            real :: rYield
            real :: rDist

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
            
            rDist = min(rDistCol, rDistBoundary)
            
            if((rDistCol .gt. rDistBoundary) .and. 
     &         (rDistBoundary .gt. 0.00001)) then 
                  

               ! THE DISTANCE TO THE NEAREST BOUNDARY IS LESS THAN THE DISTANCE 
               ! TO THE  NEXT COLLISION,  WE TRANSPORT  THE NEUTRON TO THE CELL 
               ! BOUNDARY. r(k+1) = r(k) + D Ω
               rX = rX + rDistBoundary * rU
               rY = rY + rDistBoundary * rV
               rZ = rZ + rDistBoundary * rW
               
                  call GetCurrentCell(rX, rY, rZ,
     &                               iCell, iLayer)
                  L_POS(TRACK_NID, 1) = C_POS(TRACK_NID, 1)
                  L_POS(TRACK_NID, 2) = C_POS(TRACK_NID, 2)
                  L_POS(TRACK_NID, 3) = C_POS(TRACK_NID, 3)      
                  C_POS(TRACK_NID, 1) = rX
                  C_POS(TRACK_NID, 2) = rY
                  C_POS(TRACK_NID, 3) = rZ
                  L_CEL(TRACK_NID, 1) = C_CEL(TRACK_NID, 1)
                  L_CEL(TRACK_NID, 2) = C_CEL(TRACK_NID, 2)
                  C_CEL(TRACK_NID, 1) = iCell
                  C_CEL(TRACK_NID, 2) = iLayer
            else
               ! IF COLLISION IS POSSIBLE, WE TRANSPORT THE NEUTRON TOWARDS THE
               ! COLLISION SITE. r(k+1) = r(k) + D_col Ω
               rX = rX + rDistCol * rU
               rY = rY + rDistCol * rV
               rZ = rZ + rDistCol * rW
               call GetCurrentCell(rX, rY, rZ,
     &                               iCell, iLayer)
  
                  L_POS(TRACK_NID, 1) = C_POS(TRACK_NID, 1)
                  L_POS(TRACK_NID, 2) = C_POS(TRACK_NID, 2)
                  L_POS(TRACK_NID, 3) = C_POS(TRACK_NID, 3)      
                  C_POS(TRACK_NID, 1) = rX
                  C_POS(TRACK_NID, 2) = rY
                  C_POS(TRACK_NID, 3) = rZ
                  L_CEL(TRACK_NID, 1) = C_CEL(TRACK_NID, 1)
                  L_CEL(TRACK_NID, 2) = C_CEL(TRACK_NID, 2)
                  C_CEL(TRACK_NID, 1) = iCell
                  C_CEL(TRACK_NID, 2) = iLayer   
               
               ! NEXT WE SCORE COLLISION K-EFF ESTIMATOR
               if (GetSigTot() .gt. 0.0) then
                  KEFF_COLEST = KEFF_COLEST + C_WGT(TRACK_NID,1) *
     &               GetNuBar() * GetSigFis() / GetSigTot()                  
               endif


               ! IF COLLISION IS POSSIBLE WE FIRST GENERATE THE FISSION SITES,
               ! IF THE TOTAL FISSION CROSS SECTION IS > 0.
               if(GetSigFis() .gt. 0.0) then
                  call GenerateFissionSites()
               endif
               
               ! HERE IN SURVIVAL BIASING, ABSORPTION REACTION IS PROHIBITED, BUT
               ! WE REDUCE THE PARTICLE WEIGHT INSTEAD. THIS ROUTINE ADJUSTS THE
               ! WEIGHT.
               ! DETERMINE WEIGHT ABSORBED IN SURVIVAL BIASING.
               if(GetSigAbs() .gt. 0.0) then
                  call Absorption()
               endif
               
               ! NEXT WE PROCEED FOR RUSSIAN ROULETTE TO DETERMINE WHETHER THE 
               ! NEUTRON SHOULD BE FURTHER TRACKED OR KO BE KILLED.
               call RussianRoulette()
               
               ! FINALLY, WE PROCESS SCATTERING REACTION.
               call ScatterNeutron(.true.)
               
            endif
            ! END OF IF-ELSE BLOCK CHECKING WHETHER COLLISION HAPPENS.
            
            ! HERE WE SCORE K-EFF TRACKLENGTH TALLY 
            KEFF_TRLEST = KEFF_TRLEST + C_WGT(TRACK_NID,1) *
     &         rDist * GetNuBar() * GetSigFis()
     
            if(KEFF_CYCLE .le. TOTAL_CYCLE) then
               ! HERE WE SCORE THE FLUX TRACKLENGTH TALLY. WE SCORE EACH TIME THE 
               ! NEUTRON MAKE A DEPARTURE.            
               call ScoreTrackLengthTally(TRACK_NID, 
     &               KEFF_WGT*real(N_KEFFHIS))
            endif
            return
         end subroutine
         
         
         subroutine Absorption()
            if(GetSigAbs() .le. 0.0) return
            if(GetSigTot() .le. 0.0) return
            ! DETERMINE WEIGHT ABSORBED IN SURVIVAL BIASING.
            C_AWT(TRACK_NID,1) = C_WGT(TRACK_NID,1) *
     &         GetSigAbs() / GetSigTot()
     
            ! ADJUST CURRENT PARTICLE WEIGHT
            C_WGT(TRACK_NID,1) = C_WGT(TRACK_NID,1) - 
     &         C_AWT(TRACK_NID,1)
            L_WGT(TRACK_NID,1) = C_WGT(TRACK_NID,1)

            ! HERE WE SCORE ABSORPTION KEFF ESTIMATOR.
            KEFF_ABSEST = KEFF_ABSEST + C_AWT(TRACK_NID,1) *
     &         GetNuBar() *  GetSigFis() / GetSigAbs()
         end subroutine
         
         subroutine RussianRoulette()
            real :: rEpsilon            
            rEpsilon = Rnd(int(TRACK_NID,8))
            if(C_WGT(TRACK_NID,1) .lt. KEFF_WCUT) then
               if(rEpsilon .lt. 
     &               (C_WGT(TRACK_NID,1)/KEFF_WSUR)) then
                  C_WGT(TRACK_NID,1) = KEFF_WSUR
                  L_WGT(TRACK_NID,1) = C_WGT(TRACK_NID,1)
               else
                  C_WGT(TRACK_NID,1) = 0.0
                  L_WGT(TRACK_NID,1) = 0.0
                  ISTAT(TRACK_NID,1) = N_LOST
                  
               endif
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
         subroutine GenerateFissionSites()
            
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

            real    :: rYield, rCosPhi, rSinPhi
            real    :: rEpsilon1, rEpsilon2
            ! HERE WE OBTAIN ALL OF THE CURRENT NEUTRON IDENTITY.
            iGroup = C_GRP(TRACK_NID, 1)
            iCell  = C_CEL(TRACK_NID, 1)
            iLayer = C_CEL(TRACK_NID, 2)
            rU     = C_DIR(TRACK_NID, 1)
            rV     = C_DIR(TRACK_NID, 2)
            rW     = C_DIR(TRACK_NID, 3)
            rYield = TXS_TABLE(iLayer, iCell, iGroup, 3)

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
           

            rNu = C_WGT(TRACK_NID,1) / KEFF_GUESS * rYield 
     &            / GetSigTot()
            rEpsilon = Rnd(int(TRACK_NID,8))
            if(rEpsilon .lt. (rNu - floor(rNu))) then
               rNu = floor(rNu)
            elseif(rEpsilon .ge. (rNu - floor(rNu))) then
              rNu = floor(rNu) + 1.0
            endif
            iNu = int(rNu)

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
            ! THEN USING THE PREPARED RANDOM NUMBER,  rEpsilon,  WE USE THE INVERSE 
            ! METHOD TO SAMPLE THE OUTGOING ENERGY GROUP FROM THE CUMMULATIVE DIST.

            rP(1) = 0.0

            do j=1, TXS_GRP, 1
               rP(j+1) = TXS_TABLE(iLayer, iCell, j, 4) + rP(j)
               if((rEpsilon .ge. rP(j)) .and.
     &            (rEpsilon .lt. rP(j+1))) then
                 iFissionEnergyGroup = j
                  goto 569
               endif
            enddo
569         continue

            ! **PART II: SAMPLING THE OUTGOING DIRECTION OF THE FISSION NEUTRON.
            ! FIRST WE SAMPLE THE SCATTERING COSINE, μ ∈ (-1, 1).
            rMu = 2.0 * Rnd(int(TRACK_NID,8))  - 1.0
            
            ! NEXT WE SAMPLE A RANDOM NUMBER, ϕ ∈ [0,2π).
            rRandomPhi = 2.0 * PI * Rnd(int(TRACK_NID,8))      
            
            rCosPhi = cos(rRandomPhi)
            rSinPhi = sin(rRandomPhi)
            ! NEXT WE CALCULATE THE OUTGOING NEUTRON DIRECTION USING (μ, ϕ).
            rUNew = rMu
            rVNew = sqrt(1.0-rMu*rMu) * rCosPhi   
            rWNew = sqrt(1.0-rMu*rMu) * rSinPhi 
     
            ! THEN WE ADD THE FISSION NEUTRON INTO  THE FISSION BANK SO THAT WE
            ! CAN USE IT LATER AS A SOURCE FOR NEXT CYCLE.
            call AddFissionNeutron(C_POS(TRACK_NID, 1),
     &                               C_POS(TRACK_NID, 2),
     &                               C_POS(TRACK_NID, 3),
     &                               rUNew, rVNew, rWNew,
     &                               iFissionEnergyGroup, KEFF_WGT,
     &                               iCell, iLayer)

            enddo
            !$OMP END PARALLEL DO
            ! END LOOP FOR EACH NEWBORN FISSION NEUTRON

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
            integer :: nTotal
            ! HERE WE SELECT A NEUTRON FOR N_KEFFHIS TIMES,SO THAT THE NEXT GE-
            ! NERATION HAS THE SAME AMOUNT OF NEUTRON SOURCE. THE SELECTED NEU-
            ! TRON IS THEN TRANSFERED TO THE MAIN NEUTRON BANK.
            !$OMP PARALLEL DO
            
            do n=1, N_KEFFHIS, 1
            
               ! HERE WE RANDOMLY SHUFFLE THE NEUTRON SELECTIONS.
1              i = int(ceiling(N_KEFFSRC * Rnd(int(BANK_SIZE-10,8))))

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
            
            real :: rEpsilon1, rEpsilon2
            real :: rMu, rPhi
            real, parameter :: PI = 3.141592654
            
            
594         rEpsilon1 = Rnd(int(TRACK_NID,8))
            rEpsilon2 = Rnd(int(TRACK_NID,8))
            
            rEpsilon1 = 2.0*rEpsilon1-1
            rEpsilon2 = 2.0*rEpsilon2-1
            if((rEpsilon1**2 + rEpsilon2**2) .gt. 1) then
               goto 594
            endif

            prU = 2.0*rEpsilon1**2 + 2.0*rEpsilon2**2 - 1.0
            prV = rEpsilon1 * sqrt(1.0 - prU**2)/
     &            sqrt(rEpsilon1**2 + rEpsilon2**2)
            prW = rEpsilon2 * sqrt(1.0 - prU**2)/
     &            sqrt(rEpsilon1**2 + rEpsilon2**2)     
   
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
     &                  rU, rV, rW, 1, KEFF_WGT, 1)    

                  enddo         
               enddo
            enddo   
            
            return
            
         end subroutine
         
         real function ShannonEntropy()
            integer :: i, j, k
            integer :: nMaxLayerID
            integer :: nMaxCellID
            real    :: rTotEntropy, rEntropy
            integer :: nTotSitesPerCell
            SRC_H(:,:) = 0.0d0
            nMaxLayerID = LayerCount()
            
            do k=1, N_KEFFSRC, 1
               

               SRC_H(FC_CEL(k,2), FC_CEL(k,1)) = 
     &            SRC_H(FC_CEL(k,2), FC_CEL(k,1)) + 1.0                  


            enddo

            if(RingCount() .eq. 6) nMaxCellID = 91
            if(RingCount() .eq. 7) nMaxCellID = 127
            rTotEntropy = 0.0
            do i=1, nMaxLayerID, 1
               do j=1, nMaxCellID, 1

                   if((SRC_H(i, j) .gt.  0) .and. 
     &            (SRC_H(i,j) .lt. real(N_KEFFSRC))) then
                        rEntropy = SRC_H(i,j)/
     &                    real(N_KEFFSRC)     
                rTotEntropy = rTotEntropy - rEntropy * 
     &                              log(rEntropy) / log(2.0)       
                   endif

               enddo
            enddo

            ShannonEntropy = rTotEntropy

            return 
            
         end function
  
 
      end module
