!                                                                       
!     Main TRIMON Program, Revision 180915-1.
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

      program Main
         use NeutronBankManager
         use TXSReader
         use RandomNumberGenerator
         use NeutronTrack
         use CriticalityCalculation
         use Tally
         use CellGeometry
         
!        Begin your code here...
         integer :: iStatus, i
         integer :: iTXSUnit = 20
         integer :: iSeed, nHistories, nCycle, nSkip, nRing  
         integer :: NFI = 121 
         integer :: tCol, tC, tL, tG 
         real    :: rStart, rFinish
         real    :: keffGuess
         real    :: tX, tY, tZ
         character(len=80) :: cTextTrv
         character(len=80) :: cTXSFileName
         logical :: lFlag
         character(len=30) :: cRunID
         character(len=30) :: tDate, tTime
                  call cpu_time(rStart)    

         ! HERE WE GENERATE RUN ID FOR SIMULATION IDENTIFICATION PURPOSE
         call date_and_time(DATE=tDate,TIME=tTime)       
         cRunID = trim(tDate)//'-'// tTime(1:6)// '-'// tTime(8:10) 
         print*
c         print*,' TRIMON-HGMC Version 1.190430.2335'
         print*,' TRIMON-PACKAGE Version 1.190430.2335'
         print*,' (c) 2018 Malaysian Nuclear Agency'
         print*,' (c) 2018 Universiti Sains Malaysia'
         print* 
         if(iargc() .eq. 1) then
            call getarg(1, cTXSFileName)
            if(cTXSFileName .eq. '-update-burnup') then
               call BurnupCalculation(cRunID)
               print*, ' BURN: Job ended.'
               stop
            elseif(cTXSFileName .eq. '-validate-txs') then
               cTXSFileName = 'xsdata.txs'
               open(unit=iTXSUnit, file=cTXSFileName, err=95)
               call TXSVALIDATE(iTXSUnit, lFlag, .false.)
               close(unit=iTXSUnit)
               goto 98
            
            elseif(cTXSFileName(1:1) .eq. '-') then
               print*, ' HGMC: Fatal error. Invalid parameter.'
               stop
            else
               print*, ' HGMC: Loading cross section data from ',
     &             trim(cTXSFileName)
               open(unit=915, file=cTXSFileName, status='old',
     &              err=989)
               close(unit=915)
               goto 43
989            print*, ' HGMC: Fatal error. The provided TXS',
     &            ' tape does not exists.'
               stop
43             continue
            endif
         else
            cTXSFileName = 'xsdata.txs'
            call system("libreader")
         endif
         


         open(unit=NFI, file='MAIN.INP', status='unknown',err=96)
         call rasearch(NFI,0,'#KRUN', 5, cTextTrv, iStatus)
         if (iStatus .gt. 0) goto 99
         read(NFI,*,end=97,err=97) iSeed, nHistories, nCycle, nSkip
         read(NFI,*,end=97,err=97) keffGuess
         print'(A,I0,2A)', '  HGMC: ', nHistories, 
     &      ' neutron histories will be simulated.'
         if(nCycle .lt. nSkip) then
            print*, ' HGMC: Warning. ',
     &         'The number of cycles is less than the',
     &         'number of skip cycles.'
            print*, ' HGMC: Total number of cycles is',
     &            ' adjusted to 2*NSKIP.'
            nCycle = 2*nSkip
         endif
         
         if(keffGuess .lt. 0.0) then
            print*, ' HGMC: Warning. The initial guess of',
     &            ' k-eff is less than zero.'
            print*, ' HGMC: The magnitude of k-eff is used instead.'
            keffGuess = abs(keffGuess)
         endif
         if(iSeed .lt. 0) then
            print*, ' HGMC: Warning. Random number generator seed is',
     &         ' less than zero.'
            print*, ' HGMC: The seed is set to 92090914.'
            iSeed = 92090914
         endif
         if(nHistories .lt. 3000) then
            print*, ' HGMC: Warning. ',
     &         'Insufficient neutron histories (N<3000)'
            print*, ' HGMC: Adjusting to 5000 histories.'
            nHistories = 5000
         elseif(nHistories .gt. BANK_SIZE) then
            print*, ' HGMC: Warning. The prescribed number of',
     &              ' histories is too large.'
            print*, ' HGMC: Adjusting to 50000 histories.'
            nHistories = 30000
         endif
         
         call rasearch(NFI,0,'#NRINGS', 7, cTextTrv, iStatus)
         if (iStatus .gt. 0) goto 92
         read(NFI,*,end=91,err=91) RING_COUNT

         call rasearch(NFI,0,'#NLAYERS', 8, cTextTrv, iStatus)
         if (iStatus .gt. 0) goto 92
         read(NFI,*,end=91,err=91) LAYER_COUNT
         
         call rasearch(NFI,0,'#DIMENSIONS', 11, cTextTrv, iStatus)
         if (iStatus .gt. 0) goto 94
         read(NFI,*,end=93,err=93) (RING_RADIUS(i),i=1,RING_COUNT) 
         read(NFI,*,end=93,err=93) FUEL_LENGTH
         read(NFI,*,end=93,err=93) REFLECTOR_THICKNESS

         call rasearch(NFI,0,'#POWER', 6, cTextTrv, iStatus)
         if (iStatus .gt. 0) goto 92
         read(NFI,*,end=91,err=91) NOMINAL_POWER
         

         
         open(unit=iTXSUnit, file=cTXSFileName, err=95)
!        HERE WE INITIALIZE CELL INDICES FOR MORE CODE PERFORMANCE
!        THIS IS TO PREVENT EXTRA LOOPS WHEN SEARCHING FOR CELL INDEX.
         call InitCellIndex()
         call InitTXSReader(iTXSUnit)
         print*, ' HGMC: Initializing.'
         if(LAYER_COUNT .ne. TXS_LAY) then
            print*, ' HGMC: Warning. Library layer count mismatch.'
            print'(A,A,I0,A)',
     &       '        Rectifying the problem by setting layer',
     &         ' count to ', TXS_LAY, '.'
            LAYER_COUNT = TXS_LAY
         endif
         if((RING_COUNT .eq. 6) .and. (TXS_CEL .eq. 92)) then
            goto 126
         else
            if((RING_COUNT .eq. 7) .and. (TXS_CEL .eq. 128)) then
               goto 126
            else
               print*, ' HGMC: Fatal Error. The incompatible library.'
               print*, '       Please consider changing the ring count.'
               stop
            endif         
         endif
126      continue         
         call InitTally()
         call TXSGetTable(iTXSUnit, istatus)
         
         call InitNeutronBank()
         call InitRnd(10000000_8,500000_8,int(iSeed,8))

         call CalculateKeff(nHistories,nCycle, nSkip, keffGuess,
     &      cRunID)
         close(unit=NFI)
         call BurnupCalculation(cRunID)
         call cpu_time(rFinish)
         print'("  HGMC: Total CPU time ",F0.3," min(s).")',
     &       (rFinish - rStart) / 60.0
         goto 98
         
! HERE ARE THE ERROR MESSSAGES...
91       print*, ' HGMC: Fatal Error. NRINGS card is not valid.'
         goto 100
92       print*, ' HGMC: Fatal Error. NRINGS card is not found.'
         goto 100
93       print*, ' HGMC: Fatal Error. DIMENSIONS card is not valid.'
         goto 100
94       print*, ' HGMC: Fatal Error. DIMENSIONS card is not found.'
         goto 100
95       print*, ' HGMC: Fatal Error. TRIGA Cross Section file',
     &           ' (.txs) is not found.'
         goto 101
96       print*, ' HGMC: Fatal Error. The input file is not found..'
         goto 101
97       print*, ' HGMC: Fatal Error. KRUN card input is not valid.'
         goto 100
99       print*, ' HGMC: Fatal Error. KRUN card is not found',
     &           ' in the input file.'
100      close(unit=NFI)
98       print*, ' HGMC: Job ended.'
       !  call system("pause")

101      continue

      end program

      subroutine rasearch(pNF, pnSkip, pcText, pnTextLen, 
     &                    pcTextRV, piFlag)
         
         character(len=*), intent(in)  :: pcText
         character(len=*), intent(out) :: pcTextRV ! Text return value (RV)
         integer,          intent(in)  :: pNF, pnSkip, pnTextLen
         integer,          intent(out) :: piFlag
         integer                         :: i
         
         piFlag = 0
         i = pnSkip
         rewind pNF
1        continue
         read(pNF,'(A,A)', end=999) pcTextRV
         if(pcText(1:pnTextLen) .ne. pcTextRV(i+1:i+pnTextLen)) goto 1
         return
999      piFlag = 1
         return
      end subroutine
      
  
