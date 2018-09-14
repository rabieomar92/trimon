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
         real :: rStart, rFinish
         integer :: iSeed, nHistories, nCycle, nSkip, nRing
         real    :: keffGuess
         integer :: NFI = 12
         character(len=80) :: cTextTrv
         character(len=80) :: cTXSFileName
         logical :: lFlag
         integer :: tCol, tC, tL, tG
         real    :: tX, tY, tZ
         if(iargc() .eq. 1) then
            call getarg(1, cTXSFileName)
            if(cTXSFileName .eq. '-update-burnup') then
               call BurnupCalculation()
               goto 98
            elseif(cTXSFileName .eq. '-validate-txs') then
               cTXSFileName = 'xsdata.txs'
               open(unit=iTXSUnit, file=cTXSFileName, err=95)
               call TXSVALIDATE(iTXSUnit, lFlag, .false.)
               close(unit=iTXSUnit)
               goto 98
            endif
       
         else
            cTXSFileName = 'xsdata.txs'
         endif
         print*, ' HGMC: Loading cross section data from ',
     &           trim(cTXSFileName)
         open(unit=NFI, file='MAIN.INP', status='old',err=96)
         call rasearch(NFI,0,'#KRUN', 5, cTextTrv, iStatus)
         if (iStatus .gt. 0) goto 99
         read(NFI,*,end=97,err=97) iSeed, nHistories, nCycle, nSkip,
     &       OUTLIER_CONTROL
         read(NFI,*,end=97,err=97) keffGuess
         
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
c         read(NFI,*,err=93) TOPREFLECTOR_THICKNESS
c         read(NFI,*,err=93) BOTREFLECTOR_THICKNESS

         call rasearch(NFI,0,'#POWER', 6, cTextTrv, iStatus)
         if (iStatus .gt. 0) goto 92
         read(NFI,*,end=91,err=91) NOMINAL_POWER
         
         close(unit=NFI)
         
         open(unit=iTXSUnit, file=cTXSFileName, err=95)
         call InitTXSReader(iTXSUnit)
         print'(A,$)', '  HGMC: Initializing tally counter.'
         call InitTally()
         print*, 'Done.'
         call TXSGetTable(iTXSUnit, istatus)
         
       
         
         print'(A,I0)', '  Number of cells: ', TXS_CEL*TXS_LAY
         print'(A,$)', '  HGMC: Initializing fission bank.'
         call InitNeutronBank()
         print*, 'Done.'
         print'(A,$)', '  HGMC: Initializing random numbers.'
         call InitRnd(10000000_8,500000_8,int(iSeed,8))
         print*, 'Done.'
         call cpu_time(rStart)    
         call CalculateKeff(nHistories,nCycle, nSkip, keffGuess)
         call BurnupCalculation()
         call cpu_time(rFinish)
         print'("  HGMC: Total CPU time ",G0.3," sec(s).")',
     &       (rFinish - rStart)
         goto 98
         
! HERE ARE THE ERROR MESSSAGES...
91       print*, ' Fatal Error: NRINGS card is not valid.'
         goto 100
92       print*, ' Fatal Error: NRINGS card is not found.'
         goto 100
93       print*, ' Fatal Error: DIMENSIONS card is not valid.'
         goto 100
94       print*, ' Fatal Error: DIMENSIONS card is not found.'
         goto 100
95       print*, ' Fatal Error: TRIGA Cross Section file',
     &           ' (.txs) is not found.'
         goto 101
96       print*, ' Fatal Error: The input file is not found.'
         goto 101
97       print*, ' Fatal Error: KRUN card input is not valid.'
         goto 100
99       print*, ' Fatal Error: KRUN card is not found',
     &           ' in the input file.'
100      close(unit=NFI)
98       print*, ' HGMC: Code execution terminated successfully!'
101   end program

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
