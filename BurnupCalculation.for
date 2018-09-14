!                                                                       
!     Fuel Burnup Calculation Module, Revision 180915-1.
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
!     This source code was created on 3/3/2018 11:26 AM by M. R. Omar.
!     Last revision date 15/9/2018.
!
      
      subroutine BurnupCalculation()
      
         integer, parameter :: NFO = 68
         integer, parameter :: NFI = 69
         integer, parameter :: NFEL = 70
         integer, parameter :: NFPPE = 71
         integer, parameter :: NFPPB = 72
         integer, parameter :: NFPPC = 73
         
  
         integer :: iaNBeta(3), i, j, k, iRing, iFlag, iCellID, n1
         integer :: iEl
         real    :: rBurnupIncrement, rPowerPerElement
         real    :: rCellPower, rCellPowerErr
         real    :: raBetaST8(30)
         real    :: raBetaST12(30)
         real    :: raBetaST20(30)
         real    :: raInput(30)
         integer :: nBeta
         character(len=4)  :: cFuelType
         character(len=80) :: cTextTrv
         character(len=4)  :: cIden(6), cElem(6)
         character(len=4)  :: cLoading(128), cSite(128)
         character(len=4)  :: cCellType
         real    :: rMassU, rEnrichment, rMassEr166, rMassEr167
         real    :: rOldBurnupMWD, rOldBurnupPercent
         real    :: rBurnupMWD, rBurnupPercent, rTotalCoreBurnup
         real    :: raCellPower(6400), raCellPowerErr(6400)
         character(len=8) :: cCellID
         open(unit=NFI, file="MAIN.INP", status='OLD', err=98)
         open(unit=NFEL, file="FUEL_INVENTORY.INP", 
     &        status='OLD', err=99)
         open(unit=NFPPE, file="PPE.OUT", status='OLD', err=102)
         open(unit=NFPPB, file="PDIST.OUT", status='OLD', err=104)
         open(unit=NFPPC, file="ZBURN.OUT", status='UNKNOWN', err=105)
         open(unit=NFO, file="FUEL_INVENTORY.OUT",
     &        status='UNKNOWN', err=103)
     


!     FIRST WE READ BETA PARAMETERS.
         call rasearch(NFI,0,'#BURNUP', 7, cTextTrv, iStatus)
         if (iStatus .gt. 0) goto 94
         read(NFI,*,err=97) rBurnupIncrement
         
         do j=1, 3, 1
         
         read(NFI,*,err=97) cFuelType, nBeta, (raInput(i) ,i=1,nBeta)
         if(cFuelType .eq. 'FE08') then
            iaNBeta(1) = nBeta
            do i=1, 30, 1
               raBetaST8(i) = 0.0
            enddo
            do i=1, nBeta, 1
               raBetaST8(i) = raInput(i)
            enddo
         elseif(cFuelType .eq. 'FE12') then
            iaNBeta(2) = nBeta
            do i=1, 30, 1
               raBetaST12(i) = 0.0
            enddo
            do i=1, nBeta, 1
               raBetaST12(i) = raInput(i)
            enddo         
         elseif(cFuelType .eq. 'FE20') then
            iaNBeta(3) = nBeta
            do i=1, 30, 1
               raBetaST20(i) = 0.0
            enddo
            do i=1, nBeta, 1
               raBetaST20(i) = raInput(i)
            enddo         
         endif
         
         enddo
         
         rewind NFI
!     NEXT WE READ NRINGS
         call rasearch(NFI,0,'#NRINGS',7, cTextTrv, iFlag)
         if (iFlag .gt. 0) then
            print*, ' Fatal Error: NRINGS card is not found!'
            goto 96
         endif
         read (NFI,*,err=100) iRing
         rewind NFI
!     NEXT WE READ THE CORE LOADING
         call rasearch(NFI,0,'#CORECONFIG',11, cTextTrv, iFlag)
         if (iFlag .gt. 0) then
             print*, ' BURNUP: Fatal error - LOADING card is not found!'
             goto 96
         endif
!     Read and write the loading of the first ring, A-1.
         read (NFI,'(A4,1x,A4)') cIden(1), cElem(1)
!     Write to the output parameter. This is for TXS file cell type identifying purpose.
         cLoading(1) = cElem(1)
         cSite(1)    = cIden(1)
         iCellID = 2
!     Read and write fuel elements loading for row 2 - 6
         do i=1, 15, 1
            read (NFI,'(A4,1x,A4,5(1x,A4,1x,A4))', end=100, err=100) 
     &           (cIden(j), cElem(j), j=1, 6)
            do k=1, 6, 1
               cLoading(iCellID) = cElem(k)
               cSite(iCellID)    = cIden(k)
               iCellID = iCellID + 1
            enddo
         enddo
!     Read and write row 7 (G ring) for reactors with 7 fuel rings.
         if (iRing .eq. 7) then
            do j=1, 6, 1
               read (NFI,'(A4,1x,A4,5(1x,A4,1x,A4))', end=100, err=100) 
     &              (cIden(i), cElem(i),i=1, 6)
               do k=1, 6, 1
                  cLoading(iCellID) = cElem(k)
                  cSite(iCellID)    = cIden(k)
                  iCellID = iCellID + 1
               enddo
            enddo
         endif

!     Read  and write reflector type
         read  (NFI,'(A4,1x,A4)', end=100, err=100) cIden(1), cElem(1)
         cLoading(iCellID) = 'G   '         
         cSite(iCellID)    = 'R   '
        
! ---------------------------------------------------------------------------------
      
        print*, ' BURNUP: Updating core burnup in FUEL_INVENTORY.OUT.' 
        write(NFO,'(A4,1X,A4,1X,7A10)') 'id','type','mU[g]', 'e[%]',
     &      'Er166[g]', 'Er167[g]','BU[MWd]','BU[%]','BUi[MWd]'
     
        rTotalCoreBurnup = 0.0
        do iEl=1, iCellID, 1
        
            
            call ReadEle(NFEL, cLoading(iEl), cCellType, rMassU, 
     &                   rEnrichment, rMassEr166, rMassEr167,
     &                   rOldBurnupMWD, rOldBurnupPercent)

            if((cCellType .eq. 'FE08') .or.
     &         (cCellType .eq. 'FE12') .or.
     &         (cCellType .eq. 'FE20')      ) then
  
            call ReadPPE(NFPPE, cSite(iEl), rPowerPerElement)
            rBurnupMWD = rPowerPerElement * rBurnupIncrement /
     &                   10.0**3
            rTotalCoreBurnup = rTotalCoreBurnup + rBurnupMWD
            
            rBurnupPercent = 0.0
            if(cCellType .eq. 'FE08') then
               do n1=1, iaNBeta(1), 1
                  rBurnupPercent = rBurnupPercent + 
     &               raBetaST8(n1)*(rBurnupMWD+rOldBurnupMWD)**n1
               enddo
            elseif(cCellType .eq. 'FE12') then
               do n1=1, iaNBeta(2), 1
                  rBurnupPercent = rBurnupPercent + 
     &               raBetaST12(n1)*(rBurnupMWD+rOldBurnupMWD)**n1               
               enddo            
            elseif(cCellType .eq. 'FE20') then
               do n1=1, iaNBeta(3), 1
                  rBurnupPercent = rBurnupPercent + 
     &               raBetaST20(n1)*(rBurnupMWD+rOldBurnupMWD)**n1               
               enddo            
            endif

            write(NFO,'(A4,1X,A4,1X,7F10.3)')
     &         cLoading(iEl), cCellType, rMassU, rEnrichment,
     &          rMassEr166,
     &         rMassEr167, rBurnupMWD+rOldBurnupMWD, rBurnupPercent, 
     &         rOldBurnupMWD
            endif
        enddo

!     NEXT WE READ THE POWER DISTRIBUTION FILE (PDIST.OUT) AND STORE
!     DATA IN raCellPower.
         do iCell=1, 6400, 1
            read(NFPPB,*,end=134) cCellID, raCellPower(iCell),
     &                            raCellPowerErr(iCell)
            write(NFPPC,'(A7,2X,G0.6,2X,G0.6)') cCellID,
     &                                raCellPower(iCell)*
     &                                rBurnupIncrement,
     &                                raCellPowerErr(iCell)*
     &                                rBurnupIncrement
         enddo
134      continue
         
         
         print'(A,G0.7,A,G0.7,A)', 
     &           '  BURNUP: Total core burnup after ', 
     &               rBurnupIncrement ,
     &           ' day(s) is ' , rTotalCoreBurnup, ' MWd.'
         print*, ' BURNUP: Finished writing FUEL_INVENTORY.OUT.'


! ---------------------------------------------------------------------------------         
        
         goto 96
94       print*, ' BURNUP: Fatal Error - BURNUP card is not found',
     &           ' in the input file!'
         goto 96
97       print*, ' BURNUP: Fatal Error - BURNUP card parameters ',
     &           'are not valid.'
         goto 96
98       print*, ' BURNUP: Fatal Error - the input file is missing!'
         goto 96
99       print*, ' BURNUP: Fatal Error - FUEL_INVENTORY.INP is missing!'
         goto 96
100      print*, ' BURNUP: Fatal Error - CORECONFIG card parameters ',
     &           ' are not valid.'
         goto 96
101      print*, ' BURNUP: Fatal Error - NRINGS card parameter is',
     &           ' invalid.'
         goto 96
102      print*, ' BURNUP: Fatal Error - PPE.OUT is missing!'
         goto 96
103      print*, ' BURNUP: Fatal Error - Cannot write',
     &           ' FUEL_INVENTORY.OUT!'
         goto 96
104      print*, ' BURNUP: Fatal Error - Cannot read',
     &           ' PDIST.OUT!'
         goto 96
105      print*, ' BURNUP: Fatal Error - Cannot write',
     &           ' ZBURN.OUT!'
         goto 96
         
96       continue

         close(unit=NFI)
         close(unit=NFEL)
         close(unit=NFO)
         close(unit=NFPPE)
         close(unit=NFPPB)
         close(unit=NFPPC)
         
      end subroutine


      
      
! -----------------------------------------------------------------------------
!                          PLEASE READ THIS CAREFULLY
! -----------------------------------------------------------------------------  
!                      S U B R O U T I N E ------ ReadEle
!
!     This subroutine enables us to read the fuel data from ELEM.INP.
!     This subroutine allows the programmer to search for the data of
!     the desired fuel element (searching its Fuel Site Tag, pcFuelTag).
!     Data includes the following:
!
!        * pNFEL - The file index specifying the ELEM.INP file.
!        * pcaElementTags - Characters of fuel element tags.
!        * pcElementType - Fuel element type (ST8, ST12, LEU, ...)
!        * prMassU - Mass of Uranium in the fuel element.
!        * prEnrichment - Fuel enrichment of the fuel element.
!        * prMassEr166 - Mass of Er-166 in the fuel.
!        * prMassEr167 - Mass of Er-167 in the fuel.
!        * prBuMW - Fuel element burnup in MW.
!        * prBuU235 - FUel element burnup of U-235.
!
! -----------------------------------------------------------------------------  
      subroutine ReadEle(pNFEL, pcFuelTag, pcCellType, prMassU, 
     &                   prEnrichment, prMassEr166, prMassEr167,
     &                   prBuMW, prBuU235)


         integer, intent(in) :: pNFEL
         character(len=4), intent(in)  :: pcFuelTag
         real   , intent(out) :: prMassU, prEnrichment, prMassEr166, 
     &                           prMassEr167, prBuMW, prBuU235

         character(len=4), intent(out) :: pcCellType
      
!        cFuelEleRead is just a character variable for reading fuel tag (6574 , 6945 , ...)
         character(len=4) :: cFuelTagRead
      
!        Initialize all parameter values.
         prMassU      = 0.0
         prEnrichment = 0.0
         prMassEr166  = 0.0
         prMassEr167  = 0.0
         prBuMW       = 0.0
         prBuU235     = 0.0
         pcCellType   = 'N/A '
         rewind pNFEL      

         call rbskip(pNFEL,1)

1        continue
         read(pNFEL,100,end=999,err=998) 
     &      cFuelTagRead, pcCellType, prMassU,
     &      prEnrichment, prMassEr166, prMassEr167,
     &      prBuMW, prBuU235
         if(cFuelTagRead .ne. pcFuelTag) goto 1
         return
998      continue
         print*, ' BURNUP: Fatal Error - An error has occured while',
     &              ' reading FUEL_INVENTORY.INP.'
         stop
999      continue
         if(trim(pcFuelTag) .eq. 'W' .or. 
     &      trim(pcFuelTag) .eq. 'G' .or. 
     &      trim(pcFuelTag) .eq. 'COOL' .or.
     &      trim(pcFuelTag) .eq. 'GRAP' .or.
     &      trim(pcFuelTag) .eq. 'BERY' .or.
     &      trim(pcFuelTag) .eq. 'CHN1' .or.
     &      trim(pcFuelTag) .eq. 'CHN2' .or.
     &      trim(pcFuelTag) .eq. 'CHN3' .or.
     &      trim(pcFuelTag) .eq. 'CHN4' .or.
     &      trim(pcFuelTag) .eq. 'TRCR') then  
            pcCellType = trim(pcFuelTag)
            return
         endif      
         print'(A,A4,A1)', ' BURNUP: Warning, missing fuel tag ',
     &           pcFuelTag, '!'
         pcCellType   = 'N/A '     
         return
100      format(A4,1X,A4,1X,6F10.0)

      end subroutine
      
      
      subroutine ReadPPE(pNFPPE, pcSite, prValue)
         character(len=4), intent(in) :: pcSite
         integer, intent(in) :: pNFPPE
         real, intent(out) :: prValue
         character(len=4) :: cSiteRead
         real :: rErrorRead
         prValue = 0.0
         
1        continue
         read(pNFPPE,*,end=999,err=1)
     &      cSiteRead, prValue, rErrorRead
         if(cSiteRead .ne. pcSite) goto 1
999      return
      end subroutine
      
      subroutine rbskip(pNF, pnSkip)
         implicit none
         integer, intent(in) :: pNF, pnSkip
         character(len=1) :: ch
         integer :: i
         if(pnSkip .le. 0) return
         do i=1,pnSkip
            read(pNF,'(A)') ch
         enddo
         return
      end subroutine