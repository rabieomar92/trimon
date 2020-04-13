!     --------------------------------------------------------
!     THIS FORTRAN MODULE FILE WAS GENERATED USING FORTRAN-IT.
!     --------------------------------------------------------
!     Generation time: 6/8/2018 4:54:26 PM
!     Author         : 
!     Device Name    : RABIEMSI
!     Project Name   : CellGeometry.fip
!     --------------------------------------------------------
      
      
!     Begin your code here:
      
      module LIBReaderM
         use ChannelList
         integer :: LIB_COUNT
         real    :: LIB_MASSU0, LIB_MASSFE08
         real    :: LIB_MASSFE12, LIB_MASSFE20
         integer :: LIB_GRP
         integer :: LIB_SLT = 20
         real    :: rTotalPower = 0.0
         integer, parameter :: LIB_FE08  = 1
         integer, parameter :: LIB_FE12  = 2
         integer, parameter :: LIB_FE20  = 3
         integer, parameter :: LIB_FLIP  = 4
         integer, parameter :: LIB_LEUR  = 5
         integer, parameter :: LIB_GRAP  = 6
         integer, parameter :: LIB_COOL  = 7
         integer, parameter :: LIB_G     = 8
         integer, parameter :: LIB_CHN1  = 9
         integer, parameter :: LIB_CHN2  = 10
         integer, parameter :: LIB_CHN3  = 11
         integer, parameter :: LIB_CHN4  = 12
!        LIB_GAMMA STORES THE SLOPE OF SIGMA-VS-BURNUP LINEAR LINE(MWd). 
!        STRUCTURE: LIB_GAMMA(INTERVAL_POINT,ENERGY_GROUP,DATA)
!        LIB_SIGMA STORES THE CROSS SECTIONS AT ZERO BURNUP.
!        STRUCTURE: LIB_SIGMA(INTERVAL_POINT,ENERGY_GROUP,DATA)
         real, allocatable :: LIB_GAMMA(:,:,:,:)
         real, allocatable :: LIB_SIGMA(:,:,:,:)
         real, allocatable :: LIB_POWER(:,:)
         integer :: LIB_INTCOUNT(20)
         
!     Declare your module variables here.
      
      contains
!     Declare your methods here.

      subroutine LIB_READ_HEADER_8(pNF, prMassU0, pnPrincipalGroup,
     &                         pnPrincipalIntCount, piStat)
      
         implicit none
         
         ! Define parameter(s) here...
         integer, intent(in)  :: pNF
         integer, intent(out) :: pnPrincipalGroup,
     &                           pnPrincipalIntCount
         integer, intent(out) :: piStat
         real   , intent(out) :: prMassU0
         integer :: ios
         character(len=80) :: cText
         character(len=10) :: cElem
         
         piStat = 0
         
         rewind pNF
91       read(pNF, '(A80)', iostat=ios, err=92) cText
         if(cText(1:1) .eq. '$') goto 91
         if ( LIB_flIsValid(8, cText) ) then
            read(cText(51:60),'(F10.3)') prMassU0
            read(cText(61:70),'(I10)') pnPrincipalGroup
            read(cText(71:80),'(I10)') pnPrincipalIntCount    
            piStat = 0
            goto 90
         else
            piStat = -1
            goto 91
         endif
         
92       piStat = -2         
90       continue
         return
      end subroutine

      logical function LIB_flIsValid(piCardType, pcTextLine)
         implicit none
         ! 0 = is the 000000000 card, the header card of rthe eactor core layer tape.
         ! 1 = is the 111111111 card, the header card of the reactor core cell record.
         ! 9 = is the 999999999 card, the footer card of the reactor core cell record.
         ! 8 = is the 888888888 card, the header card of the joined tape.
         ! -1 = is the data record row with 8 x F10.7 TXS elements.
         ! -2 = is the data record row with 2 x F10.7 TXS elements.
         
         integer, intent(in) :: piCardType
         character(len=80), intent(in) :: pcTextLine
         character(len=10) :: caElem(8)
         
         caElem(1) = pcTextLine(1:10)
         caElem(2) = pcTextLine(11:20)
         caElem(3) = pcTextLine(21:30)
         caElem(4) = pcTextLine(31:40)
         caElem(5) = pcTextLine(41:50)
         caElem(6) = pcTextLine(51:60)
         caElem(7) = pcTextLine(61:70)
         caElem(8) = pcTextLine(71:80)
         
         LIB_flIsValid = .false.
         
         select case (piCardType)
            case (-8)
               if( (LIB_flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(4)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(5)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(6)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(7)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(8)) .eqv. .true.)) then
                  
                  LIB_flIsValid = .true.
                  return
                  
               endif
               
            case(-1)
               if( (LIB_flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (caElem(2) .eq. '          ') .and.
     &             (caElem(3) .eq. '          ') .and.
     &             (caElem(4) .eq. '          ') .and.
     &             (caElem(5) .eq. '          ') .and.
     &             (caElem(6) .eq. '          ') .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  LIB_flIsValid = .true.
                  return
                  
               endif
            case(-2)
               if( (LIB_flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (caElem(3) .eq. '          ') .and.
     &             (caElem(4) .eq. '          ') .and.
     &             (caElem(5) .eq. '          ') .and.
     &             (caElem(6) .eq. '          ') .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  LIB_flIsValid = .true.
                  return
                  
               endif
               
            case(-3)
               if( (LIB_flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (caElem(4) .eq. '          ') .and.
     &             (caElem(5) .eq. '          ') .and.
     &             (caElem(6) .eq. '          ') .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  LIB_flIsValid = .true.
                  return
                  
               endif
            case(-4)
               if( (LIB_flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(4)) .eqv. .true.) .and.
     &             (caElem(5) .eq. '          ') .and.
     &             (caElem(6) .eq. '          ') .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  LIB_flIsValid = .true.
                  return
                  
               endif
            case(-5)
               if( (LIB_flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(4)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(5)) .eqv. .true.) .and.
     &             (caElem(6) .eq. '          ') .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  LIB_flIsValid = .true.
                  return
                  
               endif
               
            case(-6)
               if( (LIB_flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(4)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(5)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(6)) .eqv. .true.) .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  LIB_flIsValid = .true.
                  return
                  
               endif
               
            case(-7)
               if( (LIB_flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(4)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(5)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(6)) .eqv. .true.) .and.
     &             (LIB_flIsReal(caElem(7)) .eqv. .true.) .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  LIB_flIsValid = .true.
                  return
                  
               endif
            case(0)
               if( (caElem(1) .eq. ' 000000000') .and.
     &             (caElem(2) .eq. ' 000000000') .and.
     &             (caElem(3) .eq. ' 000000000') .and.
     &             (caElem(4) .eq. ' 000000000') .and.
     &             (caElem(5) .eq. ' 000000000') .and.
     &             (LIB_flIsInteger( caElem(6) ) .eqv. .true.) .and.
     &             (LIB_flIsInteger( caElem(7) ) .eqv. .true.) .and.
     &             (LIB_flIsInteger( caElem(8) ) .eqv. .true.) ) then

                  LIB_flIsValid = .true.
                  return
                  
               endif
             
            case(1)
               if( (caElem(1) .eq. ' 111111111') .and.
     &             (caElem(2) .eq. ' 111111111') .and.
     &             (caElem(3) .eq. ' 111111111') .and.
     &             (caElem(4) .eq. ' 111111111') .and.
     &             (caElem(5) .eq. ' 111111111') .and.
     &             (caElem(6) .eq. ' 111111111') .and.
     &             (LIB_flIsInteger( caElem(7) ) .eqv. .true.) .and.
     &             (LIB_flIsInteger( caElem(8) ) .eqv. .true.) ) then         

                  LIB_flIsValid = .true.
                  return
             
               endif  
             
            case(8)
               if( (caElem(1) .eq. ' 888888888') .and.
     &             (caElem(2) .eq. ' 888888888') .and.
     &             (caElem(3) .eq. ' 888888888') .and.
     &             (caElem(4) .eq. ' 888888888') .and.
     &             (caElem(5) .eq. ' 888888888') .and.
c     &             (LIB_flIsReal( caElem(6) ) .eqv. .true.) .and.
     &             (LIB_flIsInteger( caElem(7) ) .eqv. .true.) .and.
     &             (LIB_flIsInteger( caElem(8) ) .eqv. .true.) ) then          

                  LIB_flIsValid = .true.
                  return
             
               endif     

            case(9)
               if( (caElem(1) .eq. ' 999999999') .and.
     &             (caElem(2) .eq. ' 999999999') .and.
     &             (caElem(3) .eq. ' 999999999') .and.
     &             (caElem(4) .eq. ' 999999999') .and.
     &             (caElem(5) .eq. ' 999999999') .and.
     &             (caElem(6) .eq. ' 999999999') .and.
     &             (LIB_flIsInteger( caElem(7) ) .eqv. .true.) .and.
     &             (LIB_flIsInteger( caElem(8) ) .eqv. .true.) ) then           

                  LIB_flIsValid = .true.
                  return
             
               endif    
            
            case default
            
         end select
         return
      end function

      logical function LIB_flIsANumber(pcChar)
         implicit none
         character, intent(in) :: pcChar
         LIB_flIsANumber = .false.
         if (pcChar .eq. '0') LIB_flIsANumber = .true.
         if (pcChar .eq. '1') LIB_flIsANumber = .true.
         if (pcChar .eq. '2') LIB_flIsANumber = .true.
         if (pcChar .eq. '3') LIB_flIsANumber = .true.
         if (pcChar .eq. '4') LIB_flIsANumber = .true.
         if (pcChar .eq. '5') LIB_flIsANumber = .true.
         if (pcChar .eq. '6') LIB_flIsANumber = .true.
         if (pcChar .eq. '7') LIB_flIsANumber = .true.
         if (pcChar .eq. '8') LIB_flIsANumber = .true.
         if (pcChar .eq. '9') LIB_flIsANumber = .true.   
         return
      end function

      logical function LIB_flIsInteger(pcTXSElement)
         implicit none
         character(len=10), intent(in) :: pcTXSElement
         integer :: i, nDot
         LIB_flIsInteger = .true.
         if ( pcTXSElement(1:1) .ne. ' ') then
               LIB_flIsInteger = .false.
               return
         endif
         do i=2, 10, 1
            if ( LIB_flIsANumber(pcTXSElement(i:i)) .eqv. .false.) then
               LIB_flIsInteger = .false.
               return
            endif
         enddo
         return
      end function

      logical function LIB_flIsReal(pcTXSElement)
         implicit none
         character(len=10), intent(in) :: pcTXSElement
         integer :: i, nDot
         
         ! A read numbered element should only contains one (1) dot (.)
         logical :: lCondition0Passed = .false.
         ! The first char may be a blank space, or a minus sign (-), or a number [0-9].
         logical :: lCondition1Passed = .false.
         ! The second char MUST be a number [0-9].
         logical :: lCondition2Passed = .false.
         ! The third char MUST be a dot.
         logical :: lCondition3Passed = .false.
         ! The fourth and so on should be a number.
         logical :: lCondition4Passed = .false.
         
         ! Checking the zero-th condition...
         nDot = 0
         do i=1, 10, 1
            if (pcTXSElement(i:i) .eq. '.') nDot = nDot + 1
         enddo
         if (nDot .eq. 1) then
            lCondition0Passed = .true.
         else
            LIB_flIsReal = .false.
            return
         endif
         
         ! Checking the first condition...
         if (LIB_flIsANumber(pcTXSElement(1:1)) .eqv. .true.) then
            lCondition1Passed = .true.
         elseif (pcTXSElement(1:1) .eq. ' ') then
            lCondition1Passed  = .true.
         elseif (pcTXSElement(1:1) .eq. '-') then
            lCondition1Passed = .true.
         else
            LIB_flIsReal = .false.
            return                  
         endif
         
         ! Checking the second condition...
         if (LIB_flIsANumber(pcTXSElement(2:2)) .eqv. .true.) then
            lCondition2Passed = .true.
         else
            LIB_flIsReal = .false.
            return                  
         endif
         
         ! Checking the third condition...
         if (pcTXSElement(3:3) .eq. '.') then
            lCondition3Passed = .true.
         else
            LIB_flIsReal = .false.
            return         
         endif
         
         ! Checking the fourth condition...
         do i=4, 10, 1
            if (LIB_flIsANumber(pcTXSElement(i:i)) .eqv. .true.) then
               lCondition4Passed = .true.
            else
               LIB_flIsReal = .false.
               return                  
            endif            
         enddo
         
         LIB_flIsReal = lCondition1Passed .and. 
     &              lCondition2Passed .and.
     &              lCondition3Passed .and.
     &              lCondition4Passed
         return
      end function
      
      subroutine INIT_LIB_READER()
         integer :: i, j, k, l
         

            allocate( LIB_GAMMA(20,1000, 16,
     &                16+4) )  
            allocate( LIB_SIGMA(20,1000, 16,
     &                16+4) )
            allocate( LIB_POWER(20,1000) )
            
            do l=1, 20, 1
            do i=1, 1000, 1
               do j=1, 16, 1
                  do k=1, 16+4, 1
                     LIB_GAMMA(l,i,j,k) = 0.0
                     LIB_SIGMA(l,i,j,k) = 0.0
                  enddo
               enddo
               LIB_POWER(l,i) = 0.0
            enddo
            enddo
            return
      end subroutine

      subroutine LIB_GET_TABLE(pNF, piSlotID, piStat)
         implicit none
         
         ! Define parameter(s) here...
         integer, intent(in) :: pNF, piSlotID
         integer, intent(out) :: piStat
         integer :: iCurrentLine = 0
         integer :: nPrincipalGroup, nPrincipalIntCount
         integer :: iGroup, iInterval, j
         integer ::  iIntervalT
         real    :: rPowerT
         integer :: iLine = 0
         character(len=80) :: cText
         character(len=10) :: iDummy
         
!        INITIALLY WE ASSUME THERE IS NO ERROR. WE SET piStat TO 0.
         piStat = 0


!        NEXT WE READ THE TXS TAPE HEADER-8. HEADER-8 CONTAINS THE INFOR-
!        MATION ABOUT THE TOTAL NUMBER OF NEUTRON GROUPS,  THE TOTAL NUM-
!        BER OF CELLS AND THE TOTAL NUMBER OF REACTOR CORE LAYERS.
!        EXAMPLE OF HEADER-8:
!        888888888 888888888 888888888 888888888   NGROUP  NCELL  NLAYER
         call LIB_READ_HEADER_8(pNF, LIB_MASSU0, nPrincipalGroup,
     &             nPrincipalIntCount, piStat)
         iLine = iLine + 1
         if(piStat .lt. 0) then
            print*, ' LIBR: Fatal error. Failed to ',
     &              'read data table. '
            print*, '           Invalid HEADER-8 format.'
            print'(A,I0,A)', '            (Line ', 
     &            iLine, ')'
            stop
         endif
         
!        CHECK IF THE DIMENSION OF THE TABLE ARRAY SIZE PROVIDED BY  THE
!        INVOKER. IF THE ARRAY SIZE IS SMALL, THEN STOP. TELL USER THERE
!        IS INTERNAL ERROR. PROGRAMMERS MUST TAKE NOTE THAT THE SIZE  OF
!        TXS_TABLE HAS TO BE (pnL,pnC,pnG,10). 
         LIB_INTCOUNT(piSlotID) = nPrincipalIntCount
         LIB_GRP = nPrincipalGroup
         write(iDummy,'(I10)') nPrincipalGroup



!        WE IMPLEMENT NESTED LOOPS HERE TO  ITERATE THROUGH EVERY SINGLE
!        REACTOR CORE CELL AVAILABLE IN THE TAPE. THE ITERATIONS ARE LI-
!        MITTED UP TO MAX NUMBER OF CORE LAYERS (nTotalLayers)  AND  MAX
!        NUMBER OF CELLS - THE PRINCIPAL CELL COUNT (nPrincipalCell).
         !$OMP PARALLEL DO



!        HERE WE READ THE POWER INTERVAL LISTING..         
         do iInterval=1, nPrincipalIntCount, 1
95          iLine = iLine + 1
            read(pNF,'(A80)') cText
            if(cText(1:1) .eq. '$') goto 95 ! THIS IS FOR IGNORING COMMENTS
            if(cText(2:10) .eq. 'INT-CHECK') then
               read(cText(11:20),'(I10)',end=999,err=999) iIntervalT
               if(iInterval .eq. iIntervalT) then
                     read(cText(21:30),'(F10.0)',end=999,err=999)
     &                      rPowerT
                     if(iInterval .gt. 1) then
                        if(rPowerT .gt. LIB_POWER(piSlotID,
     &                      iInterval-1)) then
                            LIB_POWER(piSlotID, iInterval) = rPowerT
                        else
                           print*, ' LIBR: Fatal error, ', 
     &                         'improper power interval sequence.'  
                           print'(A,I0,A)', '              (Line ', 
     &                                      iLine, ')'
                           stop
                        endif
                     else
                         LIB_POWER(piSlotID,iInterval) = rPowerT
                     endif
               else
                  print*, ' LIBR: Fatal error, ', 
     &                    'library is corrupted.'
                  print'(A,I0,A)', '              (Line ', 
     &                                      iLine, ')'
                  stop
               endif
            else
               print*, ' LIBR: Fatal error, ',
     &                  ' library is corrupted.'
               print'(A,I0,A)', '              (Line ', 
     &                    iLine, ')'
               stop
            endif
         enddo
         
         do iInterval=1, nPrincipalIntCount, 1
         
903         iLine = iLine + 1
            read(pNF,'(A80)',end=999,err=999) cText
            if(cText(1:1) .eq. '$') goto 903 ! THIS IS FOR IGNORING COMMENTS
            if(LIB_flIsValid(1, cText)) then
               read(cText(71:80),'(I10)') iIntervalT
               if(iInterval .eq. iIntervalT) then
                  goto 904
               else
                  goto 903
               endif
            else
               goto 903
            endif            

904      continue

!        FOR EACH GROUP DEFINED FOR THE CELL, WE READ ALL COLUMNS AND
!        ASSIGN THEM TO THE TABLE MATRIX.    
         do iGroup=1, nPrincipalGroup, 1
            iLine = iLine + 1
            read(pNF,'(8F10.7)',end=999,err=999) 
     &               (LIB_GAMMA(piSlotID,iInterval, iGroup, j),
     &                j=1, nPrincipalGroup+4)
         enddo
         
         do iGroup=1, nPrincipalGroup, 1
            iLine = iLine + 1
            read(pNF,'(8F10.7)',end=999,err=999) 
     &            (LIB_SIGMA(piSlotID,iInterval, iGroup, j),
     &             j=1, nPrincipalGroup+4)
         enddo
         
         enddo ! END POWER INTERVAL LOOP
         
         return

999   print*,         ' LIBR: Fatal error. The provided',
     &                ' tape is corrupted!'
      print'(A,I0,A)', '            (Line ', iLine, ')'
      stop
      end subroutine      
      
      real function LIB_GETELEMENT(piSlot, piColumn, 
     &         piLayer, piCell, piGrp)
  

         integer, intent(in) :: piSlot, piColumn
         integer, intent(in) :: piLayer, piCell, piGrp
         integer :: iFuelID, nCell, iType, iInterval
         integer :: iFuelIDT, iRingNumber, iSiteNumber
         real :: rBUSlope, rSigma0, rSigma, rPower, rBurnup
         real :: rRPFF, rZPFF, rPFF, rTotalMass, rM, rC0
         character(len=4) :: caSiteID
         integer :: i, j, k
         caSiteID = CHL_SITEID(piCell)
         if(CHL_NRINGS .eq. 6) then
            nCell = 92
         else
            nCell = 128
         endif
         if((piSlot .lt. 1) .or. (piSlot .gt. 20)) then
            LIB_GETELEMENT = NaN
            return
         endif
         if((piColumn .lt. 1) .or. (piColumn .gt. (LIB_GRP+4))) then
            LIB_GETELEMENT = NaN
            return
         endif
         if((piLayer .lt. 1) .or. (piLayer .gt. CHL_NLAYERS)) then
            LIB_GETELEMENT = NaN
            return
         endif
         if(CHL_NRINGS .eq. 6) then
            if((piCell .lt. 1) .or. (piCell .gt. 92)) then
               LIB_GETELEMENT = NaN
               return
            endif
         elseif(CHL_NRINGS .eq. 7) then
            if((piCell .lt. 1) .or. (piCell .gt. 128)) then
               LIB_GETELEMENT = NaN
               return
            endif
         else
            LIB_GETELEMENT = NaN
            return         
         endif
         
         if(CHL_FUELID(piCell) .eq. 'GRAP') then
            rSigma = LIB_SIGMA(LIB_GRAP,1, piGrp, piColumn)
            goto 99
         endif
         if(CHL_FUELID(piCell) .eq. 'COOL') then
            rSigma = LIB_SIGMA(LIB_COOL,1, piGrp, piColumn)
            goto 99         
         endif
         if(CHL_FUELID(piCell) .eq. 'BERY') then
            rSigma = LIB_SIGMA(LIB_BERY,1, piGrp, piColumn)
            goto 99         
         endif
         if(CHL_FUELID(piCell) .eq. 'G   ') then
            rSigma = LIB_SIGMA(LIB_G,1, piGrp, piColumn)
            goto 99         
         endif
         if(CHL_FUELID(piCell) .eq. 'CHN1') then
            rSigma = LIB_SIGMA(LIB_CHN1,1, piGrp, piColumn)
            goto 99         
         endif
         if(CHL_FUELID(piCell) .eq. 'CHN2') then
            rSigma = LIB_SIGMA(LIB_CHN2,1, piGrp, piColumn)
            goto 99         
         endif
         if(CHL_FUELID(piCell) .eq. 'CHN3') then
            rSigma = LIB_SIGMA(LIB_CHN3,1, piGrp, piColumn)
            goto 99         
         endif
         if(CHL_FUELID(piCell) .eq. 'CHN4') then
            rSigma = LIB_SIGMA(LIB_CHN4,1, piGrp, piColumn)
            goto 99         
         endif        
         
         read(CHL_FUELID(piCell),'(I4)') iFuelID
         
         ! HERE WE OBTAIN PFF
         if(CHL_SITEID(piCell)(1:1) .eq. 'A') then
            rRPFF = CHL_RPFF(1)
         elseif(CHL_SITEID(piCell)(1:1) .eq. 'B') then
            rRPFF = CHL_RPFF(2)
         elseif(CHL_SITEID(piCell)(1:1) .eq. 'C') then
            rRPFF = CHL_RPFF(3)
         elseif(CHL_SITEID(piCell)(1:1) .eq. 'D') then
            rRPFF = CHL_RPFF(4)
         elseif(CHL_SITEID(piCell)(1:1) .eq. 'E') then
            rRPFF = CHL_RPFF(5)
         elseif(CHL_SITEID(piCell)(1:1) .eq. 'F') then
            rRPFF = CHL_RPFF(6)
         elseif(CHL_SITEID(piCell)(1:1) .eq. 'G') then
            rRPFF = CHL_RPFF(7)
         endif
         
         rZPFF = CHL_ZPFF(piLayer)
         rPFF = rRPFF * rZPFF
         
         k = 0
         rTotalMass = 0.0
         do i=1, nCell, 1
            read(CHL_FUELID(i),'(I4)',err=663) iFuelIDT
            rTotalMass = rTotalMass + CHL_FUELMASS(iFuelIDT) * 0.199 
            k = k + 1
663         continue
         enddo
         rPower = rPFF * (CHL_FUELMASS(iFuelID) * 0.199 / rTotalMass) *
     &         CHL_COREPOWER
         
     
         iType = CHL_FUELTYPE(iFuelID)
         
         do i=1, LIB_INTCOUNT(iType)-1, 1
            if((rPower .ge. LIB_POWER(iType,i)) .and.
     &         (rPower .lt. LIB_POWER(iType,i+1))) then
               iInterval = i+1
               goto 70
            endif
         enddo
         
         
70       continue
         
         rS1 = LIB_SIGMA(iType,iInterval+1,piGrp,piColumn)
         rS0 = LIB_SIGMA(iType,iInterval,piGrp,piColumn)
         rP1 = LIB_POWER(iType,i+1)
         rP0 = LIB_POWER(iType,i)
         rM = (rS1 - rS0) / (rP1 - rP0)
         rC0 = rS1 - rM * rP1

         rSigma0 = rM * rPower + rC0
         
         if(caSiteID(1:1) .eq. 'A') then
            iRingNumber = CHL_RINGA
         elseif(caSiteID(1:1) .eq. 'B') then
            iRingNumber = CHL_RINGB
         elseif(caSiteID(1:1) .eq. 'C') then
            iRingNumber = CHL_RINGC
         elseif(caSiteID(1:1) .eq. 'D') then
            iRingNumber = CHL_RINGD
         elseif(caSiteID(1:1) .eq. 'E') then
            iRingNumber = CHL_RINGE
         elseif(caSiteID(1:1) .eq. 'F') then
            iRingNumber = CHL_RINGF
         elseif(caSiteID(1:1) .eq. 'G') then
            iRingNumber = CHL_RINGG
         endif
         
         read(caSiteID(3:4), '(I2)') iSiteNumber
         
         rBurnup = CHL_LOCALBURNUP(iRingNumber, iSiteNumber,
     &                             piLayer)
         rBurnup = CHL_FUELBURNUPPERCENT(iFuelID)
         rBurnup = (100.0-rBurnup)*0.01*CHL_FUELMASS(iFuelID)*0.199
         rGamma = LIB_GAMMA(iType,iInterval,piGrp,piColumn)
         if(rGamma .ne. rGamma) then
            rSigma = rSigma0         
         else
            if(iType .eq. LIB_FE08) then
               rC0 = 1.0 - rGamma*LIB_MASSFE08*0.199
            elseif(iType .eq. LIB_FE12) then
               rC0 = 1.0 - rGamma*LIB_MASSFE12*0.199
            elseif(iType .eq. LIB_FE12) then
               rC0 = 1.0 - rGamma*LIB_MASSFE20*0.199
            endif
            rSigma = (rC0 + rGamma * rBurnup) * rSigma0
         endif

99       LIB_GETELEMENT = rSigma
         return
      end function
       
      end module
