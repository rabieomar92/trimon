
      module ChannelList
         integer, parameter :: CHL_FE08  = 1
         integer, parameter :: CHL_FE12  = 2
         integer, parameter :: CHL_FE20  = 3
         integer, parameter :: CHL_FLIP  = 4
         integer, parameter :: CHL_LEUR  = 5
         
         integer, parameter :: CHL_RINGA = 21
         integer, parameter :: CHL_RINGB = 22
         integer, parameter :: CHL_RINGC = 23
         integer, parameter :: CHL_RINGD = 24
         integer, parameter :: CHL_RINGE = 25
         integer, parameter :: CHL_RINGF = 26
         integer, parameter :: CHL_RINGG = 27
         
         logical :: CHL_ISLOCALBURNUP = .false.
         character(len=4) :: CHL_FUELID(128)
         character(len=4) :: CHL_SITEID(128)
         real, allocatable :: CHL_RPFF(:)
         real, allocatable :: CHL_ZPFF(:)
         
         real    :: CHL_FUELBURNUP(1:9999)
         real    :: CHL_LOCALBURNUP(30,40,30)
         integer :: CHL_FUELINVENTORY(1:9999)
         real    :: CHL_FUELBURNUPPERCENT(1:9999)
         real    :: CHL_FUELMASS(1:9999)
         integer :: CHL_FUELTYPE(1:9999)
         integer :: CHL_NRINGS, CHL_NLAYERS
         real    :: CHL_COREPOWER = 0.01
      contains

         subroutine CHL_ReadMainInput(plReadLocalBurnup)
         implicit none
            logical, intent(in) :: plReadLocalBurnup
            character(len=100) :: cText
            integer :: nRings, nLayers, iLine
            logical :: lNRINGSFound = .false.
            logical :: lNLAYERSFound = .false.
            logical :: lPFFFound = .false.
            logical :: lPowerFound = .false.
            real    :: rPower
            integer :: i,j,k,l, n
            open(unit=101, file='MAIN.INP',status='OLD', err=911)
            open(unit=102, file='FUEL_INVENTORY.INP',
     &           status='OLD',err=913)
            if(plReadLocalBurnup .eqv. .true.) then
               open(unit=103, file='ZBURN.OUT', 
     &              status='OLD', err=917)
               CHL_ISLOCALBURNUP = .true.
               call CHL_ReadLocalBurnup(103)
               goto 918
            else
               goto 917
            endif
913         print*, ' LIBR: Fatal error, fuel inventory',
     &              ' file is missing!'
            stop
            ! HERE WE READ THE FUEL INVENTORY LIST
917         continue  
            CHL_ISLOCALBURNUP = .false.
918         continue
            ! HERE WE READ THE FUEL INVENTORY LIST            
            call CHL_ReadFuelInventory(102)
            ! WE SEARCH FOR NRINGS CARD
            do iLine=1, 1000, 1
               read(101,'(A100)',end=912,err=912) cText
               if(cText(1:7) .eq. '#NRINGS') then           
                  read(101,*,err=912,end=912) nRings
                  CHL_NRINGS = nRings
                  allocate(CHL_RPFF(nRings))
                  lNRINGSFound = .true.
                  goto 91
               endif
            enddo
            
91          continue

            if(lNRINGSFound .eqv. .false.) then
               print*, ' LIBR: NRINGS card is missing',
     &                 ' in MAIN.INP!'
               stop
            else
               rewind 101
            endif
             
            ! WE SEARCH FOR NLAYERS CARD
            do iLine=1, 1000, 1
               read(101,'(A100)',end=912,err=912) cText
               if(cText(1:8) .eq. '#NLAYERS') then
                  read(101,*,err=912,end=912) nLayers
                  allocate(CHL_ZPFF(nLayers))
                  CHL_NLAYERS = nLayers
                  lNLAYERSFound = .true.
                  goto 92
               endif
            enddo

92          continue
            
            if(lNLAYERSFound .eqv. .false.) then
               print*, ' LIBR: NLAYERS card is missing',
     &                 ' in MAIN.INP!'
               stop
            else
               rewind 101
            endif
            
            ! WE SEARCH FOR POWER CARD
            do iLine=1, 1000, 1
               read(101,'(A100)',end=912,err=912) cText
               if(cText(1:6) .eq. '#POWER') then
                  read(101,*,err=912,end=912) rPower
                  CHL_COREPOWER = rPower
                  lPOWERFound = .true.
                  goto 96
               endif
            enddo

96          continue
            
            if(lPOWERFound .eqv. .false.) then
               print*, ' LIBR: POWER card is missing',
     &                 ' in MAIN.INP!'
               stop
            else
               rewind 101
            endif
            
            ! WE SEARCH FOR PFF CARD
            do iLine=1, 1000, 1
               read(101,'(A100)',end=912,err=912) cText
               if(cText(1:4) .eq. '#PFF') then
                  read(101,*,err=915,end=915) 
     &               (CHL_RPFF(i), i=1, nRings)
                  read(101,*,err=915,end=915) 
     &               (CHL_ZPFF(i), i=1, nLayers)                 
                  lPFFFound = .true.
                  goto 98
               endif
            enddo

915         print*, ' LIBR: Invalid PFF card format.'
            stop
98          continue

            
            if(lPFFFound .eqv. .false.) then
               print*, ' LIBR: PFF card is missing',
     &                 ' in MAIN.INP!'
               stop
            else
               rewind 101
            endif 
            
            ! WE SEARCH FOR CORECONFIG CARD
            do iLine=1, 1000, 1
               read(101,'(A100)', end=912, err=912) cText
               if(cText(1:11) .eq. '#CORECONFIG') then
                  n = 1
                  read(101,'(A4,1X,A4)',err=13, end=13)
     &                 CHL_SITEID(n), CHL_FUELID(n)
                       ! VALIDATE
                       if(CHL_IsValidSiteID(CHL_SITEID(n))
     &                      .eqv. .false.) then
                          goto 13
                       endif
                       if(CHL_IsValidFuelID(CHL_FUELID(n))
     &                      .eqv. .false.) then
                          goto 13
                       endif
                  do i=1, 15, 1
                     read(101,'(A4,1X,A4,5(1X,A4,1X,A4))',
     &                    err=13, end=13)
     &                  (CHL_SITEID(n+j), CHL_FUELID(n+j),j=1, 6)
                     ! VALIDATE
                     do j=1, 6, 1
                        if(CHL_IsValidSiteID(CHL_SITEID(n+j))
     &                      .eqv. .false.) then
                           goto 13
                        endif
                       if(CHL_IsValidFuelID(CHL_FUELID(n+j))
     &                      .eqv. .false.) then
                          goto 13
                       endif
                     enddo
                     n = n + 6
                  enddo
                  if(nRings .eq. 6) then
                     n = n + 1
                     CHL_SITEID(n) = 'REFL'
                     CHL_FUELID(n) = 'G   '
                  elseif(nRIngs .eq. 7) then
                     do j=1, 6, 1
                        read(101,'(A4,1X,A4,5(1X,A4,1X,A4))',
     &                       err=13, end=13)
     &                      (CHL_SITEID(n+i), CHL_FUELID(n+i),i=1,6)
                        do i=1, 6, 1
                           if(CHL_IsValidSiteID(CHL_SITEID(n+i))
     &                         .eqv. .false.) then
                              goto 13
                           endif
                           if(CHL_IsValidFuelID(CHL_FUELID(n+i))
     &                      .eqv. .false.) then
                              goto 13
                           endif
                        enddo
                        n = n + 6
                     enddo
                     n = n + 1
                     CHL_SITEID(n) = 'REFL'
                     CHL_FUELID(n) = 'G   '                
                  endif
                  
                  ! Compute total U235 Mass
                  goto 93
               endif               
            enddo
            close(unit = 101)
            close(unit = 102)
12          print*, ' LIBR: Missing CORECONFIG card.'
            stop
13          print*, ' LIBR: Invalid CORECONFIG card pattern.'
            stop
93          continue
            close(unit = 101)
            close(unit = 102)
            return
911         print*, ' LIBR: Fatal error. Failed to open MAIN.INP'
            stop
912         print*, ' LIBR: Fatal error. MAIN.INP is corrupted.'
            stop
         end subroutine
         
!        ------------------------------------------------------------
!        SITE ID VALIDATION.
!        ------------------------------------------------------------
         logical function CHL_IsValidSiteID(pcSiteID)
            character(len=4) :: pcSiteID
            integer :: iDummy
            CHL_IsValidSiteID = .false.
            if((pcSiteID(1:1) .eq. 'A') .or.
     &         (pcSiteID(1:1) .eq. 'B') .or.
     &         (pcSiteID(1:1) .eq. 'C') .or.
     &         (pcSiteID(1:1) .eq. 'D') .or.
     &         (pcSiteID(1:1) .eq. 'E') .or.
     &         (pcSiteID(1:1) .eq. 'F') .or.
     &         (pcSiteID(1:1) .eq. 'G') ) then
               CHL_IsValidSiteID = .true.
            else
               CHL_IsValidSiteID = .false.
               return
            endif
            
            if(pcSiteID(2:2) .eq. '-') then
               CHL_IsValidSiteID = .true.
            else
               CHL_IsValidSiteID = .false.
               return
            endif
            
            read(pcSiteID(3:4),'(I2)',err=15,end=15) iDummy
            if(pcSiteID(1:1) .eq. 'A') then
               if(iDummy .eq. 1) then
                  CHL_IsValidSiteID = .true.
               else
                  CHL_IsValidSiteID = .false.
                  return
               endif
            endif
            if(pcSiteID(1:1) .eq. 'B') then
               if((iDummy .ge. 1) .and. (iDummy .le. 6)) then
                  CHL_IsValidSiteID = .true.
               else
                  CHL_IsValidSiteID = .false.
                  return
               endif
            endif  
            if(pcSiteID(1:1) .eq. 'C') then
               if((iDummy .ge. 1) .and. (iDummy .le. 12)) then
                  CHL_IsValidSiteID = .true.
               else
                  CHL_IsValidSiteID = .false.
                  return
               endif
            endif   
            if(pcSiteID(1:1) .eq. 'D') then
               if((iDummy .ge. 1) .and. (iDummy .le. 18)) then
                  CHL_IsValidSiteID = .true.
               else
                  CHL_IsValidSiteID = .false.
                  return
               endif
            endif
            if(pcSiteID(1:1) .eq. 'E') then
               if((iDummy .ge. 1) .and. (iDummy .le. 24)) then
                  CHL_IsValidSiteID = .true.
               else
                  CHL_IsValidSiteID = .false.
                  return
               endif
            endif 
            if(pcSiteID(1:1) .eq. 'F') then
               if((iDummy .ge. 1) .and. (iDummy .le. 30)) then
                  CHL_IsValidSiteID = .true.
               else
                  CHL_IsValidSiteID = .false.
                  return
               endif
            endif 
            if(pcSiteID(1:1) .eq. 'G') then
               if((iDummy .ge. 1) .and. (iDummy .le. 36)) then
                  CHL_IsValidSiteID = .true.
               else
                  CHL_IsValidSiteID = .false.
                  return
               endif
            endif
            
            return
15          CHL_IsValidSiteID = .false.
            return
         end function
         
         logical function CHL_IsValidFuelID(pcFuelID)
            character(len=4), intent(in) :: pcFuelID
            integer :: iDummy
            CHL_IsValidFuelID = .false.
            
            read(pcFuelID,'(I4)',err=16,end=16) iDummy
            if((iDummy .gt. 0) .and. (iDummy .le. 9999)) then
               CHL_IsValidFuelID = .true.
               return
            endif
            
16          CHL_IsValidFuelID = .false.
            if((pcFuelID .eq. 'CHN1') .or.
     &         (pcFuelID .eq. 'CHN2') .or.
     &         (pcFuelID .eq. 'CHN3') .or.
     &         (pcFuelID .eq. 'CHN4') .or.
     &         (pcFuelID .eq. 'GRAP') .or.
     &         (pcFuelID .eq. 'COOL') .or.
     &         (pcFuelID .eq. 'BERY') ) then
               CHL_IsValidFuelID = .true.
               return
            else
               CHL_IsValidFuelID = .false.
               return
            endif
            return
         end function
         
         subroutine CHL_ReadFuelInventory(pNF)
            integer, intent(in) :: pNF
            character(len=4) :: caFuelID
            character(len=4) :: caFuelType
            real :: rDummy(6)
            integer :: iFuelID
            integer :: iDummy, i
            character(len=80) :: caDummy
            read(pNF,'(A80)') caDummy
            read(pNF,'(A80)') caDummy
201         read(pNF,*,err=200,end=202)
     &         caFuelID, caFuelType, (rDummy(i),i=1,6)

            if(CHL_IsValidFuelID(caFuelID) .eqv. .true.) then
               read(caFuelID,'(I4)',err=201,end=201) iFuelID
               if(caFuelType .eq. 'FE08') then
                  CHL_FUELINVENTORY(iFuelID) = CHL_FE08
                  CHL_FUELBURNUP(iFuelID)    = rDummy(5)
                  CHL_FUELBURNUPPERCENT(iFuelID) = rDummy(6)
                  CHL_FUELMASS(iFuelID) = rDummy(1)
                  CHL_FUELTYPE(iFuelID) = CHL_FE08 
               elseif(caFuelType .eq. 'FE12') then
                  CHL_FUELINVENTORY(iFuelID) = CHL_FE12
                  CHL_FUELBURNUP(iFuelID)    = rDummy(5)
                  CHL_FUELBURNUPPERCENT(iFuelID) = rDummy(6)
                  CHL_FUELMASS(iFuelID) = rDummy(1)
                  CHL_FUELTYPE(iFuelID) = CHL_FE12
               elseif(caFuelType .eq. 'FE20') then
                  CHL_FUELINVENTORY(iFuelID) = CHL_FE20
                  CHL_FUELBURNUP(iFuelID)    = rDummy(5)
                  CHL_FUELBURNUPPERCENT(iFuelID) = rDummy(6)
                  CHL_FUELMASS(iFuelID) = rDummy(1)
                  CHL_FUELTYPE(iFuelID) = CHL_FE20
               elseif(caFuelType .eq. 'FLIP') then
                  CHL_FUELINVENTORY(iFuelID) = CHL_FLIP
                  CHL_FUELBURNUP(iFuelID)    = rDummy(5)
                  CHL_FUELBURNUPPERCENT(iFuelID) = rDummy(6)
                  CHL_FUELMASS(iFuelID) = rDummy(1)
                  CHL_FUELTYPE(iFuelID) = CHL_FLIP
               elseif(caFuelType .eq. 'LEUR') then
                  CHL_FUELINVENTORY(iFuelID) = CHL_LEUR
                  CHL_FUELBURNUP(iFuelID)    = rDummy(5)
                  CHL_FUELBURNUPPERCENT(iFuelID) = rDummy(6)
                  CHL_FUELMASS(iFuelID) = rDummy(1)
                  CHL_FUELTYPE(iFuelID) = CHL_LEUR
               endif
               
               goto 201
            else
               goto 201
            endif

200         print*, ' LIBR: Invalid fuel inventory table format!'
            close(unit=pNF)
            stop
            
202         close(unit=pNF)
            return
         end subroutine
         
         subroutine CHL_ReadLocalBurnup(pNF)
            integer, intent(in) :: pNF
            character(len=4) :: caSiteID
            character(len=1) :: cDummy
            character(len=2) :: caLayerID
            real :: rDummy(2)
            integer :: iDummy, i
            integer :: iLayerID, iSiteNumber, iRingNumber
            
            read(pNF,*) cDummy
361         read(pNF,'(A4,A1,A2,G10.6)',err=35,end=32) 
     &         caSiteID, cDummy, caLayerID,
     &         (rDummy(i), i=1, 2)
            goto 423

35          backspace pNF
            read(pNF,'(A4,A1,A2,G10.6E2)',err=31,end=32) 
     &         caSiteID, cDummy, caLayerID,
     &         (rDummy(i), i=1, 2)          
423         continue     
            if(CHL_IsValidSiteID(caSiteID) .eqv. .true.) then
               read(caLayerID,'(I2)',err=31,end=31) iLayerID
               if((iLayerID .gt. 0) .and. (iLayerID .lt. 50)) then
                  if(caSiteID(1:1) .eq. 'A') then
                     iRingNumber = CHL_RINGA
                     read(caSiteID(3:4),'(I2)') iSiteNumber
                  elseif(caSiteID(1:1) .eq. 'B') then
                     iRingNumber = CHL_RINGB
                     read(caSiteID(3:4),'(I2)') iSiteNumber
                  elseif(caSiteID(1:1) .eq. 'C') then
                     iRingNumber = CHL_RINGC
                     read(caSiteID(3:4),'(I2)') iSiteNumber
                  elseif(caSiteID(1:1) .eq. 'D') then
                     iRingNumber = CHL_RINGD
                     read(caSiteID(3:4),'(I2)') iSiteNumber
                  elseif(caSiteID(1:1) .eq. 'E') then
                     iRingNumber = CHL_RINGE
                     read(caSiteID(3:4),'(I2)') iSiteNumber
                  elseif(caSiteID(1:1) .eq. 'F') then
                     iRingNumber = CHL_RINGF
                     read(caSiteID(3:4),'(I2)') iSiteNumber
                  elseif(caSiteID(1:1) .eq. 'G') then
                     iRingNumber = CHL_RINGG
                     read(caSiteID(3:4),'(I2)') iSiteNumber
                  endif
                  
                  CHL_LOCALBURNUP(iRingNumber,iSiteNumber,iLayerID) =
     &               rDummy(1) * 0.001
                  goto 361
               else
                  goto 31
               endif
            else
               goto 31
            endif
     
     
     
31          print*, ' LIBR: Fatal error. Invalid local ',
     &        'burnup table format!'
            stop
32          continue          

            return
         end subroutine
      end module
