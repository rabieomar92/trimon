!                                                                       
!     Neutron Bank Manager Module, Revision 180915-1.
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

      module NeutronBankManager
         
         use CellGeometry
         use TXSReader
         use RandomNumberGenerator
        
         implicit none 
         
         real   , allocatable :: C_DIR(:,:)  ! CURRENT NEUTRON DIRECTION
         real   , allocatable :: L_DIR(:,:)  ! LAST NEUTRON DIRECTION
         real   , allocatable :: C_POS(:,:)  ! CURRENT NEUTRON POSITION
         real   , allocatable :: L_POS(:,:)  ! LAST NEUTRON POSITION
         real   , allocatable :: L_COL(:,:)  ! LAST COLLISION SITE
         real   , allocatable :: B_POS(:,:)  ! BORN POSITION
         integer, allocatable :: L_RXN(:,:)  ! LAST REACTION TYPE (0-ABSORPTION, 1-SCATTER, 2-FISSION)
         integer, allocatable :: C_CEL(:,:)  ! CURRENT REACTOR CORE CELL
         integer, allocatable :: L_CEL(:,:)  ! LAST REACTOR CORE CELL
         integer, allocatable :: B_CEL(:,:)  ! BORN CELL
         integer, allocatable :: C_GRP(:,:)  ! CURRENT NEUTRON ENERGY GROUP
         integer, allocatable :: L_GRP(:,:)  ! LAST NEUTRON ENERGY GROUP
         integer, allocatable :: NCHLD(:,:)  ! NUMBER OF SECONDARY CHILD NEUTRONS PRODUCED
         integer, allocatable :: PARNT(:,:)  ! PARENT NEUTRON ID, (-1 FOR NO PARENT)
         integer, allocatable :: ISTAT(:,:)  ! STATUS ((>4)-EMPTY, 0-WAIT_QUEUE,
                                             ! 1-KILLED, 2-ESCAPED, 3-RANDOM STRIDE OVERFLOW)
         real   , allocatable :: C_WGT(:,:)  ! NEUTRON WEIGHT.
         real   , allocatable :: L_WGT(:,:)  ! LAST NEUTRON WEIGHT
         
         ! SPECIAL FISSION BANK
         
         real   , allocatable :: FC_DIR(:,:)  ! CURRENT NEUTRON DIRECTION
         real   , allocatable :: FL_DIR(:,:)  ! LAST NEUTRON DIRECTION
         real   , allocatable :: FC_POS(:,:)  ! CURRENT NEUTRON POSITION
         real   , allocatable :: FL_POS(:,:)  ! LAST NEUTRON POSITION
         real   , allocatable :: FL_COL(:,:)  ! LAST COLLISION SITE
         real   , allocatable :: FB_POS(:,:)  ! BORN POSITION
         integer, allocatable :: FL_RXN(:,:)  ! LAST REACTION TYPE (0-ABSORPTION, 1-SCATTER, 2-FISSION)
         integer, allocatable :: FC_CEL(:,:)  ! CURRENT REACTOR CORE CELL
         integer, allocatable :: FL_CEL(:,:)  ! LAST REACTOR CORE CELL
         integer, allocatable :: FB_CEL(:,:)  ! BORN CELL
         integer, allocatable :: FC_GRP(:,:)  ! CURRENT NEUTRON ENERGY GROUP
         integer, allocatable :: FL_GRP(:,:)  ! LAST NEUTRON ENERGY GROUP
         integer, allocatable :: FNCHLD(:,:)  ! NUMBER OF SECONDARY CHILD NEUTRONS PRODUCED
         integer, allocatable :: FPARNT(:,:)  ! PARENT NEUTRON ID, (-1 FOR NO PARENT)
         integer, allocatable :: FISTAT(:,:)  ! STATUS ((>4)-EMPTY, 0-WAIT_QUEUE,
                                              ! 1-KILLED, 2-ESCAPED, 3-RANDOM STRIDE OVERFLOW)
         real   , allocatable :: FC_WGT(:,:)  ! NEUTRON WEIGHT.
         real   , allocatable :: FL_WGT(:,:)  ! LAST NEUTRON WEIGHT
         integer :: N_COUNT   ! TOTAL NUMBER OF NEUTRONS IN THE BANK
         integer :: N_SOURCE  ! NUMBER OF SOURCE NEUTRONS
         integer :: N_GROUP   ! NUMBER OF NEUTRON GROUPS
         integer :: N_KEFFSRC ! NUMBER OF K-EFF SOURCE COUNT
         integer :: N_KEFFHIS ! NUMBER OF K-EFF HISTORY TRACK
         integer, parameter :: BANK_SIZE = 5000000
         ! Neutron reaction parameters
         integer, parameter :: RXN_FISSION = 3
         integer, parameter :: RXN_CAPTURE = 1
         integer, parameter :: RXN_SCATTER = 2
         ! Neutron status parameters
         integer, parameter :: N_WAIT = 0 ! WAITING/PENDING HISTORY SIMULATION
         integer, parameter :: N_DEAD = 1 ! DEAD NEUTRON.
         integer, parameter :: N_LOST = 2 ! NEUTRON LOST OUTSIDE THE REACTOR.
         integer, parameter :: N_FLOW = 3 ! NEUTRON KILLED DUE TO RANDOM NUMBER STRIDE OVERFLOW.
         integer, parameter :: N_FISS = 4 ! NEUTRON CAPTURED BY FUEL NUCLEUS.
        
      contains
   
         subroutine InitNeutronBank()
            ! OBTAIN THE NEUTRON GROUP COUNT FROM TXSREADER. DURING INITIALISATION OF
            ! OF TXSREADER VIA InitTXSReader(pNF), THE NUMBER OF NEUTRON GROUP WAS READ
            ! AND IT IS STORED IN TXS_GRO
            N_GROUP = TXS_GRP        
            N_COUNT = 0
            N_SOURCE = 0
            N_KEFFSRC = 0
            allocate( C_DIR(BANK_SIZE, 3) )
            allocate( L_DIR(BANK_SIZE, 3) )
            allocate( C_POS(BANK_SIZE, 3) )
            allocate( L_POS(BANK_SIZE, 3) )
            allocate( L_COL(BANK_SIZE, 3) )
            allocate( B_POS(BANK_SIZE, 3) )
            allocate( L_RXN(BANK_SIZE, 1) )
            allocate( C_CEL(BANK_SIZE, 2) )
            allocate( L_CEL(BANK_SIZE, 2) )
            allocate( B_CEL(BANK_SIZE, 2) )
            allocate( C_GRP(BANK_SIZE, 1) )
            allocate( L_GRP(BANK_SIZE, 1) )
            allocate( NCHLD(BANK_SIZE, 1) )
            allocate( PARNT(BANK_SIZE, 1) )
            allocate( ISTAT(BANK_SIZE, 1) )
            allocate( C_WGT(BANK_SIZE, 1) )
            allocate( L_WGT(BANK_SIZE, 1) )
            ! INITIALIZE FISSION BANK
            allocate( FC_DIR(BANK_SIZE, 3) )
            allocate( FL_DIR(BANK_SIZE, 3) )
            allocate( FC_POS(BANK_SIZE, 3) )
            allocate( FL_POS(BANK_SIZE, 3) )
            allocate( FL_COL(BANK_SIZE, 3) )
            allocate( FB_POS(BANK_SIZE, 3) )
            allocate( FL_RXN(BANK_SIZE, 1) )
            allocate( FC_CEL(BANK_SIZE, 2) )
            allocate( FL_CEL(BANK_SIZE, 2) )
            allocate( FB_CEL(BANK_SIZE, 2) )
            allocate( FC_GRP(BANK_SIZE, 1) )
            allocate( FL_GRP(BANK_SIZE, 1) )
            allocate( FNCHLD(BANK_SIZE, 1) )
            allocate( FPARNT(BANK_SIZE, 1) )
            allocate( FISTAT(BANK_SIZE, 1) )
            allocate( FC_WGT(BANK_SIZE, 1) )
            allocate( FL_WGT(BANK_SIZE, 1) )
            return
            
         end subroutine
         
         
         subroutine ClearNeutronBank()
            N_COUNT  = 0
            N_SOURCE  = 0
            return
         end subroutine
         
         subroutine ClearFissionBank()
            N_KEFFSRC = 0
            return
         end subroutine
         
         subroutine SetKeffHistoryCount(pnHistoryCount)
            integer, intent(in) :: pnHistoryCount
            if(pnHistoryCount .gt. 0) N_KEFFHIS = pnHistoryCount
            return
         end subroutine
         
         
         subroutine AddPrimaryNeutron(prX, prY, prZ, 
     &                                prU, prV, prW,
     &                                piG, prWgt)
 
            implicit none
            
            real   , intent(in) :: prX, prY, prZ, prU, prV, prW, prWgt
            integer, intent(in) :: piG
            
            integer :: iCell, iLayer
            if(N_COUNT .gt. BANK_SIZE) then
               print*, ' Fatal Error: Neutron bank is full!'
               return
            endif
            
            call GetTargetCell(prX, prY, prZ,
     &                         prU, prV, prW,
     &                         iCell, iLayer)
     
            N_COUNT = N_COUNT + 1
            N_SOURCE = N_SOURCE + 1
            
            C_DIR(N_COUNT, 1) = prU
            C_DIR(N_COUNT, 2) = prV
            C_DIR(N_COUNT, 3) = prW
            L_DIR(N_COUNT, 1) = prU
            L_DIR(N_COUNT, 2) = prV
            L_DIR(N_COUNT, 3) = prW
            C_POS(N_COUNT, 1) = prX
            C_POS(N_COUNT, 2) = prY
            C_POS(N_COUNT, 3) = prZ
            L_POS(N_COUNT, 1) = prX
            L_POS(N_COUNT, 2) = prY
            L_POS(N_COUNT, 3) = prZ
            B_POS(N_COUNT, 1) = prX
            B_POS(N_COUNT, 2) = prY
            B_POS(N_COUNT, 3) = prZ
            C_CEL(N_COUNT, 1) = iCell
            C_CEL(N_COUNT, 2) = iLayer
            L_CEL(N_COUNT, 1) = iCell
            L_CEL(N_COUNT, 2) = iLayer
            B_CEL(N_COUNT, 1) = iCell
            B_CEL(N_COUNT, 2) = iLayer
            C_GRP(N_COUNT, 1) = piG
            L_GRP(N_COUNT, 1) = piG
            NCHLD(N_COUNT, 1) = 0         ! 0 - NO CHILD NEUTRONS
            PARNT(N_COUNT, 1) = -1        ! -1 NO PARENT
            ISTAT(N_COUNT, 1) = 0         ! 0 - WAIT_QUEUE
            C_WGT(N_COUNT, 1) = prWgt
            L_WGT(N_COUNT, 1) = prWgt
            return
            
         end subroutine
         
         subroutine AddSecondaryNeutron(prX, prY, prZ, 
     &                                  prU, prV, prW, 
     &                                  piG, piParent, piLastRxn, prWgt)

            real   , intent(in) :: prX, prY, prZ, prU, prV, prW, prWgt
            integer, intent(in) :: piG, piParent, piLastRxn
            
            integer :: iCell, iLayer
            
            if(N_COUNT .gt. BANK_SIZE) then
               print*, ' Fatal Error: Neutron bank is full!'
               return
            endif
            
            call GetTargetCell(prX, prY, prZ,
     &                         prU, prV, prW,
     &                         iCell, iLayer)
     
            N_COUNT = N_COUNT + 1
            
            C_DIR(N_COUNT, 1) = prU
            C_DIR(N_COUNT, 2) = prV
            C_DIR(N_COUNT, 3) = prW
            L_DIR(N_COUNT, 1) = prU
            L_DIR(N_COUNT, 2) = prV
            L_DIR(N_COUNT, 3) = prW
            C_POS(N_COUNT, 1) = prX
            C_POS(N_COUNT, 2) = prY
            C_POS(N_COUNT, 3) = prZ
            L_POS(N_COUNT, 1) = prX
            L_POS(N_COUNT, 2) = prY
            L_POS(N_COUNT, 3) = prZ
            B_POS(N_COUNT, 1) = prX
            B_POS(N_COUNT, 2) = prY
            B_POS(N_COUNT, 3) = prZ
            C_CEL(N_COUNT, 1) = iCell
            C_CEL(N_COUNT, 2) = iLayer
            B_CEL(N_COUNT, 1) = iCell
            B_CEL(N_COUNT, 2) = iLayer
            C_GRP(N_COUNT, 1) = piG
            L_GRP(N_COUNT, 1) = piG
            NCHLD(N_COUNT, 1) = 0         ! 0 - NO CHILD NEUTRONS
            PARNT(N_COUNT, 1) = piParent  
            ISTAT(N_COUNT, 1) = 0         ! 0 - WAIT_QUEUE     
            L_COL(N_COUNT, 1) = prX
            L_COL(N_COUNT, 2) = prY
            L_COL(N_COUNT, 3) = prZ
            L_RXN(N_COUNT, 1) = piLastRxn
            C_WGT(N_COUNT, 1) = prWgt
            L_WGT(N_COUNT, 1) = prWgt
            return
            
         end subroutine

         subroutine AddFissionNeutron(prX, prY, prZ,
     &                                prU, prV, prW,
     &                                piG, prWgt)
 
            implicit none
            
            real   , intent(in) :: prX, prY, prZ, prU, prV, prW, prWgt
            integer, intent(in) :: piG
            
            integer :: iCell, iLayer
            if(N_COUNT .gt. BANK_SIZE) then
               print*, ' Fatal Error: Neutron bank is full!'
               return
            endif
            
            call GetTargetCell(prX, prY, prZ,
     &                         prU, prV, prW,
     &                         iCell, iLayer)
     
            N_KEFFSRC = N_KEFFSRC + 1
            
            FC_DIR(N_KEFFSRC, 1) = prU
            FC_DIR(N_KEFFSRC, 2) = prV
            FC_DIR(N_KEFFSRC, 3) = prW
            FL_DIR(N_KEFFSRC, 1) = prU
            FL_DIR(N_KEFFSRC, 2) = prV
            FL_DIR(N_KEFFSRC, 3) = prW
            FC_POS(N_KEFFSRC, 1) = prX
            FC_POS(N_KEFFSRC, 2) = prY
            FC_POS(N_KEFFSRC, 3) = prZ
            FL_POS(N_KEFFSRC, 1) = prX
            FL_POS(N_KEFFSRC, 2) = prY
            FL_POS(N_KEFFSRC, 3) = prZ
            FB_POS(N_KEFFSRC, 1) = prX
            FB_POS(N_KEFFSRC, 2) = prY
            FB_POS(N_KEFFSRC, 3) = prZ
            FC_CEL(N_KEFFSRC, 1) = iCell
            FC_CEL(N_KEFFSRC, 2) = iLayer
            FL_CEL(N_KEFFSRC, 1) = iCell
            FL_CEL(N_KEFFSRC, 2) = iLayer
            FB_CEL(N_KEFFSRC, 1) = iCell
            FB_CEL(N_KEFFSRC, 2) = iLayer
            FC_GRP(N_KEFFSRC, 1) = piG
            FL_GRP(N_KEFFSRC, 1) = piG
            FNCHLD(N_KEFFSRC, 1) = 0         ! 0 - NO CHILD NEUTRONS
            FPARNT(N_KEFFSRC, 1) = -1        ! -1 NO PARENT
            FISTAT(N_KEFFSRC, 1) = N_WAIT         ! 0 - WAIT_QUEUE
            FC_WGT(N_KEFFSRC, 1) = prWgt  
            FL_WGT(N_KEFFSRC, 1) = prWgt
            return
            
         end subroutine
         

         subroutine KillNeutron(piNID, piReason, piStat)
            integer, intent(in) :: piNID, piReason
            integer, intent(out) :: piStat
            
            if((piNID .gt. 0) .and. (piNID .le. N_COUNT)) then
               if(ISTAT(piNID, 1) .ne. N_WAIT) then
                  ! Case killing a killed neutron.
                  piStat = 1
                  return
               else
                  ISTAT(piNID, 1) = piReason
                  piStat = 0
                  return
               endif
            else
               ! Case invalid NID.
               piStat = 2
               return
            endif
         end subroutine
         
         integer function PrimaryCount()
            PrimaryCount = N_SOURCE
            return
         end function
         
         integer function SecondaryCount()
            SecondaryCount = N_COUNT - N_SOURCE
            return
         end function
         
         integer function TotalNeutronCount()
            TotalNeutronCount = N_COUNT
            return
         end function
         
         subroutine SetNeutronDirection(piNID, prU, prV, prW, piStat)
            integer, intent(in) :: piNID
            real   , intent(in) :: prU, prV, prW
            integer, intent(out) :: piStat
            if((piNID .gt. 0) .and. (piNID .le. N_COUNT)) then
               if(ISTAT(piNID, 1) .ne. N_WAIT) then
                  ! Case if the neutron is already killed.
                  piStat = 1
                  return
               else
                  L_DIR(piNID, 1) = C_DIR(piNID, 1)
                  L_DIR(piNID, 2) = C_DIR(piNID, 2)
                  L_DIR(piNID, 3) = C_DIR(piNID, 3)      
                  C_DIR(piNID, 1) = prU
                  C_DIR(piNID, 2) = prV
                  C_DIR(piNID, 3) = prW
                  piStat = 0
                  return
               endif
            else
               ! Case if invalid NID is specified.
               piStat = 2
               return
            endif
         end subroutine
         
         subroutine SetNeutronPosition(piNID, prX, prY, prZ, piStat)
            integer, intent(in) :: piNID
            real   , intent(in) :: prX, prY, prZ
            integer, intent(out) :: piStat  
            
            integer :: iCell, iLayer
            

            
            if((piNID .gt. 0) .and. (piNID .le. N_COUNT)) then
               if(ISTAT(piNID, 1) .ne. N_WAIT) then
                  ! Case if the neutron is not online.
                  piStat = 1
                  return
               else
               
                  call GetTargetCell(prX, prY, prZ, 
     &                               C_DIR(piNID,1),
     &                               C_DIR(piNID,1),
     &                               C_DIR(piNID,1),
     &                               iCell, iLayer)
                  if((C_DIR(piNID,1) .gt. TXS_CEL) .or. 
     &               (C_DIR(piNID,1) .lt. -2     )) then
                     C_DIR(piNID,1) = -2
                  endif
                  if((C_DIR(piNID,2) .gt. TXS_LAY) .or. 
     &               (C_DIR(piNID,2) .lt. -2     )) then
                     C_DIR(piNID,2) = -2
                  endif
                  L_POS(piNID, 1) = C_POS(piNID, 1)
                  L_POS(piNID, 2) = C_POS(piNID, 2)
                  L_POS(piNID, 3) = C_POS(piNID, 3)      
                  C_POS(piNID, 1) = prX
                  C_POS(piNID, 2) = prY
                  C_POS(piNID, 3) = prZ
                  L_CEL(piNID, 1) = C_CEL(piNID, 1)
                  L_CEL(piNID, 2) = C_CEL(piNID, 2)
                  C_CEL(piNID, 1) = iCell
                  C_CEL(piNID, 2) = iLayer
                  
                  piStat = 0
                  return
               endif
            else
               ! Case if invalid NID is specified.
               piStat = 2
               return
            endif     
               
         end subroutine
         
         subroutine SetNewEnergyGroup(piNID, piG, piStat)
            integer, intent(in) :: piNID
            integer, intent(in) :: piG
            integer, intent(out) :: piStat
            
            if((piNID .gt. 0) .and. (piNID .le. N_COUNT)) then
               if(ISTAT(piNID, 1) .ne. N_WAIT) then
                  ! Case if the neutron is not online.
                  piStat = 1
                  return
               else
                  L_GRP(piNID, 1) = C_GRP(piNID, 1)
                  C_GRP(piNID, 1) = piG
                  piStat = 0
                  return
               endif
            else
               ! Case if invalid NID is specified.
               piStat = 2
               return
            endif                
         end subroutine
         
         subroutine SetLastCollisionSite(piNID, prX, prY, prZ,
     &                                   piRxn, piStat)
            integer, intent(in) :: piNID, piRxn
            real   , intent(in) :: prX, prY, prZ
            integer, intent(out) :: piStat
            
            if((piNID .gt. 0) .and. (piNID .le. N_COUNT)) then
               if(ISTAT(piNID, 1) .ne. N_WAIT) then
                  ! Case if the neutron is not online.
                  piStat = 1
                  return
               else
                  L_COL(piNID, 1) = prX
                  L_COL(piNID, 2) = prY
                  L_COL(piNID, 3) = prZ
                  L_RXN(piNID, 1) = piRxn
                  piStat = 0
                  return
               endif
            else
               ! Case if invalid NID is specified.
               piStat = 2
               return
            endif     
            
         end subroutine
        
      end module
