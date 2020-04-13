!                                                                       
!     Neutron Bank Manager Module, Revision 200413-1.
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
!     Last revision date 13/4/2020. Neutron bank modified to accomodate
!     the survival biasing technique.
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
         integer, allocatable :: ISGEN(:,:)  ! STARTING NEUTRON GENERATION (FOR SCO)
         real   , allocatable :: B_DIR(:,:)  ! BORN DIRECTION                                    
         real   , allocatable :: C_WGT(:,:)  ! NEUTRON WEIGHT.
         real   , allocatable :: L_WGT(:,:)  ! LAST NEUTRON WEIGHT
         real   , allocatable :: C_AWT(:,:)  ! ABSORBED WEIGHT FOR SURVIVAL BIASING.
         
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
         real   , allocatable :: FB_DIR(:,:)  ! BORN DIRECTION
         real   , allocatable :: FC_AWT(:,:)  ! ABSORBED WEIGHT FOR SURVIVAL BIASING.
         
       
         
         integer :: N_COUNT  =0  ! TOTAL NUMBER OF NEUTRONS IN THE BANK
         integer :: N_SOURCE =0 ! NUMBER OF SOURCE NEUTRONS
         integer :: N_GROUP   ! NUMBER OF NEUTRON GROUPS
         integer :: N_KEFFSRC ! NUMBER OF K-EFF SOURCE COUNT
         integer :: N_KEFFHIS ! NUMBER OF K-EFF HISTORY TRACK
         integer :: N_KEFFDLY
         integer, parameter :: BANK_SIZE = 6000000_8
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
            allocate( B_DIR(BANK_SIZE, 3) )
            allocate( ISGEN(BANK_SIZE, 1) )
            allocate( C_AWT(BANK_SIZE, 1) )
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
            allocate( FB_DIR(BANK_SIZE, 3) )
            allocate( FC_AWT(BANK_SIZE, 1) )
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
     &                                piG, prWgt, piGen)
 
            implicit none
            
            real   , intent(in) :: prX, prY, prZ, prU, prV, prW, prWgt
            integer, intent(in) :: piG, piGen
            
            integer :: iCell, iLayer
            if(N_COUNT .gt. BANK_SIZE) then
               print*, ' HGMC: Fatal error. Primary bank is full!'
               stop
               return
            endif
            
            call GetCurrentCell(prX, prY, prZ,
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
            B_DIR(N_COUNT, 1) = prU
            B_DIR(N_COUNT, 2) = prV
            B_DIR(N_COUNT, 3) = prW            
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
            C_AWT(N_COUNT, 1) = 0.0
            ISGEN(N_COUNT, 1) = piGen
            return
            
         end subroutine


         subroutine AddFissionNeutron(prX, prY, prZ,
     &                                prU, prV, prW,
     &                                piG, prWgt, piCell,
     &                                piLayer)
 
            implicit none
            
            real   , intent(in) :: prX, prY, prZ, prU, prV, prW, prWgt
            integer, intent(in) :: piG, piCell, piLayer
            
            integer :: iCell, iLayer
            if(N_KEFFSRC .gt. BANK_SIZE) then
               print*, ' HGMC: Fatal error. Fission bank is full!'
               stop
               return
            endif
            
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
            FB_DIR(N_KEFFSRC, 1) = prU
            FB_DIR(N_KEFFSRC, 2) = prV
            FB_DIR(N_KEFFSRC, 3) = prW
            FC_CEL(N_KEFFSRC, 1) = piCell
            FC_CEL(N_KEFFSRC, 2) = piLayer
            FL_CEL(N_KEFFSRC, 1) = piCell
            FL_CEL(N_KEFFSRC, 2) = piLayer
            FB_CEL(N_KEFFSRC, 1) = piCell
            FB_CEL(N_KEFFSRC, 2) = piLayer
            FC_GRP(N_KEFFSRC, 1) = piG
            FL_GRP(N_KEFFSRC, 1) = piG
            FNCHLD(N_KEFFSRC, 1) = 0         ! 0 - NO CHILD NEUTRONS
            FPARNT(N_KEFFSRC, 1) = -1        ! -1 NO PARENT
            FISTAT(N_KEFFSRC, 1) = N_WAIT         ! 0 - WAIT_QUEUE
            FC_WGT(N_KEFFSRC, 1) = prWgt  
            FL_WGT(N_KEFFSRC, 1) = prWgt
            FC_AWT(N_KEFFSRC, 1) = 0.0
            return
            
         end subroutine

         integer function PrimaryCount()
            PrimaryCount = N_SOURCE
            return
         end function
         
         integer function TotalNeutronCount()
            TotalNeutronCount = N_SOURCE
            return
         end function
         

        
      end module
