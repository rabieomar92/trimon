!                                                                       
!     Monte Carlo Tally Scoring Module, Revision 180915-1.
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
!
!     THIS MODULE CONTAINS THE NECESSARY METHODS THAT ENABLE THE CODE USER TO
!     DO MONTE CARLO TALLYING FOR OBTAINING NEUTRON FLUX/POWER
!
      
      module Tally
         
         use CellGeometry
         use TXSReader
         use NeutronBankManager
         
         implicit none
         
         ! HERE WE DEFINE CELL TRACK LENGTH TALLY
         real, allocatable :: TLY_TLVAL(:,:,:)
         real, allocatable :: TLY_TLAVG(:,:,:)
         real, allocatable :: TLY_TLSQA(:,:,:)
         real, allocatable :: TLY_TLRER(:,:,:)
         ! DEFINING TRACK LENGTH FISSION CYCLE AVERAGE VALUE

         real :: TLY_ABSEST = 0.0
         real :: TLY_COLEST = 0.0
         real :: TLY_TLEST  = 0.0
         ! THIS CORRESPONDS TO THE NUMBER OF TIMES ScoreTrackLengthTally IS
         ! INVOKED.
         real :: TLY_CYCLE = 0
         
      contains
         
!        ----------------------------------------------------------------------
!        THIS SUBROUTINE INITIALIZES TALLY TABLES BY SETTING THEIR DIMENSIONS,
!        AND INITIALIZING THEIR ELEMENT VALUES TO ZERO.
!        ----------------------------------------------------------------------
         subroutine InitTally()
            integer :: i, j, k
            allocate(TLY_TLVAL(TXS_CEL,
     &                         TXS_LAY, TXS_GRP))
            allocate(TLY_TLAVG(TXS_CEL,
     &                         TXS_LAY, TXS_GRP))
            allocate(TLY_TLSQA(TXS_CEL,
     &                         TXS_LAY, TXS_GRP))
            allocate(TLY_TLRER(TXS_CEL,
     &               TXS_LAY, TXS_GRP))
            
       
            
            do i=1, TXS_CEL, 1
               do j=1, TXS_LAY, 1
                  do k=1, TXS_GRP, 1
                  
                     ! FOR TRACK LENGTH TALLY
                     TLY_TLVAL(i, j, k) = 0.0
                     TLY_TLAVG(i, j, k) = 0.0
                     TLY_TLSQA(i, j, k) = 0.0
                     TLY_TLRER(i, j, k) = 0.0

                     TLY_ABSEST = 0.0
                     TLY_TLEST = 0.0                     
                  enddo
               enddo
            enddo            
            return
            
         end subroutine

!        ----------------------------------------------------------------------     
!        THIS  SUBROUTINE  RESETS THE  TALLY  TABLES  BY  SETTING THEIR ELEMENT 
!        VALUES TO ZERO. THE PURPOSE IS TO CLEAR THE UNUSED TALLY DATA FROM THE
!        PREVIOUS FISSION CYCLE CALCULATION.
!        ----------------------------------------------------------------------
         subroutine ResetTally()
         
            integer i, j, k 
            
            do i=1, TXS_CEL, 1
               do j=1, TXS_LAY, 1
                  do k=1, TXS_GRP, 1
                  
                     ! FOR TRACK LENGTH TALLY
                     TLY_TLVAL(i, j, k) = 0.0
                     
                     TLY_ABSEST = 0.0
                     TLY_TLEST = 0.0                     
                  enddo
               enddo
            enddo
        
            return
            
         end subroutine
         
!        ----------------------------------------------------------------------
!        FOR EACH NEUTRON COLLISIONS, NEUTRON FLUX TALLIES ARE RECORDED VIA 
!        THIS SUBROUTINE. THE NEUTRON FLUX IS ESTIMATED BY USING TRACK LENGTH
!        TALLYING METHOD.
!
!        GIVEN W AS NEUTRON WEIGHT, T AS TOTAL SUM OF THE DISTANCE TRAVELLED 
!        BY ALL NEUTRONS WITHIN CELL (I,K,L), V AS THE VOLUME OF THE CELL AND
!        ΣW IS THE TOTAL WEIGHT OF . THEN THE FLUX IS GIVEN BY
!        ----------------------------------------------------------------------
         subroutine ScoreTrackLengthTally(piTrackNID, prTotalWeight)
         
            integer, intent(in) :: piTrackNID
            real   , intent(in) :: prTotalWeight
            
            real    :: rX, rY, rZ      ! CURRENT NEUTRON POSITION
            real    :: rX0, rY0, rZ0   ! LAST NEUTRON POSITION
            real    :: rTrackLength
            integer :: iCell, iLayer   ! CURRENT CELL INDICES
            integer :: iCell0, iLayer0 ! LAST CELL INDICES
            integer :: iGroup          ! CURRENT NEUTRON ENERGY GROUP
            real    :: rWeight
            real    :: rNuSigFis
            iGroup  = L_GRP(piTrackNID, 1)
            rX      = C_POS(piTrackNID, 1)
            rY      = C_POS(piTrackNID, 2)
            rZ      = C_POS(piTrackNID, 3)
            rX0     = L_POS(piTrackNID, 1)
            rY0     = L_POS(piTrackNID, 2)
            rZ0     = L_POS(piTrackNID, 3)
            rWeight = C_WGT(piTrackNID, 1)
            iCell   = C_CEL(piTrackNID, 1)
            iLayer  = C_CEL(piTrackNID, 2)
            iCell0  = L_CEL(piTrackNID, 1)
            iLayer0 = L_CEL(piTrackNID, 2)
            if((iCell .ne. iCell0) .or. (iLayer .ne. iLayer0)) then
               return
            endif
            
            if((iCell .gt. 0) .and. (iLayer .gt. 0) .and.
     &         (iCell .le. TXS_CEL) .and.
     &         (iLayer .le. TXS_LAY) .and.
     &         (iGroup .le. TXS_GRP) .and. (iGroup .gt. 0)) then
     
               rTrackLength = sqrt( ( rX - rX0 )**2 +
     &                              ( rY - rY0 )**2 + 
     &                              ( rZ - rZ0 )**2   ) 
               
               ! APPLYING NEUTRON WEIGHT
               rTrackLength = rWeight * rTrackLength / prTotalWeight
               
               ! ADDING SUM TO THE PREVIOUS VALUES.
               TLY_TLVAL(iCell, iLayer, iGroup) = 
     &            TLY_TLVAL(iCell, iLayer, iGroup) + rTrackLength
               TLY_TLEST = TLY_TLEST + rTrackLength 
               return
            endif
            
            return
            
         end subroutine
         
         
 
         subroutine ScoreKeffAbsEst(piTrackNID)
            integer, intent(in) :: piTrackNID
            
            real :: rWeight, rNuSigFis, rSigAbs
            integer :: iCell, iLayer, iGroup
            
            rWeight = C_WGT(piTrackNID, 1)
            iCell   = C_CEL(piTrackNID, 1)
            iLayer  = C_CEL(piTrackNID, 2)
            iGroup  = C_GRP(piTrackNID, 1)
            rNuSigFis = TXS_TABLE(iLayer, iCell, iGroup, 3)
            rSigAbs = TallyGetSigAbs(piTrackNID)
            TLY_ABSEST = TLY_ABSEST + nint(rNuSigFis*rWeight/rSigAbs)

            return
         end subroutine
         
         subroutine ScoreKeffColEst(piTrackNID)
            integer, intent(in) :: piTrackNID
            
            real :: rWeight, rSigTot, rNuSigFis
            integer :: iCell, iLayer, iGroup
            
            rWeight = C_WGT(piTrackNID, 1)
            iCell   = C_CEL(piTrackNID, 1)
            iLayer  = C_CEL(piTrackNID, 2)
            iGroup  = C_GRP(piTrackNID, 1)
            rNuSigFis = TXS_TABLE(iLayer, iCell, iGroup, 3)
            rSigTot = TallyGetSigTot(piTrackNID)
            TLY_COLEST = TLY_COLEST + rWeight * rNuSigFis / rSigTot
         end subroutine
         
         subroutine UpdateTallyAverages()
            
            integer :: i, j, k

            TLY_CYCLE = TLY_CYCLE + 1
            
          
            ! ESTIMATING THE FLUX DISTRIBUTION IN A CELL
            do i=1, TXS_CEL, 1
            do j=1, TXS_LAY, 1
            do k=1, TXS_GRP, 1
            
               if(TLY_CYCLE .eq. 1) then
               
                  TLY_TLAVG(i, j, k) = TLY_TLVAL(i, j, k)
                  TLY_TLSQA(i, j, k) = TLY_TLVAL(i, j, k)**2.0

               else
               
                  TLY_TLAVG(i, j, k) = (1.0 - 1.0 / 
     &               real(TLY_CYCLE)) * TLY_TLAVG(i, j, k) +
     &               TLY_TLVAL(i, j, k) / real(TLY_CYCLE)
                  TLY_TLSQA(i, j, k) = (1.0 - 1.0 / 
     &               real(TLY_CYCLE)) * TLY_TLSQA(i, j, k) +
     &               TLY_TLVAL(i, j, k)**2.0 / real(TLY_CYCLE)    
                  TLY_TLRER(i, j, k) = ( 1.0 / TLY_TLAVG(i, j, k) ) *
     &               sqrt((TLY_TLSQA(i, j, k)-TLY_TLAVG(i, j, k)**2) /
     &                      real(TLY_CYCLE-1) ) 

               endif
               

            enddo
            enddo
            enddo

            return
         end subroutine
         
         
         subroutine WriteFluxToFile(pcRunID)
            
         
            integer, parameter :: NFO = 36
            character(len=30)  :: cFileName = 'FLUX.OUT'
            integer :: i, iLayer, iGroup, iCell
            character(len=4) :: cSite
            character(len=2) :: cGrp
            character(len=30), intent(in) :: pcRunID
            open(unit=NFO, file=cFileName, status='UNKNOWN', err=99)
            goto 100
99          print*, ' Error: Could not write flux data file.'
            return
100         continue

            write(NFO,'(2A)') 'Run ID: ', trim(pcRunID)
            do iCell=1, TXS_CEL-1, 1
            do iLayer=1, TXS_LAY, 1
               call CellSite(iCell , cSite)
               write(cGrp,'(I2)') TXS_GRP*2
               write(NFO,'(A4,A1,I2.2,A2,' //
     &                    adjustr(cGrp) // 'ES12.5)')
     &            adjustl(cSite), '-', iLayer, '',
     &            (TLY_TLAVG(iCell, iLayer, i),
     &             TLY_TLRER(iCell, iLayer, i),
     &             i=1, TXS_GRP)
               enddo
            enddo
            close(NFO)
            print*, ' TRIMON: Finished writing flux to FLUX.OUT.'
         end subroutine
         
         subroutine WritePowerDistToFile(prNominalPower,pcRunID)
            
            real, intent(in) :: prNominalPower
            integer, parameter :: NFO = 39
            real, parameter :: rNu = 2.45
            real, parameter :: rC  = 3.15E-11
            real :: rPowerPerCell, rFCore, rNormFactor
            real :: rPowerPerCellError
            character(len=30)  :: cFileName = 'PDIST.OUT'
            integer :: i, iLayer, iGroup, iCell
            character(len=4) :: cSite
            character(len=2) :: cGrp
            character(len=30), intent(in) :: pcRunID
            open(unit=NFO, file=cFileName, status='UNKNOWN', err=99)
            write(NFO,'(2A)') 'Run ID: ', trim(pcRunID)
            goto 100
99          print*, ' Warning: Could not write power',
     &              ' distribution data file.'
            return
100         continue

!           HERE WE CALCULATE THE NORMALISATION FACTOR BY
!           INTEGRATING ALL TRACK LENGTH TALLIES OVER THE 
!           WHOLE CORE VOLUME.
            rFCore = 0.0
            do iCell=1, TXS_CEL, 1
            do iLayer=1, TXS_LAY, 1
            do iGroup=1, TXS_GRP, 1
            
               rFCore = rFCore + 
     &               TLY_TLAVG(iCell, iLayer, iGroup) *
     &               TXS_TABLE(iLayer, iCell, iGroup, 5)
     
            enddo
            enddo
            enddo
            
            rNormFactor = prNominalPower / rFCore
            

            do iCell=1, TXS_CEL-1, 1
            do iLayer=1, TXS_LAY, 1
            
               ! OBTAIN THE CELL SITE LABEL.
               call CellSite(iCell , cSite)
               
               rPowerPerCell = 0.0
               rPowerPerCellError = 0.0
               do iGroup=1, TXS_GRP, 1
                  ! HERE WE ACCUMULATE ALL POWER CONTRIBUTED BY
                  ! ALL NEUTRON ENERGY GROUPS.
                  rPowerPerCell = rPowerPerCell +
     &            TLY_TLAVG(iCell, iLayer, iGroup) *
     &            TXS_TABLE(iLayer, iCell, iGroup, 5)
                  
                  ! ALSO, WE ACCUMULATE THE ERROR.
                  rPowerPerCellError = rPowerPerCellError +
     &            TLY_TLRER(iCell, iLayer, iGroup) *
     &            TLY_TLAVG(iCell, iLayer, iGroup) *
     &            TXS_TABLE(iLayer, iCell, iGroup, 5)
               enddo
               
               rPowerPerCell = rPowerPerCell * rNormFactor
               rPowerPerCellError = rPowerPerCellError * rNormFactor
               
               write(NFO,'(A4,A1,I2.2,A2,2ES12.5)')
     &            adjustl(cSite), '-', iLayer, '',
     &            rPowerperCell, rPowerPerCellError
               enddo
            enddo
            close(NFO)
            print*, ' TRIMON: Finished writing power',
     &              ' distribution to PDIST.OUT.'
         end subroutine
         
         subroutine WritePowerToFile(prNominalPower,pcRunID)
            real, intent(in) :: prNominalPower
            character(len=30), intent(in) :: pcRunID
            ! Define parameter(s) here...
            integer, parameter :: NFP = 13
            character(len=30)  :: cFileName = 'PPE.OUT'
            integer :: i, iLayer, iGroup, iCell
            character(len=4) :: cSite
            character(len=2) :: cGrp
            real :: rPower, rPowerError
            real, parameter :: rNu = 2.45
            real, parameter :: rC  = 3.15E-11
            real :: rFCore      = 0.0
            real :: rNormFactor = 0.0
            open(unit=NFP, file=cFileName, status='UNKNOWN', err=99)

            write(NFP,'(2A)') 'Run ID: ', trim(pcRunID)          
            do iCell=1, TXS_CEL, 1
            do iLayer=1, TXS_LAY, 1
            do iGroup=1, TXS_GRP, 1
            
               rFCore = rFCore + 
     &               TLY_TLAVG(iCell, iLayer, iGroup) *
     &               TXS_TABLE(iLayer, iCell, iGroup, 5)
     
            enddo
            enddo
            enddo
            
            rNormFactor = prNominalPower / ( ( rC / rNu ) * rFCore )
            
            do iCell=1, TXS_CEL-1, 1

               call CellSite(iCell , cSite)
               
               rPower = 0.0
               rPowerError = 0.0
               
               do iLayer=1, TXS_LAY, 1
               do iGroup=1, TXS_GRP, 1
               
                  rPower = rPower + 
     &               TLY_TLAVG(iCell, iLayer, iGroup) *
     &               TXS_TABLE(iLayer, iCell, iGroup, 5)
                  rPowerError = rPowerError +
     &               TLY_TLRER(iCell, iLayer, iGroup)**2 *
     &               TLY_TLAVG(iCell, iLayer, iGroup)**2 *
     &               TXS_TABLE(iLayer, iCell, iGroup, 5)
                  
               enddo
               enddo
               
               rPower = rPower * ( rC / rNu ) * rNormFactor
               rPowerError =  ( rC / rNu ) * rNormFactor * 
     &                        sqrt(rPowerError)
               
               write(NFP,'(A4,A2,2ES12.5)')
     &            adjustl(cSite), '', 
     &            rPower, rPowerError

            
            enddo
            
            close(NFP)
            print*, ' TRIMON: Finished writing power-per',
     &              '-element to PPE.OUT.'
            
99          return
         end subroutine

        
         real function TallyGetSigAbs(piTrackNID)
            integer, intent(in) :: piTrackNID
            real :: rSigAbs
            integer :: nGroup
            integer :: iCell, iLayer

            ! OBTAIN CURRENT NEUTRON ENERGY GROUP AND CURRENT REACTOR CORE CELL
            nGroup = C_GRP(piTrackNID, 1)
            iCell  = C_CEL(piTrackNID, 1)
            iLayer  = C_CEL(piTrackNID, 2)
            ! IF THE NEUTRON IS CURRENTLY IN THE GRAPHITE REFLECTOR THEN SELECT 
            ! THE UNIQUE CELL ID FOR GRAPHITE REFLECTOR.
            if((iLayer .eq. -1) .and. (iCell .eq. -1)) then
               rSigAbs = TXS_TABLE(iLayer, REFLECTOR_CELLID, nGroup, 2)
            else
               rSigAbs = TXS_TABLE(iLayer, iCell, nGroup, 2)            
            endif

            ! IN CASE IF THE NEUTRON CROSS SECTION READ IS LESS THAN ZERO, WE
            ! CORRECT IT BY SETTING IT TO ZERO.
            if(rSigAbs .lt. 0) rSigAbs = 0.0
            TallyGetSigAbs = rSigAbs
            
            return
            
         end function
         
!        ----------------------------------------------------------------------
!        THIS FUNCTION CALCULATES THE TOTAL NEUTRON CROSS SECTION.
!        ----------------------------------------------------------------------
         real function TallyGetSigTot(piTrackNID)
            integer, intent(in) :: piTrackNID
            
            real    :: rSigTot
            integer :: nGroup
            integer :: iCell, iLayer
            integer :: i 
            real    :: rSigTemp
            
            ! WE OBTAIN THE CURRENT CELL ID  AND LAYER ID.  PLUS  WE OBTAIN THE 
            ! CURRENT NEUTRON ENERGY GROUP.
            nGroup = C_GRP(piTrackNID, 1)    
            iCell  = C_CEL(piTrackNID, 1)
            iLayer = C_CEL(piTrackNID, 2)
            rSigTot = 0.0
            
            ! ADD THE ABSORPTION CROSS SECTION (COLUMN-2) THIS INCLUDES FISSION.
            if((iCell .gt. TXS_CEL) .or. (iCell .lt. -1)) then
               TallyGetSigTot = 0.0
               return
            endif
            if((iLayer .gt. TXS_LAY) .or. (iLayer .lt. -1)) then
               TallyGetSigTot = 0.0
               return
            endif
            if((nGroup .gt. TXS_GRP) .or. (nGroup .lt. 1)) then
               TallyGetSigTot = 0.0
               return
            endif
            rSigTemp = TXS_TABLE(iLayer, iCell, nGroup, 2)
 
            if(rSigTemp .lt. 0) rSigTemp = 0.0
            rSigTot = rSigTot + rSigTemp
 
            ! ADD THE GROUP SCATTERING CROSS SECTIONS.
            do i=1, TXS_GRP, 1
               rSigTemp = TXS_TABLE(iLayer, iCell, nGroup, i+5)
               if(rSigTemp .lt. 0) rSigTemp = 0.0 
               rSigTot = rSigTot + rSigTemp
            enddo           
            TallyGetSigTot = rSigTot
  
            return
         end function
       
      end module
