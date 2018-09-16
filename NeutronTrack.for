!                                                                       
!     Neutron Tracking Module, Revision 180915-1.
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
!     This module defines the necessary methods of basic neutron tracking  pur-
!     pose. This includes the procedure on  choosing reaction  types,  changing 
!     neutron flight path, predicting distance to the next collision  and more.
!     It should be note that this code was written will full care. It is impor-
!     tant  to  emphasize  that  this code has to be genuine and sophisticated.
!     Otherwise,  the reactor  calculation  may leads to  error and may be mis-
!     leading.  Misleading  reactor  calculations could  even  lead  to reactor 
!     accidents!
!
!     Good luck!

      module NeutronTrack
      
         use RandomNumberGenerator
         use CellGeometry
         use NeutronBankManager
         use TXSReader
         use Tally
         
         implicit none
         
         ! CURRENT  NEUTRON ID  FOR  TRACKING  PURPOSE.   IF THE TRACK_ID VALUE
         ! CHANGED, ALL  METHODS IN THIS  MODULE  WILL REFER THE NEW VALUE.
         integer :: TRACK_NID
      contains
         
!        ----------------------------------------------------------------------
!        THIS SUBROUTINE ALLOWS THE PROGRAMMER TO CHANGE THE CURRENT NEUTRON ID
!        TO  BE  TRACKED  AND  IT  ALSO  VALIDATES  THE   NEUTRON ID,  CHECKING 
!        CHECKING  WHETHER  THE NEUTRON ID IS LISTED IN THE NEUTRON BANK.
!        ----------------------------------------------------------------------
         subroutine BeginTracking(piNID)
            integer, intent(in) :: piNID
            integer :: i
            if((piNID .gt. 0) .and. 
     &         (piNID .le. TotalNeutronCount())) then
               if(ISTAT(piNID, 1) .ne. N_WAIT) then
!                  print*, ' Warning: Cannot track neutron ', piNID, '.'
                  return
               else
                  TRACK_NID = piNID
               endif
            endif

         end subroutine
         
!        ----------------------------------------------------------------------
!        THIS FUNCTION SAMPLES THE DISTANCE TO THE NEXT COLLISION SITE. THE COL
!        LISION DISTANCE IS SAMPLED BASED ON THE EXPONENTIAL DISTRIBUTION.
!        ---------------------------------------------------------------------- 
         real function DistanceToNextCollision()
            real :: rEpsilon, rSigTot
            
            ! FIRST WE OBTAIN THE TOTAL NEUTRON CROSS SECTION.
            rSigTot = GetSigTot()
         
            ! WE CALCULATE THE DISTANCE TO NEXT COLLISION USING
            !                          ɛ
            !              D = - log ----
            !                         Σt
            !
            rEpsilon = Rnd(int(TRACK_NID,8))
            DistanceToNextCollision = - log(rEpsilon) / rSigTot
            
            return
         end function
         
!        ----------------------------------------------------------------------
!        THIS FUNCTION RETURNS THE DISTANCE  TO THE NEAREST  CELL BOUNDARY FROM 
!        THE CURRENT POSITION  OF THE CURRENT  NEUTRON  BEING TRACKED.  IT ALSO 
!        CORRECTS  THE  FLOATING  POINT  PRECISION ISSUES  THAT  CAUSES NEUTRON 
!        TRANSPORT TO BE TRAPPED BETWEEN CELL BOUNDARIES. 
!        ----------------------------------------------------------------------
         real function DistanceToNearestBoundary()
            real :: rX, rY, rZ
            real :: rU, rV, rW
            real :: rDistance
            
            ! HERE WE OBTAIN THE CURRENT POSITION AND DIRECTION OF THE NEUTRON.
            rX = C_POS(TRACK_NID, 1)
            rY = C_POS(TRACK_NID, 2)
            rZ = C_POS(TRACK_NID, 3)
            rU = C_DIR(TRACK_NID, 1)
            rV = C_DIR(TRACK_NID, 2)
            rW = C_DIR(TRACK_NID, 3)
            
            ! WE CALCULATE THE DISTANCE TO BOUNDARY USING THE FUNCTION PROVIDED
            ! BY THE GEOMETRY MODULE.
            rDistance = DistanceToBoundary(rX, rY, rZ, rU, rV, rW)
            
            ! HERE  WE  CORRECT  FLOATING  POINT  ISSUES  THAT  CAUSES  NEUTRON 
            ! BOUNDARY  TRAP.   
            if(rDistance .lt. 1E-05) then
               rDistance = rDistance + 1E-05
            endif
            
            DistanceToNearestBoundary = rDistance
            
            return
            
         end function
         
!        ----------------------------------------------------------------------
!        THIS FUNCTION CALCULATES THE TOTAL NEUTRON CROSS SECTION.
!        ----------------------------------------------------------------------
         real function GetSigTot()
            real    :: rSigTot
            integer :: nGroup
            integer :: iCell, iLayer
            integer :: i 
            real    :: rSigTemp
            
            ! WE OBTAIN THE CURRENT CELL ID  AND LAYER ID.  PLUS  WE OBTAIN THE 
            ! CURRENT NEUTRON ENERGY GROUP.
            nGroup = C_GRP(TRACK_NID, 1)    
            iCell  = C_CEL(TRACK_NID, 1)
            iLayer = C_CEL(TRACK_NID, 2)
            rSigTot = 0.0
            
            ! ADD THE ABSORPTION CROSS SECTION (COLUMN-2) THIS INCLUDES FISSION.
            if((iCell .gt. TXS_CEL) .or. (iCell .lt. -1)) then
               GetSigTot = 0.0
               return
            endif
            if((iLayer .gt. TXS_LAY) .or. (iLayer .lt. -1)) then
               GetSigTot = 0.0
               return
            endif
            if((nGroup .gt. TXS_GRP) .or. (nGroup .lt. 1)) then
               GetSigTot = 0.0
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
            GetSigTot = rSigTot
  
            return
         end function
         
!        ----------------------------------------------------------------------
!        THIS FUNCTION RETURNS THE ABSORPTION CROSS  SECTIONS  AS A FUNCTION OF
!        CURRENT NEUTRON ENERGY  GROUP AND  ALSO CURRENT CELL THAT CONTAINS THE 
!        NEUTRON.
!        ----------------------------------------------------------------------
         real function GetSigAbs()
            real :: rSigAbs
            integer :: nGroup
            integer :: iCell, iLayer
            
            
            ! OBTAIN CURRENT NEUTRON ENERGY GROUP AND CURRENT REACTOR CORE CELL
            nGroup = C_GRP(TRACK_NID, 1)
            iCell  = C_CEL(TRACK_NID, 1)
            iLayer = C_CEL(TRACK_NID, 2)
           
            
            ! IF THE NEUTRON IS CURRENTLY IN THE GRAPHITE REFLECTOR THEN SELECT 
            ! THE UNIQUE CELL ID FOR GRAPHITE REFLECTOR.
            if((iLayer .eq. -1) .and. (iCell .eq. -1)) then
               rSigAbs = TXS_TABLE(iLayer, REFLECTOR_CELLID, nGroup, 2)
            else
               if((iLayer .gt. 0) .and. (iCell .gt. 0)) then
                  rSigAbs = TXS_TABLE(iLayer, iCell, nGroup, 2)
               else
                  rSigAbs = 0.0
               endif
                        
            endif

            ! IN CASE IF THE NEUTRON CROSS SECTION READ IS LESS THAN ZERO, WE
            ! CORRECT IT BY SETTING IT TO ZERO.
            if(rSigAbs .lt. 0) rSigAbs = 0.0
            GetSigAbs = rSigAbs
            
            return
            
         end function
         
!        ----------------------------------------------------------------------
!        THIS FUNCTION RETURNS THE AVERAGE NUMBER OF  FISSION NEUTRONS PRODUCED
!        PER FISSION AS A FUNCTION OF THE ENERGY  GROUP  OF THE CURRENT NEUTRON
!        THAT IS BEING TRACKED.
!        ----------------------------------------------------------------------
         real function GetNuBar()
            ! DECLARING INTEGER THAT STORES CURRENT NEUTRON ENERGY GROUP
            integer :: iGroup
            
            ! OBTAIN THE CURRENT ENERGY GROUP FROM THE NEUTRON BANK.
            iGroup   = C_GRP(TRACK_NID, 1)
            
            ! WE BEGIN CALCULATE THE NU-BAR USING THE CORRELATED FUNCTION GIVEN
            ! BY TRKOV, INSTITUTE OF JOZEF STEFAN, 1996.
            GetNuBar = 2.55 - 0.11 * float(iGroup - 1) / 
     &                 float(max(1, TXS_GRP - 1))
     
            return
            
         end function
         
!        ----------------------------------------------------------------------
!        THIS FUNCTION RETURNS THE FISSION CROSS SECTION AS  A  FUNCTION OF THE 
!        CURRENT NEUTRON ENERGY GROUP, CURRENT  CELL ID  AND  CURRENT LAYER ID.
!        ----------------------------------------------------------------------
         real function GetSigFis()
         
            real    :: rSigFis
            integer :: iGroup
            integer :: iCell, iLayer

            ! HERE WE OBTAIN THE CURRENT ENERGY GROUP AND  CURRENT  CELL ID AND
            ! LAYER ID FROM THE NEUTRON BANK.
            iGroup = C_GRP(TRACK_NID, 1)
            iCell  = C_CEL(TRACK_NID, 1)
            iLayer = C_CEL(TRACK_NID, 2)
            
            ! IMPORTANT NOTE:  IN THE TXS TAPE,  THE  FISSION CROSS  SECTION IS
            ! GIVEN IN TERMS OF FISSION YIELD.  FISSION YIELD IS DEFINED AS THE 
            ! PRODUCT  OF  FISSION  CROSS SECTION  AND  THE  AVERAGE  NUMBER OF
            ! FISSION NEUTRON.  HENCE, TO OBTAIN  THE FISSION CROSS SECTION, WE
            ! HAVE TO  DIVIDE  THE YIELD  (OBTAINED FROM THE 3-RD COLUMN OF THE 
            ! TXS TAPE)  WITH THE AVERAGE  NUMBER OF FISSION NEUTRON.
            !
            ! IF THE NEUTRON IS CURRENTLY INSIDE THE GRAPHITE REFLECTOR (iLayer
            ! = -1 AND iCell = -1), THEN WE  OBTAIN  THE DATA WHICH IS SPECIFIC 
            ! FOR THE REFLECTOR. FOR THE REFLECTOR,  THE CELL ID MAY BE DIFFER-
            ! RENT, THUS IT IS MANDATORY TO SET THE CELL ID, iCell = REFLECTOR_
            ! CELLID. REFLECTOR_CELLID  IS DEFINED  IN THE  GEOMETRY MODULE, IT 
            ! STORES THE CELL ID OF THE GRAPHITE REFLECTOR.
            if((iLayer .eq. -1) .and. (iCell .eq. -1)) then
               rSigFis = TXS_TABLE(iLayer, REFLECTOR_CELLID,
     &                   iGroup, 3) / GetNuBar()
            else
            

                  rSigFis = TXS_TABLE(iLayer, iCell, iGroup, 3)
     &                   / GetNuBar()                  

       
            endif

            ! IN CASE IF OUR THE  DATA OBTAINED IS LESS THAN ZERO,  WE JUST SET
            ! IT TO ZERO. 
            if(rSigFis .lt. 0) rSigFis = 0.0
            GetSigFis = rSigFis
            
            return
            
         end function
         
!        ----------------------------------------------------------------------
!        THIS FUNCTION RETURNS THE SCATTERING CROSS  SECTION  AS  A FUNCTION OF
!        INCOMING NEUTRON ENERGY GROUP AND THE OUTGOING ENERGY GROUP.
!        ----------------------------------------------------------------------
         real function GetSigSca(piIncomingGroup, piOutgoingGroup)
            integer, intent(in) :: piIncomingGroup
            integer, intent(in) :: piOutGoingGroup
            
            real    :: rSigSca
            integer :: nGroup
            integer :: iCell, iLayer
            
            if(piIncomingGroup .gt. TXS_GRP) then
               GetSigSca = 0.0
               return
            endif
            ! OBTAIN THE CURRENT NEUTRON ENERGY GROUP, AND CELL INDICES.
            nGroup = C_GRP(TRACK_NID, 1)
            iCell  = C_CEL(TRACK_NID, 1)
            iLayer = C_CEL(TRACK_NID, 2)
            
            ! SINCE THE SCATTERING CROSS SECTION OF THE FIRST GROUP IS LOCATED 
            ! AT THE 6-TH COLUMN, WE HAVE TO OFFSET THE COLUMN TO +5.
            !
            ! IF THE NEUTRON IS CURRENTLY INSIDE THE GRAPHITE REFLECTOR (iLayer
            ! = -1 AND iCell = -1), THEN WE  OBTAIN  THE DATA WHICH IS SPECIFIC 
            ! FOR THE REFLECTOR. FOR THE REFLECTOR,  THE CELL ID MAY BE DIFFER-
            ! RENT, THUS IT IS MANDATORY TO SET THE CELL ID, iCell = REFLECTOR_
            ! CELLID. REFLECTOR_CELLID  IS DEFINED  IN THE  GEOMETRY MODULE, IT 
            ! STORES THE CELL ID OF THE GRAPHITE REFLECTOR.
            if((iLayer .eq. -1) .and. (iCell .eq. -1)) then

               rSigSca = TXS_TABLE(iLayer, REFLECTOR_CELLID,
     &                   piIncomingGroup, piOutgoingGroup+5)

            else
               

                  rSigSca = TXS_TABLE(iLayer, iCell,
     &                   piIncomingGroup, piOutgoingGroup+5)               


     
            endif

            ! IN CASE IF THE NEUTRON CROSS SECTION READ IS LESS THAN ZERO, WE
            ! CORRECT IT BY SETTING IT TO ZERO.            
            if(rSigSca .lt. 0) rSigSca = 0.0
            GetSigSca = rSigSca
            
            return
            
         end function

!        ----------------------------------------------------------------------
!        THIS FUNCTION RETURNS THE TOTAL SCATTERING CROSS SECTION AS A FUNCTION
!        OF CURRENT NEUTRON ENERGY GROUP AND CURRENT CELL INDICES.
!        ----------------------------------------------------------------------
         real function GetSigTotSca()
         
            real    :: rSigSca
            integer :: iGroup
            integer :: iCell, iLayer
            integer :: j
            

            ! OBTAIN THE CURRENT NEUTRON ENERGY GROUP, AND CELL INDICES.
            iGroup = C_GRP(TRACK_NID, 1)
            iCell  = C_CEL(TRACK_NID, 1)
            iLayer = C_CEL(TRACK_NID, 2)

            if(iGroup .gt. TXS_GRP) then
               GetSigTotSca = 0.0
               return
            endif            
            ! SINCE THE SCATTERING CROSS SECTION OF THE FIRST GROUP  IS LOCATED 
            ! AT THE 6-TH COLUMN,  WE HAVE TO OFFSET THE  COLUMN  TO j+5.  THIS
            ! TIME WE HAVE TO SUM ALL OF THE INDIVIDUAL GROUP SCATTERING  CROSS
            ! SECTIONS, WHICH INCLUDES THE  SUMMATION  FROM GROUP-1 TO GROUP-N.
            ! IN  TXSREADER  MODULE,  THE  TOTAL  NUMBER OF  NEUTRON  GROUPS IS
            ! STORED IN TXS_GRP.  SO WE SUM FROM GROUP 1  TO  TXS_GRP.
            !
            ! IF THE NEUTRON IS CURRENTLY INSIDE THE GRAPHITE REFLECTOR (iLayer
            ! = -1 AND iCell = -1), THEN WE  OBTAIN  THE DATA WHICH IS SPECIFIC 
            ! FOR THE REFLECTOR. FOR THE REFLECTOR,  THE CELL ID MAY BE DIFFER-
            ! RENT, THUS IT IS MANDATORY TO SET THE CELL ID, iCell = REFLECTOR_
            ! CELLID. REFLECTOR_CELLID  IS DEFINED  IN THE  GEOMETRY MODULE, IT 
            ! STORES THE CELL ID OF THE GRAPHITE REFLECTOR.      
            rSigSca = 0.0
            if((iLayer .eq. -1) .and. (iCell .eq. -1)) then

               do j=1, TXS_GRP, 1
                  rSigSca = rSigSca + 
     &                      TXS_TABLE(iLayer, REFLECTOR_CELLID,
     &                                iGroup, j+5)
               enddo
               
            else
               

                  do j=1, TXS_GRP, 1
                     rSigSca = rSigSca + 
     &                      TXS_TABLE(iLayer, iCell, iGroup, j+5)  
                  enddo           

                     
            endif

            ! IN CASE IF THE NEUTRON CROSS SECTION READ  IS LESS  THAN ZERO, WE
            ! CORRECT IT BY SETTING IT TO ZERO.  
            if(rSigSca .lt. 0) rSigSca = 0.0
            GetSigTotSca = rSigSca
            
            return
            
         end function
      
!        ----------------------------------------------------------------------
!        THIS FUNCTION WILL DECIDE WHETHER THE CURRENT   NEUTRON AT A COLLISION
!        SITE WILL FISSION OR CAPTURED.
!        ----------------------------------------------------------------------
         logical function WillFission()
         
            real :: rSigNeutronDie
            real :: rSigFis
            real :: rSigAbs
            real :: rP(3)
            real :: rEpsilon
            
            ! WE BEGIN SELECTING EITHER ABSORB OR FISSION BY TOSSING A RANDOM 
            ! NUMBER, ε. 
            rEpsilon = Rnd(int(TRACK_NID, 8))
            
            ! OBTAIN THE FISSION CROSS SECTION AND THE ABSORPTION CROSS SECTION
            ! FOR CURRENT REACTOR CORE CELL.
            rSigFis = GetSigFis()
            rSigAbs = GetSigAbs()
            
            ! CALCULATING THE CAPTURE CROSS SECTION.
            rSigNeutronDie = rSigAbs - rSigFis
            
            ! THEN WE PREPARE THE CUMMULATIVE DIST. FUNC., rP(i)            
            rP(1) = 0.0
            rP(2) = (rSigNeutronDie / rSigAbs)  + rP(1)
            rP(3) = (rSigFis / rSigAbs)         + rP(2)
            
            ! THEN BY USING THE INVERSION METHOD,  WE SELECT WHICH REACTION 
            ! REGIME CORRESPONDS TO THE ε VALUE. IF rP(1) ≤ ε < rP(2) THEN
            ! THE ABSORPTION REACTION IS POSSIBLE.
            !
            !             CUMMULATIVE DIST. FUNC. 
            !         
            !         ^
            !         |                      rP(3)
            !      1.0|.......................x 
            !         |                       : XXXXXXXXX
            !         |         rP(2)         : XXXXXXXXX
            !      0.7|...........x           : XXXXXXXXX
            !         |           :           : XXXXXXXXX
            !  ε----->|           :           : XXXXXXXXX
            !         |           :           : XXXXXXXXX
            !         |           :           : XXXXXXXXX
            !         |rP(1)      :           : XXXXXXXXX
            !       --+-----------+-----------+-----------> rxn
            !         |  ABSORB   |  FISSION  | 


            if((rEpsilon .ge. rP(1)) .and. 
     &         (rEpsilon .lt. rP(2))) then
     
               ! CASE FOR PURE ABSORPTION/CAPTURE (DIE) REACTION
               WillFission = .false.
               return
               
            elseif((rEpsilon .ge. rP(2)) .and.
     &             (rEpsilon .lt. rP(3))) then
     
               ! CASE FOR FISSION REACTION
               WillFission = .true.
               return      
               
            endif   

            return
            
         end function
         
!        ----------------------------------------------------------------------
!        THIS FUNCTION WILL DECIDE WHETHER THE CURRENT   NEUTRON AT A COLLISION
!        SITE WILL SCATTER OR ABSORBED.
!        ----------------------------------------------------------------------
         logical function WillScatter()
            real    :: rSigTot
            real    :: rSigAbs
            real    :: rP(3)
            integer :: i, k
            real :: rEpsilon
            
            ! WE BEGIN SELECTING EITHER ABSORB  OR FISSION  BY TOSSING A RANDOM
            ! NUMBER, ε. 
            rEpsilon = Rnd(int(TRACK_NID,8))
            
            ! OBTAIN THE  TOTAL  CROSS SECTION AND THE ABSORPTION CROSS SECTION
            ! FOR CURRENT REACTOR CORE CELL.
            rSigTot = GetSigTot()
            rSigAbs = GetSigAbs()
            
            ! THEN WE PREPARE THE CUMMULATIVE DIST. FUNC., rP(i)             
            rP(1) = 0.0
            rP(2) = rSigAbs   / rSigTot + rP(1)   
            rP(3) = 1.0
            
            ! THEN BY  USING  THE  INVERSION METHOD,  WE SELECT  WHICH REACTION 
            ! REGIME  CORRESPONDS  TO  THE ε  VALUE.  IF rP(1) ≤ ε < rP(2) THEN
            ! THE  ABSORPTION  REACTION IS POSSIBLE.  IF rP(2) ≤ ε < rP(3) THEN
            ! THEN SCATTERING REACTION IS POSSIBLE.
            !
            !             CUMMULATIVE DIST. FUNC. 
            !         
            !         ^
            !         |                      rP(3)
            !      1.0|.......................x 
            !         |                       : XXXXXXXXX
            !         |         rP(2)         : XXXXXXXXX
            !      0.7|...........x           : XXXXXXXXX
            !         |           :           : XXXXXXXXX
            !  ε----->|           :           : XXXXXXXXX
            !         |           :           : XXXXXXXXX
            !         |           :           : XXXXXXXXX
            !        rP(1)        :           : XXXXXXXXX
            !       --x-----------+-----------+-----------> rxn
            !         |  ABSORB   |  SCATTER  | 
            
            if((rEpsilon .ge. rP(1)) .and. 
     &         (rEpsilon .lt. rP(2))) then
               ! CASE FOR ABSORPTION REACTION
               WillScatter = .false.
               return
            elseif((rEpsilon .ge. rP(2)) .and.
     &             (rEpsilon .lt. rP(3))) then
               ! CASE FOR SCATTERING REACTION
               WillScatter = .true.
               return      
            endif   
            
            return
         end function
         
!        ----------------------------------------------------------------------
!        THIS SUBROUTINE HANDLES THE SCATTERING  PROCESS OF THE CURRENT NEUTON.
!        THIS SUBROUTINE WILL AUTOMATICALLY TRANSLATES  THE POSITION AND DIREC-
!        TION OF THE NEUTRON AS WELL AS THE OUTGOING ENERGY OF THE NEUTRON.
!        ----------------------------------------------------------------------
         subroutine ScatterNeutron()
            
            ! rMu IS THE SCATTERING COSINE, μ.
            real :: rMu
            
            ! rRandPhi IS A RANDOM ANGLE, φ ∈ [0,2π).
            real :: rRandPhi
            
            ! THIS IS THE CURRENT INCIDENT NEUTRON DIRECTION, Ω(u,v,w).
            real :: rU, rV, rW
            real :: rX, rY, rZ
            
            ! THIS IS THE CURRENT OUTGOING NEUTRON DIRECTION, Ω'(u,v,w).
            real :: rUNew, rVNew, rWNew
            
            ! A DUMMY VARIABLE. DONT WORRY ABOUT THIS.
            integer :: iStatus
            
            ! DECLARING THE INCOMING AND OUTGOING NEUTRON ENERGY GROUP.
            integer :: iIncomingGroup, iOutgoingGroup, j
            
            ! DECLARING THE CUMM. PROB. DIST. TABLE.
            real :: rPGroup(TXS_GRP+1)
            
            ! DECLARING A VARIABLE TO STORE A RANDOM NUMBER.
            real :: rEpsilon
            integer :: iCell, iLayer
            real, parameter :: PI = 3.141592654
            
            ! OBTAIN CURRENT NEUTRON DIRECTION AND POSITION.
            rU = C_DIR(TRACK_NID, 1)
            rV = C_DIR(TRACK_NID, 2)
            rW = C_DIR(TRACK_NID, 3)
            rX = C_POS(TRACK_NID, 1)
            rY = C_POS(TRACK_NID, 2)
            rZ = C_POS(TRACK_NID, 3)
            
            ! GET THE MESH CELL INDICES OF THE CELL WHERE THE NEUTRON SUBSIDES.
            iCell = C_CEL(TRACK_NID, 1)
            iLayer = C_CEL(TRACK_NID, 2)
            
            ! WE PREPARE A RANDOM NUMBER TO SAMPLE FROM SCATTERING  CUMM. DIST.
            rEpsilon = Rnd(int(TRACK_NID,8))
            
            ! OBTAIN THE CURRENT NEUTRON GROUP AS THE INCOMING GROUP.
            iIncomingGroup = C_GRP(TRACK_NID, 1)
            
            ! WE PREPARE THE CUMMULATIVE DIST FUNC FOR THE GROUP SCATTERING.
            rPGroup(1) = 0.0
            do j=1, TXS_GRP, 1
               rPGroup(j+1) = GetSigSca(iIncomingGroup, j)
     &                        / GetSigTotSca()
     &                        + rPGroup(j)
            enddo
            
            ! WE BEGIN GAMBLING IN THE CASINO. THE WE KNOW THE FATE OF THE NEU-
            ! TRON. HERE WE OBTAIN THE OUTGOING ENERGY GROUP  AS  THE RESULT OF
            ! THE RANDOM NUMBER THAT WE TOSSED EARLIER.
            do j=1, TXS_GRP, 1
               if((rEpsilon .ge. rPGroup(j  )) .and.
     &            (rEpsilon .lt. rPGroup(j+1))) then
                  iOutgoingGroup = j
               endif
            enddo
            
            ! IF CURRENT CELL IS THE REFLECTOR,  THE SCATTERING IS ISOTROPIC IN
            ! THE CENTRE  OF MASS  FRAME  OF REFERENCE.  IF NOT, WE JUST ASSUME 
            ! ISOTROPIC SCATTERING IN LAB FRAME OF REFERRENCE.
            rMu = 2.0 * Rnd(int(TRACK_NID,8)) - 1.0
            
            ! IF CURRENT CELL IS  THE  REFLECTOR CELL, USE THE FORMULA GIVEN IN 
            ! DUDERSTAT & HAMILTON PAGE 134.
C            if((iCell .eq. -1) .and. (iLayer .eq. -1)) then
C               rMu = (1.0 + 12.0*rMu) /
C     &               sqrt(12.0**2 + 2.0*12.0*rMu + 1.0)
C            endif
            
c            if(TXS_CELLTYPE(iCell) .eq. '        LW') then
c               rMu = (1.0 + 2.0*rMu) /
c     &               sqrt(2.0**2 + 2.0*2.0*rMu + 1.0)            
c            endif
            ! HERE WE GENERATE  A RANDOM ANGLE ϕ ∈ [0,2π) TO CHANGE THE NEUTRON
            ! OUTGOING ANGLE.
            rRandPhi = 2.0 * PI * Rnd(int(TRACK_NID,8))      
            
            ! USING THE SAMPLED SCATTERING COSINE AND THE RANDOM ANGLE, WE CAL-
            ! CULATE THE NEW NEUTRON DIRECTION COMPONENT, Ω'(u,v,w).
            rUNew = rMu * rU + (sqrt(1.0 - rMu**2) * 
     &              (rU * rW * cos(rRandPhi) - rV * sin(rRandPhi))) 
     &              / sqrt(1.0 - rW**2)

            rVNew = rMu * rV + (sqrt(1.0 - rMu**2) * 
     &              (rV * rW * cos(rRandPhi) + rU * sin(rRandPhi))) 
     &              / sqrt(1. - rW**2)     

            rWNew = rMu * rW - sqrt(1.0 - rMu**2) *
     &              sqrt(1. - rW**2) * cos(rRandPhi)
     
            ! ENABLE THIS CODE IF WE ASSUME THAT THE REACTOR CORE BOUNDARY IS 
            ! 100% REFLECTIVE. DISABLE THIS CODE SECTION IF WE WANT NATURAL
            ! REFLECTOR BEHAVIOUR INSTEAD.
c            if((iCell .eq. -1) .and. (iLayer .eq. -1)) then
c               rUNew = rU - 2.0*(rX*rU+rY*rV)*rX /
c     &                 RingRadius(RingCount())**2
c               rVNew = rU - 2.0*(rY*rV+rX*rU)*rY /
c     &                 RingRadius(RingCount())**2
c               rWNew = rW
c            endif

            ! REQUEST THE NEUTRON BANK TO CHANGE THE NEUTRON DIRECTION.
            call SetNeutronDirection(TRACK_NID, rUNew, rVNew, rWNew,
     &                               iStatus)
            ! REQUEST THE NEUTRON BANK TO CHANGE THE NEUTRON ENERGY. 
            call SetNewEnergyGroup(TRACK_NID, iOutgoingGroup,
     &                             iStatus)
!            print'(A,I2,A,I2,A,I4,I4,A)',
!     &      'Scatter from group ', iIncomingGroup, ' to group ',
!     &              iOutgoingGroup, ' in cell ',
!     &               C_CEL(TRACK_NID,1), C_CEL(TRACK_NID, 2), '.'
            return
         end subroutine
         
!        ----------------------------------------------------------------------
!        NOTE: THIS METHOD IS DEPRICATED BUT RETAINED FOR FUTURE REFERENCE.
!        SEE THE NEW METHOD THAT HANDLES  FISSION REACTION,  KeffFission() IN 
!        CriticalityCalculation MODULE.
!        ----------------------------------------------------------------------         
         subroutine FissionReaction()
            real :: rP(TXS_GRP+1)
            integer :: iGroup
            real    :: rEpsilon
            integer :: iFissionEnergyGroup
            integer :: iStatus
            real    :: rRandomPhi
            real    :: rMu
            integer :: iCell, iLayer, j
            real    :: rU, rV, rW
            real    :: rUNew, rVNew, rWNew
            real, parameter :: PI = 3.141592654
            real    :: rINu, rNu
            integer :: iNu, k
            
            iGroup = C_GRP(TRACK_NID, 1)
            iCell  = C_CEL(TRACK_NID, 1)
            iLayer = C_CEL(TRACK_NID, 2)
            rU     = C_DIR(TRACK_NID, 1)
            rV     = C_DIR(TRACK_NID, 2)
            rW     = C_DIR(TRACK_NID, 3)
            
!           PREPARE THE FISSION NEUTRON COUNT
            rEpsilon = Rnd(int(TRACK_NID, 8))
            rNu = GetNuBar()
            rINu = floor(rNu)
            if(rEpsilon .le. (rNu - rINu)) then
               iNu = int(rINu) + 1
            elseif(rEpsilon .gt. (rNu - rINu)) then
               iNu = int(rINu)
            endif
            
!           LOOP FOR EACH CHILD NEUTRON (NOT APPLICABLE FOR K-EFF EIGENVALUE
!           CALCULATION!
            do k=1, iNu, 1
            
!           PREPARE THE OUTGOING NEUTRON ENERGY GROUP. 
            rEpsilon = Rnd(int(TRACK_NID, 8))     
            rP(1) = 0.0
            do j=1, TXS_GRP, 1
               rP(j+1) = TXS_TABLE(iLayer, iCell, j, 4) + rP(j)
            enddo
            
            do j=1, TXS_GRP, 1
               if((rEpsilon .ge. rP(j)) .and.
     &            (rEpsilon .lt. rP(j+1))) then
                  iFissionEnergyGroup = j
               endif
            enddo
            
            rMu = 2.0 * Rnd(int(TRACK_NID,8))  - 1.0
            rRandomPhi = 2.0 * PI * Rnd(int(TRACK_NID,8))      
            
            rUNew = rMu * rU + (sqrt(1.0 - rMu**2) * 
     &              (rU * rW * cos(rRandomPhi) - rV * sin(rRandomPhi))) 
     &              / sqrt(1.0 - rW**2)

            rVNew = rMu * rV + (sqrt(1.0 - rMu**2) * 
     &              (rV * rW * cos(rRandomPhi) + rU * sin(rRandomPhi))) 
     &              / sqrt(1. - rW**2)     

            rWNew = rMu * rW - sqrt(1.0 - rMu**2) *
     &              sqrt(1. - rW**2) * cos(rRandomPhi)
            call AddSecondaryNeutron(C_POS(TRACK_NID, 1),
     &                               C_POS(TRACK_NID, 2),
     &                               C_POS(TRACK_NID, 3),
     &                               rUNew, rVNew, rWNew,
     &                               iFissionEnergyGroup,
     &                               TRACK_NID, RXN_FISSION, 1.0)
         enddo
!        END LOOP FOR EACH NEW BORN NEUTRON

            call KillNeutron(TXS_GRP, N_FISS, iStatus)
            return
         end subroutine
         
         
         subroutine AbsorbNeutron()
            integer :: iStatus
            call KillNeutron(TRACK_NID, N_DEAD, iStatus)
            return
         end subroutine
         
         subroutine NonAnalogAbsorption()
            real :: rSigTot, rSigAbs
            real :: rWeight, rWeightNew
            real :: rWeightCutOff  = 0.25
            real :: rWeightSurvive = 1.0
            integer :: iStatus
            real :: rEpsilon
            real :: rPKilled                 ! THE PROBABILITY OF GETTING KILL-
                                             ! ED IF THE WEIGHT IS BELOW CUTOFF
                                             
                                             
            rSigTot = GetSigTot()
            rSigAbs = GetSigAbs()
            rWeight = C_WGT(TRACK_NID, 1)
            rWeightNew = (1.0 - rSigAbs / rSigTot) * rWeight
            rPKilled = 1.0 - rWeight / rWeightSurvive
            if(rWeight .lt. rWeightCutOff) then
               rEpsilon = Rnd(int(TRACK_NID, 8))
               if(rEpsilon .lt. rPKilled) then
                  call KillNeutron(TRACK_NID, N_DEAD, iStatus)
               elseif((rEpsilon .ge. rPKilled) .and. 
     &                (rEpsilon .lt. 1.0)) then
                  
                  L_WGT(TRACK_NID, 1) = C_WGT(TRACK_NID, 1)
                  C_WGT(TRACK_NID, 1) = rWeightSurvive
               endif
            else
            
               L_WGT(TRACK_NID, 1) = C_WGT(TRACK_NID, 1)
               C_WGT(TRACK_NID, 1) = rWeightNew
            endif
            
         end subroutine
         
!        ----------------------------------------------------------------------
!        NOTE: THIS METHOD IS DEPRICATED BUT RETAINED FOR FUTURE REFERENCE.
!        SEE THE NEW METHOD THAT HANDLES NEUTRON TRANSPORT,  KeffTransport() IN 
!        CriticalityCalculation MODULE.
!        ----------------------------------------------------------------------             
         subroutine TransportNeutron()
            real :: rDistCol
            real :: rDistBoundary
            real :: rX, rY, rZ
            real :: rU, rV, rW
            
            integer :: iStatus
            
            
            rDistCol = DistanceToNextCollision()
            rDistBoundary = DistanceToNearestBoundary()
            rX = C_POS(TRACK_NID, 1)
            rY = C_POS(TRACK_NID, 2)
            rZ = C_POS(TRACK_NID, 3)
            rU = C_DIR(TRACK_NID, 1)
            rV = C_DIR(TRACK_NID, 2)
            rW = C_DIR(TRACK_NID, 3)
            
            if(rDistCol .gt. rDistBoundary) then  
!               print*, rDistBoundary
               rX = rX + rDistBoundary * rU
               rY = rY + rDistBoundary * rV
               rZ = rZ + rDistBoundary * rW
            else
               if(WillScatter() .eqv. .true.) then
                  call ScatterNeutron()
               else
                  if(WillFission() .eqv. .true.) then
                     call FissionReaction()
                  else 
                     call AbsorbNeutron()
                  endif
               endif
               rX = rX + rDistCol * rU
               rY = rY + rDistCol * rV
               rZ = rZ + rDistCol * rW
            endif
            call SetNeutronPosition(TRACK_NID, rX, rY, rZ, iStatus)
            ! call ScoreTrackLengthTally(TRACK_NID)
            return
         end subroutine


         
      end module
