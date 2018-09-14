!                                                                       
!     Core Unit Cell Geometry Definition, Revision 180915-1.
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
!     This source code was created on 1/10/2017 11:15 PM by M. R. Omar.
!     Last revision date 15/9/2018.
!
!
!     Methods:
!
!     (1) SUBROUTINE DISTANCETOBOUNDARY(PRX, PR, PRZ, PRU, PRV, PRW, PRDISTANCE)
!         - RETRIEVES THE DISTANCE TO THE NEAREST BOUNDARY ALONG  THE DIRECTION
!           (PRU,PRV,PRW)  OF THE NEUTRONS  FLIGHT PATH.  (PRX,PRY,PRZ)  IS THE 
!           LOCATION OF THE NEUTRON.
!
!     (2) SUBROUTINE GETTARGETCELL(PRX, PRY, PRZ, PRU, PRV, PRW, PICELLID, PILAYERID) 
!         - RETRIEVES THE IDS OF THE CELL THAT CONTAINS THE NEUTRON AT POSITION 
!           (PRX,PRY,PRZ).  IF THE  NEUTRON  IS  LOCATED  AT  A  CELL INTERFACE, 
!           THIS  SUBROUTINE  RETURNS THE IDS OF THE CELL WHICH THE NEUTRONS IS 
!           HEADING TO. 
! 
!     (3) SUBROUTINE GETCURRENTCELL(PRX, PRY, PRZ, PICELLID, PILAYERID) 
!         - RETRIEVES  THE  CURRENT  CELL THAT CONTAINS THE NEUTRON AT POSITION 
!           (PRX,PRY,PRZ).  RETURN VALUES ARE  THE CELL ID  AND  THE LAYER ID OF 
!           THE CELL. 
!           RETURN  VALUES:  PICELLID  (1-127,-1,-2);  PILAYERID  (1-10,-1,-2) 
!           -2  MEANS  OUTSIDE OF  THE  REACTOR.  -1 MEANS WITHIN  THE  GRAPHITE 
!            REFLECTOR.
! 
!     (4) SUBROUTINE RADIUSTHETA(PICELLID, PRRADIUS, PRTHETA) 
!         - RETRIEVES  THE  OUTER RADIUS AND  THE ANGLE  OF  THE SPECIFIED CELL. 
!           THE CELL  IS  IDENTIFIED BY USING  THE  CELL ID.  APPLICABLE  TO ALL 
!           CELL IN ONE REACTOR CORE LAYER. 
! 
!     (5) REAL FUNCTION LAYERHEIGHT(PILAYERID) 
!         - THIS  FUNCTION RETURNS THE LAYER HEIGHT  OF THE SPECIFIED LAYER. THE
!           LAYER IS IDENTIFIED USING LAYER ID. 
! 
!     (6) REAL FUNCTION RINGRADIUS(PICELLID) 
!         - THIS  FUNCTION RETURNS  THE RADIUS OF THE  SPECIFIED  CELL.  APPLICA- 
!           BLE TO ANY LAYERS. THE CELL IS IDENTIFIED USING CELL ID. 
! 
!     (7) SUBROUTINE CELLIDFROMSITE(PCCELLSITE, PICELLID) 
!         - RETRIEVES THE CELL ID OF THE SPECIFIED CELL SITE NAME B-01, B-02...
! 
!     (8) SUBROUTINE CELLINDEX(PICELLID, PII, PIK) 
!         - RETRIEVES THE CELL INDICES PII AND PIK FROM CELL ID. 
! 
!     (9) SUBROUTINE CELLSITE(PICELLID, PCCELLSITE) 
!         - RETRIEVES THE CELL SITE NAME FROM CELL ID. 
! 
!    (10) REAL FUNCTION THETAFROMINDEX(PII, PIK) 
!         - THIS FUNCTION RETURNS THE ANGLE OF THE SPECIFIED CELL.  THE  CELL IS 
!           IDENTIFIED BY USING THE CELL INDICES PII AND PIK. 
      
      module CellGeometry
      
      implicit none
      
!     PLEASE MODIFY THESE PARAMETERS TO THE ACTUAL TRIGA CORE DIMENSIONS. WARNING,
!     DO NOT MODIFY OTHER PARTS OF THIS CODE.
         integer :: RING_COUNT          = 6
         integer :: REFLECTOR_CELLID    = 92
         integer :: ZERO_TRUNC          = 1E-13
         real    :: REFLECTOR_THICKNESS = 30.353
         real    :: TOPREFLECTOR_THICKNESS = 0.0
         real    :: BOTREFLECTOR_THICKNESS = 0.0
         real    :: FUEL_LENGTH         = 38.1
         integer :: LAYER_COUNT         = 20
         real :: RING_RADIUS(7)   = 
     &           (/ 9.08, 18.16, 27.24, 36.32, 45.40, 54.48, 0.0/)
c         real :: LAYER_HEIGHT(LAYER_COUNT) =
c     &           (/ 1*FUEL_LENGTH / real(LAYER_COUNT),
c     &              2*FUEL_LENGTH / real(LAYER_COUNT),
c     &              3*FUEL_LENGTH / real(LAYER_COUNT),
c     &              4*FUEL_LENGTH / real(LAYER_COUNT),
c     &              5*FUEL_LENGTH / real(LAYER_COUNT),
c     &              6*FUEL_LENGTH / real(LAYER_COUNT),
c     &              7*FUEL_LENGTH / real(LAYER_COUNT),
c     &              8*FUEL_LENGTH / real(LAYER_COUNT),
c     &              9*FUEL_LENGTH / real(LAYER_COUNT),
c     &             10*FUEL_LENGTH / real(LAYER_COUNT) /)
     

     
      contains

     
      real function LayerHeight(piLayerID)
         implicit none
         integer, intent(in) :: piLayerID
         LayerHeight = 0.0
         if((piLayerID .le. LAYER_COUNT) .and. (piLayerID .gt. 0)) then
            LayerHeight = real(piLayerID) * FUEL_LENGTH / 
     &                    real(LAYER_COUNT)
         endif
         return
      end function
   
      real function RingRadius(piRingID)
         implicit none
         integer, intent(in) :: piRingID
!         real :: rRingRadius(7) = (/ 3.149, 6.298, 9.447, 
!     &                              12.596, 15.745, 18.894,
!     &                              22.043 /)

         RingRadius = 0.0
         if(piRingID .gt. 0 .and. piRingID .le. RING_COUNT) then
            RingRadius = RING_RADIUS(piRingID)
            return
         endif
         return
      end function
      
      real function ReflectorThickness()
         ReflectorThickness = REFLECTOR_THICKNESS
         return
      end function
      
      real function TReflectorThickness()
         TReflectorThickness = TOPREFLECTOR_THICKNESS
         return
      end function
      
      real function BReflectorThickness()
         BReflectorThickness = BOTREFLECTOR_THICKNESS
         return
      end function
      
      integer function RingCount()
         RingCount = RING_COUNT
         return
      end function
      
      integer function LayerCount()
         LayerCount = LAYER_COUNT
         return
      end function

      subroutine GetCurrentCell(prX, prY, prZ, piCellID, piLayerID)

         implicit none

         real,    intent(in)  :: prX, prY, prZ
         integer, intent(out) :: piCellID, piLayerID
     
         integer :: nRing(7) = (/ 1, 6, 12, 18, 24, 30, 36 /)
         real    :: rRadius, rTheta, rZ, rTemp1, rTemp2
         integer :: iI, iK, iL, iID, i, j, k, m, iQuadrant

         piCellID  = -2
         piLayerID = -2

         rRadius = sqrt(prX**2 + prY**2)

!        ---------------------------------------------------------------
!        CENTRAL THIMBLE CASE. SET CELL ID TO 1, AND PROCEED TO LAYER ID 
!        CHECKING.  
!        ---------------------------------------------------------------    
         if(rRadius .eq. 0.0) then
            iID = 1
            piCellID = 1
            goto 2
         endif
!        ----------------------------------------------------------------
!        INDICATES THAT THE PARTICLE IS WITHIN REFLECTOR. RETURNS (-1,-1)
!        IF THE PARTICLE IS OUTSIDE THE REFLECTOR & CORE, RETURNS (-2,-2)
!        ----------------------------------------------------------------

c         if((prZ .lt. 0.0) .or. 
c     &      (prZ .gt. LayerHeight(LayerCount()))) then
c               piCellID  = -2
c               piLayerID = -2        
c               return            
c         endif
c         if((rRadius .ge. RingRadius(RingCount())) .and.
c     &      (rRadius .le. (RingRadius(RingCount()) + 
c     &                     ReflectorThickness()))   ) then
c            piCellID  = -1
c            if((prZ .ge. 0.0) .and. 
c     &         (prZ .le. LayerHeight(LayerCount()))) then
c               piLayerID = -1
c               return
c            else
c               piCellID  = -2
c               piLayerID = -2        
c               return
c            endif
c         elseif(rRadius .gt. (RingRadius(RingCount()) + 
c     &          ReflectorThickness())) then
c            piCellID  = -2
c            piLayerID = -2  
c            return
c         endif
         
         if((prZ .lt. -BReflectorThickness()) .or. 
     &      (prZ .gt. (LayerHeight(LayerCount()) + 
     &                 TReflectorThickness()))) then
               piCellID  = -2
               piLayerID = -2        
               return            
         endif
         if((rRadius .ge. RingRadius(RingCount())) .and.
     &      (rRadius .le. (RingRadius(RingCount()) + 
     &                     ReflectorThickness()))   ) then
            piCellID  = -1
            if((prZ .ge. -BReflectorThickness()) .and. 
     &         (prZ .le. (LayerHeight(LayerCount()) + 
     &                 TReflectorThickness()))) then
               piLayerID = -1
               return
            else
               piCellID  = -2
               piLayerID = -2        
               return
            endif
         elseif(rRadius .lt. RingRadius(RingCount())) then
            if((prZ .gt. LayerHeight(LayerCount())) .and.
     &         (prZ .lt. (LayerHeight(LayerCount()) + 
     &                    TReflectorThickness()      ))) then
               piCellID  = -1
               piLayerID = -1
               return
            elseif((prZ .gt. -BReflectorThickness()) .and.
     &             (prZ .lt. 0.0)) then
               piCellID  = -1
               piLayerID = -1    
               return
            endif
         elseif(rRadius .gt. (RingRadius(RingCount()) + 
     &          ReflectorThickness())) then
            piCellID  = -2
            piLayerID = -2  
            return
         endif

  
!        ----------------------------------------------------------------
!        INDICATES THAT THE PARTICLE IS WITHIN THE REACTOR CORE STD CELL.
!        ----------------------------------------------------------------         
         i = 1
         do i=1, RingCount(), 1
            if(i .eq. 1) then
               if((rRadius .lt. RingRadius(i)) .and.
     &            (rRadius .ge. 0.0)) then
                  iI = 1
                  goto 1
               endif
            else
               if((rRadius .le. RingRadius(i)) .and.
     &             (rRadius .ge. RingRadius(i-1))) then
                  iI = i
                  goto 1
               endif
            endif
         enddo
1        continue

!        CLASSIFY prX AND prY INTO THE APPROPRIATE QUADRANT.
         if((prX .lt. 0) .and. (prY .ge. 0)) iQuadrant = 1
         if((prX .ge. 0) .and. (prY .ge. 0)) iQuadrant = 2
         if((prX .ge. 0) .and. (prY .lt. 0)) iQuadrant = 3
         if((prX .lt. 0) .and. (prY .lt. 0)) iQuadrant = 4         
         select case (iQuadrant)
            case (1)
               rTheta = atan(abs(prY/prX))
            case (2)
               rTheta = 3.141592654 - atan(abs(prY/prX))
            case (3)
               rTheta = 3.141592654 + atan(abs(prY/prX))
            case (4)
               rTheta = 2.0*3.141592654 - atan(abs(prY/prX))
         end select
         do j=1, 127, 1
            call CellIndex(j,i,k)
!           IF THE ITERATION HAS ARRIVED AT RING iI...
            if(i .eq. iI) then
               if(k .eq. 1) then
                  rTemp1 = ThetaFromIndex(i, k)
                  rTemp2 = ThetaFromIndex(i, nRing(iI))
                  if((rTheta .le. rTemp1) .and. (rTheta .ge. 0.0)) then
                     iID = j
                     goto 2
                  elseif((rTheta .le. (2.0*3.141592654)) 
     &               .and. (rTheta .ge. rTemp2)) then
                     iID = j
                     goto 2                  
                  endif
               else
                  rTemp1 = ThetaFromIndex(i, k)
                  rTemp2 = ThetaFromIndex(i, k-1)       
                  if((rTheta .le. rTemp1) .and.
     &             (rTheta .ge. rTemp2)) then
                     iID = j
                     goto 2
                  endif       
               endif
            endif
         enddo
2        continue
         i = 1
         do i=1, LAYER_COUNT, 1
            if(i .eq. 1) then
               if((prZ .le. LayerHeight(i)) .and. 
     &            (prZ .ge. 0.0)) then
                  iL = 1
                  goto 3
               endif
            else
               if((prZ .le. LayerHeight(i)) .and.
     &            (prZ .gt. LayerHeight(i-1))) then
                  iL = i
                  goto 3
               endif
            endif
         enddo
3        continue      
         
         piCellID  = iID
         piLayerID = iL
         
         return

      end subroutine


!     THIS SUBROUTINE RETURNS THE R-COORDINATE AND THETA-COORDINATE OF A CELL
!     WITH CELL-ID = piCellID.      
      subroutine RadiusTheta(piCellID , prRadius, prTheta)
         implicit none
         integer, intent(in) :: piCellID 
         real, intent(out) :: prRadius, prTheta
         integer :: i,k
         prRadius = 0.0
         prTheta = 0.0
         if ((piCellID .gt. 127) .or. (piCellID .lt. 1)) return
         prRadius = RingRadius(piCellID )
         call CellIndex(piCellID , i, k)
         prTheta = ThetaFromIndex(i,k)
         return
      end subroutine

      
      real function CellRadius(piCellID)
         implicit none
         integer, intent(in) :: piCellID 
         integer :: iI, iK
         CellRadius = 0.0
         if ((piCellID .gt. 127) .or. (piCellID .lt. 1)) return
         call CellIndex(piCellID , iI, iK)
         if(iI .gt. 0 .and. iI .le. RingCount()) then
            CellRadius = RingRadius(iI)
            return
         elseif(iI .gt. RingCount()) then
            CellRadius = RingRadius(RingCount())
            return
         else
            return
         endif
      end function

!     THIS FUNCTION  RETURNS CELL ID FROM SITE NAME. FOR EXAMPLE, IF THE 
!     PROGRAMMER REQUESTS FOR THE CELL-ID OF B-01, THIS FUNCTION WILL
!     RETURN 2 (CELL ID FOR B-01 IS 2).
      subroutine CellIDFromSite(pcCellSite, piCellID)
         character(len=4), intent(in) :: pcCellSite
         integer,          intent(out) :: piCellID
         integer :: nCellPerRing(7) = (/ 1, 6, 12, 18, 24, 30, 36 /)
         integer :: iDummy, iReturnValue
         
         read(pcCellSite(3:4), '(I2)') iDummy
            
         if( pcCellSite(1:1) .eq. 'A') then
            iReturnValue = iDummy
         elseif( pcCellSite(1:1) .eq. 'B') then
            iReturnValue = iDummy + 1
         elseif( pcCellSite(1:1) .eq. 'C') then
            iReturnValue = iDummy + 1 + 6         
         elseif( pcCellSite(1:1) .eq. 'D') then
            iReturnValue = iDummy + 1 + 6 + 12
         elseif( pcCellSite(1:1) .eq. 'E') then
            iReturnValue = iDummy + 1 + 6 + 12 + 18     
         elseif( pcCellSite(1:1) .eq. 'F') then
            iReturnValue = iDummy + 1 + 6 + 12 + 18 + 24    
         elseif( pcCellSite(1:1) .eq. 'G') then
            iReturnValue = iDummy + 1 + 6 + 12 + 18 + 24 + 30
         endif
         
         piCellID = iReturnValue
         return
      end subroutine

!     THIS SUBROUTINE RETURNS THE CELL INDICES (I,K) FROM THE CELL-ID.
!     IF THE CELL ID IS 3 THEN IT RETURNS (2,2) (WHICH IS B-02).      
      subroutine CellIndex(piCellID , piI, piK)
      
         implicit none
         integer, intent(in) :: piCellID 
         integer, intent(out) :: piI, piK
         integer :: NRing(7) = (/ 1, 6, 12, 18, 24, 30, 36 /)
         integer :: r
         integer :: iBase
         
         iBase = 0
         do r=1, 7, 1
            iBase = iBase + NRing(r)
            if(iBase .ge. piCellID ) then
               piI = r
               goto 1
            endif
         enddo
1        continue
         piK = piCellID 
         do r=1, piI-1, 1
            piK = piK - NRing(r)
         enddo
         
         return
      end subroutine

!     THIS SUBROUTINE ALLOWS THE PROGRAMMERS TO RETRIEVE THE CELL SITE
!     NAME FROM CELL ID. FOR EXAMPLE, IF THE CELL ID IS 2 IT RETURNS 
!     B-01.
      subroutine CellSite(piCellID , pcCellTag)
      
         implicit none
         integer, intent(in) :: piCellID 
         character(len=4), intent(out) :: pcCellTag
         integer :: NRing(7) = (/ 1, 6, 12, 18, 24, 30, 36 /)
         integer :: r, iI, iK
         integer :: iBase
         
         pcCellTag = '    '
         IF ((piCellID .gt. 127) .or. (piCellID .lt. 1)) return
         
         iBase = 0
         do r=1, 7, 1
            iBase = iBase + NRing(r)
            if(iBase .ge. piCellID ) then
               iI = r
               goto 1
            endif
         enddo
1        continue
         iK = piCellID 
         do r=1, iI-1, 1
            iK = iK - NRing(r)
         enddo
         
         if(iI .eq. 1) pcCellTag(1:2) = 'A-'
         if(iI .eq. 2) pcCellTag(1:2) = 'B-'
         if(iI .eq. 3) pcCellTag(1:2) = 'C-'
         if(iI .eq. 4) pcCellTag(1:2) = 'D-'
         if(iI .eq. 5) pcCellTag(1:2) = 'E-'
         if(iI .eq. 6) pcCellTag(1:2) = 'F-'
         if(iI .eq. 7) pcCellTag(1:2) = 'G-'  
         write(pcCellTag(3:4),'(I2.2)') iK
         
         return
      end subroutine

      
      real function ThetaFromIndex(piI, piK)
         integer, intent(in) :: piI, piK
         if((piI .lt. 2) .or. (piI .gt. 7)) then
            ThetaFromIndex = 0.0
            return
         endif
         ThetaFromIndex = (3.141592654/6.0)*(2.0*real(piK)-1.0)
     &                         / (real(piI) - 1.0)
         return
      end function 

!     -------------------------------------------------------------------------
!     SUBROUTINE DistanceToBoundary
!
!     THIS ROUTINE CALCULATES THE DISTANCE TO THE NEAREST BOUNDARY OF THE  CORE 
!     CELL ALONG THE FLIGHT DIRECTION OF THE NEUTRON.
!
!     WARNING NOTICE:  THIS ROUTINE  IS  PART OF  THE  HEART OF THE MONTE CARLO
!     CODE.  PLEASE, BEFORE ATTEMPTING  TO  MODIFY THIS CODE  PLEASE  MAKE SURE
!     THIS  ROUTINE  IS  ACCURATE  AND  AS EFFICIENT  AS  POSSIBLE.  OTHERWISE,
!     THE WHOLE MONTE CARLO SIMULATION IS NOTHING MORE THAN A TRASH BIN.
!
!     WE DEFINE SURFACES THAT BOUNDS THE CELL AS FOLLOWING:
!     THERE ARE THREE DIFFERENT TYPES OF SURFACE THAT BOUNDS THE CELL.
!     FIRST  - THE PLANAR PLANES PERPENDICULAR TO THE Z-AXIS. THIS PLANE IS PA-
!              RALLEL  TO THE X-Y PLANE.  HENCE  THE CALCULATION IS DIRECT. SEE
!              DOCUMENTATION APPENDICES FOR THE EQUATION.
!     SECOND - THE PLANAR PLANES PARALLEL TO THE Z-AXIS. SEE FIGURE 1.
!     THIRD  - THE CYLINDRICAL CURVED SURFACE CENTERED ALONG THE Z-AXIS.
!              SEE FIGURE 1.
!
!         CYLINDRICAL _________
!          SURFACE B           \
!                               \
!                      __--------------__
!                   .**                  **.
!                   \        T H E        /
!      PLANAR   _____\      C E L L      /______  PLANAR
!     SURFACE B       \                 /        SURFACE A
!                      \               /
!                       \             /
!                        \__--***--__/
!                             \
!                              \_________ CYLINDRICAL
!                                         SURFACE A
!
!                 FIGURE 1: TOP VIEW OF THE CORE CELL.
!
!     GOOD LUCK, M. R. OMAR, OCTOBER 2017.
!     -------------------------------------------------------------------------
      real function DistanceToBoundary(prX, prY, prZ,    ! Neutrons position.
     &                              prU, prV, prW)       ! Distance to the nearest boundary.
         
         implicit none
         
         real, intent(in) :: prX, prY, prZ, prU, prV, prW
         real :: prDistance
         real, parameter :: PI = 3.141592654
         real, parameter :: INFINITY = 1.0/0.0

!        Declare cell indices I, K, cell ID and layer ID.
         integer :: iI, iK, iC, iL
!        Declare variables for Planar Planes Perpendicular to Z-Axis
         real :: rDistanceA1, rDistanceB1
!        Declare variables for Planar Planes Parallel to Z-Axis
         real :: rDistanceA2, rDistanceB2
         real :: rX1, rX2, rX3, rY1, rY2, rY3, rZ1, rZ2, rZ3
         real :: rAp, rBp, rCp
!        Declare variables for Cylindrical Surfaces
         real :: rDistanceA3, rDistanceB3, rTemp1, rTemp2
         real :: rAc, rCc, rKc
!        Declare variables for surface boundary correction.
         real :: rX, rY, rZ, j
         integer :: iC2, iL2, i
         logical :: lCorrectionDone = .false.
         
         rX = prX
         rY = prY
         rZ = prZ
         lCorrectionDone = .false.
!        ---------------------------------------------------------------------
!        Determine the current cell that contains the neutron.
!        HERE  WE  MAKE  SOME  BOUNDARY CROSSING CORRECTION. THIS WAS TEDIOUS. 
!        WHEN  I  CODE  THIS  PART,  ONLY GOD AND I KNOW WHAT AM I DOING.  BUT
!        AT THE MOMENT YOU READ THIS,  ONLY GOD  KNOWS WHAT  I  HAVE DONE.  SO
!        YOU AND I HAVE TO DISSECT THIS CODE CAREFULLY.
!
!        THIS PART IS IMPORTANT FOR CALCULATION ACCURACY. OUR PROBLEM IS THAT, 
!        WHEN A NEUTRON IS LOCATED EXACTLY AT THE CELL BOUNDARY,  AT A CERTAIN
!        NEUTRON  DIRECTION  THE CELL  CROSSING OCCURS OUTSIDE THE CELL.  THIS 
!        SEVERE ERROR CAN BE SEEN WHEN A NEUTRON IS  LOCATED AT  THE  BOUNDARY
!        IS MOVING AWAY FROM THE  CELL.  WE CALL THIS AS A BLIND-SPOT PROBLEM.
! 
!        THE  SOLUTION IS BY PEEKING ALONG THE DIRECTION OF NEUTRON FLIGHT AND
!        SEE WHETHER THE FLIGHT IS STILL IN THE SAME CELL. IF ITS NOT IN THE SAME
!        CELL, WE REGISTER THE NEW CELL AS THE CURRENT CELL THAT  CONTAINS THE
!        NEUTRON. WE PEEK FROM  A DISTANCE OF  10E-10, 10E-9, ..., 10E-5 ALONG
!        THE NEUTRONS FLIGHT DIRECTION (U,V,W)
!        ---------------------------------------------------------------------
2        continue
         call GetCurrentCell(rX, rY, rZ, iC, iL)
         j = -10.0
         do i=1,5,1
            call GetCurrentCell(rX+10.0**j*prU,
     &                          rY+10.0**j*prV,
     &                          rZ+10.0**j*prW, iC2, iL2)
            j = j + 1.0
            if((iC2 .ne. iC) .or. (iL2 .ne. iL)) then
               iC = iC2
               iL = iL2
               goto 3
            endif
         enddo
3        continue
!        Determine the cell index of the specified cell ID (iC).
         if(iC .gt. 0) then
            call CellIndex(iC, iI, iK)                    
         endif

!        -----------------------------------------------------------------------         
!        IF THE NEUTRON IS LOCATED WITHIN THE GRAPHITE REFLECTOR, iL=-1 & iC=-1.
!        -----------------------------------------------------------------------
   
         if(iC .eq. -1 .and. iL .eq. -1) then

!        PROCESSING THE INNER CYLINDRICAL SURFACE OF THE GRAPHITE REFLECTOR.
            rAc = prU**2 + prV**2
            rKc = rX*prU + rY*prV
            rCc = rX**2 + rY**2 - RingRadius(RingCount())**2
            rTemp1 = (-rKc + sqrt(rKc**2 - rAc*rCc)) / rAc
            rTemp2 = (-rKc - sqrt(rKc**2 - rAc*rCc)) / rAc
         
!        We begin analyze which solution is the real solution.
            if(rCc .le. 0.0) then
!              Invalid condition here. 
               rDistanceA3 = INFINITY
            else
               if((rTemp1 .gt. 0.0) .and. 
     &            (rTemp2 .gt. 0.0)) then
                  rDistanceA3 = min(rTemp1, rTemp2)
               else
!                 Invalid condition here. 
                  rDistanceA3 = INFINITY
               endif
            endif 

!        PROCESSING THE OUTER CYLINDRICAL SURFACE OF THE GRAPHITE REFLECTOR.            
            rAc = prU**2 + prV**2
            rKc = rX*prU + rY*prV
            rCc = rX**2 + rY**2 - (RingRadius(RingCount()) + 
     &            ReflectorThickness())**2
            rTemp1 = (-rKc + sqrt(rKc**2 - rAc*rCc)) / rAc
            rTemp2 = (-rKc - sqrt(rKc**2 - rAc*rCc)) / rAc    
            if(rCc .gt. 0.0) then
!              Invalid condition here. 
               rDistanceB3 = INFINITY
            else
               if(rTemp1 .gt. 0.0) then
                  rDistanceB3 = rTemp1
               elseif(rTemp2 .gt. 0.0) then
                  rDistanceB3 = rTemp2 
               else
                  rDistanceB3 = INFINITY
               endif
            endif 

!       PROCESSING THE TOP AND BOTTOM SURFACE OF THE REFLECTOR.  
            if(abs(prW) .gt. 0.0) then
               rDistanceA1 = (LayerHeight(LayerCount()) - rZ) /  prW
               rDistanceB1 = (0.0 - rZ) / prW
               if(rDistanceA1 .le. 0.0) rDistanceA1 = INFINITY
               if(rDistanceB1 .le. 0.0) rDistanceB1 = INFINITY
            else
               rDistanceA1 = INFINITY
               rDistanceB1 = INFINITY
            endif            
            
            prDistance = min(rDistanceA1, rDistanceB1,
     &                       rDistanceA3, rDistanceB3)   
!           NOTICE: HERE WE CORRECT FLOATING POINT ISSUES THAT CAUSES NEUTRON BOUNDARY TRAP.   
            if(prDistance .lt. ZERO_TRUNC) then
               prDistance = prDistance + ZERO_TRUNC
               DistanceToBoundary = prDistance
               return
            endif
            DistanceToBoundary = prDistance
            return     
!           NOTICE: END NEUTRON BOUNDARY TRAP CORRECTION
         endif

!        -----------------------------------------------------------
!        IF THE NEUTRON IS LOCATED WITHIN THE CENTRAL THIMBLE, iC=1.
!        -----------------------------------------------------------
         if(iC .eq. 1) then
            rAc = prU**2 + prV**2
            rKc = rX*prU + rY*prV
            rCc = rX**2 + rY**2 - RingRadius(1)**2
            rTemp1 = (-rKc + sqrt(rKc**2 - rAc*rCc)) / rAc
            rTemp2 = (-rKc - sqrt(rKc**2 - rAc*rCc)) / rAc
!        We begin analyze which solution is the real solution.
            if(rCc .gt. 0.0) then
               rDistanceB3 = INFINITY
            else
               if(rTemp1 .gt. 0) then
                  rDistanceB3 = rTemp1
               elseif(rTemp2 .gt. 0) then
                  rDistanceB3 = rTemp2 
               else
                  rDistanceB3 = INFINITY
               endif
            endif
            if(abs(prW) .gt. 0) then
               rDistanceA1 = (LayerHeight(iL)   - rZ) /  prW
               rDistanceB1 = (LayerHeight(iL-1) - rZ) / prW
               if(rDistanceA1 .le. 0.0) rDistanceA1 = INFINITY
               if(rDistanceB1 .le. 0.0) rDistanceB1 = INFINITY               
            else
               rDistanceA1 = INFINITY
               rDistanceB1 = INFINITY
            endif
            prDistance = min(rDistanceB3, rDistanceA1, rDistanceB1)
!           NOTICE: HERE WE CORRECT FLOATING POINT ISSUES THAT CAUSES NEUTRON 
!           BOUNDARY TRAP.   
            if(prDistance .lt. ZERO_TRUNC) then
               prDistance = prDistance + ZERO_TRUNC
               DistanceToBoundary = prDistance
               return
            endif
            DistanceToBoundary = prDistance
            return     
!           NOTICE: END NEUTRON BOUNDARY TRAP CORRECTION
         endif

!        --------------------------------------------------------
!        IF THE NEUTRON IS LOCATED WITHIN THE STANDARD CORE CELL.
!        --------------------------------------------------------
!        We begin calculating for the planar surfaces perpendicular to z-axis.
         if(abs(prW) .gt. 0.0) then
            rDistanceA1 = (LayerHeight(iL)   - rZ) /  prW
            rDistanceB1 = (LayerHeight(iL-1) - rZ) / prW
            if(rDistanceA1 .le. 0.0) rDistanceA1 = INFINITY
            if(rDistanceB1 .le. 0.0) rDistanceB1 = INFINITY
         else
            rDistanceA1 = INFINITY
            rDistanceB1 = INFINITY
         endif
         
!        We begin calculating for the planar surface parallel to z-axis.
!        Calculating for Surface A
!        Here we define all known points (cell vertices that span on the plane). 
         rX1 = - RingRadius(iI) * cos(ThetaFromIndex(iI, iK))
         rY1 =   RingRadius(iI) * sin(ThetaFromIndex(iI, iK))
         rZ1 =   LayerHeight(iL)
         
         rX2 = - RingRadius(iI-1) * cos(ThetaFromIndex(iI, iK))
         rY2 =   RingRadius(iI-1) * sin(ThetaFromIndex(iI, iK))
         rZ2 =   LayerHeight(iL)
         
         rX3 = - RingRadius(iI) * cos(ThetaFromIndex(iI, iK))
         rY3 =   RingRadius(iI) * sin(ThetaFromIndex(iI, iK))
         rZ3 =   LayerHeight(iL-1)
!        Then we calculate the coefficients A, B, C from the formula.         
         rAp = rY2*rZ3 + rY3*rZ1 + rY1*rZ2 - rY2*rZ1 - rY1*rZ3 - rY3*rZ2
         rBp = rZ2*rX3 + rZ3*rX1 + rZ1*rX2 - rZ2*rX1 - rZ1*rX3 - rZ3*rX2
         rCp = - rAp*rX1 - rBp*rY1
!        We calculate the distance to the surface using the equation:
!
!                             (C - Ax - By)   
!                         d = -------------
!                                Au + Bv   
!        for a neutron at (x,y,z) moving in (u,v,w) direction.

         rDistanceA2 = (rCp - rAp*rX - rBp*rY) / (rAp*prU + rBp*prV)
         if(rDistanceA2 .lt. 0.0) then
            rDistanceA2 = INFINITY
         endif
         
!        We begin calculating for Surface B.
         rX1 = - RingRadius(iI) * cos(ThetaFromIndex(iI, iK-1))
         rY1 =   RingRadius(iI) * sin(ThetaFromIndex(iI, iK-1))
         rZ1 =   LayerHeight(iL)
         
         rX2 = - RingRadius(iI-1) * cos(ThetaFromIndex(iI, iK-1))
         rY2 =   RingRadius(iI-1) * sin(ThetaFromIndex(iI, iK-1))
         rZ2 =   LayerHeight(iL)
         
         rX3 = - RingRadius(iI) * cos(ThetaFromIndex(iI, iK-1))
         rY3 =   RingRadius(iI) * sin(ThetaFromIndex(iI, iK-1))
         rZ3 =   LayerHeight(iL-1)
         
         rAp = rY2*rZ3 + rY3*rZ1 + rY1*rZ2 - rY2*rZ1 - rY1*rZ3 - rY3*rZ2
         rBp = rZ2*rX3 + rZ3*rX1 + rZ1*rX2 - rZ2*rX1 - rZ1*rX3 - rZ3*rX2
         rCp = - rAp*rX1 - rBp*rY1

         rDistanceB2 = (rCp - rAp*rX - rBp*rY) / (rAp*prU + rBp*prV)
         if(rDistanceB2 .lt. 0.0) then
            rDistanceB2 = INFINITY
         endif
!        We begin calculating for the cylindrical surfaces.
         if(prU .eq. 0.0 .and. prV .eq. 0.0) then
!           The neutrons flight will never cross the cylindrical surfaces.
!           Skip calculation. Goto 1.
            rDistanceA3 = INFINITY
            rDistanceB3 = INFINITY
            goto 1
         endif
!        Calculating for Surface A 
         rAc = prU**2 + prV**2
         rKc = rX*prU + rY*prV
         rCc = rX**2 + rY**2 - RingRadius(iI-1)**2
         rTemp1 = (-rKc + sqrt(rKc**2 - rAc*rCc)) / rAc
         rTemp2 = (-rKc - sqrt(rKc**2 - rAc*rCc)) / rAc  

!        We begin analyze which solution is the real solution.
         if(rCc .le. 0.0) then
!           Invalid condition here. 
            rDistanceA3 = INFINITY
         else
            if((rTemp1 .gt. 0.0) .and. 
     &         (rTemp2 .gt. 0.0)) then
               rDistanceA3 = min(rTemp1, rTemp2)
            else
!              Invalid condition here. 
               rDistanceA3 = INFINITY              
            endif
         endif
         
!        Calculating for Surface B
         rAc = prU**2 + prV**2
         rKc = rX*prU + rY*prV
         rCc = rX**2 + rY**2 - RingRadius(iI)**2       
         rTemp1 = (-rKc + sqrt(rKc**2 - rAc*rCc)) / rAc
         rTemp2 = (-rKc - sqrt(rKc**2 - rAc*rCc)) / rAc  
         
!        We begin analyze which solution is the true solution.
         if(rCc .gt. 0.0) then
!           Invalid condition here. 
            rDistanceB3 = INFINITY
         else
            if(rTemp1 .gt. 0.0) then
               rDistanceB3 = rTemp1
            elseif(rTemp2 .gt. 0.0) then
               rDistanceB3 = rTemp2 
            else
               rDistanceB3 = INFINITY
            endif
         endif         
1        continue


         prDistance = min(rDistanceA1, rDistanceB1,
     &                    rDistanceA2, rDistanceB2,
     &                    rDistanceA3, rDistanceB3)
!     NOTICE: HERE WE CORRECT FLOATING POINT ISSUES THAT CAUSES NEUTRON BOUNDARY TRAP.   
         if(prDistance .lt. ZERO_TRUNC) then
            prDistance = prDistance + ZERO_TRUNC
            DistanceToBoundary = prDistance
            return
         endif
         DistanceToBoundary = prDistance
         return
!        NOTICE: END NEUTRON BOUNDARY TRAP CORRECTION
      end function
      
      subroutine GetTargetCell(prX, prY, prZ, 
     &                         prU, prV, prW, 
     &                         piCellID, piLayerID)
         real, intent(in) :: prX, prY, prZ
         real, intent(in) :: prU, prV, prW
         integer, intent(out) :: piCellID, piLayerID
         integer :: iLTemp, iCTemp, i
         real :: e
         
         call GetCurrentCell(prX, prY, prZ, piCellID, piLayerID)
         e = -10.0
         do i=1, 5, 1
            call GetCurrentCell(prX+10.0**e*prU,
     &                          prY+10.0**e*prV,
     &                          prZ+10.0**e*prW, iCTemp, iLTemp)
            e = e + 1.0
            if((iCTemp .ne. piCellID) .or. 
     &         (iLTemp .ne. piLayerID)) then
               piCellID = iCTemp
               piLayerID = iLTemp
               return
            endif
         enddo
         return
      end subroutine

      real function CellVolume(piCellID, piLayerID)
      
         integer, intent(in) :: piCellID
         integer, intent(in) :: piLayerID
         real, parameter :: PI = 3.141592654
         
         integer :: iI, iK, iL
         real    :: rRi, rRi1
         real    :: rThk, rThk1
         real    :: rHl, rHl1
         
         
         call CellIndex(piCellID , iI, iK)
         iL = piLayerID
         
         if((piCellID .eq. 1) .and. (piLayerID .eq. 1)) then
         
            rRi = RingRadius(1)
            rHl = LayerHeight(1)
            
            CellVolume = PI * rRi**2 * rHl
            return
            
         elseif((piCellID .gt. 1) .and. (piCellID .lt. 128) .and.
     &          (piLayerID .gt. 0) .and. 
     &          (piLayerID .lt. LayerCount())) then
         
            rRi   = RingRadius(iI)
            rRi1  = RingRadius(iI-1)
            rThk  = ThetaFromIndex(iI, iK)
            rThk1 = ThetaFromIndex(iI, iK-1)
            rHl   = LayerHeight(iL)
            rHl1  = LayerHeight(iL-1)
            CellVolume = 0.5 * (rThk-rThk1) * (rHl-rHl1) *
     &                  (rRi**2-rRi1**2)
            return
         endif
      end function
      end module
      
      
      
