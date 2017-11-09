

      subroutine GetNearestBoundaryDistance(prX, prY, prZ,    ! Neutrons position.
     &                                      prU, prV, prW,    ! Neutrons direction.
     &                                      prDistance)       ! Distance to the nearest boundary.
         
         implicit none
         
         real, intent(in)  :: prX, prY, prZ, prU, prV, prW
         real, intent(out) :: prDistance
         real, parameter :: PI = 3.141592654
         real, parameter :: INFINITY = 1.0/0.0
         real :: CellRadius, RingRadius, ThetaFromIndex, LayerHeight
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
         
!        Determine the current cell that contains the neutron.
         call GetCurrentCell(prX, prY, prZ, iC, iL)
!        Determine the cell index of the specified cell ID (iC).
         call CellIndex(iC, iI, iK)
!        We begin calculating for the planar surfaces perpendicular to z-axis.
         if(abs(prW) .gt. 0.0) then
            rDistanceA1 = (LayerHeight(iL)   - prZ) /  prW
            rDistanceB1 = (LayerHeight(iL-1) - prZ) / prW
         else
            rDistanceA1 = INFINITY
            rDistanceB1 = INFINITY
         endif
         
!        We begin calculating for the planar surface parallel to z-axis.
!        Calculating for Surface A

         rX1 = - RingRadius(iI) * cos(ThetaFromIndex(iI, iK))
         rY1 =   RingRadius(iI) * sin(ThetaFromIndex(iI, iK))
         rZ1 =   LayerHeight(iL)
         
         rX2 = - RingRadius(iI-1) * cos(ThetaFromIndex(iI, iK))
         rY2 =   RingRadius(iI-1) * sin(ThetaFromIndex(iI, iK))
         rZ2 =   LayerHeight(iL)
         
         rX3 = - RingRadius(iI) * cos(ThetaFromIndex(iI, iK))
         rY3 =   RingRadius(iI) * sin(ThetaFromIndex(iI, iK))
         rZ3 =   LayerHeight(iL-1)
         
         rAp = rY2*rZ3 + rY3*rZ1 + rY1*rZ2 - rY2*rZ1 - rY1*rZ3 - rY3*rZ2
         rBp = rZ2*rX3 + rZ3*rX1 + rZ1*rX2 - rZ2*rX1 - rZ1*rX3 - rZ3*rX2
         rCp = - rAp*rX1 - rBp*rY1
       
         rDistanceA2 = (rCp - rAp*prX - rBp*prY) / (rAp*prU + rBp*prV)
         if(rDistanceA2 .lt. 0.0) then
            rDistanceA2 = INFINITY
         endif
!        Calculating for Surface B
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

         rDistanceB2 = (rCp - rAp*prX - rBp*prY) / (rAp*prU + rBp*prV)
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
         rAc = prU**2 + prV*2
         rKc = prX*prU + prY*prV
         rCc = prX**2 + prY**2 - RingRadius(iI-1)**2

!        Calculate temporary solution-1
         rTemp1 = (-rKc + sqrt(rKc**2 - rAc*rCc)) / rAc
         rTemp2 = (-rKc - sqrt(rKc**2 - rAc*rCc)) / rAc  

!        We begin analyze which solution is the real solution.
         if(rCc .le. 0.0) then
!           Invalid condition here. 
            rDistanceA3 = INFINITY
         else
            if((rTemp1 .ge. 0.0) .and. (rTemp2 .ge. 0.0)) then
               rDistanceA3 = min(rTemp1, rTemp2)
            else
!              Invalid condition here. 
               rDistanceA3 = INFINITY              
            endif
         endif
         
!        Calculating for Surface B
         rAc = prU**2 + prV*2
         rKc = prX*prU + prY*prV
         rCc = prX**2 + prY**2 - RingRadius(iI)**2
        
!        Calculate temporary solution-1
         rTemp1 = (-rKc + sqrt(rKc**2 - rAc*rCc)) / rAc
         rTemp2 = (-rKc - sqrt(rKc**2 - rAc*rCc)) / rAc  
         
!        We begin analyze which solution is the real solution.
         if(rCc .gt. 0.0) then
!           Invalid condition here. 
            rDistanceB3 = INFINITY
         else
            if(rTemp1 .ge. 0.0) then
               rDistanceB3 = rTemp1
            elseif(rTemp1 .ge. 0.0) then
               rDistanceB3 = rTemp2 
            else
               rDistanceB3 = INFINITY
            endif
         endif         
1        continue


         prDistance = min(rDistanceA1, rDistanceB1,
     &                    rDistanceA2, rDistanceB2,
     &                    rDistanceA3, rDistanceB3)
         print*, rDistanceA1, rDistanceB1,
     &                    rDistanceA2, rDistanceB2,
     &                    rDistanceA3, rDistanceB3
         return
      end subroutine
      