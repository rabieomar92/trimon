!     --------------------------------------------------------
!     THIS FORTRAN SOURCE FILE WAS GENERATED USING FORTRAN-IT.
!     --------------------------------------------------------
!     Generation time: 6/8/2018 4:52:11 PM
!     Author         : 
!     Device Name    : RABIEMSI
!     Project Name   : CellGeometry.fip
!     --------------------------------------------------------
      
      
      program LibReader
         use LIBReaderM
         use ChannelList
!        Begin your code here...
         integer :: iStat, i, j, k
         integer :: iInt, iGrp, iType
         integer :: iLay, iCel, nCell, iFuelID
         
c         print*
c         print*, ' TRIMON-LIBREADER Version 1.181013.0058'
c         print*, ' (c) 2018 Malaysian Nuclear Agency'
c         print*, ' (c) 2018 Universiti Sains Malaysia'
c         print*
         open(unit=154, file='library\FE08.LIB',status='OLD', err=41)
         open(unit=155, file='library\FE12.LIB',status='OLD', err=42)
         open(unit=156, file='library\FE20.LIB',status='OLD', err=43)
         open(unit=157, file='library\GRAP.LIB',status='OLD', err=44)
c         open(unit=158, file='library\BERY.LIB',status='OLD', err=45)
         open(unit=159, file='library\COOL.LIB',status='OLD', err=46)
         open(unit=160, file='library\CHN1.LIB',status='OLD', err=47)
         open(unit=161, file='library\CHN2.LIB',status='OLD', err=48)
         open(unit=162, file='library\CHN3.LIB',status='OLD', err=49)
         open(unit=163, file='library\CHN4.LIB',status='OLD', err=50)
         open(unit=164, file='library\REFL.LIB',status='OLD', err=51)
         
         goto 52
41       print*, ' LIBR: Fatal error. FE08.LIB is missing.'     
         stop
42       print*, ' LIBR: Fatal error. FE12.LIB is missing.'     
         stop
43       print*, ' LIBR: Fatal error. FE20.LIB is missing.'     
         stop
44       print*, ' LIBR: Fatal error. GRAP.LIB is missing.'     
         stop
45       print*, ' LIBR: Fatal error. BERY.LIB is missing.'     
         stop
46       print*, ' LIBR: Fatal error. COOL.LIB is missing.'     
         stop
47       print*, ' LIBR: Fatal error. CHN1.LIB is missing.'     
         stop
48       print*, ' LIBR: Fatal error. CHN2.LIB is missing.'     
         stop
49       print*, ' LIBR: Fatal error. CHN3.LIB is missing.'     
         stop
50       print*, ' LIBR: Fatal error. CHN4.LIB is missing.'     
         stop
51       print*, ' LIBR: Fatal error. REFL.LIB is missing.'     
         stop 
52       continue        
         call INIT_LIB_READER()
         print*, ' LIBR: Reading FE08 element data.'
         call LIB_GET_TABLE(154, LIB_FE08, iStat)
         LIB_MASSFE08 = LIB_MASSU0
         print*, ' LIBR: Reading FE12 element data.'
         call LIB_GET_TABLE(155, LIB_FE12, iStat)
         LIB_MASSFE12 = LIB_MASSU0
         print*, ' LIBR: Reading FE20 element data.'
         call LIB_GET_TABLE(156, LIB_FE20, iStat)
         LIB_MASSFE20 = LIB_MASSU0
         print*, ' LIBR: Reading GRAP element data.'
         call LIB_GET_TABLE(157, LIB_GRAP, iStat)
c         print*, ' LIBR: Reading BERY element data.'
c         call LIB_GET_TABLE(158, LIB_BERY, iStat)
         print*, ' LIBR: Reading COOL element data.'
         call LIB_GET_TABLE(159, LIB_COOL, iStat)
         print*, ' LIBR: Reading CHN1 element data.'
         call LIB_GET_TABLE(160, LIB_CHN1, iStat)
         print*, ' LIBR: Reading CHN2 element data.'
         call LIB_GET_TABLE(161, LIB_CHN2, iStat)
         print*, ' LIBR: Reading CHN3 element data.'
         call LIB_GET_TABLE(162, LIB_CHN3, iStat)        
         print*, ' LIBR: Reading CHN4 element data.'
         call LIB_GET_TABLE(163, LIB_CHN4, iStat)  
         print*, ' LIBR: Reading REFL element data.'
         call LIB_GET_TABLE(164, LIB_G, iStat)
         close(unit=154)
         close(unit=155)
         close(unit=156)
         close(unit=157)
         close(unit=158)
         close(unit=159)         
         close(unit=160)
         close(unit=161)         
         close(unit=162)
         close(unit=163)
         close(unit=164)
         
         
         
         call CHL_ReadMainInput(.true.)
         print'(A,F0.3,A)', '  LIBR: The reactor is operating at ',
     &      CHL_COREPOWER, 'kW power.'
         print'(A,$)','  LIBR: Preparing cross section table.'
         open(unit=95, file='XSDATA.TXS',status='UNKNOWN')
!        HERE WE PRINT HEADER-8
         if(CHL_NRINGS .eq. 6) then
            write(95,'(5A10,3I10.9)') (' 888888888', i=1,5),
     &                             LIB_GRP, 92, CHL_NLAYERS
            nCell = 92
         elseif(CHL_NRINGS .eq. 7) then
            write(95,'(5A10,3I10.9)') (' 888888888', i=1,5),
     &                             LIB_GRP, 128, CHL_NLAYERS 
            nCell = 128
         endif
!        HERE WE BEGIN LOOP FOR EACH CORE LAYERS
         do iLay=1, CHL_NLAYERS, 1
         !  HERE WE WRITE HEADER-0 FOR EACH LAYER
            write(95,'(6A10,2I10.9)') (' 000000000', i=1,6),
     &                              nCell, iLay
         !  HERE WE LOOP FOR EACH CELL WITHIN A CORE LAYER.
            do iCel=1, nCell, 1
            !  HERE WE WRITE HEADER-1 FOR EACH CELL
               write(95,'(A10,2I10.9,A10)') ' 111111111',
     &                            iCel, LIB_GRP, CHL_FUELID(iCel)             

     
!           ---------------------------------------------------------
!           HERE IS THE CRUCIAL PART. WE PRINT THE CORRECTED XS HERE.
!           ---------------------------------------------------------
     
            if(CHL_FUELID(iCel) .eq. 'GRAP') then
            
               do iGrp=1, LIB_GRP, 1
                  write(95,'(8F10.7)') 0.0, 
     &            LIB_GETELEMENT(LIB_GRAP,1,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_GRAP,2,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_GRAP,3,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_GRAP,4,iLay,iCel,iGrp),
     &           (LIB_GETELEMENT(LIB_GRAP,4+k,iLay,iCel,iGrp),
     &            k=1, LIB_GRP),
     &            0.0
               enddo
               goto 162
               
            elseif(CHL_FUELID(iCel) .eq. 'COOL') then

               do iGrp=1, LIB_GRP, 1
                  write(95,'(8F10.7)') 0.0, 
     &            LIB_GETELEMENT(LIB_COOL,1,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_COOL,2,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_COOL,3,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_COOL,4,iLay,iCel,iGrp),
     &           (LIB_GETELEMENT(LIB_COOL,4+k,iLay,iCel,iGrp),
     &            k=1, LIB_GRP),
     &            0.0
               enddo
               goto 162
               
            elseif(CHL_FUELID(iCel) .eq. 'BERY') then
            
               do iGrp=1, LIB_GRP, 1
                  write(95,'(8F10.7)') 0.0, 
     &            LIB_GETELEMENT(LIB_BERY,1,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_BERY,2,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_BERY,3,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_BERY,4,iLay,iCel,iGrp),
     &           (LIB_GETELEMENT(LIB_BERY,4+k,iLay,iCel,iGrp),
     &            k=1, LIB_GRP),
     &            0.0
               enddo
               goto 162
               
            elseif(CHL_FUELID(iCel) .eq. 'G   ') then

               do iGrp=1, LIB_GRP, 1
                  write(95,'(8F10.7)') 0.0, 
     &            LIB_GETELEMENT(LIB_G,1,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_G,2,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_G,3,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_G,4,iLay,iCel,iGrp),
     &           (LIB_GETELEMENT(LIB_G,4+k,iLay,iCel,iGrp),
     &            k=1, LIB_GRP),
     &            0.0
               enddo
               goto 162
               
            elseif(CHL_FUELID(iCel) .eq. 'CHN1') then

                do iGrp=1, LIB_GRP, 1
                  write(95,'(8F10.7)') 0.0, 
     &            LIB_GETELEMENT(LIB_CHN1,1,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN1,2,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN1,3,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN1,4,iLay,iCel,iGrp),
     &           (LIB_GETELEMENT(LIB_CHN1,4+k,iLay,iCel,iGrp),
     &            k=1, LIB_GRP),
     &            0.0
               enddo
               goto 162
               
            elseif(CHL_FUELID(iCel) .eq. 'CHN2') then
      
               do iGrp=1, LIB_GRP, 1
                  write(95,'(8F10.7)') 0.0, 
     &            LIB_GETELEMENT(LIB_CHN2,1,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN2,2,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN2,3,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN2,4,iLay,iCel,iGrp),
     &           (LIB_GETELEMENT(LIB_CHN2,4+k,iLay,iCel,iGrp),
     &            k=1, LIB_GRP),
     &            0.0
               enddo
               goto 162
               
            elseif(CHL_FUELID(iCel) .eq. 'CHN3') then

               do iGrp=1, LIB_GRP, 1
                  write(95,'(8F10.7)') 0.0, 
     &            LIB_GETELEMENT(LIB_CHN3,1,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN3,2,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN3,3,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN3,4,iLay,iCel,iGrp),
     &           (LIB_GETELEMENT(LIB_CHN3,4+k,iLay,iCel,iGrp),
     &            k=1, LIB_GRP),
     &            0.0
               enddo
               goto 162
               
            elseif(CHL_FUELID(iCel) .eq. 'CHN4') then

               do iGrp=1, LIB_GRP, 1
                  write(95,'(8F10.7)') 0.0, 
     &            LIB_GETELEMENT(LIB_CHN4,1,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN4,2,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN4,3,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(LIB_CHN4,4,iLay,iCel,iGrp),
     &           (LIB_GETELEMENT(LIB_CHN4,4+k,iLay,iCel,iGrp),
     &            k=1, LIB_GRP),
     &            0.0
               enddo
               goto 162
            
            else
            
            !  ELSE (DESPITE ALL NON-FUEL CELL CHECK)
               read(CHL_FUELID(iCel),'(I4)') iFuelID
               iType = CHL_FUELTYPE(iFuelID)          
               do iGrp=1, LIB_GRP, 1
                  write(95,'(8F10.7)') 0.0, 
     &            LIB_GETELEMENT(iType,1,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(iType,2,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(iType,3,iLay,iCel,iGrp),
     &            LIB_GETELEMENT(iType,4,iLay,iCel,iGrp),
     &           (LIB_GETELEMENT(iType,4+k,iLay,iCel,iGrp),
     &            k=1, LIB_GRP),
     &            0.0
               enddo    
               goto 162
               
            endif
            

            
 162        continue           
            enddo ! FINISH CELL LOOP
            write(95,'(6A10,2I10.9)') (' 999999999', i=1,6),
     &                              nCell, iLay
         enddo ! FINISH LAYER LOOP
         
         close(unit=54)
         close(unit=55)
         close(unit=95)
         print*,'Done.'
         
      end program

