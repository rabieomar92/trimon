!
!
!     ████████╗██╗  ██╗███████╗██████╗ ███████╗ █████╗ ██████╗ ███████╗██████╗ 
!     ╚══██╔══╝╚██╗██╔╝██╔════╝██╔══██╗██╔════╝██╔══██╗██╔══██╗██╔════╝██╔══██╗
!        ██║    ╚███╔╝ ███████╗██████╔╝█████╗  ███████║██║  ██║█████╗  ██████╔╝
!        ██║    ██╔██╗ ╚════██║██╔══██╗██╔══╝  ██╔══██║██║  ██║██╔══╝  ██╔══██╗
!        ██║   ██╔╝ ██╗███████║██║  ██║███████╗██║  ██║██████╔╝███████╗██║  ██║
!        ╚═╝   ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═════╝ ╚══════╝╚═╝  ╚═╝
!                                                                    
!     MODULE TXSREADER
!
!     THIS LIBRARY WAS WRITTEN BY M. R. OMAR. OCTOBER 2017.
!     WHEN I WROTE THIS CODE, ONLY GOD  AND  I UNDERSTOOD WHAT I WAS DOING. BUT
!     NOW, ONLY GOD KNOWS... HOWEVER, I WROTE SOME NOTES ON THIS CODE:
!
!     THE PURPOSE OF THIS LIBRARY IS TO  ENABLE   PROGRAMMERS  TO  RETRIEVE  DATA
!     FROM THE TRIGA CROSS SECTION TAPE (.TXS).  A  TXS TAPE  CONSISTS OF SEVERAL
!     NEUTRON ENERGY GROUP  CONSTANTS  TROUGHOUT TRIGA  MARK-II CORE.  A TXS TAPE
!     IS GENERATED  BY  THE WIMS-INTEGRATED TRIGA CORE HOMOGENISATION (W-I-T-C-H)
!     CODE.  THE DATA IS  REPRESENTED AS A MASTER TABLE WHICH CONSISTS OF SEVERAL 
!     FIX COLUMNS AND  SEVERAL  DYNAMIC COLUMNS. THERE ARE  SIX  (6) FIX COLUMNS:
!
!              COLUMN-1     NEUTRON DIFFUSION COEFFICIENT    (CM)
!              COLUMN-2     NEUTRON ABSORPTION CROSS SECTION (1/CM)
!              COLUMN-3     FISSION YIELD (1/CM)
!              COLUMN-4     NEUTRON FISSION GROUP FRACTION
!              COLUMN-5     FISSION CROSS SECTION (1/CM)
!              COLUMN-N     EDH FACTOR
!
!     THERE  ARE  N-6  DYNAMIC COLUMNS  AND  THE NUMBER OF COLUMNS DEPENDS ON THE
!     NUMBER OF NEUTRON GROUPS DEFINED IN THE NEUTRON CROSS SECTION DATA TAPE.
!     THE DYNAMIC COLUMNS STORES THE SCATTERING CROSS SECTION MATRIX. THE COLUMNS
!     START WITH THE SCATTERING FROM GROUP-1 TO GROUP-G, WHERE G = N-6.
!     
!     THE PRINCIPAL NEUTRON GROUP IS THE NUMBER OF NEUTRON  GROUPS REPRESENTED BY
!     THE WHOLE DATA IN THE TAPE.  THIS MEANS THAT ALL OF THE DATA DEFINED IN THE
!     TAPE MUST REFLECT TO THE SAME NEUTRON GROUPS COUNT.
!
!     THE PRINCIPAL CELL COUNT IS THE NUMBER OF REACTOR CORE CELL PER CORE LAYERS
!     REPRESENTED BY THE WHOLE DATA IN THE TAPE.  IF THE PRINCIPAL CELL  COUNT IS 
!     92,  THIS MEANS THAT  THERE MUST BE 92 CELLS  IN EACH  LAYER OF THE REACTOR 
!     CORE.
!
!
!               T R I G A   R E A C T O R   C O R E       **EACH CORE LAYER IS
!                                                           DENOTED WITH INTEGER
!              _..-------------------------------.._        ID. I.E. 1,2,3,10,239
!           .-*                                     *-.        |
!           |.                                       .|        |
!           | -.._________________________________..-*|  LAYER-1 (92-CELLS)
!           |                                         |
!            *-.._________________________________..-*|  LAYER-2 (92-CELLS)
!           |                                         |  
!            *-.._________________________________..-*|  LAYER-3 (92-CELLS)
!           |                                         |
!            *-.._________________________________..-*|  LAYER-4 (92-CELLS)
!           |                                         |
!            *-.._________________________________..-*|  LAYER-5 (92-CELLS)
!           |                                         |
!            *-.._________________________________..-*|  LAYER-6 (92-CELLS)
!           |                                         |
!            *-.._________________________________..-*|  LAYER-7 (92-CELLS) 
!           |                                         |
!            *-.._________________________________..-*|  LAYER-8 (92-CELLS) 
!           |                                         |
!            *-.._________________________________..-*|  LAYER-9 (92-CELLS)
!           |                                         |
!            *-.._________________________________..-*|  LAYER-10(92-CELLS) 
!           |                                         |
!            *-.._________________________________..-* 
!
!     
!     A TXS TAPE HAS THE FOLLOWING STRUCTURE WHICH IS ALSO KNOWN AS MASTER TABLE...
!
!     [              HEADER-8                ] _______________ HEADER-8 DEFINES THE
!     [              HEADER-0                ] _______         PRINCIPAL GROUP COUNT,
!     [              HEADER-1                ]        |        PRINCIPAL CELL COUNT, AND
!     [--------------------------------------]        |        THE PRINCIPAL NUMBER OF CORE LAYERS.
!     [             CELL-TABLE               ]____    |
!     [--------------------------------------]    |   |_______ HEADER-0 INDICATES THE BEGINNING
!     [              FOOTER-9                ]    |            OF LAYER DATA.
!     [              HEADER-1                ]    |
!     [--------------------------------------]    |___________ CELL TABLE CONSISTS OF VARIABLE
!     [             CELL-TABLE               ]                 COLUMNS THAT CONTAINS REAL NUMBER
!     [--------------------------------------]                 WITH 8F10.7 FORMAT.
!     [              FOOTER-9                ]
!     [              HEADER-1                ]
!     [--------------------------------------]
!     [             CELL-TABLE               ]
!     [--------------------------------------]
!     [              FOOTER-9                ]________________ FOOTER-9 INDICATES THE END OF THE
!     [              HEADER-1                ]                 LAYER DATA.
!     [--------------------------------------]
!     [             CELL-TABLE               ]
!     [--------------------------------------]
!     [              FOOTER-9                ]
!     [              HEADER-1                ]
!     [--------------------------------------]
!     [             CELL-TABLE               ]
!     [--------------------------------------]
!     [              FOOTER-9                ]
!     :                 :                    :
!     :                 :                    :
!     [         FOOTER-9 (LAST LAYER)        ]
!
!
!     THERE ARE MANY METHODS DEFINED IN THIS LIBRARY. FORTUNATELY, THERE
!     ARE ONLY TWO (2) IMPORTANT METHODS:
!     
!     (A) TXSGETTABLE(pNF, pnL, pnC, pnG, piStat)
!         THIS METHOD READS THE MASTER TABLE AND STORE DATA IN AN ARRAY
!         (TXS_TABLE)
!
!     (B) TXSVALIDATE(pNF, plFlag, plQuiet)
!         THIS METHOD CHECKS WHETHER A TXS FILE IS CORRUPTED. A CORRUPTED
!         TXS FILE HAS FORMAT ERRORS.
!
!     SO YOU DONT HAVE TO WORRY ABOUT OTHER METHODS. THEY ARE JUST ASSIST-
!     ING METHODS.
!
!     GOOD LUCK! M.R.OMAR, OCTOBER 2017.
!
!     READ_HEADER_8 SUBROUTINE ALLOWS THE PROGRAMMER TO READ THE
!     HEADER-8 OF THE OPENED TXS TAPE. THIS SUBROUTINE RETRIEVES
!     ALL OF THE PRINCIPAL COUNTS FOR NEUTRON  GROUP,  CELL  AND 
!     CORE LAYERS.
!        pNF IS THE UNIT NUMBER OF THE OPENED TXS FILE.
!        pnPrincipalGroup IS THE PRINCIPAL GROUP COUNT.
!        pnPrincipalCell IS THE PRICIPAL CELL COUNT.
!        pnLayer IS THE TOTAL NUMBER OF LAYERS OF THE REACTOR CORE.
!        piStat IS A FLAG INDICATING THE STATUS OF THE READ PROCESS.

      module TXSReader
      
         implicit none
         
         integer :: TXS_CEL
         integer :: TXS_LAY
         integer :: TXS_GRP
         
         real, allocatable :: TXS_TABLE(:,:,:,:)
      contains
      
      subroutine InitTXSReader(pNF)
         integer, intent(in) :: pNF
         integer :: iStatus
         call READ_HEADER_8(pNF, TXS_GRP, TXS_CEL, TXS_LAY, iStatus)
         if(iStatus .eq. 0) then
            allocate( TXS_TABLE(TXS_LAY, TXS_CEL,
     &                          TXS_GRP, TXS_GRP+6) )            
         endif
      end subroutine
      
      subroutine READ_HEADER_8(pNF, pnPrincipalGroup,
     &                         pnPrincipalCell, pnLayer, piStat)
      
         implicit none
         
         ! Define parameter(s) here...
         integer, intent(in)  :: pNF
         integer, intent(out) :: pnPrincipalGroup,
     &                           pnPrincipalCell,
     &                           pnLayer
         integer, intent(out) :: piStat
         integer :: ios
         character(len=80) :: cText
         character(len=10) :: cElem
         
         piStat = 0
         
         rewind pNF
91       read(pNF, '(A80)', iostat=ios, err=92) cText
         
         if ( flIsValid(8, cText) ) then
            read(cText(51:60),'(I10)') pnPrincipalGroup
            read(cText(61:70),'(I10)') pnPrincipalCell
            read(cText(71:80),'(I10)') pnLayer
            
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
      
            
      subroutine SEARCH_CORE_LAYER(pNF, pnLayer, piStat)
      
         implicit none
         
         ! Define parameter(s) here...
         integer, intent(in) :: pNF, pnLayer
         integer, intent(out) :: piStat
         integer :: ios, nLayer, iLine, nPrincipalGroup, nPricipalCell,
     &              nTotalLayers
         character(len=80) :: cText
         character(len=10) :: cElem
         integer :: nOffset

         call READ_HEADER_8(pNF, nPrincipalGroup, nPricipalCell,
     &                      nTotalLayers, piStat)
         if(piStat .lt. 0) return         
         
         nOffSet = (nPricipalCell + nPricipalCell * nPrincipalGroup * 2)
     &             * (pnLayer - 1) + 1
     
         piStat = 0
         iLine = 0
!        At the beginning, we rewind the tape and put the cursor at the zeroth line
!        of the TXS tape.
         rewind pNF
         call FSEEK(pNF, nOffSet, 0)
!        We read the next line of the TXS tape.
91       read(pNF, '(A80)', iostat=ios, err=92) cText
         iLine = iLine + 1
!        Then we check if the current line has reaches the end of tape (eof). If yes,
!        we end the search and set piStat to -1 to denotes search failure.
         if(ios .ne. 0) then
            piStat = -1
            goto 90
         endif
         goto 93
92       cText = ''
93       continue
         if ( flIsValid(0, cText) ) then
            read(cText(71:80),'(I10)') nLayer            
            if( nLayer .eq. pnLayer ) then
               piStat = iLine
               goto 90
            else
               piStat = -1
               goto 91
            endif
         else
            piStat = -1
            goto 91
         endif
         
90       continue
         return
      end subroutine

      
      
      subroutine SEARCH_CORE_CELL(pNF, pnLayer, pnCell, piStat)
         ! Define parameter(s) here...
         integer, intent(in) :: pNF, pnLayer, pnCell
         integer, intent(out) :: piStat
         integer :: ios, nCell, iLine
         character(len=80) :: cText
         character(len=10) :: cElem
               
         call SEARCH_CORE_LAYER(pNF, pnLayer, piStat)
         if(piStat .eq. -1) return

         iLine = piStat 
         piStat = 0
91       read(pNF, '(A80)', iostat=ios, err=92) cText
         iLine = iLine + 1
!        Now we check if eof is encountered. If eof is encountered, we exit
!        this subroutine by setting the search status (piStat) to -1 to denote
!        a search failure.
         if(ios .ne. 0) then
            piStat = -1
            goto 90
         endif
         goto 93
         
!        If we have trouble during reading the line, we handle this
!        by setting the line as a blank line.
92       cText = ''

93       continue
         
         if(flIsValid(1, cText)) then
!           If a HEADER-1 is found we check whether the HEADER-1 is the one that
!           we are looking for. If the HEADER-1 defines the neutron cross section
!           table for the desired cell, then we exit this subroutine (goto 90) 
!           and set the search status (piStat) to current line number (iLine) to 
!           denote a search success.
            read(cText(11:20),'(I10)') nCell
            if(nCell .eq. pnCell) then
               piStat = iLine
               goto 90
            else
!              Seems that the current cursor has not arrived at the FOOTER-9
!              of the current core layer, we continue search for HEADER-1 in
!              the next line. We continue read the next line and see if we 
!              could find HEADER-1. (goto 91)
               piStat = -1
               goto 91
            endif
         elseif(flIsValid(9, cText)) then
!           If the cursor has arrived at the FOOTER-9 of the current layer,
!           The the cross section table HEADER-1 is not found in the tape.
!           We end the search by exiting this subroutine (goto 90).
            piStat = -1
            goto 90
         else
            piStat = -1
            goto 91
         endif
         
90       continue
         return
      end subroutine

!     THIS METHOD IS DEPRICATED BECAUSE THE ALGORITHM FOR RETRIEVING DATA
!     FROM THE MASTER TABLE IS SIGNIFICANTLY SLOW. WE KEEP THIS FOR FUTURE
!     REFERENCE.    
      subroutine GET_DATA_ELEM(pNF, piLayer, piCell, piGroup,
     &                          piCol, prValue, piStat)
     
         implicit none
         
         ! Define parameter(s) here...
         integer, intent(in) :: pNF, piLayer, piCell, piCol, piGroup
         real   , intent(out) :: prValue
         integer, intent(out) :: piStat
         
         integer :: nPrincipalGroup, nPricipalCell, nTotalLayers
         integer :: iG
         character(len=80) :: cTextR1, cTextR2
         
         piStat = 0
         prValue = 0.0000000
         
         if( (piCol .lt. 1) .or. (piCol .gt. 10) ) then
            piStat = -2
            return
         endif
         
         call READ_HEADER_8(pNF, nPrincipalGroup, nPricipalCell,
     &                      nTotalLayers, piStat)
         if(piStat .lt. 0) return
         
         call SEARCH_CORE_CELL(pNF, piLayer, piCell, piStat)
         if(piStat .lt. 0) return
         
         do iG=1, nPrincipalGroup, 1
 
            if(iG .ne. piGroup) then

               read(pNF,'(A80)')
               read(pNF,'(A80)')
               piStat = piStat + 2
               
            elseif(iG .eq. piGroup) then
            
               read(pNF,'(A80)') cTextR1
               read(pNF,'(A80)') cTextR2
               if(piCol .le. 8) piStat = piStat + 1
               if(piCol .gt. 8) piStat = piStat + 2
               
               select case (piCol)
                  case (1)
                     read(cTextR1(1:10),'(F10.7)') prValue
                     return
                  case (2)
                     read(cTextR1(11:20),'(F10.7)') prValue
                     return                     
                  case (3)
                     read(cTextR1(21:30),'(F10.7)') prValue
                     return                  
                  case (4)
                     read(cTextR1(31:40),'(F10.7)') prValue
                     return                  
                  case (5)
                     read(cTextR1(41:50),'(F10.7)') prValue
                     return                  
                  case (6)
                     read(cTextR1(51:60),'(F10.7)') prValue
                     return                  
                  case (7)
                     read(cTextR1(61:70),'(F10.7)') prValue
                     return                  
                  case (8)
                     read(cTextR1(71:80),'(F10.7)') prValue
                     return                  
                  case (9)
                     read(cTextR2(1:10),'(F10.7)') prValue
                     return                  
                  case (10)
                     read(cTextR2(11:20),'(F10.7)') prValue
                     return                  
                  case default
                     piStat = -2
                     return
               end select
            endif      
         enddo

         return
      end subroutine
      

!     --------------------------------------------------------------
!     SUBROUTINE TXSGETTABLE(pNF, piStat)
!     --------------------------------------------------------------
!     THIS SUBROUTINE WAS WRITTEN BY M.R.OMAR ON 16TH OCTOBER 2017.
!     THIS SUBROUTINE ENABLES THE PROGRAMMER TO OBTAIN THE NEUTRON
!     CROSS SECTION DATA FROM A TXS TAPE. THE DATA WILL BE STORED
!     IN A TABLE ARRAY (TXS_TABLE).
!
!     TXS_TABLE CAN BE VISUALISED AS FOLLOWS:
!
!         T A B L E   TXS_TABLE(TXS_LAY,TXS_CEL,TXS_GRP,nColumn)
!      ______________________________________________________
!     |        |       |        |          |     |           |
!     | ILAYER | ICELL | IGROUP | COLUMN-1 | ... | COLUMN-N  |
!     |________|_______|________|__________|_____|___________|
!     |    1   |   1   |   1    | 0.112121 | ... | 1.232232  |
!     |    1   |   1   |   2    | 0.332321 | ... | 1.923221  |
!     :    :   :   :   :   :    :     :    :  :  :     :     :
!     |   10   |   92  |   4    | 0.465534 | ... | 1.223222  |
!     |________|_______|________|__________|_____|___________|
!
!     COLUMN-1 IS THE GROUP-AVERAGED DIFFUSION CONSTANT. (cm)
!     COLUMN-2 IS THE GROUP-AVERAGED ABSORPTION CROSS SECTION. (1/cm)
!     COLUMN-3 IS THE GROUP-AVERAGED FISSION YIELD (1/cm)
!     COLUMN-4 IS NEUTRON SPECTRUM GROUP FRACTION.
!     COLUMN-5 IS THE GROUP-AVERAGED FISSION CROSS SECTION (1/cm)
!     COLUMN-6 IS THE SCATTERING CROSS SECTION FROM (1 -> IGROUP)
!     COLUMN-7 IS THE SCATTERING CROSS SECTION FROM (2 -> IGROUP)
!                                :
!                                :
!     COLUMN-(N-1) IS THE SCATTERING CROSS SECTION FROM (N -> IGROUP)
!     COLUMN-N IS THE EDH FACTOR.
!
      subroutine TXSGETTABLE(pNF, piStat)
     
         implicit none
         
         ! Define parameter(s) here...
         integer, intent(in) :: pNF
         integer, intent(out) :: piStat
         logical :: lTapeValid
         integer :: nPrincipalGroup, nPricipalCell, nTotalLayers
         integer :: iG, iLay, iCel, j
         character(len=80) :: cTextR1, cTextR2
         character(len=10) :: iDummy
         
!        INITIALLY WE ASSUME THERE IS NO ERROR. WE SET piStat TO 0.
         piStat = 0

!        BEFORE WE PROCEED READING THE TAPE,  WE VERIFY WHETHER THE TAPE
!        IS VALID OR NOT. INVALID TAPE CAN CAUSE ERRORS WHEN READING THE
!        CORRUPTED TAPE.
         call TXSVALIDATE(pNF, lTapeValid, .true.)
         if(lTapeValid) print*, ' TAPE-READ: The tape is valid.'
         if(lTapeValid .eqv. .false.) then
            print*, ' TAPE-READ: The tape is not valid. Abort reading.'
            piStat = -1
            return
         endif
         
         print*, ' TAPE-READ: Loading data from the TXS tape.',
     &           ' Please wait.'
     
!        NEXT WE READ THE TXS TAPE HEADER-8. HEADER-8 CONTAINS THE INFOR-
!        MATION ABOUT THE TOTAL NUMBER OF NEUTRON GROUPS,  THE TOTAL NUM-
!        BER OF CELLS AND THE TOTAL NUMBER OF REACTOR CORE LAYERS.
!        EXAMPLE OF HEADER-8:
!        888888888 888888888 888888888 888888888   NGROUP  NCELL  NLAYER
         call READ_HEADER_8(pNF, nPrincipalGroup, nPricipalCell,
     &                      nTotalLayers, piStat)
         if(piStat .lt. 0) return
         
!        CHECK IF THE DIMENSION OF THE TABLE ARRAY SIZE PROVIDED BY  THE
!        INVOKER. IF THE ARRAY SIZE IS SMALL, THEN STOP. TELL USER THERE
!        IS INTERNAL ERROR. PROGRAMMERS MUST TAKE NOTE THAT THE SIZE  OF
!        TXS_TABLE HAS TO BE (pnL,pnC,pnG,10). 

         write(iDummy,'(I10)') nPrincipalGroup
         print*, ' TAPE-READ: The PRINCIPAL neutron group count is ',
     &       adjustl(iDummy//'.')
         write(iDummy,'(I10)') nPricipalCell
         print*, ' TAPE-READ: The PRINCIPAL cell count is ',
     &       adjustl(iDummy//'.')
         write(iDummy,'(I10)') nTotalLayers
         print*, ' TAPE-READ: The total number of core layers is ',
     &       adjustl(iDummy//'.')

!        WE IMPLEMENT NESTED LOOPS HERE TO  ITERATE THROUGH EVERY SINGLE
!        REACTOR CORE CELL AVAILABLE IN THE TAPE. THE ITERATIONS ARE LI-
!        MITTED UP TO MAX NUMBER OF CORE LAYERS (nTotalLayers)  AND  MAX
!        NUMBER OF CELLS - THE PRINCIPAL CELL COUNT (nPrincipalCell).
         do iLay=1, nTotalLayers, 1
         
            write(iDummy,'(I10)') iLay         
            print*, ' TAPE-READ: Loading data for reactor core layer-',
     &           adjustl(iDummy//'.')
     
         do iCel=1, nPricipalCell, 1
         
!        WE BEGIN INVOKE THE SEARCH_CORE_CELL METHOD TO  GET  THE  CORE 
!        CELL.         
         call SEARCH_CORE_CELL(pNF, iLay, iCel, piStat)
         if(piStat .lt. 0) return

!        FOR EACH GROUP DEFINED FOR THE CELL, WE READ ALL COLUMNS AND
!        ASSIGN THEM TO THE TABLE MATRIX.    
         do iG=1, nPrincipalGroup, 1
            read(pNF,'(8F10.7)') (TXS_TABLE(iLay, iCel, iG, j),
     &                            j=1, nPrincipalGroup+6)
         enddo
         
         enddo
         enddo
         print*, ' TXS-READ: Finish reading neutron cross section tape.'
         return
      end subroutine

      
      
!     THE METHODS BELOW THIS SECTION ALLOWS THE PROGRAMMERS TO VALIDATE
!     A TXS TAPE. PROGRAMMERS CAN UTILISE THESE METHODS TO CHECK FOR 
!     ERRORS AND INVALID FORMATTING IN THE TXS TAPE.
!     THE PROGRAMMER IS ADVISED NOT TO MODIFY  OTHER METHODS EXCEPT THE
!     FOLLOWING:
!
!        TXSVALIDATE(pNF, plFlag, plQuiet)
!
!     TXSVALIDATE ALLOWS THE PROGRAMMER TO VALIDATE THE PROGRAMMERS TO
!     VALIDATE THE TXS TAPE. pNF IS THE UNIT NUMBER OF THE TXS TAPE,
!     plFlag IS A FLAG INDICATING IF A TAPE IS CORRUPTED/IS INCORRECT
!     IN FORMAT. plQuiet ALLOWS THE METHOD TO SURPRESS THE VALIDATION
!     WARNING AND ERROR MESSAGES.
!
!     GOOD LUCK! REGARDS, M.R. OMAR.

      logical function flIsANumber(pcChar)
         implicit none
         character, intent(in) :: pcChar
         flIsANumber = .false.
         if (pcChar .eq. '0') flIsANumber = .true.
         if (pcChar .eq. '1') flIsANumber = .true.
         if (pcChar .eq. '2') flIsANumber = .true.
         if (pcChar .eq. '3') flIsANumber = .true.
         if (pcChar .eq. '4') flIsANumber = .true.
         if (pcChar .eq. '5') flIsANumber = .true.
         if (pcChar .eq. '6') flIsANumber = .true.
         if (pcChar .eq. '7') flIsANumber = .true.
         if (pcChar .eq. '8') flIsANumber = .true.
         if (pcChar .eq. '9') flIsANumber = .true.   
         return
      end function

      logical function flIsInteger(pcTXSElement)
         implicit none
         character(len=10), intent(in) :: pcTXSElement
         integer :: i, nDot
         flIsInteger = .true.
         if ( pcTXSElement(1:1) .ne. ' ') then
               flIsInteger = .false.
               return
         endif
         do i=2, 10, 1
            if ( flIsANumber(pcTXSElement(i:i)) .eqv. .false.) then
               flIsInteger = .false.
               return
            endif
         enddo
         return
      end function

      logical function flIsReal(pcTXSElement)
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
            flIsReal = .false.
            return
         endif
         
         ! Checking the first condition...
         if (flIsANumber(pcTXSElement(1:1)) .eqv. .true.) then
            lCondition1Passed = .true.
         elseif (pcTXSElement(1:1) .eq. ' ') then
            lCondition1Passed  = .true.
         elseif (pcTXSElement(1:1) .eq. '-') then
            lCondition1Passed = .true.
         else
            flIsReal = .false.
            return                  
         endif
         
         ! Checking the second condition...
         if (flIsANumber(pcTXSElement(2:2)) .eqv. .true.) then
            lCondition2Passed = .true.
         else
            flIsReal = .false.
            return                  
         endif
         
         ! Checking the third condition...
         if (pcTXSElement(3:3) .eq. '.') then
            lCondition3Passed = .true.
         else
            flIsReal = .false.
            return         
         endif
         
         ! Checking the fourth condition...
         do i=4, 10, 1
            if (flIsANumber(pcTXSElement(i:i)) .eqv. .true.) then
               lCondition4Passed = .true.
            else
               flIsReal = .false.
               return                  
            endif            
         enddo
         
         flIsReal = lCondition1Passed .and. 
     &              lCondition2Passed .and.
     &              lCondition3Passed .and.
     &              lCondition4Passed
         return
      end function

      logical function flIsValid(piCardType, pcTextLine)
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
         
         flIsValid = .false.
         
         select case (piCardType)
            case (-8)
               if( (flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (flIsReal(caElem(4)) .eqv. .true.) .and.
     &             (flIsReal(caElem(5)) .eqv. .true.) .and.
     &             (flIsReal(caElem(6)) .eqv. .true.) .and.
     &             (flIsReal(caElem(7)) .eqv. .true.) .and.
     &             (flIsReal(caElem(8)) .eqv. .true.)) then
                  
                  flIsValid = .true.
                  return
                  
               endif
               
            case(-1)
               if( (flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (caElem(2) .eq. '          ') .and.
     &             (caElem(3) .eq. '          ') .and.
     &             (caElem(4) .eq. '          ') .and.
     &             (caElem(5) .eq. '          ') .and.
     &             (caElem(6) .eq. '          ') .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  flIsValid = .true.
                  return
                  
               endif
            case(-2)
               if( (flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (caElem(3) .eq. '          ') .and.
     &             (caElem(4) .eq. '          ') .and.
     &             (caElem(5) .eq. '          ') .and.
     &             (caElem(6) .eq. '          ') .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  flIsValid = .true.
                  return
                  
               endif
               
            case(-3)
               if( (flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (caElem(4) .eq. '          ') .and.
     &             (caElem(5) .eq. '          ') .and.
     &             (caElem(6) .eq. '          ') .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  flIsValid = .true.
                  return
                  
               endif
            case(-4)
               if( (flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (flIsReal(caElem(4)) .eqv. .true.) .and.
     &             (caElem(5) .eq. '          ') .and.
     &             (caElem(6) .eq. '          ') .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  flIsValid = .true.
                  return
                  
               endif
            case(-5)
               if( (flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (flIsReal(caElem(4)) .eqv. .true.) .and.
     &             (flIsReal(caElem(5)) .eqv. .true.) .and.
     &             (caElem(6) .eq. '          ') .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  flIsValid = .true.
                  return
                  
               endif
               
            case(-6)
               if( (flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (flIsReal(caElem(4)) .eqv. .true.) .and.
     &             (flIsReal(caElem(5)) .eqv. .true.) .and.
     &             (flIsReal(caElem(6)) .eqv. .true.) .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  flIsValid = .true.
                  return
                  
               endif
               
            case(-7)
               if( (flIsReal(caElem(1)) .eqv. .true.) .and.
     &             (flIsReal(caElem(2)) .eqv. .true.) .and.
     &             (flIsReal(caElem(3)) .eqv. .true.) .and.
     &             (flIsReal(caElem(4)) .eqv. .true.) .and.
     &             (flIsReal(caElem(5)) .eqv. .true.) .and.
     &             (flIsReal(caElem(6)) .eqv. .true.) .and.
     &             (flIsReal(caElem(7)) .eqv. .true.) .and.
     &             (caElem(8) .eq. '          ') ) then
               
                  flIsValid = .true.
                  return
                  
               endif
            case(0)
               if( (caElem(1) .eq. ' 000000000') .and.
     &             (caElem(2) .eq. ' 000000000') .and.
     &             (caElem(3) .eq. ' 000000000') .and.
     &             (caElem(4) .eq. ' 000000000') .and.
     &             (caElem(5) .eq. ' 000000000') .and.
     &             (caElem(6) .eq. ' 000000000') .and.
     &             (flIsInteger( caElem(7) ) .eqv. .true.) .and.
     &             (flIsInteger( caElem(8) ) .eqv. .true.) ) then

                  flIsValid = .true.
                  return
                  
               endif
             
            case(1)
               if( (caElem(1) .eq. ' 111111111') .and.
     &             (caElem(4) .eq. '          ') .and.
     &             (caElem(5) .eq. '          ') .and.
     &             (caElem(6) .eq. '          ') .and.
     &             (caElem(7) .eq. '          ') .and.
     &             (caElem(8) .eq. '          ') .and.
     &             (flIsInteger( caElem(2) ) .eqv. .true.) .and.
     &             (flIsInteger( caElem(3) ) .eqv. .true.) ) then            

                  flIsValid = .true.
                  return
             
               endif
             
            case(8)
               if( (caElem(1) .eq. ' 888888888') .and.
     &             (caElem(2) .eq. ' 888888888') .and.
     &             (caElem(3) .eq. ' 888888888') .and.
     &             (caElem(4) .eq. ' 888888888') .and.
     &             (caElem(5) .eq. ' 888888888') .and.
     &             (flIsInteger( caElem(6) ) .eqv. .true.) .and.
     &             (flIsInteger( caElem(7) ) .eqv. .true.) .and.
     &             (flIsInteger( caElem(8) ) .eqv. .true.) ) then            

                  flIsValid = .true.
                  return
             
               endif     

            case(9)
               if( (caElem(1) .eq. ' 999999999') .and.
     &             (caElem(2) .eq. ' 999999999') .and.
     &             (caElem(3) .eq. ' 999999999') .and.
     &             (caElem(4) .eq. ' 999999999') .and.
     &             (caElem(5) .eq. ' 999999999') .and.
     &             (caElem(6) .eq. ' 999999999') .and.
     &             (flIsInteger( caElem(7) ) .eqv. .true.) .and.
     &             (flIsInteger( caElem(8) ) .eqv. .true.) ) then            

                  flIsValid = .true.
                  return
             
               endif    
            
            case default
            
         end select
         return
      end function
     
      
      
      subroutine TXSVALIDATE(pNF, plFlag, plQuiet)
         
!        Begin your code here...
         implicit none
         integer, intent(in) :: pNF
         logical, intent(in) :: plQuiet
         logical, intent(out) :: plFlag
         character(len=80) :: cArg
         character(len=80) :: cText
         integer :: ieof = 0
         integer :: nLinePerRow, jLine, nFULL, iTYPE
         character(len=10) :: cNCell, cNGroup, cNLayer, cLayerID,
     &                        cCellID, cLine, cPGroup
         integer :: nGroup, iRec, iLine, nPrincipalGroup,
     &              iCell, iGrp, iLay, nLayer, nWarning
         
         plFlag = .true.
   
         iLine = 0
         write(cLine,'(I10)') iLine
!        First we open the joined TXS file.
         rewind pNF
         goto 92
91       if(plQuiet) goto 70
         print*, ' FATAL ERROR. Fail to load TRIGA Cross ',
     &           'Section (TXS) file: ', trim(cArg)
70       continue
         plFlag = .false.
         return
92       continue
         if(plQuiet) goto 71
         print*
         print*, ' TAPE-READ: Validating the TXS tape.'
71       continue
!        Then we read the first line of the file.      
94       read(pNF, '(A80)', iostat=ieof, err=93) cText
         iLine = iLine + 1
         write(cLine,'(I10)') iLine
!        Check whether the tape is already reaches the end of file.
         if(ieof .ne. 0) then
            if(plQuiet) goto 72
            print*, ' TAPE ERROR at line ', adjustl(cLine//'.')
            print*, ' --- The tape content is empty.'
72          continue
            plFlag = .false.
            return
         endif

!        Search for the 8s header. If not found -> stop.
         if(flIsValid(8, cText) .eqv. .false.) then
            if(plQuiet) goto 73
            print*, ' TAPE ERROR at line ', adjustl(cLine//'.')
            print*, ' HEADER-8 of the tape is missing.'
73          continue
            plFlag = .false.
            return
         else
!           If 8s header is found, set the principal number of core cell per layer,
!           cNCell and the number of core cell layer, cNLayer.
            cNCell = cText(61:70)
            if(plQuiet) goto 74
            print*, ' TAPE-READ: The principal number of ',
     &              'core cells per layer is ',
     &               adjustl(cNCell // '.')
74          continue
            cNLayer = cText(71:80)
            if(plQuiet) goto 75
            print*, ' TAPE-READ: The number of core layers is ',
     &               adjustl(cNLayer // '.')
75          continue
            read(cText(51:60),'(I10)') nGroup
            nPrincipalGroup = nGroup
            write(cPGroup,'(I10)') nGroup
            if(plQuiet) goto 76
            print'(A,A10)', '  TAPE-READ: The principal number of ' //
     &      'macroscopic neutron energy groups is ',
     &       adjustl(cPGroup // '.')
76          continue
!           We now calculate the number of line per cross section table row.
            if(floor(real(nPrincipalGroup+6)/real(8)) .lt. 
     &         (real(nPrincipalGroup+6)/real(8))) then
               nLinePerRow = int(floor(real(nPrincipalGroup+6)/real(8)))
     &                       + 1
            else
               nLinePerRow = int(floor(real(nPrincipalGroup+6)/real(8)))
            endif
            nFULL = int(floor(real(nPrincipalGroup+6)/8.0))
            iTYPE = - ((nPrincipalGroup+6) - nFULL * 8)
         endif
         
         
         iLay = 0
96       read(pNF, '(A80)', iostat=ieof, err=93) cText
         iLine = iLine + 1
         write(cLine,'(I10)') iLine
         if(ieof .ne. 0) then
            read(cNLayer,'(I10)') nLayer
            if (iLay .eq. nLayer) then
               if(plQuiet) goto 77
               print*, ' TAPE-READ: Validation complete.'
77             continue
               plFlag = .true.
               return
            else
               if(plQuiet) goto 78
               print*, ' TAPE ERROR at line ', adjustl(cLine//'.')
               print*, ' --- Part of the tape is missing while',
     &         ' reading HEADER-1'
               print*, '     Layer ID: ', cLayerID
               print*, '     Cell    : ', cCellID
78             continue
               plFlag = .false.
               return
            endif
         endif

102      iLay = iLay + 1         
         if(flIsValid(0, cText) .eqv. .false.) then
            iLay = iLay - 1
            if ( cText(1:10) .eq. ' 000000000') then
               if(plQuiet) goto 79
               print*, ' TAPE ERROR at line ', adjustl(cLine//'.')
               print*, ' --- Invalid HEADER-0 format is detected',
     &                 ' at line ', adjustl(cLine//'.')
               print*, '     Last Layer ID : ', cLayerID
               print*, '     Last Cell  ID : ', cCellID
79             continue
               plFlag = .false.
               return
            endif
            goto 96
         else
            if (cText(61:70) .ne. cNCell) then
               if(plQuiet) goto 80
               print*, ' TAPE ERROR at line ', adjustl(cLine//'.')
               print*, ' --- Cell/layer count mismatch.'
80             continue
               plFlag = .false.
               return
            else
!           We obtain the layer ID that we are processing.
               cLayerId = cText(71:80)
101            iCell = 0
97             read(pNF, '(A80)', iostat=ieof, err=93) cText
               iLine = iLine + 1
               write(cLine,'(I10)') iLine
               if(ieof .ne. 0) then
                  if(plQuiet) goto 81
                  print*, ' TAPE ERROR at line ', adjustl(cLine//'.')
                  print*, ' --- Part ',
     &            'of the tape is missing after reading HEADER-1.'
                  print*, '      Last Layer ID : ', cLayerID
                  print*, '      Last Cell  ID : ', cCellID
81                continue
                  plFlag = .false.
                  return
               endif
               
!              CHECKING HEADER-1 VALIDITY. IF HEADER-1 IS FOUND WITH INCORRECT FORMAT,
!              THEN THE TAPE  IS  CONSIDERED  INVALID.  IF HEADER-1 IS  FOUND WITH THE 
!              CORRECT FORMAT, WE READ THE NEXT LINE AND SEARCH FOR A FOLLOWING  8-COL
!              DATA RECORD. IF TGE 8-COL RECORD IS NOT FOUND,  THE TAPE IS  CONSIDERED
!              NOT VALID - NO COMMENTS/BLACK SPACES ARE ALLOWED AFTER HEADER-1.

               if(flIsValid(1, cText) .eqv. .false.) then
!                 IF HEADER-1 IS NOT FOUND, SEARCH FOR FOOTER-9.
                  if(flIsValid(9, cText) .eqv. .false.) then
!                    IF ITS NOT A FOOTER-9, CHECK IF THE LINE IS AN INVALID HEADER-1.                    
                     if ( cText(1:10) .eq. ' 111111111') then
!                       IF IT IS AN INVALID HEADER-1, THEN THE TAPE IS INVALID. END.
                        if(plQuiet) goto 82
                        print*, ' TAPE ERROR at line ',
     &                  adjustl(cLine//'.')
                        print*, ' --- Invalid HEADER-1 format is ',
     &                 'detected at line ', adjustl(cLine//'.')
                        print*, '     Last Layer ID : ', cLayerID
                        print*, '     Last Cell  ID : ', cCellID
82                      continue     
                        plFlag = .false.
                        return
                     endif
!                    IF IT IS NEITHER A HEADER-1, OR AN INVALID HEADER-1 OR A FOOTER-9,
!                    CONTINUE READ THE NEXT LINE SEARCH FOR HEADER-1 (GOTO LABEL 97)
!                    WE IGNORE AND PROCEED BECAUSE COMMENTS ARE ALLOWED AFTER A HEADER-0.
                     goto 97
                  else
!                    IF FOOTER-9 IS FOUND, WE CHECK WHETHER THE FOOTER MATCHES WITH THE
!                    CURRENT REACTOR CORE LAYER HEADER-0. IF NOT, THEN THE TAPE IS NOT
!                    VALID. IF YES, WE READ THE NEXT LINE AND SEARCH FOR A NEW HEADER-0
!                    FOR THE NEXT REACTOR CORE LAYER (GOTO LABEL 96).
                     if ((cText(61:70) .eq. cNCell) .and.
     &                   (cText(71:80) .eq. cLayerId)) then
                        goto 96
                     else 
                        if(plQuiet) goto 83
                        print*, ' TAPE ERROR at line ', 
     &                     adjustl(cLine//'.')
                        print*, ' --- Cell/layer ID mismatch.'
                        print*, '      Last Layer ID : ', cLayerID
                        print*, '      Last Cell  ID : ', cCellID
83                      continue                        
                        plFlag = .false.
                        return
                     endif
                  endif
               else
               
!                 IF HEADER-1 IS FOUND AND VALID, WE EXTRACT THE CELL ID AND THE NEUTRON 
!                 GROUP COUNT FROM HEADER-1.
                  iCell = iCell + 1
                  cCellId = cText(11:20)
                  cNGroup = cText(21:30)
                  read(cText(21:30),'(I10)') nGroup
!                 IF THE HEADER-1 IS NOT THE FIRST HEADER-1 FROM THE FIRST REACTOR CORE LAYER,
!                 WE CHECK IF THE NEUTRON GROUP DEFINED BY THE HEADER-1 IS CONSISTENT WITH
!                 THE PRINCIPAL NUMBER OF NEUTRON ENERGY GROUPS. THE PRINCIPAL NUMBER WAS 
!                 OBTAINED FROM THE 6-TH ELEMENT OF HEADER-8 DEFINED AT THE FIRST LINE OF 
!                 THE TXS TAPE.
                  if ((iCell .ne. 1) .or. (iLay .ne. 1)) then
                     if (nPrincipalGroup .ne. nGroup) then
                        if(plQuiet) goto 84
                        print*, ' TAPE ERROR at line ', 
     &                     adjustl(cLine//'.')
                        print*, ' --- Neutron energy group count' //
     &                  ' mismatch.'
                        print*, '     Layer ID: ', cLayerID
                        print*, '     Cell    : ', cCellID
84                      continue                 
                        plFlag = .false.
                        return                     
                     endif
                  endif
                  iGrp = 0
!                 WE READ THE NEXT LINE. THE NEXT LINE AFTER HEADER-1 MUST BE A 8-COL DATA
!                 RECORD, WHICH IS A PART OF THE ROW OF A NEUTRON CROSS SECTION TABLE 
!                 THAT BELONGS TO THE LAST HEADER-1 READ. IF IT IS NOT A 8-COL DATA RECORD,
!                 THEN THE TAPE IS NOT VALID.
                  jLine = 0
98                read(pNF, '(A80)', iostat=ieof, err=93) cText
                  if(ieof .ne. 0) then
                     if(plQuiet) goto 85
                     print*, ' TAPE ERROR at line ', adjustl(cLine//'.')
                     print*, ' --- Part ',
     &               'of the tape is missing after a 8-column record.'
                     print*, '     Layer ID: ', cLayerID
                     print*, '     Cell    : ', cCellID
85                   continue                     
                     plFlag = .false.
                     return
                  endif
                  iLine = iLine + 1
                  jLine = jLine + 1
                  write(cLine,'(I10)') iLine
                  if((jLine .le. nFULL) .and. 
     &               (flIsValid(-8, cText) .eqv. .false.)) then
                     if(plQuiet) goto 86
                     print*, ' TAPE ERROR at line ', adjustl(cLine//'.')
                     print*, ' --- Part ',
     &               'of the tape is missing after a 8-column record.'
                     print*, '     Layer ID: ', cLayerID
                     print*, '     Cell    : ', cCellID
86                   continue                     
                     plFlag = .false.
                     return
                  else
!                    IF THE LINE AFTER THE HEADER-1 IS A 8-COL DATA, THEN WE PROCEED TO
!                    READ THE NEXT LINE. THE NEXT LINE AFTER A 8-COL DATA MUST BE A 2-COL
!                    DATA RECORD. IF A N-COL DATA RECORD IS NOT FOUND, THE THE TAPE IS 
!                    NOT VALID.
!                    NOTE: COMMENTS ARE NOT ALLOWED AFTER HEADER-1, 8-COL DATA RECORD AND
!                          2-COL DATA RECORD. HOWEVER, COMMENTS ARE ALLOWED AFTER THE 2-COL
!                          DATA RECORD OF THE LAST ROWS OF THE NEUTRON CROSS SECTION TABLE
!                          OF THE LAST TABLE-1 READ.
100                  iGrp = iGrp + 1
99                   read(pNF, '(A80)', iostat=ieof, err=93) cText   
                     if(ieof .ne. 0) then
                        plFlag = .false.
                        return
                     endif
                     iLine = iLine + 1
                     write(cLine,'(I10)') iLine
                     if(flIsValid(iTYPE, cText) .eqv. .false.) then
                     if(plQuiet) goto 87
                     print*, ' TAPE ERROR at line ', adjustl(cLine//'.')
                     print*, ' --- Part ',
     &               'of the tape is missing after a ',
     &               'N-column record.'
                     print*, '     Layer ID: ', cLayerID
                     print*, '     Cell    : ', cCellID
87                   continue                     
                     plFlag = .false.
                     return
                     else 
                        if(iGrp .eq. nGroup) then
                           goto 97
                        endif
                        goto 98
                     endif
                  endif
               endif
            endif
         endif
93       continue
95       continue
         plFlag = .false.
         return
      end subroutine
      
      
      end module