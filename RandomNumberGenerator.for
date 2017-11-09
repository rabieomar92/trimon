!     
!     Module RandomNumberGenerator. Author M. R. Omar, October 2017. Revision
!     171102-1.
!     (c) Universiti Sains Malaysia
!     (c) Agensi Nuklear Malaysia
!
!
!     Methods:
!
!     (1) SUBROUTINE INITRND(PNHIST, PNSTRIDE, PISEED0)
!         - THIS SUBROUTINE ALLOWS THE PROGRAMMER TO INITIALIZE THE RANDOM NUM-
!           BER GENERATOR.  PNHIST IS THE NUMBER OF PARTICLE HISTORY,  PNSTRIDE 
!           IS THE STRIDE  LENGTH  (PER PARTICLE HISTORY)  AND  PISEED0  IS THE 
!           STARTING SEED.
!
!     (2) REAL FUNCTION RND(PIHIST)
!         - THIS FUNCTION RETURNS  THE NEXT  SEQUENCE OF THE  RANDOM NUMBER FOR
!           THE SPECIFIED HISTORY CHANNEL (PIHIST). 
!      
!     (3) INTEGER FUNCTION FUTURESEED(PNN, PISEED)
!         - THIS FUNCTION SKIPS PNN NUMBER  OF SEQUENCE  AND RETURNS (PNN+1)-TH 
!           RANDOM NUMBER SEQUENCE. REFER MANUAL FOR MORE INFORMATION.
!
!     Notes:
!
!     THE RANDOM NUMBER GENERATION PROCESS  IS A MARKOV  PROCESS,  WHERE RANDOM 
!     NUMBERS ARE GENERATED SEQUENTIALLY BASED ON THEIR PREVIOUS  RANDOM NUMBER
!     SEQUENCE. THUS,  WE CAN DELIVER  A  HIGH SPEED RANDOM NUMBER  CALCULATION
!     WITHOUT WASTING COMPUTER MEMORY.
!
!     GOOD LUCK. BY M. R. OMAR.

      module RandomNumberGenerator
         
         implicit none

!        NUMBER OF PARTICLE HISTORY CHANNELS.
         integer(8) :: N_HIST = 999999_8
!        STRIDE LENGTH
         integer(8), parameter :: N_STR = 152917_8
!        RANDOM NUMBER SEQUENCE ARRAY. 
         integer(8), allocatable :: iaSeed(:)
         integer(8), parameter :: I_MUL = 2806196910506780709_8
         integer(8), parameter :: I_ADD = 1_8
         integer(8), parameter :: N_BIT = 63_8
         integer(8), parameter :: I_MOD = ibset(0_8, N_BIT)
         integer(8), parameter :: I_MSK = not(I_MOD)
         real(8)   , parameter :: R_NOR = 2._8**(-N_BIT)
         
      contains

!     -------------------------------------------------------------------------
!     THIS FUNCTION GENERATES A PSEUDO-RANDOM NUMBER USING A LINEAR CONGRUENT
!     GENERATOR FOR THE SPECIFIED HISTORY CHANNEL.
!     -------------------------------------------------------------------------

      real(8) function Rnd(piHist)
         integer(8), intent(in) :: piHist
         iaSeed(piHist) = iand(I_MUL * iaSeed(piHist) + I_ADD, I_MSK)
         Rnd = iaSeed(piHist) * R_NOR
         return
      end function
      
!     -------------------------------------------------------------------------
!     THIS SUB SETS UP THE RANDOM NUMBER GENERATOR, DETERMINING THE SEED AND
!     THE VALUES OF LINEAR CONGRUENTIAL CONSTANTS.
!     -------------------------------------------------------------------------

      subroutine InitRnd(pnHist, pnStride, piSeed0)
         integer(8), intent(in) :: pnHist
         integer(8), intent(in) :: pnStride
         integer(8), intent(in) :: piSeed0
         integer :: i  
         
         N_HIST = pnHist
         allocate(iaSeed(N_HIST))
         
         do i=1, N_HIST, 1
            iaSeed(i) = FutureSeed(i*pnStride, piSeed0)
         enddo
         
         return
         
      end subroutine
      
      integer(8) function FutureSeed(pnN, piSeed)
         integer(8), intent(in) :: pnN
         integer(8), intent(in) :: piSeed
         integer(8) :: nSkip
         integer(8) :: iG
         integer(8) :: iC
         integer(8) :: iGNew
         integer(8) :: iCNew
         
         nSkip = pnN
         do while(nSkip .lt. 0_8)
            nSkip = nSkip + I_MOD
         enddo
         
         nSkip = iand(nSkip, I_MSK)
         
         iG = I_MUL
         iC = I_ADD
         iGNew = 1
         iCNew = 0
         do while(nSkip .gt. 0_8)
            if(btest(nSkip,0)) then
               iGNew = iand(iGNew*iG, I_MSK)
               iCNew = iand(iCNew*iG + iC, I_MSK)
            endif
            iC = iand((iG+1)*iC, I_MSK)
            iG = iand(iG*iG, I_MSK)
            nSkip = ishft(nSkip, -1)
         enddo
         
         FutureSeed = iand(iGNew*piSeed + iCNew, I_MSK)
         
         return
         
      end function
      
      end module
