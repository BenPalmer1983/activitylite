PROGRAM activitylite



!---------------------------------------------
!  Calculated the number of atoms of each isotope in a decay chain
!  Follows only one branch, gives total activity of this branch
!  Production rate of parent isotope
!  Input starting number of atoms, half life (s) and branching factor for each isotope
!
!  Example input file (chain used by :http://wordpress.mrreid.org/)
!
!  0.0D0    1.0D0    1.0D1   10   1           ! parent production rate, tStart, tEnd, data points, calcType
!  Rd-209   3.96D0       6.02D6               ! Parent  (1)  t1/2  N(t) 
!  Po-215   2D-3         0.0D0   1.0D0        ! Child A (2)  t1/2  N(t)  B(1,2)
!  Pb-211   2166.0D0     0.0D0   1.0D0        ! Child B (3)  t1/2  N(t)  B(1,2)  
!  Bi-211   128.0D0      0.0D0   1.0D0        ! Child C (4)  t1/2  N(t)  B(1,2) 
!  Po-211   0.516        0.0D0   1.0D0        ! Child D (5)  t1/2  N(t)  B(1,2) 
!  Pb-207   -1           0.0D0   1.0D0        ! Child E (6)  t1/2  N(t)  B(1,2)
!
!
!  Compile using the make.sh script
!
!  Example usage:
!  
!  ./activitylite.x chain.in
!
!
!






! Force declaration of all variables
  Implicit None
! Data types  
  Integer, Parameter :: SingleReal = Selected_Real_Kind(6,37)         ! single real, 6 decimal precision, exponent range 37
  Integer, Parameter :: DoubleReal = Selected_Real_Kind(15,307)       ! double real, 15 decimal precision, exponent range 307
  Integer, Parameter :: QuadrupoleReal = Selected_Real_Kind(33,4931)  ! quadrupole real
  Integer, Parameter :: TinyInteger = Selected_Int_Kind(1)            ! tiny integer    1 byte
  Integer, Parameter :: SmallInteger = Selected_Int_Kind(4)           ! small integer    4 bytes -2E31 to 2E31-1
  Integer, Parameter :: StandardInteger = Selected_Int_Kind(8)        ! standard integer 8 bytes -2E63 to 2E63-1
  Integer, Parameter :: LongInteger = Selected_Int_Kind(12)           ! long integer
  Integer, Parameter :: VeryLongInteger = Selected_Int_Kind(32)       ! very long integer
! Constants
  Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000
  Real(kind=QuadrupoleReal), Parameter :: lnTwoQ = 0.6931471805599453094172321214581765680755D0

  Call runCalc() 

  Contains
 
  Subroutine runCalc()
! force declaration of all variables
    Implicit None
! declare variables
    Real(kind=DoubleReal), Parameter :: lnTwo = 0.693147180559945D0
    Character(len=255) :: inputFile
    Character(len=255) :: fileRow, rowTemp
    Character(len=64) :: bufferA, bufferB, bufferC, bufferD, bufferE
    Character(len=16), Dimension(1:100) :: labelArray
    Integer(kind=StandardInteger) :: i, j, ios, isotopeCount, tIncs, calcType
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChange
    Real(kind=DoubleReal) :: time, w, tStart, tEnd, activity
! Get file name
    call get_command_argument(1,inputFile)
! count isotopes    
    Open(UNIT=1,FILE=inputFile)
    isotopeCount = -1
    Do i=1,maxFileRows
      Read(1,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then  ! break out
        EXIT
      End If
      rowTemp = RemoveSpaces(fileRow)
      If(rowTemp(1:1).ne." ")Then
        isotopeCount = isotopeCount + 1
      End If
    End Do
    Close(1)
! Allocate arrays
    Allocate(decayDataArray(1:isotopeCount,1:6))    
    Allocate(isotopeChange(1:isotopeCount,1:12))
! Read input   
    Print *,"Input Decay Chain"
    Print *,"--------------------------------------------------------------------------------------"
    Open(UNIT=1,FILE=inputFile)
    isotopeCount = 0
    Do i=1,maxFileRows
      Read(1,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then  ! break out
        EXIT
      End If
      If(i.eq.1)Then
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE
        Read(bufferA,*) w
        Read(bufferB,*) tStart
        Read(bufferC,*) tEnd
        Read(bufferD,*) tIncs
        Read(bufferE,*) calcType
        Print *,"Parent production rate:    ",w
        Print *,"Print activity start time: ",tStart
        Print *,"Print activity end time:   ",tEnd
        Print *,"Data points:               ",tIncs        
      Else
        rowTemp = RemoveSpaces(fileRow)
        If(rowTemp(1:1).ne." ")Then
          isotopeCount = isotopeCount + 1
          If(isotopeCount.eq.1)Then
            Read(fileRow,*) bufferA, bufferB, bufferC
            Read(bufferB,*) decayDataArray(isotopeCount,3)  ! Half life
            Read(bufferC,*) decayDataArray(isotopeCount,2)  ! Start atoms
            decayDataArray(isotopeCount,4) = 1.0D0 
            labelArray(isotopeCount) = bufferA(1:16)
            print *,labelArray(isotopeCount),decayDataArray(isotopeCount,3),&
            decayDataArray(isotopeCount,2),decayDataArray(isotopeCount,4)
          Else 
            Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
            Read(bufferB,*) decayDataArray(isotopeCount,3)  ! Half life
            Read(bufferC,*) decayDataArray(isotopeCount,2)  ! Start atoms
            Read(bufferD,*) decayDataArray(isotopeCount,4)  ! Branching factor    
            labelArray(isotopeCount) = bufferA(1:16)
            print *,labelArray(isotopeCount),decayDataArray(isotopeCount,3),&
            decayDataArray(isotopeCount,2),decayDataArray(isotopeCount,4)
          End If
        End If
      End If  
    End Do
    Close(1) 
    Print *,"--------------------------------------------------------------------------------------"
! string format
    bufferB = BlankString(bufferB)
    write(bufferB,"(I4)") isotopeCount
    bufferB = RemoveSpaces(bufferB)
    bufferB = "("//trim(bufferB)//"(A16,A1))"
! ----    
    bufferC = BlankString(bufferC)
    write(bufferC,"(I4)") isotopeCount
    bufferC = RemoveSpaces(bufferC)
    bufferC = "("//trim(bufferC)//"(E14.6,A3))"
! print out    
    rowTemp = BlankString(rowTemp)
    write(rowTemp, trim(bufferB)) (labelArray(j)," ", j=1,isotopeCount)
    rowTemp = "Time/s         "//trim(rowTemp)//"                       "
    rowTemp = rowTemp(1:(15+isotopeCount*17))//"Activity/Bq"
    print *,rowTemp
    Do i=1,tIncs
      time = tStart+(i-1)*((tEnd-tStart)/(tIncs-1))
      isotopeChange = CalcIsotopeAmount(100.0D0,decayDataArray,time)   
      If(calcType.eq.1)Then  ! Print out analytic result (and numeric where analytic stops)
        activity = 0.0D0
        Do j=1,isotopeCount
          If(isotopeChange(j,7).gt.0)Then
            activity = activity + isotopeChange(j,4)*isotopeChange(j,8)
          End If
        End Do
        bufferA = BlankString(bufferA)
        rowTemp = BlankString(rowTemp)
        write(bufferA,"(E14.6)") time
        write(bufferB,"(E14.6)") activity
        write(rowTemp, trim(bufferC)) (isotopeChange(j,4),"   ", j=1,isotopeCount)
        print *,bufferA(1:15),rowTemp(1:(isotopeCount*17)),bufferB
      End If
      If(calcType.eq.2)Then  ! Print out numeric only
        activity = 0.0D0
        Do j=1,isotopeCount
          If(isotopeChange(j,7).gt.0)Then
            activity = activity + isotopeChange(j,12)*isotopeChange(j,8)
          End If
        End Do
        bufferA = BlankString(bufferA)
        rowTemp = BlankString(rowTemp)
        write(bufferA,"(E14.6)") time
        write(bufferB,"(E14.6)") activity
        write(rowTemp, trim(bufferC)) (isotopeChange(j,12),"   ", j=1,isotopeCount)
        print *,bufferA(1:15),rowTemp(1:(isotopeCount*17)),bufferB
      End If
      
    End Do
  End Subroutine runCalc
  
  
  
! ------------------------------------------------------------------------!
! Decay Functions
! ------------------------------------------------------------------------!

  Function CalcIsotopeAmount(w,decayDataArray,t) RESULT (isotopeChange)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j,decaySteps,decayStepCounter, noChanges
    Real(kind=DoubleReal) :: halfLifeChange, randNumber, w, t
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChange
    Real(kind=DoubleReal) :: stableLimit
! Quadrupole Reals    
    Real(kind=QuadrupoleReal) :: resultQ, resultGS, tQ, tempQ
    Real(kind=QuadrupoleReal), Dimension(1:20) :: L ! Lambda
    Real(kind=QuadrupoleReal), Dimension(1:20) :: N ! Starting number of atoms
    Real(kind=QuadrupoleReal), Dimension(1:20) :: E ! Exp
    Real(kind=QuadrupoleReal), Dimension(1:19) :: B ! Exp
! -------------------------------------------------    
! decaySteps really means decay isotopes in chain (steps = decaySteps-1)
! ------------------------------------------------- 
! Input decay chain array:   
! decayDataArray(i,1) !Tally key
! decayDataArray(i,2) !No. Atoms
! decayDataArray(i,3) !Half life
! decayDataArray(i,4) !branching factor
! decayDataArray(i,5) !isotope Z
! decayDataArray(i,6) !isotope A   
!-------------------------------------------------    
! Output decay chain array:   
! isotopeChange(i,1)    !Tally key
! isotopeChange(i,2)    !Change in isotope amount
! isotopeChange(i,3)    !Start amount
! isotopeChange(i,4)    !End amount
! isotopeChange(i,5)    !Isotope Z
! isotopeChange(i,6)    !Isotope A
! isotopeChange(i,7)    !T1/2
! isotopeChange(i,8)    !Decay constant
! isotopeChange(i,9)    !Branching factor
! isotopeChange(i,10)   !Parent production rate  
! isotopeChange(i,11)   !Time    
! isotopeChange(i,12)   !GS End    
! -------------------------------------------------   
! Init variables
   tQ = t
! -------------------------------------------------   
! Alter decay chain
! -------------------------------------------------  
! - If dTime * decay constant lt 1.0D-14 then assume stable for purposes of simulation
    decayStepCounter = 0
    Do i=1,size(decayDataArray,1)
      stableLimit = (log(2.0D0)/decayDataArray(i,3))*t
      decayStepCounter = decayStepCounter + 1
      If(stableLimit.lt.1.0D-14)Then
        decayDataArray(i,3) = -1    !set as stable
        Exit
      End If
    End Do    
! Resize array
    decayDataArray = ArraySize2DDouble(decayDataArray,decayStepCounter)
! -------------------------------------------------   
! Set stable isotope decay constant very small to avoid infinity error   
! -------------------------------------------------   
    Do i=1,size(decayDataArray,1)      
      If(decayDataArray(i,3).eq.(-1))Then
        decayDataArray(i,3) = 1.0D100
      End If
    End Do   
! -------------------------------------------------   
! Break same decay constants by ~1E-3% to avoid singularities
! -------------------------------------------------   
    noChanges = 0
    Do While(noChanges.eq.0)
      noChanges = 1
      Do i=1,size(decayDataArray,1)      
        Do j=1,size(decayDataArray,1)      
          If(i.ne.j)Then
            If(decayDataArray(i,3).eq.decayDataArray(j,3))Then
              Call RANDOM_NUMBER(randNumber)
              halfLifeChange = 0.1D0+randNumber*0.9D0
              halfLifeChange = decayDataArray(i,3)*1D-5*halfLifeChange
              decayDataArray(i,3) = decayDataArray(i,3)+halfLifeChange
              decayDataArray(j,3) = decayDataArray(j,3)-halfLifeChange
              noChanges = 0
            End If
          End If
        End Do 
      End Do   
    End Do
! set decay steps/isotopes
    decaySteps = size(decayDataArray,1)
! allocate isotopeChange array
    Allocate(isotopeChange(1:decaySteps,1:12))
! Fill with starting data
    Do i=1,decaySteps
      isotopeChange(i,1) = decayDataArray(i,1)
      isotopeChange(i,2) = 0.0D0          !default no change
      isotopeChange(i,3) = decayDataArray(i,2)
      isotopeChange(i,4) = decayDataArray(i,2)    !default no change
      isotopeChange(i,5) = decayDataArray(i,5)
      isotopeChange(i,6) = decayDataArray(i,6)
      isotopeChange(i,7) = decayDataArray(i,3)
      isotopeChange(i,8) = log(2.0D0)/decayDataArray(i,3)
      isotopeChange(i,9) = decayDataArray(i,4)
      isotopeChange(i,10) = w
      isotopeChange(i,11) = t
      isotopeChange(i,12) = 0.0D0          !default no change
    End Do
! Store lambda starting atom number data
    Do i=1,decaySteps
      If(decayDataArray(i,3).gt.9.9D99)Then
        L(i) = 0.0D0
      Else
        L(i) = lnTwoQ/isotopeChange(i,7)
      End If
      N(i) = isotopeChange(i,3)
      tempQ = -1.0D0*L(i)*tQ
      E(i) = exp(tempQ)
      B(i) = decayDataArray(i,4)
    End Do
!    
! nP -> nA -> nB -> nC -> nD ...
!    
!Set starting variables 
    If(decaySteps.ge.1)Then    
! calc nP
      resultQ = (w/L(1))*(1-E(1))+N(1)*E(1)
      resultGS = CalcIsotopeAmountGS(tQ,1,isotopeChange)
      If(ISNAN(resultQ))Then ! solve numerically
        isotopeChange(1,4) = dble(resultGS)
        isotopeChange(1,12) = dble(resultGS)  
      Else  
        isotopeChange(1,4) = dble(resultQ)  
        isotopeChange(1,12) = dble(resultGS)    
      End If
    End If  
    If(decaySteps.ge.2)Then      
! calc nA     
      resultQ = B(2)*L(1)*w*(1.0D0/(L(1)*L(2))+E(1)/(L(1)*(L(1)-L(2)))-&
                E(2)/(L(2)*(L(1)-L(2))))+&
                B(2)*L(1)*N(1)*(E(1)/(L(2)-L(1))+E(2)/(L(1)-L(2)))+&
                N(2)*E(2)            
      resultGS = CalcIsotopeAmountGS(tQ,2,isotopeChange)
      If(ISNAN(resultQ))Then ! solve numerically
        isotopeChange(2,4) = dble(resultGS) 
        isotopeChange(2,12) = dble(resultGS)  
      Else  
        isotopeChange(2,4) = dble(resultQ)  
        isotopeChange(2,12) = dble(resultGS)      
      End If
    End If  
    If(decaySteps.ge.3)Then   
! child B terms
      resultQ = &
      w*B(2)*B(3)*L(1)*L(2)*&
      (1.0D0/(L(1)*L(2)*L(3))-&
      E(1)/(L(1)*(L(1)-L(2))*(L(1)-L(3)))+&
      E(2)/(L(2)*(L(1)-L(2))*(L(2)-L(3)))+&
      E(3)/(L(3)*(L(1)-L(3))*(L(3)-L(2))))+&
      B(2)*B(3)*L(1)*L(2)*N(1)*&
      (E(1)/((L(1)-L(2))*(L(1)-L(3)))-&
      E(2)/((L(1)-L(2))*(L(2)-L(3)))-&
      E(3)/((L(1)-L(3))*(L(3)-L(2))))+&
      B(3)*L(2)*N(2)*&          
      (E(1)/(L(2)-L(1))+E(2)/(L(1)-L(2)))+&
      N(3)*E(3)
      resultGS = CalcIsotopeAmountGS(tQ,3,isotopeChange)
      If(ISNAN(resultQ))Then ! solve numerically
        isotopeChange(3,4) = dble(resultGS)   
        isotopeChange(3,12) = dble(resultGS)  
      Else  
        isotopeChange(3,4) = dble(resultQ) 
        isotopeChange(3,12) = dble(resultGS)       
      End If
    End If    
! Numeric inverse laplace for remainder    
    If(decaySteps.ge.4)Then  
      Do i=4,decaySteps
        resultGS = CalcIsotopeAmountGS(tQ,i,isotopeChange)
        isotopeChange(i,4) = dble(resultGS)  
        isotopeChange(i,12) = dble(resultGS)  
      End Do
    End If
! Adjust the isotope values    
    Do i=1,decaySteps 
      If(isotopeChange(i,4).lt.0.0D0)Then
        isotopeChange(i,4) = 0.0D0
      End If
    End Do 
! Store changes in isotope amounts
    Do i=1,size(isotopeChange,1)
      isotopeChange(i,2) = isotopeChange(i,4) - isotopeChange(i,3)
    End Do
  End Function CalcIsotopeAmount
  
  
  Function CalcIsotopeAmountGS(t,isotopeStep,isotopeChangeIn) RESULT (output)
! Force declaration of all variables
    Implicit None
! Declare variables  
    Integer(kind=StandardInteger) :: i, isotopeStep, M, k
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChangeIn
    Real(kind=QuadrupoleReal), Dimension(1:50) :: weightingQ    
    Real(kind=QuadrupoleReal), Dimension(1:20) :: L ! Lambda
    Real(kind=QuadrupoleReal), Dimension(1:20) :: N ! Starting number of atoms
    Real(kind=QuadrupoleReal), Dimension(1:20) :: B ! Starting number of atoms
    Real(kind=QuadrupoleReal) :: kQ, w, t, ft, s, FS, output 
!-------------------------------------------------    
! Output decay chain array:   
! isotopeChange(i,1)    !Tally key
! isotopeChange(i,2)    !Change in isotope amount
! isotopeChange(i,3)    !Start amount
! isotopeChange(i,4)    !End amount
! isotopeChange(i,5)    !Isotope Z
! isotopeChange(i,6)    !Isotope A
! isotopeChange(i,7)    !T1/2
! isotopeChange(i,8)    !Decay constant
! isotopeChange(i,9)    !Branching factor
! isotopeChange(i,10)   !Parent production rate   
! ------------------------------------------------- 
! Init variables
    M = 8    
    weightingQ = GaverStehfestWeightingQ(M,weightingQ)
    w = isotopeChangeIn(1,10)
    output = 0.0D0
! Adjust the isotope values    
    Do i=1,isotopeStep 
      If(isotopeChangeIn(i,4).lt.0.0D0)Then
        isotopeChangeIn(i,4) = 0.0D0
      End If
    End Do 
! Store lambda starting atom number data
    Do i=1,isotopeStep
      L(i) = lnTwoQ/isotopeChangeIn(i,7)
      N(i) = isotopeChangeIn(i,3)
      If(i.eq.1)Then
        B(i) = 1.0D0
      Else
        B(i) = isotopeChangeIn(i,9)
      End If
    End Do
! Perform calculation
    ft = 0.0D0
    Do k=1,2*M
      kQ = 1.0D0 * k
      s = (kQ*lnTwoQ)/t
! -----------------------  
      FS = (1.0D0/(s+L(1)))*(w/s+N(1))
      Do i=2,isotopeStep
        FS = (1.0D0/(s+L(i)))*(B(i)*L(i-1)*FS+N(2))
      End Do
      !FS = (1.0D0/(s+L(1)))*(w/s+N(1))
! -----------------------  
      ft = ft + weightingQ(k)*FS   
    End Do
    ft = (lnTwoQ/t)*ft 
    output = Dble(ft)
    !isotopeChangeOut(isotopeStep,4) = Dble(ft)
  End Function CalcIsotopeAmountGS
  
  
  
  
  Function GaverStehfestWeighting(N, weightingIn) RESULT (weighting)
! Force declaration of all variables
    Implicit None
! Declare variables      
    Integer(kind=StandardInteger) :: N
    Integer(kind=StandardInteger) :: j, k, jStart, jEnd
    Real(kind=DoubleReal) :: factor, wSum
    !Real(kind=DoubleReal), Dimension(1:2*N) :: weighting
    Real(kind=DoubleReal), Dimension(:) :: weightingIn
    Real(kind=DoubleReal), Dimension(1:size(weightingIn)) :: weighting
! Init array
    weighting = 0.0D0
! k loop
    Do k=1,2*N
      factor = (-1)**(k+N)/(1.0D0*FactorialDP(N))
      jStart = Floor((k+1)/2.0D0)
      jEnd = min(k,N)
      wSum = 0.0D0
! j loop
      Do j=jStart,jEnd
        wSum = wSum + 1.0D0*(j**(N+1))*BinomialCoefficientDP(N,j)*&
        BinomialCoefficientDP(2*j,j)*BinomialCoefficientDP(j,k-j)
      End Do
      weighting(k) = factor*wSum
    End Do
  End Function GaverStehfestWeighting  
  
  Function GaverStehfestWeightingQ(N, weightingIn) RESULT (weighting)
! Force declaration of all variables
    Implicit None
! Declare variables      
    Integer(kind=StandardInteger) :: N
    Integer(kind=StandardInteger) :: j, k, jStart, jEnd
    Real(kind=QuadrupoleReal) :: factor, wSum
    Real(kind=QuadrupoleReal), Dimension(:) :: weightingIn
    Real(kind=QuadrupoleReal), Dimension(1:size(weightingIn)) :: weighting
! Init array
    weighting = 0.0D0
! k loop
    Do k=1,2*N
      factor = (-1)**(k+N)/(1.0D0*FactorialQ(N))
      jStart = Floor((k+1)/2.0D0)
      jEnd = min(k,N)
      wSum = 0.0D0
! j loop
      Do j=jStart,jEnd
        wSum = wSum + 1.0D0*(j**(N+1))*BinomialCoefficientQ(N,j)*&
        BinomialCoefficientQ(2*j,j)*BinomialCoefficientQ(j,k-j)
      End Do
      weighting(k) = factor*wSum
    End Do
  End Function GaverStehfestWeightingQ
   
  Function Factorial(input) RESULT (output)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,input
    Integer(kind=StandardInteger) :: output
! calculate factorial
    output = 1
    Do i=1,input
      output = i * output
    End Do
  End Function Factorial

  Function BinomialCoefficient(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: c,n,k
! calculate factorial
    c = Factorial(n)/(Factorial(n-k)*Factorial(k))
  End Function BinomialCoefficient
  
  Function FactorialDP(input) RESULT (output)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,input
    Integer(kind=VeryLongInteger) :: tempInt
    Real(kind=DoubleReal) :: output
! calculate factorial
    tempInt = 1
    Do i=1,input
      tempInt = i * tempInt
    End Do
    output = 1.0D0*tempInt
  End Function FactorialDP

  Function BinomialCoefficientDP(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: n,k
    Real(kind=DoubleReal) :: c, nDP, nkDP, kDP
! calculate factorial
    nDP = FactorialDP(n)
    nkDP = Factorial(n-k)
    kDP = FactorialDP(k)
    c = 1.0D0*nDP/(nkDP*kDP)
  End Function BinomialCoefficientDP 
  
  Function FactorialQ(input) RESULT (output)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,input
    Integer(kind=VeryLongInteger) :: tempInt
    Real(kind=QuadrupoleReal) :: tempQ
    Real(kind=QuadrupoleReal) :: output
! calculate factorial
    tempInt = 1
    tempQ = 1
    Do i=1,input
      If(i.le.33)Then
        tempInt = i * tempInt
        tempQ = 1.0D0*tempInt
      End If  
      If(i.eq.34)Then
        tempQ = 1.0D0*i*tempInt 
      End If   
      If(i.ge.35)Then
        tempQ = 1.0D0*i*tempQ 
      End If       
    End Do
    output = tempQ
  End Function FactorialQ

  Function BinomialCoefficientQ(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: n,k
    Real(kind=QuadrupoleReal) :: c, nDP, nkDP, kDP
! calculate factorial
    nDP = FactorialDP(n)
    nkDP = Factorial(n-k)
    kDP = FactorialDP(k)
    c = 1.0D0*nDP/(nkDP*kDP)
  End Function BinomialCoefficientQ
  
  Function ArraySize2DDouble (inputArray,arraySizeHeight,arraySizeWidthIn) &
    RESULT (outputArray)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i, j
    Integer(kind=StandardInteger) :: arraySizeHeight
    Integer(kind=StandardInteger), optional :: arraySizeWidthIn
    Integer(kind=StandardInteger) :: arraySizeWidth
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: outputArray
! catch optional width
    If(present(arraySizeWidthIn))Then
      arraySizeWidth = arraySizeWidthIn
    Else
      arraySizeWidth = size(inputArray,2)
    End If
! Allocate output array
    Allocate(outputArray(1:arraySizeHeight,1:arraySizeWidth))
! transfer data
    Do i=1,arraySizeHeight
      Do j=1,arraySizeWidth
        If(i.le.size(inputArray,1).and.j.le.size(inputArray,2))Then
          outputArray(i,j) = inputArray(i,j)
        Else
          outputArray(i,j) = 0.0D0
        End If
      End Do
    End Do
  End Function ArraySize2DDouble
  

    Function BlankString (input) RESULT (output)
      Character(*), INTENT(IN) :: input
      Character(Len(input)) :: output
      Integer(kind=StandardInteger) :: i
      Do i=1,Len(input)
        output(i:i) = " "
      End Do
    End Function BlankString
  
    Function RemoveSpaces (input) RESULT (output)
      CHARACTER(*), INTENT(IN) :: input
      CHARACTER(LEN(input)) :: outputTemp
      CHARACTER(LEN(input)) :: output
! -- Local variables
      INTEGER :: i,j
! -- Copy input string
      outputTemp = input
! Blank output
      Do i = 1, LEN( outputTemp )
        output( i:i ) = " "
      End Do
! transfer outputtemp to output without spaces
      j = 0
      Do i = 1, LEN( outputTemp )
        If(outputTemp( i:i ).ne." ")Then
          j = j + 1
          output( j:j ) = outputTemp( i:i )
        End If
      End Do
    End Function RemoveSpaces

End Program activitylite