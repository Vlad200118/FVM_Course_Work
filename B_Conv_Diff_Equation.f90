Subroutine Calc_Conv_Dif_T (NI, NJ, V, T, CellVolume, IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter, SCHEME)
    Implicit none
    integer :: NI, NJ, I, J
    integer :: Iter
    Real :: Re, Pr, SCHEME, CFL, VNM
    Real :: maxRes, Eps, dx,dy, dTimeConv, dTimeDiff
    Real :: TLeft, TRight
    Real :: nu, H, U, coefDiff
    Real,Dimension(NI-1,NJ-1):: CellVolume
    Real,Dimension(NI,NJ-1,2) ::IFaceCenter, IFaceVector
    Real,Dimension(NI-1,NJ,2) ::JFaceCenter, JFaceVector
    Real,Dimension(0:NI,0:NJ,2) ::CellCenter
    Real,Dimension(0:NI,0:NJ,2) :: V
    Real,Dimension(0:NI,0:NJ) :: Res, T

    Real,allocatable,dimension(:,:) ::  deltaTime
    Real,allocatable,dimension(:,:) :: DivTV
    Real,allocatable,dimension(:,:,:)::GradT
    Real,allocatable,dimension(:,:):: LaplacianT

    allocate(DivTV( 0:NI, 0:NJ), GradT( 0:NI, 0:NJ, 2), LaplacianT(0:NI, 0:NJ))
    allocate(deltaTime(0:NI, 0:NJ))

    TLeft = 0.0
    TRight = 1.0

    U = 1.0
    H = 1.0

    Re = 100
    Pr = 0.7
    maxRes = 1.0
    Eps = 1.E-7

    Iter = 0
    T(:,:) = (TLeft + TRight)/2.0

    CFL = 0.4
    VNM = 0.4
    nu = U * H / Re
    coefDiff = nu / Pr

    DO J =1, NJ-1
      DO I=1, NI-1
        dx = abs((IFaceCenter(I+1,J,1) - IFaceCenter(I,J,1)))
        dy = abs((JFaceCenter(I,J+1,2) - JFaceCenter(I,J,2)))
        dTimeConv = CFL * min (dx/Norm2(V(I,J,:)), dy/Norm2(V(I,J,:)))
        dTimeDiff = VNM * min(dx*dx, dy*dy)/2.0/coefDiff
        deltaTime(I,J) = dTimeConv*dTimeDiff/(dTimeConv+dTimeDiff)
      END DO
    END DO
    write(*,*)maxval(deltaTime)
!pause
    DO WHILE(maxRes .GE. Eps)

        Iter = Iter + 1

        T(0,0:NJ)=TLeft
        T(NI,0:NJ)=TRight
        T(0:NI,0)=T(0:NI,1)
        T(0:NI,NJ)=T(0:NI,NJ-1)

        Call B_CalcGradientGGI(NI,NJ,T,CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter, GradT)
        Call B_CalcLaplacian(NI,NJ,T,GradT,CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter, LaplacianT)
        Call B_CalcDivergancePV(NI,NJ,V,T,GradT,CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter,DivTV,SCHEME)

        Res(1:NI-1,1:NJ-1) = DivTV(1:NI-1,1:NJ-1) - LaplacianT(1:NI-1,1:NJ-1)/Re/Pr
        T(1:NI-1,1:NJ-1) = T(1:NI-1,1:NJ-1) - Res(1:NI-1,1:NJ-1) * deltaTime (1:NI-1,1:NJ-1)


        maxRes = maxval(abs(Res(1:NI-1, 1:NJ-1)))
        if (MOD(Iter, 500) .EQ. 0) then
            Write(*,*) Iter, maxRes
        end if
    END DO

    Deallocate(DivTV, GradT, LaplacianT, deltaTime)
End subroutine
