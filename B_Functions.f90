Function Pressure(X,Y)
    implicit none
    real :: Pressure, X, Y
    Pressure = x**1.0+y**1.0
End Function

Subroutine CalcInitV(I, J, NI, NJ, CellCenter, V)
    implicit none
    integer:: I, J, NI, NJ
    Real,Dimension(0:NI,0:NJ,2):: V
    Real,Dimension(0:NI,0:NJ,2):: CellCenter
    V(I,J,1) = 0.0 + 1.0*(CellCenter(I, J, 1)**0.0)*(CellCenter(I, J, 2)**0.0)
    V(I,J,2) = 0.0 + 1.0*(CellCenter(I, J, 2)**0.0)*(CellCenter(I, J, 1)**0.0)
end subroutine

Function LinearInterpolation(TargetPoint, PointIn, PointOut, funIn, funOut)
    implicit none
    real, dimension(2) :: TargetPoint, PointIn, PointOut
    real::LenIn, LenOut, funIn, funOut, LinearInterpolation

    LenIn  = Norm2(PointIn(:) - TargetPoint(:))!¬нутри €чейки
    LenOut = Norm2(PointOut(:) - TargetPoint(:))!—наружи €чейки €чейки

    LinearInterpolation = (funIn*LenOut + funOut*LenIn)/(LenOut + LenIn)
End Function


Subroutine CalcGradPExact(NI, NJ, CellCenter, GradPExact)
    implicit none
    integer:: I, J, NI, NJ
    Real,Dimension(0:NI,0:NJ,2):: GradPExact
    Real,Dimension(0:NI,0:NJ,2):: CellCenter
    DO J = 1, NJ-1
        DO I = 1, NI-1
            GradPExact(I,J,1) = 0.0 + 3.0*CellCenter(I, J, 1)**2.0
            GradPExact(I,J,2) = 0.0 + 3.0*CellCenter(I, J, 2)**2.0
        END DO
    END DO
end subroutine

Subroutine CalcDivVExact(NI, NJ, CellCenter, DivVExact)
    implicit none
    integer:: I, J, NI, NJ
    Real,Dimension(0:NI,0:NJ):: DivVExact
    Real,Dimension(0:NI,0:NJ,2):: CellCenter
    DO J = 1, NJ-1
        DO I = 1, NI-1
            DivVExact(I,J) = 0.0 + 1.0*(CellCenter(I, J, 1)**0.0)*(CellCenter(I, J, 1)**0.0) + &
             1.0*(CellCenter(I, J, 2)**0.0)*(CellCenter(I, J, 2)**0.0)
        END DO
    END DO
end subroutine

Subroutine CalcRotVExact(NI, NJ, CellCenter, RotVExact)
    implicit none
    integer:: I, J, NI, NJ
    Real,Dimension(0:NI,0:NJ):: RotVExact
    Real,Dimension(0:NI,0:NJ,2):: CellCenter
    DO J = 1, NJ-1
        DO I = 1, NI-1
            RotVExact(I,J) = 0.0 + 3.0*(CellCenter(I, J, 1)**2.0 )+ &
            0.0*(CellCenter(I, J, 2)**0.0)
        END DO
    END DO
end subroutine

Subroutine CalcLaplacianPExact(NI, NJ, CellCenter, LaplacianPExact)
    implicit none
    integer:: I, J, NI, NJ
    Real,Dimension(0:NI,0:NJ):: LaplacianPExact
    Real,Dimension(0:NI,0:NJ,2):: CellCenter
    DO I = 1, NI - 1
        DO J = 1, NJ - 1
            LaplacianPExact(I,J) = 0.0 + 12.0*(CellCenter(I, J, 1)**2.0 )+ &
            12.0*(CellCenter(I, J, 2)**2.0)
        END DO
    END DO
end subroutine



Subroutine CalcErrorVector(NI, NJ, Var, ExactVar, ErrorVar)
    implicit none
    integer:: I, J, NI, NJ
    Real,Dimension(0:NI,0:NJ,2):: Var, ExactVar, ErrorVar
    DO  J = 1, NJ-1
      DO  I = 1,NI-1
        ErrorVar(I,J,:) = abs(ExactVar(I,J,:) - Var(I,J,:))/abs(ExactVar(I,J,:))
      ENDDO
  ENDDO

end subroutine

Subroutine CalcError(NI, NJ, Var, ExactVar, ErrorVar)
    implicit none
    integer:: I, J, NI, NJ
    Real,Dimension(0:NI,0:NJ):: Var, ExactVar, ErrorVar
    DO  J = 1, NJ-1
      DO  I = 1,NI-1
        ErrorVar(I,J) = abs(ExactVar(I,J) - Var(I,J))/ abs(ExactVar(I,J))
      ENDDO
  ENDDO

end subroutine

