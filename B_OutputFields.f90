Subroutine B_OutputFields(IO,NI,NJ,X,Y,P,V,GradP,GradPExact, GradPError,DivV, DivVExact,DivVError, &
                            RotV, RotVExact, RotVError, LaplacianP, LaplacianPExact, LaplacianPError)
  integer :: NI, NJ
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P
  Real,Dimension(0:NI,0:NJ):: DivV, DivVExact, DivVError
  Real,Dimension(0:NI,0:NJ):: RotV, RotVExact, RotVError
  Real,Dimension(0:NI,0:NJ):: LaplacianP, LaplacianPExact, LaplacianPError
  Real,Dimension(0:NI,0:NJ,2):: V
  Real,Dimension(0:NI,0:NJ,2):: GradP, GradPExact, GradPError

  Write(IO,*) 'VARIABLES = "X", "Y", "P","V_x","V_y","gradP_x", "gradP_y","gradP_xExact","gradP_yExact","gradP_xError", "gradP_yError","DivV", "DivVExact", "DivVError", "RotV", "RotVExact", "RotVError", "LaplacianP", "LaplacianPExact", "LaplacianPError"'
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-30]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ) 
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)

  Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)

  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,2)

  Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F14.7)') GradPExact(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') GradPExact(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F14.7)') GradPError(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') GradPError(1:NI-1,1:NJ-1,2)

  Write(IO,'(100F14.7)') DivV(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') DivVExact(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') DivVError(1:NI-1,1:NJ-1)

  Write(IO,'(100F14.7)') RotV(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') RotVExact(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') RotVError(1:NI-1,1:NJ-1)

  Write(IO,'(100F14.7)') LaplacianP(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') LaplacianPExact(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') LaplacianPError(1:NI-1,1:NJ-1)

End Subroutine 
