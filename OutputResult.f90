Subroutine B_OutputResult(IO,NI,NJ,X,Y,P,V,GradP,DivV,RotV, LaplacianP)
  integer :: NI, NJ
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P
  Real,Dimension(0:NI,0:NJ):: DivV, RotV, LaplacianP
  Real,Dimension(0:NI,0:NJ,2):: V
  Real,Dimension(0:NI,0:NJ,2):: GradP

  Write(IO,*) 'VARIABLES = "X", "Y", "P","V_x","V_y","gradP_x", "gradP_y","DivV", "RotV", "LaplacianP"'
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-30]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ)
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)

  Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)

  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,2)

  Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,2)

  Write(IO,'(100F14.7)') DivV(1:NI-1,1:NJ-1)

  Write(IO,'(100F14.7)') RotV(1:NI-1,1:NJ-1)

  Write(IO,'(100F14.7)') LaplacianP(1:NI-1,1:NJ-1)


End Subroutine


Subroutine B_OutputResultConv_Diff(IO,NI,NJ,X,Y,P,V,T)
  integer :: NI, NJ
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::T, P
  Real,Dimension(0:NI,0:NJ,2):: V

  Write(IO,*) 'VARIABLES = "X", "Y", "P","V_x","V_y", "T"'
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-30]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ)
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)

  Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)

  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,2)

  Write(IO,'(100F14.7)') T(1:NI-1,1:NJ-1)

End Subroutine
