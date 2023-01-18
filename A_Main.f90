Program Main

  character(*), parameter:: InputFile='input.txt', OutputFile='data.plt' ! names of input and output files
  character MeshFile*30, SolutionFile*30, ResultFile*30, str       ! name of file with computational mesh
  integer, parameter:: IO = 12 ! input-output unit
  integer :: SolutionFileRead, SolveConvDiff, GGI
  Real :: SCHEME

  real,allocatable,dimension(:,:):: X,Y,CellVolume ! scalar arrays
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector ! vector arrays

  real,allocatable,dimension(:,:,:)::V
  real,allocatable,dimension(:,:)::P, T

  real,allocatable,dimension(:,:) :: DivV, DivVExact, DivVError
  real,allocatable,dimension(:,:) :: RotV, RotVExact, RotVError
  real,allocatable,dimension(:,:,:)::GradP, GradPExact, GradPError!
  real,allocatable,dimension(:,:):: LaplacianP, LaplacianPExact, LaplacianPError

!===  READ INPUT FILE ===
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  READ(IO,*) MeshFile  ! read name of file with computational mesh
  READ(IO,*) SolutionFileRead
  READ(IO,*) SolutionFile  ! read name of file with computational mesh
  READ(IO,*) ResultFile
  READ(IO,*) GGI
  READ(IO,*) SCHEME
  READ(IO,*) SolveConvDiff
  CLOSE(IO)

!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ
  WRITE(*,*) 'NI, NJ = ',NI,NJ
  WRITE(*,*) 'GGI = ', GGI
  WRITE(*,*) 'Divergance SCHEME = ', SCHEME

!=== ALLOCATE ALL ARRAYS ===
  WRITE(*,*) 'Allocate arrays'       
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(T(0:NI,0:NJ))
  allocate(V(0:NI,0:NJ,2))   ! Velocity
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes    
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter( NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector( NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter( NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector( NI-1,NJ,2)) ! Face Vectors for I-faces

  allocate(DivV( 0:NI,0:NJ), DivVExact( 0:NI,0:NJ), DivVError( 0:NI,0:NJ))
  allocate(RotV( 0:NI,0:NJ), RotVExact( 0:NI,0:NJ), RotVError( 0:NI,0:NJ))
  allocate(GradP( 0:NI,0:NJ,2), GradPExact( 0:NI,0:NJ,2), GradPError( 0:NI,0:NJ,2))
  allocate(LaplacianP( 0:NI,0:NJ), LaplacianPExact( 0:NI,0:NJ), LaplacianPError( 0:NI,0:NJ))

!===  READ GRID ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile
  IF (SolutionFileRead .EQ. 1) THEN
    READ(IO,*) ((X(I,J),Y(I,J),rtmp,I=1,NI),J=1,NJ)
  ELSE
    READ(IO,*) ((X(I,J),Y(I,J),I=1,NI),J=1,NJ)
  END IF
  CLOSE(IO)

!=== CALCULATE METRIC ===
  WRITE(*,*) 'Calculate metric'       
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector) 
  
!=== INITIATE FIELDS ===
  WRITE(*,*) 'Initiate fields'       
  DO  J = 0,NJ
    DO  I = 0,NI
      Call CalcInitV(I, J, NI, NJ, CellCenter, V);
      P(I,J) = Pressure(CellCenter(I,J,1),CellCenter(I,J,2))
    ENDDO
  ENDDO

!=== READ SOLUTION ===
    IF (SolutionFileRead .EQ. 1) THEN
    WRITE(*,*) 'Read solution'
    OPEN (IO, FILE = SolutionFile)
    READ(IO,*) str
    READ(IO,*) str
    READ(IO,*) ((rtmp,rtmp, V(I,J,1),V(I,J,2),rtmp,P(I,J),rtmp,rtmp, I=0, NI), J=0, NJ)
    CLOSE(IO)
    END IF

  WRITE(*,*) 'Calculate derivatives'

!=== CALCULATE GRADIENT ===
  WRITE(*,*) 'Calculate gradient'

  If (GGI .EQ. 0) THEN
    Call B_CalcGradientGG(NI,NJ,P,CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter, GradP)
  ELSE
    Call B_CalcGradientGGI(NI,NJ,P,CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter, GradP)
  END IF
  Call CalcGradPExact(NI, NJ, CellCenter, GradPExact)
  Call CalcErrorVector(NI, NJ, GradP, GradPExact, GradPError)

!=== CALCULATE LAPLACIAN ===
  WRITE(*,*) 'Calculate Laplacian'
  Call B_CalcLaplacian(NI,NJ,P,GradP,CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter, LaplacianP)

  Call CalcLaplacianPExact(NI, NJ, CellCenter, LaplacianPExact)
  Call CalcError(NI,NJ, LaplacianP, LaplacianPExact, LaplacianPError)

!=== CALCULATE DIVERGANCE ===
  WRITE(*,*) 'Calculate Divergance'
  !Call B_CalcDiverganceV(NI,NJ,V,CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter,DivV)

  Call B_CalcDivergancePV(NI,NJ,V,P,GradP,CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter,DivV,SCHEME)

  Call CalcDivVExact(NI, NJ, CellCenter, DivVExact)
  Call CalcError (NI, NJ, DivV, DivVExact, DivVError)

!=== CALCULATE ROTOR ===
  WRITE(*,*) 'Calculate Rotor'
  Call B_CalcRotorV(NI,NJ,V,CellVolume,IFaceVector,JFaceVector, IFaceCenter,JFaceCenter,CellCenter,RotV)

  Call CalcRotVExact(NI, NJ, CellCenter, RotVExact)
  Call CalcError (NI, NJ, RotV, RotVExact, RotVError)


  if (SolveConvDiff .EQ. 1) THEN
  WRITE(*,*) 'Conv_Dif'
  Call Calc_Conv_Dif_T (NI, NJ, V, T, CellVolume, IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter, SCHEME)

  WRITE(*,*) 'Output fields to file: ', ResultFile
  Open(IO,FILE=ResultFile)
  Call B_OutputResultConv_Diff(IO,NI,NJ,X,Y,P,V,T)
  Close(IO)

  End if

 if(SolutionFileRead .EQ. 0 .AND.SolveConvDiff .EQ. 0) then
!=== OUTPUT FIELDS ===
  WRITE(*,*) 'Output fields to file: ', ResultFile
  Open(IO,FILE=ResultFile)
  Call B_OutputFields(IO,NI,NJ,X,Y,P,V,GradP,GradPExact, GradPError,DivV, DivVExact,DivVError, RotV, RotVExact, RotVError, LaplacianP, LaplacianPExact, LaplacianPError)
  Close(IO)
End if

 if(SolutionFileRead .EQ. 1 .AND.SolveConvDiff .EQ. 0) then
  WRITE(*,*) 'Output fields to file: ', ResultFile
  Open(IO,FILE=ResultFile)
  Call B_OutputResult(IO,NI,NJ,X,Y,P,V,GradP,DivV,RotV, LaplacianP)
  Close(IO)
End if
  Deallocate (X,Y,CellVolume,IFaceVector,JFaceVector, IFaceCenter,JFaceCenter,CellCenter)
  Deallocate (P, V, T)
  Deallocate (GradP,GradPExact,GradPError)


END PROGRAM Main
