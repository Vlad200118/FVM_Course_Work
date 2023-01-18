SUBROUTINE B_CalcLaplacian(NI,NJ,P,GradP, CellVolume,IFaceVector,JFaceVector,IFaceCenter, &
                            JFaceCenter, CellCenter, LaplacianP)
    Implicit none
    Integer:: NI, NJ ,I, J, IFace, I_neighbor, J_neighbor
    Real :: DerivativeP, length
    Real :: LinearInterpolation
    Real,Dimension(NI-1,NJ-1):: CellVolume
    Real,Dimension(0:NI,0:NJ)::P, LaplacianP
    Real,Dimension(NI,NJ-1,2) ::IFaceCenter, IFaceVector
    Real,Dimension(NI-1,NJ,2) ::JFaceCenter, JFaceVector
    Real,Dimension(0:NI,0:NJ,2) ::CellCenter, GradP
    Real,Dimension(2)::FaceCenter, faceOrt, FaceVector, LROrt, GradLinearInterpolation
    Integer :: Bound


    LaplacianP(:,:) = 0.0
    DO J=1, NJ-1
        DO I=1, NI-1
        DO IFace=1, 4 !Цикл по соседям
        !DerivativeP = 0
        Bound = 0
            select case (IFace)
                case(1)!сосед снизу
                    I_neighbor = I
                    J_neighbor = J-1
                    FaceCenter(:) =  JFaceCenter(I,J,:)
                    FaceVector(:) = -JFaceVector(I,J,:)!минус из-за внешней нормали
                    IF(J .EQ. 1) THEN
                        Bound = 1
                    END IF

                case(2)!сосед сверху
                    I_neighbor = I
                    J_neighbor = J+1
                    FaceCenter(:) = JFaceCenter(I,J+1,:)
                    FaceVector(:) = JFaceVector(I,J+1,:)
                    IF(J .EQ. (NJ-1)) THEN
                        Bound = 1
                    END IF

                case(3)!сосед слева
                    I_neighbor = I-1
                    J_neighbor = J
                    FaceCenter(:) =  IFaceCenter(I,J,:)
                    FaceVector(:) = -IFaceVector(I,J,:)!минус из-за внешней нормали
                    IF(I .EQ. 1) THEN
                        Bound = 1
                    END IF
                case(4)!сосед справа
                    I_neighbor = I+1
                    J_neighbor = J
                    FaceCenter(:) = IFaceCenter(I+1,J,:)
                    FaceVector(:) = IFaceVector(I+1,J,:)
                    IF(I .EQ. (NI-1) ) THEN
                        Bound = 1
                    END IF
            end select
            !! ВЫЧИСЛЕНИЕ У ГРАНИЦ
                faceOrt(:) = FaceVector(:)/Norm2(FaceVector(:))

                IF (Bound .EQ. 1) THEN
                    length = Norm2(FaceCenter(:) - CellCenter(I, J, :))
                    LROrt(:) = (FaceCenter(:) - CellCenter(I, J, :))/ length
                    DerivativeP = 5.0/3.0*(P(I_neighbor, J_neighbor)-P(I,J)) / length - 2.0/3.0*DOT_PRODUCT(faceOrt(:),GradP(I,J,:))
                    GradLinearInterpolation(:) = GradP(I,J,:)
                ELSE
                    length = Norm2(CellCenter(I_neighbor, J_neighbor,:) - CellCenter(I, J, :))
                !!Вектор между LR
                    LROrt(:) = (CellCenter(I_neighbor, J_neighbor,:) - CellCenter(I, J, :))/ length

                    GradLinearInterpolation(1) = LinearInterpolation(FaceCenter(:), CellCenter(I,J,:), &
                                   CellCenter(I_neighbor,J_neighbor,:), GradP(I,J,1), GradP(I_neighbor,J_neighbor,1))

                    GradLinearInterpolation(2) = LinearInterpolation(FaceCenter(:), CellCenter(I,J,:), &
                                   CellCenter(I_neighbor,J_neighbor,:), GradP(I,J,2), GradP(I_neighbor,J_neighbor,2))

                    DerivativeP = (P(I_neighbor, J_neighbor) - P(I, J))/length
                    !DerivativeP = DerivativeP + DOT_PRODUCT(faceOrt(:)-LROrt(:),GradLinearInterpolation(:))
                END IF
                    DerivativeP = DerivativeP + DOT_PRODUCT(faceOrt(:)-LROrt(:),GradLinearInterpolation(:))
                    LaplacianP(I, J) = LaplacianP(I, J) + DerivativeP * Norm2(FaceVector(:))
                END DO

                LaplacianP(I, J) = LaplacianP(I, J)/CellVolume(I,J)
        END DO
    END DO
End Subroutine


