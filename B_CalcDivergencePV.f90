Subroutine B_CalcDivergancePV(NI,NJ,V,P,GradP,CellVolume,IFaceVector,JFaceVector,IFaceCenter, &
                            JFaceCenter,CellCenter, DivPV, SCHEME)
    Implicit none
    Integer :: NI, NJ ,I, J, IFace, I_neighbor, J_neighbor
    Real :: SCHEME
    Real :: LinearInterpolation, PFace, PLinearInterpolation, PFOU, PSOU, GFace
    Real, Dimension (2) :: VLinearInterpolation
    Real,Dimension(NI-1,NJ-1) :: CellVolume
    Real,Dimension(0:NI,0:NJ,2) :: V, GradP
    Real,Dimension(0:NI,0:NJ) :: DivPV, P
    Real,Dimension(NI,NJ-1,2) :: IFaceCenter, IFaceVector
    Real,Dimension(NI-1,NJ,2) :: JFaceCenter, JFaceVector
    Real,Dimension(0:NI,0:NJ,2) :: CellCenter
    Real,Dimension(2) :: FaceCenter, FaceVector

    DivPV(:,:) = 0.0
    PFace = 0.0
    DO J =1, NJ-1
        DO I=1, NI-1
            DO IFace=1, 4 !Цикл по соседям
                select case (IFace)
                case(1)!сосед снизу
                    I_neighbor = I
                    J_neighbor = J-1
                    FaceCenter(:) =  JFaceCenter(I,J,:)
                    FaceVector(:) = -JFaceVector(I,J,:)!минус из-за внешней нормали

                case(2)!сосед сверху
                    I_neighbor = I
                    J_neighbor = J+1
                    FaceCenter(:) = JFaceCenter(I,J+1,:)
                    FaceVector(:) = JFaceVector(I,J+1,:)

                case(3)!сосед слева
                    I_neighbor = I-1
                    J_neighbor = J
                    FaceCenter(:) =  IFaceCenter(I,J,:)
                    FaceVector(:) = -IFaceVector(I,J,:)!минус из-за внешней нормали

                case(4)!сосед справа
                    I_neighbor = I+1
                    J_neighbor = J
                    FaceCenter(:) = IFaceCenter(I+1,J,:)
                    FaceVector(:) = IFaceVector(I+1,J,:)

                end select

                    VLinearInterpolation(1) = LinearInterpolation(FaceCenter(:), CellCenter(I,J,:), &
                                   CellCenter(I_neighbor,J_neighbor,:),V(I,J,1), V(I_neighbor,J_neighbor,1))!V_linear_x

                    VLinearInterpolation(2) = LinearInterpolation(FaceCenter(:), CellCenter(I,J,:), &
                                   CellCenter(I_neighbor,J_neighbor,:),V(I,J,2), V(I_neighbor,J_neighbor,2))!V_linear_y

                    PLinearInterpolation = LinearInterpolation(FaceCenter(:), CellCenter(I,J,:), &
                                   CellCenter(I_neighbor,J_neighbor,:),P(I,J), P(I_neighbor,J_neighbor))

                    GFace = DOT_PRODUCT(VLinearInterpolation(:), FaceVector(:))

                    IF (SCHEME .EQ. 2.0) then
                        IF (GFace .GE. 0) then
                            PFOU = P(I,J)
                        else
                            PFOU = P(I_neighbor,J_neighbor)
                        end if
                        PFace = PFOU
                    end if

                    IF (SCHEME .LE. 1.0) then
                        IF (GFace .GE. 0) then
                            PSOU = P(I,J) + DOT_PRODUCT(GradP(I,J,:),FaceCenter(:) - CellCenter(I,J,:))
                        else
                            PSOU = P(I_neighbor,J_neighbor) + DOT_PRODUCT(GradP(I_neighbor,J_neighbor,:),FaceCenter(:) - CellCenter(I_neighbor,J_neighbor,:))
                        end if
                        PFace = (1.0-SCHEME)*PLinearInterpolation + SCHEME * PSOU
                    end if


                    DivPV(I,J) = DivPV(I,J) + PFace*GFace

            END DO
            DivPV(I,J)=DivPV(I,J)/CellVolume(I,J)

        END DO
    END DO
End Subroutine
