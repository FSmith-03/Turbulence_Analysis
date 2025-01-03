! This automatically generated Fortran wrapper file allows codes
! written in Fortran to be called directly from C and translates all
! C-style arguments into expected Fortran-style arguments (with
! assumed size, local type declarations, etc.).


SUBROUTINE C_MATRIX_ROTATE(THETAX, THETAY, THETAZ, R_DIM_1, R_DIM_2, R) BIND(C)
  USE ISO_FORTRAN_ENV, ONLY: INT64
  IMPLICIT NONE
  REAL, INTENT(IN) :: THETAX
  REAL, INTENT(IN) :: THETAY
  REAL, INTENT(IN) :: THETAZ
  INTEGER(KIND=SELECTED_INT_KIND(18)), INTENT(OUT) :: R_DIM_1
  INTEGER(KIND=SELECTED_INT_KIND(18)), INTENT(OUT) :: R_DIM_2
  REAL, ALLOCATABLE, SAVE, DIMENSION(:,:) :: R_LOCAL
  INTEGER(KIND=INT64), INTENT(OUT) :: R

  INTERFACE
    SUBROUTINE MATRIX_ROTATE(THETAX, THETAY, THETAZ, R)
      IMPLICIT NONE
      REAL, INTENT(IN) :: THETAX
      REAL, INTENT(IN) :: THETAY
      REAL, INTENT(IN) :: THETAZ
      REAL, INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: R
    END SUBROUTINE MATRIX_ROTATE
  END INTERFACE

  CALL MATRIX_ROTATE(THETAX, THETAY, THETAZ, R_LOCAL)
  
  R_DIM_1 = SIZE(R_LOCAL,1)
  R_DIM_2 = SIZE(R_LOCAL,2)
  R = LOC(R_LOCAL(1,1))
END SUBROUTINE C_MATRIX_ROTATE


SUBROUTINE C_MAIN_CALCULATION(INPUT_INTS_DIM_1, INPUT_INTS, A_LIST_DIM_1, A_LIST_DIM_2, A_LIST, THETA_LIST_DIM_1, THETA_LIST_DIM_2,&
& THETA_LIST) BIND(C)
  IMPLICIT NONE
  INTEGER(KIND=SELECTED_INT_KIND(18)), INTENT(IN) :: INPUT_INTS_DIM_1
  INTEGER, INTENT(IN), DIMENSION(INPUT_INTS_DIM_1) :: INPUT_INTS
  INTEGER(KIND=SELECTED_INT_KIND(18)), INTENT(IN) :: A_LIST_DIM_1
  INTEGER(KIND=SELECTED_INT_KIND(18)), INTENT(IN) :: A_LIST_DIM_2
  REAL, INTENT(IN), DIMENSION(A_LIST_DIM_1,A_LIST_DIM_2) :: A_LIST
  INTEGER(KIND=SELECTED_INT_KIND(18)), INTENT(IN) :: THETA_LIST_DIM_1
  INTEGER(KIND=SELECTED_INT_KIND(18)), INTENT(IN) :: THETA_LIST_DIM_2
  REAL, INTENT(IN), DIMENSION(THETA_LIST_DIM_1,THETA_LIST_DIM_2) :: THETA_LIST

  INTERFACE
    SUBROUTINE MAIN_CALCULATION(INPUT_INTS, A_LIST, THETA_LIST)
      IMPLICIT NONE
      INTEGER, INTENT(IN), DIMENSION(:) :: INPUT_INTS
      REAL, INTENT(IN), DIMENSION(:,:) :: A_LIST
      REAL, INTENT(IN), DIMENSION(:,:) :: THETA_LIST
    END SUBROUTINE MAIN_CALCULATION
  END INTERFACE

  CALL MAIN_CALCULATION(INPUT_INTS, A_LIST, THETA_LIST)
END SUBROUTINE C_MAIN_CALCULATION


SUBROUTINE C_SENSOR_LINE_GENERATOR(X_BOUNDARY, NXF, POS_VECTOR_DIM_1, POS_VECTOR_DIM_2, POS_VECTOR) BIND(C)
  USE ISO_C_BINDING, ONLY: C_BOOL
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: X_BOUNDARY
  INTEGER, INTENT(IN) :: NXF
  INTEGER(KIND=SELECTED_INT_KIND(18)), INTENT(IN) :: POS_VECTOR_DIM_1
  INTEGER(KIND=SELECTED_INT_KIND(18)), INTENT(IN) :: POS_VECTOR_DIM_2
  REAL, DIMENSION(POS_VECTOR_DIM_1,POS_VECTOR_DIM_2) :: POS_VECTOR

  INTERFACE
    FUNCTION SENSOR_LINE_GENERATOR(X_BOUNDARY, NXF) RESULT(POS_VECTOR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: X_BOUNDARY
      INTEGER, INTENT(IN) :: NXF
      REAL, DIMENSION(3,NXF) :: POS_VECTOR
    END FUNCTION SENSOR_LINE_GENERATOR
  END INTERFACE

  POS_VECTOR = SENSOR_LINE_GENERATOR(X_BOUNDARY, NXF)
END SUBROUTINE C_SENSOR_LINE_GENERATOR

