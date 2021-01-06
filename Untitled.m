z&Kx ~(";~&S~o:NX)
ConPIzX Ah&l, APPRh2, APPROXS, APPROX4
APPROXl(T) * (T * .5641896) / (.5 + (T * T))
APPROX2(T,U) = (T * (1.410474 + W.5641896))/ (.75 + (U *(S.+U)))
APPROXS(T) - ( 16.4955 + T * (20.20933 + T * (11.96482 +
-I T * (3.778987 + 0,5642236*T))))
. / ( 16.4955 + T * (38.82363 + T *
-l (39.27121 + T * (21.69274 + T * (6.699398 + T))))) 
162 F. ScHRElBR
C
C
C
100
C
C
C
200
C
C
300
C
C
C
C
C
C
400
C
C
APPROX4(T,U) =(T * (36183.31 - U * (3321.99 - U * (1540.767 - u
-I *(219.031 - U *(35.7668 - U *(1.320522 - U * -56419))))))
-I / ($2066.6 - U * (24322.8 - U'* (9022.23 - U * (21¶6.18
_I - U * (364.219 - U * (61.5704 - U * (1.34144 - U))))))))
IF (Y.GT.15.) TKEN
--~~~__-_-__~--__-~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_~~~~~________
all point8 dre in region I
DO 100, I=O,NX
T = CMPLX(Y,-X(1))
PRBFCT(1) = APPROXl(T)
_-_~_~----__~---__~_~~~~~--~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_______
ELSE IF (Y.LT.15. .AND. Y.QE.5.5) TRBN
_~~_~_-~~____--__~_~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_______
points are in region I or region II
DO 200, I=O,NX
T = CMPLX(Y,-X(1))
S = ABS(X(1)) + Y
IF (S .QE. 15.) TREN
PRBFCT(1) = APPROXl(T)
ELSE
=T*T
kBFCTCI) - APPROX2(T,U)
END IF
COHTINUE
~~_~_~~~__~~~_______~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_____
ELSE IF (Y.LT.5.5 .AND. Y.QT.0.75) TREN
_____~______________~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_~_____
DO 300, IIo,NX
T = cHPLX(Y,-X(1))
s - ABS(X(1)) + Y
IF (S .GE. 15.) THBN
PRBFCZf(1) = APPROXl(T)
ELSE IF (S.LT.5.5) THEN
PRBFCT(1) = APPROXs(T)
ELSE
U =T*T
PRBFCT(1) = APPROXZ(T,U)
END IF
CONTINUE
____________________~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ELSE
_______~_____~_~____~~~~~~~~~~~~~~~~~~--~-----~~~~---~~~~~~~~~~
Da 400, I=O,NX
T = CRPLX(Y,-X(1))
AX = ABS(X(1))
S =AX+Y
IF (S .GB. 15.) THEN
region I
PRBFCT(I)= APPROXl(T)
ELSE IF (S.LT.15.0 .AND. S.GE.5.5) THEN
region II
U=T*T
PRBFCT(I)= APPROXZ(T,U)
ELSE IF (S.LT.5.5 .AND. Y.GE.(O.l95*AX-0.176)) THEN
region III
PRBFCT(I)= APPROX3(T)
ELSE
region IV
U =T*T
PRBFCT(I)= CEXP(U) - APPROX4(T,U)
END IF
cotmNUR
______~_~_~~~_~_~~~_~~~~~~~~~~~~~~~~~~~-~~~~~~~~~~~-----~~~~~~~
END IF
REnJRN
END