RF<a,c,d>:=RationalFunctionField(Rationals(),3); P:=PolynomialRing(RF);
PR<X0,Z0,Y0,T0,X1,Z1,Y1,T1>:=PolynomialRing(RF,8);

//((X2:Z2),(Y2:Z2)):=((X0:Z0),(Y0:Z0))-((X1:Z1),(Y1:Z1))
X2:=(X0*Z1-Z0*X1)*(T0*T1+d*Y0*Y1); Z2:=(Z0*Z1+a*X0*X1)*(T0*T1-d*Y0*Y1);
Y2:=(Z0*Z1+a*X0*X1)*(Y0*T1-T0*Y1); T2:=(Z0*Z1-a*X0*X1)*(T0*T1+d*Y0*Y1);
//((X3:Z3),(Y3:Z3)):=((X0:Z0),(Y0:Z0))+((X1:Z1),(Y1:Z1))
X3:=(X0*Z1+Z0*X1)*(T0*T1-d*Y0*Y1); Z3:=(Z0*Z1-a*X0*X1)*(T0*T1+d*Y0*Y1);
Y3:=(Z0*Z1-a*X0*X1)*(Y0*T1+T0*Y1); T3:=(Z0*Z1+a*X0*X1)*(T0*T1-d*Y0*Y1);

//Projective differential addition formulas without precomputation.
X3d:=Z2*(X0^2*Z1^2-Z0^2*X1^2); Z3d:=X2*(Z0^2*Z1^2-a^2*X0^2*X1^2);
Y3d:=T2*(Y0^2*T1^2-T0^2*Y1^2); T3d:=Y2*(T0^2*T1^2-d^2*Y0^2*Y1^2);

//Check modulo the quotient relations.
QR<X0,Z0,Y0,T0,X1,Z1,Y1,T1>:=RingOfFractions(quo<PR|
    Y0*T0*(Z0^2+a*X0^2)-c*X0*Z0*(T0^2+d*Y0^2),
    Y1*T1*(Z1^2+a*X1^2)-c*X1*Z1*(T1^2+d*Y1^2)
>);
QR!(Y2*T2*(Z2^2+a*X2^2)-c*X2*Z2*(T2^2+d*Y2^2));
QR!(Y3*T3*(Z3^2+a*X3^2)-c*X3*Z3*(T3^2+d*Y3^2));
QR!(X3/Z3-X3d/Z3d);
QR!(Y3/T3-Y3d/T3d);

RF<c,d>:=RationalFunctionField(Rationals(),2); P:=PolynomialRing(RF);
PR<X0,Z0,Y0,T0,X1,Z1,Y1,T1>:=PolynomialRing(RF,8);

//Assumption.
a:=1;

//((X2:Z2),(Y2:Z2)):=((X0:Z0),(Y0:Z0))-((X1:Z1),(Y1:Z1))
X2:=(X0*Z1-Z0*X1)*(T0*T1+d*Y0*Y1); Z2:=(Z0*Z1+X0*X1)*(T0*T1-d*Y0*Y1);
Y2:=(Z0*Z1+X0*X1)*(Y0*T1-T0*Y1); T2:=(Z0*Z1-X0*X1)*(T0*T1+d*Y0*Y1);
//((X3:Z3),(Y3:Z3)):=((X0:Z0),(Y0:Z0))+((X1:Z1),(Y1:Z1))
X3:=(X0*Z1+Z0*X1)*(T0*T1-d*Y0*Y1); Z3:=(Z0*Z1-X0*X1)*(T0*T1+d*Y0*Y1);
Y3:=(Z0*Z1-X0*X1)*(Y0*T1+T0*Y1); T3:=(Z0*Z1+X0*X1)*(T0*T1-d*Y0*Y1);

//Precomputation.
kappadash:=(X0^2-Z0^2)/(X0^2+Z0^2);

//Projective differential addition formulas with precomputation.
X3d:=Z2*((X1^2-Z1^2)-kappadash*(X1^2+Z1^2));
Z3d:=X2*((X1^2-Z1^2)+kappadash*(X1^2+Z1^2));

//Check modulo the quotient relations.
QR<X0,Z0,Y0,T0,X1,Z1,Y1,T1>:=RingOfFractions(quo<PR|
    Y0*T0*(Z0^2+X0^2)-c*X0*Z0*(T0^2+d*Y0^2),
    Y1*T1*(Z1^2+X1^2)-c*X1*Z1*(T1^2+d*Y1^2)
>);
QR!(Y2*T2*(Z2^2+X2^2)-c*X2*Z2*(T2^2+d*Y2^2));
QR!(Y3*T3*(Z3^2+X3^2)-c*X3*Z3*(T3^2+d*Y3^2));
QR!(X3/Z3-X3d/Z3d);

