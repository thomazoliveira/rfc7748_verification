RF<a,d>:=RationalFunctionField(Rationals(),2);
P:=PolynomialRing(RF); K<r>:=quo<P|P.1^2-d/a>;
PR<X0,Z0,Y0,T0,X1,Z1,Y1,T1>:=PolynomialRing(K,8);

//((X2:Z2),(Y2:Z2)):=((X0:Z0),(Y0:Z0))-((X1:Z1),(Y1:Z1))
X2:=X0*Y1*T0*Z1-Y0*X1*Z0*T1; Z2:=Z0*T1*T0*Z1-d*X0*X1*Y0*Y1;
Y2:=Y0*Y1*Z0*Z1+a*X0*X1*T0*T1; T2:=Z0*T1*T0*Z1+d*X0*X1*Y0*Y1;
//((X3:Z3),(Y3:Z3)):=((X0:Z0),(Y0:Z0))+((X1:Z1),(Y1:Z1))
X3:=X0*Y1*T0*Z1+Y0*X1*Z0*T1; Z3:=Z0*T1*T0*Z1+d*X0*X1*Y0*Y1;
Y3:=Y0*Y1*Z0*Z1-a*X0*X1*T0*T1; T3:=Z0*T1*T0*Z1-d*X0*X1*Y0*Y1;

//Assumed input coordinates.
rYY0:=r*Y0^2; TT0:=T0^2; rYY1:=r*Y1^2; TT1:=T1^2; rYY2:=r*Y2^2; TT2:=T2^2;
//Precomputation.
epsilon:=(1-r)/(1+r)*(TT0+rYY0)/(TT0-rYY0);

//Projective differential addition formulas with precomputation.
rYY3:= TT2*((TT1-rYY1)-epsilon*(TT1+rYY1))^2;
TT3:=rYY2*((TT1-rYY1)+epsilon*(TT1+rYY1))^2;

//Check modulo the quotient relations.
QR<X0,Z0,Y0,T0,X1,Z1,Y1,T1>:=RingOfFractions(quo<PR|
    a*X0^2*T0^2+Y0^2*Z0^2-Z0^2*T0^2-d*X0^2*Y0^2,
    a*X1^2*T1^2+Y1^2*Z1^2-Z1^2*T1^2-d*X1^2*Y1^2
>);
QR!(a*X3^2*T3^2+Y3^2*Z3^2-Z3^2*T3^2-d*X3^2*Y3^2);
QR!(r*Y3^2/T3^2-rYY3/TT3);
