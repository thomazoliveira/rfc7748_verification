RF<A,B>:=RationalFunctionField(Rationals(),2); P:=PolynomialRing(RF);
PR<y0,y1,x0,x1>:=PolynomialRing(RF,4);

//(x3,y3):=(x0,y0)+(x1,y1) Source: Explicit Formulas Database (EFD).
x3:=B*(y0-y1)^2/(x0-x1)^2-A-x0-x1;
y3:=(2*x0+x1+A)*(y0-y1)/(x0-x1)-B*(y0-y1)^3/(x0-x1)^3-y0;
//(x2,y2):=(x0,y0)-(x1,y1)
x2:=B*(y0+y1)^2/(x0-x1)^2-A-x0-x1;
y2:=(2*x0+x1+A)*(y0+y1)/(x0-x1)-B*(y0+y1)^3/(x0-x1)^3-y0;
//Proposed alternative formulas.
x3t:=(1-x0*x1)/(x0-x1)*(x0*y1-y0*x1)/(x0*y1+y0*x1);
y3t:=(y0*(x1^2-1)-y1*(x0^2-1))/(x0-x1)^2*(x0*y1-y0*x1)/(x0*y1+y0*x1);
//Differential addition formulas.
x3d:=(1-x0*x1)^2/(x2*(x0-x1)^2);
y3d:=(y0^2*(1-x1^2)^2-y1^2*(1-x0^2)^2)/(y2*(x0-x1)^4);

//Check modulo the quotient relations.
QR<x0,y0,x1,y1>:=RingOfFractions(quo<PR|
    B*y0^2-(x0^3+A*x0^2+x0),B*y1^2-(x1^3+A*x1^2+x1)
>);
QR!(B*y2^2-(x2^3+A*x2^2+x2)); QR!(B*y3^2-(x3^3+A*x3^2+x3));
QR!(x3-x3d); QR!(y3-y3d); QR!(x3-x3t); QR!(y3-y3t);

RF<A,B>:=RationalFunctionField(Rationals(),2); P:=PolynomialRing(RF);
PR<X0,Z0,Y0,T0,X1,Z1,Y1,T1>:=PolynomialRing(RF,8);

//((X2:Z2),(Y2:Z2)):=((X0:Z0),(Y0:Z0))-((X1:Z1),(Y1:Z1))
X2:=(Z0*Z1-X0*X1)*(X0*T0*Y1*Z1+Y0*Z0*X1*T1);
Z2:=(X0*Z1-Z0*X1)*(X0*T0*Y1*Z1-Y0*Z0*X1*T1);
Y2:=(X0*T0*Y1*Z1+Y0*Z0*X1*T1)*(Y0*Z0^2*T1*(X1^2-Z1^2)+T0*(X0^2-Z0^2)*Y1*Z1^2);
T2:=(X0*T0*Y1*Z1-Y0*Z0*X1*T1)*T0*T1*(X0*Z1-Z0*X1)^2;
//((X3:Z3),(Y3:Z3)):=((X0:Z0),(Y0:Z0))+((X1:Z1),(Y1:Z1))
X3:=(Z0*Z1-X0*X1)*(X0*T0*Y1*Z1-Y0*Z0*X1*T1);
Z3:=(X0*Z1-Z0*X1)*(X0*T0*Y1*Z1+Y0*Z0*X1*T1);
Y3:=(X0*T0*Y1*Z1-Y0*Z0*X1*T1)*(Y0*Z0^2*T1*(X1^2-Z1^2)-T0*(X0^2-Z0^2)*Y1*Z1^2);
T3:=(X0*T0*Y1*Z1+Y0*Z0*X1*T1)*T0*T1*(X0*Z1-Z0*X1)^2;

//Precomputation.
mu:=(X0+Z0)/(X0-Z0);

//Projective differential addition formulas with precomputation.
X3d:=Z2*((X1+Z1)+mu*(X1-Z1))^2;
Z3d:=X2*((X1+Z1)-mu*(X1-Z1))^2;

//Check modulo the quotient relations.
QR<X0,Z0,Y0,T0,X1,Z1,Y1,T1>:=RingOfFractions(quo<PR|
    B*Y0^2*Z0^3-(X0^3*T0^2+A*X0^2*T0^2*Z0+X0*T0^2*Z0^2),
    B*Y1^2*Z1^3-(X1^3*T1^2+A*X1^2*T1^2*Z1+X1*T1^2*Z1^2)
>);
QR!(B*Y2^2*Z2^3-(X2^3*T2^2+A*X2^2*T2^2*Z2+X2*T2^2*Z2^2));
QR!(B*Y3^2*Z3^3-(X3^3*T3^2+A*X3^2*T3^2*Z3+X3*T3^2*Z3^2));
QR!(X3/Z3-X3d/Z3d);

