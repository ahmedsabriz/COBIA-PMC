
// BZZMATH: Release 7.0

//	============================< POLY5.HPP >===================================
// * Examples: c:\bzzmath\exampled\dxpoly5.cpp											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-2001	Date Written.


#ifndef BZZ_POLY5_HPP
#define BZZ_POLY5_HPP

//	============================================================================
//	==========================< class BzzPoly5 >================================
//	============================================================================

class BzzPoly5
{
private:
	double a,b,c,d,e,f;
	double x0,y0,yp0,x1,y1,yp1,x2,y2,yp2;
	double x10,y10x10,x20,x21,y21x21;

public:
	// default
	BzzPoly5(void)
		{
		a = b = c = d = e = f = 0.;
		}

	BzzPoly5(double xx0,double yy0,double yyp0,
		double xx1,double yy1,double yyp1);

	void operator()(double xx0,double yy0,double yyp0,
		double xx1,double yy1,double yyp1);
	void operator()(double xx2,double yy2,double yyp2);
	void operator()(double xx2,double yy2);
	double operator()(double x);
};


#endif // BZZ_POLY5_HPP

