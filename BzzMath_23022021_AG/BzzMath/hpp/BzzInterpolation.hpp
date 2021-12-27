// BZZMATH: Release 7.0

//	=========================< BzzInterpolation.hpp >===========================
//	* BzzInterpolation Class																	*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 14,15)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: BzzMath\Examples\BzzMathBasic\Interpolation\Interpolation.cpp	*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	07-1989	Date Written
//	04-1994	Revision

//	============================================================================
//	******	BzzInterpolation constructors:												*
//	* BzzInterpolation p; // default															*
//	* BzzInterpolation p(x,y);																	*
//	============================================================================

#ifndef BZZ_INTERPOLATION_HPP
#define BZZ_INTERPOLATION_HPP

//	============================================================================
//	==========================< class BzzInterpolation >========================
//	============================================================================
class BzzInterpolation : public BzzBaseClass
{
private:
	friend void Delete(BzzInterpolation* a);
	int n, numPoints;
	double* a, * x, * y, * d, * c, * w, * xz, * yz;
	void Initialize(int nn);
	void Delete(void);
	char parLagrange;
	char parNewton;
	void NewtonParameters(void);
	void LagrangeParameters(void);
	void Ordering(double z);

public:
	BzzInterpolation(void);
	BzzInterpolation(BzzVector& xx, BzzVector& yy);
	double Polynomial(double zz);
	double Newton(double zz);
	double Lagrange(double zz);
	BzzVector Neville(double zz);

	double Rational(double zz);
	//	double Thiele(double zz);
	BzzVector BulirschStoer(double zz);

	void operator()(BzzVector& xx, BzzVector& yy);

	~BzzInterpolation(void);
	virtual void ObjectBzzPrint(void);
};

#endif // BZZ_INTERPOLATION_HPP