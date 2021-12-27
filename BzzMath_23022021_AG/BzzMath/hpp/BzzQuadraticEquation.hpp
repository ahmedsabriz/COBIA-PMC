// BZZMATH: Release 7.0

//	========================< BzzQuadraticEquation.hpp >==========================
//	* Class BzzQuadraticEquation for quadratic equations									*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 3)					*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// *																									*
// * Examples: BzzMath\examples\BzzMathBasic\QuadraticEquation\					*
//					Quadratic\Quadratic.cpp										*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	04-1996	Date Written

// Solve the quadratic equation:
// ax2 + bx + c = 0

////////////////// Release 5.0
// 08-2003	Added operator(a,b,c).
// 09_2008 Added GetStationaryPoint function

#ifndef BZZ_QUADRATIC_DOUBLE_HPP
#define BZZ_QUADRATIC_DOUBLE_HPP

//	============================================================================
//	====================< class BzzQuadraticEquation >============================
//	============================================================================

class BzzQuadraticEquation : public BzzBaseClass
{
private:
	enum
	{
		NOT_DEFINED = 0,
		IMPOSSIBLE = 1,
		IDENTITY = 2,
		ONE_REAL = 3,
		ONE_REAL_NEGATIVE_OVERFLOW = 4,
		ONE_REAL_POSITIVE_OVERFLOW = 5,
		TWO_REAL = 6,
		TWO_REAL_FIRST_NEGATIVE_OVERFLOW = 7,
		TWO_REAL_SECOND_POSITIVE_OVERFLOW = 8,
		TWO_REAL_FIRST_NEGATIVE_OVERFLOW_SECOND_POSITIVE_OVERFLOW = 9,
		COMPLEX = 10
	};
	int quadraticSolutionType;

	double a, b, c, x1, x2;
public:
	BzzQuadraticEquation(double aa, double bb, double cc)
	{
		(*this)(aa, bb, cc);
	}
	void operator()(double aa, double bb, double cc);
	BzzQuadraticEquation(void)
	{
		a = b = c = 0.;
		quadraticSolutionType = NOT_DEFINED;
	}
	virtual void ObjectBzzPrint(void);
	int QuadraticSolutionType(void)
	{
		return quadraticSolutionType;
	}
	char GetX1(double* xx1);
	char GetX2(double* xx2);
	char GetReal(double* xx1);
	char GetImaginary(double* xx2);
	// 1 OK; 0 No stationary point or overflaw;
	char GetStationaryPoint(double* stationaryPoint);
};

#endif // BZZ_QUADRATIC_DOUBLE_HPP