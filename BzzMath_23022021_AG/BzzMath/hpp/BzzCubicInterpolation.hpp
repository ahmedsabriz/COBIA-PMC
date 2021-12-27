// BZZMATH: Release 7.0

//	=====================< BzzCubicInterpolation.cpp >====================
//	* classes BzzCubicBase, BzzCubicSpline,								*
//	* BzzCubicHermite, BzzCubicSmooth										*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 14, 15)			*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: BzzMath\Examples\BzzMathBasic\CubicInterpolation\					*
// *				CubicHermite\CubicHermite.cpp								*
// * Examples: BzzMath\Examples\BzzMathBasic\CubicInterpolation\					*
// *				CubicSmooth\CubicSmooth.cpp								*
// * Examples: BzzMath\Examples\BzzMathBasic\CubicInterpolation\					*
// *				CubicSpline\CubicSpline.cpp								*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	04-1992	Date Written.
//	01-2000	Added operator () for initializing.

////////////////// Release 5.0
//	10-2003	Added GetXMin and  GetXMax functions

// BZZMATH: Release 7.0
//	10-2011	Added BzzFiniteElementPointsForHermitePolynomial class
//	10-2011	Added BzzPiecewiseHermitePolynomialInterpolation class
//	10-2011	Added BzzPiecewiseLinearInterpolation class
//	10-2011	Added BzzPiecewiseLinearInterpolation class
//	01-2012	Added BzzFiniteElementPoints class
//	01-2012	Added BzzPiecewisePolynomialInterpolation class

#ifndef BZZ_CUBIC_DOUBLE_HPP
#define BZZ_CUBIC_DOUBLE_HPP

//	============================================================================
//	========================< class BzzCubicBase >========================
//	============================================================================

class BzzCubicBase
{
protected:
	static const char* const BZZ_ERROR;
	int numPoints;
	BzzVector	x,
		y,
		a, b, c; // a(x-xi)3 + ... + d

// generic constructor
	BzzCubicBase(BzzVector& xx, BzzVector& yy);
	void Initialize(BzzVector& xx, BzzVector& yy);

public:
	//	============================================================================
	//	**************************< constructors >**********************************
	//	============================================================================

		// default
	BzzCubicBase(void)
	{
		numPoints = 0;
	}
	// copy-initializer

//	============================================================================
//	**************************< destructor >************************************
//	============================================================================

	~BzzCubicBase(void) {};

	//	============================================================================
	//	*******************************< Functions >********************************
	//	============================================================================

	double Get(double xx);
	double Get(double xx, char ider);
	double BzzIntegral(void);
	double GetXMin(void)
	{
		return x[1];
	}
	double GetXMax(void)
	{
		return x[numPoints];
	}
};

//	============================================================================
//	=========================< class BzzCubicSpline >=====================
//	============================================================================
class BzzCubicSpline : public BzzCubicBase
{
	//enum SplineType
	//		{
	//		NATURAL,
	//		ASSIGNED_SECOND_DERIVATIVE,
	//		SECOND_DERIVATIVE_ADJACENT,
	//		ASSIGNED_FIRST_DERIVATIVE,
	//		PERIODIC
	//		};
private:
	SplineType splineType;
	double d1, dn;
	double y21, y2n;
	double alfa, beta;

	void GetCoefficients(void);
	void Initialize(BzzVector& xx, BzzVector& yy,
		SplineType sp = SECOND_DERIVATIVE_ADJACENT);
	void Initialize(BzzVector& xx, BzzVector& yy, SplineType sp,
		double d1, double d2);

public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzCubicSpline(void) {}
	// copy-initializer
	// NATURAL, PERIODIC,SECOND_DERIVATIVE_ADJACENT(alfa = beta = 1)
	BzzCubicSpline(BzzVector& xx, BzzVector& yy,
		SplineType sp = SECOND_DERIVATIVE_ADJACENT);

	// ASSIGNED_SECOND_DERIVATIVE
	// SECOND_DERIVATIVE_ADJACENT
	// ASSIGNED_FIRST_DERIVATIVE
	BzzCubicSpline(BzzVector& xx, BzzVector& yy, SplineType sp,
		double d1, double d2);

	void operator () (BzzVector& xx, BzzVector& yy,
		SplineType sp = SECOND_DERIVATIVE_ADJACENT);
	void operator () (BzzVector& xx, BzzVector& yy, SplineType sp,
		double d1, double d2);
	double operator () (double xx) { return Get(xx); }
	double operator () (double xx, char ider) { return Get(xx, ider); };

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzCubicSpline(void) {};

	//	============================================================================
	//	*************************< assignment operator >****************************
	//	============================================================================
	BzzCubicSpline& operator = (const BzzCubicSpline& rval);
};

//	============================================================================
//	===========================< class BzzCubicHermite >==================
//	============================================================================
class BzzCubicHermite : public BzzCubicBase
{
private:
	void GetCoefficients(void);
	void Initialize(BzzVector& xx, BzzVector& yy,
		BzzVector& yy1);
	void Initialize(BzzVector& xx, BzzVector& yy);
	BzzVector z, zx;

public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzCubicHermite(void) {};
	// copy-initializer

	BzzCubicHermite(BzzVector& xx, BzzVector& yy,
		BzzVector& yy1);
	BzzCubicHermite(BzzVector& xx, BzzVector& yy);

	void operator () (BzzVector& xx, BzzVector& yy,
		BzzVector& yy1);
	void operator () (BzzVector& xx, BzzVector& yy);

	double operator () (double xx) { return Get(xx); }
	double operator () (double xx, char ider) { return Get(xx, ider); };
	void GetMiddlePoints(BzzVector* zz);
	void GetFirstDerivativeInMiddlePoints(BzzVector* d2y);
	void GetValueInMiddlePoints(BzzVector* yy);
	void GetSecondDerivativeInSupportPoints(BzzVector* d2y);
	void GetSecondDerivativeInMiddlePoints(BzzVector* d2y);

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzCubicHermite(void) {};

	//	=========================================================================
//	*************************< assignment operator >****************************
//	============================================================================
	BzzCubicHermite& operator = (const BzzCubicHermite& r);
};

//	============================================================================
//	===========================< class BzzCubicSmooth >===================
//	============================================================================
class BzzCubicSmooth : public BzzCubicBase
{
private:
	double smooth; // smmoth = 0. linear smmoth = 1. CubicHermite
	void GetCoefficients(void);

public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzCubicSmooth(void) {};
	// copy-initializer

	BzzCubicSmooth(BzzVector& xx, BzzVector& yy);
	BzzCubicSmooth(BzzVector& xx, BzzVector& yy, double sm);

	void operator () (BzzVector& xx, BzzVector& yy);
	void operator () (BzzVector& xx, BzzVector& yy, double sm);

	double operator () (double xx) { return Get(xx); }
	double operator () (double xx, char ider) { return Get(xx, ider); };

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzCubicSmooth(void) {};

	//	============================================================================
	//	*************************< assignment operator >****************************
	//	============================================================================
	BzzCubicSmooth& operator = (const BzzCubicSmooth& r);
};

//	============================================================================
//	========< class BzzFiniteElementPointsForHermitePolynomial >================
//	============================================================================
class BzzFiniteElementPointsForHermitePolynomial : public  BzzBaseClass
{
	friend class BzzPiecewiseHermitePolynomialInterpolation;
private:
	int	polynomiumOrder,
		numExteriorPoints,
		numElements,
		numInteriorPointsInEachElement,
		numInteriorPoints,
		numSupportPointsForValues,
		numSupportPointsForDerivatives,
		numSupportPointsInEachElement,
		numSupportPoints,
		numConnectionPoints;

	BzzVector	supportPointsForValues,
		supportPointsForDerivatives,
		supportPoints,
		deltaExteriorPoints, uDelta, uDelta2,
		aux, auxx,
		lambda, tau;

	BzzMatrix	A, B, C;
	//	BzzMatrix W; // credo di no
	//	BzzFactorizedGauss G; // credo di no

public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzFiniteElementPointsForHermitePolynomial(void);
	// copy-initializer
	void operator()(int pOrder, int nElements, double xA, double xB);
	void SetExteriorPoints(int pOrder, BzzVector& x);
	void operator()(int pOrder, BzzVector& x);
	void operator()(int pOrder, BzzVector& xx, BzzVectorInt& ni);

	//	============================================================================
	//	*********************************< Functios >*******************************
	//	============================================================================
	int GetPolynomiumOrder(void) { return polynomiumOrder; }
	int GetNumExteriorPoints(void) { return numExteriorPoints; }
	int GetNumElements(void) { return numElements; }
	int GetNumInteriorPointsInEachElement(void) { return numInteriorPointsInEachElement; }
	int GetNumInteriorPoints(void) { return numInteriorPoints; }
	int GetNumSupportPointsForValues(void) { return numSupportPointsForValues; }
	int GetNumSupportPointsForDerivatives(void) { return numSupportPointsForDerivatives; }
	int GetNumSupportPointsInEachElement(void) { return numSupportPointsInEachElement; }
	int GetNumSupportPoints(void) { return numSupportPoints; }
	int GetNumConnectionPoints(void) { return numConnectionPoints; }

	//	void GetCoefficients(BzzVector *bb);
	//	void SetExteriorPoints(int n,BzzVector &x);
	void GetSupportPointsForValues(BzzVector* zz)
	{
		*zz = supportPointsForValues;
	}
	void GetSupportPointsForDerivatives(BzzVector* zz)
	{
		*zz = supportPointsForDerivatives;
	}
	void GetDeltaExteriorPoints(BzzVector* dd)
	{
		*dd = deltaExteriorPoints;
	}
	void GetSupportPoints(BzzVector* zz)
	{
		*zz = supportPoints;
	}

	//	void GetValueInSupportPoints(BzzVector *yy);
	//	void GetFirstDerivativeInSupportPoints(BzzVector *dd1);
	//	void GetSecondDerivativeInSupportPoints(BzzVector *dd2);

	void GetConnectionPoints(BzzVector* xx);

	//	void GetValueInConnectionPoints(BzzVector *yy);
	//	void GetFirstDerivativeInConnectionPoints(BzzVector *dd1);
	//	void GetSecondDerivativeInConnectionPoints(BzzVector *dd2);

	//	void GetValueInMiddlePoints(BzzVector *yy);
	//	void GetSecondDerivativeInMiddlePoints(BzzVector *d2y);
	// double operator()(double xx,int ii);

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzFiniteElementPointsForHermitePolynomial(void) {};

	//	=========================================================================
	//	*************************< assignment operator >****************************
	//	============================================================================
	//BzzFiniteElementPointsForHermitePolynomial &operator = (const BzzFiniteElementPointsForHermitePolynomial &r);

	virtual void ObjectBzzPrint(void) {}
};

//	============================================================================
//	==========================< class BzzFiniteElementPoints >==================
//	============================================================================
class BzzFiniteElementPoints : public  BzzBaseClass
{
	friend class BzzPiecewisePolynomialInterpolation;
private:
	int	polynomiumOrder,
		numExteriorPoints,
		numElements,
		numInteriorPointsInEachElement,
		numInteriorPoints,
		numPointsInEachElement,
		numSupportPointsForValues,
		numSupportPointsForDerivatives,
		numSupportPointsInEachElement,
		numSupportPoints,
		numConnectionPoints,
		twoOverlappedPoints;

	BzzVector	supportPointsForValues,
		supportPointsForDerivatives,
		exteriorPoints,
		supportPoints,
		supportPointsLegendre,
		connectionPoints,
		deltaExteriorPoints, uDelta, uDelta2,
		aux, auxx,
		lambda, tau;

	BzzMatrix	A, B, C;
	//	BzzMatrix W; // credo di no
	//	BzzFactorizedGauss G; // credo di no

public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzFiniteElementPoints(void);
	// copy-initializer
	void SetTwoOverlappedPoints(void) { twoOverlappedPoints = 1; }
	/////////////////
	// One point
	void SetExteriorPoints(int pOrder, BzzVector& x);
	// Automatic selection
	void operator()(int pOrder, int nElements, double xA, double xB);
	// Imposed selection
	void operator()(int pOrder, BzzVector& x);
	// Imposed variable selection
	void operator()(int pOrder, BzzVector& xx, BzzVectorInt& ni);
	/////////////////
	// Two points
	void SetExteriorPoints(int pOrder, int nElements, double xA, double xB);

	//	============================================================================
	//	*********************************< Functios >*******************************
	//	============================================================================
	int GetTwoOverlappedPoints(void) { return twoOverlappedPoints; }
	int GetPolynomiumOrder(void) { return polynomiumOrder; }
	int GetNumExteriorPoints(void) { return numExteriorPoints; }
	int GetNumElements(void) { return numElements; }
	int GetNumConnectionPoints(void) { return numConnectionPoints; }
	int GetNumInteriorPointsInEachElement(void) { return numInteriorPointsInEachElement; }
	int GetNumPointsInEachElement(void) { return numPointsInEachElement; }
	int GetNumSupportPointsForValues(void) { return numSupportPointsForValues; }
	int GetNumSupportPoints(void) { return numSupportPoints; }
	int GetNumInteriorPoints(void) { return numInteriorPoints; }
	int GetNumSupportPointsInEachElement(void) { return numSupportPointsInEachElement; }

	//	void GetCoefficients(BzzVector *bb);
	//	void SetExteriorPoints(int n,BzzVector &x);
	int GetNumSupportPointsForDerivatives(void) { return numSupportPointsForDerivatives; }
	void GetSupportPointsForValues(BzzVector* zz)
	{
		*zz = supportPointsForValues;
	}
	void GetExteriorPoints(BzzVector* zz)
	{
		*zz = exteriorPoints;
	}
	void GetDeltaExteriorPoints(BzzVector* dd)
	{
		*dd = deltaExteriorPoints;
	}
	void GetSupportPoints(BzzVector* zz)
	{
		*zz = supportPoints;
	}
	void GetSupportPointsLegendre(BzzVector* zz)
	{
		*zz = supportPointsLegendre;
	}
	void GetConnectionPoints(BzzVector* zz)
	{
		*zz = connectionPoints;
	}

	// Eliminare
	void GetSupportPointsForDerivatives(BzzVector* zz)
	{
		*zz = supportPointsForDerivatives;
	}

	//	void GetValueInSupportPoints(BzzVector *yy);
	//	void GetFirstDerivativeInSupportPoints(BzzVector *dd1);
	//	void GetSecondDerivativeInSupportPoints(BzzVector *dd2);

	//	void GetValueInConnectionPoints(BzzVector *yy);
	//	void GetFirstDerivativeInConnectionPoints(BzzVector *dd1);
	//	void GetSecondDerivativeInConnectionPoints(BzzVector *dd2);

	//	void GetValueInMiddlePoints(BzzVector *yy);
	//	void GetSecondDerivativeInMiddlePoints(BzzVector *d2y);
	// double operator()(double xx,int ii);

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzFiniteElementPoints(void) {};

	//	=========================================================================
	//	*************************< assignment operator >****************************
	//	============================================================================
	//BzzFiniteElementPoints &operator = (const BzzFiniteElementPoints &r);

	virtual void ObjectBzzPrint(void) {}
};

//	============================================================================
//	===========< class BzzPiecewiseHermitePolynomialInterpolation >==============
//	============================================================================
class BzzPiecewiseHermitePolynomialInterpolation : public  BzzBaseClass
{
	friend class BzzFiniteElementPointsForHermitePolynomial;
private:

	//	int	polynomiumOrder,
	//			numExteriorPoints,
	//			numElements,
	//			numInteriorPointsInEachElement,
	//			numInteriorPoints,
	//			numSupportPointsForValues,
	//			numSupportPointsForDerivatives,
	//			numSupportPointsInEachElement,
	//			numSupportPoints;

	//	BzzVector	supportPointsForValues,
	//					supportPointsForDerivatives,
	//					supportPoints,
	//					deltaExteriorPoints,uDelta,uDelta2;

	BzzFiniteElementPointsForHermitePolynomial fep;

	BzzVector	y, d1,
		aux, auxx,
		lambda, tau,
		a, b, c, d, e, f, g,
		c2, d3, e4, f5, g6, d6, e12, f20, g30, w; // numElements

//	BzzMatrix	A,B,C;
	BzzMatrix W;
	BzzFactorizedGauss G;

public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzPiecewiseHermitePolynomialInterpolation(void);
	// copy-initializer
	void operator()(BzzFiniteElementPointsForHermitePolynomial& fep, BzzVector& yy, BzzVector& dd, int coefficients = 0);

	//	BzzPiecewiseHermitePolynomialInterpolation(BzzVector &xx,BzzVector &yy);

	//	void operator () (BzzVector &xx,BzzVector &yy,
	//			BzzVector &yy1);
	//	void operator () (BzzVector &xx,BzzVector &yy);

	//	double operator () (double xx){return Get(xx);}
	//	double operator () (double xx,char ider){return Get(xx,ider);};

	//	============================================================================
	//	*********************************< Functions >*******************************
	//	============================================================================

	//	void GetCoefficients(BzzVector *bb);
	//	void SetExteriorPoints(int n,BzzVector &x);
	//	void GetSupportPointsForValues(BzzVector *zz)
	//		{*zz = supportPointsForValues;}
	//	void GetSupportPointsForDerivatives(BzzVector *zz)
	//		{*zz = supportPointsForDerivatives;}
	//	void GetDeltaExteriorPoints(BzzVector *dd)
	//		{*dd = deltaExteriorPoints;}

	//	void GetSupportPoints(BzzVector *zz)
	//		{* zz = supportPoints;}
	void GetValueInSupportPoints(BzzFiniteElementPointsForHermitePolynomial& fep, BzzVector* yy);
	void GetFirstDerivativeInSupportPoints(BzzFiniteElementPointsForHermitePolynomial& fep, BzzVector* dd1);
	void GetSecondDerivativeInSupportPoints(BzzFiniteElementPointsForHermitePolynomial& fep, BzzVector* dd2);

	//	void GetConnectionPoints(BzzVector *xx);
	void GetValueInConnectionPoints(BzzFiniteElementPointsForHermitePolynomial& fep, BzzVector* yy);
	void GetFirstDerivativeInConnectionPoints(BzzFiniteElementPointsForHermitePolynomial& fep, BzzVector* dd1);
	void GetSecondDerivativeInConnectionPoints(BzzFiniteElementPointsForHermitePolynomial& fep, BzzVector* dd2);

	void GetSecondDerivativeInLeftPoint(BzzFiniteElementPointsForHermitePolynomial& fep, double* dd2);
	void GetSecondDerivativeInRightPoint(BzzFiniteElementPointsForHermitePolynomial& fep, double* dd2);
	//	void GetValueInMiddlePoints(BzzVector *yy);
	//	void GetSecondDerivativeInMiddlePoints(BzzVector *d2y);
	double operator()(BzzFiniteElementPointsForHermitePolynomial& fep, double xx, int ii = 0);
	void operator()(BzzFiniteElementPointsForHermitePolynomial& fep, double xx, double* yy, double* dd);
	void operator()(BzzFiniteElementPointsForHermitePolynomial& fep, double xx, double* yy, double* dd, double* d2);

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzPiecewiseHermitePolynomialInterpolation(void) {};

	//	=========================================================================
	//	*************************< assignment operator >****************************
	//	============================================================================
	//BzzPiecewiseHermitePolynomialInterpolation &operator = (const BzzPiecewiseHermitePolynomialInterpolation &r);

	virtual void ObjectBzzPrint(void) {}
};

//	============================================================================
//	==================< class BzzPiecewisePolynomialInterpolation >==============
//	============================================================================
class BzzPiecewisePolynomialInterpolation : public  BzzBaseClass
{
	friend class BzzFiniteElementPoints;
private:
	int twoPoint;

	//	int	polynomiumOrder,
	//			numExteriorPoints,
	//			numElements,
	//			numInteriorPointsInEachElement,
	//			numInteriorPoints,
	//			numSupportPointsForValues,
	//			numSupportPointsForDerivatives,
	//			numSupportPointsInEachElement,
	//			numSupportPoints;

	//	BzzVector	supportPointsForValues,
	//					supportPointsForDerivatives,
	//					supportPoints,
	//					deltaExteriorPoints,uDelta,uDelta2;

	//	BzzFiniteElementPoints fep;

	BzzVector	y, d1, d2;
	//					aux,auxx,
	//					lambda,tau,
	//					a,b,c,d,e,f,g,
	//					c2,d3,e4,f5,g6,d6,e12,f20,g30,w; // numElements

	//	BzzMatrix	A,B,C;
	BzzMatrix Y, D1, D2;
	//	BzzFactorizedGauss G;

public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzPiecewisePolynomialInterpolation(void);
	// copy-initializer
	void operator()(BzzFiniteElementPoints& fep, BzzVector& yy, int coefficients = 0);

	//	BzzPiecewisePolynomialInterpolation(BzzVector &xx,BzzVector &yy);

	//	void operator () (BzzVector &xx,BzzVector &yy,
	//			BzzVector &yy1);
	//	void operator () (BzzVector &xx,BzzVector &yy);

	//	double operator () (double xx){return Get(xx);}
	//	double operator () (double xx,char ider){return Get(xx,ider);};

	//	============================================================================
	//	*********************************< Functions >*******************************
	//	============================================================================

	//	void GetCoefficients(BzzVector *bb);
	//	void SetExteriorPoints(int n,BzzVector &x);
	//	void GetSupportPointsForValues(BzzVector *zz)
	//		{*zz = supportPointsForValues;}
	//	void GetSupportPointsForDerivatives(BzzVector *zz)
	//		{*zz = supportPointsForDerivatives;}
	//	void GetDeltaExteriorPoints(BzzVector *dd)
	//		{*dd = deltaExteriorPoints;}

	//	void GetSupportPoints(BzzVector *zz)
	//		{* zz = supportPoints;}
	void GetFirstDerivativeInSupportPoints(BzzFiniteElementPoints& fep, BzzVector* dd1);
	void GetSecondDerivativeInSupportPoints(BzzFiniteElementPoints& fep, BzzVector* dd2);

	//	void GetConnectionPoints(BzzVector *xx);
	//	void GetValueInConnectionPoints(BzzFiniteElementPointsForHermitePolynomial &fep,BzzVector *yy);
	//	void GetFirstDerivativeInConnectionPoints(BzzFiniteElementPointsForHermitePolynomial &fep,BzzVector *dd1);
	//	void GetSecondDerivativeInConnectionPoints(BzzFiniteElementPointsForHermitePolynomial &fep,BzzVector *dd2);

	//	void GetSecondDerivativeInLeftPoint(BzzFiniteElementPointsForHermitePolynomial &fep,double *dd2);
	//	void GetSecondDerivativeInRightPoint(BzzFiniteElementPointsForHermitePolynomial &fep,double *dd2);
	//	void GetValueInMiddlePoints(BzzVector *yy);
	//	void GetSecondDerivativeInMiddlePoints(BzzVector *d2y);
	//	double operator()(BzzFiniteElementPointsForHermitePolynomial &fep,double xx,int ii = 0);
	//	void operator()(BzzFiniteElementPointsForHermitePolynomial &fep,double xx,double *yy,double *dd);
	//	void operator()(BzzFiniteElementPointsForHermitePolynomial &fep,double xx,double *yy,double *dd,double *d2);

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzPiecewisePolynomialInterpolation(void) {};

	//	=========================================================================
	//	*************************< assignment operator >****************************
	//	============================================================================
	//BzzPiecewisePolynomialInterpolation &operator = (const BzzPiecewisePolynomialInterpolation &r);

	virtual void ObjectBzzPrint(void) {}
};

//	============================================================================
//	==================< class BzzPiecewiseLinearInterpolation >==================
//	============================================================================
class BzzPiecewiseLinearInterpolation : public  BzzBaseClass
{
private:
	int	numElements,
		numSupportPoints,
		numInterpolators;
	BzzVector	supportPoints,
		y, b, dx;
	BzzMatrix Y, B;
	//	BzzVectorInt envelope;
public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzPiecewiseLinearInterpolation(void)
	{
		numElements = numSupportPoints = 0;
	}
	// copy-initializer
	void operator()(BzzVector& x, BzzVector& yy);
	double operator()(double xx, int i = 0);
	void operator()(double xx, double* v, double* d);

	void operator()(BzzVector* x, BzzMatrix* YY);
	void GetValues(double xx, BzzVector* yy);
	void GetFirstDerivatives(double xx, BzzVector* dd);
	void operator()(double xx, BzzVector* v, BzzVector* d);

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzPiecewiseLinearInterpolation(void) {};

	//	============================================================================
	//	*******************************< Functions >********************************
	//	============================================================================
	void GetEnvelope(void);
	void GetBvpAbscissas(int numA, BzzVector* ab);

	//	============================================================================
	//	***********************************< Print >********************************
	//	============================================================================
	virtual void ObjectBzzPrint(void);
};

//	============================================================================
//	=====================< class BzzPiecewiseInterpolation >=====================
//	============================================================================
class BzzPiecewiseInterpolation : public  BzzBaseClass
{
private:

	int	numElements,
		numSupportPoints;
	BzzVector supportPoints, y;

public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzPiecewiseInterpolation(void)
	{
		numElements = numSupportPoints = 0;
	}
	// copy-initializer
	void operator()(BzzVector& x, BzzVector& yy);
	double operator()(double xx);

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzPiecewiseInterpolation(void) {};

	//	============================================================================
	//	***********************************< Print >********************************
	//	============================================================================
	virtual void ObjectBzzPrint(void) {}
};

#endif // BZZ_CUBIC_DOUBLE_HPP