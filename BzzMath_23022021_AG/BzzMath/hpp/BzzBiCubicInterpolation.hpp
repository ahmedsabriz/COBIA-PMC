// BZZMATH: Release 7.0

//	=======================< BzzBiCubicInterpolation.cpp >======================
//	* classes BzzBiCubicBase, BzzBiCubicSpline, BzzBiCubicHermite, 				*
//	*			BzzBiCubicSmooth																	*
// * Examples: BzzMath\Examples\BzzMathBasic\BiCubicInterpolation\				*
// *				BiCubicHermite\BiCubicHermite.cpp										*
// * Examples: BzzMath\Examples\BzzMathBasic\BiCubicInterpolation					*
// *				\CubicSmooth\BiCubicSmooth.cpp											*
// * Examples: BzzMath\Examples\BzzMathBasic\BiCubicInterpolation\				*
// *				BiCubicSpline\BiCubicSpline.cpp											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	10-2003	Date Written.

////////////////// Release 5.0

#ifndef BZZ_BICUBIC_DOUBLE_HPP
#define BZZ_BICUBIC_DOUBLE_HPP

//	============================================================================
//	=======================< class BzzBiCubicBase >=============================
//	============================================================================

class BzzBiCubicBase
{
protected:
	static const char* const BZZ_ERROR;
	int numRows, numColumns;
	BzzVector	x1,
		x2,
		aux;

	BzzMatrix Y;

	BzzVectorArray a1, Y1;

	double x1Min, x2Min, x1Max, x2Max;

	//	void Initialize(BzzVector &xx,BzzVector &yy);

public:
	//	============================================================================
	//	**************************< constructors >**********************************
	//	============================================================================

		// default
	BzzBiCubicBase(void)
	{
		numRows = 0;numColumns = 0;
	}
	// copy-initializer

//	============================================================================
//	**************************< destructor >************************************
//	============================================================================

	~BzzBiCubicBase(void) {};

	//	============================================================================
	//	*******************************< Functions >********************************
	//	============================================================================
};

//	============================================================================
//	=========================< class BzzBiCubicSpline >===========================
//	============================================================================
class BzzBiCubicSpline : public BzzBiCubicBase
{
private:
	enum DataTypeSpline
	{
		MATRIX,
		VECTOR_ARRAY
	}dataType;

	BzzCubicSpline* s1, s2;

public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzBiCubicSpline(void) { s1 = 0; }
	// copy-initializer
	// NATURAL, PERIODIC,SECOND_DERIVATIVE_ADJACENT(alfa = beta = 1)
	BzzBiCubicSpline(BzzVector& xx1, BzzVector& xx2, BzzMatrix& YY);
	void operator () (BzzVector& xx1, BzzVector& xx2, BzzMatrix& YY);

	BzzBiCubicSpline(BzzVectorArray& xx1, BzzVector& xx2, BzzVectorArray& YY);
	void operator () (BzzVectorArray& xx1, BzzVector& xx2, BzzVectorArray& YY);

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzBiCubicSpline(void);

	//	============================================================================
	//	****************************< operator() > *********************************
	//	============================================================================
	double operator() (double x1, double x2);
	void Stretch(int rows, int columns, BzzVector* z1,
		BzzVector* z2, BzzMatrix* Z);
};

//	============================================================================
//	===========================< class BzzBiCubicHermite >========================
//	============================================================================
class BzzBiCubicHermite : public BzzBiCubicBase
{
private:
	enum DataTypeHermite
	{
		MATRIX,
		VECTOR_ARRAY
	}dataType;

	BzzCubicHermite* s1, s2;

public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzBiCubicHermite(void) { s1 = 0; }
	// copy-initializer
	// NATURAL, PERIODIC,SECOND_DERIVATIVE_ADJACENT(alfa = beta = 1)
	BzzBiCubicHermite(BzzVector& xx1, BzzVector& xx2, BzzMatrix& YY);
	void operator () (BzzVector& xx1, BzzVector& xx2, BzzMatrix& YY);

	BzzBiCubicHermite(BzzVectorArray& xx1, BzzVector& xx2, BzzVectorArray& YY);
	void operator () (BzzVectorArray& xx1, BzzVector& xx2, BzzVectorArray& YY);

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzBiCubicHermite(void);

	//	============================================================================
	//	****************************< operator() > *********************************
	//	============================================================================
	double operator() (double x1, double x2);
	void Stretch(int rows, int columns, BzzVector* z1,
		BzzVector* z2, BzzMatrix* Z);
};

//	============================================================================
//	===========================< class BzzBiCubicSmooth >========================
//	============================================================================
class BzzBiCubicSmooth : public BzzBiCubicBase
{
private:
	enum DataTypeSmooth
	{
		MATRIX,
		VECTOR_ARRAY
	}dataType;

	BzzCubicSmooth* s1, s2;
	double smooth; // smmoth = 0. linear smmoth = 1. CubicHermite

public:
	//	============================================================================
	//	*****************************< constructors > ******************************
	//	============================================================================
		// default
	BzzBiCubicSmooth(void) { s1 = 0; }
	// copy-initializer

	BzzBiCubicSmooth(BzzVector& xx1, BzzVector& xx2, BzzMatrix& YY);
	void operator () (BzzVector& xx1, BzzVector& xx2, BzzMatrix& YY);
	BzzBiCubicSmooth(BzzVector& xx1, BzzVector& xx2, BzzMatrix& YY, double sm);
	void operator () (BzzVector& xx1, BzzVector& xx2, BzzMatrix& YY, double sm);

	BzzBiCubicSmooth(BzzVectorArray& xx1, BzzVector& xx2, BzzVectorArray& YY);
	void operator () (BzzVectorArray& xx1, BzzVector& xx2, BzzVectorArray& YY);
	BzzBiCubicSmooth(BzzVectorArray& xx1, BzzVector& xx2, BzzVectorArray& YY, double sm);
	void operator () (BzzVectorArray& xx1, BzzVector& xx2, BzzVectorArray& YY, double sm);

	//	============================================================================
	//	****************************< destructor > *********************************
	//	============================================================================
	~BzzBiCubicSmooth(void);

	//	============================================================================
	//	****************************< operator() > *********************************
	//	============================================================================
	double operator() (double x1, double x2);
	void Stretch(int rows, int columns, BzzVector* z1,
		BzzVector* z2, BzzMatrix* Z);
};

#endif // BZZ_BICUBIC_DOUBLE_HPP