// BZZMATH: Release 7.0

//	======================< BzzFunctionRoot.hpp >========================
//	* Class BzzFunctionRoot for root finding										*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-1993	Date Written.
//	09-2009	Added BzzFunctionRootMP class.

#ifndef ZERO_DOUBLE_HPP
#define ZERO_DOUBLE_HPP

//	============================================================================
//	======================< class BzzFunctionRoot >=============================
//	============================================================================

class BzzFunctionRoot : public BzzBaseClass
{
private:
	enum MethodBzzFunctionRoot
	{
		INTERPOLATION,
		BISECTION
	};

	static const char* const BZZ_ERROR;
	static const int ITMAX, NPMAX;
	static const double TOL_REL;
	static const double TOL_ABS;
	static const double TOL_ABS_Y;

	char control,
		stop;
	int iterationStatus;
	double tLeft,
		yLeft,
		tRight,
		yRight,
		tNew,
		yNew,
		tSolution,
		ySolution,
		* p, * q, * prev, * dp,
		tExt, yExt,
		tTollAbs, tTollRel,
		yTollAbs;

	double (*y)(double t);

	int iter, iterTotal, niter;
	int numPoints;
	MethodBzzFunctionRoot methodBzzFunctionRoot, methodBzzFunctionRootUsed;
	char speed;

	void Initialize(void);
	void Deinitialize(void);
	char ControlBzzFunctionRoot(void);
	char ControlNewBzzFunctionRoot(void);
	char ZeroSolve(int ni);

public:
	//	============================================================================
	//	**************************< constructors >**********************************
	//	============================================================================
		// default constructor
	BzzFunctionRoot(void);

	// copy constructor
	BzzFunctionRoot(const BzzFunctionRoot& rval);

	// BzzFunctionRoot z(tL,tR,yzer);
	BzzFunctionRoot(double tL, double tR, double (*yy)(double t));
	void operator()(double tL, double tR, double (*yy)(double t));

	// BzzFunctionRoot z(tL,tR,yzer,yL,yR);
	BzzFunctionRoot(double tL, double tR, double (*yy)(double t),
		double yL, double yR);
	void operator()(double tL, double tR, double (*yy)(double t),
		double yL, double yR);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzFunctionRoot(void);

	//	============================================================================
	//	*********************< Non-modifying access functions >*********************
	//	============================================================================
	int IterationCounter(void)
	{
		return iterTotal;
	}
	double TLeft(void)
	{
		return tLeft;
	}
	double YLeft(void)
	{
		return yLeft;
	}
	double TRight(void)
	{
		return tRight;
	}
	double YRight(void)
	{
		return yRight;
	}
	double TSolution(void)
	{
		return tSolution;
	}
	double YSolution(void)
	{
		return ySolution;
	}

	//	============================================================================
	//	*******************< Modifying access functions >***************************
	//	============================================================================
	char operator()(void)
	{
		return (*this)(100);
	}
	char operator()(int ni);

	// return == -2 The function values have the same sign
	// return < 0 Search failed before reaching the solution
	// return = 0 Found problems
	// return = 1 Maximum number of iterations
	// return = 2 (tRight - tLeft) <=
	//				tTolAbs + EPS*(fabs(tLeft) + fabs(tRight))
	// return = 3 fabs(ySolution) <= yTollAbs

	void SetTolRel(double tr);
	void SetTolAbs(double ta);
	void SetTolY(double yt);

	//	============================================================================
	//	======================< Non-modifying functions >===========================
	//	============================================================================

	//	*****************************< BzzPrint >***********************************
	virtual void ObjectBzzPrint(void);
};

//	============================================================================
//	=========================< class BzzGetFunctionRoot >=======================
//	============================================================================

class BzzGetFunctionRoot
{
public:
	double (*y)(double t);
	double GetFunctionValue(double t)
	{
		return y(t);
	}
};

//	============================================================================
//	====================< class BzzFunctionRootMP >=============================
//	============================================================================

class BzzFunctionRootMP : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static const double TOL_REL;
	static const double TOL_ABS;
	static const double TOL_ABS_Y;

	int numThreads;
	BzzGetFunctionRoot* ptrThreads;

	char	control,
		stop;
	double	tLeft,
		yLeft,
		tRight,
		yRight,
		tSolution,
		ySolution,
		tTolAbs, tTolRel,
		yTolAbs;

	double (*y)(double t);

	int iter, iterTotal, niter, functionCallCounter;
	BzzVector tv, yv, tMP, yMP;
	BzzVectorInt iSort;

	//	void Initialize(void);
	//	void Deinitialize(void);
	char ControlBzzFunctionRoot(void);
	//	char ControlNewBzzFunctionRoot(void);
	//	char ZeroSolve(int ni);

public:
	//	============================================================================
	//	**************************< constructors >**********************************
	//	============================================================================
		// default constructor
	BzzFunctionRootMP(void);

	// copy constructor
	BzzFunctionRootMP(const BzzFunctionRootMP& rval);

	// BzzFunctionRootMP z(tL,tR,yzer);
	BzzFunctionRootMP(double tL, double tR, double (*yy)(double t));
	void operator()(double tL, double tR, double (*yy)(double t));

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzFunctionRootMP(void);

	//	============================================================================
	//	*********************< Non-modifying access functions >*********************
	//	============================================================================
	int IterationCounter(void)
	{
		return iterTotal;
	}
	int FunctionCallCounter(void)
	{
		return functionCallCounter;
	}
	double TLeft(void)
	{
		return tLeft;
	}
	double YLeft(void)
	{
		return yLeft;
	}
	double TRight(void)
	{
		return tRight;
	}
	double YRight(void)
	{
		return yRight;
	}
	double TSolution(void)
	{
		return tSolution;
	}
	double YSolution(void)
	{
		return ySolution;
	}

	//	============================================================================
	//	*******************< Modifying access functions >***************************
	//	============================================================================
	char operator()(void)
	{
		return (*this)(100);
	}
	char operator()(int ni);

	// return < 0 Search failed before reaching the solution
	// return = 0 Found problems
	// return = 1 Maximum number of iterations
	// return = 2 (tRight - tLeft) <=
	//				tTolAbs + EPS*(fabs(tLeft) + fabs(tRight))
	// return = 3 fabs(ySolution) <= yTolAbs

	void SetTolRel(double tr);
	void SetTolAbs(double ta);
	void SetTolY(double yt);

	//	============================================================================
	//	======================< Non-modifying functions >===========================
	//	============================================================================

	//	*****************************< BzzPrint >***********************************
	virtual void ObjectBzzPrint(void);
};

#endif // ZERO_DOUBLE_HPP