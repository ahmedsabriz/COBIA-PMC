// BZZMATH: Release 7.0

//	===================< BzzMinimizationMonoObject.hpp >========================
//	* Class BzzMinimizationMonoObject for monodimensional minimization			*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 21, 22)			*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	04-2006	Date Written.

#ifndef MINIMUM_DOUBLE_MONO_OBJECT_HPP
#define MINIMUM_DOUBLE_MONO_OBJECT_HPP

//	============================================================================
//	====================< class BzzMyMinimizationMonoObject >===================
//	============================================================================
class BzzMyMinimizationMonoObject : public BzzBaseClass
{
public:
	virtual double GetFunctionValue(double t) { return 0.; }
	virtual void ObjectBzzPrint(void) {};
};

//	============================================================================
//	===================< class BzzMinimizationMonoObject >================
//	============================================================================

class BzzMinimizationMonoObject : public BzzBaseClass
{
private:
	enum MethodBzzMinimizationMonoObject
	{
		PARABOLIC_INTERPOLATION,
		GOLDEN_SECTION
	};

	enum BracketingStatus
	{
		BRACKET,
		NO_BRACKET
	}bracketingStatus;

	static const char* const BZZ_ERROR;
	static const int ITMAX, NPMAX;
	static const double GOLDEN;
	static const double GOLDEN_MINUS_ONE;
	static const double TOL_REL;
	static const double TOL_ABS;
	static const double TOL_ABS_F;

	BzzMyMinimizationMonoObject* ptrObject;

	char control,
		out,
		speed,
		fifht,
		stop;
	double lamda, dt21, df21, dt31, df31, dt32, tFifth;

	double tLeft, fLeft,
		tRight, fRight,
		tSolution, fSolution, // the best
		t2, f2, // the second one
		t3, f3, // the third
		tNew, fNew,
		tExt, fExt,
		tTollRel, tTollAbs, fToll,
		delta,
		dxm;

	int iter,
		iterTotal,
		niter,
		numPoints,
		numParabolic,
		iParabolic,
		iterationStatus;

	MethodBzzMinimizationMonoObject methodBzzMinimization,
		methodBzzMinimizationUsed;

	void Initialize(void);
	void Deinitialize(void);

	char SearchBzzMinimization(void);
	char BracketBzzMinimization(void);
	char ControlBzzMinimization(void);
	char ControlNewBzzMinimization(void);
	char MinimumMonoSolve(int ni);

public:
	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================
		// default constructor
	BzzMinimizationMonoObject(void);

	// copy constructor
	BzzMinimizationMonoObject(const BzzMinimizationMonoObject& rval);

	// first constructor BzzMinimizationMonoObject z(tL,tR,fmin);
	BzzMinimizationMonoObject(double tL, double tR, BzzMyMinimizationMonoObject* ptrO);
	void operator()(double tL, double tR, BzzMyMinimizationMonoObject* ptrO);

	// second constructor BzzMinimizationMonoObject z(tL,tR,fmin,yL,yR);
	BzzMinimizationMonoObject(double tL, double tR, BzzMyMinimizationMonoObject* ptrO,
		double yL, double yR);
	void operator()(double tL, double tR, BzzMyMinimizationMonoObject* ptrO,
		double yL, double yR);

	// third constructor BzzMinimizationMonoObject z(t0,fmin);
	BzzMinimizationMonoObject(double t0, BzzMyMinimizationMonoObject* ptrO);
	void operator()(double t0, BzzMyMinimizationMonoObject* ptrO);

	// fourth constructor BzzMinimizationMonoObject z(t0,fmin,y0);
	BzzMinimizationMonoObject(double t0, BzzMyMinimizationMonoObject* ptrO, double y0);
	void operator()(double t0, BzzMyMinimizationMonoObject* ptrO, double y0);

	// fifth constructor BzzMinimizationMonoObject z(t0,tL,tR,fmin);
	BzzMinimizationMonoObject(double to, double tL, double tR, BzzMyMinimizationMonoObject* ptrO);
	void operator()(double t0, double tL, double tR, BzzMyMinimizationMonoObject* ptrO);

	void SetTolRel(double tr);
	void SetTolAbs(double ta);
	void SetTolF(double ft);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzMinimizationMonoObject(void) {};

	//	============================================================================
	//	******************< Non-modifying access functions >************************
	//	============================================================================
	int IterationCounter(void)
	{
		return iterTotal;
	}
	double TLeft(void)
	{
		return tLeft;
	}
	double FLeft(void)
	{
		return fLeft;
	}
	double TRight(void)
	{
		return tRight;
	}
	double FRight(void)
	{
		return fRight;
	}
	double TSolution(void)
	{
		return tSolution;
	}
	double FSolution(void)
	{
		return fSolution;
	}

	//	============================================================================
	//	****************************< Modifying access functions >******************
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
	// return = 3 fabs(fSolution) <= fTollAbs
	// return = 4 The function is constant in the three best points
	// return = 5 The number of parabolic interpolations is equal to fixed value

	void SetParabolicMaxNumber(int nn)
	{
		if (nn > 0)numParabolic = nn;
	}

	//	============================================================================
	//	*******************************< Non-modifying functions >******************
	//	============================================================================
	virtual void ObjectBzzPrint(void);
};

#endif // MINIMUM_DOUBLE_MONO_OBJECT_HPP