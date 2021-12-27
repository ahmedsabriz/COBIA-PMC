// BZZMATH: Release 7.0

//	===================< BzzMinimizationMono.hpp >========================
//	* Class BzzMinimizationMono for monodimensional minimization			*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 21, 22)			*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	07-1997	Date Written.
//	09-2009	Added BzzMinimizationMonoMP class.

#ifndef MINIMUM_DOUBLE_MONO_HPP
#define MINIMUM_DOUBLE_MONO_HPP

//	============================================================================
//	=========================< class BzzMinimizationMono >=====================
//	============================================================================

class BzzMinimizationMono : public BzzBaseClass
{
private:
	enum MethodBzzMinimizationMono
	{
		PARABOLIC_INTERPOLATION,
		GOLDEN_SECTION
	};

	enum MonoMonoOrMulti
	{
		MONO,
		MULTI_MONO,
		MULTI_PI
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

	char control,
		out,
		speed,
		stop;
	double lamda, dt21, df21, dt31, df31, dt32;

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

	double (*ptrF)(double t);
	double (*ptrMulti)(BzzVector& x);
	double F(double t);
	MonoMonoOrMulti monoMonoMulti;
	int iMono;
	BzzVector xMulti, x0, pi, xL, xU;
	void FindBounds(void);
	int numVariables;

	int iter,
		iterTotal,
		niter,
		numPoints,
		iterationStatus;

	MethodBzzMinimizationMono methodBzzMinimization,
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
	BzzMinimizationMono(void);

	// copy constructor
	BzzMinimizationMono(const BzzMinimizationMono& rval);

	// first constructor BzzMinimizationMono z(tL,tR,fmin);
	BzzMinimizationMono(double tL, double tR, double (*y)(double t));
	void operator()(double tL, double tR, double (*y)(double t));
	BzzMinimizationMono(double (*y)(double t), double tL, double tR);
	void operator()(double (*y)(double t), double tL, double tR);

	// second constructor BzzMinimizationMono z(tL,tR,fmin,yL,yR);
	BzzMinimizationMono(double tL, double tR, double (*y)(double t),
		double yL, double yR);
	void operator()(double tL, double tR, double (*y)(double t),
		double yL, double yR);

	// third constructor BzzMinimizationMono z(t0,fmin);
	BzzMinimizationMono(double t0, double (*y)(double t));
	void operator()(double t0, double (*y)(double t));

	// fourth constructor BzzMinimizationMono z(t0,fmin,y0);
	BzzMinimizationMono(double t0, double (*y)(double t), double y0);
	void operator()(double t0, double (*y)(double t), double y0);

	//	==========================< Multi Dimensional Space >========================
		// first constructor BzzMinimizationMono z(i,x0,fmin,tmin,tmax);
	BzzMinimizationMono(int i, BzzVector& x00, double (*fMulti)(BzzVector& x),
		double tL, double tR);
	void operator()(int i, BzzVector& x00, double (*fMulti)(BzzVector& x),
		double tL, double tR);

	// second constructor BzzMinimizationMono z(i,x0,F0,fmin,tmin,tmax);
	BzzMinimizationMono(int i, BzzVector& x00, double F00,
		double (*fMulti)(BzzVector& x), double tL, double tR);
	void operator()(int i, BzzVector& x00, double F00,
		double (*fMulti)(BzzVector& x), double tL, double tR);

	// third constructor BzzMinimizationMono z(i,x0,fmin,xL,xU);
	BzzMinimizationMono(int i, BzzVector& x00,
		double (*fMulti)(BzzVector& x),
		BzzVector& xL, BzzVector& xU);
	void operator()(int i, BzzVector& x00, double (*fMulti)(BzzVector& x),
		BzzVector& xL, BzzVector& xU);

	// fourth constructor BzzMinimizationMono z(i,x0,F0,fmin,tmin,tmax);
	BzzMinimizationMono(int i, BzzVector& x00, double F00,
		double (*fMulti)(BzzVector& x),
		BzzVector& xL, BzzVector& xU);
	void operator()(int i, BzzVector& x00, double F00,
		double (*fMulti)(BzzVector& x),
		BzzVector& xL, BzzVector& xU);

	// fifth constructor BzzMinimizationMono z(pi,x0,fmin,xL,xU);
	BzzMinimizationMono(BzzVector& pi, BzzVector& x00,
		double (*fMulti)(BzzVector& x),
		BzzVector& xL, BzzVector& xU);
	void operator()(BzzVector& pi, BzzVector& x00,
		double (*fMulti)(BzzVector& x),
		BzzVector& xL, BzzVector& xU);

	// fourth constructor BzzMinimizationMono z(pi,x0,F0,fmin,xL,xU);
	BzzMinimizationMono(BzzVector& pi, BzzVector& x00, double F00,
		double (*fMulti)(BzzVector& x),
		BzzVector& xL, BzzVector& xU);
	void operator()(BzzVector& pi, BzzVector& x00, double F00,
		double (*fMulti)(BzzVector& x),
		BzzVector& xL, BzzVector& xU);

	void SetTolRel(double tr);
	void SetTolAbs(double ta);
	void SetTolF(double ft);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzMinimizationMono(void) {};

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

	//	============================================================================
	//	*******************************< Non-modifying functions >******************
	//	============================================================================
	virtual void ObjectBzzPrint(void);
};

//	============================================================================
//	=========================< class BzzGetMinimizationMono >=======================
//	============================================================================

class BzzGetMinimizationMono
{
public:
	//	int monoMulti,iMono;
	//	BzzVector x0,xMulti,pi;
	double (*y)(double t);
	//	double (*ptrMulti)(BzzVector &xMulti);
	double GetFunctionValue(double t)
	{
		//		switch(monoMulti)
		//			{
		//			case 1:
		return y(t);
		//				break;
		//			case 2:
		//				xMulti[iMono] = t;
		//				return ptrMulti(xMulti);
		//				break;
		//			case 3:
		//				xMulti = x0 + t * pi;
		//				return ptrMulti(xMulti);
		//				break;
		//			}
		//		return 1.;
	}
};

//	============================================================================
//	====================< class BzzMinimizationMonoMP >=============================
//	============================================================================

class BzzMinimizationMonoMP : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static const double TOL_REL;
	static const double TOL_ABS;
	static const double TOL_ABS_Y;

	enum MonoMonoOrMulti
	{
		MONO,
		MULTI_MONO,
		MULTI_PI
	};

	int numThreads;
	BzzGetMinimizationMono* ptrThreads;
	int MaxDt(double d1, double d2, double d3);
	char	control;
	double	tLeft,
		yLeft,
		tRight,
		yRight,
		tSolution,
		ySolution,
		tTolAbs, tTolRel,
		yTolAbs;

	double (*y)(double t);

	MonoMonoOrMulti monoMonoMulti;

	int iter, iterTotal, niter, functionCallCounter;
	BzzVector tv, yv, tMP, yMP;
	BzzVectorInt iSort;

	//	void Initialize(void);
	//	void Deinitialize(void);
	char ControlBzzMinimizationMono(void);
	//	char ControlNewBzzMinimizationMono(void);
	//	char ZeroSolve(int ni);

public:
	//	============================================================================
	//	**************************< constructors >**********************************
	//	============================================================================
		// default constructor
	BzzMinimizationMonoMP(void);

	// copy constructor
	BzzMinimizationMonoMP(const BzzMinimizationMonoMP& rval);

	// FIRST
	// BzzMinimizationMonoMP z(tL,tR,yzer);
	BzzMinimizationMonoMP(double tL, double tR, double (*yy)(double t));
	void operator()(double tL, double tR, double (*yy)(double t));

	// SECOND
	// constructor BzzMinimizationMono z(tL,tR,fmin,t0,y0);
	BzzMinimizationMonoMP(double tL, double tR, double (*y)(double t),
		double t0, double y0);
	void operator()(double tL, double tR, double (*y)(double t), double t0, double y0);

	// constructor BzzMinimizationMono m(i,x0,F0,fmin,xL,xR);
//	BzzMinimizationMonoMP(int i,BzzVector &x00,double F00,
//		double (*fMulti)(BzzVector &x),
//		BzzVector &xL,BzzVector &xU);
//	void operator()(int i,BzzVector &x00,double F00,
//		double (*fMulti)(BzzVector &x),
//		BzzVector &xL,BzzVector &xU);
//	void operator()(double tL,double tR,double (*fMulti)(BzzVector &x),double t0,double y0);

	// constructor BzzMinimizationMono z(pi,x0,F0,fmin,xL,xU);
//	BzzMinimizationMonoMP(BzzVector &pi,BzzVector &x00,double F00,
//		double (*fMulti)(BzzVector &x),
//		BzzVector &xL,BzzVector &xU);
//	void operator()(BzzVector &pi,BzzVector &x00,double F00,
//		double (*fMulti)(BzzVector &x),
//		BzzVector &xL,BzzVector &xU);

//	============================================================================
//	*****************************< destructor >*********************************
//	============================================================================
	~BzzMinimizationMonoMP(void);

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
	double FSolution(void)
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

#endif // MINIMUM_DOUBLE_MONO_HPP