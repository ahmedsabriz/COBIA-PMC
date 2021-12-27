// BZZMATH: Release 7.0

//	============================================================================
//	=============< class BzzMinimizationMultiVeryRobustObject >=================
//	============================================================================
//	* Min F = F(x)																					*
//	* x >= L (lower bounds)																		*
//	* x <= U (upper bounds)																		*
// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-2012

#ifndef MINIMIZATION_MULTI_VERY_ROBUST_OBJECT
#define MINIMIZATION_MULTI_VERY_ROBUST_OBJECT

//=============================================================================================
class BzzMyMinimizationVeryRobustObject : public BzzMyMinimizationRobustObject
{
public:
	BzzVector x;
	int iVar;
	virtual double GetFunctionValue(BzzVector& x) { return 0.; }
	virtual double GetFunctionValue(double t);//{return 0.;}
//		{
//		x[iVar] = t;
//		printf("\nt %e",t);getchar();
//		return GetFunctionValue(x);
//		}
	virtual void ObjectBzzPrint(void) {};
};

//	=============< class BzzMinimizationMultiVeryRobustObject >=================

class BzzMinimizationMultiVeryRobustObject : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static const int MAX_ITER;
	static const double X_TOL_REL;
	static const double X_TOL_ABS;
	static const double F_TOL_REL;
	static const double F_TOL_ABS;
	int	numVariables,
		iter,
		iterFeasible,
		iterTotal,
		iterTotalFeasible,
		nIter,
		maxIter,
		printTasks,
		printSubTasks,
		numThreads,
		numT;

	double xTolRel, xTolAbs, fTolRel, fTolAbs;

	char control;
	double	f0,
		fSolution;

	BzzVector	x00,
		xSolution,
		xL,
		xU,
		dx;

	BzzMinimizationMonoVeryRobustObject* mono;
	BzzMinimizationRobustObject* multi;
	BzzVector* x0;
	double* y0;
	//	double (*ptrFun)(BzzVector &x);
	BzzMyMinimizationVeryRobustObject* ptrObject;

public:

	//	============================================================================
	//	********************************< constructors >****************************
	//	============================================================================
		// default constructor
	BzzMinimizationMultiVeryRobustObject(void);

	// copy constructor
	BzzMinimizationMultiVeryRobustObject(const BzzMinimizationTwoVeryRobust& rval);

	void operator()(BzzVector& x00,
		//		double (*yy)(BzzVector &x),
		BzzMyMinimizationVeryRobustObject* ptro,
		BzzVector& xxL, BzzVector& xxU);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzMinimizationMultiVeryRobustObject(void);

	//	============================================================================
	//	******************< Non-modifying access functions >************************
	//	============================================================================
	double GetXSolution(BzzVector* xS)
	{
		*xS = xSolution;return fSolution;
	}
	void GetSolutions(BzzMatrix* Xs, BzzVector* fS);

	int GetIterTotalFeasible(void)
	{
		return iterTotalFeasible;
	}
	int GetIterTotal(void)
	{
		return iterTotal;
	}

	//	============================================================================
	//	****************************< Modifying access functions >******************
	//	============================================================================
	// return control = 0 iterations > maxIter
	// return control = 1 Solution found solutions with the required precision
	char operator() (void);

	void SetTolRelX(double tr) { xTolRel = tr; } // default = 1.e-3
	void SetTolAbsX(double ta) { xTolAbs = ta; } // default = 1.e-3
	void SetTolRelF(double fr) { fTolRel = fr; } // default = .9;
	void SetTolAbsF(double fa) { fTolAbs = fa; } // default = .01;
	void SetMaxIter(int mxit) { maxIter = mxit; } // default = 500

	void StopTaskPrint(void) { printTasks = 0; }
	void StopSubTasksPrint(void) { printSubTasks = 0; }
	void SetTasksPrint(void) { printTasks = 1; }
	void SetSubTasksPrint(int psb = 1)
	{
		if (psb > 0)
			printSubTasks = psb;
		else
			printSubTasks = 1;
	}

	//	============================================================================
	//	*******************************< Non-modifying functions >******************
	//	============================================================================
	virtual void ObjectBzzPrint(void);
};

#endif // MINIMIZATION_MULTI_VERY_ROBUST_OBJECT