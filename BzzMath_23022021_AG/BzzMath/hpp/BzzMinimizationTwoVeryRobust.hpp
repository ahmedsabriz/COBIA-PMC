// BZZMATH: Release 7.0

//	===================< BzzMinimizationTwoVeryRobust.hpp >=========================
//	* Class BzzMinimizationTwoVeryRobust for bidimensional minimization				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	09-2009	Date Written.

#ifndef MINIMIZATION_TWO_VERY_ROBUST
#define MINIMIZATION_TWO_VERY_ROBUST

//	============================================================================
//	===========< class BzzMyMinimizationTwoMonoRobustIntern >================
//	============================================================================
class BzzMyMinimizationTwoMonoRobustIntern : public BzzMyMinimizationRobustObject
{
	double (*ptrF)(BzzVector& x);
	int iMonoIntern;
	int iMonoExtern;
	BzzVector x;
public:
	virtual double GetFunctionValue(double t);
	virtual void ObjectBzzPrint(void);
	void operator()(int ivi, int ive, BzzVector& xx,
		double (*ptr)(BzzVector& x))
	{
		iMonoIntern = ivi;iMonoExtern = ive;x = xx;ptrF = ptr;
	}
	void operator()(double t)
	{
		x[iMonoExtern] = t;
	}
};

//	============================================================================
//	===========< class BzzMyMinimizationTwoMonoRobustExtern >================
//	============================================================================
//      BzzMinimizationMonoVeryRobustObject
class BzzMyMinimizationTwoMonoRobustExtern : public BzzMyMinimizationRobustObject
{
	int iMonoExtern;
	double tMinIntern, tMaxIntern;
	BzzMyMinimizationTwoMonoRobustIntern myIntern;
	BzzMinimizationMonoVeryRobustObject monoRobustInternObject;
public:
	int iterF, iterT;
	virtual double GetFunctionValue(double t);
	virtual void ObjectBzzPrint(void);
	BzzMyMinimizationTwoMonoRobustExtern(int ivi, int ive, BzzVector& xx,
		double (*ptr)(BzzVector& x),
		double tMi, double tMa)
	{
		iterF = iterT = 1;
		tMinIntern = tMi;
		tMaxIntern = tMa;
		iMonoExtern = ive;
		myIntern(ivi, ive, xx, ptr);
	}
};

//	============================================================================
//	===================< class BzzMinimizationTwoVeryRobust >===============
//	============================================================================
//	============================================================================
//	************************< Twodimensional Minimum >************************
//	============================================================================
//	*************************< Unconstrained Minimization >*********************
//	* Min F = F(x)																					*
//	* x >= L (lower bounds)																		*
//	* x <= U (upper bounds)																		*

class BzzMinimizationTwoVeryRobust : public BzzBaseClass
{
	friend class BzzMyMinimizationMultiMonoNonRobustIntern;
	friend class BzzMyMinimizationMultiMonoNonRobustExtern;
	friend class BzzMyMinimizationMultiMonoRobustIntern;
	friend class BzzMyMinimizationMultiMonoRobustExtern;
private:
	static const char* const BZZ_ERROR;
	static const int MAX_ITER;
	static const double TOL_REL;
	static const double TOL_ABS;
	int	numVariables,
		iter,
		iterFeasible,
		iterTotal,
		iterTotalFeasible,
		nIter,
		maxIter,
		numGlobalUnfesible,
		iMonoIntern, iMonoExtern,
		//			plotTask,
		printTasks,
		printSubTasks,
		extraPointsNumber,
		plot;

	char control;
	double	f0,
		fSolution;

	BzzVector	x0,
		x,
		xSolution,
		xL,
		xU;

	double (*ptrF)(BzzVector& x);

public:
	void SetExtraPointsNumberForMultipleMinimaSearch(int extra)
	{
		if (extra > 10)
			extraPointsNumber = extra;
	}
	void SetPlot(void)
	{
		plot = 1;
	}
	//	============================================================================
	//	********************************< constructors >****************************
	//	============================================================================
		// default constructor
	BzzMinimizationTwoVeryRobust(void);

	// copy constructor
	BzzMinimizationTwoVeryRobust(const BzzMinimizationTwoVeryRobust& rval);

	BzzMinimizationTwoVeryRobust(BzzVector& x00, int i, int j, double FF0,
		double (*yy)(BzzVector& x),
		BzzVector& xxL, BzzVector& xxU);
	void operator()(BzzVector& x00, int i, int j, double FF0,
		double (*yy)(BzzVector& x),
		BzzVector& xxL, BzzVector& xxU);

	BzzMinimizationTwoVeryRobust(BzzVector& x00, double FF0,
		double (*yy)(BzzVector& x),
		BzzVector& xxL, BzzVector& xxU);
	void operator()(BzzVector& x00, double FF0,
		double (*yy)(BzzVector& x),
		BzzVector& xxL, BzzVector& xxU);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzMinimizationTwoVeryRobust(void)
	{}

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
	// return control < 0 Solution not yet founded with the required precision
	//			(itertions >= required max iteration)
	// return control = 0 Founded solution with a value < to setTollF required
	// return control > 0 Founded solutions with the required precision
	char operator() (void);

	//	void SetTolRel(double tr){control = -1; tTollRel = fabs(tr);}
	//	void SetTolAbs(double ta){control = -1; tTollAbs = fabs(ta);}
	//	void SetTolF(double yt = TOL_ABS_Y){control = -1; yToll = yt;}
	void SetMaxIter(int mxit) { maxIter = mxit; }

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
	//	void SetPlot(void){plotTask = 1;}

	//	============================================================================
	//	*******************************< Non-modifying functions >******************
	//	============================================================================
	virtual void ObjectBzzPrint(void);
};

#endif // MINIMIZATION_TWO_VERY_ROBUST