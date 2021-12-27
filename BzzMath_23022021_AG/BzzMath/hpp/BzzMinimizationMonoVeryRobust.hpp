// BZZMATH: Release 7.0

//	==================< BzzMinimizationMonoVeryRobust.hpp >=======================
//	* Class BzzMinimizationMonoVeryRobust for minimization in double precision	*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	04-2003	Date Written.
//	05-2008	Added GetUnfeasibleGlobal function.
//	05-2008	Added ResetUnfeasibleGlobal function.

#ifndef MINIMIZATION_ROBUST_MONO_HPP
#define MINIMIZATION_ROBUST_MONO_HPP

//	============================================================================
//	=========================< class BzzMinimizationMonoVeryRobust >=====================
//	============================================================================

class BzzMinimizationMonoVeryRobust : public BzzBaseClass
{
private:
	enum FunctionTypeBzzMinimization
	{
		MONO_MONO,
		MONO_MULTI_MONO,
		MONO_MULTI_MULTI,
		MONO_SCANNING,
		LINEAR_MULTI,
		UNCONSTRAINED_MULTI,
	}functionType;

	enum UnconstrainedMethodBzzMinimization
	{
		START,
		MOVE,
		RE_START,
		//		TWO_SIMPLEX,
		MULTI_SIMPLEX,
		ROB,
		ROS,
		ROD,
		EXIT
	}unconstrainedMethodBzzMinimization;

	static const char* const BZZ_ERROR;
	static const int MAX_ITER;
	static const double TOL_REL;
	static const double TOL_ABS;
	static const double TOL_REL_Y;
	static const double TOL_ABS_Y;
	static const double DT_RATIO;
	static const double DT_MIN;
	static const double U_F_RATIO;

	char	control, unimodal, unfeasibleGlobal;

	double	t0, y0,
		t1, y1,
		t2, y2,
		t3, y3,
		tMin, tMax,
		tBreak,
		tNew, yNew,
		tUnfeasible,
		tSolution, ySolution,
		lambda,
		lambdaMax,
		dt21, dy21,
		dt31, dy31,
		dt32, dy32,
		tTollRel, tTollAbs, yToll,
		unfeasibleFeasibleRatio, oneMunfeasibleFeasibleRatio,
		dtToll;

	BzzVector t, y, tCriticInf, tCriticSup;
	BzzVector tAux, yAux, tAux1;
	BzzVector x0, pi, pi1, x1, x2, x3;
	//	BzzVector piA,piB; // for TWO_TWO_MULTI_MULTI

	double (*yMono)(double t);
	double (*yMulti)(BzzVector& x);
	//	void (*fMulti)(BzzVector &x,BzzVector &f);
	double (*FMulti)(BzzVector& x);
	double Y(double t);
	//	double Phiw(BzzVector &x);

	int	iVar,
		iter,
		iterFeasible,
		iterTotal,
		iterTotalInConstructor,
		iterFeasibleInConstructor,
		iterTotalFeasible,
		iterMultiRos,
		nIter,
		maxIter,
		numPoints,
		nCritic,
		iVar1, iVar2,
		simplexIterations;

	int	printTasks,
		printSubTasks,
		plotTask;

	void InitializeMono(void);
	void InitializeScanningMono(void) { BzzError("TODO InitializeScanningMono"); }
	void Deinitialize(void);

	void FindtMax(double tt1, double ttUnfeasible);
	void FindtMin(double ttUnfeasible, double ttn);
	void FindtCritic(double tt1, double ttUnfeasible, double ttn);
	char MonoSearch(void);
	char UnconstrainedMultiSearch(void);
	void MonoPrint(void);
	void MonoPrintNle(void);
	void FindBounds(BzzVector& xx, BzzVector& pp, double* tmi, double* tma);
	void UnconstrainedMultiPrint(void);
	void UnconstrainedMultiStart(void);
	void UnconstrainedMultiSimplex(void);
	void UnconstrainedMultiMove(void);
	void UnconstrainedMultiRob(void);
	void UnconstrainedMultiRos(void);
	void UnconstrainedMultiRod(void);

	// Multidimensional
	double	f0,
		fSolution,
		fP,
		fM,
		fI,
		fR,
		fF,
		percMove,
		percDMI,
		dMI,
		dMF,
		dIF,
		tNew1;

	int	numVariables,
		numVertices,
		numVerticesMinusOne,
		oldCase,
		countLeftRight;

	BzzVector			x,
		xL, // lower bound
		xU, // upper bound
		x0M,
		xSolution,
		xP,
		xM,
		xI,
		xR,
		xF,
		weights,
		functionsx,
		//							functions0,
		//							functionsSolution,
		//							functionsP,
		//							functionsM,
		//							functionsI,
		//							functionsR,
		//							functionsF,
		aux;

	BzzVector YS;
	BzzMatrix XS;
	BzzVector tS1;
	BzzVector piOld;

	BzzVectorInt sorted;
	BzzMinimizationMonoVeryRobust* mono;

public:
	//	============================================================================
	//	********************************< constructors >****************************
	//	============================================================================
		// default constructor
	BzzMinimizationMonoVeryRobust(void);

	// copy constructor
	BzzMinimizationMonoVeryRobust(const BzzMinimizationMonoVeryRobust& rval);

	//	============================================================================
	//	************************< Monodimensional Minimum >*************************
	//	============================================================================
		//	MONO_MONO
	BzzMinimizationMonoVeryRobust(double tt0, double yy0,
		double (*y)(double t),
		double tmin, double tmax);
	void operator()(double tt0, double yy0,
		double (*y)(double t),
		double tmin, double tmax);

	BzzMinimizationMonoVeryRobust(double tt0,
		double (*y)(double t),
		double tmin, double tmax);
	void operator()(double tt0,
		double (*y)(double t),
		double tmin, double tmax);

	BzzMinimizationMonoVeryRobust(double (*y)(double t),
		double tmin, double tmax);
	void operator()(double (*y)(double t),
		double tmin, double tmax);

	//	MONO_MULTI_MONO
	BzzMinimizationMonoVeryRobust(int iv, BzzVector& xx0, double yy0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	void operator()(int iv, BzzVector& xx0, double yy0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	BzzMinimizationMonoVeryRobust(int iv, BzzVector& xx0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	void operator()(int iv, BzzVector& xx0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);

	//	MONO_MULTI_MULTI
	BzzMinimizationMonoVeryRobust(BzzVector& ppi, BzzVector& xx0, double yy0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	void operator()(BzzVector& ppi, BzzVector& xx0, double yy0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	BzzMinimizationMonoVeryRobust(BzzVector& ppi, BzzVector& xx0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	void operator()(BzzVector& ppi, BzzVector& xx0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	BzzMinimizationMonoVeryRobust(BzzVector& ppi, BzzVector& xx0, double yy0,
		double (*y)(BzzVector& x),
		BzzVector& xxL, BzzVector& xxU);
	void operator()(BzzVector& ppi, BzzVector& xx0, double yy0,
		double (*y)(BzzVector& x),
		BzzVector& xxL, BzzVector& xxU);
	BzzMinimizationMonoVeryRobust(BzzVector& ppi, BzzVector& xx0,
		double (*y)(BzzVector& x),
		BzzVector& xxL, BzzVector& xxU);
	void operator()(BzzVector& ppi, BzzVector& xx0,
		double (*y)(BzzVector& x),
		BzzVector& xxL, BzzVector& xxU);

	//	MONO_SCANNING
	BzzMinimizationMonoVeryRobust(double tt0, double yy0,
		double tmin, double tmax);
	void operator()(double tt0, double yy0,
		double tmin, double tmax);

	//	============================================================================
	//	************************< Multidimensional Minimum >************************
	//	============================================================================
	//	*************************< Unconstrained Minimization >*********************
	//	* Min F = F(x)																					*
	//	* x >= L (lower bounds)																		*
	//	* x <= U (upper bounds)																		*

		//	UNCONSTRAINED_MULTI,
	BzzMinimizationMonoVeryRobust(BzzVector& x00, double FF0,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU);
	void operator()(BzzVector& x00, double FF0,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU);

	BzzMinimizationMonoVeryRobust(BzzVector& x00,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU);
	void operator()(BzzVector& x00,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzMinimizationMonoVeryRobust(void)
	{
		delete mono;
	}

	//	============================================================================
	//	******************< Non-modifying access functions >************************
	//	============================================================================
	double GetTSolution(void)
	{
		return tSolution;
	}
	double GetYSolution(void)
	{
		return ySolution;
	}
	double GetFSolution(void)
	{
		return fSolution;
	}
	void GetSolutions(BzzVector* tSolutions, BzzVector* ySolutions);
	void GetXSolution(BzzVector* xS)
	{
		*xS = xSolution;
	}
	void GetXSolutions(BzzMatrix* Xs);
	void GetSolutions(BzzMatrix* Xs, BzzVector* tSolutions,
		BzzVector* ySolutions);
	void GetXLateralToSolution(int lateral, BzzVector* xLeft, double* fLeft,
		BzzVector* xRight, double* fRight);

	double GetCurvature(void) { return lambdaMax; }
	int GetIterTotalFeasible(void)
	{
		return iterTotalFeasible;
	}
	int GetIterTotal(void)
	{
		return iterTotal;
	}
	int GetUnfeasibleGlobal(void)
	{
		return unfeasibleGlobal;
	}
	void ResetUnfeasibleGlobal(void)
	{
		unfeasibleGlobal = 0;
	}

	//	============================================================================
	//	****************************< Modifying access functions >******************
	//	============================================================================
	// return control < 0 Solution not yet founded with the required precision
	//			(itertions >= required max iteration)
	// return control = 0 Founded solution with a value < to setTollF required
	// return control > 0 Founded solutions with the required precision
	char operator() (void)
	{
		return (*this)(maxIter);
	}
	char operator()(int ni);

	double Scanning(double* yyy, int* unf, char* cc);
	void SetTolRel(double tr) { control = -1; tTollRel = fabs(tr); }
	void SetTolAbs(double ta) { control = -1; tTollAbs = fabs(ta); }
	void SetTolF(double yt = TOL_ABS_Y) { control = -1; yToll = yt; }
	void SetMaxIter(int mxit) { maxIter = mxit; }

	// standard value 10
	void SetExtraPointsNumberForMultipleMinimaSearch(int numP)
	{
		control = -1;
		if (numP <= 1)
			dtToll = 1.;
		else
			dtToll = 1. / double(numP);
	}

	void SetUnfeasibleFeasibleRatio(double rat)
	{
		control = -1;
		if (rat < .001)
			rat = .001;
		else if (rat > .5)
			rat = .5;
		unfeasibleFeasibleRatio = rat;
		oneMunfeasibleFeasibleRatio = 1. - unfeasibleFeasibleRatio;
	}

	void Unimodal(void) { unimodal = 1; }

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
	void SetPlot(void) { plotTask = 1; }

	//	============================================================================
	//	*******************************< Non-modifying functions >******************
	//	============================================================================
	virtual void ObjectBzzPrint(void);
	void BzzPlot(const char* title = "Monodimensional Minimum",
		const char* abscissa = "t", const char* ordinate = "y");
	void MultiPlot(const char* title = "Bidimensional Minimum",
		const char* abscissa = "x1", const char* ordinate = "x2");
};

#endif // MINIMIZATION_ROBUST_MONO_HPP