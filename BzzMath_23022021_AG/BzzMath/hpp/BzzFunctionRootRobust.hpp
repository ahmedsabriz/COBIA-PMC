// BZZMATH: Release 7.0

//	==================< BzzFunctionRootRobust.hpp >=======================
//	* Class BzzFunctionRootRobust for minimization in double precision	*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	04-2004	Date Written.

#ifndef FUNCTION_ROOT_DOUBLE_ROBUST_HPP
#define FUNCTION_ROOT_DOUBLE_ROBUST_HPP

//	============================================================================
//	=========================< class BzzFunctionRootRobust >=====================
//	============================================================================

class BzzFunctionRootRobust : public BzzBaseClass
{
private:
	/*
		enum UnconstrainedMethodBzzMinimization
			{
			START,
			MOVE,
			RE_START,
			TWO_SIMPLEX,
			MULTI_SIMPLEX,
			ROB,
			ROS,
			ROD,
			EXIT
			}unconstrainedMethodBzzMinimization;
	*/
	static const char* const BZZ_ERROR;
	static const int MAX_ITER;
	static const double TOL_REL;
	static const double TOL_ABS;
	static const double TOL_REL_Y;
	static const double TOL_ABS_Y;
	static const double DT_RATIO;
	static const double DT_MIN;
	static const double U_F_RATIO;

	char	control, unimodal;

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
	BzzVector x0, pi, pi1, x, x1, x2, x3;
	BzzVector piA, piB; // for TWO_TWO_MULTI_MULTI

	double (*yMono)(double t);
	double (*yMulti)(BzzVector& x);
	void (*fMulti)(BzzVector& x, BzzVector& f);
	double (*FMulti)(BzzVector& x);
	double Y(double t);
	double Phiw(BzzVector& x);

	int	iVar,
		iter,
		iterFeasible,
		iterTotal,
		iterTotalInConstructor,
		iterFeasibleInConstructor,
		iterTotalFeasible,
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
	/*
		BzzVector	xL, // lower bound
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
	*/
	BzzVector YS;
	BzzMatrix XS;
	BzzVector tS1;

	BzzVectorInt sorted;
	BzzFunctionRootRobust* mono;

public:
	//	============================================================================
	//	********************************< constructors >****************************
	//	============================================================================
		// default constructor
	BzzFunctionRootRobust(void);

	// copy constructor
	BzzFunctionRootRobust(const BzzFunctionRootRobust& rval);

	//	============================================================================
	//	************************< Monodimensional Minimum >*************************
	//	============================================================================
		//	MONO_MONO
	BzzFunctionRootRobust(double tt0, double yy0,
		double (*y)(double t),
		double tmin, double tmax);
	void operator()(double tt0, double yy0,
		double (*y)(double t),
		double tmin, double tmax);

	BzzFunctionRootRobust(double tt0,
		double (*y)(double t),
		double tmin, double tmax);
	void operator()(double tt0,
		double (*y)(double t),
		double tmin, double tmax);

	BzzFunctionRootRobust(double (*y)(double t),
		double tmin, double tmax);
	void operator()(double (*y)(double t),
		double tmin, double tmax);

	//	MONO_MULTI_MONO
	BzzFunctionRootRobust(int iv, BzzVector& xx0, double yy0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	void operator()(int iv, BzzVector& xx0, double yy0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	BzzFunctionRootRobust(int iv, BzzVector& xx0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	void operator()(int iv, BzzVector& xx0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);

	//	MONO_MULTI_MONO_NLE
	BzzFunctionRootRobust(int iv, BzzVector& xx0, double yy0,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w,
		double tmin, double tmax);
	void operator()(int iv, BzzVector& xx0, double yy0,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w,
		double tmin, double tmax);
	BzzFunctionRootRobust(int iv, BzzVector& xx0,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w,
		double tmin, double tmax);
	void operator()(int iv, BzzVector& xx0,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w,
		double tmin, double tmax);

	//	void MonoMultiMono(int iv,BzzVector &xx0,double yy0,
	//		double (*y)(BzzVector &x),
	//		double tmin,double tmax);

		//	MONO_MULTI_MULTI
	BzzFunctionRootRobust(BzzVector& ppi, BzzVector& xx0, double yy0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	void operator()(BzzVector& ppi, BzzVector& xx0, double yy0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	BzzFunctionRootRobust(BzzVector& ppi, BzzVector& xx0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	void operator()(BzzVector& ppi, BzzVector& xx0,
		double (*y)(BzzVector& x),
		double tmin, double tmax);

	//	MONO_MULTI_MULTI_NLE
	BzzFunctionRootRobust(BzzVector& ppi, BzzVector& xx0, double yy0,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w,
		double tmin, double tmax);
	void operator()(BzzVector& ppi, BzzVector& xx0, double yy0,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w,
		double tmin, double tmax);
	BzzFunctionRootRobust(BzzVector& ppi, BzzVector& xx0,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w,
		double tmin, double tmax);
	void operator()(BzzVector& ppi, BzzVector& xx0,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w,
		double tmin, double tmax);

	// TWO_MONO_MONO
	BzzFunctionRootRobust(BzzVector& xx1, double yy1,
		BzzVector& xx2, double yy2,
		double (*y)(BzzVector& x),
		double tmin, double tmax);
	void operator()(BzzVector& xx1, double yy1,
		BzzVector& xx2, double yy2,
		double (*y)(BzzVector& x),
		double tmin, double tmax);

	// TWO_MONO_MONO_NLE
	BzzFunctionRootRobust(BzzVector& xx1, double yy1,
		BzzVector& xx2, double yy2,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w,
		double tmin, double tmax);
	void operator()(BzzVector& xx1, double yy1,
		BzzVector& xx2, double yy2,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w,
		double tmin, double tmax);

	//	THREE_MONO_MULTI_MULTI
	BzzFunctionRootRobust(BzzVector& xx1, double yy1,
		BzzVector& xx2, double yy2,
		BzzVector& xx3, double yy3,
		double (*yy)(BzzVector& x));
	void operator()(BzzVector& xx1, double yy1,
		BzzVector& xx2, double yy2,
		BzzVector& xx3, double yy3,
		double (*yy)(BzzVector& x));

	//	THREE_MONO_MULTI_MULTI_NLE
	BzzFunctionRootRobust(BzzVector& xx1, double yy1,
		BzzVector& xx2, double yy2,
		BzzVector& xx3, double yy3,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w);
	void operator()(BzzVector& xx1, double yy1,
		BzzVector& xx2, double yy2,
		BzzVector& xx3, double yy3,
		void (*y)(BzzVector& x, BzzVector& f), BzzVector& w);

	// TWO_TWO_MONO_MULTI
	BzzFunctionRootRobust(int iv1, int iv2, BzzVector& x00,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU);
	void operator()(int iv1, int iv2, BzzVector& x00,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU);

	// TWO_TWO_MULTI_MULTI
	BzzFunctionRootRobust(BzzVector& ppi1, BzzVector& ppi2,
		BzzVector& x00,
		double (*yy)(BzzVector& x),
		BzzVector* xxU);
	void operator()(BzzVector& ppi1, BzzVector& ppi2,
		BzzVector& x00,
		double (*yy)(BzzVector& x),
		BzzVector* xxU);

	//	MONO_SCANNING
	BzzFunctionRootRobust(double tt0, double yy0,
		double tmin, double tmax);
	void operator()(double tt0, double yy0,
		double tmin, double tmax);

	//	============================================================================
	//	************************< Multidimensional Minimum >************************
	//	============================================================================

	//	*****************************< Linear Programming >*************************
	//	* Min F = sTx																					*
	//	* Ex = e	(mE linear equations)															*
	//	* Dx >= d (mD linear disequations)														*
	//	* x >= L (lower bounds)																		*
	//	* x <= U (upper bounds)																		*

		//	LINEAR_MULTI,
	BzzFunctionRootRobust(BzzVector& x00, BzzVector* ss,
		BzzMatrixSparse* EE, BzzVector* ee,
		BzzMatrixSparse* DD, BzzVector* dd,
		BzzVectorInt* iiE, BzzVector* xxE,
		BzzVectorInt* iiL, BzzVector* xxL,
		BzzVectorInt* iiU, BzzVector* xxU);

	// constructor FILE .BZZ
	BzzFunctionRootRobust(char* fileLP);

	//	*************************< Unconstrained Minimization >*********************
	//	* Min F = F(x)																					*
	//	* x >= L (lower bounds)																		*
	//	* x <= U (upper bounds)																		*

		//	UNCONSTRAINED_MULTI,
	BzzFunctionRootRobust(BzzVector& x00, double FF0,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU);
	void operator()(BzzVector& x00, double FF0,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU);

	//	***************************< Constrained Minimization >*********************
	//	* Min F = F(x)																					*
	//	* h(x) = 0	(nh non linear equations)													*
	//	* g(x) >= 0	(ng non linear disequations)												*
	//	* Ex = e	(mE linear equations)															*
	//	* Dx >= d (mD linear disequations)														*
	//	* x >= L (lower bounds)																		*
	//	* x <= U (upper bounds)																		*

		//	CONSTRAINED_MULTI
	BzzFunctionRootRobust(BzzVector& x00, double (*y)(BzzVector& x),
		int nh, void (*h)(BzzVector& x, BzzVector& h),
		int ng, void (*g)(BzzVector& x, BzzVector& g),
		BzzMatrixSparse* EE, BzzVector* ee,
		BzzMatrixSparse* DD, BzzVector* dd,
		BzzVectorInt* iiE, BzzVector* xxE,
		BzzVectorInt* iiL, BzzVector* xxL,
		BzzVectorInt* iiU, BzzVector* xxU);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzFunctionRootRobust(void)
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
	void GetSolutions(BzzVector* tSolutions, BzzVector* ySolutions);

	//	double GetCurvature(void){return lambdaMax;}
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
	// return control = 0 Founded solution with a value < to setTollY required
	// return control > 0 Founded solutions with the required precision
	char operator() (void)
	{
		return (*this)(maxIter);
	}
	char operator()(int ni);

	//	double Scanning(double *yyy,int *unf,char *cc);
	void SetTolRel(double tr);
	void SetTolAbs(double ta);
	void SetTolY(double yt);
	void SetMaxIter(int mxit);

	// standard value 100
	void SetExtraPointsNumberForMultipleRootSearch(int numP)
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

	void OnlyOneRoot(void) { unimodal = 1; }

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
	//void BzzPlot(char *title = "Monodimensional Minimum",
	//	char *abscissa = "t",char *ordinate = "y");
	void BzzPlot(const char* title = "Monodimensional Minimum",
		const char* abscissa = "t", const char* ordinate = "y");
	//void MultiPlot(char *title = "Bidimensional Minimum",
	//	char *abscissa = "x1",char *ordinate = "x2");
	void MultiPlot(const char* title = "Bidimensional Minimum",
		const char* abscissa = "x1", const char* ordinate = "x2");
};

#endif // FUNCTION_ROOT_DOUBLE_ROBUST_HPP