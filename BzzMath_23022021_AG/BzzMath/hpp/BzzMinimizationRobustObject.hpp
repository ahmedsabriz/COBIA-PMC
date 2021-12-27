// BZZMATH: Release 7.0

//	=======< BzzMinimizationRobustObject.hpp >========================
//	* The class BzzMinimizationRobustObject is used for unconstrained*
// * optimization																					*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-2008	Date Written.

//	============================================================================
//	******	constructors:																		*
//	* BzzMinimizationRobustObject m;m(x0,FunMin);													*
//	* BzzMinimizationRobustObject(x0,FunMin);														*
//	****************************************************************************
//	****** Functions for minimization:														*
//	* char control = m();																		*
//	* char control = m(niter);																	*
//	****************************************************************************
//	* Other access functions:																	*
//	* int numIterations = m.TotalIterations();											*
//	* int numIterations = m.IterationCounter();											*
//	* double F = m.GetBzzMinimumF();															*
//	* double F = m.GetSolution(&x);															*
//	* m.BzzPrint("Comments");																	*
//	* m.BzzMessage("Comments");																*
//	****************************************************************************
//	* Functions for assigning tollerances:													*
//	****************************************************************************
//	* Function for save partial results:													*
//	* Save("SAVE.TXT");																			*
//	****************************************************************************
//	* Restart: m(x0);																				*
//	* Restart: m.Restart(x0);																	*
//	****************************************************************************

#ifndef BZZ_SIMPLEX_PLUS_OBJECT_DOUBLE_HPP
#define BZZ_SIMPLEX_PLUS_OBJECT_DOUBLE_HPP

//	============================================================================
//	===========< class BzzMyMinimizationRobustObject >===============
//	============================================================================
//class BzzMyMinimizationRobustObject : public BzzBaseClass
//	{
//public:
//	virtual double GetFunctionValue(BzzVector &x){return 0.;}
//	virtual void ObjectBzzPrint(void){};
//	};

//	============================================================================
//	=================< class BzzMinimizationRobustObject >==================
//	============================================================================
class BzzMinimizationRobustObject : public BzzBaseClass
{
	friend class BzzMinimumMultiRobust;
private:
	enum MethodStatusBzzSimplexPlus
	{
		EXIT,
		START,
		SIMPLEX_ACTION
	}methodStatus;

	enum FunctionTypeBzzSimplex
	{
		TWO_IV1_IV2,
		TWO_PI1_PI2,
		UNCONSTRAINED_MULTI,
		CONSTRAINED_MULTI
	}functionSimplexType;

	static const char* const BZZ_ERROR;

	FILE* bzzFileSave;
	char saveResults;

	int	numVariables,
		numVertices,
		niter, iter, iterTotal, iterTotalFeasible, iterRidge,
		numIterations,
		iterSimplexCollapsed,
		printTasks,
		printSubTasks;

	char	control,
		controlCollapsed,
		stop;

	BzzVector	x0,	// start
		xi,	// current
		x,		// generic
		xMin, // x in minimum
		xCollapsed,
		h,
		* v,	 // vertices
		f,		// function value in v
		d, // control callapse
		vB,	 // base
		vR,	 // reflection
		vE,	 // expansion
		vC,	 // contraction
		vA,    //
		vP,	 // projected vertex
		dvPvB,
		dxa,
		aux,
		xL, xU; // boundary on variables

	BzzVector z0, zi, piA, piB;
	int iVar1, iVar2;

	double	fi,
		f0,
		fMin,
		//				fRidge,
		fCollapsed,
		dfMinMax,
		normdbx,
		normB,
		fx,
		xTollRel,
		xTollAbs,
		fTollRel,
		fTollAbs,
		fSaved;

	double	fB,	 // base
		fR,	 // reflection
		fE,	 // expansion
		fA,
		fP,
		fC;	 // contraction
	double yToll;

	BzzVectorInt	sortedV;

	//	double (*ptrFun)(BzzVector &x);
	BzzMyMinimizationRobustObject* ptrObject;

	double MinFunction(BzzVector& x);

	void CommonInitializations(BzzVector& x00);
	void Start(void);
	void Baricenter(void);
	void SimplexAction(void);
	char MinimumSolve(void);

public:
	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================
		// default constructor
	BzzMinimizationRobustObject(void);

	// copy constructor
	BzzMinimizationRobustObject
	(const BzzMinimizationRobustObject& rval);

	// SimplexPlus
	BzzMinimizationRobustObject
	(BzzVector& x00, BzzMyMinimizationRobustObject* ptro);
	void operator ()
		(BzzVector& x00, BzzMyMinimizationRobustObject* ptro);

	BzzMinimizationRobustObject
	(BzzVector& x00, double ff, BzzMyMinimizationRobustObject* ptro);
	void operator ()
		(BzzVector& x00, double ff, BzzMyMinimizationRobustObject* ptro);

	BzzMinimizationRobustObject
	(BzzVector& x00, BzzMyMinimizationRobustObject* ptro,
		BzzVector* xxL, BzzVector* xxU);
	void operator ()
		(BzzVector& x00, BzzMyMinimizationRobustObject* ptro,
			BzzVector* xxL, BzzVector* xxU);
	BzzMinimizationRobustObject
	(BzzVector& x00, BzzMyMinimizationRobustObject* ptro,
		BzzVector& xxL, BzzVector& xxU);
	void operator ()
		(BzzVector& x00, BzzMyMinimizationRobustObject* ptro,
			BzzVector& xxL, BzzVector& xxU);

	BzzMinimizationRobustObject
	(BzzVector& x00, double ff, BzzMyMinimizationRobustObject* ptro,
		BzzVector* xxL, BzzVector* xxU);
	void operator ()
		(BzzVector& x00, double ff, BzzMyMinimizationRobustObject* ptro,
			BzzVector* xxL, BzzVector* xxU);

	BzzMinimizationRobustObject
	(BzzVector& x00, double ff, BzzMyMinimizationRobustObject* ptro,
		BzzVector& xxL, BzzVector& xxU);
	void operator ()
		(BzzVector& x00, double ff, BzzMyMinimizationRobustObject* ptro,
			BzzVector& xxL, BzzVector& xxU);

	//TWO_IV1_IV2
	// (xS,fS,yy,iv1,tL1,tU1,iv2,tL2,tU2)
	BzzMinimizationRobustObject
	(BzzVector& xxS, double ffS, BzzMyMinimizationRobustObject* ptro,
		int iv1, double ttL1, double ttU1,
		int iv2, double ttL2, double ttU2);
	void operator ()
		(BzzVector& xxS, double ffS, BzzMyMinimizationRobustObject* ptro,
			int iv1, double ttL1, double ttU1,
			int iv2, double ttL2, double ttU2);

	//TWO_PI1,PI2
	// (xS,fS,yy,pi1,pi2,)
	BzzMinimizationRobustObject
	(BzzVector& xxS, double ffS, BzzMyMinimizationRobustObject* ptro,
		BzzVector& pi1, BzzVector& pi2,
		BzzVector* xxL, BzzVector* xxU);
	void operator ()
		(BzzVector& xxS, double ffS, BzzMyMinimizationRobustObject* ptro,
			BzzVector& pi1, BzzVector& pi2,
			BzzVector* xxL, BzzVector* xxU);

	void operator ()
		(BzzVector& x00, double ff);

	void operator ()(BzzVector& x00);
	void Restart(BzzVector& x0);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzMinimizationRobustObject(void);

	//	============================================================================
	//	********************< Non-modifying access functions >**********************
	//	============================================================================
	double GetSolution(BzzVector* xx)
	{
		(*xx) = xMin;
		return fMin;
	}

	double GetBzzMinimumF(void) { return fMin; }
	double GetMinimumF(void) { return fMin; }
	int TotalIterations(void) { return iterTotal; }
	int IterationCounter(void) { return iterTotal; }
	int GetIterTotalFeasible(void)
	{
		return iterTotalFeasible;
	}
	int GetIterTotal(void)
	{
		return iterTotal;
	}
	virtual void ObjectBzzPrint(void);
	void Save(char* sv);
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
	void SetTolF(double yt = -BZZ_BIG) { control = -1; yToll = yt; }
	void SetTolAbsF(double tollAbsF)
	{
		control = -1;
		fTollAbs = tollAbsF;
	}
	void SetTolRelF(double tollRelF)
	{
		control = -1;
		if (tollRelF > 10. * BZZ_MACH_EPS_DOUBLE)
			fTollRel = tollRelF;
		else
			fTollRel = 10. * BZZ_MACH_EPS_DOUBLE;
	}
	void SetTolAbsX(double tollAbsX)
	{
		control = -1;
		xTollAbs = tollAbsX;
	}
	void SetTolRelX(double tollRelX)
	{
		control = -1;
		if (tollRelX > 10. * BZZ_MACH_EPS_DOUBLE)
			xTollRel = tollRelX;
		else
			xTollRel = 10. * BZZ_MACH_EPS_DOUBLE;
	}

	//	============================================================================
	//	**********************< Modifying functions >*******************************
	//	============================================================================
	char operator ()(void)
	{
		return (*this)(0);
	}
	char operator ()(int ni);

	//	if(control < 0)
	//		::BzzPrint("\nSearch failed before reaching"
	//		" the solution\n");
	//	if(control == -3)
	//		::BzzPrint("\nUnfeasible points\n");
	//	if(control == 0)
	//		::BzzPrint("\nFound problems\n");
	//	if(control == 1)
	//		::BzzPrint("\nMaximum number of iterations\n");
	//	if(control == 2)
	//		::BzzPrint("\nSatisfied test on dx\n");
	//	if(control == 3)
	//		::BzzPrint("\nfSolution <= fToll\n");
	//	if(control == 4)
	//		::BzzPrint("\nSatisfied test on the function values for simplex vertices");
};

#endif // BZZ_SIMPLEX_PLUS_OBJECT_DOUBLE_HPP