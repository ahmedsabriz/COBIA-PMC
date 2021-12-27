// BZZMATH: Release 7.0
//	=====================< BzzMinimizationLargeSparseNewton.hpp >========================
// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	12-2012	Date Written.

#ifndef BZZ_MINIMIZATION_LARGE_SPARSE_NEWTON
#define BZZ_MINIMIZATION_LARGE_SPARSE_NEWTON

//	============================================================================
//	================< class BzzMinimizationLargeSparseNewton >==================
//	============================================================================
class BzzMinimizationLargeSparseNewton : public BzzBaseClass
{
private:
	enum LargeMethodType
	{
		CG_METHOD,
		PCG_METHOD
	}largeMethodType;

	enum LargeMethodSelected
	{
		CG_CG_METHOD,
		CG_PCG_METHOD,
		PCG_PCG_METHOD
	}selectedMethod;

	static const char* const BZZ_ERROR;
	int	numVariables,
		numThreads,
		printTasks,
		printSubTasks,
		functionCalls, totalFunctionCalls,
		maxFunctionCalls, maxTotalFunctionCalls,
		internalIterations, maxInternalIterations,
		externalIterations, externalTotalIterations,
		maxExternalIterations, maxExternalTotalIterations,
		hessianEvaluation,
		convergenceSpeed,
		conditionerBand;

	char	control,
		stop;

	double (*ptrFun)(BzzVector& x);

	BzzMatrixSparseSymmetricLocked G;
	BzzFactorizedSymmetricBand M;
	BzzMatrixSymmetricBand B;

	BzzVector	xi, xmi, x, xOpt,
		x0,
		gi;

	BzzVector xL, xU;

	double	fi, fmi, f, fOpt, f1, f2,
		f0,
		xTolAbs,
		xTolRel;

	void CGCG(void);
	void CGPCG(void);
	void PCGPCG(void);
	BzzVector dCG, zCG, qCG, pCG, vCG, wCG, gCG;

public:
	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================
		// default constructor
	BzzMinimizationLargeSparseNewton(void);

	BzzMinimizationLargeSparseNewton
	(BzzVector* x00, double (*ptr)(BzzVector& x));
	void operator()(BzzVector* x00, double (*ptr)(BzzVector& x));

	BzzMinimizationLargeSparseNewton
	(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG, double (*ptr)(BzzVector& x));
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		double (*ptr)(BzzVector& x));

	//	============================================================================
	//	********************< Non-modifying access functions >**********************
	//	============================================================================
	void SetMaxExternalIterations(int me)
	{
		if (me <= 0)
			BzzError("SetMaxExternalIterations <= 0");
		maxExternalIterations = me;
	}

	void SetMaxExternalTotalIterations(int me)
	{
		if (me <= 0)
			BzzError("SetMaxExternalTotalIterations <= 0");
		maxExternalTotalIterations = me;
	}

	virtual void ObjectBzzPrint(void);

	//	============================================================================
	//	**********************< Modifying functions >*******************************
	//	============================================================================
	char operator ()(void);
};

#endif // BZZ_MINIMIZATION_LARGE_SPARSE_NEWTON