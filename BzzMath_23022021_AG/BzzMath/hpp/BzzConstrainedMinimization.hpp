// BZZMATH: Release 7.0

//	==================< BzzConstrainedMinimization.hpp >=======================
//	* Class BzzConstrainedMinimization	*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	04-2013	Date Written.

#ifndef CONTRAINED_MINIMIZATION_HPP
#define CONTRAINED_MINIMIZATION_HPP

//	============================================================================
//	===================< class BzzConstrainedMinimization >=================
//	============================================================================
//==================
//	Min F = F(x)
//	xL <= x <= xU
//	Ex = e
//	dL <= Dx <= dU
//	h(x) = 0
//	pL <= p(x) <= pU
//==================
class BzzConstrainedMinimization;
class BzzMyMinimizationConstrained : public BzzMyMinimizationVeryRobustObject
{
	friend class BzzConstrainedMinimization;
private:
	BzzConstrainedMinimization* ptrConstrained;
	double weightLinearEquality, weightLinearInequalityL, weightLinearInequalityU;
	double weightNonlinearEquality, weightNonlinearInequalityL, weightNonlinearInequalityU;
public:
	~BzzMyMinimizationConstrained(void) {};
	BzzMyMinimizationConstrained(void)
	{
		ptrConstrained = 0;
	}
	BzzMyMinimizationConstrained(BzzConstrainedMinimization* ptr)
	{
		ptrConstrained = ptr;
	}
	void operator()(BzzConstrainedMinimization* ptr)
	{
		ptrConstrained = ptr;
	}
	void SetWeights(double w1, double w2, double w3, double w4, double w5, double w6)
	{
		weightLinearEquality = w1;
		weightLinearInequalityL = w2;
		weightLinearInequalityU = w3;
		weightNonlinearEquality = w4;
		weightNonlinearInequalityL = w5;
		weightNonlinearInequalityU = w6;
	}

	virtual double GetFunctionValue(BzzVector& x);
	virtual void ObjectBzzPrint(void) {};
};

class BzzMyConstrainedMinimizationMonoObject : public BzzMyMinimizationMonoObject
{
	friend class BzzConstrainedMinimization;
	//	BzzVector x0,dx;
	BzzConstrainedMinimization* ptrConstrained;
	double weightLinearEquality, weightLinearInequalityL, weightLinearInequalityU;
	double weightNonlinearEquality, weightNonlinearInequalityL, weightNonlinearInequalityU;
	//	BzzVector xS,dx;
public:
	~BzzMyConstrainedMinimizationMonoObject(void) {};
	BzzMyConstrainedMinimizationMonoObject(void)
	{
		ptrConstrained = 0;
	}
	BzzMyConstrainedMinimizationMonoObject(BzzConstrainedMinimization* ptr)
	{
		ptrConstrained = ptr;
	}
	void operator()(BzzConstrainedMinimization* ptr)
	{
		ptrConstrained = ptr;
	}
	void SetWeights(double w1, double w2, double w3, double w4, double w5, double w6)
	{
		weightLinearEquality = 100. * w1;
		weightLinearInequalityL = 100. * w2;
		weightLinearInequalityU = 100. * w3;
		weightNonlinearEquality = 100. * w4;
		weightNonlinearInequalityL = 100. * w5;
		weightNonlinearInequalityU = 100. * w6;
	}
	virtual double GetFunctionValue(double t);
	virtual void ObjectBzzPrint(void) {};
};

class BzzConstrainedMinimization : public BzzBaseClass
{
	friend class BzzMyMinimizationConstrained;
	friend class BzzMyConstrainedMinimizationMonoObject;
private:
	BzzMinimizationMultiVeryRobustObject m;
	BzzMinimizationMonoObject mono;
	BzzMyMinimizationConstrained myptr;
	BzzMyConstrainedMinimizationMonoObject myMonoptr;
	enum BzzHessianType
	{
		NO_Hessian,
		Hessian_Dense,
		Hessian_Sparse
	}hessianType;
	enum BzzMatricesType
	{
		NO_Matrices,
		Matrices_Dense,
		Matrices_Sparse
	}matricesType;

	enum BzzConstrainedMinimizationType
	{
		LinearFunction_Bound,
		LinearFunction_Bound_Linear,
		SparseLinearFunction_Bound_Linear,
		LinearFunction_Bound_Linear_Equality,
		SparseLinearFunction_Bound_Linear_Equality,
		LinearFunction_Bound_Linear_Equality_Inequality,
		SparseLinearFunction_Bound_Linear_Equality_Inequality,

		QuadraticFunction_Bound,
		SparseQuadraticFunction_Bound,
		QuadraticFunction_Bound_Linear,
		SparseQuadraticFunction_Bound_Linear,
		QuadraticFunction_Bound_Linear_Equality,
		SparseQuadraticFunction_Bound_Linear_Equality,
		QuadraticFunction_Bound_Linear_Equality_Inequality,
		SparseQuadraticFunction_Bound_Linear_Equality_Inequality,

		GenericFunction_Bound,
		SoftGenericFunction_Bound, // Soft viene data la struttura dell'Hessiano
		GenericFunction_Bound_Linear,
		SoftGenericFunction_Bound_Linear,
		SparseGenericFunction_Bound_Linear,
		SoftSparseGenericFunction_Bound_Linear,
		GenericFunction_Bound_Linear_Equality,
		SoftGenericFunction_Bound_Linear_Equality,
		SparseGenericFunction_Bound_Linear_Equality,
		SoftSparseGenericFunction_Bound_Linear_Equality,
		GenericFunction_Bound_Linear_Equality_Inequality,
		SoftGenericFunction_Bound_Linear_Equality_Inequality,
		SparseGenericFunction_Bound_Linear_Equality_Inequality,
		SoftSparseGenericFunction_Bound_Linear_Equality_Inequality
	}constrainedMinimizationType;

	enum BzzConstrainedFunctionType
	{
		LinearFunction,
		QuadraticFunction,
		SoftGenericFunction,
		GenericFunction
	}constrainedFunctionType;

	static const char* const BZZ_ERROR;
	static const int MAX_ITER;
	//	static const double TOL_REL;
	//	static const double TOL_ABS;
	int	numVariables,
		numLinearEquality,
		numLinearInequality,
		numNonLinearEquality,
		numNonLinearInequality,
		kktSize,
		atticSize,
		control,
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

	double	f0,
		fSolution,
		fx;

	BzzVector	x0,
		xSolution,
		x,
		g, // gradient
		s,
		e, ee,
		dL,
		d,
		dU,
		xL,
		xU,
		h,
		p,
		pL,
		pU,
		dx,
		xKKT,
		bKKT,
		aux1,
		aux2;

	BzzMatrix	E,
		D,
		A,
		KKT;

	BzzMatrixSparseLocked	ES,
		DS,
		AS,
		KKTS,
		HS,
		PS, PPS;
	BzzVectorInt nlH, nlP;

	double cq;
	BzzMatrixSymmetric G;
	BzzMatrixSparseSymmetricLocked GS;

	double w1, w2, w3, w4, w5, w6;
	BzzVectorInt	dActive, // -1 0 1
		pActive, // -1 0 1
		vActive; // -1 0 1
	int	numActiveD,
		numActiveP, numNonActiveP,
		numActiveV;
	BzzVectorInt	iActiveD,
		iActiveP, iNonActiveP, iAP,
		iActiveV;

	BzzVector lambdaE, lambdaD, lambdaH, lambdaP, lambdaV;
	BzzNonLinearSystemUtilities nlsuH, nlsuP;
	BzzNonLinearSystemUtilities nlsuHA, nlsuPA;

	int iterKKT;
	int ControlActiveConstraints(void);
	int ControlActiveConstraintsForAttic(void);
	int BuildKKTSystemDense(void);
	int BuildKKTSystemSparse(void);
	int BuildAtticSystemDense(void);
	int BuildAtticSystemSparse(void);

	double (*ptrFun)(BzzVector& x);
	void (*ptrEquality)(BzzVector& x, BzzVector& f);
	void (*ptrInequality)(BzzVector& x, BzzVector& f);

	//	BzzMinimizationMonoVeryRobust *mono;
	//	BzzMinimizationRobust *multi;
	//	BzzVector *x0;
	//	double *y0;

	int LinearFunctionBound(void);
	int LinearFunctionBoundLinear(void);
	int SparseLinearFunctionBoundLinear(void);
	int LinearFunctionBoundLinearEquality(void);
	int SparseLinearFunctionBoundLinearEquality(void);
	int LinearFunctionBoundLinearEqualityInequality(void);
	int SparseLinearFunctionBoundLinearEqualityInequality(void);

	int QuadraticFunctionBound(void);
	int SparseQuadraticFunctionBound(void);
	int QuadraticFunctionBoundLinear(void);
	int SparseQuadraticFunctionBoundLinear(void);
	int QuadraticFunctionBoundLinearEquality(void);
	int SparseQuadraticFunctionBoundLinearEquality(void);
	int QuadraticFunctionBoundLinearEqualityInequality(void);
	int SparseQuadraticFunctionBoundLinearEqualityInequality(void);

	int GenericFunctionBound(void);
	int SoftGenericFunctionBound(void);
	int GenericFunctionBoundLinear(void);
	int SoftGenericFunctionBoundLinear(void);
	int SparseGenericFunctionBoundLinear(void);
	int SoftSparseGenericFunctionBoundLinear(void);
	int GenericFunctionBoundLinearEquality(void);
	int SoftGenericFunctionBoundLinearEquality(void);
	int SparseGenericFunctionBoundLinearEquality(void);
	int SoftSparseGenericFunctionBoundLinearEquality(void);
	int GenericFunctionBoundLinearEqualityInequality(void);
	int SoftGenericFunctionBoundLinearEqualityInequality(void);
	int SparseGenericFunctionBoundLinearEqualityInequality(void);
	int SoftSparseGenericFunctionBoundLinearEqualityInequality(void);

public:
	//////////////////////////////////////////////////////
	// TODO
	// void SetInviolableLowerBounds(BzzVector &ilb);
	// void SetInviolableUpperBounds(BzzVector &iub);
	//////////////////////////////////////////////////////

	//	============================================================================
	//	********************************< constructors >****************************
	//	============================================================================
		// default constructor
	BzzConstrainedMinimization(void);

	BzzConstrainedMinimization(const BzzMinimizationTwoVeryRobust& rval)	// copy constructor
	{
		BzzError("Copy Constructor");
	}

	// =====================< Linear function >========================
	// LinearFunction_Bound
	void operator()(BzzVector* x00, BzzVector* ss,
		BzzVector* xxL, BzzVector* xxU);

	//	=========================< LinearFunction_Bound_Linear >=========================
	void operator()(BzzVector* x00, BzzVector* ss,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU);

	//	===================< SparseLinearFunction_Bound_Linear >=========================
	void operator()(BzzVector* x00, BzzVector* ss,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU);

	//	=====================< LinearFunction_Bound_Linear_Equality >====================
	void operator()(BzzVector* x00, BzzVector* ss,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h));

	//	================< SparseLinearFunction_Bound_Linear_Equality >===================
	void operator()(BzzVector* x00, BzzVector* ss,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzMatrixSparseLocked* HH, BzzVectorInt* nnlH);

	//	===========< LinearFunction_Bound_Linear_Equality_Inequality >===================
	void operator()(BzzVector* x00, BzzVector* ss,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzVector* ppL, BzzVector* ppU, void (*pp)(BzzVector& x, BzzVector& p));
	//	========< SparseLinearFunction_Bound_Linear_Equality_Inequality >================
	void operator()(BzzVector* x00, BzzVector* ss,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzMatrixSparseLocked* HH, BzzVectorInt* nnlH,
		BzzVector* ppL, BzzVector* ppU, void (*pp)(BzzVector& x, BzzVector& p),
		BzzMatrixSparseLocked* PP, BzzVectorInt* nnlP);

	// =====================< Quadratic Function >========================
//	=========================< QuadraticFunction_Bound >=============================
	void operator()(BzzVector* x00, BzzMatrixSymmetric* GG, BzzVector* ss, double cc,
		BzzVector* xxL, BzzVector* xxU);
	//	=======================< SparseQuadraticFunction_Bound >=========================
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		BzzVector* ss, double cc, BzzVector* xxL, BzzVector* xxU);

	//	=========================< QuadraticFunction_Bound_Linear >======================
	void operator()(BzzVector* x00, BzzMatrixSymmetric* GG, BzzVector* ss, double cc,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU);

	//	================< SparseQuadraticFunction_Bound_Linear >=========================
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		BzzVector* ss, double cc,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU);

	//	===============< QuadraticFunction_Bound_Linear_Equality >=======================
	void operator()(BzzVector* x00, BzzMatrixSymmetric* GG, BzzVector* ss, double cc,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h));

	//	===========< SparseQuadraticFunction_Bound_Linear_Equality >=====================
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		BzzVector* ss, double cc,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzMatrixSparseLocked* HH, BzzVectorInt* nnlH);

	//	=============< QuadraticFunction_Bound_Linear_Equality_Inequality >==============
	void operator()(BzzVector* x00, BzzMatrixSymmetric* GG, BzzVector* ss, double cc,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzVector* ppL, BzzVector* ppU, void (*pp)(BzzVector& x, BzzVector& p));

	//	==============< SparseQuadraticFunction_Bound_Linear_Equality_Inequality >=======
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		BzzVector* ss, double cc,
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzMatrixSparseLocked* HH, BzzVectorInt* nnlH,
		BzzVector* ppL, BzzVector* ppU, void (*pp)(BzzVector& x, BzzVector& p),
		BzzMatrixSparseLocked* PP, BzzVectorInt* nnlP);

	// =====================< Generic Function >========================
	//
//	=========================< GenericFunction_Bound >===============================
	void operator()(BzzVector* x00, double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU);
	//	=========================< SoftGenericFunction_Bound >===========================
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU);

	//	=========================< GenericFunction_Bound_Linear >========================
	void operator()(BzzVector* x00, double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU);

	//	=========================< SoftFunction_Bound_Linear >========================
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU);

	//	==================< SparseGenericFunction_Bound_Linear >=========================
	void operator()(BzzVector* x00,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU);

	//	==================< SoftSparseGenericFunction_Bound_Linear >=========================
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU);

	//	==================< GenericFunction_Bound_Linear_Equality >======================
	void operator()(BzzVector* x00, double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h));

	//	==================< SoftGenericFunction_Bound_Linear_Equality >===================
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h));

	//	=============< SparseGenericFunction_Bound_Linear_Equality >=====================
	void operator()(BzzVector* x00,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzMatrixSparseLocked* HH, BzzVectorInt* nnlH);

	//	=============< SoftSparseGenericFunction_Bound_Linear_Equality >=====================
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzMatrixSparseLocked* HH, BzzVectorInt* nnlH);

	//	================< GenericFunction_Bound_Linear_Equality_Inequality >=============
	void operator()(BzzVector* x00, double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzVector* ppL, BzzVector* ppU, void (*pp)(BzzVector& x, BzzVector& p));

	//	================< SoftGenericFunction_Bound_Linear_Equality_Inequality >=============
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrix* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrix* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzVector* ppL, BzzVector* ppU, void (*pp)(BzzVector& x, BzzVector& p));

	//	==============< SparseGenericFunction_Bound_Linear_Equality_Inequality >=========
	void operator()(BzzVector* x00,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzMatrixSparseLocked* HH, BzzVectorInt* nnlH,
		BzzVector* ppL, BzzVector* ppU, void (*pp)(BzzVector& x, BzzVector& p),
		BzzMatrixSparseLocked* PP, BzzVectorInt* nnlP);

	//	==============< SoftSparseGenericFunction_Bound_Linear_Equality_Inequality >=========
	void operator()(BzzVector* x00, BzzMatrixSparseSymmetricLocked* GG,
		double (*yy)(BzzVector& x),
		BzzVector* xxL, BzzVector* xxU,
		BzzMatrixSparseLocked* EE, BzzVector* ee,
		BzzVector* ddL, BzzMatrixSparseLocked* DD, BzzVector* ddU,
		int nnH, void (*hh)(BzzVector& x, BzzVector& h),
		BzzMatrixSparseLocked* HH, BzzVectorInt* nnlH,
		BzzVector* ppL, BzzVector* ppU, void (*pp)(BzzVector& x, BzzVector& p),
		BzzMatrixSparseLocked* PP, BzzVectorInt* nnlP);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzConstrainedMinimization(void);

	//	============================================================================
	//	******************< Non-modifying access functions >************************
	//	============================================================================
	double GetXSolution(BzzVector* xS)
	{
		*xS = xSolution;return fSolution;
	}
	void GetSolutions(BzzMatrix* Xs, BzzVector* fS);

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
	int operator() (void);

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

#endif // CONTRAINED_MINIMIZATION_HPP