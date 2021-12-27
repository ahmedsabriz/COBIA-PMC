
// BZZMATH: Release 7.0

//	=====================< BzzQuadraticProgramming.hpp >========================
//	* BzzMatrixSparseSymmetricLocked: Class for sparse and symmetrical matrices*
// ============================================================================
// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	12-2012	Date Written.															 

#ifndef BZZ_QUADRATIC_PROGRAMMING
#define BZZ_QUADRATIC_PROGRAMMING

class BzzQuadraticProgramming : public BzzBaseClass
	{
private:
	enum BzzQuadraticProgrammingType
		{
		DENSE_BOUND,
		SPARSE_BOUND,
		DENSE_BOUND_EQUALITY,
		SPARSE_BOUND_EQUALITY,
		DENSE_BOUND_EQUALITY_INEQUALITY,
		SPARSE_BOUND_EQUALITY_INEQUALITY,
		}quadraticProgrammingType;

	static const char *const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	int	numVariables,
			numEquality,
			numInequality,
			numReducedVariables,
			numActiveBound,
			numActiveInequality,
			feasible, // 0 non ancora feasible 1 in xFeasible, fFeasible
			control, // 1 solution 0 nt yet -1 unfeasible
			printTasks,
			printSubTasks;

	//Rows,numColumns,numElements;
	int whoAmI;

// small problems	
	BzzMatrixSymmetric G,R;
	BzzFactorizedSymmetric S;
	BzzMatrix AA;
	BzzVector bb,z,y;
	BzzFactorizedLQ A;
	BzzMatrix N;
	BzzFactorizedGauss GG;

// large problems
	BzzMatrixSparseSymmetricLocked GS,RS;
	BzzVector dS,dSr;
	BzzVectorInt rS,rSr;
	BzzVectorInt cS,cSr;
	BzzVector lS,lSr;

	BzzMatrix E,D;
	BzzMatrixSparseLocked ES,DS;

	BzzVector	s,r,
					x,
					xL,
					xU,
					xFeasible,
					dxUL,epsUL,
					x0,
					xSolution,
					xA,xB,xC,
					g,
					e,aa,
					d,dL,dU,dd,
					lambda;

	BzzVectorInt	lowerBound,
						upperBound,
						boundedVariables,
						lambdaVariables,
						activeBound, // -1 lower 0 passive 1 upper 
						passiveBound, // 0 can be active 1 must be passive
						iActiveBound, // dimensioned numActiveBound; iActiveBound[i] = bound active
						activeInequality, // -1 lower 0 passive 1 upper 
						iActiveInequality; // dimensioned numActiveInequality; iActiveInequality[i] = Inequality active
	
	double f0,fSolution,fFeasible,tA,tB,tC,tP,fA,fB,fC,alfa;

	void ControlBound(BzzVector *xx);
	void GetGradient(BzzVector &x);
	void GetDirection(void);
	void GetParabolicSolution(void);
	double GetFunction(BzzVector &x);
	int FindConstraint(int *con,double *t); 
	void GetReducedKKT(void);

public:
//	============================================================================
//	*****************************< constructors >*******************************
//	============================================================================
	// default constructor
	BzzQuadraticProgramming(void);

// quadraticProgrammingType = DENSE_BOUND
	void operator()(BzzMatrixSymmetric *GG,BzzVector *ss,BzzVector *ll,BzzVector *uu);
	void operator()(BzzVector *x00,BzzMatrixSymmetric *GG,BzzVector *ss,BzzVector *ll,BzzVector *uu);
// quadraticProgrammingType = SPARSE_BOUND
	void operator()(BzzMatrixSparseSymmetricLocked *GGS,BzzVector *ss,BzzVector *ll,BzzVector *uu);
	void operator()(BzzVector *x00,BzzMatrixSparseSymmetricLocked *GGS,BzzVector *ss,BzzVector *ll,BzzVector *uu);

// quadraticProgrammingType = DENSE_BOUND_EQUALITY
	void operator()(BzzMatrixSymmetric *GG,BzzVector *ss,BzzVector *ll,BzzVector *uu,
		BzzMatrix *EE,BzzVector *ee);
	void operator()(BzzVector *x00,BzzMatrixSymmetric *GG,BzzVector *ss,BzzVector *ll,BzzVector *uu,
		BzzMatrix *EE,BzzVector *ee);
	void DenseBoundEquality(void);
// quadraticProgrammingType = SPARSE_BOUND_EQUALITY
	void operator()(BzzMatrixSparseSymmetricLocked *GGS,BzzVector *ss,BzzVector *ll,BzzVector *uu,
		BzzMatrixSparseLocked *EES,BzzVector *ee);
	void operator()(BzzVector *x00,BzzMatrixSparseSymmetricLocked *GGS,BzzVector *ss,BzzVector *ll,BzzVector *uu,
		BzzMatrixSparseLocked *EES,BzzVector *ee);

// quadraticProgrammingType = DENSE_BOUND_EQUALITY_INEQUALITY
	void operator()(BzzMatrixSymmetric *GG,BzzVector *ss,BzzVector *ll,BzzVector *uu,
		BzzMatrix *EE,BzzVector *ee,BzzMatrix *DD,BzzVector *ddL,BzzVector *ddU);
	void operator()(BzzVector *x00,BzzMatrixSymmetric *GG,BzzVector *ss,BzzVector *ll,BzzVector *uu,
		BzzMatrix *EE,BzzVector *ee,BzzMatrix *DD,BzzVector *ddL,BzzVector *ddU);
	void DenseBoundEqualityInequality(void);
// quadraticProgrammingType = SPARSE_BOUND_EQUALITY_INEQUALITY
	void operator()(BzzMatrixSparseSymmetricLocked *GGS,BzzVector *ss,BzzVector *ll,BzzVector *uu,
		BzzMatrixSparseLocked *EES,BzzVector *ee,
		BzzMatrixSparseLocked *DDS,BzzVector *ddL,BzzVector *ddU);
	void operator()(BzzVector *x00,BzzMatrixSparseSymmetricLocked *GGS,BzzVector *ss,BzzVector *ll,BzzVector *uu,
		BzzMatrixSparseLocked *EES,BzzVector *ee,
		BzzMatrixSparseLocked *DDS,BzzVector *ddL,BzzVector *ddU);
	
	void operator()(void);

//	============================================================================
//	=========================< Non-modifying functions >========================
//	============================================================================
	void SetTasksPrint(void){printTasks = 1;}
	void SetSubTasksPrint(int psb = 1)
		{
		if(psb > 0)
			printSubTasks = psb;
		else 
			printSubTasks = 1;
		}
//	*******************************< BzzPrint >*********************************
	virtual void ObjectBzzPrint(void);
};

#endif // BZZ_QUADRATIC_PROGRAMMING
