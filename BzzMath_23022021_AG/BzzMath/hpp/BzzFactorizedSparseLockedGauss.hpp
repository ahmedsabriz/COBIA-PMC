// BZZMATH: Release 7.0

// ======================================================================================
//	========================< BzzFactorizedSparseLockedGauss.hpp >========================
// Sparse Linear Systems
// ======================================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	12-2013	Date Written

#ifndef BZZ_FACTORIZED_SPARSE_LOCKED_GAUSS_HPP
#define BZZ_FACTORIZED_SPARSE_LOCKED_GAUSS_HPP
class BzzFactorizedSparseLockedGauss : public BzzBaseClass
{
	friend void Solve
	(BzzFactorizedSparseLockedGauss* A, BzzVector* bx);
	friend void TransposeSolve
	(BzzFactorizedSparseLockedGauss* A, BzzVector* bx);
private:
	enum FactorizationStatus
	{
		UNFACTORIZED,
		FACTORIZED
	}factorizationStatus;

	enum BalancingStatus
	{
		NO_BALANCE,
		STARTING_BALANCE,
		VERIFY_BALANCE,
		NON_MODIFY_BALANCE
	}balancingStatus;

	enum OrderingStatus
	{
		NO_ORDERING,
		ORDERING,
		FORCE_ORDERING
	}orderingStatus;

	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	int numRows, numColumns, numElements, numThreads;
	double density;
	int whoAmI;
	int originalLowerBand, originalUpperBand;
	int lowerBand, upperBand;
	BzzMatrixSparseLocked A;
	BzzMatrixCoefficientsLocked B;
	BzzMatrixBand D;
	BzzMatrix M;
	BzzFactorizedBandGauss C;
	BzzVector bx;

	int getTime;
	double	startTime,
		elapsedTime,
		readingTime,
		orderingTime,
		buildBandMatrixTime,
		factorizationTime,
		solutionTime,
		totalTime;
	BzzVectorInt rB, cB;
	BzzVectorInt newRowsOrder;
	BzzVectorInt newColumnsOrder;
	void BuildBand(int iStart, int iEnd);
	int factorized;

public:
	//	============================================================================
	//	*************************< constructors >***********************************
	//	============================================================================
		// default constructor
	BzzFactorizedSparseLockedGauss(void);
	void operator()(BzzMatrix& AA, BzzVector& bb);
	void operator()(BzzMatrixSparseLocked* AA, BzzVector* bb);
	void operator()(BzzVector* bb);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedSparseLockedGauss(void) { countInScope--; }

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void) { BzzError("TODO PRINT"); }
	void GetTime(void) { getTime = 1; }
	void NoOrdering(void) { orderingStatus = NO_ORDERING; }
	void Ordering(void) { orderingStatus = ORDERING; }
	void ForceOrdering(void) { orderingStatus = FORCE_ORDERING; }

	//	============================================================================
	//	========================< Linear System Functions >=========================
	//	============================================================================
	//	void Solve(BzzVector *x);
	//	void Solve(BzzVector &b,BzzVector *x);
	//	double SystemConditioning(void);
};

#endif // BZZ_FACTORIZED_SPARSE_LOCKED_GAUSS_HPP