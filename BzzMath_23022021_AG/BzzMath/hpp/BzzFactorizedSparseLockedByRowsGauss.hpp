// BZZMATH: Release 7.0

///////////////////////////////////////
// TODO:
// ConditionNumber
///////////////////////////////////////

//	================< BzzFactorizedSparseLockedByRowsGauss.hpp >========================
//	* BzzFactorizedSparseLockedByRowsGauss: Class for solution of square linear			*
// * sparse systems																				*
// * Description:																					*
// * Examples: BzzMath\examples\BzzMathBasic\LinearSystems\							*
// *				FactorizedSparseLockedGauss\FactorizedSparseLockedGauss.cpp			*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	04-2002	Date Written.
//	03-2012	Added BzzFactorizedGaussAttic class.
//	04-2012	Added BzzFactorizedGaussSparseAttic class.

//	============================================================================
//	******* Functions for solving square linear sparse systems:						*
//	* Solve(&A,b,&x);																				*
//	* Solve(&A,&bx);																				*
//	****************************************************************************

#ifndef BZZ_FACTORIZED_DOUBLE_SPARSE_LOCKED_BYROWS_GAUSS_HPP
#define BZZ_FACTORIZED_DOUBLE_SPARSE_LOCKED_BYROWS_GAUSS_HPP

class BzzMatrixSparse;
struct ElementBzzMatrixSparse;

//	============================================================================
//	==============< class BzzFactorizedSparseLockedByRowsGauss >================
//	============================================================================
class BzzFactorizedSparseLockedByRowsGauss : public BzzBaseClass
{
	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================
	friend void Solve
	(BzzFactorizedSparseLockedByRowsGauss* A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedSparseLockedByRowsGauss* A, BzzVector& bx, BzzVector* x);

	friend void Solve
	(BzzFactorizedSparseLockedByRowsGauss* A, BzzMatrix* bx);
	friend void Solve
	(BzzFactorizedSparseLockedByRowsGauss* A, BzzMatrix& bx, BzzMatrix* x);

	//friend void TransposeSolve
	//		 (BzzFactorizedSparseLockedByRowsGauss *A,BzzVector *bx);
	//friend void TransposeSolve
	//		 (BzzFactorizedSparseLockedByRowsGauss *A,BzzVector &bx,BzzVector *x);

private:
	enum LockedSparseFactorizationStatus
	{
		MATRIX_NON_AVAILABLE,
		UNFACTORIZED,
		FACTORIZED
	}factorizationStatus;

	static const char* const BZZ_ERROR;
	static const char* const BZZ_ERR_DESTROYED;
	static int count; // for whoAmI
	static int countInScope;

	BzzMatrixSparse S, O;
	int numRows, numColumns;

	BzzVector* l, * r, * a;
	BzzVector norm1;
	BzzVector w;
	BzzVectorInt* il, * ir, * ia;
	BzzVectorInt sizea, sizel, sizer;
	BzzVectorInt indx;

	char singular;	// singular = 1 when singular
	char ordering;
	char pivoting;
	int whoAmI;
	BzzVector norm;
	//	BzzVectorInt ordRows,ordColumns,pivotRows;
	BzzVectorInt ordRows, ordColumns;

	void Factorization(void);
	void Solution(BzzVector* bx);
	void Solution(BzzMatrix* bx);

	//	void TransposeSolution(BzzVector *bx);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedSparseLockedByRowsGauss(void);

	// copy constructor
//	BzzFactorizedSparseLockedByRowsGauss(const BzzFactorizedSparseLockedByRowsGauss &rval);

		//	from BzzMatrixSparse
	BzzFactorizedSparseLockedByRowsGauss
	(BzzMatrixSparse& rval);

	// FILE ASCII
	BzzFactorizedSparseLockedByRowsGauss(char* filematrix);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedSparseLockedByRowsGauss(void);
	void Deinitialize(void);

	//	============================================================================
	//	*****************************< Access functions >***************************
	//	============================================================================
		// number of rows
	int Rows(void) const
	{
		return numRows;
	}

	// number of columns
	int Columns(void) const
	{
		return numColumns;
	}

	int WhoAmI(void) const { return whoAmI; }

	//	============================================================================
	//	***************************< assignment operators >*************************
	//	============================================================================
	void operator =
		(BzzMatrixSparse& rval);

	// transforms a BzzMatrixSparse in BzzFactorizedSparseLockedByRowsGauss
	// BzzMatrixSparse is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixSparse* lval, BzzFactorizedSparseLockedByRowsGauss* rval);

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
	//	void BzzPrintExistingElements(void);
	//	void BzzPrintStructure(void);
	//	int CountElements(void);
	//	double Condition(void);

	//	============================================================================
	//	=====================< Modifying Functions >================================
	//	============================================================================
	void NonOrdering(void)
	{
		if (factorizationStatus == UNFACTORIZED)
			ordering = 0;
	}
	void NonPivoting(void)
	{
		if (factorizationStatus == UNFACTORIZED)
			pivoting = 0;
	}
	void Ordering(void)
	{
		if (factorizationStatus == UNFACTORIZED)
			ordering = 1;
	}
	void Pivoting(void)
	{
		BzzError("TODO Pivoting in BzzFactorizedSparseLockedByRowsGauss");
		if (factorizationStatus == UNFACTORIZED)
			pivoting = 1;
	}
};

//	============================================================================
//	===================< class BzzFactorizedGaussAttic >========================
//	============================================================================
// Che cosa DEVE FARE:
// dato un sistema sottodimensionato Ax = b
// trova le xd e le xu Ad.xd = b - Au.xu
// e date le x trova le xd
//
//
//
//
//

class BzzFactorizedGaussAttic : public BzzBaseClass
{
private:
	enum BzzGaussAtticFactorizationStatus
	{
		MATRIX_NON_AVAILABLE,
		UNFACTORIZED,
		FACTORIZED
	}bzzGaussAtticFactorizationStatus;

	enum BzzGaussAtticSparsity
	{
		DENSE,
		//		DENSE_SPARSE,
		//		SPARSE_DENSE,
		SPARSE
	}bzzGaussAtticSparsity;

	enum BzzGaussAtticBalancing
	{
		NON_AVAILABLE_BALANCE,
		STARTING_BALANCE,
		NON_MODIFY_BALANCE,
		VERIFY_BALANCE
	}bzzGaussAtticBalancing;

	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	static int const MAX_ROWS;
	static int const MAX_ROWS_COLUMNS;
	static double const ZERO_FRACTION;
	static double const EPS_PIVOT;

	double epsPivot;

	int	numOriginalRows, numRows,
		numColumns,
		noBalancing,
		ordering,
		pivoting,
		whoAmI;

	BzzVectorInt	ipColumns,
		jDependent,
		jIndependent,
		deletedEquations,
		linearlyIndependentEquations;

	BzzVector	bForBalancing,
		b,
		bLast,
		rowsCoefficients, columnsCoefficients,
		rad, urad;
	BzzMatrix	originalDenseA,
		modifiedDenseA,
		A,
		factorizedDenseAd,
		denseAu, denseAd,
		denseAm1dAu,
		denseN;
private:
	void DenseToSparse(void);
	void SparseToDense(void);
	void DenseBalancing(void);
	void SparseBalancing(void);
	void DenseFactorization(void);
	void DenseFactorization(BzzVector& weightX);
	BzzVector weightX;
	void SparseFactorization(void);
	void DenseSolution(BzzVector* x);
	void SparseSolution(BzzVector* x);
	void DenseTransposeSolution(BzzVector& c, BzzVector* y);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedGaussAttic(void);
	void operator()(BzzMatrix& AA, BzzVector& bb);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedGaussAttic(void) {};
	//	void Deinitialize(void);

	//	============================================================================
	//	*****************************< Access functions >***************************
	//	============================================================================
	void SetEpsPivot(double ep)
	{
		epsPivot = ep;
	}
	// number of rows
	int OriginalRows(void) { return numOriginalRows; }

	int Rows(void) { return numRows; }

	// number of columns
	int Columns(void) { return numColumns; }

	int WhoAmI(void) { return whoAmI; }

	void SetNoBalancing(void)
	{
		if (bzzGaussAtticFactorizationStatus != FACTORIZED)
			noBalancing = 1; // non bilancia
	}

	void GetDependentVariables(BzzVectorInt* jDep)
	{
		if (bzzGaussAtticFactorizationStatus == FACTORIZED)
			*jDep = jDependent;
	}
	void GetIndependentVariables(BzzVectorInt* jIndep)
	{
		if (bzzGaussAtticFactorizationStatus == FACTORIZED)
			*jIndep = jIndependent;
	}

	void GetLinearlyDependentEquations(BzzVectorInt* iDep)
	{
		if (bzzGaussAtticFactorizationStatus == FACTORIZED)
			*iDep = deletedEquations;
	}

	void GetLinearlyDependentAndIndependentEquations(BzzVectorInt* iDep, BzzVectorInt* iIndep);

	int GetNumberOfLinearlyIndependentEquations(void)
	{
		if (bzzGaussAtticFactorizationStatus == FACTORIZED)
			return numRows;
	}

	void GetLinearlyIndependentEquations(BzzVectorInt* iIndep);

	void GetMatricesAdAu(BzzMatrix* Ad, BzzMatrix* Au);
	void GetModifiedMatrixA(BzzMatrix* modifiedMatrixA);
	void GetMatrixAm1dAu(BzzMatrix* Am1dAu);
	void GetMatrixN(BzzMatrix* N);
	void GetFactorizedMatrix(BzzMatrix* AA);

	//	============================================================================
	//	***************************< assignment operators >*************************
	//	============================================================================

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
	//	void Get

	//	============================================================================
	//	============================< Modifying functions	=========================
	//	============================================================================
	//	void DeleteRow(int row);
	//	void AddRow(BzzVector &r);

	//	void NonOrdering(void)
	//		{
	//		if(factorizationStatus == UNFACTORIZED)
	//			ordering = 0;
	//		}
	//	void NonPivoting(void)
	//		{
	//		if(factorizationStatus == UNFACTORIZED)
	//			pivoting = 0;
	//		}
	//	void Ordering(void)
	//		{
	//		if(factorizationStatus == UNFACTORIZED)
	//			ordering = 1;
	//		}
	//	void Pivoting(void)
	//		{
	//		BzzError("TODO Pivoting in BzzFactorizedGaussAttic");
	//		if(factorizationStatus == UNFACTORIZED)
	//			pivoting = 1;
	//		}

	//	============================================================================
	//	========================< Linear System Functions >=========================
	//	============================================================================
	void Solve(BzzVector* x);
	void Solve(BzzVector& b, BzzVector* x);
	void Solve(BzzVector* x, BzzVector& weightX);
	void Solve(BzzVector& b, BzzVector* x, BzzVector& weightX);

	void Factorize(BzzVector* x);
	void Factorize(BzzVector& b, BzzVector* x);
	void Factorize(BzzVector* x, BzzVector& weightX);
	void Factorize(BzzVector& b, BzzVector* x, BzzVector& weightX);

	void TransposeSolve(BzzVector& c, BzzVector* y);
	double SystemConditioning(void);
};

//	============================================================================
//	===================< class BzzFactorizedGaussSparseAttic >========================
//	============================================================================
// Che cosa DEVE FARE:
// dato un sistema sottodimensionato Ax = b
// trova le xd e le xu Ad.xd = b - Au.xu
// e date le x trova le xd
//
//
//
//
//

class BzzFactorizedGaussSparseAttic : public BzzBaseClass
{
private:
	enum BzzGaussAtticSparseFactorizationStatus
	{
		MATRIX_NON_AVAILABLE,
		UNFACTORIZED,
		FACTORIZED
	}bzzGaussAtticSparseFactorizationStatus;

	enum BzzGaussAtticSparseBalancing
	{
		NON_AVAILABLE_BALANCE,
		STARTING_BALANCE,
		NON_MODIFY_BALANCE,
		VERIFY_BALANCE
	}bzzGaussAtticSparseBalancing;

	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	static int const MAX_ROWS;
	static int const MAX_ROWS_COLUMNS;
	static double const ZERO_FRACTION;
	static double const EPS_PIVOT;

	double epsPivot;

	int	numOriginalEquations, numEquations,
		numOriginalVariables, numVariables,
		numColumns,
		numOriginalElements, numElements,
		noBalancing,
		whoAmI,
		numConstraints;

	FILE* fileGaussSparseAttic;

	BzzVectorInt	rOriginal, r,
		cOriginal, c;

	//						ipColumns,
	//						jDependent,
	//						jIndependent,
	//						deletedEquations,
	//						linearlyIndependentEquations;

	BzzVector	vOriginal, v,
		b,
		xOriginal, x,
		//rowsCoefficients,columnsCoefficients,
		rad, urad;

	void InitializeGaussSparseAttic(void);
	/*
		BzzMatrix	originalDenseA,
						modifiedDenseA,
						A,
						factorizedDenseAd,
						denseAu,denseAd,
						denseAm1dAu,
						denseN;

		void DenseToSparse(void);
		void SparseToDense(void);
		void DenseBalancing(void);
		void SparseBalancing(void);
		void DenseFactorization(void);
		void SparseFactorization(void);
		void DenseSolution(BzzVector *x);
		void SparseSolution(BzzVector *x);
	*/
public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedGaussSparseAttic(void);
	void operator()(BzzMatrix* EE, BzzVector* ee);
	void operator()(int m, int n, BzzVectorInt* rr, BzzVectorInt* cc, BzzVector* vv,
		BzzVector* bb);
	void operator()(char* fileLP);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedGaussSparseAttic(void) {};
	//	void Deinitialize(void);
	/*
	//	============================================================================
	//	*****************************< Access functions >***************************
	//	============================================================================
		// number of rows
		int OriginalRows(void){return numOriginalRows;}

		int Rows(void){return numRows;}

		// number of columns
		int Columns(void){return numColumns;}

		int WhoAmI(void){return whoAmI;}

		void SetNoBalancing(void)
			{
			if(bzzGaussAtticFactorizationStatus != FACTORIZED)
				noBalancing = 1; // non bilancia
			}

		void GetDependentVariables(BzzVectorInt *jDep)
			{
			if(bzzGaussAtticFactorizationStatus == FACTORIZED)
				*jDep = jDependent;
			}
		void GetIndependentVariables(BzzVectorInt *jIndep)
			{
			if(bzzGaussAtticFactorizationStatus == FACTORIZED)
				*jIndep = jIndependent;
			}

		void GetLinearlyDependentEquations(BzzVectorInt *iDep)
			{
			if(bzzGaussAtticFactorizationStatus == FACTORIZED)
				*iDep = deletedEquations;
			}

		void GetLinearlyDependentAndIndependentEquations(BzzVectorInt *iDep,BzzVectorInt *iIndep);

		int GetNumberOfLinearlyIndependentEquations(void)
			{
			if(bzzGaussAtticFactorizationStatus == FACTORIZED)
				return numRows;
			}

		void GetLinearlyIndependentEquations(BzzVectorInt *iIndep);

		void GetMatricesAdAu(BzzMatrix *Ad,BzzMatrix *Au);
		void GetModifiedMatrixA(BzzMatrix *modifiedMatrixA);
		void GetMatrixAm1dAu(BzzMatrix *Am1dAu);
		void GetMatrixN(BzzMatrix *N);

	//	============================================================================
	//	***************************< assignment operators >*************************
	//	============================================================================

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================
	*/
	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
	/*
	//	============================================================================
	//	============================< Modifying functions	=========================
	//	============================================================================
	//	void DeleteRow(int row);
	//	void AddRow(BzzVector &r);

	//	void NonOrdering(void)
	//		{
	//		if(factorizationStatus == UNFACTORIZED)
	//			ordering = 0;
	//		}
	//	void NonPivoting(void)
	//		{
	//		if(factorizationStatus == UNFACTORIZED)
	//			pivoting = 0;
	//		}
	//	void Ordering(void)
	//		{
	//		if(factorizationStatus == UNFACTORIZED)
	//			ordering = 1;
	//		}
	//	void Pivoting(void)
	//		{
	//		BzzError("TODO Pivoting in BzzFactorizedGaussSparseAttic");
	//		if(factorizationStatus == UNFACTORIZED)
	//			pivoting = 1;
	//		}

	//	============================================================================
	//	========================< Linear System Functions >=========================
	//	============================================================================
		void Solve(BzzVector *x);
		void Solve(BzzVector &b,BzzVector *x);
		double SystemConditioning(void);
	*/
};

#endif // BZZ_FACTORIZED_DOUBLE_SPARSE_LOCKED_BYROWS_GAUSS_HPP