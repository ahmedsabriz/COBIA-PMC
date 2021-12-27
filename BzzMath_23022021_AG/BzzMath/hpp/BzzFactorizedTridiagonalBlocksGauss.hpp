// BZZMATH: Release 7.0
/*
Linear systems with tridiagonal block structure
Example:
Example:
*xxxxx.........
x*xxxx.........
xx*xxx.........
xxx*xxxxx......
xxxx*xxxx......
xxxxx*xxx......
...xxx*xxxxx...
...xxxx*xxxx...
...xxxxx*xxx...
......xxx*xxxxx
......xxxx*xxxx
......xxxxx*xxx
.........xxx*xx
.........xxxx*x
.........xxxxx*
*/

//	==========< BzzFactorizedTridiagonalBlocksGauss.hpp >====================
//	* BzzFactorizedTridiagonalBlocksGauss: Class for solution of square linear			*
// * systems with tridiagonal block structure											*
// * Examples: c:\bzzmath\examples\ExFactTridiagonalBlocks.cpp						*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	08-2002	Date Written.

////////////////// Release 5.0
//	08-2003	Added TransposeSolve.
//	08-2003	Added ConditionNumber.
//	08-2003	Added Determinant.

////////////////// Release 6.0
//	07-2009	BzzFactorizedGauss parallelized wth OpenMP.

//	============================================================================
//	******* Functions for solving square tridiagonal block systems:				*
//	* Solve(&A,b,&x);																				*
//	* Solve(&A,&bx);																				*
//	* Solve(&A,B,&X);																				*
//	* Solve(&A,&BX);																				*
//	* TransposeSolve(&A,b,&x);																	*
//	* TransposeSolve(&A,&bx);																	*
//	* TransposeSolve(&A,B,&X);																	*
//	* TransposeSolve(&A,&BX);																	*
//	****************************************************************************

#ifndef BZZ_FACTORIZED_DOUBLE_TRIDIAGONAL_BLOCK_HPP
#define BZZ_FACTORIZED_DOUBLE_TRIDIAGONAL_BLOCK_HPP

//	============================================================================
//	=====================< class BzzFactorizedTridiagonalBlocksGauss >==========================
//	============================================================================
class BzzFactorizedTridiagonalBlocksGauss : public BzzBaseClass
{
	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================
	friend void Solve
	(BzzFactorizedTridiagonalBlocksGauss* A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedTridiagonalBlocksGauss* A, BzzVector& bx, BzzVector* x);

	friend void Solve
	(BzzFactorizedTridiagonalBlocksGauss* A, BzzMatrix* BX);
	friend void Solve
	(BzzFactorizedTridiagonalBlocksGauss* A, BzzMatrix& B, BzzMatrix* X);

	friend void TransposeSolve
	(BzzFactorizedTridiagonalBlocksGauss* A, BzzVector* bx);
	friend void TransposeSolve
	(BzzFactorizedTridiagonalBlocksGauss* A, BzzVector& bx, BzzVector* x);

	friend void TransposeSolve
	(BzzFactorizedTridiagonalBlocksGauss* A, BzzMatrix* BX);
	friend void TransposeSolve
	(BzzFactorizedTridiagonalBlocksGauss* A, BzzMatrix& B, BzzMatrix* X);

protected:
	enum OdebFactorizationStatus
	{
		MATRIX_NON_AVAILABLE,
		UNFACTORIZED,
		FACTORIZED
	}factorizationStatus;

	static const char* const BZZ_ERROR;
	static const char* const BZZ_ERR_DESTROYED;
	static int count; // for whoAmI
	static int countInScope;

	int	numRows,
		numColumns,
		size,
		blockDimensions,
		numBlockDiagonal;

	double** matrix;

	double norm;
	BzzVectorInt lowerBand, upperBand, pivotBand, indx;
	char singular;	// singular = 1 when singular
	int whoAmI;

	void Factorization(void);
	void Solution(BzzVector* bx);
	void TransposeSolution(BzzVector* bx);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedTridiagonalBlocksGauss(void);

	// copy constructor
//	BzzFactorizedTridiagonalBlocksGauss(const BzzFactorizedTridiagonalBlocksGauss &rval);

	//
	BzzFactorizedTridiagonalBlocksGauss(int numV, int bloDim);
	void SetDimensions(int numV, int bloDim);
	// assigns and receives vector values with control
	double& operator()(int i, int j);
	// assigns and receives vector values without control
	double* operator [] (int i)
	{
		return matrix[i];
	}

	// FILE ASCII
//	BzzFactorizedTridiagonalBlocksGauss(char *filematrix);

//	============================================================================
//	********************************< destructor >******************************
//	============================================================================
	~BzzFactorizedTridiagonalBlocksGauss(void) { Delete(this); }
	friend void Delete(BzzFactorizedTridiagonalBlocksGauss* rval);

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

	int BlockDimensions(void) const
	{
		return blockDimensions;
	}

	int BumBlockDiagonal(void) const
	{
		return numBlockDiagonal;
	}

	int WhoAmI(void) const { return whoAmI; }

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================
	BzzVectorInt GetLowerElementInRows(void);
	void GetLowerElementInRows(BzzVectorInt* low);
	BzzVectorInt GetUpperElementInRows(void);
	void GetUpperElementInRows(BzzVectorInt* up);

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
	void BzzPrintStructure(void);
	int CountElements(void)
	{
		return size - 1;
	}
	double ConditionNumber(void);
	double Determinant(void);

	//	============================================================================
	//	==========================< Modifying functions	==========================
	//	============================================================================
	void SetToZeroAndReset(void);

	//	============================================================================
	//	***********************< assignment operators >*****************************
	//	============================================================================
		// transforms a BzzMatrix in Factorized
		// creates an independent matrix
	BzzFactorizedTridiagonalBlocksGauss& operator = (BzzMatrixTridiagonalBlocks& A);

	// transforms a BzzMatrix in Factorized
	// BzzMatrix is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixTridiagonalBlocks* A, BzzFactorizedTridiagonalBlocksGauss* B);

	//	============================================================================
	//	*********************< Special function for ODE >***************************
	//	============================================================================
	friend void ProductForOde(double ch, BzzFactorizedTridiagonalBlocksGauss* A);

	//	============================================================================
	//	*********************< Special function for DAE >***************************
	//	============================================================================
	friend void ProductForDae(double ch, double c, BzzVectorInt& iDer,
		BzzFactorizedTridiagonalBlocksGauss* A);
};

#endif // BZZ_FACTORIZED_DOUBLE_TRIDIAGONAL_BLOCK_HPP