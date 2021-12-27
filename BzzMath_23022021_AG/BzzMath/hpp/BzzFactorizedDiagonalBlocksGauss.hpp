// BZZMATH: Release 7.0
/*
Linear systems with diagonal block structure
Example:
*xx.............
x*x.............
xx*.............
...*xxx.........
...x*xx.........
...xx*x.........
...xxx*.........
.......*xx......
.......x*x......
.......xx*......
..........*xx...
..........x*x...
..........xx*...
.............*xx
.............x*x
.............xx*
*/

//	==================< BzzFactorizedDiagonalBlocksGauss.hpp >===============
//	* BzzFactorizedDiagonalBlocksGauss: Class for solution of square linear	*
// * systems with diagonal block structure												*
// * Examples: BzzMath\Examples\BzzMathBasic\LinearSystems\							*
// * FactorizedDiagonalBlocksdGauss\FactorizedDiagonalBlocksdGauss.cpp	*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	09-2003	Date Written.

////////////////// Release 5.0

//	============================================================================
//	******* Functions for solving square linear diagonal block systems:			*
//	* Solve(&A,b,&x);																				*
//	* Solve(&A,&bx);																				*
//	* Solve(&A,B,&X);																				*
//	* Solve(&A,&BX);																				*
//	* TransposeSolve(&A,b,&x);																	*
//	* TransposeSolve(&A,&bx);																	*
//	* TransposeSolve(&A,B,&X);																	*
//	* TransposeSolve(&A,&BX);																	*
//	****************************************************************************

#ifndef BZZ_FACTORIZED_DOUBLE_DIAGONAL_BLOCK_HPP
#define BZZ_FACTORIZED_DOUBLE_DIAGONAL_BLOCK_HPP

//	============================================================================
//	=====================< class BzzFactorizedDiagonalBlocksGauss >==========================
//	============================================================================
class BzzFactorizedDiagonalBlocksGauss : public BzzBaseClass
{
	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================
	friend void Solve
	(BzzFactorizedDiagonalBlocksGauss* A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedDiagonalBlocksGauss* A, BzzVector& bx, BzzVector* x);

	friend void Solve
	(BzzFactorizedDiagonalBlocksGauss* A, BzzMatrix* BX);
	friend void Solve
	(BzzFactorizedDiagonalBlocksGauss* A, BzzMatrix& B, BzzMatrix* X);

	friend void TransposeSolve
	(BzzFactorizedDiagonalBlocksGauss* A, BzzVector* bx);
	friend void TransposeSolve
	(BzzFactorizedDiagonalBlocksGauss* A, BzzVector& bx, BzzVector* x);

	friend void TransposeSolve
	(BzzFactorizedDiagonalBlocksGauss* A, BzzMatrix* BX);
	friend void TransposeSolve
	(BzzFactorizedDiagonalBlocksGauss* A, BzzMatrix& B, BzzMatrix* X);

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
	BzzVectorInt lowerBand, upperBand, pivotBand, indx, blocksDimension, noPivot;
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
	BzzFactorizedDiagonalBlocksGauss(void);

	// copy constructor
//	BzzFactorizedDiagonalBlocksGauss(const BzzFactorizedDiagonalBlocksGauss &rval);

	//
	BzzFactorizedDiagonalBlocksGauss(int numV, int bloDim);
	void SetDimensions(int numV, int bloDim);
	BzzFactorizedDiagonalBlocksGauss(BzzVectorInt& diag);
	void SetDimensions(BzzVectorInt& diag);

	// assigns and receives vector values with control
	double& operator()(int i, int j);
	// assigns and receives vector values without control
	double* operator [] (int i)
	{
		return matrix[i];
	}

	// FILE ASCII
//	BzzFactorizedDiagonalBlocksGauss(char *filematrix);

//	============================================================================
//	********************************< destructor >******************************
//	============================================================================
	~BzzFactorizedDiagonalBlocksGauss(void) { Delete(this); }
	friend void Delete(BzzFactorizedDiagonalBlocksGauss* rval);

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
		// transforms a Matrix in Factorized
		// creates an independent matrix
	BzzFactorizedDiagonalBlocksGauss& operator = (BzzMatrixDiagonalBlocks& A);

	// transforms a Matrix in Factorized
	// Matrix is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixDiagonalBlocks* A, BzzFactorizedDiagonalBlocksGauss* B);
};

#endif // BZZ_FACTORIZED_DOUBLE_DIAGONAL_BLOCK_HPP