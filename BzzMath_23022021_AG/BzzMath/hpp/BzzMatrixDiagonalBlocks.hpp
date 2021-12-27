// BZZMATH: Release 7.0
/*
Matrices with diagonal block structure
Example:
*.........
.*xxx.....
.x*xx.....
.xx*x.....
.xxx*.....
.....*xx..
.....x*x..
.....xx*..
........*x
........x*

Example:
*xxx........
x*xx........
....xxxx....
....xxxx....
........xxxx
........xxxx

Example:
*x..
x*..
xx..
xx..
..xx
..xx
..xx
..xx

*/

//	=======================< BzzMatrixDiagonalBlocks.hpp >=================
//	* BzzMatrixDiagonalBlocks: matrices with tridiagonal block structure	*
// * Examples: BzzMath\examples\BzzMathBasic\LinearAlgebra\							*
// *				MatrixDiagonalBlocks\MatrixDiagonalBlocks.cpp			*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris

////////////////// Release 5.0
//	09-2003	Date Written.
//	02-2005	Added DiagonalProduct function.

//	============================================================================

#ifndef BZZ_MATRIX_DOUBLE_DIAGONAL_BLOCK_HPP
#define BZZ_MATRIX_DOUBLE_DIAGONAL_BLOCK_HPP

//	============================================================================
//	=====================< class BzzMatrixDiagonalBlocks >==========================
//	============================================================================
class BzzMatrixDiagonalBlocks : public BzzBaseClass
{
	friend void DiagonalProduct(char* file, BzzVector& x, BzzVector* y);
	friend void Product(double v, BzzMatrixDiagonalBlocks* A);
	friend void Product(double v, BzzMatrixDiagonalBlocks& B,
		BzzMatrixDiagonalBlocks* A);
	friend void Product(BzzMatrixDiagonalBlocks& A, BzzVector& x,
		BzzVector* y);
	friend void Product(BzzMatrixDiagonalBlocks& A,
		BzzMatrixDiagonalBlocks& B,
		BzzMatrixDiagonalBlocks* C);
	friend void TProduct(BzzMatrixDiagonalBlocks& A, BzzVector& x,
		BzzVector* y);
	friend class BzzFactorizedDiagonalBlocksGauss;
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixDiagonalBlocks* A, BzzFactorizedDiagonalBlocksGauss* B);
	friend void IProduct(BzzFactorizedDiagonalBlocksGauss* A,
		BzzMatrixDiagonalBlocks& B,
		BzzMatrixDiagonalBlocks* C);
	friend class BzzDaeSparse;
	friend class BzzOdeSparseStiff;
	friend class BzzDaeSparseObject;
	friend class BzzOdeSparseStiffObject;

private:
	static const char* const BZZ_ERROR;
	static const char* const BZZ_ERR_DESTROYED;
	static int count; // for whoAmI
	static int countInScope;

	int	numRows,
		numColumns,
		size,
		blockDimensions,
		blockDimensionsColumn,
		numBlockDiagonal;

	double** matrix;

	BzzVectorInt lowerBand, upperBand, blocksDimension;
	int whoAmI;
public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzMatrixDiagonalBlocks(void);

	// copy constructor
	BzzMatrixDiagonalBlocks(BzzMatrixDiagonalBlocks& rval);

	//
	BzzMatrixDiagonalBlocks(int numV, int bloDim);
	void SetDimensions(int numV, int bloDim);

	BzzMatrixDiagonalBlocks(BzzVectorInt& diag);
	void SetDimensions(BzzVectorInt& diag);

	BzzMatrixDiagonalBlocks(int numR, int numC, int bloDimR, int bloDimC);
	void SetDimensions(int numR, int numC, int bloDimR, int bloDimC);

	// assigns and receives vector values with control
	double& operator()(int i, int j);
	// assigns and receives vector values without control
	double* operator [] (int i)
	{
		return matrix[i];
	}

	// FILE ASCII
//	BzzMatrixDiagonalBlocks(char *filematrix);

//	============================================================================
//	********************************< destructor >******************************
//	============================================================================
	~BzzMatrixDiagonalBlocks(void) { Delete(this); }
	friend void Delete(BzzMatrixDiagonalBlocks* rval);

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
	int BlockRowDimensions(void) const
	{
		return blockDimensions;
	}
	int BlockColumnDimensions(void) const
	{
		return blockDimensionsColumn;
	}
	int NumBlockDiagonal(void) const
	{
		return numBlockDiagonal;
	}

	int WhoAmI(void) const { return whoAmI; }

	//	============================================================================
	//	***************************< assignment operators >*************************
	//	============================================================================
	BzzMatrixDiagonalBlocks& operator =
		(BzzMatrixDiagonalBlocks& rval);

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================
	BzzVectorInt GetLowerElementInRows(void);
	void GetLowerElementInRows(BzzVectorInt* low);
	BzzVectorInt GetUpperElementInRows(void);
	void GetUpperElementInRows(BzzVectorInt* up);
	BzzVectorInt GetBlocksDimensions(void)
	{
		return blocksDimension;
	}

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
	void BzzPrintStructure(void);
	int CountElements(void)
	{
		return size - 1;
	}
	//	double Condition(void);

	//	============================================================================
	//	==========================< Modifying functions	==========================
	//	============================================================================
	void SetToZeroAndReset(void);
};

#endif // BZZ_MATRIX_DOUBLE_DIAGONAL_BLOCK_HPP