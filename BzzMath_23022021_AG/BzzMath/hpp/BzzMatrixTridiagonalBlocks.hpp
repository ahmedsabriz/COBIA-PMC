// BZZMATH: Release 7.0
/*
Matrices with tridiagonal block structure
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

//	========================< BzzMatrixTridiagonalBlocks.hpp >===================
//	* BzzMatrixTridiagonalBlocks: matrices with tridiagonal block structure		*
// * Examples: c:\bzzmath\examples\ExMatrixTridiagonalBlocks.cpp					*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	11-2002	Date Written.

////////////////// Release 5.0
//	08-2003	Added GetLowerElementInRows and GetUpperElementInRows.
//	08-2003	Added TProduct.

//	============================================================================

#ifndef BZZ_MATRIX_DOUBLE_TRIDIAGONAL_BLOCK_HPP
#define BZZ_MATRIX_DOUBLE_TRIDIAGONAL_BLOCK_HPP

//	============================================================================
//	=====================< class BzzMatrixTridiagonalBlocks >==========================
//	============================================================================
class BzzMatrixTridiagonalBlocks : public BzzBaseClass
{
	friend void Product(BzzMatrixTridiagonalBlocks& A, BzzVector& x,
		BzzVector* y);
	friend void TProduct(BzzMatrixTridiagonalBlocks& A, BzzVector& x,
		BzzVector* y);
	friend class BzzFactorizedTridiagonalBlocksGauss;
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixTridiagonalBlocks* A, BzzFactorizedTridiagonalBlocksGauss* B);
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
		numBlockDiagonal;

	double** matrix;

	BzzVectorInt lowerBand, upperBand;
	char singular;	// singular = 1 when singular
	int whoAmI;

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzMatrixTridiagonalBlocks(void);

	// copy constructor
//	BzzMatrixTridiagonalBlocks(const BzzMatrixTridiagonalBlocks &rval);

	//
	BzzMatrixTridiagonalBlocks(int numV, int bloDim);
	void SetDimensions(int numV, int bloDim);
	// assigns and receives vector values with control
	double& operator()(int i, int j);
	// assigns and receives vector values without control
	double* operator [] (int i)
	{
		return matrix[i];
	}

	// FILE ASCII
//	BzzMatrixTridiagonalBlocks(char *filematrix);

//	============================================================================
//	********************************< destructor >******************************
//	============================================================================
	~BzzMatrixTridiagonalBlocks(void) { Delete(this); }
	friend void Delete(BzzMatrixTridiagonalBlocks* rval);

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
	int NumBlockDiagonal(void) const
	{
		return numBlockDiagonal;
	}

	int WhoAmI(void) const { return whoAmI; }

	//	============================================================================
	//	****************************< assignment operators >************************
	//	============================================================================
	BzzMatrixTridiagonalBlocks& operator =
		(BzzMatrixTridiagonalBlocks& rval);

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
	//	double Condition(void);

	//	============================================================================
	//	==========================< Modifying functions	==========================
	//	============================================================================
	void SetToZeroAndReset(void);
};

#endif // BZZ_MATRIX_DOUBLE_TRIDIAGONAL_BLOCK_HPP