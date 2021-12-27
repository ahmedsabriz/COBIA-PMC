// BZZMATH: Release 7.0

//	==========================< BzzMatrixIntBand.hpp >==========================
//	* BzzMatrixIntBand: Class for operations with banded int matrices				*
// * Examples: BzzMath\Examples\BzzMathBasic\IntegerAlgebra\						*
// *				MatrixIntBand\MatrixIntBand.cpp													*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	12-1999	Date Written.

#ifndef BZZ_MATRIX_INT_BAND_HPP
#define BZZ_MATRIX_INT_BAND_HPP

class BzzVectorInt;

//	============================================================================
//	==========================< class BzzMatrixIntBand >===========================
//	============================================================================

class BzzMatrixIntBand : public BzzBaseClass
{
	friend class BzzSave;
	friend class BzzLoad;

private:
	static const char* const BZZ_ERROR;
	static int count; // per whoAmI
	int numRows, numColumns;
	int size;
	int	lowerBand,
		upperBand,
		totalBand;

	int** matrix;

	int whoAmI;

	// initialising constructors
	void Initialize(int m, int n, int low, int up);

	// re-initialising
	void ReInitialize(int m, int n, int low, int up);
	friend void ChangeDimensions(int m, int n, int low, int up, BzzMatrixIntBand* B);

	// deinitialisation
	void Deinitialize(void);

	// preparing assignments
	void PrepCopy(int rRows, int rColumns);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzMatrixIntBand(void);

	// copy constructor
	BzzMatrixIntBand(const BzzMatrixIntBand& rval);

	// sizing constructor
	BzzMatrixIntBand(int rows, int columns, int low, int up);

	// initialises from formatted File
	// BzzMatrixIntBand A("MAT.DAT");
	BzzMatrixIntBand(char* filematrix);

	// initialises from binary File
	// BzzMatrixIntBand A('*',"MAT.BIN"); See Save
	BzzMatrixIntBand(char, char* filematrix);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzMatrixIntBand(void);

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

	// assigns and receives values with control
	int& operator ()
		(int row, int column);

	// assigns and receives values without control
	int* operator [] (int r)
	{
		return matrix[r];
	}

	int LowerBand(void) { return lowerBand; }
	int UpperBand(void) { return upperBand; }
	//	============================================================================
	//	***************************< assignment operators >*************************
	//	============================================================================
	BzzMatrixIntBand& operator =
		(BzzMatrixIntBand& rval);

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);

	//	********************************< Save >************************************
	void Save(char* filematrix); // formatted
	void Save(char, char* filematrix);// binary

	//	*****************************< Max and Min >********************************
	//to have the position Max(im,jc)
	int Max(int* imax = 0, int* jmax = 0);
	int MaxAbs(int* imax = 0, int* jmax = 0);

	int Min(int* imin = 0, int* jmin = 0);
	int MinAbs(int* imin = 0, int* jmin = 0);

	//	============================================================================
	//	=====================< Modifying Functions >================================
	//	============================================================================
	friend void Delete(BzzMatrixIntBand* A);
	void SetToZeroAndReset(void);// eliminates all elements
	friend void Swap(BzzMatrixIntBand* A, BzzMatrixIntBand* B);
	void SetDiagonal(int j, BzzVectorInt& rval);
	void SetDiagonal(int j, int xf);
	// recovery of formatted Save
	friend void Load(BzzMatrixIntBand* A, char* filemat);
	// recovery of binary Save
	friend void Load(BzzMatrixIntBand* A, char, char* filemat);
};

#endif // BZZ_MATRIX_INT_BAND_HPP