// BZZMATH: Release 7.0

//	=========================< BzzMatrixInt.hpp >===============================
//	* Class BzzMatrixInt Class for operations int matrices							*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 3)					*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\BzzMathBasic\IntegerAlgebra\					*
// *				MatrixInt\MatrixInnt.cpp													*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	02-1997	Date Written.
//	07-1998	Added AppendRow, AppendColumn functions.

////////////////// Release 4.0
//	10-2001	Modifyed AppendRow, AppendColumn functions.
//	10-2001	Modifyed InsertRow, InsertColumn functions.

////////////////// Release 5.0
// 10-2003	Added MaxColumn,MaxAbsColumn,MinColumn,MinAbsColumn functions.
// 10-2003	Added MaxRow,MaxAbsRow,MinRow,MinAbsRow functions.
// 11-2003	Added SetToZero function.
//	07-2004	Added AppendRowsFromMatrix function.
//	01-2005	Added PrintByRows and PrintTranspose functions.
//	01-2005	Added MaxRows,MaxAbsRows,MaxColumns,MaxAbsColumns functions.
//	01-2005	Added MinRows,MinAbsRows,MinColumns,MinAbsColumns functions.
//	02-2006	Added CopyDataFromVector function.
//	02-2004	Added UseMatrixAsVector function.
//	02-2004	Added UseMatrixRowAsVector function.
//	02-2006	Added GetRowHandle function.
//	10-2006	Added ReorderByColumns function.
//	12-2006	Added constructor from file as for MatrixSparse.
//	06-2008	Added SwapMaxMinOrder function.
//	06-2008	Added AntiTranspose function.
//	06-2008	Added HorizontalReflection and VerticalReflection functions.

////////////////// Release 6.0
//	10-2008	Added MoveRowFromkToj function.
//	10-2008	Added MoveColumnromkToj function.

//	============================================================================
//	******	BzzMatrixInt constructors:														*
//	* BzzMatrixInt A; // default																*
//	* BzzMatrixInt A = B; // copy-initializer												*
//	* BzzMatrixInt A(m,n); // sizes and places at 0										*
//	* BzzMatrixInt A(2,3,1,2,3,4,5,6); // matrix 2X3									*
//	* int x[6]={1,2,3,4,5,6};																	*
//	* BzzMatrixInt A(2,3,x); // from array													*
//	* BzzMatrixInt A("MAT.DAT"); // Formatted	file										*
//	* BzzMatrixInt A("MAT.MTX",'*'); // Formatted file as for Sparse				*
//	* BzzMatrixInt A('*',MAT.BIN"); // Binary file										*
//	* BzzMatrixInt A(m,n,B); // submatrix of B											*
//	* BzzMatrixInt A(m,n,i,j,B);// submatrix of B from i,j							*
//	****************************************************************************
//	***** Access functions:																		*
//	* i = A.Rows(); // numRows																	*
//	* i = A.Columns(); // numColumns															*
//	* who = A.WhoAmI();																			*
//	* xf = A.GetValue(i,j);																		*
//	* v = GetRow(i);																				*
//	* v = GetColumn(j);																			*
//	* v = GetDiagonal(j); // j = 0 principal positive at right						*
//	* xf = A(i,j);																					*
//	* A(i,j) = xf;																					*
//	* xf = A[i][j];																				*
//	* A[i][j] = xf;																				*
//	* A.SetValue(i,j,xf);																		*
//	* A.SetRow(i,v);																				*
//	* A.SetRow(i,xf);																				*
//	* A.SetColumn(j,v);																			*
//	* A.SetColumn(j,xf);																			*
//	* A.SetDiagonal(j,xf); // j = 0 principal positive at right						*
//	* A.SetDiagonal(j,v);  // j = 0 principal positive at right						*
//	* A.SetBzzMatrix(xf);																		*
//	*	int count = BzzMatrixInt::ObjectCount();											*
//	*	int countInScope = BzzMatrixInt::ObjectCountInScope();						*
//	****************************************************************************
//	***** Assignments:																			*
//	* A = B;	// BzzMatrixInt																	*
//	* A = v;	// BzzVectorInt																	*
//	* A = 3;	// int																				*
//	****************************************************************************
//	*****	Operators for composing matrices:												*
//	* A = B&&C; // Adds C beneath B															*
//	* A = B||C; // Adds C onto the side of B												*
//	****************************************************************************
//	*****	Operators for tests:																	*
//	* if(A == B)																					*
//	* if(A != B)																					*
//	****************************************************************************
//	* Other functions:																			*
//	* A.BzzPrint("Comment"); // always														*
//	* A.BzzMessage("Comment"); // only if bzzMessageYesNo == 'Y'					*
//	* xf = A.Max(&imax,&jmax);																	*
//	* xf = A.MaxAbs(&imax,&jmax);																*
//	* xf = A.Min(&imin,&jmin);																	*
//	* xf = A.MinAbs(&imin,&jmin);																*
//	* xf = A.MaxColumn(j,&imax);																*
//	* xf = A.MaxAbsColumn(j,&imax);															*
//	* xf = A.MinColumn(j,&imin);																*
//	* xf = A.MinAbsColumn(j,&imin)															*
//	* xf = A.MaxRow(i,&jmax);																	*
//	* xf = A.MaxAbsRow(i,&jmax);																*
//	* xf = A.MinRow(i,&jmin);																	*
//	* xf = A.MinAbsRow(i,&jmin)																*
//	* A.Save("MAT.DAT");	// formatted														*
//	* A.Save('*',"MAT.DAT");	// unformatted												*
//	* Load(&A,"MAT.DAT");																		*
//	* Load(&A,'*',"MAT.DAT");																	*
//	* A.InsertRow(i,v);																			*
//	* A.InsertRow(i,xf);																			*
//	* A.InsertColumn(j,v);																		*
//	* A.InsertColumn(j,xf);																		*
//	* A.AppendRow(v);																				*
//	* A.AppendRow(xf);																			*
//	* A.AppendRowsFromMatrix(iv,B);															*
//	* A.AppendColumn(v);																			*
//	* A.AppendColumn(xf);																		*
//	* A.DeleteRow(i);																				*
//	* A.DeleteRows(iv);																			*
//	* A.DeleteColumn(j);																			*
//	* A.DeleteColumns(jw);																		*
//	* Delete(&A);																					*
//	* A.SetToZero();																				*
//	* ChangeDimensions(newr,newc,&A);														*
//	* Swap(&A,&B);																					*
//	* A.SwapRows(i,j);																			*
//	* A.SwapColumns(i,j);																		*
//	* SwapMaxMinOrder(&A);																		*
//	* Transpose(&A);																				*
//	* AntiTranspose(&A);																			*
//	* HorizontalReflection(&A);																*
//	* VerticalReflection(&A);																	*
//	* A.UseMatrixAsVector(&a);																	*
//	* A.UseMatrixRowAsVector(row,&a);														*
//	****************************************************************************

#ifndef BZZ_MATRIX_INT_HPP
#define BZZ_MATRIX_INT_HPP

// preventive statements
class BzzVectorInt;

//	============================================================================
//	============================< class BzzMatrixInt >==========================
//	============================================================================
class BzzMatrixInt : public BzzBaseClass
{
	friend class BzzVectorInt;
	friend class BzzSave;
	friend class BzzLoad;

private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	int** matrix;
	int numRows, numColumns;
	int size;
	int whoAmI;
	char shadow;

	// initialisation constructors
	void Initialize(int mrows, int ncolumns);

	// deinitialisation
	void Deinitialize(void);

	// prepares assignments
	void PrepCopy(int rRows, int rColumns);

	// private constructor BzzMatrixInt A type('*',3,5);
	BzzMatrixInt(char, int mrows, int columns);

public:
	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default // BzzMatrixInt A;
	BzzMatrixInt(void);

	// copy-initializer // BzzMatrixInt A = B;
#if BZZ_COMPILER == 101
	BzzMatrixInt(const BzzMatrixInt& rval);
#else
	BzzMatrixInt(BzzMatrixInt& rval);
#endif

	// sizes and initialises at 0 //BzzMatrixInt A(3,5);
	BzzMatrixInt(int rows, int columns);

	// sizes and initialises
	// BzzMatrixInt A(2,3,1,2,3,4,5,6);
	BzzMatrixInt(int rows, int columns, int a11, ...);

	// initialises from array
	// int w[6]={1,2,3,4,5,6}; BzzMatrixInt A(3,5,w);
	BzzMatrixInt(int rows, int columns, int* initvalues);

	// initialises from formatted File
	// BzzMatrixInt A("MATINT.DAT");
	BzzMatrixInt(char* filematrix);

	// initialises from binary File
	// BzzMatrixInt A('*',"MATINT.BIN"); See Save
	BzzMatrixInt(char, char* filematrix);

	// makes a submatrix of rows,columns
	BzzMatrixInt(int rows, int columns, const BzzMatrixInt& rval);

	// likewise starting from irow,jcol
	BzzMatrixInt(int rows, int columns, int irow, int jcol,
		const BzzMatrixInt& rval);

	// from BzzVectorInt
	BzzMatrixInt(const BzzVectorInt& rval);

	//	============================================================================
	//	****************************< destructor >**********************************
	//	============================================================================
	~BzzMatrixInt(void);

	//	============================================================================
	//	********************< Non-modifying access functions >**********************
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
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

	// receives the values with control
	int GetValue(int row, int col) const;

	//takes the row i of a matrix
	BzzVectorInt GetRow(int i) const;
	void GetRow(int i, BzzVectorInt* r);

	//takes the column j of a matrix
	BzzVectorInt GetColumn(int j) const;
	void GetColumn(int i, BzzVectorInt* r);

	// takes the diagonal i of a squared matrix
	// i = 0 principal
	// i > 0 right diagonal
	// i < 0 left diagonal
	BzzVectorInt GetDiagonal(int i) const;

	//	============================================================================
	//	*********************< Modifying access functions >*************************
	//	============================================================================
		// assigns and receives vector values with control
	int& operator () (int row, int col);

	// assigns and receives vector values without control
	int* operator [] (int r)
	{
		return matrix[r];
	}

	// assigns values with control
	void SetValue(int row, int column, int val);

	//substitutes column i with BzzVectorInt rval
	void SetRow(int i, const BzzVectorInt& rval);

	//substitutes column i with xf
	void SetRow(int i, const int& xf);

	// substitutes column j with BzzVectorInt rval
	void SetColumn(int j, const BzzVectorInt& rval);

	// substitutes column j with xf
	void SetColumn(int j, const int& xf);

	// substitutes diagonal j with Vector rval
	void SetDiagonal(int j, BzzVectorInt& rval);

	// substitutes diagonal j with xf
	void SetDiagonal(int j, const int& xf);

	// set all coefficients to xf
	void SetBzzMatrix(const int& xf);

	//	============================================================================
	//	**********************< assignment operators >******************************
	//	============================================================================
	BzzMatrixInt& operator =
		(const BzzMatrixInt& rval);

	BzzMatrixInt& operator =
		(const BzzVectorInt& rval);

	void operator = (int c);
	//	============================================================================
	//	*********************< operators for composing matrices >*******************
	//	============================================================================

			//adds a matrix beneath another
	friend BzzMatrixInt operator &&
		(const BzzMatrixInt& lval, const BzzMatrixInt& rval);

	//add a matrix to the side of another
	friend BzzMatrixInt operator ||
		(const BzzMatrixInt& lval, const BzzMatrixInt& rval);

	//	============================================================================
	//	***************************< test operators >*******************************
	//	============================================================================
	friend char operator ==
		(const BzzMatrixInt& lval, const BzzMatrixInt& rval);

	friend char operator !=
		(const BzzMatrixInt& lval, const BzzMatrixInt& rval);

	//	============================================================================
	//	=======================< Non-modifying functions >==========================
	//	============================================================================

	//	******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);
	void PrintByRows(void);
	void PrintTranspose(void);

	//	********************************< Save >************************************
	void Save(char* filematrix); // formatted
	void Save(char, char* filematrix);// binary

//	*****************************< Max and Min >********************************
	//to have the position Max(im,jc)
	int Max(int* imax = 0, int* jmax = 0);
	int MaxAbs(int* imax = 0, int* jmax = 0);

	int Min(int* imin = 0, int* jmin = 0);
	int MinAbs(int* imin = 0, int* jmin = 0);

	//	***************************< Max and Min >*******************************
		//to have the position Max(im,jc)
	int MaxColumn(int j, int* imax = 0);
	int MaxAbsColumn(int j, int* imax = 0);

	int MinColumn(int j, int* imin = 0);
	int MinAbsColumn(int j, int* imin = 0);

	//	******************************< Max and Min >*******************************
		//to have the position Max(im,jc)
	int MaxRow(int i, int* jmax = 0);
	int MaxAbsRow(int i, int* jmax = 0);

	int MinRow(int i, int* jmin = 0);
	int MinAbsRow(int i, int* jmin = 0);

	void MaxRows(BzzVectorInt* m, BzzVectorInt* im);
	void MaxAbsRows(BzzVectorInt* m, BzzVectorInt* im);
	void MaxColumns(BzzVectorInt* m, BzzVectorInt* im);
	void MaxAbsColumns(BzzVectorInt* m, BzzVectorInt* im);

	void MinRows(BzzVectorInt* m, BzzVectorInt* im);
	void MinAbsRows(BzzVectorInt* m, BzzVectorInt* im);
	void MinColumns(BzzVectorInt* m, BzzVectorInt* im);
	void MinAbsColumns(BzzVectorInt* m, BzzVectorInt* im);

	//	============================================================================
	//	===========================< Modifying Functions >==========================
	//	============================================================================
	void InsertRow(int i, int x);
	void InsertRow(int i, BzzVectorInt& v);
	void InsertColumn(int j, int x);
	void InsertColumn(int j, BzzVectorInt& v);
	void AppendRow(int x);
	void AppendRow(BzzVectorInt& v);
	void AppendRowsFromMatrix(BzzVectorInt& iv, BzzMatrixInt& B);
	void AppendColumn(int x);
	void AppendColumn(BzzVectorInt& v);
	void DeleteRow(int i);
	void DeleteRows(BzzVectorInt& iv);
	void DeleteColumn(int j);
	void DeleteColumns(BzzVectorInt& iv);

	friend void Delete(BzzMatrixInt* A); // eliminate BzzMatrixInt
	void SetToZero(void) // Set all coefficients to zero
	{
		memset(matrix[0], 0, size * sizeof(int));
	}

	friend void ChangeDimensions(int rows, int columns,
		BzzMatrixInt* result, char);

	// recovery of formatted Save
	friend void Load(BzzMatrixInt* A, char* filevector);

	// recovery of binary Save
	friend void Load(BzzMatrixInt* A, char, char* filevector);

	friend void Transpose(BzzMatrixInt* A);
	friend void AntiTranspose(BzzMatrixInt* A);
	friend void HorizontalReflection(BzzMatrixInt* A);
	friend void VerticalReflection(BzzMatrixInt* A);

	friend void SwapMaxMinOrder(BzzMatrixInt* A);
	friend void Swap(BzzMatrixInt* lval, BzzMatrixInt* rval);
	void SwapRows(int i, int j);
	void SwapColumns(int i, int j);
	void MoveRowFromkToj(int k, int j);
	void MoveColumnFromkToj(int k, int j);

	friend void ReorderByRows(BzzMatrixInt* A, BzzVectorInt& iS);
	friend void ReorderByColumns(BzzMatrixInt* A, BzzVectorInt& iS);
	friend void CopyDataFromVector(BzzMatrixInt* lval, BzzVectorInt& rval);
	void UseMatrixAsVector(BzzVectorInt* a);
	void UseMatrixRowAsVector(int i, BzzVectorInt* a);
	int* GetHandle(void) { return matrix[0] + 1; }
	int* GetRowHandle(int i) { return matrix[i] + 1; }
};

// Friend functions with default arguments
void ChangeDimensions(int rows, int columns,
	BzzMatrixInt* result, char zero = 0);


#endif // BZZ_MATRIX_INT_HPP