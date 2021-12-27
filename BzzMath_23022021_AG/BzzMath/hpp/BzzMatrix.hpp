// BZZMATH: Release 7.0

//	=====================< BzzMatrix.hpp >================================
//	* Class BzzMatrix for operations between matrices and vectors			*
//	* in double precision																		*
// * Description:																					*
// *					Dal Fortan al C++	(Capitolo 12)										*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *							and																	*
// *					Scientific C++	(Chapter 12)											*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *																									*
// * Examples: c:\bzzmath\examples\BzzMathBasic\LinearAlgebra\						*
// *				Matrix\Matrix.cpp												*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-1991	Date Written.
//	11-1992	English version.
//	01-1994	Modified Max, Min, MaxAbs, MinAbs.
//	03-1994	Added BzzBaseClass for BzzPrint and BzzMessage.
//	09-1994	Added shadow variable for returning object.
// 09-1994	Added functions ObjectCount, ObjectCountInScope.
// 09-1994	Added functions IProduct, ProductI, ITProduct,ProductIT.
//	11-1994	Conversion to double precision done.
//	02-1995	Added the Bzz prefix to the names of the classes.
//	02-1995	Added ProductT functions.
//	07-1996	The function Product has been improved.
//	04-1997	The function Product has been further enhanced.
//	05-1997	Added BzzPowInt for BzzMatrix.
//	07-1998	Added AppendRow, AppendColumn functions.
//	12-1998	Added Stretch function.
//	12-1998	Added new version for GetRow, GetColumn, GetDiagonale.
//	03-1999	Added ElementByElementProduct.
// 03-1999	Added GetHandle function.
// 05-1999	Added GetRowsNorm2 function.

////////////////// Release 4.0
//	10-2000	Added Balance function.
//	09-2001	Added NormalizeColumn function.
//	09-2001	Added CenterAndNormalizeColumn function.
//	09-2001	Added NormalizeColumns function.
//	09-2001	Added CenterAndNormalizeColumns function.
//	10-2001	Modifyed AppendRow, AppendColumn functions.
//	10-2001	Modifyed InsertRow, InsertColumn functions.
//	04-2003	The function Product has been improved.
//	04-2003	The function TProduct has been improved.
//	04-2003	The function ProductT has been improved.
// 04-2003	Added GetColumnsNorm2 function.

////////////////// Release 5.0
// 07-2003	Added HessenbergByElimination function.
// 07-2003	Added Balance function.
// 10-2003	Added MaxColumn,MaxAbsColumn,MinColumn,MinAbsColumn functions.
// 10-2003	Added MaxRow,MaxAbsRow,MinRow,MinAbsRow functions.
// 10-2003	Added CenterColumn,CenterColumns functions.
// 10-2003	Added GetPointsInMultiGrid function.
// 10-2003	Added GetNewPointMinimumDistanceFromExistingPoints function.
// 10-2003	Added GetMaximumDistance function.
// 11-2003	Added SetToZero function.
//	12-2003	Added NormalizeRow and NormalizeRows functions.
//	07-2004	Added AppendRowsFromMatrix function.
//	11-2004	Added GetRowHandle function.
//	11-2004	Added GetSumElements function.
//	11-2004	Added GetSumAbsElements function.
//	11-2004	Added UseMatrixAsVector function.
//	11-2004	Added RowsElementsProduct and ColumnsElementsProduct function.
//	11-2004	Added ReorderByRows function.
// 12-2004	Added GetRowsSum and GetRowsSumAbs functions.
//	12-2004	Added UseMatrixRowAsVector function.
//	01-2005	Added MaxRows,MaxAbsRows,MaxColumns,MaxAbsColumns functions.
//	01-2005	Added MinRows,MinAbsRows,MinColumns,MinAbsColumns functions.
//	01-2005	Added StretchColumns, StretchRows functions.
//	02-2005	Added CopyDataFromVector function.
//	10-2006	Added ReorderByColumns function.
//	12-2006	Added constructor from file as for MatrixSparse.
//	12-2006	Added GetColumnsSum and GetColumnsSumAbs functions.
//	06-2008	Added SwapMaxMinOrder function.
//	06-2008	Added AntiTranspose function.
//	06-2008	Added HorizontalReflection and VerticalReflection functions.
//	09-2008	Added MoveRowFromkToj function.
//	09-2008	Added MoveColumnromkToj function.

////////////////// Release 6.0
//	12-2008	Added BzzBalance functions.
//	03-2009	Added GetHatMatrix function.
//	03-2009	Added GetSecludedObservations function.
//	04-2009	Added GetClusters function.
//	04-2009	Added GetGoodExperiment function.
// 05-2009	Added GetMatrixHandle function.
//	05-2009	Added TProductT function.
//	05-2009	Added OpenMP for Product, TProduct, ProductT, TProductT functions.
//	09-2009	Added GetExperimentWithDifferentIndependentVariables functions.
// 10-2009	Added ProductForSelectedColumn function.
// 10-2009	Added ProductForSelectedRow function.
// 12-2011	Added BvpNormalizeRows function.

////////////////// Release 7.0
// 05-2013	Added ReorderByRowsAndColumns function.
// 06-2014  Added ProductForSelectedRows function.
// 06-2014  Added ProductForSelectedColumns function.
// 06-2014  Added ProductForSelectedRowsAndColumns function.

//	============================================================================
//	******	BzzMatrix constructors:													*
//	* BzzMatrix A; // default															*
//	* BzzMatrix A = B; // copy-initializer											*
//	* BzzMatrix A(m,n); // sizes and places at 0									*
//	* BzzMatrix A(2,3,1.,2.,3.,4.,5.,6.); // matrix 2X3						*
//	* double x[6]={1.,2.,3.,4.,5.,6.};														*
//	* BzzMatrix A(2,3,x); // from array												*
//	* BzzMatrix A("MAT.DAT"); // Formatted	file									*
//	* BzzMatrix A("MAT.MTX",'*'); // Formatted file as for Sparse			*
//	* BzzMatrix A('*',"MAT.BIN"); // Binary file									*
//	* BzzMatrix A(m,n,B); // submatrix of B										*
//	* BzzMatrix A(m,n,i,j,B);// submatrix of B from i,j						*
//	* BzzMatrix A = v; // from BzzVector									*
//	* BzzMatrix A = D; // from BzzMatrixDiagonal							*
//	* BzzMatrix A = R;#if BZZ_COMPILER == 101
// from BzzMatrixRight								*
//	* BzzMatrix A = L; // from BzzMatrixLeft								*
//	* BzzMatrix A = R; // from BzzMatrixSymmetric						*
//	****************************************************************************
//	***** Access functions:																		*
//	* i = A.Rows(); // numRows																	*
//	* i = A.Columns(); // numColumns															*
//	* who = A.WhoAmI();																			*
//	* xf = A.GetValue(i,j);																		*
//	* v = A.GetRow(i);																			*
//	* A.GetRow(i,&v);																				*
//	* v = A.GetColumn(j);																		*
//	* A.GetColumn(j,&v);																			*
//	* v = A.GetDiagonal(j); // j = 0 principal positive at right					*
//	* A.GetDiagonal(j,&v); // j = 0 principal positive at right						*
//	* xf = A(i,j);																					*
//	* A(i,j) = xf;																					*
//	* xf = A[i][j];																				*
//	* A[i][j] = xf;																				*
//	* A.SetValue(i,j,xf);																		*
//	* A.SetRow(i,xf);																				*
//	* A.SetColumn(j,v);																			*
//	* A.SetColumn(j,xf);																			*
//	* A.SetDiagonal(j,xf); // j = 0 principal positif at right						*
//	* A.SetDiagonal(j,v);  // j = 0 principal positif at right						*
//	* A.SetBzzMatrix(xf);																		*
//	* int count = BzzMatrix::ObjectCount();										*
//	* int countInScope = BzzMatrix::ObjectCountInScope();						*
//	* A.UseMatrixAsVector(&a);																	*
//	* A.UseMatrixRowAsVector(row,&a);														*
//	****************************************************************************
//	***** Assignments:																			*
//	* A = B;	// BzzMatrix																*
//	* A = v;	// BzzeVector																*
//	* A = D;	// BzzMatrixDiagonal														*
//	* A = R;	// BzzMatrixRight															*
//	* A = L;	// BzzMatrixLeft															*
//	* A = S;	// BzzMatrixSymmetric													*
//	****************************************************************************
//	***** Implemented operations :															*
//	* Sum(A,B,&C); // C = A + B;																*
//	* C = A + B;	// C = A + B;																*
//	* Sum(A,B,&A); // A = A + B;																*
//	* Sum(&A,B);	// A = A + B;																*
//	* A += B;		// A = A + B;																*
//	* Sum(A,B,&B); // B = A + B;																*
//	* Sum(A,&B);	// B = A + B;																*
//	* Sum(A,A,&A); // A = A + A;																*
//	* Sum(&A);		// A = A + A;																*
//	* Difference(A,B,&C); // C = A - B;														*
//	* C = A - B;			 // C = A - B;														*
//	* Difference(A,B,&A); // A = A - B;														*
//	* Difference(&A,B);	// A = A - B;														*
//	* A -= B;				 // A = A - B;														*
//	* Difference(A,B,&B); // B = A - B;														*
//	* Difference(A,&B);	// B = A - B;														*
//	* Difference(A,A,&A); // A = A - A;														*
//	* Difference(&A);		// A = A - A;														*
//	* Minus(A,&B);		 // B = -A;																*
//	* B = -A;				// B = -A;															*
//	* Minus(&A);			// A = -A;															*
//	* Product(A,B,&C);	// C = A*B;															*
//	* C = A*B;																						*
//	* Product(A,B,&A);	// A = A*B;															*
//	* Product(&A,B);		// A = A*B;															*
//	* A *= B;				// A = A*B;															*
//	* Product(A,B,&B);	// B = A*B;															*
//	* Product(A,&B);		// B = A*B;															*
//	* Product(A,A,&A);	// A = A*A;															*
//	* Product(&A);		 // A = A*A;															*
//	* Product(A,x,&y);	// y = A*x;															*
//	* y = A*x;																						*
//	* Product(A,x,&x);	// x = A*x;															*
//	* Product(A,&x);		// x = A*x;															*
//	* Product(3.,A,&B);	// B = 3.*A;														*
//	* B = 3.*A;																						*
//	* Product(3.,&A);	 // A = 3.*A;															*
//	* A.ProductForSelectedColumn(j,c,&v);													*
//	* A.ProductForSelectedRow(j,c,&v);														*
//	* A *= 3.;				// A = 3.*A;														*
//	* TProduct(A,B,&C);	// C = ATB;															*
//	* C = A%B;				// C = ATB;															*
//	* TProduct(A,B,&A);	// A = ATB;															*
//	* TProduct(&A,B);	 // A = ATB;															*
//	* A %= B;				// A = ATB;															*
//	* TProduct(A,B,&B);	// B = ATB;															*
//	* TProduct(A,&B);	 // B = ATB;															*
//	* TProduct(A,A,&A);	// A = ATA;															*
//	* TProduct(&A);		// A = ATA;															*
//	* A.SelfTProduct(&S);		// S = ATA;													*
//	* A %= A;				// A = ATA;															*
//	* y = A%x;				// y = ATx;															*
//	* ProductT(x,y,&A);	// A = xyT;															*
//	* A = x->*y;			// A = xyT;															*
//	* ProductT(A,B,&C);	// C = ABT															*
//	* ProductT(A,&B);		// B = ABT															*
//	* ProductT(&A,B);		// A = ABT															*
//	* ProductT(&A);		// A = AAT															*
//	* A.SelfProductT(&S);		// S = AAT;													*
//	* IProduct(A,B,&C);	// C = A-1B															*
//	* IProduct(&A,B);		// A = A-1B															*
//	* IProduct(A,&B);		// B = A-1B															*
//	* ProductI(A,B,&C);	// C = AB-1															*
//	* ProductI(&A,B);		// A = AB-1															*
//	* ProductI(A,&B);		// B = AB-1															*
//	* ITProduct(A,B,&C);	// C = (A-T)B														*
//	* ITProduct(&A,B);	// A = (A-T)B														*
//	* ITProduct(A,&B);	// B = (A-T)B														*
//	* ProductIT(A,B,&C);	// C = A(B-T)														*
//	* ProductIT(&A,B);	// A = A(B-T)														*
//	* ProductIT(A,&B);	// B = A(B-T)														*
//	* IProductT(A,B,&C);	// C = (A-1)BT														*
//	* IProduct(&A,B);		// A = (A-1)BT														*
//	* IProduct(A,&B);		// B = A-1BT														*
//	* TProductI(A,B,&C);	// C = (AT)B-1														*
//	* TProductI(&A,B);	// A = (AT)B-1														*
//	* TProductI(A,&B);	// B = (AT)B-1														*
//	* ElementByElementProduct(A,B,&C);														*
//	* A.RowsElementsProduct(v);																*
//	* A.ColumnsElementsProduct(v);															*
//	* Division(B,3.,&A); // A = B/3.;														*
//	* A = B/3.;																						*
//	* Division(&A,3.);	// A = A/3.;														*
//	* A /= 3.;				// A = A/3.															*
//	* BzzPowInt(&A,5);// A = A^5																*
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
//	* xf = A.NormT();																				*
//	* xf = A.NormR();																				*
//	* xf = A.NormC();																				*
//	* xf = A.NormF();																				*
//	* xf = A.NormI();																				*
//	* xf = A.Norm1();																				*
//	* xf = A.Norm2();																				*
//	* A.Save("MAT.DAT");	// formatted														*
//	* A.Save('*',"MAT.DAT");	// unformatted												*
//	* Load(&A,"MAT.DAT");																		*
//	* Load(&A,'*',"MAT.DAT");																	*
//	* Load(&A,'*','*',"SPAR.DAT");															*
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
// * A.Stretch(m,n);																				*
// * StretchRows(xOld,AOld,xNew,&ANew);													*
// * StretchRows(xOld,AOld,xNew,&ANew,.9);												*
// * StretchColumns(xOld,AOld,xNew,&ANew);												*
// * StretchColumns(xOld,AOld,xNew,&ANew,.9);											*
//	* Swap(&A,&B);																					*
//	* A.SwapRows(i,j);																			*
//	* A.SwapColumns(i,j);																		*
//	* SwapMaxMinOrder(&A);																		*
//	* Transpose(&A);																				*
//	* AntiTranspose(&A);																			*
//	* HorizontalReflection(&A);																*
//	* VerticalReflection(&A);																	*
//	* SumRankOne(u,vT,&A);		// A = A + uvT;											*
//	* SumRankOne(3.,u,vT,&A);	// A = A + 3.*uvT;										*
//	* SumRankOne(u,vT,3.,&A);	// A = A + uvT/3.;										*
// * double *matrix = A.GetHandle();														*
// * double **m = GetMatrixHandle();														*
// * A.GetRowsNorm2(&norm2R);																	*
// * A.GetColumnsNorm2(&norm2C);																*
// * A.GetRowsSum(&sumR);																		*
// * A.GetRowsSumAbs(&sumAbsR);																*
// * A.GetColumnsSum(&sumC);																	*
// * A.GetColumnsSumAbs(&sumAbsC);															*
// * A.Balance(&coeffRows,&coeffColumns);													*
// * norm = A.NormalizeRow(j);																*
// * NormalizeRows(&A);																			*
// * NormalizeRows(&A,&norm);																	*
// * norm = A.NormalizeColumn(j);															*
// * A.CenterAndNormalizeColumn(j,&mean,&norm);											*
// * NormalizeColumns(&A);																		*
// * NormalizeColumns(&A,&norm);																*
// * CenterAndNormalizeColumns(&A);															*
// * GetPointsInMultiGrid(numPointsForEachAxis,xMin,xMax,&X);						*
// * double minimumDistance = GetNewPointMinimumDistanceFromExistingPoints		*
// * 			(newPoint,ExistingPoints);														*
// * double maximumDistance = GetMaximumDistance(A,&i,&j);							*
// * ReorderByRows(&A,iS);																		*
// * ReorderByColumns(&A,iS);																	*
// * ReorderByRowsAndColumns(&A,iRows,iColumns);										*
//	****************************************************************************

#ifndef BZZ_MATRIX_DOUBLE_HPP
#define BZZ_MATRIX_DOUBLE_HPP

// preventive statements
class BzzVectorInt;
class BzzVector;
class BzzMatrix;
class BzzMatrixDiagonal;
class BzzMatrixLeft;
class BzzMatrixRight;
class BzzMatrixSymmetric;

class BzzMatrixSparse;

//	============================================================================
//	========================< class BzzMatrix >===========================
//	============================================================================
class BzzMatrix : public BzzBaseClass
{
	friend class BzzVector;
	friend class BzzMatrixDiagonal;
	friend class BzzMatrixRight;
	friend class BzzMatrixLeft;
	friend class BzzMatrixSymmetric;
	friend class BzzMatrixSparse;
	friend class BzzMatrixBand;
	friend class BzzFactorized;
	friend class BzzFactorizedPLR;
	friend class BzzFactorizedGauss;
	friend class BzzFactorizedQRLQ;
	friend class BzzFactorizedQR;
	friend class BzzFactorizedLQ;
	friend class BzzFactorizedSVD;
	friend class BzzFactorizedSymmetricBase;
	friend class BzzFactorizedBandGauss;
	friend class BzzSave;
	friend class BzzLoad;

private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	double** matrix;
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

	// for assigning and initialising from BzzMatrixDiagonal
	void CopyDiagonal(const BzzMatrixDiagonal& rval);

	// for assigning and initialising from BzzMatrixRight
	void CopyRight(const BzzMatrixRight& rval);

	// for assigning and initialising from BzzMatrixLeft
	void CopyLeft(const BzzMatrixLeft& rval);

	// for assigning and initialising from BzzMatrixSymmetric
	void CopySymm(const BzzMatrixSymmetric& rval);

	// private constructor BzzMatrix A type('*',3,5);
	BzzMatrix(char, int mrows, int columns);

	// for Jacobian
	BzzVector xJacobian, fJacobian;

public:
	//	============================================================================
	//	**************************< constructors >**********************************
	//	============================================================================
		// default // BzzMatrix A;
	BzzMatrix(void);

	// copy-initializer // BzzMatrix A = B;

#if BZZ_COMPILER == 101
	BzzMatrix(const BzzMatrix& rval);
#else
	BzzMatrix(BzzMatrix& rval);
#endif

	// sizes and initialises at 0 //BzzMatrix A(3,5);
	BzzMatrix(int rows, int columns);

	// sizes and initialises
	// BzzMatrix A(2,3,1.,2.,3.,4.,5.,6.);
	BzzMatrix(int rows, int columns, double a11, ...);

	// initialises from array
	// double w[6]={1.,2.,3.,4.,5.6}; BzzMatrix A(3,5,w);
	BzzMatrix(int rows, int columns, double* initvalues);

	// initialises from formatted File
	// BzzMatrix A("MAT.DAT");
	BzzMatrix(char* filematrix);

	// initialises from formatted File  as for Sparse
	// BzzMatrix A("MAT.MTX",'*');
	BzzMatrix(char* filematrix, char);

	// initialises from binary File
	// BzzMatrix A('*',"MAT.BIN"); See Save
	BzzMatrix(char, char* filematrix);

	// makes a submatrix of rows,columns
	BzzMatrix(int rows, int columns, const BzzMatrix& rval);

	// likewise starting from irow,jcol
	BzzMatrix(int rows, int columns, int irow, int jcol,
		const BzzMatrix& rval);

	// from BzzVector
	BzzMatrix(const BzzVector& rval);

	// from BzzMatrixDiagonal
	BzzMatrix(const BzzMatrixDiagonal& rval);

	// from BzzMatrixRight
	BzzMatrix(const BzzMatrixRight& rval);

	// from BzzMatrixLeft
	BzzMatrix(const BzzMatrixLeft& rval);

	// from BzzMatrixSymmetric
	BzzMatrix(const BzzMatrixSymmetric& rval);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzMatrix(void);

	//	============================================================================
	//	**********************< Non-modifying access functions >********************
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
	double GetValue(int row, int col) const;

	//takes the row i of a matrix
	BzzVector GetRow(int i) const;
	void GetRow(int i, BzzVector* v);

	//takes the column j of a matrix
	BzzVector GetColumn(int j) const;
	void GetColumn(int j, BzzVector* v);

	// takes the diagonal i of a matrix
	// i = 0 principal
	// i > 0 right diagonal
	// i < 0 left diagonal
	BzzVector GetDiagonal(int i) const;
	void GetDiagonal(int i, BzzVector* v);

	//	============================================================================
	//	*********************< Modifying access functions >*************************
	//	============================================================================
		// assigns and receives vector values with control
	double& operator () (int row, int col);

	// assigns and receives vector values without control
	double* operator [] (int r)
	{
		return matrix[r];
	}

	// assigns values with control
	void SetValue(int row, int column, double val);

	//substitutes column i with rval
	void SetRow(int i, const BzzVector& rval);

	//substitutes column i with xf
	void SetRow(int i, const double& xf);

	// substitutes column j with rval
	void SetColumn(int j, const BzzVector& rval);

	// substitutes column j with xf
	void SetColumn(int j, const double& xf);

	// substitutes diagonal j with Vector rval
	void SetDiagonal(int j, BzzVector& rval);

	// substitutes diagolnal j with xf
	void SetDiagonal(int j, double xf);

	// set all coefficients to xf
	void SetBzzMatrix(const double& xf);

	//	============================================================================
	//	**********************< assignment operators >******************************
	//	============================================================================
	BzzMatrix& operator =
		(const BzzMatrix& rval);

	BzzMatrix& operator =
		(const BzzVector& rval);

	BzzMatrix& operator =
		(const BzzMatrixDiagonal& rval);

	BzzMatrix& operator =
		(const BzzMatrixLeft& rval);

	BzzMatrix& operator =
		(const BzzMatrixRight& rval);

	BzzMatrix& operator =
		(const BzzMatrixSymmetric& rval);

	void operator = (double c);

	// transforms a BzzMatrix in Factorized
	// BzzMatrix gets destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrix* lval, BzzFactorized* rval);
	void UseMatrixAsVector(BzzVector* a);
	void UseMatrixRowAsVector(int i, BzzVector* a);

	//	============================================================================
	//	*********************< operators for composing matrices >*******************
	//	============================================================================

			//adds a matrix beneath another
	friend BzzMatrix operator &&
		(const BzzMatrix& lval, const BzzMatrix& rval);

	//add a matrix to the side of another
	friend BzzMatrix operator ||
		(const BzzMatrix& lval, const BzzMatrix& rval);

	//	============================================================================
	//	***************************<  test operators >******************************
	//	============================================================================
	friend char operator ==
		(const BzzMatrix& lval, const BzzMatrix& rval);

	friend char operator !=
		(const BzzMatrix& lval, const BzzMatrix& rval);

	//	============================================================================
	//	===============================< OPERATIONS >===============================
	//	============================================================================

	//	============================================================================
	//	********************************< Sum >*************************************
	//	============================================================================
		 // Sum(A,B,&C); C = A + B;
	friend void Sum
	(const BzzMatrix& lval, const BzzMatrix& rval,
		BzzMatrix* result);

	// C = A + B;
	friend BzzMatrix operator +
		(const BzzMatrix& lval, const BzzMatrix& rval);

	// Sum(&A,B); A = A + B;
	friend void Sum
	(BzzMatrix* lvalAndResult, const BzzMatrix& rval);

	// A += B; A = A + B;
	BzzMatrix& operator +=
		(const BzzMatrix& rval);

	// Sum(B,&A); A = B + A;
	friend void Sum
	(const BzzMatrix& lval, BzzMatrix* rvalAndResult);

	// Sum(&A); A = A + A;
	friend void Sum
	(BzzMatrix* lvalRvalAndResult);

	//	============================================================================
	//	*****************************< Difference >*********************************
	//	============================================================================
		// Difference(A,B,&C); C = A - B;
	friend void Difference
	(const BzzMatrix& lval, const BzzMatrix& rval,
		BzzMatrix* result);

	// C = A - B;
	friend BzzMatrix operator -
		(const BzzMatrix& lval, const BzzMatrix& rval);

	// Difference(&A,B); A = A - B;
	friend void Difference
	(BzzMatrix* lvalAndResult, const BzzMatrix& rval);

	// A -= B; A = A - B;
	BzzMatrix& operator -=
		(const BzzMatrix& rval);

	// Difference(B,&A); A = B - A;
	friend void Difference
	(const BzzMatrix& lval, BzzMatrix* rvalAndResult);

	// Difference(&A); A = A - A;
	friend void Difference
	(BzzMatrix* lvalRvalAndResult);

	//	============================================================================
	//	*******************************< Minus >************************************
	//	============================================================================
	friend void Minus
	(const BzzMatrix& rval, BzzMatrix* result);

	friend BzzMatrix operator -
		(const BzzMatrix& rval);

	friend void Minus
	(BzzMatrix* rvalAndResult);

	//	============================================================================
	//	******************************< Product >***********************************
	//	============================================================================
	friend void Product // C=A*B;
	(BzzMatrix& lval, BzzMatrix& rval,
		BzzMatrix* result);

	friend BzzMatrix operator *		// C=A*B;
		(BzzMatrix& lval, BzzMatrix& rval);

	friend void Product		 // A = A * B;
	(BzzMatrix* lvalAndResult, BzzMatrix& rval);

	BzzMatrix& operator *=			// A = A * B;
		(BzzMatrix& rval);

	friend void Product		 // B = A * B;
	(BzzMatrix& lval, BzzMatrix* rvalAndResult);

	friend void Product		 // A = A * A;
	(BzzMatrix* lvalRvalAndResult);

	friend void Product		// y =A*x;
	(const BzzMatrix& lval, const BzzVector& rval,
		BzzVector* result);

	friend BzzVector operator *	 // y = A*x;
		(const BzzMatrix& lval, const BzzVector& rval);

	friend void Product		 // A = 3.*B;
	(double lval, const BzzMatrix& rval, BzzMatrix* result);

	friend BzzMatrix operator *		// A = 3.*B;
		(double lval, const BzzMatrix& rval);

	friend void Product		 // A = 3.* A;
	(double lval, BzzMatrix* rvalAndResult);

	BzzMatrix& operator *=	 // A = 3.*A;
		(double rval);

	friend void Product(BzzMatrix& lval, BzzMatrixBand& rval,
		BzzMatrix* result);

	friend void Product(BzzMatrixBand& lval, BzzMatrix& rval,
		BzzMatrix* result);
	void ProductForSelectedColumn(int j, double c, BzzVector* v);
	void ProductForSelectedRow(int j, double c, BzzVector* v);

	friend void ProductForSelectedRows(BzzMatrix& A, BzzVector& x,
		BzzVectorInt& rows, BzzVector* v);
	friend void ProductForSelectedColumns(BzzMatrix& A, BzzVector& x,
		BzzVectorInt& col, BzzVector* v);
	friend void ProductForSelectedRowsAndColumns(BzzMatrix& A, BzzVector& x,
		BzzVectorInt& rows, BzzVectorInt& col, BzzVector* v);

	//	============================================================================
	//	*****************************< TProduct >***********************************
	//	============================================================================
	friend void TProduct // C=ATB;
	(BzzMatrix& lval, BzzMatrix& rval,
		BzzMatrix* result);

	friend BzzMatrix operator %
		(BzzMatrix& lval, BzzMatrix& rval); // ATB

		// TProduct(A,x,&y); y =ATx; y =A%x;
	friend void TProduct(BzzMatrix& lval,
		BzzVector& rval, BzzVector* result);

	friend BzzVector operator %
		(BzzMatrix& lval, BzzVector& rval);

	friend void TProduct		 // A = ATB;
	(BzzMatrix* lvalAndResult, BzzMatrix& rval);

	BzzMatrix& operator %=
		(BzzMatrix& rval);

	friend void TProduct		 // B = ATB;
	(BzzMatrix& lval, BzzMatrix* rvalAndResult);

	friend void TProduct		 // A = A % A;
	(BzzMatrix* lvalRvalAndResult);

	void SelfTProduct		 // S = A % A;
	(BzzMatrixSymmetric* S);

	// TProduct(A,&x); x = ATx;
	friend void TProduct
	(BzzMatrix& lval, BzzVector* rvalAndresult);

	//	============================================================================
	//	*****************************< ProductT >***********************************
	//	============================================================================
		// ProductT(x,y,&A); A = xyT; A = x->*y;
	friend void ProductT(const BzzVector& lval,
		const BzzVector& rval, BzzMatrix* result);

	// A = x->*y; A = xyT;
	friend BzzMatrix operator ->*
		(const BzzVector& lval, const BzzVector& rval);

	friend void ProductT		 // C = ABT ;
	(BzzMatrix& lval, BzzMatrix& rval,
		BzzMatrix* result);

	friend void ProductT		 // A = ABT ;
	(BzzMatrix* lvalAndResult, BzzMatrix& rval);

	friend void ProductT		 // B = ABT ;
	(BzzMatrix& lval, BzzMatrix* rvalAndResult);

	friend void ProductT		 // A = AAT ;
	(BzzMatrix* lvalRvalAndResult);

	void SelfProductT		 // S = AAT;
	(BzzMatrixSymmetric* S);

	//	============================================================================
//	*****************************< TProductT >***********************************
//	============================================================================
	friend void TProductT // C=ATBT;
	(BzzMatrix& lval, BzzMatrix& rval,
		BzzMatrix* result);

	friend void IProduct		// C = A-1B
	(BzzMatrix& A, BzzMatrix& B, BzzMatrix* C);

	friend void IProduct		// A = A-1B
	(BzzMatrix* A, BzzMatrix& B);

	friend void IProduct		// B = A-1B
	(BzzMatrix& A, BzzMatrix* B);

	friend void ProductI		// C = AB-1
	(BzzMatrix& A, BzzMatrix& B, BzzMatrix* C);

	friend void ProductI		// A = AB-1
	(BzzMatrix* A, BzzMatrix& B);

	friend void ProductI		// B = AB-1
	(BzzMatrix& A, BzzMatrix* B);

	friend void ITProduct		// C = (A-T)B
	(BzzMatrix& A, BzzMatrix& B, BzzMatrix* C);

	friend void ITProduct		// A = (A-T)B
	(BzzMatrix* A, BzzMatrix& B);

	friend void ITProduct		// B = (A-T)B
	(BzzMatrix& A, BzzMatrix* B);

	friend void ProductIT		// C = A(B-T)
	(BzzMatrix& A, BzzMatrix& B, BzzMatrix* C);

	friend void ProductIT		// A = A(B-T)
	(BzzMatrix* A, BzzMatrix& B);

	friend void ProductIT		// B = A(B-T)
	(BzzMatrix& A, BzzMatrix* B);

	friend void IProductT		// C = (A-1)BT
	(BzzMatrix& A, BzzMatrix& B, BzzMatrix* C);

	friend void IProductT		// A = (A-1)BT
	(BzzMatrix* A, BzzMatrix& B);

	friend void IProductT		// B = A-1BT
	(BzzMatrix& A, BzzMatrix* B);

	friend void TProductI		// C = (AT)B-1
	(BzzMatrix& A, BzzMatrix& B, BzzMatrix* C);

	friend void TProductI		// A = (AT)B-1
	(BzzMatrix* A, BzzMatrix& B);

	friend void TProductI		// B = (AT)B-1
	(BzzMatrix& A, BzzMatrix* B);
	friend void ElementByElementProduct(BzzMatrix& A, BzzMatrix& B,
		BzzMatrix* C);
	void RowsElementsProduct(BzzVector& v);
	void ColumnsElementsProduct(BzzVector& v);

	//	============================================================================
	//	******************************< Division >**********************************
	//	============================================================================
	friend void Division	// A =B/3.;
	(const BzzMatrix& lval, double rval, BzzMatrix* result);

	friend BzzMatrix operator /		// A= B/3.;
		(const BzzMatrix& lval, double rval);

	friend void Division			 // A = A/3.;
	(BzzMatrix* lvalAndResult, double rval);

	BzzMatrix& operator /=				// A /= 3.;
		(double rval);

	//	============================================================================
	//	=========================< Non-modifying functions >========================
	//	============================================================================

	//	*****************************< BzzPrint >***********************************
	virtual void ObjectBzzPrint(void);

	//	********************************< Save >************************************
	void Save(char* filematrix); // formatted
	void Save(char, char* filematrix);// binary

//	******************************< Max and Min >*******************************
	//to have the position Max(im,jc)
	double Max(int* imax = 0, int* jmax = 0);
	double MaxAbs(int* imax = 0, int* jmax = 0);

	double Min(int* imin = 0, int* jmin = 0);
	double MinAbs(int* imin = 0, int* jmin = 0);

	//	***************************< Max and Min >*******************************
		//to have the position Max(im,jc)
	double MaxColumn(int j, int* imax = 0);
	double MaxAbsColumn(int j, int* imax = 0);

	double MinColumn(int j, int* imin = 0);
	double MinAbsColumn(int j, int* imin = 0);

	//	******************************< Max and Min >*******************************
		//to have the position Max(im,jc)
	double MaxRow(int i, int* jmax = 0);
	double MaxAbsRow(int i, int* jmax = 0);

	double MinRow(int i, int* jmin = 0);
	double MinAbsRow(int i, int* jmin = 0);

	void MaxRows(BzzVector* m, BzzVectorInt* im);
	void MaxAbsRows(BzzVector* m, BzzVectorInt* im);
	void MaxColumns(BzzVector* m, BzzVectorInt* im);
	void MaxAbsColumns(BzzVector* m, BzzVectorInt* im);

	void MinRows(BzzVector* m, BzzVectorInt* im);
	void MinAbsRows(BzzVector* m, BzzVectorInt* im);
	void MinColumns(BzzVector* m, BzzVectorInt* im);
	void MinAbsColumns(BzzVector* m, BzzVectorInt* im);

	//	********************************< Norms >***********************************
	double NormT(void);
	double NormR(void);
	double NormC(void);
	double NormF(void);
	double NormI(void) { return NormR(); }
	double Norm1(void) { return NormC(); }
	double Norm2(void);
	double GetSumElements(void);
	double GetSumAbsElements(void);

	//	============================================================================
	//	========================< Modifying Functions >=============================
	//	============================================================================
	void InsertRow(int i, double x);
	void InsertRow(int i, BzzVector& v);
	void InsertColumn(int j, double x);
	void InsertColumn(int j, BzzVector& v);
	void AppendRow(double x);
	void AppendRow(BzzVector& v);
	void AppendRowsFromMatrix(BzzVectorInt& iv, BzzMatrix& B);
	void AppendColumn(double x);
	void AppendColumn(BzzVector& v);
	void DeleteRow(int i);
	void DeleteRows(BzzVectorInt& iv);
	void DeleteColumn(int j);
	void DeleteColumns(BzzVectorInt& iv);
	friend void Delete(BzzMatrix* A); // eliminate BzzMatrix
	void SetToZero(void) // Set all coefficients to zero
	{
		memset(matrix[0], 0, size * sizeof(double));
	}
	friend void ChangeDimensions(int rows, int columns,
		BzzMatrix* result, char);
	void Stretch(int rows, int columns);
	friend void StretchRows(BzzVector& xOld, BzzMatrix& AOld,
		BzzVector& xNew, BzzMatrix* ANew, double);
	friend void StretchColumns(BzzVector& xOld, BzzMatrix& AOld,
		BzzVector& xNew, BzzMatrix* ANew, double);

	// recovery of formatted Save
	friend void Load(BzzMatrix* A, char* filevector);

	// recovery of binary Save
	friend void Load(BzzMatrix* A, char, char* filevector);

	// recovery of formatted Save from sparse matrix
	friend void Load(BzzMatrix* A, char, char, char* filematrix);

	friend void Transpose(BzzMatrix* A);
	friend void AntiTranspose(BzzMatrix* A);
	friend void HorizontalReflection(BzzMatrix* A);
	friend void VerticalReflection(BzzMatrix* A);

	friend char Inverse(BzzMatrix* A);	// Gauss Jordan

	friend void SumRankOne(const BzzVector& u,
		const BzzVector& vT, BzzMatrix* result);
	friend void SumRankOne(double product, const BzzVector& u,
		const BzzVector& vT, BzzMatrix* result);
	friend void SumRankOne(const BzzVector& u,
		const BzzVector& vT, double divisor, BzzMatrix* result);

	friend void SwapMaxMinOrder(BzzMatrix* A);
	friend void Swap(BzzMatrix* lval, BzzMatrix* rval);
	void SwapRows(int i, int j);
	void SwapColumns(int i, int j);
	void MoveRowFromkToj(int k, int j);
	void MoveColumnFromkToj(int k, int j);

	friend void CopyDataFromVector(BzzMatrix* lval, BzzVector& rval);

	friend void Solve(const BzzMatrixRight& R, BzzMatrix* BX);
	friend void Solve(const BzzMatrixLeft& L, BzzMatrix* BX);
	friend void BzzPowInt(BzzMatrix* A, int n);
	double** GetMatrixHandle(void) { return matrix; }
	double* GetHandle(void) { return matrix[0] + 1; }
	double* GetRowHandle(int i) { return matrix[i] + 1; }
	void GetRowsNorm2(BzzVector* norm2R);
	void GetColumnsNorm2(BzzVector* norm2C);
	void GetRowsSum(BzzVector* sumR);
	void GetRowsSumAbs(BzzVector* sumAbsR);
	void GetColumnsSum(BzzVector* sumC);
	void GetColumnsSumAbs(BzzVector* sumAbsC);
	void Balance(BzzVector* coeffRows, BzzVector* coeffColumns);
	void BzzBalance(BzzVector* coeffRows, BzzVector* coeffColumns);
	void BzzBalance(BzzVector& b, BzzVector* coeffRows, BzzVector* coeffColumns);
	friend void HessenbergByElimination(BzzMatrix* A);
	double NormalizeRow(int j);
	friend void NormalizeRows(BzzMatrix* A);
	friend void NormalizeRows(BzzMatrix* A, BzzVector* norm);
	double NormalizeColumn(int j);
	double CenterColumn(int j);
	void CenterAndNormalizeColumn(int j, double* mean, double* norm);
	friend void NormalizeColumns(BzzMatrix* A);
	friend void NormalizeColumns(BzzMatrix* A, BzzVector* norm);
	friend void CenterColumns(BzzMatrix* A);
	friend void CenterColumns(BzzMatrix* A, BzzVector* mean);
	friend void CenterAndNormalizeColumns(BzzMatrix* A);
	friend void CenterAndNormalizeColumns(BzzMatrix* A,
		BzzVector* mean, BzzVector* norm);
	void GetSecludedObservations(BzzVectorInt* ix = 0, BzzMatrix* X = 0,
		BzzVectorInt* ip = 0, BzzMatrix* P = 0, BzzVectorInt* ir = 0, BzzMatrix* R = 0);
	void GetClusters(BzzMatrixInt* I, BzzMatrix* C);
	friend void GetPointsInMultiGrid(BzzVectorInt& numPointsForEachAxis,
		BzzVector& xMin, BzzVector& xMax, BzzMatrix* XX);
	friend void GetNewPointMinimumDistanceFromExistingPoints(BzzVector& newPoint,
		BzzMatrix& ExistingPoints, double* minimumDistance, int* existing);
	friend void GetNewPointsMinimumDistanceFromExistingPoints(BzzMatrix& newPoints,
		BzzMatrix& existingPoints, BzzVector* minimumDistance,
		BzzVectorInt* newP, BzzVectorInt* existingP);
	friend double GetMaximumDistance(BzzMatrix& A, int* i, int* j);
	void GetGoodExperiment(BzzVector& xMin, BzzVector& xMax, BzzVector* x);
	void GetExperimentWithDifferentIndependentVariables(BzzVector& xMin, BzzVector& xMax, BzzVector* x);
	friend void ReorderByRows(BzzMatrix* A, BzzVectorInt& iS);
	friend void ReorderByColumns(BzzMatrix* A, BzzVectorInt& iS);
	friend void ReorderByRowsAndColumns(BzzMatrix* A, BzzVectorInt& iRows, BzzVectorInt& iColumns);
	void GetHatMatrix(BzzMatrix* H);
	void BvpNormalizeRows(void);
};

// ============================================================================
// ======================< class BzzOrthogonalization >========================
// ============================================================================

class BzzOrthogonalization
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

public:
	int BzzModifiedGramSchmidtRows(BzzMatrix* A);
	int BzzHouseholderColumns(BzzMatrix* A);
	int BzzHouseholderRows(BzzMatrix* A);
	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor BzzOrthogonalization o;
	BzzOrthogonalization(void) {};

	// ============================================================================
	// *****************************< destructor >	*******************************
	// ============================================================================
	~BzzOrthogonalization(void) {};
};


void ChangeDimensions(int rows, int columns,
	BzzMatrix* result, char zero = 0);
void Stretch(int rows, int columns);
void StretchRows(BzzVector& xOld, BzzMatrix& AOld,
	BzzVector& xNew, BzzMatrix* ANew, double smooth = .5);
void StretchColumns(BzzVector& xOld, BzzMatrix& AOld,
	BzzVector& xNew, BzzMatrix* ANew, double smooth = .5);


#endif // BZZ_MATRIX_DOUBLE_HPP