// TODO: OrderingLinearProgramming

// BZZMATH: Release 7.0

//	=======================< BzzMatrixSparse.hpp >========================
//	* Class BzzMatrixSparse for operations between sparse matrices and	*
// * vectors in double precision																*
// * Description:																					*
// *					Dal Fortan al C++	(Capitolo 13)										*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *							and																	*
// *					Scientific C++	(Chapter 13)											*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *																									*
// * Examples: BzzMath\examples\BzzMathBasic\LinearAlgebra\							*
// *				MatrixSparse\MatrixSparse.cpp								*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	04-1991	Date Written.
//	11-1992	English version.
//	03-1994	Added BzzBaseClass for BzzPrint and BzzMessage.
// 09-1994	Added functions Existing, CountElements, ObjectCount,
//				ObjectCountInScope.
//	11-1994	Conversion to double precision done.
//	02-1994	Added the Bzz prefix to the names of the classes.
//	04-1997	Added OrderingGauss.
//	05-1997	Added GetRow, GetColumn.
//	05-1997	Added BzzPrintStructure.
//	07-1997	Added TransposeSolveLeft, TransposeSolveRight.
//	10-1997	Added OrderingLinearProgramming.
//	10-1997	Added AddOneColumn, AddOneRow, AddNColumns, AddNRows.
//	10-1997	Added DeleteRow, DeleteColumn.
// 05-1999	Added GetRowsNorm2 function.
// 05-1999	Added GetRowsNorm2X function.
//	05-1999	Added DeleteLastNRows.
//	07-1999	Added Append, Merge, SplitByColumns, SplitByRows.
//	10-1999	Added GetRowsExistence.

/////////// Release 4.0
//	01-2000	Added BeginScanning.
//	04-2000	Added ReleaseElement,ReleaseAllElementsInRow,
//				ReleaseAllElementsInColumn,ReleaseAllElementsInMatrix.
//	04-2000	Added InsertStoredElement.
//	04-2000	Added class BzzMatrixSparseStore.
//	06-2000	Added GetNumEquationsForEachVariable.
//	06-2000	Added GetNumVariablesForEachEquation.
// 09-2000	Added GetStartingElementInRow.
// 09-2000	Added OrderingLQ.
// 10-2000	Added MoveRowsFromFirstMatrixAndAppendToSecond.
//	11-2000	Added GetNumVariablesForEachEquationAndNumEquationsForEachVariable.
//	11-2000	Added GetEquationsContainingSingletons.
//	11-2000	Added GetLinearlyDependentCouplesOfEquations.
//	11-2000	Added GetLinearlyDependentCouplesOfEquations,
//				GetPositiveLinearlyDependentCouplesOfEquations,
//				GetNegativeLinearlyDependentCouplesOfEquations.
//	07-2001	Added SetDiagonal.
//	04-2002	Added new version of OrderingGauss.
//	04-2002	Added GetRowsNorm1.

/////////// Release 5.0
//	08-2003	Added assignement from BzzMatrix.
//	09-2003	Added CountElementsOutsideLowerBand.
//	09-2003	Added CountElementsOutsideUpperBand.
//	09-2003	Added CountElementsOutsideBands.
//	10-2004	Added OrderingGaussWithoutPivotimg.
//	10-2004	Added GaussFillingUp

/////////// Release 6.0
//	03-2005	Added SwapRows and SwapColumns functions.
//	03-2005	Added SwapElements function.
//	03-2005	Added Transpose function.
//	12-2006	Added GetRowsSum.
//	12-2006	Added GetColumnsSum.

//	============================================================================
//	******	BzzMatrixSparse constructors:											*
//	* BzzMatrixSparse A; // default													*
//	* BzzMatrixSparse A = B; // copy-initializer									*
//	* BzzMatrixSparse A(m,n); // sizes												*
//	* BzzMatrixSparse A("SPAR.DAT"); // Formatted	file						*
//	****************************************************************************
//	***** Access functions:																		*
//	* i = A.Rows(); // numRows																	*
//	* i = A.Columns(); // numColumns															*
//	* who = A.WhoAmI();																			*
//	* xf = A(i,j);																					*
//	* xf = A.GetValue(i,j);																		*
//	* r = A.GetRow(i);																			*
//	* c = A.GetColumn(i);																		*
//	* A(i,j) = xf;																					*
//	* A.SetValue(i,j,xf);																		*
//	* char t = A.Existing(i,j);																*
//	* double *ptrVal = A.Scanning(&i,&j,&val);											*
//	* A.BeginScanning();																			*
//	* int countElements = A.CountElements();												*
//	* int countElementsOutsideLowerBand = A.CountElementsOutsideLowerBand(10);	*
//	* int countElementsOutsideUpperrBand =													*
//	*			A.CountElementsOutsideUpperrBand(15);										*
//	* A.CountElementsOutsideBands(10,15,&low,&up);										*
// * A.BzzPrintExistingElements();															*
// * A.BzzPrintStructure();																	*
//	* int count = BzzMatrixSparse::ObjectCount();								*
//	* int countInScope = BzzMatrixSparse::ObjectCountInScope();				*
//	****************************************************************************
//	****************************************************************************
//	***** Assignments:																			*
//	* A = B;	// BzzMatrixSparse														*
//	* A = B;	// B BzzMatrix																*
//	****************************************************************************
//	***** Implemented operations :															*
//	* C = A + B;	// C = A + B																*
//	* A += B;		// A = A + B																*
//	* C = A - B;	// C = A - B																*
//	* A -= B;		// A = A - B																*
//	* Product(A,B,&C);	// C = AB															*
//	* Product(A,x,&y);	// y = Ax															*
//	* y = A*x;;		// y = Ax																	*
// * TProduct(A,x,&y);		// y = ATx;														*
// * TProduct(A,&x);			// x = ATx;														*
//	* y = A%x;;		// y = ATx																	*
//	* Product(3.,&A);	// y = Ax																*
//	****************************************************************************
//	***** BzzMessage and BzzPrint																*
//	* A.BzzPrint("Comment"); // always														*
//	* A.BzzMessage("Comment"); // only if bzzMessageYesNo == 'Y'					*
//	****************************************************************************
//	***** Save and Load																			*
//	* A.Save("spar.dat");		// formatted												*
//	* A.Save('*',"spar.bin");	// non formatted											*
//	* Load(&A,"spar.dat");		// formatted												*
//	* Load(&A,'*',"spar.bin");	// non formatted											*
//	****************************************************************************
//	***** Remove, Delete, Clean, ChangeDimensions and Swap							*
//	* A.RemoveElement(3,151);																	*
//	* A.RemoveAllElementsInRow(row);															*
//	* A.RemoveAllElementsInColumn(column);													*
//	* A.RemoveAllElementsInMatrix();		// eliminates all elements					*
//	* A.CleanMatrix(eps);	// eliminates those <= eps									*
//	* A.DeleteRow(5);																				*
//	* A.DeleteLastNRows(3);																		*
//	* A.DeleteColumn(9);																			*
//	* Delete(&A);																					*
//	* Transpose(&A);																				*
//	* ChangeDimensions(newr,newc,&A);														*
//	* A.AddOneColumn();																			*
//	* A.AddNColumns();																			*
//	* A.AddOneRow();																				*
//	* A.AddNRows();																				*
//	* Swap(&A,&B);																					*
//	* A.SwapRows(123,435);																		*
//	* A.SwapColumns(215,7);																		*
//	* A.SwapElements(215,7,415,212);															*
//	****************************************************************************
//	***** SetDiagonal																				*
//	* A.SetDiagonal(0,5.);	// Principal													*
//	* A.SetDiagonal(-5,2.);	// Lower															*
//	* A.SetDiagonal(8,-3.);	// Upper															*
//	****************************************************************************
//	***** Norms																						*
//	* double normi = A.NormI();																*
//	****************************************************************************
//	***** Functions for tiangular linear systems											*
//	* A.SolveRight(&bx);																			*
//	* A.SolveLeft(&bx);																			*
//	* A.TransposeSolveRight(&bx);																*
//	* A.TransposeSolveLeft(&bx);																*
//	****************************************************************************

#ifndef BZZ_SPARSE_DOUBLE_HPP
#define BZZ_SPARSE_DOUBLE_HPP

//	============================================================================
//	========================< class BzzMatrixSparse >=====================
//	============================================================================

class BzzFactorizedSparseGauss;

struct ElementBzzMatrixSparse
{
	int column;
	double value;
	ElementBzzMatrixSparse* next;
};

class BzzMatrixSparseStore
{
public:
	ElementBzzMatrixSparse* elementStore;
	BzzMatrixSparseStore(void);
	~BzzMatrixSparseStore(void);
	friend void Delete(BzzMatrixSparseStore* store);
	int CountElements(void);
};

class BzzMatrixSparse : public BzzBaseClass
{
	friend class BzzFactorizedSparseGauss;
	friend class BzzLinearProgramming;
	friend class BzzFactorizedSparseLQ;
	friend class BzzFactorizedSparseLockedLQ;
	friend class BzzMatrixCoefficientsExistence;
	friend class BzzFactorizedBandGauss;
	friend class BzzMatrixBand;
	friend class BzzMatrixSparseLockedByRows;
	friend class BzzMatrixSparseLockedByColumns;
	friend class BzzFactorizedSparseUnspecified;
	friend class BzzSave;
	friend class BzzLoad;
	friend class BzzLinearProgrammingAttic;

	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixSparse* lval, BzzFactorizedSparseGauss* rval);
	friend void ReplaceBzzMatrixSparseWithBzzMatrixSparseLocked
	(BzzMatrixSparse* S, BzzMatrixSparseLockedByRows* P);
	friend void ReplaceBzzMatrixSparseWithMinusBzzMatrixSparseLocked
	(BzzMatrixSparse* S, BzzMatrixSparseLockedByRows* P);
	friend void ReplaceBzzMatrixSparseWithBzzMatrixSparseLocked
	(BzzMatrixSparse* S, BzzMatrixSparseLockedByColumns* P);
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixSparse* lval, BzzFactorizedSparseLQ* rval);
	friend void AnalyzeCompatibilities(BzzMatrixSparse& E, BzzVector& e,
		BzzVector& L, BzzVector& U, BzzVectorInt* ie);
	friend void AnalyzeEqualityDegeneration(BzzMatrixSparse& E, BzzVector& e,
		BzzVector& L, BzzVector& U, BzzVector& x,
		BzzVectorInt* ie, BzzVectorInt* ive);
	friend void AnalyzeInequalityDegeneration(BzzMatrixSparse& E, BzzVector& e,
		BzzVector& L, BzzVector& U, BzzVector& x,
		BzzVectorInt* ie, BzzVectorInt* ive);
	friend char FindVariablesForDegeneration(BzzMatrixSparse& E, BzzVector& e,
		BzzMatrixSparse& D, BzzVector& d,
		BzzVector& L, BzzVector& U, BzzVector& x,
		BzzVectorInt* iv);

private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	int numRows, numColumns;
	int whoAmI;
	ElementBzzMatrixSparse** elementRow;
	ElementBzzMatrixSparse* elemScanning;

	int rowScanning;

	int	lowerBand,
		upperBand;

	double& InsertElement
	(ElementBzzMatrixSparse* elem, int row,
		int column, int first);
	void Initialize(int rows, int columns);
	void Copy(const BzzMatrixSparse& rval);

public:
	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default constructor
	BzzMatrixSparse(void);

	// copy constructor
	BzzMatrixSparse(const BzzMatrixSparse& rval);

	// sizing constructor
	BzzMatrixSparse(int rows, int columns);
	// initialises from formatted File

	// BzzMatrixSparse A("SPAR.DAT");
	BzzMatrixSparse(char* filematrix);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzMatrixSparse(void);

	//	============================================================================
	//	****************************< Access functions >****************************
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

	// assigns and receives vector values with control
	double& operator ()
		(int row, int column);

	void SetValue(int i, int j, double val);
	double GetValue(int i, int j);
	BzzVector GetRow(int i);
	BzzVector GetColumn(int i);
	char Existing(int i, int j);
	double* Scanning(int* i, int* j, double* val);
	ElementBzzMatrixSparse* GetStartingElementInRow(int i);
	void BeginScanning(void);
	int CountElements(void);
	int CountElementsOutsideLowerBand(int lowerB);
	int CountElementsOutsideUpperBand(int upB);
	void CountElementsOutsideBands(int lowerB, int upB, int* low, int* up);
	void CountElementsOutsideBands(int* l, int* u, BzzVectorInt* low, BzzVectorInt* up);
	void BzzPrintExistingElements(void);
	void BzzPrintStructure(void);
	void GetRowsExistence(BzzVectorInt* iRows);
	void GetNumEquationsForEachVariable
	(BzzVectorInt* numEquationsForEachVariable);
	void GetNumVariablesForEachEquation
	(BzzVectorInt* numVariablesForEachEquation);
	void GetNumVariablesForEachEquationAndNumEquationsForEachVariable
	(BzzVectorInt* numVariablesForEachEquation,
		BzzVectorInt* numEquationsForEachVariable);
	void GetEquationsContainingSingletons
	(BzzVectorInt* numVariablesForEachEquation,
		BzzVectorInt* numEquationsForEachVariable,
		BzzVectorInt* equationsContainingSingletons);
	void GetLinearlyDependentCouplesOfEquations
	(BzzVectorInt* firstEquation, BzzVectorInt* secondEquation,
		double cosAlfa = .99999);
	void GetPositiveLinearlyDependentCouplesOfEquations
	(BzzVectorInt* firstEquation, BzzVectorInt* secondEquation,
		double cosAlfa = .99999);
	void GetNegativeLinearlyDependentCouplesOfEquations
	(BzzVectorInt* firstEquation, BzzVectorInt* secondEquation,
		double cosAlfa = -.99999);
	void DecoupleEqualitySystem(BzzVector* e, double cosAlfa = .999);

	//	============================================================================
	//	**************************< assignment operators >**************************
	//	============================================================================
	BzzMatrixSparse& operator =
		(const BzzMatrixSparse& rval);
	BzzMatrixSparse& operator =
		(BzzMatrix& rval);

	//	============================================================================
	//	=============================< OPERATIONS >=================================
	//	============================================================================

	//	============================================================================
	//	*********************************< Sum >************************************
	//	============================================================================

	friend BzzMatrixSparse operator +
		(BzzMatrixSparse& lval,
			BzzMatrixSparse& rval);

	BzzMatrixSparse& operator +=
		(BzzMatrixSparse& rval);

	//	============================================================================
	//	******************************< Difference >********************************
	//	============================================================================

	friend BzzMatrixSparse operator -
		(BzzMatrixSparse& lval,
			BzzMatrixSparse& rval);

	BzzMatrixSparse& operator -=
		(BzzMatrixSparse& rval);

	//	============================================================================
	//	********************************< Product >*********************************
	//	============================================================================
		// Product(A,B,&C);
	friend void Product(BzzMatrixSparse& lval,
		BzzMatrixSparse& rval, BzzMatrixSparse* result);

	// Product(A,B,&C);
	friend void Product(BzzMatrixSparse& lval,
		BzzMatrixSparse& rval, BzzMatrix* result);

	// Product(A,B,&C); // A, C BzzMatrix
	friend void Product(BzzMatrix& lval,
		BzzMatrixSparse& rval, BzzMatrix* result);

	// Product(A,B,&C); // B, C BzzMatrix
	friend void Product(BzzMatrixSparse& lval,
		BzzMatrix& rval, BzzMatrix* result);

	// Product(A,B,&C); // A BzzMatrix
	friend void Product(BzzMatrix& lval,
		BzzMatrixSparse& rval, BzzMatrixSparse* result);

	// Product(A,B,&C); // B BzzMatrix
	friend void Product(BzzMatrixSparse& lval,
		BzzMatrix& rval, BzzMatrixSparse* result);

	// Product(A,x,&y);
	friend void Product
	(BzzMatrixSparse& lval, BzzVector& rval,
		BzzVector* result);

	// Product(A,&x);
	friend void Product(BzzMatrixSparse& lval,
		BzzVector* rvalAndResult);

	// y = A*x;
	friend BzzVector operator *
		(BzzMatrixSparse& lval, BzzVector& rval);

	// Product(3.,&A);
	friend void Product(double c, BzzMatrixSparse* result);

	//	============================================================================
	//	*******************************< TProduct >*********************************
	//	============================================================================
		// TProduct(A,B,&C); // A,B,C BzzMatrixSparse
	friend void TProduct(BzzMatrixSparse& lval,
		BzzMatrixSparse& rval, BzzMatrixSparse* result);

	// TProduct(A,x,&y);
	friend void TProduct
	(BzzMatrixSparse& lval, BzzVector& rval,
		BzzVector* result);

	// TProduct(A,&x);
	friend void TProduct(BzzMatrixSparse& lval,
		BzzVector* rvalAndResult);

	// y = A%x;
	friend BzzVector operator %
		(BzzMatrixSparse& lval, BzzVector& rval);

	//	============================================================================
	//	*******************************< ProductT >*********************************
	//	============================================================================
		// ProductT(A,B,&C); // A,B,C BzzMatrixSparse
	friend void ProductT(BzzMatrixSparse& lval,
		BzzMatrixSparse& rval, BzzMatrixSparse* result);

	//	============================================================================
	//	========================< Non-modifying functions >=========================
	//	============================================================================

	//	*******************************< BzzPrint >*********************************
	virtual void ObjectBzzPrint(void);
	void	Save(char* filematrix); // formatted
	void	Save(char, char* filematrix); // non formatted
	void FindBands(int* low, int* up);

	//	============================================================================
	//	===========================< Modifying Functions >==========================
	//	============================================================================
	void RemoveElement(int row, int column);
	void RemoveAllElementsInRow(int row);
	void RemoveAllElementsInColumn(int column);
	friend void Delete(BzzMatrixSparse* A);
	friend void Transpose(BzzMatrixSparse* A);
	void RemoveAllElementsInMatrix(void);// eliminates all elements
	void CleanMatrix(double eps); // eliminates those <= eps
	void ReleaseElement(int row, int column, BzzMatrixSparseStore& store);
	void ReleaseAllElementsInRow(int row, BzzMatrixSparseStore& store);
	void ReleaseAllElementsInColumn(int column, BzzMatrixSparseStore& store);
	void ReleaseAllElementsInMatrix(BzzMatrixSparseStore& store);
	void InsertStoredElement(int row, int column, double value,
		BzzMatrixSparseStore& store);
	void InsertStoredElement(int row, int column, double value,
		ElementBzzMatrixSparse* elem);
	void AddOneColumn(void) { numColumns++; };
	void AddNColumns(int n) { numColumns += n; }
	void AddOneRow(void);
	void AddNRows(int n);
	void DeleteColumn(int column);
	void DeleteRow(int row);
	void DeleteLastNRows(int n);
	void Append(BzzMatrixSparse* A);
	void Merge(BzzMatrixSparse* A);
	void SplitByRows(BzzMatrixSparse* A, int m);
	void SplitByColumns(BzzMatrixSparse* A, int n);
	friend void MoveRowsFromFirstMatrixAndAppendToSecond(BzzVectorInt& iRows,
		BzzMatrixSparse* A, BzzMatrixSparse* B);
	friend void ChangeDimensions(int nr, int nc,
		BzzMatrixSparse* m);
	friend void Swap(BzzMatrixSparse* lval,
		BzzMatrixSparse* rval);
	friend void Load(BzzMatrixSparse* A, char* filevector);
	friend void Load(BzzMatrixSparse* A, char, char* filevector);
	void SetDiagonal(int i, double c);
	void SwapRows(int i, int j);
	void SwapColumns(int i, int j);
	void SwapElements(int i1, int j1, int i2, int j2);

	//	============================================================================
	//	==============================< Norms >=====================================
	//	============================================================================
	double NormI(void);
	void GetRowsNorm2(BzzVector* norm2R);
	void GetRowsNorm1(BzzVector* norm1R);
	void GetColumnsNorm2(BzzVector* norm2C);
	void GetRowsNorm2X(BzzVector& x, BzzVector* norm2X);
	void GetRowsSum(BzzVector* sumR);
	void GetColumnsSum(BzzVector* sumC);

	//	============================================================================
	//	======================< System solution functions >=========================
	//	============================================================================
	void SolveRight(BzzVector* bx);
	void SolveLeft(BzzVector* bx);
	void SolveLeft(BzzVector* bx, char);

	void TransposeSolveRight(BzzVector* bx);
	void TransposeSolveLeft(BzzVector* bx);

	void TransposeSolveRight(BzzVector* bx, char);
	void TransposeSolveLeft(BzzVector* bx, char);

	friend void Solve(BzzMatrixSparse* A, BzzVector* bx);
	double ConditionNumberRight(void);
	double ConditionNumberLeft(void);
	void OrderingGauss(BzzMatrixSparse* A, BzzVectorInt* ordRows,
		BzzVectorInt* ordColumns);
	void OrderingGauss(BzzMatrix* A, BzzVectorInt* ordRows,
		BzzVectorInt* ordColumns);
	void OrderingGauss(BzzVectorInt* ordRows, BzzVectorInt* ordColumns);
	void OrderingGaussWithoutPivoting(BzzMatrixSparse* A, BzzVectorInt* ordRows,
		BzzVectorInt* ordColumns);
	void OrderingGaussWithoutPivoting(BzzMatrix* A, BzzVectorInt* ordRows,
		BzzVectorInt* ordColumns);
	void OrderingGaussWithoutPivoting(BzzVectorInt* ordRows, BzzVectorInt* ordColumns);
	int GaussFillingUp(void);

	void OrderingLinearProgramming(BzzVectorInt* ordRows,
		BzzVectorInt* ordColumns);
	void OrderingLQ(BzzVectorInt* ordRows, BzzVectorInt* ordColumns);
	int ChoiceLinearProgrammingVariables(BzzVector& unorm2V,
		BzzVectorInt* varForFactorization);
};

//	============================================================================
//	=============================< Modifying functions >========================
//	============================================================================
/*
	void SwapRCV(BzzVectorInt *rr,BzzVectorInt *cc,BzzVector *vv);
	void ReorderingRows(BzzVectorInt &nro);
	void ReorderingColumns(BzzVectorInt &nco);
	void ReorderingRowsAndColumns(BzzVectorInt &nro,BzzVectorInt &nco);

	void CleanRows(BzzVectorInt *dr);
	void CleanColumns(BzzVectorInt *dc);
	void CleanColumnsAndEvaluateResiduals(BzzVectorInt *dc,BzzVector *vc,
		BzzVector *res);
	void EvaluateResidualsForSelectedColumns(BzzVectorInt *dc,BzzVector *vc,
		BzzVector *res);
	void EvaluateResidualsAndChangeSignForSelectedColumns(BzzVectorInt *dc,BzzVector *vc,
		BzzVector *res);
//	void CleanRowsAndColumns(BzzVectorInt *dr,BzzVectorInt *dc);
	void DeleteRows(BzzVectorInt &dr);
	void DeleteRows(BzzVectorInt *dr);
	void DeleteColumns(BzzVectorInt *dc);
//	void DeleteRowsAndColumns(BzzVectorInt *dr,BzzVectorInt *dc);
*/

#endif // BZZ_SPARSE_DOUBLE_HPP