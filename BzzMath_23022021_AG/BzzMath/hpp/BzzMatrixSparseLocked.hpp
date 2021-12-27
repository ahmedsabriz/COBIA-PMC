// BZZMATH: Release 7.0

//	==============< BzzMatrixSparseLocked.hpp >========================
//	* BzzMatrixSparseLocked: Class for sparse matrices*
// ============================================================================
// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	11-2010	Date Written.
//	11-2012	Inserted RemoveEquationsAndVariables function.
//	11-2012	Inserted RemoveEquations function.
//	11-2012	Inserted RemoveVariables function.
// 02-2013  Modified

#ifndef BZZ_MATRIX_SPARSE_LOCKED
#define BZZ_MATRIX_SPARSE_LOCKED

class BzzMatrixSparseLocked : public BzzBaseClass
{
	friend class BzzSave;
	friend class BzzLoad;
	friend class BzzMatrixCoefficientsLocked;
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	static const double HESSIAN_TOLERANCE;
	static const double DX_TOLERANCE;
	int numRows, numColumns;
	int whoAmI;
	int findDependence;

	// For Sum
	BzzVector aux, aux1, aux2; // for Jacobian
	BzzVectorInt iaux; // also for FindDependence
	BzzVectorInt jaux; // also for FindDependence
	BzzVectorInt rx, ry, rz, re, rw;
	BzzVector vx, vy, vz, ve, vw;
	BzzVectorInt nonLinearities; // for Jacobian and Hessians

	// for VariablesAndEquationsOrdering
	BzzVectorInt naux;

	// for FindDependence
	int	maxElementsInRow, iMaxRow,
		minElementsInRow, iMinRow,
		maxElementsInColumn, iMaxColumn,
		minElementsInColumn, iMinColumn,
		lowerBand,
		upperBand,
		numNullRows,
		numNullColumns,
		numSingletonRows,
		numSingletonColumns,
		numAssignedVariablesForSingletons,
		reduced,
		reducedLowerBand,
		reducedUpperBand,
		numElementsUpPrincipalDiagonal,
		numElementsDownPrincipalDiagonal,
		numJacobianCalls, // 0 first call
		numHessiansCalls; // 0 first call

	int BnumGroup;
	BzzVectorInt BnumVariablesInGroup;
	BzzVectorInt BgroupForEachVariable;
	BzzVectorIntArray BvariablesInGroup;
	BzzVectorInt	BnumVariablesInEachEquation,
		BnumEquationsForEachVariable,
		BnumLinearVariablesInEachEquation;
	BzzVectorIntArray BcI;

	double hessiansTolerance, dxTolerance;

	// for Jacobian
	BzzVector xJacobian, fJacobian;
	BzzVectorInt	variablesNonPresent;
	// for hessians
	BzzVector h;
	BzzVector v0, vp;
	int numElementsForHessians;
	BzzVectorInt hessian, rh, ch;
	BzzVector vh;

	//	BzzVector x0,x,h,fp,fa,fm;
	//	double f0;
	//	double (*ptrFun)(BzzVector &x);

	////////////////////////////////////////////
public:
	int	numGroup;
	int numElements;
	//  1. 0. 3. 0. 0. 6.
	//  0. 7. 0. 0. 8. 0.
	//  2. 0. 0. 0. 5. 0.
	//  0. 0. 0. 0. 0. 0.
	// Locked
	BzzVector v; // v(7,1.,3.,6.,8.,5.,2.,7.)
	BzzVectorInt r, // r(7,1,1,1,2,3,3,2)
		c; // c(7,1,3,6,5,5,1,2)
private:
	void BuildRCVFromRows(void);
	void BuildRCVFromColumns(void);

public:
	// LockedByRows
	BzzVectorIntArray rI;
	// rI(1)(3,1,3,6)
	// rI(2)(2,5,2)
	// rI(3)(2,1,5)
	// rI(4)(0)
	BzzVectorIntArray rNL; // contiene le colonne c[k] + 10000000 * nonLinearities[k]
	BzzVectorArray rV;
	// rV(1)(3,1.,3.,6.)
	// rV(2)(2,8.,7.)
	// rV(3)(2,2.,5.)
	// rV(4)(0)
	BzzVectorIntArray rIL; // contiene le colonne lineari in quella riga
	BzzVectorArray rVL; // contiene il valore dello Jacobiano lineare

private:
	void BuildRowsFromRCV(void); // r,c,v->rI,rV

public:
	// LockedByColumns
	BzzVectorIntArray cI;
	// cI(1)(2,1,3)
	// cI(2)(1,2)
	// cI(3)(1,1)
	// cI(4)(0)
	// cI(5)(2,3,2)
	// cI(6)(1,1)
	BzzVectorArray cV;
	// cV(1)(2,1.,2.)
	// cV(2)(1,7.)
	// cV(3)(1,3.)
	// cV(4)(0)
	// cV(5)(2,5.,8.)
	// cV(6)(1,6.)
private:
	void BuildColumnsFromRCV(void); // r,c,v->cI,cV

	void Initialize(void);
	BzzVectorInt auxi; // serve per i Build
//	BzzHessians *hessians;
	BzzVectorInt rO;
	BzzVectorInt cO;
	BzzVector vO;

public:
	int numLinear, numQuadratic, numNonlinear;
	BzzVectorInt newRowsOrderBeforeAutoSwap, newColumnsOrderBeforeAutoSwap;

	BzzVectorInt sortNonlinearities;
	BzzVectorInt	numVariablesInEachEquation,
		numEquationsForEachVariable,
		numLinearVariablesInEachEquation;
	BzzVectorInt newVariablesNumber;
	int numHessiansElements;
	BzzVectorInt hHessian, rHessian, cHessian;
	BzzVector vHessian;
	void SetHessiansTolerance(double tol)
	{
		hessiansTolerance = tol;
	}
	void SetDxTolerance(double tol)
	{
		dxTolerance = tol;
	}
	void SetNonlinearities(BzzVectorInt* nl)
		//		{Swap(&nonLinearities,nl);}
	{
		(*this)(numRows, numColumns, &r, &c, &v, nl);
	}
	// for FindDependence
// dependence[var][j] var = 1,numColumns; j = 1,numEqWithSuchVariable[var]
//	BzzVectorIntArray dependence; // � uguale a cI

 // variableInGroup[group][i] group = 1,numGroup; i = 1,numVariablesInGroup[group]
	BzzVectorIntArray variablesInGroup;

	// numVariablesInGroup[group] group = 1,numGroup;
	BzzVectorInt numVariablesInGroup;

	BzzVectorInt groupForEachVariable;

	//	int	**dependence,
	//			*numEqWithSuchVariable, >> numEquationsForEachVariable
	//		 	**variableInGroup,
	//			numGroup,
	//			*numVariablesInGroup;

	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzMatrixSparseLocked(void);

	// dimensions constructor
	BzzMatrixSparseLocked(int nr, int nc);
	void operator()(int nr, int nc);
	void operator()(BzzVectorInt* rr, BzzVectorInt* cc,
		BzzVector* vv);

	// copy constructor
	BzzMatrixSparseLocked(const BzzMatrixSparseLocked& rval);

	// constructor
	BzzMatrixSparseLocked(int nr, int nc, BzzVectorInt* rr, BzzVectorInt* cc,
		BzzVector* vv);
	void operator()(int nr, int nc, BzzVectorInt* rr, BzzVectorInt* cc,
		BzzVector* vv);

	// constructor
	BzzMatrixSparseLocked(int nr, int nc, BzzVectorInt* rr, BzzVectorInt* cc,
		BzzVector* vv, BzzVectorInt* nl);
	void operator()(int nr, int nc, BzzVectorInt* rr, BzzVectorInt* cc,
		BzzVector* vv, BzzVectorInt* nl);

	// band matrix
	BzzMatrixSparseLocked(int nr, int nc, int low, int up);
	void operator()(int nr, int nc, int low, int up);

	// BzzMatrixSparseLocked A("SPAR.DAT");
	BzzMatrixSparseLocked(char* filematrix);

	// return 1 only if numElements * max < numRows * numColumns
	// in this case the matrix is replaced

	friend int ReplaceBzzMatrixWithBzzMatrixSparseLocked(int max,
		BzzMatrix* A, BzzMatrixSparseLocked* B);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzMatrixSparseLocked(void);

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
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }
	int GetGroupNumber(void) { return numGroup; }

	//	============================================================================
	//	=========================< Non-modifying functions >========================
	//	============================================================================

	//	*******************************< BzzPrint >*********************************

	virtual void ObjectBzzPrint(void);
	void BzzPrintByColumns(void);
	void BzzPrintStructureByRows(void);
	void BzzPrintStructureByColumns(void);
	void PlotStructure(void);

	void GetNumVariablesInEachEquation(BzzVectorInt* ve)
	{
		*ve = numVariablesInEachEquation;
	}
	void GetNumEquationsForEachVariable(BzzVectorInt* ev)
	{
		*ev = numEquationsForEachVariable;
	}
	void GetNumVariablesForEachEquationAndNumEquationsForEachVariable(BzzVectorInt* ve,
		BzzVectorInt* ev)
	{
		*ve = numVariablesInEachEquation;*ev = numEquationsForEachVariable;
	}
	//	double MaxAbsColumnForSelectedRows(int j,BzzVectorInt &row,int *i);
	//	void BzzColumnBalanceForSelectedRows(BzzVectorInt &row);
	//	void GetRCV(BzzVectorInt *rr,BzzVectorInt *cc,BzzVector *vv);
	//
	void	Save(char* filematrix); // formatted
//
//	void	Save(char,char *filematrix); // non formatted
	void FindBands(int* low, int* up);
	void FindReducedBands(int reduced, int* low, int* up);
	void FindDependence(void);
	void PrintDependence(void);
	void BuildMatrixByRows(BzzVectorInt& rowsB, BzzMatrixSparseLocked* B);
	void BuildJacobian(BzzVector& x, BzzVector& f,
		void SystemName(BzzVector& x, BzzVector& f));
	void JacobianUpdate(BzzVector& x1, BzzVector& f1,
		void SystemName(BzzVector& x, BzzVector& f));
	//	void JacobianUpdate(BzzVector &x1,BzzVector &f1,
	//		void SystemName(BzzVector &x,BzzVector &f),
	//		BzzMatrixSparseLocked &toModify);
	void BuildJacobianAndHessians(BzzVector& x, BzzVector& f,
		void SystemName(BzzVector& x, BzzVector& f));
	void BuildJacobianAndHessiansUpdate(BzzVector& x, BzzVector& f,
		void SystemName(BzzVector& x, BzzVector& f));
	void JacobianAndHessiansUpdate(BzzVector& x, BzzVector& f,
		void SystemName(BzzVector& x, BzzVector& f));

	void operator =
		(BzzMatrixSparseLocked& B);
	friend void Swap(BzzMatrixSparseLocked* A, BzzMatrixSparseLocked* B);

	int GetRowsWithSingletons(BzzVectorInt* rS);

	//	============================================================================
	//	=============================< Modifying functions >========================
	//	============================================================================
	void ReorderingRows(BzzVectorInt& nro);
	void ReorderingColumns(BzzVectorInt& nco);
	void ReorderingRowsAndColumns(BzzVectorInt& nro, BzzVectorInt& nco);
	void VariablesAndEquationsOrdering(BzzVectorInt* newRowsOrder,
		BzzVectorInt* newColumnsOrder);
	//	void FindDoublets(int version);
	void FindDoublets(int version, BzzVector* rhs);

	void FindBaricenter(BzzVectorInt& equationType, BzzVector& rhs, BzzVector& lowerBounds,
		BzzVector& upperBounds, BzzVector* baricenter);
	void RemoveEquations(BzzVectorInt& rr);
	void SetEquation(int eq, BzzVectorInt* cc, BzzVector* vv);
	void RemoveVariables(BzzVectorInt& cc);
	void SetColumn(int col, BzzVectorInt* rr, BzzVector* vv);
	void RemoveEquationsAndVariables(BzzVectorInt& rr, BzzVectorInt& cc);

	void SwapRows(int i, int j);
	void SwapColumns(int i, int j);

	void MoveRowFromkToj(int k, int j);
	void MoveColumnFromkToj(int i, int j);

	void SortVariablesInRows(void);
	void SortEquationsInColumns(void);

	friend void Transpose(BzzMatrixSparseLocked* A);
	void ModifyColumnsElements(BzzMatrixSparseLocked& B);

	void Append(BzzMatrixSparseLocked& A, BzzMatrixSparseLocked& B);
	void Append(BzzMatrixSparseLocked& A, BzzMatrixSparseLocked& B, BzzVectorInt& rowsB);
	void Append(BzzMatrixSparseLocked& A, BzzVectorInt& rowsA, BzzMatrixSparseLocked& B);
	void Append(BzzMatrixSparseLocked& A, BzzVectorInt& rowsA, BzzMatrixSparseLocked& B,
		BzzVectorInt& rowsB);

	void Append(BzzMatrixSparseLocked& A, BzzVector& a,
		BzzMatrixSparseLocked& B, BzzVector& b, BzzVector* c);
	void Append(BzzMatrixSparseLocked& A, BzzVector& a,
		BzzMatrixSparseLocked& B, BzzVector& b, BzzVectorInt& rowsB, BzzVector* c);
	void Append(BzzMatrixSparseLocked& A, BzzVector& a,
		BzzVectorInt& rowsA, BzzMatrixSparseLocked& B, BzzVector& b, BzzVector* c);
	void Append(BzzMatrixSparseLocked& A, BzzVector& a,
		BzzVectorInt& rowsA, BzzMatrixSparseLocked& B, BzzVector& b, BzzVectorInt& rowsB, BzzVector* c);
	void SplitByRows(BzzVectorInt& rowsB, BzzMatrixSparseLocked* B);
	void SplitByColumns(BzzVectorInt& columnsB, BzzMatrixSparseLocked* B);
	void SplitByRowsAndColumns(BzzVectorInt& rowsB, BzzVectorInt& columnsB,
		BzzMatrixSparseLocked* B);

	//	============================================================================
	//	*************************< Algebraic Operations >***************************
	//	============================================================================

	//	============================================================================
	//	********************************< Product >*********************************
	//	============================================================================
	friend void Product(double c, BzzMatrixSparseLocked& A, BzzMatrixSparseLocked* B);
	friend void Product(double c, BzzMatrixSparseLocked* A);

	friend void Product(BzzMatrixSparseLocked& SL, BzzVector& x, BzzVector* y);
	friend void Product(BzzMatrixSparseLocked& SL, BzzVectorSparse& x, BzzVector* y);

	friend void ProductForSelectedColumns(BzzMatrixSparseLocked& SL, BzzVector& x,
		BzzVectorInt& col, BzzVector* y);
	friend void ProductForSelectedColumns(BzzMatrixSparseLocked& SL, BzzVectorSparse& x,
		BzzVectorInt& col, BzzVector* y);

	friend void ProductForSelectedRows(BzzMatrixSparseLocked& SL, BzzVector& x,
		BzzVectorInt& rows, BzzVector* y);
	friend void ProductForSelectedRows(BzzMatrixSparseLocked& SL, BzzVectorSparse& x,
		BzzVectorInt& rows, BzzVector* y);

	friend void ProductForSelectedRowsAndColumns(BzzMatrixSparseLocked& SL, BzzVector& x,
		BzzVectorInt& rows, BzzVectorInt& col, BzzVector* y);
	//	friend void ProductForSelectedRowsAndColumns(BzzMatrixSparseLocked &SL,BzzVectorSparse &x,
	//			BzzVectorInt &rows,BzzVectorInt &col,BzzVector *y);

	friend void Product(BzzMatrixSparseLocked& S1, BzzMatrixSparseLocked& S2,
		BzzMatrix* A);
	friend void Product(BzzMatrixSparseLocked& S1, BzzMatrixSparseLocked& S2,
		BzzMatrixSparseLocked* S3);

	//	============================================================================
	//	*******************************< ProductT >*********************************
	//	============================================================================
	friend void ProductT(BzzVectorSparse& x, BzzVectorSparse& y,
		BzzMatrixSparseLocked* SL);

	//	============================================================================
	//	***********************************< Sum >**********************************
	//	============================================================================
	// A, B, C stessi r e c. Senza controllo
	friend void PerfectSum(BzzMatrixSparseLocked& A, BzzMatrixSparseLocked& B,
		BzzMatrixSparseLocked* C);

	//=====
	// TODO
	friend void Sum(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
		BzzMatrixSparseLocked& A, BzzMatrixSparseLocked& B,
		BzzMatrixSparseLocked* C);
	//friend void Sum(BzzMatrixSparseLocked &A,BzzMatrixSparseLocked &B,
	//		BzzMatrixSparseLocked *C);
	//=====

	//	============================================================================
	//	*******************************< Difference >*******************************
	//	============================================================================
	// A, B, C stessi r e c. Senza controllo
	friend void PerfectDifference(BzzMatrixSparseLocked& A, BzzMatrixSparseLocked& B,
		BzzMatrixSparseLocked* C);

	//	============================================================================
	//	=================< Functions for Gauss Factorization >======================
	//	============================================================================
	void FindPivot(void);
};

#endif // BZZ_MATRIX_SPARSE_LOCKED