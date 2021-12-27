//=====================================================================================
//=======================================================================================
//===================================================================================
//	==============< BzzMatrixCoefficientsLocked.hpp >========================
//	* BzzMatrixCoefficientsLocked: Class for sparse matrices analysis*
// ============================================================================
// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	04-2013	Date Written.

#ifndef BZZ_MATRIX_COEFFICIENTS_LOCKED
#define BZZ_MATRIX_COEFFICIENTS_LOCKED

class BzzMatrixCoefficientsLocked : public BzzBaseClass
{
	friend class BzzSave;
	friend class BzzLoad;
	friend class BzzMatrixSparseLocked;
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	int numRows, numColumns, numElements;
	int whoAmI;
	int findDependence;

	BzzVectorInt	variablesNonPresent;
	// For Sum
	BzzVector aux, aux1; // for Jacobian
	BzzVectorInt iaux; // also for FindDependence
	BzzVectorInt jaux; // also for FindDependence
	BzzVectorInt rz;
	BzzVector vz;

	// for VariablesAndEquationsOrdering
	BzzVectorInt naux, numVariablesInOneEqualityOnly;
	int numXE1Variables, numXE1Equations;
	int noSingleton;

	// for FindDependence
	int	numGroup,
		maxElementsInRow, iMaxRow,
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
		reducedLowerBand,
		reducedUpperBand,
		numElementsUpPrincipalDiagonal,
		numElementsDownPrincipalDiagonal;
	/*
		// for Jacobian
		BzzVector xJacobian,fJacobian;
		BzzVectorInt	variablesNonPresent;
		// for hessians
		BzzVector h;
		BzzVector v0,vp;
		int numElementsForHessians;
		BzzVectorInt hessian,rh,ch;
		BzzVector vh;

	//	BzzVector x0,x,h,fp,fa,fm;
	//	double f0;
	//	double (*ptrFun)(BzzVector &x);

	*/
	////////////////////////////////////////////

		//  1. 0. 3. 0. 0. 6.
		//  0. 7. 0. 0. 8. 0.
		//  2. 0. 0. 0. 5. 0.
		//  0. 0. 0. 0. 0. 0.
		// Locked
	BzzVectorInt r, // r(7,1,1,1,2,3,3,2)
		c; // c(7,1,3,6,5,5,1,2)
	void BuildRCFromRows(void);
	void BuildRCFromColumns(void);

	// LockedByRows
	BzzVectorIntArray rI;
	void BuildRowsFromRC(void); // r,c,v->rI,rV

	// LockedByColumns
	BzzVectorIntArray cI;
	void BuildColumnsFromRC(void); // r,c,v->cI,cV

	void Initialize(void);
	BzzVectorInt auxi; // serve per i Build
//	BzzHessians *hessians;

public:
	BzzVectorInt newRowsOrderBeforeAutoSwap, newColumnsOrderBeforeAutoSwap;
	BzzVectorInt	numVariablesInEachEquation,
		numEquationsForEachVariable;
	BzzVectorInt newVariablesNumber;
	BzzVectorIntArray variablesInGroup;
	BzzVectorInt numVariablesInGroup;
	BzzVectorInt groupForEachVariable;

	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzMatrixCoefficientsLocked(void);

	// dimensions constructor
	BzzMatrixCoefficientsLocked(int nr, int nc);

	// copy constructor
	BzzMatrixCoefficientsLocked(const BzzMatrixCoefficientsLocked& rval);

	// constructor
	BzzMatrixCoefficientsLocked(int nr, int nc, BzzVectorInt* rr, BzzVectorInt* cc);

	void operator()(int nr, int nc);
	void operator()(int nr, int nc, BzzVectorInt* rr, BzzVectorInt* cc);
	void operator()(BzzVectorInt* rr, BzzVectorInt* cc);

	// BzzMatrixCoefficientsLocked A("SPAR.DAT");
	BzzMatrixCoefficientsLocked(char* filematrix);

	// return 1 only if numElements * max < numRows * numColumns
	// in this case the matrix is replaced

//	friend int ReplaceBzzMatrixWithBzzMatrixCoefficientsLocked(int max,
//		BzzMatrix *A,BzzMatrixCoefficientsLocked *B);

// ============================================================================
// *****************************< destructor >*********************************
// ============================================================================
	~BzzMatrixCoefficientsLocked(void);

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
	void SetNoSingleton(void)
	{
		noSingleton = 1;
	}
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
	int GetNumSingletonRows(void)
	{
		return numSingletonRows;
	}
	void FindDependence(void);
	void PrintDependence(void);
	/*
		void BuildJacobian(BzzVector &x,BzzVector &f,
			void SystemName(BzzVector &x,BzzVector &f));
		void JacobianUpdate(BzzVector &x1,BzzVector &f1,
			void SystemName(BzzVector &x,BzzVector &f));
	//	void JacobianUpdate(BzzVector &x1,BzzVector &f1,
	//		void SystemName(BzzVector &x,BzzVector &f),
	//		BzzMatrixCoefficientsLocked &toModify);
		void BuildJacobianAndHessians(BzzVector &x,BzzVector &f,
			void SystemName(BzzVector &x,BzzVector &f));
		void BuildJacobianAndHessiansUpdate(BzzVector &x,BzzVector &f,
			void SystemName(BzzVector &x,BzzVector &f));
		void JacobianAndHessiansUpdate(BzzVector &x,BzzVector &f,
			void SystemName(BzzVector &x,BzzVector &f));
	*/
	void operator =
		(BzzMatrixCoefficientsLocked& B);
	void operator =
		(BzzMatrixSparseLocked& B);
	friend void Swap(BzzMatrixCoefficientsLocked* A, BzzMatrixCoefficientsLocked* B);

	//	============================================================================
	//	=============================< Modifying functions >========================
	//	============================================================================
	void ReorderingRows(BzzVectorInt& nro);
	void ReorderingColumns(BzzVectorInt& nco);
	void ReorderingRowsAndColumns(BzzVectorInt& nro, BzzVectorInt& nco);
	void VariablesAndEquationsOrdering(BzzVectorInt* newRowsOrder,
		BzzVectorInt* newColumnsOrder);

	void RemoveEquations(BzzVectorInt& rr);
	void SetEquation(int eq, BzzVectorInt* cc);
	void RemoveVariables(BzzVectorInt& cc);
	void SetColumn(int col, BzzVectorInt* rr);
	void RemoveEquationsAndVariables(BzzVectorInt& rr, BzzVectorInt& cc);

	void SwapRows(int i, int j);
	void SwapColumns(int i, int j);

	void MoveRowFromkToj(int k, int j);
	void MoveColumnFromkToj(int i, int j);

	void SortVariablesInRows(void);
	void SortEquationsInColumns(void);

	friend void Transpose(BzzMatrixCoefficientsLocked* A);
	//	void ModifyColumnsElements(BzzMatrixCoefficientsLocked &B);
	void Append(BzzMatrixCoefficientsLocked& A, BzzMatrixCoefficientsLocked& B);
	void Append(BzzMatrixCoefficientsLocked& A, BzzMatrixCoefficientsLocked& B, BzzVectorInt& rowsB);
	void Append(BzzMatrixCoefficientsLocked& A, BzzVectorInt& rowsA, BzzMatrixCoefficientsLocked& B);
	void Append(BzzMatrixCoefficientsLocked& A, BzzVectorInt& rowsA, BzzMatrixCoefficientsLocked& B,
		BzzVectorInt& rowsB);
	void SplitByRows(BzzVectorInt& rowsB, BzzMatrixCoefficientsLocked* B);
	void SplitByColumns(BzzVectorInt& columnsB, BzzMatrixCoefficientsLocked* B);
	void SplitByRowsAndColumns(BzzVectorInt& rowsB, BzzVectorInt& columnsB,
		BzzMatrixCoefficientsLocked* B);
};

#endif // BZZ_MATRIX_COEFFICIENTS_LOCKED