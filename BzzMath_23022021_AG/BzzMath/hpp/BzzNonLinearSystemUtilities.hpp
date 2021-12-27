// BZZMATH: Release 7.0

//	=======================< BzzNonLinearSystemUtilities.hpp >=====================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	04-2013	Date Written

#ifndef NON_LINEAR_SYSTEM_UTILITIES_HPP
#define NON_LINEAR_SYSTEM_UTILITIES_HPP

//	============================================================================
//	=================< class BzzNonLinearSystemJacobian >=======================
//	============================================================================

class BzzNonLinearSystemJacobianAndHessians
{
private:
	double zerDer, eta2, hj, xh, xdh, hInv;
	BzzVector aaa, bbb, ccc, fp;
	BzzVectorInt xType;
public:
	void (*ptrSystem)(BzzVector& xx, BzzVector& ff);
	void GetJacobian(int jStart, int jEnd, BzzVector& xL,
		BzzVector& xU, BzzVectorInt& xType, BzzVector& xi, BzzVector& fi, BzzMatrix& Ji);
	void UpdateJacobian(int jStart, int jEnd, BzzVectorInt& xType, BzzVector& dxi,
		BzzVector& dfi, BzzMatrix& Ji);
	void GetJacobianAndDiagonalHessians(int jacobian, int jStart, int jEnd, BzzVector& xL,
		BzzVector& xU, BzzVectorInt& xType, BzzVector& xi, BzzVector& fi, BzzMatrix& Ji,
		BzzVector& h, BzzMatrix& Fp, BzzMatrixSymmetric* G);
	void GetLowerDiagonalHessians(int jStart, int jEnd, BzzVector& xL,
		BzzVector& xU, BzzVectorInt& xType, BzzVector& xi, BzzVector& fi, BzzMatrix& Ji,
		BzzVector& h, BzzMatrix& Fp, BzzMatrixSymmetric* G);
};

//	============================================================================
//	==================< class BzzNonLinearSystemUtilities >=====================
//	============================================================================

class BzzNonLinearSystemUtilities : public BzzBaseClass
{
	static const char* const BZZ_ERROR;

	int	numVariables,
		numEquations,
		numThread,
		starting,
		jacobian,
		hessians,
		denseSparse; // 0 default 1 dense 2 sparse

	BzzVector	xL,
		xU,
		xJacobian,
		fJacobian;

	// openMP
	BzzVectorInt iStart, iEnd;

	// dense
	BzzVectorInt xType, xTypeOriginal;
	BzzNonLinearSystemJacobianAndHessians* getJac;

	// sparse
	BzzVectorIntArray xTypeArrayOriginal;
	BzzVectorInt rOriginal, cOriginal, raux, caux;
	BzzVector vOriginal, vaux;
	int numElementsOriginal;

	void (*ptrSystem)(BzzVector& xx, BzzVector& ff);
	void Initialize(void);

public:
	// dense
	BzzVector dx, df;
	BzzVector h;
	BzzMatrix J, Fp;
	BzzMatrixSymmetric* G;

	// sparse
	BzzMatrixSparseLocked JS;

	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================

		// default constructor
	BzzNonLinearSystemUtilities(void);

	// operator()
	void operator()(int numE, BzzVector& xxL, BzzVector& xxU, BzzVectorInt& xType,
		void (*ptrS)(BzzVector& xx, BzzVector& ff));
	void operator()(int numE, BzzVector& xxL, BzzVector& xxU,
		BzzVectorIntArray& xTypeArray, void (*ptrS)(BzzVector& xx, BzzVector& ff));

	//	============================================================================
	//	*********************************< destructor >*****************************
	//	============================================================================
	~BzzNonLinearSystemUtilities(void)
	{
		if (numEquations != 0)
#if BZZ_COMPILER != 0
			delete[] G;
#else
			//delete[numEquations + 1] G;
			delete[] G;
#endif

		if (numThread != 0)
#if BZZ_COMPILER != 0
			delete[] getJac;
#else
			//delete[numThread] getJac;
			delete[] getJac;
#endif
	}

	//	============================================================================
	//	******************************< Functions >*********************************
	//	============================================================================
	void SetDenseMatrix(void)
	{
		denseSparse = 1;
	}
	void SetSparseMatrix(void)
	{
		denseSparse = 2;
	}

	void SetFixedVariables(BzzVectorInt& fV);
	void BuildJacobian(BzzVector& xi, BzzVector& fi);
	int UpdateJacobian(BzzVector& xi, BzzVector& fi); // 0 not updated 1 updated 2 rebuilded
	void BuildJacobianAndDiagonalHessians(BzzVector& xi, BzzVector& fi);
	void BuildLowerDiagonalHessians(BzzVector& xi, BzzVector& fi);
	void BuildJacobianAndHessians(BzzVector& xi, BzzVector& fi);

	virtual void ObjectBzzPrint(void);
};

#endif // NON_LINEAR_SYSTEM_UTILITIES_HPP