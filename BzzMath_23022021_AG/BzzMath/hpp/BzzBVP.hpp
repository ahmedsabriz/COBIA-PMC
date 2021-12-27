// BZZMATH: Release 7.0

//	==========================< BzzBVP.hpp >==========================
//	* BzzBVP class: class for BVP				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	01-2014	Date Written

#ifndef BZZ_BVP_HPP
#define BZZ_BVP_HPP

//	===============================================================
//	======================< class BzzBVP >=========================
//	===============================================================
class BzzBVP : public BzzBaseClass
{
public:
	int	numEquations,
		numVariables,
		polynomiumOrder, originalPolynomiumOrder,
		numElements,
		numExteriorPoints,
		numSupportPointsForValues,
		numSupportPointsForDerivatives,
		numSupportPoints,
		numSupportPointsInEachElement,
		numConnectionPoints,
		iterCritical,
		maxNumCritical,
		iter;

	BzzFiniteElementPointsForHermitePolynomial fep;  // 3, 4, 5, 6

	BzzPiecewiseHermitePolynomialInterpolation* phpi;

	BzzPiecewiseLinearInterpolation pLin;

	BzzVector	xL, // punti iniziali in Linear
		sv, // SupportPointsForValues
		sd, // SupportPointsForDerivatives
		s, // SupportPoints for residuals
		deltaExteriorPoints,
		xc, // SupportPoints for contiguous points
		residuals, maxResiduals, vaResiduals, d1Residuals, d2Residuals,
		aux1, aux2;

	BzzMatrix	YL, // valori di Y iniziali in Linear
					// variabili del sistema non lineare
		Y, // valori di y per i polinomi Dimensions(numEquations x numSupportPointsForValues)
		D1, // valori di y' per i polinomi Dimensions(numEquations x numSupportPointsForDerivatives)

		// valori per i residui
		SY, // valori di y Dimensions(numSupportPoints x numVariables)
		SD1, // valori di y' Dimensions(numSupportPoints x numVariables)
		SD2; // valori di y" Dimensions(numSupportPoints x numVariables)

	BzzVector	y, // riga di Y (numSupportPointsForValues)
		d1, // riga di D1 Dimensions(numSupportPointsForDerivatives)

		sy, // riga di SY Dimensions(numVariables)
		sd1, // riga di SD1 Dimensions(numVariables)
		sd2; // riga di SD2 Dimensions(numVariables)

	BzzVector d2Left, d2Right;

	int initialFinal, fixedPolynomialOrder;
	int numVariablesNLS;
	BzzNonLinearSystemSparseObject nls;
	BzzVector z0, z; // Non Linear System variables

	BzzVectorInt	equationValuePresence, // ci sono o no y nei residui
		equationFirstDerivativePresence, // ci sono o no derivate prime
		equationSecondDerivativePresence, // ci sono o no derivate seconde
		leftBoundarySecondDerivativePresence, // ci sono o no derivate seconde
		rightBoundarySecondDerivativePresence, // ci sono o no derivate seconde
		numInternalPoints, numOriginalInternalPoints;

	int presenceOfSecondDerivativeInLeft; // 0 nessuna
	int presenceOfSecondDerivativeInRight; // 0 nessuna
	int constantParameter; // 1 se c'è 0 se non c'è

	BzzVector	weights,
		criticalParametersSoft,
		criticalParametersRequired,
		criticalParameters,
		criticalParametersOK, // valore dei parametri per la soluzione trovata
		scaleForPlot;
	//					zz;
	//	BzzMatrix ZZ;
		// bvpOK	< 0 bvp failed
		// bvpOK	== 0 bvp solved with softer critical parameters
		// bvpOK	== 1 bvp solved with the required critical parameters
	int criticalParametersPresent, solutionOK, bvpOK, nnxx;
	int criticalParOK; // 1 se giunge a soluzione per criticalParametersOK

	void(*ptrResidualsName)(double x, BzzVector& y,
		BzzVector& d1, BzzVector& d2, BzzVector& cP,
		BzzVector& r);
	void(*ptrLeftBoundaryName)(double x, BzzVector& y,
		BzzVector& d1, BzzVector& d2, BzzVector& cP,
		BzzVector& r, double* ri);
	void(*ptrRightBoundaryName)(double x, BzzVector& y,
		BzzVector& d1, BzzVector& d2, BzzVector& cP,
		BzzVector& r, double* rf);
	int ConstantParameter(int numIter);

public:
	//=================================< Functions >================================
	void SetPolynomiumOrder(int po)
	{
		fixedPolynomialOrder = 1;
		if (po < 3)
		{
			originalPolynomiumOrder = polynomiumOrder = 3;
		}
		else if (po > 6)
		{
			originalPolynomiumOrder = polynomiumOrder = 6;
		}
		else
		{
			originalPolynomiumOrder = polynomiumOrder = po;
		}
	}
	void SetScaleForPlot(BzzVector& sc)
	{
		scaleForPlot = sc;
	}

	int GetSolution(BzzVector* xV, BzzMatrix* yV)
	{
		if (bvpOK == 1)
		{
			*xV = sv;
			*yV = Y;
		}
		else if (bvpOK == 0)
		{
			BzzError("TODO bvpOK == 0");
		}
		return bvpOK;
	}
	//=========================< constructor >============================
	BzzBVP(void);
	//==================< Inizializer >================================
	void operator()(int nE, int nV,
		BzzVectorInt& evp, BzzVectorInt& ed1p, BzzVectorInt& ed2p,
		BzzVectorInt& ld2p, BzzVectorInt& rd2p,
		void (*ptrRes)(double x, BzzVector& y,
			BzzVector& d1, BzzVector& d2, BzzVector& cP,
			BzzVector& r),
		void (*ptrLeft)(double x, BzzVector& y,
			BzzVector& d1, BzzVector& d2, BzzVector& cP,
			BzzVector& r, double* ri),
		void (*ptrRigh)(double x, BzzVector& y,
			BzzVector& d1, BzzVector& d2, BzzVector& cP,
			BzzVector& r, double* rf),
		BzzVector& cPS, BzzVector& cPR,
		BzzVector& xx, BzzMatrix& YY, BzzVectorInt& nip, int initialFinal);

	int operator()(int numIter = 10);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzBVP(void);

	//==================================< BzzPrint >================================
	void ObjectBzzPrint(void);
};

//=============================================================================================
class MyNonLinearSystemObjectBVP : public BzzMyNonLinearSystemSparseObject
{
	friend class BzzBVP;
private:
	BzzBVP* ptrBVP;
public:
	MyNonLinearSystemObjectBVP(BzzBVP* ptr)
	{
		ptrBVP = ptr;
	}
	virtual void GetResiduals(BzzVector& x, BzzVector& f);
	virtual void ObjectBzzPrint(void);
};

//=============================================================================================
class MyNonLinearSystemObjectBVPConstant : public BzzMyNonLinearSystemSparseObject
{
	friend class BzzBVP;
private:
	BzzBVP* ptrBVP;
public:
	MyNonLinearSystemObjectBVPConstant(BzzBVP* ptr)
	{
		ptrBVP = ptr;
	}
	virtual void GetResiduals(BzzVector& x, BzzVector& f);
	virtual void ObjectBzzPrint(void);
};

#endif // BZZ_BVP_HPP

/*
// BZZMATH: Release 7.0

//	==========================< BzzBVP.hpp >==========================
//	* BzzBVP class: class for BVP				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	01-2014	Date Written

#ifndef BZZ_BVP_HPP
#define BZZ_BVP_HPP

//	===============================================================
//	======================< class BzzBVP >=========================
//	===============================================================
class BzzBVP : public BzzBaseClass
	{
public:
	int	numEquations,
			numVariables,
			polynomiumOrder,
			numElements,
			numExteriorPoints,
			numSupportPointsForValues,
			numSupportPointsForDerivatives,
			numSupportPoints,
			numSupportPointsInEachElement,
			numConnectionPoints;

	BzzFiniteElementPointsForHermitePolynomial fep;  // 3, 4, 5, 6

	BzzPiecewiseHermitePolynomialInterpolation *phpi;

	BzzPiecewiseLinearInterpolation pLin;

	BzzVector	xL, // punti iniziali in Linear
					sv, // SupportPointsForValues
					sd, // SupportPointsForDerivatives
					s, // SupportPoints for residuals
					deltaExteriorPoints,
					xc; // SupportPoints for contiguous points

	BzzMatrix	YL, // valori di Y iniziali in Linear
					Y, // valori di y nei numSupportPointsForValues

					// salvo il primo e l'ultimo sono i valori in xC
					YY, // valori di y nei numSupportPointsForDerivatives
					D1, // valori di y' nei numSupportPointsForDerivatives
					D2, // valori di y" nei numSupportPointsForDerivatives

					SY, // valori di y nei numSupportPoints
					SD1, // valori di y' nei numSupportPoints
					SD2; // valori di y" nei numSupportPoints

// credo sia meglio!!!!
// serve per UseMatrixRowAsVector
// per avere le y e le d1 di una variabile
//	BzzVector	*y, // riga di Y(numEquations,numSupportPointsForValues) // !!!!!!!!
//					yy, // riga di YY // non dovrebbe servire
//					*d1, // riga di D1
//					d2, // riga di D2 // non dovrebbe servire
// forse non così
	BzzVector	y, // riga di Y(numSupportPointsForValues,numEquations)
					yy, // riga di YY
					d1, // riga di D1
					d2, // riga di D2 // non dovrebbe servire

// credo sia meglio!!!!
// serve per UseMatrixRowAsVector
					sy, // riga di SY(numSupportPoints,numEquations); // !!!!!
					sd1, // riga di SD1
					sd2; // riga di SD2

	int initialFinal;
	int numVariablesNLS;
	BzzNonLinearSystemSparseObject nls;
	BzzVector z0,z; // Non Linear System variables

	BzzVectorInt	equationValuePresence, // ci sono o no y nei residui
						equationFirstDerivativePresence, // ci sono o no derivate prime
						equationSecondDerivativePresence, // ci sono o no derivate seconde
						leftBoundarySecondDerivativePresence, // ci sono o no derivate seconde
						rightBoundarySecondDerivativePresence, // ci sono o no derivate seconde
						numInternalPoints;

	int presenceOfSecondDerivativeInLeft; // 0 nessuna
	int presenceOfSecondDerivativeInRight; // 0 nessuna
	int constantParameter; // 1 se c'è 0 se non c'è

	BzzVector	criticalParametersSoft,
					criticalParametersRequired,
					criticalParameters,
					criticalParametersOK; // valore dei parametri per la soluzione trovata
	int criticalParOK; // 1 se giunge a soluzione per criticalParametersOK

	void(*ptrResidualsName)(double x, BzzVector &y,
	  BzzVector &d1, BzzVector &d2,BzzVector &cP,
		BzzVector &r);
	void(*ptrLeftBoundaryName)(double x, BzzVector &y,
	  BzzVector &d1, BzzVector &d2,BzzVector &cP,
		BzzVector &r);
	void(*ptrRightBoundaryName)(double x, BzzVector &y,
	  BzzVector &d1, BzzVector &d2,BzzVector &cP,
		BzzVector &r);

public:
//=========================< constructor >============================
	BzzBVP(void);
//==================< Inizializer >================================
	void operator()(int nE,int nV,
		BzzVectorInt &evp,BzzVectorInt &ed1p,BzzVectorInt &ed2p,
		BzzVectorInt &ld2p,BzzVectorInt &rd2p,
		void (*ptrRes)(double x,BzzVector &y,
			BzzVector &d1,BzzVector &d2,BzzVector &cP,
			BzzVector &r),
		void (*ptrLef)(double x,BzzVector &y,
			BzzVector &d1,BzzVector &d2,BzzVector &cP,
			BzzVector &r),
		void (*ptrRigh)(double x,BzzVector &y,
			BzzVector &d1,BzzVector &d2,BzzVector &cP,
			BzzVector &r),
		BzzVector &cPS,BzzVector &cPR,
		BzzVector &xx,BzzMatrix &YY,BzzVectorInt &nip,int initialFinal);

//	============================================================================
//	******************************< destructor >********************************
//	============================================================================
   ~BzzBVP(void);

//==================< BzzPrint >================================
	void ObjectBzzPrint(void){};
	};

//=============================================================================================
class MyNonLinearSystemObjectBVP : public BzzMyNonLinearSystemSparseObject
	{
friend class BzzBVP;
private:
	BzzBVP *ptrBVP;
public:
	MyNonLinearSystemObjectBVP(BzzBVP *ptr)
		{ptrBVP = ptr;}
	virtual void GetResiduals(BzzVector &x,BzzVector &f);
	virtual void ObjectBzzPrint(void);
	};

#endif // BZZ_BVP_HPP

*/