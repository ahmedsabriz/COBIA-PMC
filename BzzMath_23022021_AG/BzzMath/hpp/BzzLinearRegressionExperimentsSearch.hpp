// BZZMATH: Release 7.0

//	===================< BzzLinearRegressionExperimentsSearch.hpp >=============
//	* BzzLinearRegressionExperimentsSearch class for linear regression	best	*
//	*		experiments search																	*
// * Examples: c:\bzzmath\examples\BzzMathAdvabced\Regressions\					*
// *	LinearRegressionExperimentsSearch\LinearRegressionExperimentsSearch.cpp	*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	10-2003	Date Written.

////////////////// Release 5.0

//	============================================================================
//	******	BzzLinearRegressionExperimentsSearch constructors:						*
//	* BzzLinearRegressionExperimentsSearch r; 											*
//	* BzzLinearRegressionExperimentsSearch r(numP,FuntionName,xMin,xMax); 		*
//	* BzzLinearRegressionExperimentsSearch r(numP,FuntionName,Xnew); 				*
//	* BzzLinearRegressionExperimentsSearch r(Xbase,numP,FuntionName,xMin,xMax);*
//	* BzzLinearRegressionExperimentsSearch r(Xbase,numP,FuntionName,Xnew); 		*
//	****************************************************************************
//	***** Access functions:																		*
//	* who = r.WhoAmI();																			*
//	* int count = BzzLinearRegressionExperimentsSearch::ObjectCount();									*
//	* int countInScope = BzzLinearRegressionExperimentsSearch::ObjectCountInScope();				*
//	****************************************************************************
//	* Functions:																					*
//	* r(numP,FuntionName,xMin,xMax); 														*
//	* r(numP,FuntionName,Xnew); 																*
//	* r(Xbase,numP,FuntionName,xMin,xMax);													*
//	* r(Xbase,numP,FuntionName,Xnew); 														*
//	* r.SelectCriterion(n);																		*
//	* r();																							*
//	* r(numExperiments);																			*
//	* r.BzzPrint("Comment");																	*
//	****************************************************************************

#ifndef LINREG_EXPERIMENTS_DOUBLE_SEARCH_HPP
#define LINREG_EXPERIMENTS_DOUBLE_SEARCH_HPP

//	============================================================================
//	============< class BzzLinearRegressionExperimentsSearch >==================
//	============================================================================

class BzzLinearRegressionExperimentsSearch : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static const char* const WARNING;

	static int count; // for whoAmI
	static int countInScope;

	BzzFactorizedSVD svd;
	double percOf10;

	void (*ptrModel)(BzzVector& x, BzzVector& f);

	BzzMatrix	Fbase, // Original data
		Xbase,
		Fnew,
		Xnew,
		Ftotal,
		Xtotal,
		Xselected, // 8 x numX
		XexperimentsSelected, // numExperimentsSelected x numX
		FbaseP,
		FnewP,
		V, // SVD
		U,	// SVD
		P;	// Principal Component

	BzzVectorInt	iSelected, // 10
		iExperimentsSelected; // numExperimentsSelected

	BzzVector	b,	  // parametres
		xMin, xMax,
		x,
		f,
		d,	  // singular values
		f1, f2, f3, f4, f5, f6, f7,
		minimumDistanceX,
		minimumDistanceF,
		minimumDistanceSVD,
		minimumDistancexSVD,
		fCriterionBase,
		fCriterionOptimum,
		bTolerance,
		bzzColumnsCorrelation;

	int	numExperimentsBase,
		numExperimentsNew,
		numExperimentsTotal,
		numExperimentsSelected,
		numParameters,
		numX,
		defaultExperimentSelected,
		svdM1ExperimentSelected,
		svdM2ExperimentSelected,
		//			intercept, //  = 1 model with intercept
		whoAmI;

	int	selectedCriterion;
	// 0 default mixed
	// 1 Max det
	// 2 Max dMin
	// 3 Max fT(F1TF1)-1f
	// 4 Min Max fT(F1TF1)-1f
	// 5 Max sum d
	// 6 Max 1/cond
	// 7 Max det / sum
	// 8 Max Minimum distance in x space Criterion
	// 9 Max Minimum distance in f space Criterion
	// 10 Max Minimum distance in Principal Components x space Criterion
	// 11 Max Minimum distance in Principal Components f space Criterion

	double dMax, dMin;

	void Initialize(void);
	void DefaultSelection(void);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
	BzzLinearRegressionExperimentsSearch(void);
	BzzLinearRegressionExperimentsSearch
	(const BzzLinearRegressionExperimentsSearch& rval);

	//	* BzzLinearRegressionExperimentsSearch r(numP,FuntionName,xMin,xMax); 		*
	BzzLinearRegressionExperimentsSearch(int numP,
		void (*ptrM)(BzzVector& x, BzzVector& f),
		BzzVector& xMin, BzzVector& xMax);
	void operator()(int numP,
		void (*ptrM)(BzzVector& x, BzzVector& f),
		BzzVector& xMin, BzzVector& xMax);

	//	* BzzLinearRegressionExperimentsSearch r(numP,FuntionName,Xnew); 				*
	BzzLinearRegressionExperimentsSearch(int numP,
		void (*ptrM)(BzzVector& x, BzzVector& f),
		BzzMatrix& Xn);
	void operator()(int numP,
		void (*ptrM)(BzzVector& x, BzzVector& f),
		BzzMatrix& Xn);

	//	* BzzLinearRegressionExperimentsSearch r(Xbase,numP,FuntionName,xMin,xMax);*
	BzzLinearRegressionExperimentsSearch(int numP, BzzMatrix& Xb,
		void (*ptrM)(BzzVector& x, BzzVector& f),
		BzzVector& xMin, BzzVector& xMax);
	void operator()(int numP, BzzMatrix& Xb,
		void (*ptrM)(BzzVector& x, BzzVector& f),
		BzzVector& xMin, BzzVector& xMax);

	//	* BzzLinearRegressionExperimentsSearch r(Xbase,numP,FuntionName,Xnew); 		*
	BzzLinearRegressionExperimentsSearch(int numP, BzzMatrix& Xb,
		void (*ptrM)(BzzVector& x, BzzVector& f),
		BzzMatrix& Xn);
	void operator()(int numP, BzzMatrix& Xb,
		void (*ptrM)(BzzVector& x, BzzVector& f),
		BzzMatrix& Xn);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzLinearRegressionExperimentsSearch(void) { countInScope--; }

	// ============================================================================
	// *******************< Non-modifying access functions >***********************
	// ============================================================================
		// number of experiments
	//	int NumExperiments(void) const {return numExperiments;}

		// number of parameters
	int NumParameters(void) const { return numParameters; }

	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }
	void GetSelectedExperiments(BzzMatrix* newE)
	{
		*newE = XexperimentsSelected;
	}

	//	============================================================================
	//	*********************< Modifying access functions >*************************
	//	============================================================================
	void SelectCriterion(int s)
	{
		if (s < 0 || s  > 11)
			s = 10;
		selectedCriterion = s;
	}

	void SetPercentageOfCriterion10(double p = .0)
	{
		if (p < 0. || p > .8)
		{
			percOf10 = 0.8;
			selectedCriterion = 0;
		}
		else
		{
			if (p >= 0.1 && selectedCriterion == 4)
				selectedCriterion = 3;
			percOf10 = p;
		}
	}

	void operator()(void);

	//	============================================================================
	//	=======================< Non-modifying functions >==========================
	//	============================================================================
	virtual void ObjectBzzPrint(void);

	//	============================================================================
	//	===========================< Modifying Functions >==========================
	//	============================================================================
};

#endif // LINREG_EXPERIMENTS_DOUBLE_SEARCH_HPP