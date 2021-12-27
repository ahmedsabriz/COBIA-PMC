
// BZZMATH: Release 7.0

//	=============================< BzzOdeb.hpp >==========================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	01-2005	Date Written															 

#ifndef ODEB_DOUBLE_HPP
#define ODEB_DOUBLE_HPP

class BzzOdebObject : public BzzBaseClass
	{
public:
	virtual void GetNonLinearFunctionsValue(int el,double x,BzzVector &y,
		BzzVector &a,BzzVector &b,BzzVector &c,
		BzzVector &d,BzzVector &e,BzzVector &f,
		BzzVector &h,BzzVector &p,BzzVector &w){};
	virtual void GetNonLinearFunctionsValue(int el,double x,BzzVector &y,
		BzzVector &b,BzzVector &c,BzzVector &e){};
	virtual void GetInitialBoundaryNonLinearFunctionsValue(double x,BzzVector &y,
		BzzVector &s,BzzVector &m,BzzVector &n,BzzVector &t){};
	virtual void GetInitialBoundaryNonLinearFunctionsValue(double x,BzzVector &y,
		BzzVector &n){};
	virtual void GetFinalBoundaryNonLinearFunctionsValue(double x,BzzVector &y,
		BzzVector &s,BzzVector &m,BzzVector &n,BzzVector &t){};
	virtual void GetFinalBoundaryNonLinearFunctionsValue(double x,BzzVector &y,
		BzzVector &n){};
	virtual void ObjectBzzPrint(void){};
	};

class BzzOdeb : public BzzBaseClass
	{
private:
	enum BzzOdebType
		{
		QUASI_LINEARIZATION,
		ORTHOGONAL_COLLOCATION
		}odebType;

	BzzOdebObject *ptrOdeb;
	static const char *const BZZ_ERROR;
	static int MAX_T_ITERATION;
	static const double E1_REL;
	static const double E1_ABS;
	static const double E2_REL;
	static const double E2_ABS;
	static const double T_MAX;
	static const double T_START;

	char first;
	int	numVariables,
			numPointsForEachElement,
			numElements,
			numSupportPoints,
			numEquations,
			numIterations,
			iterTotal,
			iterT;
//			numtDiDj;

	double absTolerance,relTolerance,tSolution;

	int	numAlfa,
			numBeta,
			numGamma,
			numDelta,
			numEpsylon,
			numPhi,
			numEta,
			numPsi,
			numOmega,
			numSigmaI,numMuI,numNuI,numTauI,
			numSigmaF,numMuF,numNuF,numTauF;

	int eigenvalue;

	BzzVector	xElements,
							tSupportPointsInEachElement,
							xSupportPoints,
							y,yD,yDelta,
							yMin,yMax;

	int minSet,maxSet;

	BzzMatrix	Y0,
							Y,
							YSolution,
							Z0,
							Equation0,
							Equation,
							EquationSolution;

	BzzVector y1,y2,yy;

	BzzVector	alfa,
							beta,
							gamma,
							delta,
							epsylon,
							phi,
							eta,
							psi,
							omega,
							sigmaI,
							muI,
							nuI,
							tauI,
							sigmaF,
							muF,
							nuF,
							tauF,
							aux1ForDerivatives,
							aux2ForDerivatives;

	BzzMatrix	malfa,
							mbeta,
							mgamma,
							mdelta,
							mepsylon,
							mphi,
							meta,
							mpsi,
							momega,
							msigmaI,
							mmuI,
							mnuI,
							mtauI,
							msigmaF,
							mmuF,
							mnuF,
							mtauF;

	double uDeltaT;
	BzzMatrix	tMatrixForFirstDerivatives,
							tMatrixForSecondDerivatives,
							*xMatrixForFirstDerivatives,
							*xMatrixForSecondDerivatives,
							*yMatrixFirstDerivatives,
							*yMatrixSecondDerivatives;
//							tDiDj,
//							*xDiDj;

	BzzVectorInt nA,nB,nC,nD,nE,nF,nH,nP,nW;
	BzzVectorInt vB,vF;
//	BzzVectorInt nSi,nMi,nNi,nTai,nSf,nMf,nNf,nTaf;
//	BzzVectorIntArray dA,dB,dC,dD,dE,dF,dH,dP,dW;
//	BzzVectorIntArray dSi,dMi,dNi,dTai,dSf,dMf,dNf,dTaf;
//	BzzMatrixInt DSi,DMi,DNi,DTai,DSf,DMf,DNf,DTaf;

//	BzzVector	f,ga,gb;
//	BzzMatrix	Fy,Fy1,Fy2,Gay,Gay1,Gay2,Gby,Gby1,Gby2;

	BzzFactorizedStaircaseGauss staircaseGauss;
	BzzMatrixStaircase staircase;
	BzzVectorInt lowerBand,upperBand;
	BzzVector w,ww;

	void ObjectBzzPrint(void);
	void GetTSupportPointsInEachElement(void);
	void GetTMatricesForDerivatives(void);
	void GetXSupportPoints(void);
	void GetXMatricesForDerivatives(void);
	void GetAllResiduals(void);
	void GetElementResiduals(int element);
	int ConstraintsControl(void); // return 1 vincoli soddisfatti 0 violati

public:
//	BzzMatrixInt Da,Db,Dc,Dd,De,Dh,Dp,Dw;
//	BzzVector DA,DB,DC,DD,DE,DH,DP,DW;

//	============================================================================
//	*****************************< constructors >*******************************
//	============================================================================
	// default
	BzzOdeb(void);

	BzzOdeb(int numP,BzzVector &x);
	void operator()(int numP,BzzVector &x);

	BzzOdeb(BzzOdebObject *oo,int numP,BzzVector &x,
			BzzMatrix &Z,
			BzzVectorInt &na,BzzVectorInt &nb,BzzVectorInt &vb,BzzVectorInt &nc,
			BzzVectorInt &nd,BzzVectorInt &ne,
			BzzVectorInt &nf,BzzVectorInt &vf,
			BzzVectorInt &nh,BzzVectorInt &np,BzzVectorInt &nw,
			BzzVectorInt &nsi,BzzVectorInt &nmi,BzzVectorInt &nni,BzzVectorInt &ntai,
			BzzVectorInt &nsf,BzzVectorInt &nmf,BzzVectorInt &nnf,BzzVectorInt &ntaf);

	void operator()(BzzOdebObject *oo,int numP,BzzVector &x,
			BzzMatrix &Z,
			BzzVectorInt &na,BzzVectorInt &nb,BzzVectorInt &vb,BzzVectorInt &nc,
			BzzVectorInt &nd,BzzVectorInt &ne,
			BzzVectorInt &nf,BzzVectorInt &vf,
			BzzVectorInt &nh,BzzVectorInt &np,BzzVectorInt &nw,
			BzzVectorInt &nsi,BzzVectorInt &nmi,BzzVectorInt &nni,BzzVectorInt &ntai,
			BzzVectorInt &nsf,BzzVectorInt &nmf,BzzVectorInt &nnf,BzzVectorInt &ntaf);

// Quasi linearization
//	BzzOdeb(BzzOdebObject *oo,int numP,BzzVector &x,
//			BzzMatrix &Z);
//	void operator()(BzzOdebObject *oo,int numP,BzzVector &x,
//			BzzMatrix &Z);
	
	// return -1 NON CONVERGE
	// return 1 OK
	int operator()(int it = 5);

//	============================================================================
//	******************************< destructor >********************************
//	============================================================================
	void Deinitialize(void);
   ~BzzOdeb(void){BzzWarning("TODO destructor");};
   
//	============================================================================
//	==============================< Functions > ================================
//	============================================================================
	void GetMatrixForFirstDerivatives(BzzMatrix *D1);
	void SetMinimumConstraints(BzzVector &ymi);
	void SetMaximumConstraints(BzzVector &yma);
	void SetAbsTolerance(double tolA = 0.)
		{absTolerance = tolA;}
	void SetRelTolerance(double tolR = 1.e-5)
		{relTolerance = tolR;}
/*
	int GetNumStep(void){return numStep;}
	int GetNumFunction(void){return numFunction;}
	int GetNumFunctionForJacobian(void){return numFunctionForJacobian;}
	int GetNumAnalyticalJacobian(void){return numAnalyticalJacobian;}
	int GetNumNumericalJacobian(void){return numNumericalJacobian;}
	int GetNumFactorization(void){return numFactorization;}
	int GetNumSolution(void){return numSolution;}
	double GetHUsed(void){return hUsedInPreviousStep;}
	double GetHInNextStep(void){return hInNextStep;}
	int GetOrderUsed(void){return orderUsed;}
	int GetOrderInNextStep(void){return orderInNextStep;}
	int GetCalculationState(void){return odeCalculationState;}
	int GetComponentWithLargestError(void);
	BzzVector GetEstimatedErrors(void);
	BzzVector GetTollRel(void);
	BzzVector GetTollAbs(void);
	void GetInitAndEndTimeStep(double *tInitStep,double *tEndStep)
		{*tInitStep = t - hUsedInPreviousStep;*tEndStep = t;}
	double GetTimeInMeshPoint(void){return t;}
	BzzVector GetYInMeshPoint(void){return z[0];}
	BzzVector GetY1InMeshPoint(void){return f;}
	int GetOdeCalculationState(void){return odeCalculationState;}	
	void BzzPrintErrorState(void);
	void StepPrint(char *sv,char *pr = "Labels");
	void StepPrint(char *sv,BzzVectorInt &isv,char *pr = "Labels");
	void StepPrint(void (*stepPrintOut)(BzzVector &y,double t));

	void SetH0(double h00);
	void SetHMin(double hm);
	void SetHMax(double hm){hMax = hm;}
	void SetMaxStep(int maxS);
	void SetTolAbs(double tA);
	void SetTolAbs(BzzVector *tA);
	void SetTolAbs(const BzzVector &tA);
	void SetTolRel(double tR);
	void SetTolRel(BzzVector *tR);
	void SetTolRel(const BzzVector &tR);
	void SetMinimumConstraints(BzzVector *yMi);
	void SetMinimumConstraints(BzzVector &yMi);	
	void SetMaximumConstraints(BzzVector *yMa);	
	void SetMaximumConstraints(BzzVector &yMa);		

	void StopTaskPrint(void){printTasks = 0;}
	void StopSubTasksPrint(void){printSubTasks = 0;}
	void SetTasksPrint(void){printTasks = 1;}
	void SetSubTasksPrint(int psb = 1)
		{
		if(psb > 0)
			printSubTasks = psb;
		else 
			printSubTasks = 1;
		}
	void StopIntegrationBeforeRecalcuatingJacobian(int nJ = 1)
		{stopIntegrationBeforeRecalcuatingJacobian = nJ;}
	BzzVector operator() (double tF);
	BzzVector Operator(double tF){return (*this)(tF);}
	BzzVector operator() (double tF,double tC);
	BzzVector Operator(double tF,double tC){return (*this)(tF,tC);}
*/
	};


#endif // ODEB_DOUBLE_HPP
