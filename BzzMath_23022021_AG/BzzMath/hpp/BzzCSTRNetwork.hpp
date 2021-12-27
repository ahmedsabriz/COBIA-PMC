// BZZMATH: Release 7.0

// =========================< BzzCSTRNetwork.HPP >=============================
// * Class BzzCSTRNetwork for *
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
/////////// Release 5.
//	11-2004	Date Written.

// ============================================================================
// ****** Constructors for BzzCSTRNetwork:												*

// ****************************************************************************
// ***** Access functions :																	*
// **
// ****************************************************************************
// ***** Assignment:																				*
// **
// ****************************************************************************
// ***** BzzPrint and BzzMessage																*
// * v.BzzPrint();																				*
// * v.BzzMessage();																				*
// **
// ****************************************************************************
// ***** Other functions:																		*
// **
// ****************************************************************************
// ****************************************************************************

#ifndef BZZ_CSTR_NETWORK_HPP
#define BZZ_CSTR_NETWORK_HPP

class BzzFactorizedDiagonalBlocksGaussAndMatrixLocked;

// ============================================================================
// ============================< class BzzCSTRNetwork >=============================
// ============================================================================

class BzzCSTRNetwork : public BzzBaseClass
{
	friend class BzzNonLinearSystemSparse;
	//friend class MyNonLinearSystemObjectCSTR;
	friend class MyOdeSystemObjectOneCSTR;
	friend class MyOdeSystemObjectAllCSTR;
private:
	static const double R_KIN;
	static const double R_CTOT;
	static const double UR_CTOT;
	static const double D1;
	static const double D2;
	static const double D3;
	static const double DSQ;

	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	// initialise constructors
//	void Initialize(int nc);

	int	whoAmI,
		printTasks,
		printSubTasks;

	int memoTemperature; // 1 salva su file 0 memorizza
	int memoNewton; // 1 salva su file 0 memorizza

	double tolRel, tolAbs;

	int	sameTemperature;
	int	numComponents, // num Rows nuij (200)
		numReactions,	// num columns nuij (3000)
		numEquil, // num reaction with equilibrium
		numCSTRReactors; // num reactors

	char** strReaction;
	char** strComponent;

	BzzVectorInt	numDir1,
		numDir2,
		numDir3,
		numDir4,
		numDir5,
		numInvTot1,
		numInvTot2,
		numInvTot3,
		numInvTot4,
		numInvTot5,
		numInvEq1,
		numInvEq2,
		numInvEq3,
		numInvEq4,
		numInvEq5,
		jDir1,
		jDir2,
		jDir3,
		jDir4,
		jDir5,
		jInvTot1,
		jInvTot2,
		jInvTot3,
		jInvTot4,
		jInvTot5,
		jInvEq1,
		jInvEq2,
		jInvEq3,
		jInvEq4,
		jInvEq5,
		jEquil,	// 0 se non c'�1 se c'�equilibrio
		reactionWithEquil, // indice j della reazione // dimensionato per jEquil=1
		jCorpoMorto, // 0,1,2,3..
		jNumEfficiency,
		* iEfficiency;
	BzzVector valDir5, valInvTot5, valInvEq5;

	BzzVector* efficiency;
	BzzVector	beta1, beta2,
		E1, E2,
		k01, k02,
		qanm,
		aFall, bFall, cFall, dFall, eFall,
		coeffM;

	BzzVector* aDH; // alta T
	BzzVector* bDH; // bassa T
	BzzVector* aDS; // alta T
	BzzVector* bDS; // bassa T
	// TODO a e b per derivate!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	BzzVector T1, T2, T3;

	BzzVector	componentDH,
		componentDS,
		reactionDH,
		reactionDS,
		sumNuij,
		molecularWeight;

	BzzMatrix	componentDHm,
		componentDSm;

	BzzVector	temperature,
		pression,
		volume,
		massRate,
		massInput,
		massOutput,
		sM,
		logTm, loguRTm, cTotm,
		feedInkRreactor;

	BzzMatrixSparse	Mg,
		Ms;

	BzzMatrixSparseLockedByRows	Ls,
		Ld;

	//	BzzFactorizedSparseGauss	Fg; // TODO Locked!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	BzzFactorizedSparseLockedByRowsGauss Fg;//,Ff;

	BzzMatrix	massFractionsInReactorsSolution, // reactors as rows
		massFractionsInReactors,
		volumeMolecularWeightT,
		feed,
		auxiliar;

	BzzVector	massFractionsInReactorsSolutionV,
		massFractionsInReactorsV,
		volumeMolecularWeightTV,
		feedV,
		auxiliarV;

	void SetSolution(BzzMatrix& XX)
	{
		massFractionsInReactorsSolution = XX;
	}

	//	BzzVectorInt cstrSequence,cstrSequenceAutoSwap,cstrExit,sumExit;
	BzzVectorInt cstrSequence, cstrSequenceAutoSwap, sumExit;
	BzzVectorInt cluster, abraCadabra;
	//	BzzVector exitInkReactor;
	BzzVectorIntArray cstrOut;
	BzzVectorIntArray cstrConnect;

	//	BzzMatrix reactionsRateInReactors; // reactors as rows

		// r[j] = rDirT[j] *(rDirC[j] + rInvC[j] / kEq[j]);
	BzzVector	r,
		k1, k2,
		rDirT,
		rDirC,
		//							kEq,
		uKeq,
		rInvC;
	BzzMatrix	uKeqm,
		k1m,
		k2m,
		logFcentm;

	// Memo for Derivatives
	// One point
//	double	wT,mT,
//				wP,mP;

	double T, P, logT, loguRT, cTot;
	double F1Stop, maxResOK, maxResOpt, resMeanOpt;
	int maxCountTutto;

	BzzVector	mc,
		mr,
		mrDirT,
		mrDirC,
		mR,
		mrInvC,
		wF,
		logFcent;//,mlogFcent;

//	void GetMolecularWeight(int i);
	void MemoTemperatureFunctions(void);
	BzzVector maRes, cRes, moRes, RRes;
	//	void ExpandCluster(BzzVectorInt &baseCluster,BzzVectorInt *newCluster);
	void ExpandCluster(BzzVectorInt& baseCluster, BzzVectorInt* newCluster,
		BzzMatrix& baseMassFractionsInReactors,
		BzzMatrix* newMassFractionsInReactors);
	//	void GenerateClusterStart(int numExpansion);
public:
	int	kReactor;
	//BzzVector feedInkRreactor;
	void ExpandCluster(int nExpansions);

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor BzzCSTRNetwork v;
	BzzCSTRNetwork(void);

	// copy-initializer
//	BzzCSTRNetwork(BzzCSTRNetwork &rval){BzzError("TODO)";}

	// from FILE
	BzzCSTRNetwork(const char* fileSt, const char* fileKin, const char* fileEq, const char* fileNetwork);
	void operator()(const char* fileSt, const char* fileKin, const char* fileEq, const char* fileNetwork);

	BzzCSTRNetwork(const char* fileNetwork, const char* first);
	void operator()(const char* fileNetwork, const char* first,
		int cicloCluster, int cicloDiffusion, int iaia, int relaxation);

	BzzCSTRNetwork(const char* fileSt, const char* fileKin, const char* fileEq);
	void operator()(const char* fileStconst, const char* fileKin, const char* fileEq);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzCSTRNetwork(void);

	// ============================================================================
	// *******************< Non-modifying access functions >***********************
	// ============================================================================
	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }
	void SetTolRel(double tolr) { tolRel = tolr; }
	void SetTolAbs(double tola) { tolAbs = tola; }
	void SetSumExit(BzzVectorInt& sumE) { sumExit = sumE; }
	void SetMemoTemperature(void) { memoTemperature = 1; } // 1 salva su file 0 memorizza
	void SetMemoNewton(void) { memoNewton = 1; } // 1 salva su file 0 memorizza
	int GetNumComponents(void) { return numComponents; }
	int GetNumReactions(void) { return numReactions; }
	int GetNumCSTRReactors(void) { return numCSTRReactors; }

	// ============================================================================
	// ******************************< Setting functions >*************************
	// ============================================================================
	void StopTaskPrint(void) { printTasks = 0; }
	void StopSubTasksPrint(void) { printSubTasks = 0; }
	void SetTasksPrint(void) { printTasks = 1; }
	void SetSubTasksPrint(int psb = 1)
	{
		if (psb > 0)
			printSubTasks = psb;
		else
			printSubTasks = 1;
	}

	// ============================================================================
	// **********************< Modifying access functions >************************
	// ============================================================================

	// ============================================================================
	// *************************< assignment operators >***************************
	// ============================================================================
	BzzCSTRNetwork& operator =
		(const BzzCSTRNetwork& rval);

	// ============================================================================
	// ====================< Non-modifying functions >=============================
	// ============================================================================

	// ******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);
	void OutputPrint(int ciclo);

	//	********************************< Save >************************************
	void Save(const char* file); // formatted
	void Save(char, const char* file);// binary

// ============================================================================
// ======================< Modifying Functions >===============================
// ============================================================================
	friend void Delete(BzzCSTRNetwork* result); // eliminates BzzCSTRNetwork

	// recovery from Save
//	friend void Load
//			(BzzCSTRNetwork *result,char *file***); // formatted
//	friend void Load
//			(BzzCSTRNetwork *result,char,char *file***);// binary

// ============================================================================
// =========================< Other functions >================================
// ============================================================================
	void operator()(BzzMatrix* massFractionsInReactors); // reactors as rows
	void operator()(BzzVectorInt& restart,
		BzzMatrix* massFractionsInReactors); // reactors as rows
	int GetFirst(void);
	void GetPreSecond(void);
	void GetSecond(void);
	int GetThird(void);
	//	void GetSecond(BzzMatrix *massFractionsInReactors);
	void GetReactionsRateInAllReactorsFromMassFractions
	(BzzMatrix& massFractionsInReactors, // reactors as rows
		BzzMatrix& reactionsRateInReactors); // reactors as rows
	void GetReactionsRateInAllReactorsFromMassFractions
	(BzzVector& massFractionsInReactorsV, // reactors as rows
		BzzVector& reactionsRateInReactors); // reactors as rows
	void GetReactionsRateFromConcentrations(BzzVector& c,
		BzzVector* R);
	void GetDerivativesC(BzzMatrix* dRC);
	void GetDiagonalMatricesForLinearizedSistem
	(BzzMatrix& massFractionsInReactors);
	void GetDiagonalMatricesForLinearizedSistem
	(char* file, BzzVector& mfV);

	//	void GetAllResiduals(BzzVector &m,BzzVector &residuals);
	void GetResiduals(BzzMatrix& massFractionsInReactors,
		BzzMatrix& residuals);
	void GetAllResiduals(BzzVector& m, BzzVector& residuals);
	void GetResiduals(BzzVector& m, BzzVector& residuals);
	void GetJacobian(BzzVector& x, BzzMatrix& JJ);

	void GetMolarFractionsFromMassFractions(BzzVector& molarFractions,
		BzzVector& massFractions, double wM);
	void GetMassFractionsFromMolarFractions(BzzVector& massFractions,
		BzzVector& molarFractions, double wM);
	void GetConcentrationsFromMolarFractions(BzzVector& concentrations,
		BzzVector& molarFractions, double cTot);
	void GetMolarFractionsFromConcentrations(BzzVector& molarFractions,
		BzzVector& concentrations, double cTot);

	double GetMeanWeightFromMolarFractions(BzzVector& molarFractions);
	double GetMeanWeightFromMassFractions(BzzVector& massFractions);
	void operator()(double T, double P, BzzVector& c, BzzVector* R);
	void SwapMassFractionsInReactors(BzzMatrix* massFractionsInReactors);
	//	void operator()(BzzVector pT,BzzVector pP,
	//		BzzMatrix &pc,BzzMatrix *pR);
	//	void ResetNumPoints(int nP);
	//	void SetBaseForDerivatives(void);

	//	void DerivativesT(BzzVector *dRT);
	//	void DerivativesT(int iP,BzzVector *dRT);

	//	void DerivativesC(int iP,BzzMatrix *dRC);
};

#endif // BZZ_CSTR_NETWORK_HPP