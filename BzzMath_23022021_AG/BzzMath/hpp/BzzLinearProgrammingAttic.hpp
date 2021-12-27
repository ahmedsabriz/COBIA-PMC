// BZZMATH: Release 7.0

//	==============< BzzLinearProgrammingAttic.HPP >=============================
//	* BzzLinearProgrammingAttic: Class for linear programming						*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	07-2000	Date Written.
//	10-2009	Revision.
//	11-2010	New version

//	============================================================================
//	* This class solve the following linear programming								*
//	* Min F = sTx																					*
//	* Ex = e	(mE linear equality)																*
//	* Dx <= c (mD linear inequality)														*
//	****************************************************************************
//	============================================================================
//	******* Constructor																			*
//	* BzzLinearProgrammingAttic lpa(&x0,&s,numE,numElementsE,&E,&e,				*
			//numD,numElementsD,&D,&c);														*
//	* BzzLinearProgrammingAttic lpa;															*
//	* lpa(&x0,&s,numE,numElementsE,&E,&e,numD,numElementsD,&D,&c);					*
//	* BzzLinearProgrammingAttic lpa("FILELP.ATC");										*
//	****************************************************************************
//	******* Functions																				*
// * lp(n);																							*
//	****************************************************************************

#ifndef BZZ_LINEAR_PROGRAMMING_ATTIC
#define BZZ_LINEAR_PROGRAMMING_ATTIC

// ============================================================================
// =========================< class BzzMonoAttic >=============================
// ============================================================================

class BzzMonoAttic : public BzzBaseClass
{
private:
	enum BzzMonoAttiSolutionType
	{
		NOT_READY = 0,
		NORMAL = 1,
		ROOF = 2,
		FLOOR = 3,
		LOWER = 4,
		UPPER = 5,
		UNFEASIBLE_LOWER = 6,
		UNFEASIBLE_UPPER = 7
	}monoAttiSolutionType;

	static const char* const BZZ_ERROR;
	static const double BZZ_BIG_ATTIC;
	static int count; // for whoAmI
	static int countInScope;
	int whoAmI;
	int i, j, k, itMin, ilowup;
	int numVariables;
	int noBestWorstVertex;
	int workingVariable, vertex, numWorkingVariable;
	double ed, t, tMin, bEq, bigAttic, f;
	BzzVectorInt activeBound;
	BzzVectorInt satisfiedBound;
	BzzVector x, lambda, eq, s, lowerBounds, upperBounds;
	BzzVector xFloor, xRoof, xLower, xUpper, dxUL, d;
	double bLower, bUpper, bFloor, bRoof,
		fLower, fUpper, fFloor, fRoof;
	BzzVector sorted;
	BzzVectorInt iSorted;

public:
	void SetBigAttic(double big)
	{
		bigAttic = fabs(big);
	}
	int WhoAmI(void) { return whoAmI; }
	BzzMonoAttic(void);
	void operator()(BzzVector* cc, BzzVector* ss,
		BzzVector* xxL, BzzVector* xxU);
	int operator()(double b);
	virtual void ObjectBzzPrint(void);
};

//	============================================================================
//	=============< class BzzVariablesInOneEqualityOnly >========================
//	============================================================================
class BzzVariablesInOneEqualityOnly : public BzzBaseClass
{
private:
	int variable;
	int equality;
	int numVariables;
	double coefficient;
	BzzVectorInt iVariable;
	BzzVector coefVariable;
public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzVariablesInOneEqualityOnly(void) {}
	~BzzVariablesInOneEqualityOnly(void) { BzzError("URCA"); }
	void InitializeObject(int v, int e, double val, BzzVectorInt* iv, BzzVector* vv);
	friend double TProduct(BzzVariablesInOneEqualityOnly& vioeo, BzzVector& xd);

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
};
//	============================================================================
//	=================< class BzzLinearProgrammingAttic >========================
//	============================================================================
class BzzLinearProgrammingAttic : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // per whoAmI
	static int countInScope;
	static const double EPS_CONSTRAINTS;
	static const double EPS_COSINE_CONSTRAINTS;
	static const int MAX_FACTORIZATION_BEFORE_DELETE;
	static const int MIN_NEGATIVE_LAMBDA_BEFORE_DELETE;
	static const double BZZ_BIG_ATTIC;

	///////////////////////////////////////////////////////////////////////////
	// Per fattorizzazione e soluzione
	int numRows, numColumns, numElements, marameo;
	int numLeft, numK, numRight;
	BzzVector diag, auxd;
	BzzVectorArray vdaR, vdaL;
	BzzVectorIntArray viaR, viaL;
	BzzVectorInt ip, ipp, ipSingular, jpSingular, numR, jp, numL;//,iCol;//,dimC;
	BzzVectorInt ipl, auxi, auxiLeft, auxijDecode;
	int* ptri;
	double* ptrd;

	void InitializeGauss(void);
	int FactorizationWithPivoting(int k, int pivot = 0);
	void SolutionWithPivoting(BzzVector* bx);
	void TransposeSolutionWithPivoting(BzzVector* bx);

	///////////////////////////////////////////////////////////////////////////

	int whoAmI;
	FILE* fileLinearProgramming;
	FILE* fileBzzTmp;

	int	verbosity, subTaskVerbosity,
		versionSelected,
		printEqualityEquationsAnalysis,
		printStartingEqualityEquationsAnalysis,
		printEqualitySingletons,
		printEqualityVariablesAnalysis,
		printStartingEqualityVariablesAnalysis,
		printInequalityEquationsAnalysis,
		printStartingInequalityEquationsAnalysis,
		printInequalityVariablesAnalysis,
		printStartingInequalityVariablesAnalysis,
		printBounds,
		printSatisfiedEquations,
		printRedundancy,
		plotEqualityEquations,
		plotEqualityVariables,
		plotInequalityEquations,
		plotInequalityVariables,
		printWarning,
		printPartialTimes,
		printEqualityInsertion,
		printInequalityInsertion,
		printObjectFunctionAnalysis,
		getEqualityInequalityInsertion, // 0 solo inequality 1 anche le equality
		getInequalitySingletonInsertion,
		getInequalityReinsertion;

	int	numVariables, numVariablesOriginal,
		numEqualityConstraints, numEqualityConstraintsOriginal, // num Equality
		numElementsE, numElementsEOriginal,
		numInequalityConstraints, numInequalityConstraintsOriginal,  // num Inequality
		numDeletedEqualityConstraints, numDeletedInequalityConstraints, numDeletedNonNegativityConstraints,
		numStartingNonNegativityConstraints,
		numElementsD, numElementsDOriginal,
		numFreeVariables, numFreeVariablesOriginal,
		numEqualitySingleton, numEqualityDoublet, numEqualityTriplet,
		numEqualityLargerThanNumberEquality,
		numInequalitySingleton, numInequalityDoublet, numInequalityTriplet,
		numSingletonInEqualityMatrixColumns,
		numSingletonInInequalityMatrixColumns,
		numVariablesNonPresentInEquality,
		numVariablesNonPresentInEqualityWithNullRHS,
		numVariablesNonPresentInInequality,
		numVariablesNonPresentInEqualityAndInInequality,
		numVariablesNonPresentInEqualityAndSingletonInInequality,
		numVariablesNonPresentInInequalityAndSingletonInEquality,
		numInequalityEqualitySingleton,
		numInequalityLowerBound,
		numInequalityUpperBound,
		numNonNegativityConstraints, // numVariables - numFreeVariables
		numNonNegativityConstraintsOriginal,
		numActiveConstraints,
		numNullRHSTotal, numNullE, numNullC,
		numVariablesInEqualityConstraints,
		numVariablesNonInEqualityConstraints,
		numSatisfiedInequalityConstraints,
		numPassiveInequalityConstraints,
		numUnsatisfiedInequalityConstraints,
		numSatisfiedNonNegativityConstraints,
		numPassiveNonNegativityConstraints,
		numUnsatisfiedNonNegativityConstraints,
		numSatisfiedEqualityConstraints,
		numUnsatisfiedEqualityConstraintsPlus,
		numUnsatisfiedEqualityConstraintsMinus,
		numInequalitySingletonInsertion,
		numFloorPossibleConstraints,
		numRoofPossibleConstraints,
		numRoofPossibleVariables,
		// floor- 3 tutti i coefficienti sono - e s <= 0 si può aumentare x e la funzione non peggiora
		numVariablesWithFloorConstraintsMinus,
		// floor+ 4 tutti i coefficienti sono + e s > 0 si può diminuire x e la funzione non peggiora
		numVariablesWithFloorConstraintsPlus,
		// roof- tutti i coefficienti sono - e s > 0 si può aumentare x ma la funzione peggiora
		numVariablesWithRoofConstraintsMinus,
		// roof+ 1 tutti i coefficienti sono + e s <= 0 si può diminuire x ma la funzione peggiora
		numVariablesWithRoofConstraintsPlus;

	int	primalDual,
		artificial, // 0 normale 1 artificiale
		imaxS, imaxE, //Singleton Equality
		iminE, iminI, iminN, // Equality,Inequality,NonNegativity
		imaxI, imaxN,
		iter, iterTotal, maxIterations,
		iterF, iterTotalF, maxFactorization;

	int	setDualSelection,// = 1 permette di selezionare il duale se meglio; 0 impone il primal; 2 impone il dual
		getBalance, // 1 fa il bilanciamento 0 Non lo fa
		getNormalization, // 1 fa la normalizzazione 0 non lo fa
		getFloorVariablesSearch,
		floorSuccessful, // 0 no floor or unsuccess 1 success
		getLinearCombinations, // 1 verifica le combinazioni lineari con LQ
		// 0 non fa la prevenzione
		// 1 previene la degeneracy default
		getDegeneracyPrevenction,
		numFactorizedEquality,
		numFactorizedLocked,
		numFactorizedHyperLocked, // numEqualitySingleton + numInequalityEqualitySingleton
		numFactorizedInequality,
		numFactorizedNonNegativity,
		numFactorizedRows,
		linearDependent,
		numReinsertInequalityNonSingleton,
		numConstraintsJ,
		numLinearlyDependentEqualityConstraints;

	double startTime, endTime, totalTime, elapsedTime, partialTime;

	double	f,
		f0,
		fMin,
		fF, fFO,
		fE,
		fA,
		epsCosineConstraints,
		tMin;

	int feasibleF, feasibleE, // 1 se feasible; 0 se non c'è; -1 se infeasible
		feasibleFounded; // 0 se non ha ancora trovato alcun punto feasible; 1 se è passato da un feasible
	int ifloorPlus, ifloorPlusNonN, ifloorMinus, iminIfloor, iminNfloor;
	int solutionFounded; // 1 founded -1 infeasible solution -2 unbounded solution
	int setPrintFile;
	int numPositiveOF, numNegativeOF, numNullOF;

	BzzVector	x0,
		x,
		w, w0,
		xMin,
		xF, xFO,
		xE,
		xFloorPlus, xFloorMinus,
		s, sOriginal, sArtificial, sa, sSorted, saSorted,
		vE, vEOriginal,
		e, eOriginal,
		vD, vDOriginal,
		c, cOriginal,
		h, d, lambda,
		pE, pD, pN,
		qE, qD, qN,
		norm;

	int numOutliersInequality, numOutliersEquality;
	BzzVectorInt iOutliersInequality, iOutliersEquality;
	BzzVectorInt outliersInequality, outliersEquality;
	BzzVectorInt	variablesNonPresentInEquality,
		variablesNonPresentInInequality,
		variablesNonPresentInEqualityAndInInequality,
		nonNegativityForVariablesNonPresentInEqualityAndInInequality,
		variablesNonPresentInEqualityWithNullRHS,
		variablesNonPresentInEqualityAndSingletonInInequality,
		variablesNonPresentInInequalityAndSingletonInEquality;
	int numOutliersInequalityVariables, numOutliersEqualityVariables;
	BzzVectorInt iOutliersInequalityVariables, iOutliersEqualityVariables;
	BzzVectorInt outliersInequalityVariables, outliersEqualityVariables;

	BzzVectorInt numColumnsJ, insertJ;
	BzzVectorIntArray rowJ, colJ;

	BzzVectorInt	rE, rEOriginal,
		cE, cEOriginal,
		rD, rDOriginal,
		cD, cDOriginal,
		iw, iw0,
		objectFunctionCoefficients, // i primi sono negativi, poi positivi e poi nulli
		variablesInEqualityConstraints,
		equalitySingleton,
		equalityDoublet,
		variableSingleton,
		inequalitySingleton,
		inequalityDoublet,
		inequalityEqualitySingleton,
		variableInequalityEqualitySingleton,
		singletonInEqualityMatrixColumns,
		singletonInInequalityMatrixColumns,
		inequalityLowerBound,
		inequalityUpperBound,
		variableWithLowerBound,
		variableWithUpperBound,
		freeVariables, freeVariablesOriginal,
		variablesWithNonNegativityConstraints, variablesWithNonNegativityConstraintsOriginal,
		numVariablesInEachInequalityConstraint,
		numInequalityConstraintsForEachVariable,
		numVariablesInEachEqualityConstraint,
		numEqualityConstraintsForEachVariable,
		numEqualityAndInequalityConstraintsForEachVariable,
		numEqualityConstraintsWithNullRHSForEachVariable,
		satisfiedEqualityConstraints,// 1, sat; 2 unsatP; 3 unsatM
		satisfiedInequalityConstraints, // 1 sat; 2 unsat; 3 passive
		satisfiedNonNegativityConstraints,//  // 1 sat; 2 unsat; 3 passive
		floorRoofPossibleConstraints,
		roofPossibleVariables,
		floorRoofVariables,  // -1 NO 1 roof+ 2 roof- 3 floor- 4 floor+
		variablesWithRoofConstraintsPlus,
		variablesWithRoofConstraintsMinus,
		variablesWithFloorConstraintsMinus,
		variablesWithFloorConstraintsPlus,
		activeEquality,
		activeInequality,
		activeNonNegativity,
		factorizedRows, // numVariables  1 equality 2 inequality 3 non negativity constraint
		iFactorizedRows, // numVariables indice del vincolo nella sua categoria
		linearCombinations, // in equality constraints
		factorizedEquality,
		factorizedInequality,
		factorizedNonNegativity,
		reinsertInequality,
		exchangeE, exchangeD,
		isSorted, isaSorted,
		aux, aux1, aux2;

	// 0 free
	// 1 non negativity and s < 0
	// 2 non negativity and s > 0
	// 3 non negativity and s = 0
	// + 10 satisfied
	// + 20 violated
	// + 30 passive
	// 102 inequality lowerbound
	// 101 inequality upperbound
	// Es 32 nonnegativity passive and s > 0
	////////////////////////////////////////
	// 0
	// 1 non negativity satisfied, s > 0
	// 2 lower bound satisfied, s > 0
	// 3 upper bound satisfied, s < 0
	BzzVectorInt variablesCharacteristics;

	// 0 unbounded
	// 1 lower only
	// 2 upper only
	// 3 different lower and upper
	// 4 singleton
	BzzVectorInt variableBoundCharacteristics;

	BzzVectorIntArray rowD;
	BzzVectorArray valRowD;
	BzzVectorIntArray colD;
	BzzVectorArray valColD;

	BzzVectorIntArray rowE;
	BzzVectorArray valRowE;
	BzzVectorIntArray rowESorted;
	BzzVectorIntArray rowEiSorted;
	BzzVectorArray valRowESorted;

	BzzVectorIntArray colE;
	BzzVectorArray valColE;

	BzzVector	equalityRowsBalance,
		inequalityRowsBalance,
		variablesBalance;

	BzzVector lowerBound, upperBound; // dimensionato numVariables
	int numLower, numUpper; // nuove varaibili
	BzzVector lower, upper; // dimensionati numLower,numUpper
	BzzVectorInt iLower, iUpper; // dimensionati numLower,numUpper
	int numSingletons;
	BzzVector singletons;
	BzzVectorInt iSingletons;
	double fSingletons; // valore della funzione con i singletons inseriti

	void Initialize(void);
	void InitializePrint(void);
	void Balance(void);
	void VariablesInEqualityConstraints(void);
	void VariablesInEachInequalityConstraint(void);
	void InequalityConstraintsForEachVariable(void);
	void OutliersInequalityVariables(void);
	void FindBounds(void);
	void FindRedundancy(void);

	void VariablesInEachEqualityConstraint(void);
	void EqualityConstraintsForEachVariable(void);
	void OutliersEqualityVariables(void);
	void CheckEqualityConstraints(void);
	void CheckInequalityConstraints(void);
	void BuildDualProblem(void);
	void FloorRoofPossibleConstraints(void);
	void VariablesWithFloorConstraints(void);
	void BuildArtificialFunction(void);
	void BuildArtificialFunctionForEquality(void); // TODO
//	void GetNewEqualityConstraintAndX(void);
	void GetNewInequalityConstraintAndX(void);
	void GetNewEqualityInequalityConstraintAndX(void); // TODO
	void GetNewInequalityConstraintAndXArtificial(void);
	void GetNewEqualityInequalityConstraintAndXArtificial(void); // TODO
	void InsertSingleton(void);
	void InsertEquality(void);
	void InequalitySingletonInsertionNumber(void);
	void InequalitySingletonInsertion(void);
	void CheckFeasibilityForInequalityConstraints(void);
	int FindFeasiblePoint(void);
	//	void BuildMatrixJ(void);
	void FirstStep(void);
	void GetObjectFunctionCoefficients(void);
	void GetLinearCombinations(void);
	void EqualityFactorizationAndKKTSolution(void);
	void EqualityFactorization(void);
	void EqualityInsertion(void);
	void InequalityFactorization(void);
	//void ArtificialFactorization(void);
	void KKTSolution(void);
	void EqualityResiduals(void);
	void SatisfiedEqualityNumber(void);
	void InequalityResiduals(void);
	void SatisfiedInequalityNumber(void);
	void AtticFeasibleVertex(int vero);
	void AtticFeasibleNonVertex(void);
	void AtticInfeasibleVertex(int vero);
	void AtticInfeasibleNonVertex(void);
	void DeleteFactorization(int row);
	void DeleteFactorization(void);
	void InsertFactorization(void);
	void VariablesCharacteristics(void);
	void CompleteInequalitySingletonInsertion(void);
	void InequalityReinsertion(void);
	void CheckSingleton(void);
	void ModifyStartingPointSingletonCheck(void);
	void ModifyStartingPointHyperSingletonCheck(void);
	////////////////////////////////////////////////////////
		// NEW FUNCTIONS and Variables
	double bigAttic;
	int assignedBounds;
	int	numXaIVariables, // variabili ausiliarie o artificiali per le Inequality
		numXaEVariables; // variabili artificiali per le Equality
//iSortedNumVariablesInEquality
	BzzVectorInt sortedNumVariablesInEachEqualityConstraint,
		sortedVariableInEachEqualityConstraint;
	BzzVectorInt	sortedNumEqualityConstraintsForEachVariable,
		sortedEqualityConstraintForEachVariable;
	BzzVectorInt sortedNumVariablesInEachInequalityConstraint,
		sortedVariableInEachInequalityConstraint;
	BzzVectorInt	sortedNumInequalityConstraintsForEachVariable,
		sortedInequalityConstraintForEachVariable;
	BzzVectorInt	sortedNumEqualityAndInequalityConstraintForEachVariable,
		sortedEqualityAndInequalityConstraintForEachVariable;

	BzzVector sortedDeltaBounds;
	BzzVectorInt iSortedDeltaBounds;

	int	numEqualityWithNullRHS,
		numEqualityWithPositiveRHS,
		numEqualityWithNegativeRHS;
	BzzVectorInt	equalityWithNullRHS,
		equalityWithPositiveRHS,
		equalityWithNegativeRHS;
	double sArtificialValue;

	int numTriangularBlocks;
	BzzVectorInt rowsInTriangularBlocks;
	BzzVectorInt columnsInTriangularBlocks;

	BzzVector	xaI,
		saI,
		lowerXaI,
		upperXaI,
		xaE,
		saE,
		lowerXaE,
		upperXaE;

	int	numXE1Variables, // numero di variabili che compaiono solo in una equality
		numXI1Variables, // numero di variabili che compaiono solo in una inequality
		numEqualityWithxE1, // numero di equality che hanno almeno unsingleton xE1
		numInequalityWithxI1; // numero di inequality che hanno almeno unsingleton xI1

	BzzVectorInt	iXE1Variables, // quali sono le variabili che compaiono solo in una equality
		iXI1Variables; // quali sono le variabili che compaiono solo in una inequality

// per ogni equality indica il numero di variabili che compaiono solo in quella equality
// es in equality 1 ci sono 3 variabili che...
	BzzVectorInt numVariablesInOneEqualityOnly; // dimensionato numEqualityConstraints

	// per ogni equality indica le variabili che compaiono solo in quella equality
	// es. nella equality 1 le tre variabili sono 25, 26, 27
	BzzVectorIntArray variablesInOneEqualityOnly; // dimensionato numVariablesInOneEqualityOnly
	BzzVectorArray coefForVariablesInOneEqualityOnly;
	BzzVector yE1Lower, yE1Upper; // bound per l'insieme delle vInOneEqOnly

	BzzVectorInt numVariablesInOneInequalityOnly; // dimensionato numInequalityConstraints
	BzzVectorIntArray variablesInOneInequalityOnly; // dimensionato numVariablesInOneInequalityOnly

	// vettore di oggetti dimensionato come il numero di variabili
	// che compaiono in una sola equality
	// per ogni oggetto del vettorevariabile che compare in una sola equality
	BzzVariablesInOneEqualityOnly* vInE;

	// forse non serve basta il segno
		// -1 se artificiale
		// 0 se vincolo attivo
		// se vincolo passivo
	//	BzzVectorInt	ixaI,
	//		iaxE;

	BzzVector yPlot;
	void BzzPlot(const char* str, BzzVector& yPlot);
	void Preprocessing(void);
	void SetDimensions(void);
	void BuildRCVEfromRow(void);
	void BuildRCVEfromCol(void);
	void BuildRowEValRowE(void);
	void BuildColEValColE(void);

	void BuildRCVDfromRow(void);
	void BuildRCVDfromCol(void);
	void BuildRowDValRowD(void);
	void BuildColDValColD(void);

	void DeleteEquality(int eq);
	void DeleteEqualities(BzzVectorInt& eq);
	void DeleteVariableInEquality(int var);
	void DeleteVariablesInEquality(BzzVectorInt& var);

	void DeleteInequality(int eq);
	void DeleteInequalities(BzzVectorInt& eq);
	void DeleteVariableInInequality(int var);
	void DeleteVariablesInInequality(BzzVectorInt& var);
	void CheckBounds(void);
	void DeleteRedundancy(void);

	void InitializeBounds(void);
	void FindSingletons(void);

	void SelectStartingPoint(void);
	void CheckResiduals(void);
	int FloorVariablesSearch(void);
	//	void AuxiliaryVariables(void);
	void ProblemAnalysis(void);
	void FindVariablesInOneEqualityOnly(void);
	void FindVariablesInOneInequalityOnly(void);
	void FindVariablesNonPresentInEqualityAndInInequality(void);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzLinearProgrammingAttic(void);

	// copy constructor
//	BzzLinearProgrammingAttic(const BzzLinearProgrammingAttic &rval);

	// constructor
	void operator()(BzzVector* xx, BzzVector* ss,
		int numE, int numElE, BzzVectorInt* rrE, BzzVectorInt* ccE, BzzVector* vvE, BzzVector* ee,
		int numD, int numElD, BzzVectorInt* rrD, BzzVectorInt* ccD, BzzVector* vvD, BzzVector* cc,
		BzzVectorInt* ffr);
	void operator()(BzzVector* xx, BzzVector* xxL, BzzVector* xxU, BzzVector* ss,
		int numE, int numElE, BzzVectorInt* rrE, BzzVectorInt* ccE, BzzVector* vvE, BzzVector* ee,
		int numD, int numElD, BzzVectorInt* rrD, BzzVectorInt* ccD, BzzVector* vvD, BzzVector* cc);

	// constructor
	void operator()(BzzVector* xx, BzzVector* ss,
		BzzMatrix* EE, BzzVector* ee, BzzMatrix* DD, BzzVector* cc,
		BzzVectorInt* ffr);
	// con bound
	void operator()(BzzVector* xx,
		BzzVector* xxL, BzzVector* xxU, BzzVector* ss,
		BzzMatrix* EE, BzzVector* ee, BzzMatrix* DD, BzzVector* cc);

	// constructor FILE .ATC
	void operator()(char* fileLP);

	int operator()(int maxI = 0, int maxF = 0);
	// return 0 non ancora
	// TODO sceglire

//	============================================================================
//	********************************< destructor >******************************
//	============================================================================
	~BzzLinearProgrammingAttic(void)
	{
		countInScope--;
		if (numEqualityConstraints != 0)
		{
//#if BZZ_COMPILER != 0
			delete[] vInE;
//#else
//			delete[numEqualityConstraints + 1] vInE;
//#endif
		}

		if (setPrintFile == 1)
		{
			fclose(fileLinearProgramming);
			bzzFileOut = fileBzzTmp;
		}
	}

	//	============================================================================
	//	*****************************< Access functions >***************************
	//	============================================================================

	int WhoAmI(void) const { return whoAmI; }
	double GetSolution(BzzVector* xs)
	{
		*xs = x;
		return s % x;
	}

	//	============================================================================
	//	***************************< assignment operators >*************************
	//	============================================================================

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================
	void SetBigAttic(double big)
	{
		bigAttic = fabs(big);
	}

	void SetNoInequalityReinsertion(void)
	{
		getInequalityReinsertion = 0;
	}
	void GetEqualityInequalityInsertion(void)
	{
		getEqualityInequalityInsertion = 1;
	}
	void SetNoInequalitySingletonInsertion(void)
	{
		getInequalitySingletonInsertion = 0;
	}

	void PrintObjectFunctionAnalysis(int i);
	void PrintObjectFunctionAnalysis(void)
	{
		printObjectFunctionAnalysis = 1;
	}

	void PrintEqualityInsertion(void)
	{
		printEqualityInsertion = 1;
	}

	void PrintInequalityInsertion(void)
	{
		printInequalityInsertion = 1;
	}

	void PrintWarning(void)
	{
		printWarning = 1;
	}

	void PrintPartialTimes(void)
	{
		printPartialTimes = 1;
	}

	void PrintStartingEqualityEquationsAnalysis(void)
	{
		printStartingEqualityEquationsAnalysis = 1;
	}

	void PrintEqualityEquationsAnalysis(int dummy);
	void PrintEqualityEquationsAnalysis(void)
	{
		printEqualityEquationsAnalysis = 1;
	}

	void PrintEqualitySingletons(void)
	{
		printEqualitySingletons = 1;
	}

	void PrintStartingInequalityEquationsAnalysis(void)
	{
		printStartingInequalityEquationsAnalysis = 1;
	}

	void PrintInequalityEquationsAnalysis(int dummy);
	void PrintInequalityEquationsAnalysis(void)
	{
		printInequalityEquationsAnalysis = 1;
	}

	void PrintBounds(void)
	{
		printBounds = 1;
	}

	void PrintSatisfiedEquations(void)
	{
		printSatisfiedEquations = 1;
	}

	void PrintRedundancy(void)
	{
		printRedundancy = 1;
	}

	void PlotEqualityEquationsStructure(void)
	{
		plotEqualityEquations = 1;
	}
	void PlotInequalityEquationsStructure(void)
	{
		plotInequalityEquations = 1;
	}
	void PlotEqualityVariablesStructure(void)
	{
		plotEqualityVariables = 1;
	}
	void PlotInequalityVariablesStructure(void)
	{
		plotInequalityVariables = 1;
	}

	void PrintStartingEqualityVariablesAnalysis(void)
	{
		printStartingEqualityVariablesAnalysis = 1;
	}

	void PrintEqualityAndInequalityVariablesAnalysis(int dummy);
	void PrintEqualityVariablesAnalysis(int dummy);
	void PrintEqualityVariablesAnalysis(void)
	{
		printEqualityVariablesAnalysis = 1;
	}

	void PrintStartingInequalityVariablesAnalysis(void)
	{
		printStartingInequalityVariablesAnalysis = 1;
	}

	void PrintInequalityVariablesAnalysis(int dummy);
	void PrintInequalityVariablesAnalysis(void)
	{
		printInequalityVariablesAnalysis = 1;
	}

	void EqualityIndipendenceControl(void)
	{
		getLinearCombinations = 1;
	}

	void SetPrintFile(char* ptrF)
	{
		setPrintFile = 1;
		if ((fileLinearProgramming = fopen(ptrF, "w")) == NULL)
			BzzError("I can't open %s", ptrF);
		fileBzzTmp = bzzFileOut;
		bzzFileOut = fileLinearProgramming;
	}

	void SetCosinePrecision(double eps)
	{
		if (eps > 1.e-12 && eps < 1.e-3)
			epsCosineConstraints = eps;
	}

	void SetSubTaskVerbosity(int verb = 0)
	{
		if (verb >= 0)
			subTaskVerbosity = verb;
		else
			subTaskVerbosity = 0;
	}

	void SetTaskVerbosity(int verb = 0)
	{
		if (verb >= 0)
			verbosity = verb;
		else
			verbosity = 1;
		switch (verbosity)
		{
		case 0:
			break;
		case 6:
			PlotEqualityEquationsStructure();
			PlotEqualityVariablesStructure();
			PlotInequalityEquationsStructure();
			PlotInequalityVariablesStructure();
		case 5:
			PrintStartingEqualityEquationsAnalysis();
			PrintStartingInequalityEquationsAnalysis();
			PrintStartingEqualityVariablesAnalysis();
			PrintStartingInequalityVariablesAnalysis();
		case 4:
			PrintEqualitySingletons();
			PrintEqualityInsertion();
			PrintInequalityInsertion();
			PrintWarning();
			PrintBounds();
		case 3:
			PrintEqualityEquationsAnalysis();
			PrintEqualityVariablesAnalysis();
			PrintInequalityEquationsAnalysis();
			PrintInequalityVariablesAnalysis();
			PrintSatisfiedEquations();
		case 2:
			PrintPartialTimes();
		case 1:
			PrintObjectFunctionAnalysis();
			break;
		default:
			break;
		}
	}

	void SetPrimalProblem(void)
	{
		setDualSelection = 0;
	}

	void SetDualProblem(void)
	{
		setDualSelection = 2;
	}

	void SetNoFloorSearch(void)
	{
		getFloorVariablesSearch = 0;
	}

	void SetNoBalance(void)
	{
		// devo decidere quale è il default
		getBalance = 0;
	}

	void SetBalance(void)
	{
		getBalance = 1;
	}

	void SetNormalization(void)
	{
		// default
		getNormalization = 1;
	}

	void SetNoNormalization(void)
	{
		getNormalization = 0;
	}

	//TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// CAMBIARE in SetNoDegeneracyPrevention
	void SetNoDegeneracyPrevention(void)
	{
		getDegeneracyPrevenction = 0;
	}

	void SetDegeneracyPrevention(void)
	{
		getDegeneracyPrevenction = 1;
	}

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);

	//	============================================================================
	//	=====================< Modifying Functions >================================
	//	============================================================================
};

#endif // BZZ_LINEAR_PROGRAMMING_DOUBLE