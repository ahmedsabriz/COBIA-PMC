
// BZZMATH: Release 7.0

// ===================================< BzzQuadraticProgramming.HPP >============================
// * Class BzzQuadraticProgramming																					*
// * Description:												 																	*
// *																																	*
// * Examples:																														* 
// * C:\BzzMath\Examples\BzzMathAdvanced\Minimization\QuadraticProgramming\QuadraticProgramming.cpp	*
// ====================================================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	01-2008	Date Written.															 

// ============================================================================
// ****** Constructors for BzzQuadraticProgramming:							*
// * BzzQuadraticProgramming qp; // default										*
// * qp(&S,&s,&E,&e,&D,&d,&xMin,&xMax);													*
// * BzzQuadraticProgramming qp(&S,&s,&E,&e,&D,&d,&xMin,&xMax);			*
// ****************************************************************************
// ***** Access functions :																	*
// **
// ****************************************************************************
// ***** Assignment:																				*
// **
// ****************************************************************************
// ***** BzzPrint																					*
// * o.BzzPrint();																				*
// **
// ****************************************************************************
// ***** Implemented operations :															*
// **
// ****************************************************************************
// ***** Other functions:																		*
// **
// ****************************************************************************
// ****************************************************************************

#ifndef BZZ_QUADRATIC_PROGRAMMING_DOUBLE
#define BZZ_QUADRATIC_PROGRAMMING_DOUBLE


// ============================================================================
// ============================< class BzzQuadraticProgramming >=============================
// ============================================================================

class BzzQuadraticProgramming : public BzzBaseClass
	{
private:
	static const char *const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;	

	// initialise constructors
//	void Initialize(int nc);

	int	whoAmI,
			printTasks,
			printSubTasks;
	int	numVariables,
			numEqualityConstraints,
			numInequalityConstraints,
			numActiveConstraints;

	double	F,
				tolAbs,
				tolRel;

//	BzzMatrixSymmetric S;
	BzzMatrix S,E,D;
	BzzVector s,e,d,xMin,xMax;
	BzzMatrix A,Q,N;
	BzzVector x,z,b,q,l,v,y;
	BzzFactorizedGauss G;
	BzzFactorizedLQ AA;

public:

// ============================================================================
// ***************************< constructors >*********************************
// ============================================================================
	// default constructor BzzQuadraticProgramming v;
	BzzQuadraticProgramming(void); 

	// copy-initializer
	BzzQuadraticProgramming(BzzQuadraticProgramming &rval)
		{BzzError("%s copy-initializer not implemented",BzzQuadraticProgramming::BZZ_ERROR);}
	
	// other constructor
	BzzQuadraticProgramming
		(BzzMatrixSymmetric *S,BzzVector *s,
		BzzMatrix *E,BzzVector *e,
		BzzMatrix *D,BzzVector *d,
		BzzVector *xMin,BzzVector *xMax);

	void operator()
		(BzzMatrixSymmetric *S,BzzVector *s,
		BzzMatrix *E,BzzVector *e,
		BzzMatrix *D,BzzVector *d,
		BzzVector *xMin,BzzVector *xMax);

// ============================================================================
// *****************************< destructor >*********************************
// ============================================================================
	~BzzQuadraticProgramming(void){};

// ============================================================================
// *******************< Non-modifying access functions >***********************
// ============================================================================
	int WhoAmI(void) const {return whoAmI;}
	static int ObjectCount(void){return count;}
	static int ObjectCountInScope(void){return countInScope;}

// ============================================================================
// ******************************< Setting functions >*************************
// ============================================================================
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

// ============================================================================
// **********************< Modifying access functions >************************
// ============================================================================

// ============================================================================
// *************************< assignment operators >***************************
// ============================================================================

// ============================================================================
// ====================< Non-modifying functions >=============================
// ============================================================================

// ******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void){};

//	********************************< Save >************************************
//	void Save(char *file***); // formatted
//	void Save(char,char *file***);// binary

// ============================================================================
// ======================< Modifying Functions >===============================
// ============================================================================

	// recovery from Save
//	friend void Load
//			(BzzQuadraticProgramming *result,char *file***); // formatted
//	friend void Load
//			(BzzQuadraticProgramming *result,char,char *file***);// binary

// ============================================================================
// =========================< Other functions >================================
// ============================================================================
	void operator()(int maxIt = 100);
};

#endif // BZZ_QUADRATIC_PROGRAMMING_DOUBLE


