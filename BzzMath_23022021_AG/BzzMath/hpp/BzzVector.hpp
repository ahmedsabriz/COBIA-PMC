// BZZMATH: Release 7.0

// =========================< BzzVector.hpp >============================
// * Class BzzVector for operations between matrices							*
// * and vectors in double precision														*
// * Description: Dal Fortan al C++	(Capitolo 12)										*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *							and																	*
// *					Scientific C++	(Chapter 12)											*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1993)														*
// *																									*
// * Examples: c:\bzzmath\examples\BzzMathBasic\LinearAlgebra\						*
// *           Vector\Vector.cpp											 	*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
// 04-1991	Date Written
// 11-1992	English version
// 01-1994	Modified Max, Min, MaxAbs, MinAbs.
// 03-1994	Added BzzBaseClass for BzzPrint and BzzMessage.
// 09-1994	Added shadow variable for returning object.
// 09-1994	Added functions ObjectCount, ObjectCountInScope.
// 09-1994	Added functions IProduct, ITProduct.
// 10-1994	Minor changes throughout.
// 11-1994	Conversion to double precision done.
// 02-1995	Added the Bzz prefix to the names of the classes.
// 12-1996	Added GetHandle function for mixed language: C++ / FORTRAN.
//	12-1996	Added DeleteOneElementEveryN.
//	01-1997	Added Stretch.
//	04-1997	Added InsertSortedVectorInSortedVector.
//	04-1997	Added InsertElementInSortedVector.
//	06-1998	Added BzzPowInt.
//	10-1998	Added IsVectorSorted.
//	01-1999	Added MedianAbs.
//	01-1999	Added LinearRegressionCoefficient.
//	02-1999	Added CenterAndNormalize.
//	03-1999	Added ElementByElementProduct.
//	03-1999	Added constructor from sub vector.

////////////////// Release 4.0
//	02-2000	Added GetSumElements.
//	03-2000	Added UnbiasedStandardDeviation.
//	03-2000	Added ResidualsUnbiasedNormalDeviate and MeanAbsoluteDeviation.
//	03-2000	Added StatisicalAnalysis.
//	03-2000	Added DeleteLastNElements.
//	06-2000	Added v == c, v > c, v >= c, v < c, v <= c.
//	02-2001	Added Linearize.
//	03-2001	Added Mad2.
//	03-2001	Added LeastMedianSquareLocation.
//	03-2001	Added TrimmedMedianDeviation.
//	04-2001	Added new version for Stretch.
//	05-2001	Added new version for GetBzzVector.
//	05-2001	Added new version for Stretch.
//	09-2001	Added new version for Sort.
//	04-2002	Added Reorder.
//	04-2002	LocateInFirstNSortedElements.
//	04-2002	Added MoveElementFromkToj.
//	04-2003	Added InsertElementInFirstNSortedElements.

////////////////// Release 5.0
//	10-2003	Added DeleteFirstNElements.
//	10-2003	Added DeleteFirstAndLastNElements.
//	10-2003	Added Center function.
//	10-2003	Added SetToZero function.
//	01-2004	Added TrimmedMean function.
//	11-2004	Added GetSumAbsElements.
//	01-2005	Added SumAsFractions and SumAsFractionsWithControl.
//	02-2005	Added UseSubVectorAsVector function.
//	02-2005	Added CopyDataFromVector function.
//	12-2006	Added MinMax function.
//	12-2006	Added AitkenAccelerator function.
//	06-2008	Added SwapMaxMinOrder function.

////////////////// Release 6.0
//	10-2008	Added NormM function.
//	10-2008	Added GetPerpendicularIntersection function.
//	11-2008	Added MixedErrorEstimation function.
//	11-2008	Added DeleteMin, DeleteMax and DelteMinMax functions.
//	11-2008	Added CleverMean function.
//	05-2009	Added CleverVariance function.
//	05-2009	Added GetCleverMeanVarianceOutliers function.
//	12-2009	Modified DeleteElement function.
//	01-2010	Added MinMaxMean function.
//	01-2010	Added ReplaceNElements function.
//	05-2010	Added BzzOutliersDetection class.
//	12-2011	Added BvpNormalize function.
//	11-2012	Added BzzVectorSparse class.
//	12-2012	Added Sum in BzzVectorSparse class.

////////////////// Release 7.0
//	04-2013	Added BuildForwardGradient function.
//	04-2013	Added BuildGradientAndHessian functions.
//	04-2013	Added BuildGradientAndPositiveHessian functions.
//	06-2013	Added SwappedReordering function.
//	06-2013	Added ReorderSubVector function.
//	06-2013	Added SwappedSubVectorReordering function.

// ============================================================================
// ****** Constructors for BzzVector:												*
// * BzzVector v; // default															*
// * BzzVector v = x; // copy-initializer											*
// * BzzVector v(n); // sizes v and set its coifficients = 0				*
// * BzzVector v(5,1.,2.,3.,4.,5.); // vector 1,2,3,4,5						*
// * double x[5]={1.,2.,3.,4.,5.};															*
// * BzzVector v(5,x); // from array												*
// * BzzVector w(5,v); // subvector of v											*
// * BzzVector w(5,3,v); // subvector of v starting from element 3		*
// * BzzVector v("VET.DAT"); // from formatted	file							*
// * BzzVector v('*',VET.BIN"); // from binary file							*
// ****************************************************************************
// ***** Access functions :																	*
// *	int i = v.Size(); // BzzVector v size										*
// *	int who = v.WhoAmI(); // BzzVector count									*
// *	double xf = v(i);	// with control													*
// *	v(i) = 4.;			// with control													*
// *	double xf = v[i];	// inline without control										*
// *	v[i] = 7.;			// inline without control										*
// *	double xf = v.GetValue(i); // with control										*
// *	v.SetValue(i,7.);			// with control											*
// *	BzzVector w = GetBzzVector(nc,i,v);// w subvector of v		*
// *	w.GetBzzVector(nc,i,v);// w subvector of v								*
// *	w.UseSubVectorAsVector(nc,i,v); // w subvector of v with the same data!!*
// *	v.SetBzzVector(i,w); // BzzVector w placed in v					*
// *	int count = BzzVector::ObjectCount();										*
// *	int countInScope = BzzVector::ObjectCountInScope();					*
// ****************************************************************************
// ***** Assignment:																				*
// *	v = w; // Assignment																		*
// *	v = 3.; // Assignment from double													*
// *	v = A; // Assignment	from BzzMatrix											*
// * w(5,v); // w subvector of v	dimensioned 5											*
// ****************************************************************************
// ***** BzzPrint and BzzMessage																*
// *	v.BzzPrint();																				*
// *	v.BzzMessage();																			*
// *	v.BzzPrint("The value of a(3) is %e",a[3]);										*
// *	v.BzzMessage("The value of a(3) is %e",a[3]);									*
// ****************************************************************************
// ***** Implemented operations :															*
// *	Sum (a,b,&c);	// c = a + b;															*
// *	c = a + b;		// c = a + b;															*
// *	Sum (a,b,&a);	// a = a + b;															*
// *	Sum (&a,b);	 // a = a + b;																*
// *	a += b;			// a = a + b;															*
// *	Sum (a,b,&b);	// b = a + b;															*
// *	Sum (b,&a);	 // a = b + a;																*
// *	Sum (a,a,&a);	// a = a + a;															*
// *	Sum (&a);		// a = a + a;															*
// *	Difference(a,b,&c);	// c = a - b;													*
// *	c = a - b;				// c = a - b;													*
// *	Difference(a,b,&a);	// a = a - b;													*
// *	Difference(&a,b);	 // a = a - b;														*
// *	a -= b;					// a = a - b;													*
// *	Difference(a,b,&b);	// b = a - b;													*
// *	Difference(b,&a);	 // a = b - a;														*
// *	Difference(a,a,&a);	// a = a - a;													*
// *	Difference(&a);		// a = a - a;													*
// *	Minus (a,&b);			// b = -a;														*
// *	b = -a;																						*
// *	Minus (&a);			 // a = -a;															*
// *	Product(A,x,&y);		// y = A*x; A BzzMatrix								*
// *	y = A*x;				 // y = A*x															*
// *	Product(A,x,&x);		// x = A*x;														*
// *	Product(A,&x);		 // x = A*x;														*
// *	Product(R,x,&y);		// y = R*x; R BzzMatrixRight								*
// *	Product(R,&x);			// x = R*x; R BzzMatrixRight								*
// *	y =R*x;																						*
// * Product(L,x,&y);		// y = L*x; L BzzMatrixLeft								*
// * Product(L,&x);			// x = L*x; L BzzMatrixLeft								*
// * Product(S,x,&y);		// y = S*x; S BzzMatrixSymmetric							*
// * Product(S,&x);			// x = S*x; S BzzMatrixSymmetric							*
// * y = S*x;																						*
// *	y = L*x;																						*
// *	Product(3.,x,&y);	 // y = 3.*x;														*
// *	y = 3.*x;																					*
// *	Product(3.,x,&x);	 // x = 3.*x;														*
// *	Product(3.,&x);		// x = 3.*x;													*
// *	x *= 3.;				 // x = 3.*x;														*
// *	TProduct(x,y,&xf);	// xf = xTy;													*
// *	xf = Dot(x,y);		 // xf = xTy;														*
// *	xf = x%y;				// xf = xTy;													*
// *	TProduct(A,x,&y);	 // y = ATx;														*
// *	y = A%x;				 // y = ATx;														*
// *	TProduct(A,x,&x);	 // x = ATx;														*
// *	TProduct(A,&x);		// x = ATx;														*
// * TProduct(L,x,&y);		// y = LTx;														*
// * TProduct(L,&x);			// x = LTx;														*
// * TProduct(R,x,&y);		// y = RTx;														*
// * TProduct(R,&x);			// x = RTx;														*
// *	ProductT(x,y,&A);	 // A = xyT;														*
// *	A = x->*y;				// A = xyT;														*
// *	IProduct(A,x,&y);		// y = A-1x;													*
// *	IProduct(A,&x);		// x = A-1x;													*
// * IProduct(D,x,&y);		// y = D-1x;													*
// * IProduct(D,&x);		// x = D-1x;														*
// * IProduct(L,x,&y);		// y = L-1x;													*
// * IProduct(L,&x);		// x = L-1x;														*
// * IProduct(R,x,&y);		// y = R-1x;													*
// * IProduct(R,&x);		// x = R-1x;														*
// * IProduct(S,x,&y);		// y = S-1x;													*
// * IProduct(S,&x);		// x = S-1x;														*
// *	ITProduct(A,x,&y);	// y = A-Tx;													*
// *	ITProduct(A,&x);		// x = A-Tx;													*
// * ITProduct(L,x,&y);	// y = L-Tx;														*
// * ITProduct(L,&x);		// x = L-Tx;													*
// * ITProduct(R,x,&y);	// y = R-Tx;														*
// * ITProduct(R,&x);		// x = R-Tx;													*
// * ElementByElementProduct(x,y,&z);														*
// *	Division(x,3.,&y);	// y = x/3.;													*
// *	y = x/3.;																					*
// *	Division(&x,3.);		// x = x/3.;													*
// *	x /= 3.;				 // x = x/3.;														*
// * SumAsFractions(x0,d,x);																	*
// * SumAsFractionsWithControl																*
// ****************************************************************************
// *****	Operators for forming vectors:													*
// *	v = w&&y; // appends one vector to another										*
// ****************************************************************************
// *****	Operators for tests:																	*
// *	if(v == w)																					*
// *	if(v != w)																					*
// * if(v == c)																					*
// * if(v > c)																						*
// * if(v >= c)																					*
// * if(v < c)																						*
// * if(v <= c)																					*
// ****************************************************************************
// *	Other functions:																			*
// *	double xf;																					*
// *	xf = v.Max();																				*
// *	xf = v.Max(&imax);																		*
// *	xf = v.MaxAbs();																			*
// *	xf = v.MaxAbs(&imax);																	*
// *	xf = v.Min();																				*
// *	xf = v.Min(&imin);																		*
// *	xf = v.MinAbs();																			*
// *	xf = v.MinAbs(&imin);																	*
// *	v.MinMax(&iMin,&min,&iMax,&max);														*
// *	v.MinMaxMean(&iMin,&min,&iMax,&max,&mean);										*
// *	xf = v.Norm1();																			*
// *	xf = v.Norm2();																			*
// *	xf = v.NormI();																			*
// *	xf = v.NormM();																			*
// *	xf = v.GetSumElements();																*
// *	xf = v.GetSumAbsElements();															*
// *	v.Save("Vet.DAT");																		*
// *	v.Save('*',"VET.BIN");																	*
// *	Load(&v,"VET.DAT");																		*
// *	Load(&v,'*',"VET.BIN");																	*
// *	Delete(&v); // eliminate a BzzVector										*
// *	v.SetToZeroDelete(); // set all elements to zero								*
// *	ChangeDimensions(newdim,&v);															*
// *	Swap(&x,&y);																				*
// *	norm = Normalize(&v); // vTv = 1														*
// *	mean = Center(&v); 																		*
// *	CenterAndNormalize(&v);																	*
// *	CenterAndNormalize(&v,&mean,&normCentered);										*
// *	Reverse(&v);																				*
// *	Sort(&v);																					*
// *	Sort(&v,&iS);																				*
// *	SwapMaxMinOrder(&v);																		*
// *	Reorder(&v,iS);																			*
// *	LocateInSortedVector(v,f);																*
// *	j = v.LocateInSortedVector(f);														*
// *	j = v.LocateInFirstNSortedElements(n,f);											*
// *	j = v.IsVectorSorted();// return 0 non sorted; 1 ascending; 2 descending*
// *	v.FirstInLastOut(f);																		*
// *	v.LastInFirstOut(f);																		*
// *	v.Append(f);																				*
// *	v.Append(w);																				*
// *	v.Insert(3,f);																				*
// *	v.Insert(3,w);																				*
// *	v.InsertSortedVectorInSortedVector(w);												*
// *	v.InsertElementInSortedVector(4.);													*
// *	v.InsertElementInFirstNSortedElements(10.,5);									*
// *	v.MoveElementFromkToj(k,j);															*
// *	v.DeleteElement(3);																		*
// *	v.DeleteLastNElements(n);																*
// *	v.DeleteFirstNElements(n);																*
// *	v.DeleteFirstAndLastNElements(n);													*
// *	v.DeleteElements(iv);																	*
// *	v.DeleteOneElementEveryN(n);															*
// *	v.DeleteMin(); v.DeleteMax(); v.DelteMinMax();									*
// *	v.Stretch(n);																				*
// *	v.Stretch(xOld,xNew);																	*
// *	v.Stretch(xOld,xNew,.8);																*
// *	Stretch(xOld,vOld,xNew,&vNew,.8);													*
// *	v.Linearize(x1,xn);																		*
// *	v.SwapElements(i,j);																		*
// *	v(n,w); // double w[n]; double *w = new double[n];								*
// *	double *vector = v.GetHandle();														*
// *	w = sin(v);																					*
// *	w = cos(v);																					*
// *	w = tan(v);																					*
// *	w = asin(v);																				*
// *	w = acos(v);																				*
// *	w = atan(v);																				*
// *	w = sinh(v);																				*
// *	w = cosh(v);																				*
// *	w = tanh(v);																				*
// *	w = sqrt(v);																				*
// *	w = log(v);																					*
// *	w = log10(v);																				*
// *	w = exp(v);																					*
// *	w = BzzPowInt(v,n);																		*
// ****************************************************************************
// ***** Functions for Statistics:															*
// *	f = Mean(v);																				*
// *	f = Mean("Vect.DAT");																	*
// *	f = Median(v);																				*
// *	f = Remedian(v);																			*
// *	f = Remedian("Vect.DAT");																*
// *	f = TrimmedMean(fraction,v);															*
// *	f = TrimmedMean(fraction,"Vect.DAT")												*
// *	f = TrimmedMean(v);																		*
// *	f = CleverMean(v);																		*
// *	f = CleverMean("Vect.dat");															*
// *	f = LeastMedianSquareLocation(v);													*
// *	f = Variance(v);																			*
// *	f = Variance(v,mean);																	*
// *	f = Mad2(v);																				*
// *	f = Mad2(v,median);																		*
// *	f = Mad2("Vect.DAT");																	*
// *	f = Mad2("Vect.DAT",median);															*
// *	f = CleverVariance(v);																	*
// *	f = CleverVariance("Vect.dat");																	*
// *	f = StandardDeviation(v);																*
// *	f = StandardDeviation(v,mean);														*
// *	f = UnbiasedStandardDeviation(v);													*
// *	f = UnbiasedStandardDeviation(v,mean);												*
// *	f = Mad(v);																					*
// *	f = Mad(v,median);																		*
// *	f = Mad("Vect.DAT");																		*
// *	f = Mad("Vect.DAT",median);															*
// *	f = MeanAbsoluteDeviation(v);															*
// *	f = Range(v);																				*
// *	w = ResidualsNormalDeviate(v);														*
// *	w = ResidualsUnbiasedNormalDeviate(v);												*
// *	w = ResidualsNormalDeviate(v,mean,std);											*
// *	w = ResidualsRobustDeviate(v);														*
// *	w = ResidualsRobustDeviate(v,median,mad);											*
// *	v.StatisticalAnalysis("Analysis");													*
// ****************************************************************************

#ifndef BZZ_VECTOR_DOUBLE_HPP
#define BZZ_VECTOR_DOUBLE_HPP

// preventive statements
class BzzVectorInt;
class BzzVector;
class BzzMatrix;
class BzzMatrixDiagonal;
class BzzMatrixLeft;
class BzzMatrixRight;
class BzzMatrixSymmetric;

class BzzMatrixSparse;
class BzzMatrixSparseLockedByRows;
class BzzMatrixSparseLockedByColumns;

// ============================================================================
// =============================< class BzzVector >============================
// ============================================================================

class BzzVector : public BzzBaseClass
{
	friend class BzzMatrix;
	friend class BzzMatrixDiagonal;
	friend class BzzMatrixRight;
	friend class BzzMatrixLeft;
	friend class BzzMatrixSymmetric;

	friend class BzzMatrixSparse;
	friend class BzzMatrixSparseLokedByRows;
	friend class BzzMatrixSparseLokedByColumns;

	friend class BzzFactorized;
	friend class BzzFactorizedPLR;
	friend class BzzFactorizedGauss;
	friend class BzzFactorizedQRLQ;
	friend class BzzFactorizedQR;
	friend class BzzFactorizedLQ;
	friend class BzzFactorizedSVD;
	friend class BzzFactorizedSymmetricBase;
	friend class BzzFactorizedSparseGauss;
	friend BzzVector LinearLeastSquaresWithLinearEqualityConstraints
	(BzzMatrix& A, BzzVector& b,
		BzzMatrix& C, BzzVector& d);
	friend BzzVector LinearLeastSquaresWithLinearEqualityConstraints
	(BzzMatrix* A, BzzVector* b,
		BzzMatrix* C, BzzVector* d);
	friend BzzVector NonNegativeLinearLeastSquares
	(BzzMatrix& A, BzzVector& b);
	friend BzzVector NonNegativeLinearLeastSquares
	(BzzMatrix* A, BzzVector* b);
	friend class BzzSave;
	friend class BzzLoad;
	friend void Load(BzzMatrixSparse* A, char, char* filematrix);
	friend void AitkenAccelerator(BzzVector& x1, BzzVector& x2, BzzVector& x3,
		BzzVector* x);
	friend void AitkenAccelerator(BzzVector& x1, BzzVector& x2, BzzVector& x3,
		BzzVector& x4, BzzVector& x5, BzzVector& x6,
		BzzVector* x);
	friend void AitkenAccelerator(BzzVector& x1, BzzVector& x2, BzzVector& x3,
		BzzVector& x4, BzzVector& x5, BzzVector& x6,
		BzzVector& x7, BzzVector& x8, BzzVector* x);
	friend void GetPerpendicularIntersection(BzzVector& x, BzzVector& a, double b,
		BzzVector* z);
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;
	static const double ETA;
	double* vector;
	int dimensions;
	int whoAmI;
	char shadow;
	char matrixAsVector;
	char subVectorAsVector;

	// initialise constructors
	void Initialize(int nc);

	// private constructor BzzVector('*',nc)
	BzzVector(char, int nc);

public:

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor BzzVector v;
	BzzVector(void);

	// copy-initializer
#if BZZ_COMPILER > 100
	BzzVector(const BzzVector& rval);
#else
	BzzVector(BzzVector& rval);
#endif
	// sizes and initialises at 0
	BzzVector(int nc);

	// sizes and initialises
	BzzVector(int nc, double v1, ...);

	// data pointer constructor
	BzzVector(int nc, double* initvalues);

	// subvector of v
	BzzVector(int nc, BzzVector& v);

	// subvector of v starting from i
	BzzVector(int nc, int i, BzzVector& v);

	// file v(FILE) constructor;
	BzzVector(char* filevector);

	// unformatted v('*',FILE) constructor saved with Save
	BzzVector(char, char* filevector);

	// ============================================================================
	// *****************************< destructor >	*******************************
	// ============================================================================
	~BzzVector(void);

	// ============================================================================
	// ***********************< Non-modifying access functions >*******************
	// ============================================================================
	int Size(void) const { return dimensions; } // dimensions

	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

	// receives the value of the vector with control
	double GetValue(int i) const;

	// subvector (similar to BzzVector(nc,ielem,BzzVector)
	friend BzzVector GetBzzVector
	(int nc, int ielem, const BzzVector& rval);

	// ============================================================================
	// ************************< Modifying access functions >**********************
	// ============================================================================
		// subvector (similar to BzzVector(nc,ielem,BzzVector)
	void GetBzzVector
	(int nc, int ielem, const BzzVector& rval);
	void UseSubVectorAsVector(int nc, int ielem, const BzzVector& rval);

	// assigns and receives vector values with control
	double& operator () (int i);

	// assigns and receives values without control
	double& operator []
	(int i) {
		return vector[i];
	}

	const double& operator []
	(int i) const {
		return vector[i];
	}

	// assigns vector values with control
	void SetValue(int i, double val);

	// modifies the vector starting from ielem with rval
	void SetBzzVector(int ielem, const BzzVector& rval);

	// ============================================================================
	// ***************************< assignment operators >*************************
	// ============================================================================
	BzzVector& operator =
		(const BzzVector& rval);

	void operator =
		(double c);

	BzzVector& operator =
		(const BzzMatrix& rval);

	//w.(n,v);
	void operator() (int n, BzzVector& v);
	void operator() (int n, int start, BzzVector& v);

	// ============================================================================
	// **********************< operators for composing vectors >*******************
	// ============================================================================
			//adds one vector to another in sequence
	friend BzzVector operator &&
		(const BzzVector& lval, const BzzVector& rval);

	// ============================================================================
	// **************************< operators for tests >***************************
	// ============================================================================
	friend char operator ==
		(const BzzVector& lval, const BzzVector& rval);

	friend char operator !=
		(const BzzVector& lval, const BzzVector& rval);

	char operator != (double c);
	char operator == (double c);
	char operator > (double c);
	char operator >= (double c);
	char operator < (double c);
	char operator <= (double c);

	// ============================================================================
	// ****************< Gradient and Hessian building >***************************
	// ============================================================================
	void BuildForwardGradient(BzzVector& x, double F, double (*ptrF)(BzzVector& x));
	void BuildGradientAndHessian(BzzVector& x, double F, double (*ptrF)(BzzVector& x),
		BzzMatrixSymmetric* G);
	void BuildGradientAndHessian(BzzVectorInt& r, BzzVectorInt& c,
		BzzVector& x, double F, double (*ptrF)(BzzVector& x),
		BzzMatrixSymmetric* G);
	void BuildGradientAndPositiveHessian(BzzVector& x, double F, double (*ptrF)(BzzVector& x),
		BzzMatrixSymmetric* G);
	void BuildGradientAndPositiveHessian(BzzVectorInt& r, BzzVectorInt& c,
		BzzVector& x, double F, double (*ptrF)(BzzVector& x),
		BzzMatrixSymmetric* G);

	// ============================================================================
	// ============================< OPERATIONS >==================================
	// ============================================================================

	// ============================================================================
	// *******************************< Sum >**************************************
	// ============================================================================
	friend void Sum(const BzzVector& lval,
		const BzzVector& rval,
		BzzVector* result); // Sum (a,b,&c); c = a + b;

	friend BzzVector operator +
		(const BzzVector& lval, const BzzVector& rval);

	// Sum (&a,b); a = a + b;
	friend void Sum(BzzVector* lvalAndResult,
		const BzzVector& rval);

	BzzVector& operator +=
		(const BzzVector& rval);

	// Sum (b,&a); a = b + a;
	friend void Sum(const BzzVector& lval,
		BzzVector* rvalAndResult);

	// Sum (&a); a = a + a;
	friend void Sum(BzzVector* lvalRvalAndResult);
	friend double SumAsFractions(BzzVector& x0, BzzVector& d,
		BzzVector* x, double, double);
	friend double SumAsFractionsWithControl(BzzVector* x0, BzzVector* d,
		BzzVector* x);

	// ============================================================================
	// *****************************< Difference >*********************************
	// ============================================================================
		// Difference(a,b,&c); c = a - b;
	friend void Difference(const BzzVector& lval,
		const BzzVector& rval, BzzVector* result);

	// c = a - b;
	friend BzzVector operator -
		(const BzzVector& lval, const BzzVector& rval);

	// Difference(&a,b); a = a - b;
	friend void Difference(BzzVector* lvalAndResult,
		const BzzVector& rval);

	// a -= b; a = a - b;
	BzzVector& operator -=
		(const BzzVector& rval);

	// Difference(b,&a); a = b - a;
	friend void Difference(const BzzVector& lval,
		BzzVector* rvalAndResult);
	// Difference(&a); a = a - a;
	friend void Difference(BzzVector* lvalRvalAndResult);

	// ============================================================================
	// *******************************< Minus >************************************
	// ============================================================================
		// Minus (a,&b); b = -a;
	friend void Minus(const BzzVector& rval,
		BzzVector* result);

	BzzVector operator -();// unary minus

	// Minus (&a); a = -a;
	friend void Minus(BzzVector* rvalAndResult);

	// ============================================================================
	// ******************************< Product >***********************************
	// ============================================================================
		// Product(A,x,&y); y = A*x;
	friend void Product
	(const BzzMatrix& lval, const BzzVector& rval,
		BzzVector* result);

	// y = A*x;
	friend BzzVector operator *
		(const BzzMatrix& lval, const BzzVector& rval);

	// Product(A,&x); x = A*x;
	friend void Product
	(const BzzMatrix& lval, BzzVector* rvalAndresult);

	friend BzzVector operator *
		(BzzMatrixDiagonal& lval, BzzVector& rval);

	// Product(L,x,&y); y = L*x;
	friend void Product
	(const BzzMatrixLeft& lval, const BzzVector& rval,
		BzzVector* result);

	// Product(L,&x); x = L*x;
	friend void Product(const BzzMatrixLeft& lval, BzzVector* rvalAndResult);

	friend BzzVector operator *	 // y = L*x;
		(const BzzMatrixLeft& lval, const BzzVector& rval);

	// Product(R,x,&y); y = R*x;
	friend void Product
	(const BzzMatrixRight& lval, const BzzVector& rval,
		BzzVector* result);

	// Product(R,&x); x = R*x;
	friend void Product
	(const BzzMatrixRight& lval, BzzVector* rvalAndResult);

	friend BzzVector operator *	 //y =R*x;
		(const BzzMatrixRight& lval, const BzzVector& rval);

	// Product(S,x,&y); y = S*x;
	friend void Product
	(BzzMatrixSymmetric& lval, const BzzVector& rval,
		BzzVector* result);

	// Product(S,&x); x = S*x;
	friend void Product
	(BzzMatrixSymmetric& lval, BzzVector* rvalAndResult);

	friend BzzVector operator *
		(BzzMatrixSymmetric& lval, const BzzVector& rval);

	friend BzzVector operator *
		(BzzMatrixSparse& lval, BzzVector& rval);

	// Product(P,x,&y); y = P*x;
	friend void Product
	(const BzzMatrixSparse& lval, const BzzVector& rval,
		BzzVector* result);

	// Product(P,&x); x = P*x;
	friend void Product(const BzzMatrixSparse& lval,
		BzzVector* rvalAndResult);

	// Product(3.,x,&y); y = 3.*x;
	friend void Product
	(double lval, const BzzVector& rval, BzzVector* result);

	friend BzzVector operator *
		(double lval, const BzzVector& rval);

	// Product(3.,&x); x = 3.*x;
	friend void Product
	(double lval, BzzVector* rvalAndResult);

	BzzVector& operator *= (double rval);

	// ============================================================================
	// ***************************< TProduct, Dot >********************************
	// ============================================================================
		// TProduct(x,y,&c); c = xTy; c = x%y; c = Dot(x,y);
	friend void TProduct(BzzVector& lval,
		BzzVector& rval, double* result);

	// c = Dot(x,y); c = xTy;
	friend double Dot(BzzVector& lval,
		BzzVector& rval);

	// c = x%y; c = xTy;
	friend double operator %
		(BzzVector& lval, BzzVector& rval);

	// TProduct(A,x,&y); y =ATx; y =A%x;
	friend void TProduct(BzzMatrix& lval,
		BzzVector& rval, BzzVector* result);

	// y = ATx;
	friend BzzVector operator %
		(BzzMatrix& lval, BzzVector& rval);

	// TProduct(A,&x); x = ATx;
	friend void TProduct
	(BzzMatrix& lval, BzzVector* rvalAndresult);

	// TProduct(L,x,&y);
	friend void TProduct
	(BzzMatrixLeft& lval, BzzVector& rval,
		BzzVector* result);

	// TProduct(L,&x);
	friend void TProduct
	(BzzMatrixLeft& lval, BzzVector* rvalAndResult);

	// TProduct(R,x,&y);
	friend void TProduct
	(BzzMatrixRight& lval, BzzVector& rval,
		BzzVector* result);

	// TProduct(R,&x);
	friend void TProduct
	(BzzMatrixRight& lval, BzzVector* rvalAndResult);

	// TProduct(P,x,&y);
	friend void TProduct
	(BzzMatrixSparse& lval, BzzVector& rval,
		BzzVector* result);

	// TProduct(P,&x);
	friend void TProduct
	(BzzMatrixSparse& lval, BzzVector* rvalAndResult);

	// y = L%x;
	friend BzzVector operator %
		(BzzMatrixLeft& lval, BzzVector& rval);

	// y = R%x;
	friend BzzVector operator %
		(BzzMatrixRight& lval, BzzVector& rval);

	// y = P%x;
	friend BzzVector operator %
		(BzzMatrixSparse& lval, BzzVector& rval);

	// ============================================================================
	// ******************************< ProductT >**********************************
	// ============================================================================
		// ProductT(x,y,&A); A = xyT; A = x->*y;
	friend void ProductT(const BzzVector& lval,
		const BzzVector& rval, BzzMatrix* result);

	friend BzzMatrix operator ->*
		(const BzzVector& lval, const BzzVector& rval);

	// ============================================================================
	// ******************************< IProduct >**********************************
	// ============================================================================
		// IProduct(A,x,&y); x = A-1y;
	friend void IProduct(BzzMatrix& A, BzzVector& x,
		BzzVector* y);

	// IProduct(A,&x); x = A-1x;
	friend void IProduct(BzzMatrix& A, BzzVector* x);

	// IProduct(D,x,&y);
	friend void IProduct
	(BzzMatrixDiagonal& lval, BzzVector& rval,
		BzzVector* result);

	// IProduct(D,&x);
	friend void IProduct
	(BzzMatrixDiagonal& rval, BzzVector* rvalAndResult);

	// IProduct(L,x,&y);
	friend void IProduct
	(const BzzMatrixLeft& lval, const BzzVector& rval,
		BzzVector* result);

	// IProduct(L,&x);
	friend void IProduct
	(const BzzMatrixLeft& lval, BzzVector* rvalAndResult);

	// IProduct(R,x,&y);
	friend void IProduct
	(const BzzMatrixRight& lval, const BzzVector& rval,
		BzzVector* result);

	// IProduct(R,&x);
	friend void IProduct
	(const BzzMatrixRight& lval, BzzVector* rvalAndResult);

	// IProduct(S,x,&y);
	friend void IProduct
	(const BzzMatrixSymmetric& lval, const BzzVector& rval,
		BzzVector* result);

	// IProduct(S,&x);
	friend void IProduct
	(const BzzMatrixSymmetric& lval, BzzVector* rvalAndResult);

	// IProduct(P,x,&y);
	friend void IProduct
	(const BzzMatrixSparse& lval, const BzzVector& rval,
		BzzVector* result);

	// IProduct(P,&x);
	friend void IProduct
	(const BzzMatrixSparse& lval, BzzVector* rvalAndResult);

	// ============================================================================
	// *****************************< ITProduct >**********************************
	// ============================================================================
		// ITProduct(A,x,&y); x = A-Ty;
	friend void ITProduct(BzzMatrix& A, BzzVector& x,
		BzzVector* y);

	// ITProduct(A,&x); x = A-Tx;
	friend void ITProduct(BzzMatrix& A, BzzVector* x);

	// ITProduct(L,x,&y);
	friend void ITProduct
	(const BzzMatrixLeft& lval, const BzzVector& rval,
		BzzVector* result);

	// ITProduct(L,&x);
	friend void ITProduct
	(const BzzMatrixLeft& lval, BzzVector* rvalAndResult);

	// ITProduct(R,x,&y);
	friend void ITProduct
	(const BzzMatrixRight& lval, const BzzVector& rval,
		BzzVector* result);

	// ITProduct(R,&x);
	friend void ITProduct
	(const BzzMatrixRight& lval, BzzVector* rvalAndResult);

	// ITProduct(P,x,&y);
	friend void ITProduct
	(const BzzMatrixSparse& lval, const BzzVector& rval,
		BzzVector* result);

	// ITProduct(P,&x);
	friend void ITProduct
	(const BzzMatrixSparse& lval, BzzVector* rvalAndResult);

	// ============================================================================
	// ***************************< ElementByElementProduct >**********************
	// ============================================================================
	friend void ElementByElementProduct(BzzVector& l, BzzVector& r,
		BzzVector* s);

	// ============================================================================
	// ******************************< Division >**********************************
	// ============================================================================

		// Division(x,3.,&y); y = x/3.;
	friend void Division
	(const BzzVector& lval, double rval, BzzVector* result);

	friend BzzVector operator /		// y = x/3.
		(const BzzVector& lval, double rval);

	// Division(&x,3.); x = x/3.;
	friend void Division
	(BzzVector* lvalAndResult, double rval);

	// x /= 3.; x = x/3.;
	BzzVector& operator /= (double rval);

	// ============================================================================
	// =======================< Non-modifying functions >==========================
	// ============================================================================

	// *****************************< BzzPrint >***********************************
	virtual void ObjectBzzPrint(void);
	void StatisticalAnalysis(const char* s = "");
	void Save(char* filevector); // formatted
	void Save(char, char* filevector);// binary

// ***************************< Max and Min >**********************************
	double Max(int* imax = 0); //to have the Max position (im)
	double MaxAbs(int* imax = 0);

	double Min(int* imin = 0);
	double MinAbs(int* imin = 0);
	void MinMax(int* iMin, double* min, int* iMax, double* max);
	void MinMaxMean(int* iMin, double* min, int* iMax, double* max, double* mean);

	// ******************************< Norms >*************************************
	double Norm1(void);
	double Norm2(void);
	double NormI(void);
	double NormM(void);
	double GetSumElements(void);
	double GetSumAbsElements(void);

	// ============================================================================
	// ===========================< Modifying Functions >==========================
	// ============================================================================
	friend void Delete(BzzVector* result); // eliminates BzzVector
	void SetToZero(void)
	{
		memset(vector, 0, (dimensions + 1) * sizeof(double));
	}
	void SetToZero(int start, int end);
	friend void ChangeDimensions(int dim, BzzVector* result, char);

	// recovery from Save
	friend void Load
	(BzzVector* result, char* filevector); // formatted
	friend void Load
	(BzzVector* result, char, char* filevector);// binary

	friend void SumRankOne(const BzzVector& u,
		const BzzVector& vT, BzzMatrix* result);
	friend void SumRankOne(double product, const BzzVector& u,
		const BzzVector& vT, BzzMatrix* result);
	friend void SumRankOne(const BzzVector& u,
		const BzzVector& vT, double divisor, BzzMatrix* result);

	//normalises (xTx=1) and returns norm2
	friend double Normalize(BzzVector* result);
	friend double Center(BzzVector* result);
	friend double CenterAndNormalize(BzzVector* result);
	friend void CenterAndNormalize(BzzVector* result,
		double* mean, double* normCentered);

	friend void Reverse(BzzVector* result);	//inverts the vector

	friend void Sort(BzzVector* result); // orders the vector
	friend void Sort(BzzVector* result, BzzVectorInt* iS); // orders the vector

	// return vector v ordered with iS
	friend void Reorder(BzzVector* v, BzzVectorInt& iS);
	// return vector v ordered with iS
	friend void SwappedReording(BzzVector* v, BzzVectorInt& iS);
	// return a portion of vector v starting from k ordered with iS
	friend void ReorderSubVector(int k, BzzVector* v, BzzVectorInt& iS);
	// return a portion of vector v starting from k ordered with iS
	friend void SwappedSubVectorReordering(int k, BzzVector* v, BzzVectorInt& iS);

	// swaps the contents of two vectors
	friend void Swap(BzzVector* lval, BzzVector* rval);
	friend void CopyDataFromVector(BzzMatrix* lval, BzzVector& rval);
	//reverse the order Max min
	friend void SwapMaxMinOrder(BzzVector* x);

	// ============================================================================
	// ==========================< Other functions >===============================
	// ============================================================================
	friend int LocateInSortedVector(BzzVector& v, double f);
	int LocateInSortedVector(double f);
	int LocateInFirstNSortedElements(int n, double f);
	int IsVectorSorted(void); // return 0 non sorted; 1 ascending; 2 descending
	void FirstInLastOut(double f);
	void LastInFirstOut(double f);
	void Append(double f);
	void Append(BzzVector& v);
	void Insert(int i, double f);
	void Insert(int i, BzzVector& v);
	int InsertElementInSortedVector(double f);
	int InsertElementInFirstNSortedElements(int nE, double f);
	void InsertSortedVectorInSortedVector(BzzVector& vj);
	void ReplaceNElements(BzzVector& b);
	void ReplaceNElements(int j, BzzVector& b);
	void MoveElementFromkToj(int k, int j);
	double DeleteElement(int i);
	void DeleteLastNElements(int n);
	void DeleteFirstNElements(int n);
	void DeleteFirstAndLastNElements(int n);
	void DeleteElements(BzzVectorInt v);
	void DeleteOneElementEveryN(int n);
	void DeleteMax(void);
	void DeleteMin(void);
	void DeleteMinMax(void);
	void Stretch(int i);
	void Stretch(BzzVector& xOld, BzzVector& xNew, double);
	friend void Stretch(BzzVector& xOld, BzzVector& vOld,
		BzzVector& xNew, BzzVector* vNew, double);
	void Linearize(double x1, double x2);
	void SwapElements(int i, int j);
	void operator () (int n, double* w);
	double* GetHandle(void) { return vector + 1; }
	void BvpNormalize(void);
	friend BzzVector sin(BzzVector& v);
	friend BzzVector cos(BzzVector& v);
	friend BzzVector tan(BzzVector& v);
	friend BzzVector asin(BzzVector& v);
	friend BzzVector acos(BzzVector& v);
	friend BzzVector atan(BzzVector& v);
	friend BzzVector sinh(BzzVector& v);
	friend BzzVector cosh(BzzVector& v);
	friend BzzVector tanh(BzzVector& v);
	friend BzzVector sqrt(BzzVector& v);
	friend BzzVector log(BzzVector& v);
	friend BzzVector log10(BzzVector& v);
	friend BzzVector exp(BzzVector& v);
	friend BzzVector BzzPowInt(BzzVector& v, int n);

	// ============================================================================
	// *******************************< Statistics >*******************************
	// ============================================================================

	friend double Mean(BzzVector& v);
	friend double Mean(char* filevector);
	friend double Median(BzzVector& v);
	friend double MedianAbs(BzzVector& v);
	friend double Remedian(BzzVector& v);
	friend double Remedian(char* filevector);
	friend double TrimmedMean(double fraction, BzzVector& v);
	friend double TrimmedMean(double fraction, char* filevector);
	friend double CleverMean(BzzVector& v);
	friend double CleverMean(char* filevector);
	//	friend double TrimmedMean(BzzVector &v);
	friend double LeastMedianSquareLocation(BzzVector& v);
	friend double Variance(BzzVector& v);
	friend double Variance(BzzVector& v, double vm);
	friend double Variance(char* filevector);
	friend double Variance(char* filevector, double vm);
	friend double Mad2(BzzVector& v);
	friend double Mad2(BzzVector& v, double md);
	friend double Mad2(char* filevector);
	friend double Mad2(char* filevector, double md);
	friend double CleverVariance(BzzVector& v);
	friend double CleverVariance(char* filevector);
	friend double StandardDeviation(BzzVector& v);
	friend double StandardDeviation(BzzVector& v, double vm);
	friend double UnbiasedStandardDeviation(BzzVector& v);
	friend double UnbiasedStandardDeviation(BzzVector& v, double vm);
	friend double Mad(BzzVector& v);
	friend double Mad(BzzVector& v, double md);
	friend double Mad(char* filevector);
	friend double Mad(char* filevector, double md);
	friend double TrimmedMedianDeviation(BzzVector& v, int h);
	friend double TrimmedMedianDeviation(BzzVector& v, int h, double md);
	friend double MeanAbsoluteDeviation(BzzVector& v);
	friend double Range(BzzVector& v);
	friend double LinearRegressionCoefficient(BzzVector& y,
		BzzVector& x);
	friend void GetCleverMeanVarianceOutliers(BzzVector& v,
		double* mean, double* variance, BzzVectorInt* iout, BzzVector* out);
	friend void GetCleverMeanVarianceOutliers(char* filevector,
		double* mean, double* variance, BzzVectorInt* iout, BzzVector* out);
	friend void GetCleverMeanVarianceOutliers(BzzVector& v,
		BzzVector* mean, BzzVector* variance, BzzVectorInt* iout, BzzVector* out);
	friend void GetCleverMeanVarianceOutliers(char* filevector,
		BzzVector* mean, BzzVector* variance, BzzVectorInt* iout, BzzVector* out);

	friend BzzVector ResidualsNormalDeviate(BzzVector& v);
	friend BzzVector ResidualsUnbiasedNormalDeviate(BzzVector& v);
	friend BzzVector ResidualsNormalDeviate(BzzVector& v,
		double vm, double std);
	friend BzzVector ResidualsRobustDeviate(BzzVector& v);
	friend BzzVector ResidualsRobustDeviate(BzzVector& v, double md,
		double mad);
	friend double MixedErrorEstimation(BzzVector& x, BzzVector& s, BzzVector* e,
		double, double);
};

// ============================================================================
// ====================< class BzzOutliersDetection >==========================
// ============================================================================

class BzzOutliersDetection : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	BzzVector vectorBase;
	int numElements, numElementsBase, numIteration, numElementsModified, maxClever;
	int numParameters;
	int whoAmI;
	BzzVector cleverMean, cleverVariance, outliers, cleverSumOfSquares;
	BzzVector max, min;
	BzzVectorInt iOutliers, iSort, iOriginal;
	BzzVectorInt imin, imax;
	double sum, sumq, cm, cv, va, standard, clsos, testValue;
	// initialise constructors
	void Initialize(void);

public:
	double CleverMean(void) { return cm; }
	double CleverVariance(void) { return cv; }
	void GetCleverMeanVarianceOutliers(BzzVector* cm, BzzVector* cv,
		BzzVectorInt* iout, BzzVector* out) {
		*cm = cleverMean;*cv = cleverVariance;
		*iout = iOutliers;*out = outliers;
	}
	void GetMinMax(BzzVectorInt* imi, BzzVector* mi, BzzVectorInt* ima, BzzVector* ma)
	{
		*imi = imin;*mi = min;*ima = imax;*ma = max;
	}
	double CleverSumOfSquares(void) { return clsos; }
	void GetCleverSumOfSquaresOutliers(BzzVector* clsos,
		BzzVectorInt* iout, BzzVector* out) {
		*clsos = cleverSumOfSquares;
		*iout = iOutliers;*out = outliers;
	}
	void SetTestValue(double ts = 2.5) { if (ts > 0.)testValue = ts; }

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor BzzOutliersDetection;
	BzzOutliersDetection(void);
	BzzOutliersDetection(BzzVector* vc, int mxmima = 500);
	void operator()(BzzVector* vc, int mxmima = 500);
	BzzOutliersDetection(char* fileVector, int mxmima = 500);
	void operator()(char* fileVector, int mxmima = 500);
	BzzOutliersDetection(int nump, BzzVector* r);
	void operator()(int nump, BzzVector* r);
	virtual void ObjectBzzPrint(void);
};

// ============================================================================
// =======================< class BzzVectorSparse >============================
// ============================================================================
class BzzMatrixSparseLocked;

class BzzVectorSparse : public BzzBaseClass
{
	friend class BzzMatrixSparseLocked;
	friend void ProductT(BzzVectorSparse& x, BzzVectorSparse& y,
		BzzMatrixSparseLocked* SL);
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	int dimensions, numElements;
	int whoAmI;
	BzzVectorInt r;
	BzzVector v;
	friend void Swap(BzzVectorSparse* x, BzzVectorSparse* y)
	{
		Swap(&x->dimensions, &y->dimensions);
		Swap(&x->numElements, &y->numElements);
		Swap(&x->r, &y->r);
		Swap(&x->v, &y->v);
	}

public:

	int WhoAmI(void) const { return whoAmI; }
	int Size(void) { return dimensions; }
	int GetElementsNumber(void) { return numElements; }
	// ============================================================================
	// *****************************< constructor >	*******************************
	// ============================================================================
		// default
	BzzVectorSparse(void)
	{
		dimensions = 0;count++;countInScope++;whoAmI = count;
	}
	void operator ()(int dim, BzzVectorInt& rr, BzzVector& vv);

	// ============================================================================
	// *****************************< destructor >	*******************************
	// ============================================================================
	~BzzVectorSparse(void) { countInScope--; }

	virtual void ObjectBzzPrint(void);

	// ============================================================================
	// ============================< OPERATIONS >==================================
	// ============================================================================
	//	TProduct >> xTy
	friend double TProduct(BzzVector& aux, BzzVectorSparse& x, BzzVectorSparse& y);
	friend double TProduct(BzzVector& x, BzzVectorSparse& y);
	friend double TProduct(BzzVectorSparse& x, BzzVector& y);
	// Product y = Ax
	// Questa � in BzzMatrixSparseLocked.hpp
	friend void Product(BzzMatrixSparseLocked& SL, BzzVector& x, BzzVector* y);
	// Questa � in BzzMatrixSparseLocked.hpp
	friend void Product(BzzMatrixSparseLocked& SL, BzzVectorSparse& x, BzzVector* y);
	friend void Product(BzzMatrix& A, BzzVectorSparse& x, BzzVector* y);
	friend void ProductT(BzzVectorSparse& x, BzzVectorSparse& y,
		BzzMatrixSparseLocked* SL);
	friend void ProductTWithoutPivot(int ip, int jp, BzzVectorSparse& x, BzzVectorSparse& y,
		BzzMatrixSparseLocked* SL);
	// Sum
	friend void Sum(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
		BzzVectorSparse& x, BzzVectorSparse& y, BzzVectorSparse* z, int);
	friend void Sum(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
		BzzVectorInt& rx, BzzVector& vx, BzzVectorInt& ry, BzzVector& vy,
		BzzVectorInt* rz, BzzVector* vz, int);
	// Difference
	friend void Difference(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
		BzzVectorSparse& x, BzzVectorSparse& y, BzzVectorSparse* z, int);
	friend void Difference(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
		BzzVectorInt& rx, BzzVector& vx, BzzVectorInt& ry, BzzVector& vy,
		BzzVectorInt* rz, BzzVector* vz, int);

	friend double GetRatioForCoupleOfElements(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
		BzzVectorInt& rx, BzzVector& vx, BzzVectorInt& ry, BzzVector& vy,
		BzzVectorInt* rz, BzzVector* vz, BzzVectorInt* re, BzzVector* ve);
	friend double GetRatioForCoupleOfElements(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
		BzzVectorInt& rx, BzzVector& vx, BzzVectorInt& ry, BzzVector& vy,
		BzzVectorInt* rz, BzzVector* vz, BzzVectorInt* re,
		BzzVectorInt* r1, BzzVector* v1, BzzVectorInt* r2, BzzVector* v2,
		double epsylon);
};

// Friend functions with default arguments
// SumAsFractions
double SumAsFractions(BzzVector& x0, BzzVector& d,
	BzzVector* x, double dper = 0., double uxper = 0.);
// ChangeDimensions
void ChangeDimensions(int dim, BzzVector* result, char zero = 0);
// Stretch
void Stretch(BzzVector& xOld, BzzVector& xNew, double smooth = .5);
void Stretch(BzzVector& xOld, BzzVector& vOld,
	BzzVector& xNew, BzzVector* vNew, double smooth = .5);
double MixedErrorEstimation(BzzVector& x, BzzVector& s, BzzVector* e,
	double alfa = 1., double beta = 1.);
// Sum
void Sum(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
	BzzVectorSparse& x, BzzVectorSparse& y, BzzVectorSparse* z, int zero = 1);
void Sum(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
	BzzVectorInt& rx, BzzVector& vx, BzzVectorInt& ry, BzzVector& vy,
	BzzVectorInt* rz, BzzVector* vz, int zero = 1);
// Difference
void Difference(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
	BzzVectorSparse& x, BzzVectorSparse& y, BzzVectorSparse* z, int zero = 1);
void Difference(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
	BzzVectorInt& rx, BzzVector& vx, BzzVectorInt& ry, BzzVector& vy,
	BzzVectorInt* rz, BzzVector* vz, int zero = 1);
// GetRatioForCoupleOfElements
double GetRatioForCoupleOfElements(BzzVectorInt& jaux, BzzVectorInt& iaux, BzzVector& aux,
	BzzVectorInt& rx, BzzVector& vx, BzzVectorInt& ry, BzzVector& vy,
	BzzVectorInt* rz, BzzVector* vz, BzzVectorInt* re,
	BzzVectorInt* r1, BzzVector* v1, BzzVectorInt* r2, BzzVector* v2,
	double epsylon = 1.e-12);


#endif // BZZ_VECTOR_DOUBLE_HPP