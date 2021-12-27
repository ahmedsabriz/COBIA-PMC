// BZZMATH: Release 7.0

//	==========================< BzzComplex.hpp >==========================
//	* Class BzzComplex for operations between complex numbers				*
// * Description:																					*
// *																									*
//	============================================================================

// Revision History (MM-YYYY)
//	05-1994	Date Written.
//	Author	Guido Buzzi-Ferraris
//	02-1995	Added the Bzz prefix to the names of the classes.

// BZZMATH: Release 7.0
// 10-2008

//	============================================================================
//	****** Constructors for BzzComplex:												*
//	* BzzComplex c; // default															*
//	* BzzComplex c = z; // copy-initializer										*
//	* BzzComplex c(3,2); // initializer												*
//	* BzzComplex c(3.F,2.F); // initializer										*
//	* BzzComplex c(3.,2.); // initializer											*
//	* BzzComplex c = 1; // conversion from int									*
//	* BzzComplex c = 1.F; // conversion from double								*
//	* BzzComplex c = 1.; // // conversion from double							*
//	****************************************************************************
//	***** Access functions (c = x +iy):														*
//	*	double x = GetReal(c); 																	*
//	*	double y = GetImaginary(c); 															*
//	*	c(2,3); // c = 2 + i3																	*
//	*	c(2.F,3.F); // c = 2 + i3																*
//	*	c(2.,3.); // c = 2 + i3																	*
//	****************************************************************************
//	***** Assignment:																				*
//	*	c = z; // from complex																	*
//	*	c = z = w; // multiple assignments													*
//	*	c = 2; // from int																		*
//	*	c = 2.F; // from double																	*
//	*	c = 2.; // from double																	*
//	****************************************************************************
//	***** BzzPrint and BzzMessage																*
//	*	c.BzzPrint();																				*
//	*	c.BzzMessage();																			*
//	*	c.BzzPrint("The real part of c is %e",c.GetReal());							*
//	*	c.BzzMessage("The real part of c is %e",c.GetReal());							*
//	****************************************************************************
//	***** c, z, w BzzComplex															*
//	***** f int or double or double															*
//	***** Implemented operations :															*
//	*	c = z + w;																					*
//	*	c = f + z;																					*
//	*	c = z + f;																					*
//	*	c += z;																						*
//	*	c += f;																						*
//	*	c = z - w;																					*
//	*	c = f - z;																					*
//	*	c = z - f;																					*
//	*	c -= z;																						*
//	*	c -= f;																						*
//	*	c = z * w;																					*
//	*	c = f * z;																					*
//	*	c = z * f;																					*
//	*	c *= z;																						*
//	*	c *= f;																						*
//	*	c = z / w;																					*
//	*	c = f / z;																					*
//	*	c = z / f;																					*
//	*	c /= z;																						*
//	*	c /= f;																						*
//	****************************************************************************
//	*****	Operators for tests:																	*
//	*	if(a == b)																					*
//	*	if(a != b)																					*
//	****************************************************************************
//	*****	Other functions:																		*
//	*	int k;																						*
//	*	double f;																					*
//	*	BzzComplex c,z;																	*
//	*	Swap(&c,&z);																				*
//	*	f = Modulus(c);																			*
//	*	f = Argument(c);																			*
//	*	z = Inverse(c); // z = 1./c															*
//	*	z = Conjugate(c);																			*
//	*	z = exp(c);																					*
//	*	z = sin(c);																					*
//	*	z = cos(c);																					*
//	*	z = tan(c);																					*
//	*	z = asin(c);																				*
//	*	z = acos(c);																				*
//	*	z = atan(c);																				*
//	*	z = sinh(c);																				*
//	*	z = cosh(c);																				*
//	*	z = tanh(c);																				*
//	*	z = sinh(c);																				*
//	*	z = cosh(c);																				*
//	*	z = tanh(c);																				*
//	*	z = log(c);																					*
//	*	z = log10(c);																				*
//	*	z = sqrt(c);																				*
//	*	z = pow(c,k);																				*
//	*	z = pow(c,f);																				*
//	*	z = pow(c,z);																				*
//	****************************************************************************

#ifndef BZZ_COMPLEX_DOUBLE_HPP
#define BZZ_COMPLEX_DOUBLE_HPP

//	============================================================================
//	========================< class BzzComplex >==========================
//	============================================================================

class BzzComplex : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static const double BZZ_LOG10;
	static int count; // for whoAmI
	int whoAmI;

	char	expxState,
		expyState,
		modulusState,
		argumentState;

	double complex[2],
		expx, cosy, siny,	// expxState = 1
		expy, cosx, sinx,	// expyState = 1
		modulus,				// modulusState = 1
		argument;			// argumentState = 1

	void Initialize(void);

public:

	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default constructor BzzComplex c;
	BzzComplex(void);

	// copy-initializer
	BzzComplex(const BzzComplex& c);

	// initializer
	BzzComplex(int r, int i);

	// initializer
	BzzComplex(float r, float i);

	// initializer
	BzzComplex(double r, double i);

	// conversion
	BzzComplex(int r);

	// conversion
	BzzComplex(float r);

	// conversion
	BzzComplex(double r);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzComplex(void) {};

	//	============================================================================
	//	****************************< Access functions >****************************
	//	============================================================================
	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }

	friend inline double GetReal(const BzzComplex& c)
	{
		return c.complex[0];
	}

	friend inline double GetImaginary(const BzzComplex& c)
	{
		return c.complex[1];
	}

	void operator()(int r, int i);

	void operator()(float r, float i);

	void operator()(double r, double i);

	//	============================================================================
	//	*************************< assignment operators >***************************
	//	============================================================================
	BzzComplex& operator =
		(const BzzComplex& c);

	//	============================================================================
	//	*************************< operators for tests >****************************
	//	============================================================================
	friend char operator ==
		(const BzzComplex& lval, const BzzComplex& rval);

	friend char operator !=
		(const BzzComplex& lval, const BzzComplex& rval);

	//	============================================================================
	//	=============================< OPERATIONS >=================================
	//	============================================================================

	//	============================================================================
	//	********************************< Sum >*************************************
	//	============================================================================
	friend const BzzComplex operator +
		(const BzzComplex& lval, const BzzComplex& rval);

	BzzComplex& operator +=
		(const BzzComplex& rval);

	//	============================================================================
	//	*****************************< Difference >*********************************
	//	============================================================================
		// c = a - b;
	friend BzzComplex operator -
		(const BzzComplex& lval, const BzzComplex& rval);

	// a -= b; a = a - b;
	BzzComplex& operator -=
		(const BzzComplex& rval);

	//	============================================================================
	//	********************************< Minus >***********************************
	//	============================================================================
	BzzComplex operator -();// unary minus

//	============================================================================
//	*******************************< Product >**********************************
//	============================================================================
	//c = a * b;
	friend BzzComplex operator *
		(const BzzComplex& lval, const BzzComplex& rval);

	// c = a * 3.F;
	friend BzzComplex operator * (const BzzComplex& lval, double x);

	// c = 3.F*a;
	friend BzzComplex operator * (double x, const BzzComplex& rval);

	// a *= b; a = a * b;
	BzzComplex& operator *= (const BzzComplex& rval);

	// c *= 3.;
	BzzComplex& operator *= (double x);
	//	============================================================================
	//	*******************************< Division >*********************************
	//	============================================================================
		// c = a / b;
	friend BzzComplex operator /
		(const BzzComplex& lval, BzzComplex& rval);

	// c = a/3.F
	friend BzzComplex operator / (const BzzComplex& lval, double x);

	// c = 3./a
	friend BzzComplex operator / (double x, BzzComplex& rval);

	// a /= b; a = a / b;
	BzzComplex& operator /=
		(const BzzComplex& rval);

	// c /= 3.;
	BzzComplex& operator /= (double x);

	//	============================================================================
	//	====================< Non-modifying functions >=============================
	//	============================================================================

	//	******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);

	//	============================================================================
	//	======================< Modifying Functions >===============================
	//	============================================================================
	friend void Swap(BzzComplex* lval, BzzComplex* rval);
	friend double Modulus(BzzComplex& c);
	friend double Argument(BzzComplex& c);
	friend BzzComplex Inverse(const BzzComplex& c);
	friend BzzComplex Conjugate(const BzzComplex& c);

	friend BzzComplex sin(BzzComplex& c);
	friend BzzComplex cos(BzzComplex& c);
	friend BzzComplex tan(BzzComplex& c);
	friend BzzComplex asin(BzzComplex& c);
	friend BzzComplex acos(BzzComplex& c);
	friend BzzComplex atan(BzzComplex& c);

	friend BzzComplex sinh(BzzComplex& c);
	friend BzzComplex cosh(BzzComplex& c);
	friend BzzComplex tanh(BzzComplex& c);
	friend BzzComplex asinh(BzzComplex& c);
	friend BzzComplex acosh(BzzComplex& c);
	friend BzzComplex atanh(BzzComplex& c);

	friend BzzComplex exp(BzzComplex& c);

	friend BzzComplex log(BzzComplex& c);
	friend BzzComplex log10(BzzComplex& c);

	friend BzzComplex sqrt(BzzComplex& c);

	friend BzzComplex pow(BzzComplex& c, int k);
	friend BzzComplex pow(BzzComplex& c, double a);
	friend BzzComplex pow(BzzComplex& c, BzzComplex& a);
};

#endif // BZZ_COMPLEX_DOUBLE_HPP