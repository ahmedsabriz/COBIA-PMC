// BZZMATH: Release 7.0
// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	03-2013	Date Written.

#ifndef BZZ_OPTNOV_NEW_HPP
#define BZZ_OPTNOV_NEW_HPP

class BzzMyMinimizationMonoOptnov : public BzzMyMinimizationMonoObject
{
	int var;
	BzzVector x;
public:
	double (*ptrMulti)(BzzVector& x);
	BzzMyMinimizationMonoOptnov(void) {}
	void operator()(int v, BzzVector* xx)
	{
		var = v;x = *xx;
	}
	virtual double GetFunctionValue(double t);
	virtual void ObjectBzzPrint(void);
};

class BzzOptnov : public BzzBaseClass
{
private:
	int numVariables;
	BzzMyMinimizationMonoOptnov funM;
	BzzMinimizationMonoObject m;
	BzzVector	d,
		x, xL, xU;
	BzzVectorInt iSort;
	double fi1, fi, f;
public:
	BzzOptnov(void) {}
	int operator()(BzzVector& xi1, double fi1, BzzVector& xi, double fi,
		double Name(BzzVector& x), BzzVector* xx, double* ff);
	virtual void ObjectBzzPrint(void) {}
};

#endif // BZZ_OPTNOV_NEW_HPP