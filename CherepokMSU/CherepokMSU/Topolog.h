#pragma once
#include <iostream>
#include <vector>
#include <fstream>


using namespace std;
class Topol;
class TopolEN;
class TopolNE;
class TopolNeN;

class Solver;

class Topol
{
public:
	vector<int> Ja;
	vector<int> Ia;
	vector<double> A;
	int Nn=0;
	int Ne=0;
	void print();

	void fillA();
	
};

class TopolEN :public Topol
{
public:
	void generator(int width, int height, int square, int triangle);
};

class TopolNE :public Topol
{
public:
	void convertNE(TopolEN& EN);
};

class TopolNeN :public Topol
{
public:
	void conevertNeN(TopolNE& NE, TopolEN& EN);
};

class Solver
{
public:
	void SLAU(Topol& topol, double eps, int maxit, vector<double>& b, vector<double>& x, int& k, double& res);
	void SpMV(Topol& topol, vector<double>& b, vector<double>& res);
	void dot(vector<double>& a, vector<double>& b, double& res);
	void axpy(vector<double>& a, vector<double>& b, double alfa, vector<double>& res);
	void inverseM(Topol& A, Topol& M);
	double norma(vector<double>& a);
};