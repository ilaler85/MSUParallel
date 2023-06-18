#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <math.h>
#include "Topolog.h"

using namespace std;

void Topol::print()
{
	cout << "Ia\n";
	for (size_t i = 0; i < Ia.size(); ++i)
		cout << Ia[i] << endl;
	
	cout << "\nJa\n";
	for (rsize_t i = 0; i < Ja.size(); ++i)
		cout << Ja[i] << endl;
}

void Topol::fillA()
{
	A.resize(Ja.size());
	for (int i = 0; i < Ia.size()-1; ++i)
	{
		double sum = 0;
		int diag = -1;
		for (int col = Ia[i]; col < Ia[i + 1]; ++col)
		{
			int j = Ja[col];
			if (i != j)
			{
				double a = cos(i*j+i+j);
				A[col] = a;
				sum += abs(a);
			}
			else
				diag = col;
		}
		if (diag != -1)
			A[diag] = 2 * sum;
	}

}	

void TopolEN::generator(int width, int height, int square, int triangle)
{
	Ia.push_back(0);
	int countNodes = 0;
	for (int i = 0; i < width; ++i)
	{
		for (int j = 0; j < height; ++j)
		{
			int set = square + triangle;
			int position = i * width + j;
			int NodeLeftUp = i * (width + 1) + j;
			int NodeRightUp = NodeLeftUp + 1;
			int NodeLeftDown = (i + 1) * (width + 1) + j;
			int NodeRightDown = NodeLeftDown + 1;
			if (position % set < square)
			{
				Ja.push_back(NodeLeftUp);
				Ja.push_back(NodeRightUp);
				Ja.push_back(NodeLeftDown);
				Ja.push_back(NodeRightDown);
				countNodes += 4;
				Ia.push_back(countNodes);
			}
			else
			{
				Ja.push_back(NodeLeftUp);
				Ja.push_back(NodeRightUp);
				Ja.push_back(NodeLeftDown);
				countNodes += 3;
				Ia.push_back(countNodes);
				Ja.push_back(NodeRightUp);
				Ja.push_back(NodeLeftDown);
				Ja.push_back(NodeRightDown);
				countNodes += 3;
				Ia.push_back(countNodes);
			}
		}
	}
	Ne = (int)Ia.size() - 1;
	Nn = (width + 1) * (height + 1);
}


void TopolNE::convertNE(TopolEN& EN)
{
	Ne = EN.Ne;
	Nn = EN.Nn;
	Ia.assign(Nn + 1, 0);
	for (int i = 0; i < Ne; ++i)
	{
		for (int j = EN.Ia[i]; j < EN.Ia[i + 1]; ++j)
		{
			++Ia[EN.Ja[j] + 1];
		}
	}

	for (int i = 1; i < Nn + 1; ++i)
	{
		Ia[i] += Ia[i - 1];
	}

	Ja.assign(Ia[Nn], 0);
	for (int i = 0; i < Ne; ++i)
	{
		for (int j = EN.Ia[i]; j < EN.Ia[i + 1]; ++j)
		{
			int node = EN.Ja[j];
			int ipos = Ja[Ia[node + 1] - 1]++;
			Ja[Ia[node] + ipos] = i;
		}
	}
}

void TopolNeN::conevertNeN(TopolNE& NE, TopolEN& EN)
{
	Nn = EN.Nn;
	Ia.reserve(Nn + 1);
	Ia.push_back(0);
	for (int node = 0; node < Nn; ++node)
	{
		set<int> sosedi;
		for (int i = NE.Ia[node]; i < NE.Ia[node + 1]; ++i)
		{
			int elem = NE.Ja[i];
			for (int j = EN.Ia[elem]; j<EN.Ia[elem+1];++j)
			{
				sosedi.insert(EN.Ja[j]);
			}
		}

		for (auto i:sosedi)
		{
			Ja.push_back(i);
		}
		Ia.push_back(Ia.back() + sosedi.size());
	}
}

void vectorPrint(vector<double>& v)
{
	for (auto d : v)
		cout << d << "  ";
	cout << endl;
}


void Solver::SLAU(Topol& topol, double eps, int maxit, vector<double>& b, vector<double>& x, int& k, double& res)
{
	int N = topol.Nn;
	cout << N << "\t" << topol.Ia.size() - 1 << endl;
	k = 0;
	double ro = 0;
	double roBack = 0;
	double betta = 0;
	double alfa = 0;
	Topol M;
	inverseM(topol, M);
	vector<double> r = b;
	vector<double> z(N);
	vector<double> p(N);
	vector<double> q(N);
	do
	{
		++k;
		SpMV(M ,r, z);
		ro = 0;

		dot(r, z, ro);
		if (k == 1)
		{
			p = z;
		}
		else
		{
			betta = ro / roBack;
			axpy(p, z, betta, p);
		}
		SpMV(topol, p, q);
		double tmp = 0;
		dot(p, q, tmp);
		alfa = ro / tmp;
		axpy(p, x, alfa, x);
		axpy(q, r, -alfa, r);
		roBack = ro;
		res = norma(r);
		cout << res <<endl;
	} while ((ro > eps * eps) && (k < maxit));
	
	
}


double Solver::norma(vector<double>& a)
{
	double result = 0;
	dot(a, a, result);
	return sqrt(result);
}

void Solver::inverseM(Topol& topol, Topol& M)
{
	for (int i = 0; i < topol.Ia.size() -1; ++i)
	{
		M.Ia.push_back(i);
		M.Ja.push_back(i);
	}
	M.Ia.push_back(topol.Ia.size() - 1);

	for (int i = 0; i < topol.Ia.size() - 1; ++i)
	{
		int diag = -1;
		for (int col = topol.Ia[i]; col < topol.Ia[i + 1]; ++col)
		{
			int j = topol.Ja[col];
			if (i == j)
			{
				M.A.push_back(1 / topol.A[col]);
				break;
			}
		}

	}
}

void Solver::SpMV(Topol& topol, vector<double>& b, vector<double>& res)
{
	
	for (int i = 0; i < topol.Ia.size()-1; ++i)
	{
		double localsum = 0;
		
		for (int j = topol.Ia[i]; j < topol.Ia[i+1]; ++j)
		{
			int stolbec = topol.Ja[j];
			localsum += b[stolbec] * topol.A[j];
		}
		res[i] = localsum;
	}
}

void Solver::axpy(vector<double>& x, vector<double>& y, double alfa, vector<double>& res)
{
	for (int i = 0; i < x.size(); ++i)
	{
		res[i] = x[i] * alfa + y[i];
	}
}

void Solver::dot(vector<double>& a, vector<double>& b, double& res)
{
	for (int i = 0; i < a.size(); ++i)
	{
		res += a[i] * b[i];
	}
}

void fillB(vector<double>& b, int width)
{
	for (int i = 0; i < width; ++i)
	{
		b.push_back(sin(i));
	}
}


int main(int args, char* argv[])
{
	int width = 100;
	int height = 100;
	int square = 5;
	int triangle = 10;

	double eps = 0.00001;
	int maxit = 1000;
	/*if (args == 1)
	{
		cout << "The input is 4 numbers n1, n2 grid sizes k4, k3 alternating square and triangular cells";
		return 1;
	}
	if (args < 5)
	{
		cout << "Input error count of data";
	}
	width = atoi(argv[1]);
	height = atoi(argv[2]);
	if (width < 1 || height < 1)
	{
		cout << "Error entering the grid size";
		return 2;
	}
	square = atoi(argv[3]);
	triangle = atoi(argv[4]);
	if (square < 0 || triangle < 0)
	{
		cout << "Grid cell alternation input error";
		return 3;
	}*/
	/*Topol test;
	
	test.Ia = {0,2,3};
	test.Ja = {0,1,1};
	test.A = { 2,1,1 };
	test.Nn = 2;
	vector<double> a1 = { 1, 2, 3 };
	vector<double> b1 = { 3, 2, 3 };
	vector<double> res (3);
	Solver sol;
	Topol Mi;
	sol.inverseM(test, Mi);
	sol.SpMV(test,a1,res);
	for (auto i : Mi.A)
		cout << i << endl;*/

	TopolNE topNE;
	TopolEN topEN;
	TopolNeN topNeN;
	topEN.generator(width, height, square, triangle);
	topNE.convertNE(topEN);
	topNeN.conevertNeN(topNE, topEN);
	topNeN.fillA();
	vector<double> b;
	fillB(b, topNeN.Ia.size()-1);
	vector<double> x(topNeN.Ia.size() - 1, 0);
	int k;
	double res;
	Solver sol;
	sol.SLAU(topNeN, eps,maxit,b,x,k,res);
	//topNeN.print();
}