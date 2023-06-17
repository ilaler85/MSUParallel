#pragma once
#include <iostream>
#include <vector>
#include <fstream>


using namespace std;
class Topol;
class TopolEN;
class TopolNE;
class TopolNeN;

class Topol
{
public:
	vector<int> Ja;
	vector<int> Ia;
	vector<int> A;
	int Nn;
	int Ne;
	void print();
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