#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
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

int main(int args, char* argv[])
{
	int width = 3;
	int height = 3;
	int square = 2;
	int triangle = 1;
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


	TopolNE topNE;
	TopolEN topEN;
	TopolNeN topNeN;
	topEN.generator(width, height, square, triangle);
	topNE.convertNE(topEN);
	topNeN.conevertNeN(topNE, topEN);
	topNeN.print();
}