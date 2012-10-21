#include "path.h"

Path::Path(void)
{
	cost = 0;
	demand = 0;
}

Path::Path(list <int> n, double c, double d)
{
	nodes = n;
	cost = c;
	demand = d;
}

Path::Path(const Path& p)
{
	nodes = p.nodes;
	cost = p.cost;
	demand = p.demand;
}


Path::~Path(void)
{
}

void Path::calc_demand(int day, double ** de_rank)
{
	demand = 0;
	list <int>::iterator it;
	for (it = nodes.begin() ; it != nodes.end() ; it++)
	{
		if (*it != 0)
			demand += de_rank[*it - 1][day];
	}
}

void Path::print(void)
{
	cout << "Path: node: ";
	for (list<int>::iterator it = nodes.begin() ; it!= nodes.end() ; it++)
	{
		cout << *it << " ";
	}
	cout << endl << "cost : " << cost << " demand : " << demand << endl;
}

bool Path::operator== ( Path &p2)  
{  
	list<int>::iterator it1, it2;

	if (nodes.size() != p2.nodes.size()) 
		return false;

	for (it1 = nodes.begin(), it2 = p2.nodes.begin() ; it1 != nodes.end() && it2 != p2.nodes.end() ; it1++, it2++)
		if (*it1 != *it2)
			return false;

	return true;
}  
  
bool Path::operator!= (Path &p2)  
{  
    list<int>::iterator it1, it2;

	if (nodes.size() != p2.nodes.size()) 
		return true;

	for (it1 = nodes.begin(), it2 = p2.nodes.begin() ; it1 != nodes.end() && it2 != p2.nodes.end() ; it1++, it2++)
		if (*it1 != *it2)
			return true;

	return false;
}  
