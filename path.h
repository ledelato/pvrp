#pragma once
#include <list>
#include <iostream>

using namespace std;

/**
 * Path represents the route of a single vehicle on a single day
 */
class Path {
public:
  list<int> nodes;
  double cost;
  double demand;
  Path(void);
  Path(list<int> n, double c, double d);
  Path(const Path& p);
  ~Path(void);
  void calc_demand(int day, double ** de_rank);
  void print(void);
  bool operator ==(Path &p2);
  bool operator !=(Path &p2);
};

