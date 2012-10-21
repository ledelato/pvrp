#include <iostream>
#include <list>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <cassert>
#include <cfloat>
#include "path.h"
#define DEBUG 0
#define MAXN 200
#define MAXV 20
#define MAXD 20
#define EP 1e-9

using namespace std;

/*parameters*/
double cap_factor = 10; /*phi: penalty factor for overcapacity*/
double div_factor = 0.015;
int p = 5;
double delta = 0.1;
double alpha = 0;
double beta = 0;
int i_total = 50; /*ita: total number of iterations*/
int init_type = 0; /*initialization type: 1 - maintain consistency, 0 - no consistency*/
double tabu_weight = 15; /*weight to adjust tabu length*/

char* name = "input/tabu/convrp4.txt";
char* routename = "output/alpha/set1/prob6/route.txt";
char* recordname = "tabu/record.txt";
char *paraname = "input/parameter/para1.txt";

int n_total = 0; /*total number of nodes*/
int d_total = 0; /*total number of days*/
int l_total = 0; /*total number of regions*/
int v_total = 0; /*total number of vehicles*/
int cap; /*vehicle capacity*/

double **arr_d; /*pointer to the array of distance matrix*/
double **new_d; /*pointer to the ordered distance matrix*/
double **co_old; /*coordinate date input - raw*/
double **co; /*ordered coordinate*/
int *arr_reg; /*pointer to the array of region information*/
double **demand; /*pointer to the array of demand matrix*/
double **de_rank; /*pointer to the ranked demand array*/
int *reg; /*region record for each node--after rank*/

double *con_weight; /*customer familiarity weights*/
double *con_weight_copy;
double *reg_weight; /*region familiarity weights*/
double *reg_weight_copy;

vector<vector<vector<vector<int> > > > neighbours;
vector<vector<int> > visit_time; //visit time by vehicles
vector<vector<int> > region_time; //region visit time by vehicles
vector<vector<int> > best_visit_time; //visit time by vehicles--for best feasible solutions
vector<vector<int> > best_region_time; //region visit time by vehicles--for best feasible solutions
vector<vector<vector<int> > > region_freq(d_total,
    vector<vector<int> >(v_total, vector<int>(l_total, 0))); /*how many nodes in the region that a vehicle visits*/

//Cache
Path * remove_cache[MAXN][MAXV][MAXD];
Path * insert_cache[MAXN][MAXV][MAXD];

//Timing
double geni_insert_time = 0;
double geni_remove_time = 0;
double cal_cost_time = 0;
double frequency_time = 0;
double cal_obj_time = 0;

//debugging print
void printroute(const list<int>& route) {
  cout << "printroute" << endl;
  list<int> temproute;
  temproute = route;
  list<int>::iterator it;

  for (it = temproute.begin(); it != temproute.end(); it++) {
    cout << *it << " ";
  }
  cout << "-100" << endl;

  //temproute.~list();
}

/*This function calculates the traveling cost along one route*/
double tcost(const Path& oldroute, double **d)
/* link<int>route: the link list of route
 double d: distance matrix*/
{
  /*check for exception*/

  Path newroute;
  newroute = oldroute;

  if (newroute.nodes.empty())
    return 0;

  //int i, j;
  /*for(i = 0; i<=n_total; i = i+1)
   {
   for (j = 0; j<n_total; j = j+1)
   {
   cout<<d[i][j]<<" ";
   }
   cout<<endl;
   }*/

  double t = 0;
  int node1 = 0;
  int node2 = 0;
  list<int>::iterator it = newroute.nodes.begin();
  ++it;
  while (it != newroute.nodes.end()) {
    --it;
    node1 = *it;
    ++it;
    node2 = *it;
    t = t + d[node1][node2];
    ++it;/*move iterator forward*/
  }
  --it;
  node1 = *it;
  it = newroute.nodes.begin();
  node2 = *it;
  t = t + d[node1][node2]; /*the distance between the last node and depot*/

  //newroute.~list();
  return t;
}

/*nearest p neighbor for node i, on day r*/
vector<int> p_nb(int p, int node, double **d, int day, int totalserve) {
  vector<int> nearnb(p, 0);
  //nearnb = new int [p];
  return nearnb;

  int index, i, j;

  if (node != 0) {
    vector<vector<double> > temp(totalserve + 3, vector<double>(5, 0));
    /*for every node except for depot*/
    index = 0;
    temp[index][0] = d[node][0];
    temp[index][1] = 0;
    index = index + 1;

    for (i = 1; i <= n_total; i = i + 1) {
      //cout<<"temp"<<endl;
      if ((demand[i - 1][day] != 0) && (i != node)) {
        temp[index][0] = new_d[node][i];
        temp[index][1] = i;
        //cout<<temp[index][0]<<" "<<temp[index][1]<<endl;
        index = index + 1;

      }
    }
    if (index != totalserve) {
      cout << "p_nb error!" << endl;
    }
    sort(temp.begin(), temp.end());

    j = 0;
    for (i = 0; i < n_total; i = i + 1) {
      if (j < p) {
        if ((temp[i][1] != 0) || (temp[i][0] != 0)) {
          nearnb[j] = (int) temp[i][1];
          //cout<<nearnb[j]<<" ";
          j = j + 1;
        }
      }
    }
    //cout<<endl;
  } else {
    vector<vector<double> > temp(totalserve + 3, vector<double>(5, 0));
    /*for depot*/
    index = 0;

    for (i = 1; i <= n_total; i = i + 1) {
      if (demand[i - 1][day] != 0) {
        temp[index][0] = d[0][i];
        temp[index][1] = i;
        index = index + 1;
      }
    }
    if (index != totalserve) {
      cout << "p_nb error!" << endl;
    }
    sort(temp.begin(), temp.end());

    j = 0;
    for (i = 0; i < n_total; i = i + 1) {
      if (j < p) {
        if ((temp[i][1] != 0) || (temp[i][0] != 0)) {
          nearnb[j] = (int) temp[i][1];
          //cout<<nearnb[j]<<" ";
          j = j + 1;
        }
      }
    }
    //cout<<endl;

  }

  //cout<<"the end of nearnb"<<endl;

  return nearnb; /*the end of nearest neighborhood searching*/
}

/*This function determines if node2 is in the p-neighborhood of node1*/
bool isinp(int node1, int node2, const vector<vector<int> >& nearnb, int p) {
  int i;
  //cout<<"isinp"<<endl;
  for (i = 0; i < p; i = i + 1) {
    //cout<<nearnb[node1][i]<<endl;
    if (nearnb[node1][i] == node2)
      return 1;
  }
  return 0;
}

/*This function determines if route contains q of p-neighborhood of node*/
bool isnear(const list<int>& oldroute, int node,
    const vector<vector<int> > &neighbor, int p, int q) {
  //printroute(route);

  list<int> route;
  route = oldroute;

  int n = 0; /*counter of already containted neighbors*/
  list<int>::iterator it;
  bool test;

  for (it = route.begin(); it != route.end(); ++it) {
    //cout<<*it<<endl;
    test = isinp(node, *it, neighbor, p);
    //cout<<"isinp:"<<test<<endl;
    if (test && (*it != 0))
      n = n + 1; /*exclude depot*/
  }
  //route.~list();

  if (n >= q)
    return 1;
  else
    return 0;
}

/*This function determines if node i is in the candidate set W(lambda)*/
bool isinw(int node, const vector<int>& candidate) {
  /*int i = 0;
   for (i = 0; i<(int)candidate.size(); i = i+1)
   {
   if( candidate[i]==node)
   return 1;
   }
   return 0;*/

  return 1; //testing
}

/*This function checks if the assignment is feasible*/
//TODO: Remove d_total and the other global variables from here.
bool isfeasible(const vector<vector<int> >& assign, double cap, double **demand,
    int d_total, int v_total, int n_total)
    /*assign: assignment table;
     cap: capacity of vehicles
     demand: demand table;
     d_total: total period length
     v_total: total number of vehicles
     n_total: total number of nodes*/
    {
  int r, i, j; /*index*/
  double *fill; /*the left capacity on vehicle*/
  fill = new double[v_total + 3];
  int temp;

  for (r = 0; r < d_total; r = r + 1) /*on each day*/
  {
    for (j = 0; j < v_total; j = j + 1) /*initialization*/
    {
      fill[j] = cap;
    }
    for (i = 1; i <= n_total; i = i + 1) /*for each node*/
    {
      temp = assign[i - 1][r] - 1;
      fill[temp] = fill[temp] - demand[i - 1][r];
      //cout<<demand[i-1][r]<<endl;
    }
    for (j = 0; j < v_total; j = j + 1) /*check for feasibility*/
    {
      if (fill[j] < 0) {
        return 0;
      }
    }
  }

  delete[] fill;
  return 1;
}

list<int> remove1(const list<int>& original_route, int j, int k, int oldnode,
    int*myints) {

  list<int>::iterator it;
  list<int> seg1, seg2, seg3;
  list<int> route;
  list<int> newroute;
  list<int> newnewroute;

  list<int> oldroute;
  oldroute = original_route;

  int length = (int) oldroute.size();

  //order old route
  it = find(oldroute.begin(), oldroute.end(), oldnode);

  route.insert(route.end(), ++it, oldroute.end());
  route.insert(route.end(), oldroute.begin(), --it);

  int n = 0;
  int *newints;
  newints = new int[length - 1];
  for (it = route.begin(); it != route.end(); it++) {
    newints[n] = *it;
    n = n + 1;
  }

  seg1.resize(j + 1);
  seg2.resize(k - j);
  seg3.resize(length - k - 2);
  copy(newints, newints + j + 1, seg1.begin());
  copy(newints + j + 1, newints + k + 1, seg2.begin());
  copy(newints + k + 1, newints + length - 1, seg3.begin());

  reverse(seg1.begin(), seg1.end());
  reverse(seg2.begin(), seg2.end());

  newroute.insert(newroute.end(), seg1.begin(), seg1.end());
  newroute.insert(newroute.end(), seg2.begin(), seg2.end());
  newroute.insert(newroute.end(), seg3.begin(), seg3.end());

  it = find(newroute.begin(), newroute.end(), 0);
  newnewroute.insert(newnewroute.end(), it, newroute.end());
  newnewroute.insert(newnewroute.end(), newroute.begin(), it);

  delete[] newints;
  return newnewroute;
}

list<int> remove2(const list<int>& original_route, int l, int j, int k,
    int oldnode, int*myints) {

  list<int>::iterator it;
  list<int> seg1, seg2, seg3, seg4;
  list<int> route;
  list<int> newroute;
  list<int> newnewroute;

  list<int> oldroute;
  oldroute = original_route;

  int length = (int) oldroute.size();

  //order old route
  it = find(oldroute.begin(), oldroute.end(), oldnode);

  route.insert(route.end(), ++it, oldroute.end());
  route.insert(route.end(), oldroute.begin(), --it);

  int n = 0;
  int *newints;
  newints = new int[length - 1];
  for (it = route.begin(); it != route.end(); it++) {
    newints[n] = *it;
    n = n + 1;
  }

  seg1.resize(l + 1);
  seg2.resize(j - l);
  seg3.resize(k - j);
  seg4.resize(length - k - 2);
  copy(newints, newints + l + 1, seg1.begin());
  copy(newints + l + 1, newints + j + 1, seg2.begin());
  copy(newints + j + 1, newints + k + 1, seg3.begin());
  copy(newints + k + 1, newints + length - 1, seg4.begin());

  reverse(seg1.begin(), seg1.end());
  reverse(seg3.begin(), seg3.end());

  newroute.insert(newroute.end(), seg3.begin(), seg3.end());
  newroute.insert(newroute.end(), seg1.begin(), seg1.end());
  newroute.insert(newroute.end(), seg2.begin(), seg2.end());
  newroute.insert(newroute.end(), seg4.begin(), seg4.end());

  it = find(newroute.begin(), newroute.end(), 0);
  newnewroute.insert(newnewroute.end(), it, newroute.end());
  newnewroute.insert(newnewroute.end(), newroute.begin(), it);

  delete[] newints;
  return newnewroute;
}

list<int> insert1(const list<int>& original_route, int idi, int idj, int idk,
    int newnode, int*myints, int size, int category) {
  int i = -1;
  for (int l = 0; l < size; l++) {
    if (myints[l] == idi) {
      i = l;
    }
  }
  assert(i!=-1);

  assert(category != -1);

  int r;

  int *newints;
  int index = 0;
  newints = new int[size];

  list<int> oldroute;

  if (category == 0) {
    //does not need to change direction
    for (r = i; r < i + size; r = r + 1) {
      newints[index] = myints[(r + size) % size];
      index = index + 1;
    }
  } else {
    //need to change direction
    for (r = i; r > i - size; r = r - 1) {
      newints[index] = myints[(r + size) % size];
      index = index + 1;
    }

  }

  list<int> newroute, newnewroute;
  list<int>::iterator it;
  list<int> seg1, seg2;

  newroute.insert(newroute.begin(), idi);
  newroute.insert(newroute.end(), newnode);
  newroute.insert(newroute.end(), idj);

  index = 1;
  while (newints[index % size] != idj) {
    seg1.insert(seg1.end(), newints[index % size]);
    index = index + 1;
  }
  seg1.reverse();
  newroute.insert(newroute.end(), seg1.begin(), seg1.end());
  index = index + 1;

  newroute.insert(newroute.end(), idk);
  while (newints[index] != idk) {
    seg2.insert(seg2.end(), newints[index]);
    index = index + 1;
  }
  seg2.reverse();
  newroute.insert(newroute.end(), seg2.begin(), seg2.end());
  index = index + 1;

  while (index < size) {
    newroute.insert(newroute.end(), newints[index]);
    index = index + 1;
  }

  it = find(newroute.begin(), newroute.end(), 0);
  newnewroute.insert(newnewroute.end(), it, newroute.end());
  newnewroute.insert(newnewroute.end(), newroute.begin(), it);

  delete[] newints;

  return newnewroute;

}

list<int> insert2(const list<int>& original_route, int i, int l, int j, int k,
    int newnode, int *myints, int size, int category) {

  int idi, idj, idk, idl;
  int prel, prek, nextj, nexti;
  idi = myints[i];
  idj = myints[j];
  idk = myints[k];
  idl = myints[l];

  int r;

  int *newints;
  int index = 0;
  newints = new int[size];

  if (category == 0) {
    //does not need to change direction
    prel = myints[(l - 1 + size) % size];
    prek = myints[(k - 1 + size) % size];
    nexti = myints[(i + 1 + size) % size];
    nextj = myints[(j + 1 + size) % size];
    for (r = i; r < i + size; r = r + 1) {
      newints[index] = myints[(r + size) % size];
      index = index + 1;
    }
  } else {
    //needs to change direction
    prel = myints[(l + 1 + size) % size];
    prek = myints[(k + 1 + size) % size];
    nexti = myints[(i - 1 + size) % size];
    nextj = myints[(j - 1 + size) % size];
    for (r = i; r > i - size; r = r - 1) {
      newints[index] = myints[(r + size) % size];
      index = index + 1;
    }

  }

  list<int> newroute, newnewroute;
  list<int>::iterator it;
  list<int> seg1, seg2, seg3, seg4;

  int flag = 0; //flag to show if operated or not

  newroute.insert(newroute.begin(), idi);
  newroute.insert(newroute.end(), newnode);

  //seg1.resize(l-i-1);
  //seg2.resize(j-l+1);
  //seg3.resize(k-j-1);
  //seg4.resize(size-k+i+1);

  index = 1;

  while (newints[index % size] != idl) {
    seg1.insert(seg1.end(), newints[index % size]);
    index = index + 1;
  }
  seg1.reverse();

  while (newints[index % size] != nextj) {
    seg2.insert(seg2.end(), newints[index % size]);
    index = index + 1;
  }
  seg2.reverse();
  newroute.insert(newroute.end(), seg2.begin(), seg2.end());

  while (newints[index % size] != idk) {
    seg3.insert(seg3.end(), newints[index % size]);
    index = index + 1;
  }
  newroute.insert(newroute.end(), seg3.begin(), seg3.end());

  newroute.insert(newroute.end(), seg1.begin(), seg1.end());

  while (index < size) {
    newroute.insert(newroute.end(), newints[index]);
    index = index + 1;
  }

  it = find(newroute.begin(), newroute.end(), 0);
  newnewroute.insert(newnewroute.end(), it, newroute.end());
  newnewroute.insert(newnewroute.end(), newroute.begin(), it);

  delete[] newints;
  return newnewroute;
}

int find_order(int *myints, vector<int> v, int size) {
  vector<int> idv(v.size(), -1);

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < v.size(); j++) {
      if (myints[i] == v[j]) {
        idv[j] = i;
      }
    }
  }

  for (int i = 0; i < idv.size(); i++)
    assert(idv[i] != -1);

  idv.push_back(idv[0]);
  int last, bump;
  last = idv[0];
  bump = 0;
  for (int i = 1; i < idv.size(); i++) {
    if (last > idv[i]) {
      bump++;
    }
    last = idv[i];
  }

  if (bump <= 1)
    return 0;

  last = idv[0];
  bump = 0;
  for (int i = 1; i < idv.size(); i++) {
    if (last < idv[i]) {
      bump++;
    }
    last = idv[i];
  }

  if (bump <= 1)
    return 1;
  return -1;
}

/*This function performs both type I and type II GENI insertion*/
/*route: route before insertion;
 newnode: new node to insert;
 neighbor: neighborhood matrix;
 p: parameter to for nearest p-neighbor;
 arr_d: distance matrix*/
Path geni_insert(const Path& original_route, int newnode, int day, int vehicle,
    int p, double **arr_d, int *myints, int size) {

  double start = clock();
  double tempbest = DBL_MAX;
  double tempcost;
  Path newroute, temproute;
  list<int>::iterator it;
  int i, j, k, l;

  Path oldroute;
  oldroute = original_route;

  //oldroute.print();

  list<int>::const_iterator it2;
  int *next = new int[n_total + 1];
  int *prev = new int[n_total + 1];
  int last = 0;
  it2 = original_route.nodes.begin();
  it2++;
  for (; it2 != original_route.nodes.end(); it2++) {
    next[last] = *it2;
    last = *it2;
  }
  next[last] = 0;

  last = 0;
  it2 = original_route.nodes.end();
  it2--;
  for (; it2 != original_route.nodes.begin(); it2--) {
    prev[last] = *it2;
    last = *it2;
  }
  prev[last] = 0;

  //mar 07 2009
  int newid = newnode;
  int nexti, nextj, nextk;
  int newi, newj, newk, newl;
  int newnexti, newnextj, newprek, newprel;
  int newcate;
  int prek, prel;
  int idi, idj, idk, idl;
  int bestorder;
  double tempsave;
//	double oldcost;
  double diff;
//	oldcost = tcost(oldroute, new_d);
  double tempsavebest = 1000;
  bool flagchange = 0;
  bool flagfeasible = 0;
  bool category;

  if (size <= 2) {
    temproute = oldroute;

    // TODO: replace this temproute.nodes.insert(temproute.nodes.end(), newnode);
    temproute.nodes.insert(temproute.nodes.begin(), newnode);
    it = find(temproute.nodes.begin(), temproute.nodes.end(), 0);
    newroute.nodes.insert(newroute.nodes.end(), it, temproute.nodes.end());
    newroute.nodes.insert(newroute.nodes.end(), temproute.nodes.begin(), it);
    // TODO: end
    newroute.cost = tcost(newroute, new_d);
    newroute.demand = oldroute.demand + de_rank[newnode - 1][day];
  } else {
    //original
    for (i = 0; i < p && i < neighbours[day][vehicle][newnode].size(); i++) {
      int vi = neighbours[day][vehicle][newnode][i];
      for (j = 0; j < p && j < neighbours[day][vehicle][newnode].size();
          j = j + 1) {
        if (i == j)
          continue;
        int vj = neighbours[day][vehicle][newnode][j];

        //TODO: Use STL set to remove duplicates, improve performance.
        vector<int> tmpv;
        for (k = 0; k < p && k < neighbours[day][vehicle][next[vi]].size();
            k = k + 1) {
          tmpv.push_back(neighbours[day][vehicle][next[vi]][k]);
        }

        for (k = 0; k < p && k < neighbours[day][vehicle][prev[vi]].size();
            k = k + 1) {
          tmpv.push_back(neighbours[day][vehicle][prev[vi]][k]);
        }

        for (k = 0; k < 2 * p && k < tmpv.size(); k = k + 1) {
          int vk = tmpv[k];
          if (vk == vj || vk == vi)
            continue;

          vector<int> nodes;
          nodes.push_back(vi);
          nodes.push_back(vj);
          nodes.push_back(vk);
          int order = find_order(myints, nodes, size);
          assert(order != -1);

          if (order == 0) {
            tempsave = new_d[vi][newnode] + new_d[newnode][vj]
                - new_d[vi][next[vi]] - new_d[vj][next[vj]]
                + new_d[next[vi]][vk] - new_d[vk][next[vk]]
                + new_d[next[vj]][next[vk]];
          } else {
            tempsave = new_d[vi][newnode] + new_d[newnode][vj]
                - new_d[vi][prev[vi]] - new_d[vj][prev[vj]]
                + new_d[prev[vi]][vk] - new_d[vk][prev[vk]]
                + new_d[prev[vj]][prev[vk]];
          }

          if (DEBUG)
            cout << newnode << " " << vi << " " << vj << " " << vk << " "
                << tempsave << endl;

          if (tempsave < tempsavebest) {

            newi = vi;
            newj = vj;
            newk = vk;
            tempsavebest = tempsave;
            bestorder = order;

            flagchange = 1;
          }

        }
      }
    }

    if (flagchange == 1) {
      if (DEBUG) {
        cout << "here" << endl;
        cout << newnode << " " << newi << " " << newj << " " << newk << " "
            << tempsavebest << endl;
        system("pause");
      }
      newroute.nodes = insert1(oldroute.nodes, newi, newj, newk, newnode,
          myints, size, bestorder);
      if (DEBUG)
        cout << "Cost : " << oldroute.cost << " " << tempsavebest << " "
            << tcost(newroute, new_d) << endl;
      //assert(abs(oldroute.cost + tempsavebest - tcost(newroute, new_d)) < EP);
      newroute.cost = oldroute.cost + tempsavebest;
      newroute.demand = oldroute.demand + de_rank[newnode - 1][day];

    }

    else {
      newroute = oldroute;
    }

    // calculate reverce maping
    vector<int> route_index(n_total + 1, -1);
    for (i = 0; i < size; i++) {
      route_index[myints[i]] = i;
    }

    tempsavebest = DBL_MAX;
    flagchange = 0;
    if (size >= 5) {
      for (i = 0; i < p && i < neighbours[day][vehicle][newnode].size(); i++) {
        int vi = neighbours[day][vehicle][newnode][i];
        int ii = route_index[vi];
        for (j = 0; j < p && j < neighbours[day][vehicle][newnode].size();
            j++) {
          int vj = neighbours[day][vehicle][newnode][j];
          int jj = route_index[vj];
          if (vi == vj)
            continue;

          set<int> tmps1;
          for (k = 0; k < p && k < neighbours[day][vehicle][next[vi]].size();
              k++)
            tmps1.insert(neighbours[day][vehicle][next[vi]][k]);
          for (k = 0; k < p && k < neighbours[day][vehicle][prev[vi]].size();
              k++)
            tmps1.insert(neighbours[day][vehicle][prev[vi]][k]);

          set<int>::iterator it;
          for (it = tmps1.begin(); it != tmps1.end(); it++) {
            int vk = *it;
            int kk = route_index[vk];
            if (vk == vj)
              continue;

            set<int> tmps2;
            for (l = 0; l < p && l < neighbours[day][vehicle][next[vj]].size();
                l++)
              tmps2.insert(neighbours[day][vehicle][next[vj]][l]);
            for (l = 0; l < p && l < neighbours[day][vehicle][prev[vj]].size();
                l++)
              tmps2.insert(neighbours[day][vehicle][prev[vj]][l]);

            set<int>::iterator it2;
            for (it2 = tmps2.begin(); it2 != tmps2.end(); it2++) {
              int vl = *it2;
              int ll = route_index[vl];
              if (vl == vi)
                continue;

              /*
               idi = myints[i];
               idj = myints[j];
               idk = myints[k];
               idl = myints[l];
               */

              flagfeasible = 0;
              if (((ii < ll) && (ll < jj) && (jj < kk))
                  || ((ll < jj) && (jj < kk) && (kk < ii))
                  || ((jj < kk) && (kk < ii) && (ii < ll))
                  || ((kk < ii) && (ii < ll) && (ll < jj))) {
                category = 0;
                flagfeasible = 1;

                nexti = myints[(ii + 1 + size) % size];
                nextj = myints[(jj + 1 + size) % size];
                prek = myints[(kk - 1 + size) % size];
                prel = myints[(ll - 1 + size) % size];

              } else if (((kk < jj) && (jj < ll) && (ll < ii))
                  || ((jj < ll) && (ll < ii) && (ii < kk))
                  || ((ll < ii) && (ii < kk) && (kk < jj))
                  || ((ii < kk) && (kk < jj) && (jj < ll))) {
                category = 1;
                flagfeasible = 1;

                nexti = myints[(ii - 1 + size) % size];
                nextj = myints[(jj - 1 + size) % size];
                prek = myints[(kk + 1 + size) % size];
                prel = myints[(ll + 1 + size) % size];

              } else {
                continue;
              }

              if ((vk != vj) && (vk != nextj) && (vl != vi) && (vl != nexti)
                  && flagfeasible) {

                tempsave = new_d[newnode][vj] + new_d[nextj][vl]
                    + new_d[prek][prel] + new_d[vk][nexti] + new_d[vi][newnode]
                    - new_d[nextj][vj] - new_d[vk][prek] - new_d[vi][nexti]
                    - new_d[prel][vl];

                if (tempsave < tempsavebest) {
                  tempsavebest = tempsave;
                  newi = ii;
                  newj = jj;
                  newk = kk;
                  newl = ll;
                  newcate = category;
                  flagchange = 1;
                }
              }
            }
          }
        }
      }
    }

    if (flagchange == 1) {
      temproute.nodes = insert2(oldroute.nodes, newi, newl, newj, newk, newnode,
          myints, size, newcate);
      if (DEBUG)
        cout << "Cost : " << oldroute.cost + tempsavebest << " "
            << tcost(temproute, new_d) << endl;
      //assert(abs((oldroute.cost + tempsavebest - tcost(temproute, new_d))) < EP);
      temproute.cost = oldroute.cost + tempsavebest;
      temproute.demand = oldroute.demand + de_rank[newnode - 1][day];

      if (temproute.cost < newroute.cost) {
        newroute = temproute;
      }
    }

  }

  geni_insert_time += clock() - start;

  delete[] next;
  delete[] prev;

  return newroute;
}

Path geni_insert_cached(const Path& original_route, int newnode, int day,
    int vehicle, int p, double **arr_d, int *myints, int size) {
  if (insert_cache[newnode][vehicle][day] == NULL) {
    Path * path = new Path();
    *path = geni_insert(original_route, newnode, day, vehicle, p, arr_d, myints,
        size);
    insert_cache[newnode][vehicle][day] = path;
  }
  return *insert_cache[newnode][vehicle][day];
}

Path geni_remove(const Path& original_route, int newnode, int day, int p,
    double **arr_d, int *myints, int size) {
  double start = clock();
  double tempbest = DBL_MAX;
  double tempcost;
  Path newroute, temproute;
  list<int>::iterator it;
  int j, k, l;

  Path oldroute;
  oldroute = original_route;
  //oldroute.print();

  /*int *myints = NULL;
   myints = new int(oldroute.size());
   int n = 0;
   for(it = oldroute.begin(); it!=oldroute.end(); it++)
   {
   myints[n] = *it;
   n = n+1;
   }

   int size = n;*/

  int idj, idk, idi, idl;
  int nexti, nextj, nextk, nextl;

  idi = myints[0];
  int i = 0;
  while (idi != newnode) {
    idi = myints[i];
    i = i + 1;
  }
  idi = myints[i - 2];
  nexti = myints[(i + size) % size];

  double tempsave, tempsavebest;
  tempsavebest = DBL_MAX;

  int newl, newj, newk;

  if (size <= 3) {
    if (size <= 2) {
      newroute.nodes.clear();
      newroute.cost = 0;
      newroute.demand = 0;
    } else {
      newroute.nodes.insert(newroute.nodes.end(), idi);
      newroute.nodes.insert(newroute.nodes.end(), nexti);
      newroute.cost = tcost(newroute, new_d);
      newroute.calc_demand(day, de_rank);
    }
  } else {
    for (j = 0; j < size - 2; j = j + 1) {
      for (k = j + 1; k < size - 2; k = k + 1) {
        idj = myints[(j + i) % size];
        nextj = myints[(j + 1 + i) % size];
        idk = myints[(k + i) % size];
        nextk = myints[(k + 1 + i) % size];

        tempsave = new_d[idi][idj] + new_d[nexti][idk] + new_d[nextj][nextk]
            - new_d[idi][newnode] - new_d[newnode][nexti] - new_d[idj][nextj]
            - new_d[idk][nextk];

        if (tempsave < tempsavebest) {
          tempsavebest = tempsave;

          newj = j;
          newk = k;
        }
      }
    }

    newroute.nodes = remove1(oldroute.nodes, newj, newk, newnode, myints);
    if (DEBUG)
      cout << "Costs : " << oldroute.cost + tempsavebest << " "
          << tcost(newroute, new_d) << endl;
    //assert(abs(tempsavebest + oldroute.cost - tcost(newroute, new_d)) < EP);
    newroute.cost = tempsavebest + oldroute.cost;
    newroute.demand = oldroute.demand - de_rank[newnode - 1][day];

    if (size > 4) {
      tempsavebest = DBL_MAX;

      for (l = 0; l < size - 2; l = l + 1) {
        for (j = l + 1; j < size - 2; j = j + 1) {
          for (k = j + 1; k < size - 2; k = k + 1) {
            idl = myints[(l + i) % size];
            nextl = myints[(l + i + 1) % size];
            idj = myints[(j + i) % size];
            nextj = myints[(j + 1 + i) % size];
            idk = myints[(k + i) % size];
            nextk = myints[(k + 1 + i) % size];

            tempsave = new_d[nexti][nextl] + new_d[idj][nextk] + new_d[idi][idk]
                + new_d[nextj][idl] - new_d[newnode][nexti] - new_d[idl][nextl]
                - new_d[idj][nextj] - new_d[idk][nextk] - new_d[idi][newnode];

            if (tempsave < tempsavebest) {
              tempsavebest = tempsave;

              newl = l;
              newj = j;
              newk = k;
            }
          }
        }
      }
      temproute.nodes = remove2(oldroute.nodes, newl, newj, newk, newnode,
          myints);
      if (DEBUG)
        cout << "Costs : " << oldroute.cost + tempsavebest << " "
            << tcost(temproute, new_d) << endl;
      //assert(abs(oldroute.cost + tempsavebest - tcost(temproute, new_d)) < EP);
      temproute.cost = oldroute.cost + tempsavebest;
      temproute.demand = oldroute.demand - de_rank[newnode - 1][day];

      if (temproute.cost < newroute.cost) {
        newroute = temproute;
      }
    }

  }

  geni_remove_time += clock() - start;

  return newroute;
}

// geni_remove_cached(i+1, vtemp-1, r, route[r][vtemp-1], neighbor[r], p, new_d, myints, n);
Path geni_remove_cached(int newnode, int routeId, int day,
    const Path& original_route, int p, double **arr_d, int *myints, int size) {
  if (remove_cache[newnode][routeId][day] == NULL) {
    Path * path = new Path();
    *path = geni_remove(original_route, newnode, day, p, arr_d, myints, size);
    remove_cache[newnode][routeId][day] = path;
  }
  return *remove_cache[newnode][routeId][day];

}

