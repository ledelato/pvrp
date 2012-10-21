/*January 19th,2009 
 The third version of Tabu search for extended PVRP model
 with functions wrapped up*/

#include "geni.h"
#include "other.h"
#include "load.h"
#include <cstdlib>
#include <climits>
#include <csignal>
#include <ctime>

int root;

struct solution_change {
  int moved_node;
  int day;
  int old_route;
  int new_route;
};
//what solutions has changed from last time

typedef vector<vector<int> > assignment;

bool cmp(int a, int b) {
  return new_d[root][a] < new_d[root][b];
} //compare two value if a<b

void calculate_neighbor(vector<vector<Path> > &route) {
  neighbours.resize(d_total);
  for (int i = 0; i < d_total; i++) {
    neighbours[i].resize(v_total);
    for (int j = 0; j < v_total; j++) {
      //neighbours[i][j].resize(n_total);

      vector<int> nodes;

      for (int k = 0; k <= n_total; k++) {
        nodes.clear();
        list<int>::iterator it;
        for (it = route[i][j].nodes.begin(); it != route[i][j].nodes.end();
            it++) {
          if (*it != k && *it != -100)
            nodes.push_back(*it);
        }
        // insert 0 in empty neighbourhood
        if (nodes.empty())
          nodes.push_back(0);

        neighbours[i][j].push_back(nodes);
        root = k;
        sort(neighbours[i][j][k].begin(), neighbours[i][j][k].end(), cmp);
      }
    }
  }
}

void fix_neighbours(solution_change c) {

  if (c.moved_node == -1)
    return;

  // TODO: to improve this for speed

  //remove
  for (int i = 0; i <= n_total; i++) {
    for (int j = 0; j < (signed) neighbours[c.day][c.old_route][i].size(); j++)
      if (neighbours[c.day][c.old_route][i][j] == c.moved_node) {
        neighbours[c.day][c.old_route][i][j] =
            neighbours[c.day][c.old_route][i][neighbours[c.day][c.old_route][i].size()
                - 1];
        neighbours[c.day][c.old_route][i].pop_back();
        root = i;
        sort(neighbours[c.day][c.old_route][i].begin(),
            neighbours[c.day][c.old_route][i].end(), cmp);
      }
  }

  //insert
  for (int i = 0; i <= n_total; i++) {
    if (i != c.moved_node) {
      neighbours[c.day][c.new_route][i].push_back(c.moved_node);
      root = i;
      sort(neighbours[c.day][c.new_route][i].begin(),
          neighbours[c.day][c.new_route][i].end(), cmp);
    }
  }

}

void print_neighbours() {
  for (int d = 0; d < d_total; d++) {
    cout << "Day " << d << endl;

    for (int r = 0; r < v_total; r++) {
      cout << "\tRoute " << r << endl;

      for (int i = 0; i < n_total; i++) {
        cout << "\t\t";
        for (int j = 0; j < p && j < neighbours[d][r][i].size(); j++) {
          cout << " (" << neighbours[d][r][i][j] << ","
              << new_d[i][neighbours[d][r][i][j]] << ") ";
        }
        cout << endl;
      }
    }
  }
}

void print_all_routes(vector<vector<Path> > &route) {
  for (int i = 0; i < d_total; i++) {
    cout << "On day " << i << ":" << endl;
    for (int j = 0; j < v_total; j++) {
      //cout<<"For vehicle "<<j <<":"<<endl;
      route[i][j].print();
    }
  }
}

solution_change find_change(assignment &old_assign, assignment &new_assign) {
  solution_change sc;
  int num_changed = 0;

  for (int i = 0; i < n_total; i++) {
    for (int d = 0; d < d_total; d++) {
      if (old_assign[i][d] != new_assign[i][d]) {
        num_changed++;
        sc.day = d;
        sc.moved_node = i + 1;
        sc.new_route = new_assign[i][d] - 1;
        sc.old_route = old_assign[i][d] - 1;
      }
    }
  }

  assert(num_changed <= 1);
  if (num_changed == 0) {
    cout << "No change found" << endl;
    sc.moved_node = -1;
  }

  return sc;
}

/**
 * Sets the visit time of each customer from each vehicle to be zero
 */
void initialize_visit_time(vector<vector<int> > &visit_time) {
  visit_time.resize(n_total);

  for (int i = 0; i < n_total; i++) {
    visit_time[i].resize(v_total);
    for (int j = 0; j < v_total; j++) {
      visit_time[i][j] = 0;
    }
  }
}

/**
 * Sets the time that each vehicle visits each region to zero
 */
void initialize_region_time(vector<vector<int> > &region_time) {
  region_time.resize(l_total);

  for (int i = 0; i < l_total; i++) {
    region_time[i].resize(v_total);
    for (int j = 0; j < v_total; j++) {
      region_time[i][j] = 0;
    }
  }
}

/*
 * Set the frequency with which each driver visits each region on each day to zero
 */

void initialize_region_freq(vector<vector<vector<int> > > &region_freq) {
  region_freq.resize(d_total);

  for (int i = 0; i < d_total; i++) {
    region_freq[i].resize(v_total);
    for (int j = 0; j < v_total; j++) {
      region_freq[i][j].resize(l_total);
      for (int k = 0; k < l_total; k++) {
        region_freq[i][j][k] = 0;
      }
    }
  }
}

void print_matrix(vector<vector<int> > &matrix, int xsize, int ysize) {
  cout << "Printing matrix..." << endl;

  for (int i = 0; i < xsize; i++) {
    cout << "Node " << i << ": ";

    for (int j = 0; j < ysize; j++) {
      cout << " " << matrix[i][j];
    }
    cout << endl;
  }
}

void calculate_visit_time(vector<vector<Path> > &route,
    vector<vector<int> > &visit_time) {

  initialize_visit_time(visit_time);
  for (int i = 0; i < d_total; i++) {

    for (int j = 0; j < v_total; j++) {
      list<int>::iterator it;
      for (it = route[i][j].nodes.begin(); it != route[i][j].nodes.end();
          it++) {
        if (*it != 0)
          visit_time[*it - 1][j] = visit_time[*it - 1][j] + 1;
      }
    }
  }
}

void calculate_region_freq(vector<vector<Path> > &route) {
  initialize_region_freq(region_freq);

  if (DEBUG) {
    for (int i = 0; i < n_total; i++) {
      cout << reg[i] << " ";
    }
    cout << endl;
  }

  for (int i = 0; i < d_total; i++) {
    for (int j = 0; j < v_total; j++) {
      list<int>::iterator it;
      for (it = route[i][j].nodes.begin(); it != route[i][j].nodes.end();
          it++) {
        if (*it != 0)
          region_freq[i][j][reg[*it - 1]] += 1;
      }
    }
  }
}

void calculate_region_time(vector<vector<Path> > &route,
    vector<vector<int> > &region_time) {
  initialize_region_time(region_time);

  for (int i = 0; i < d_total; i++) {
    for (int j = 0; j < v_total; j++) {
      for (int r = 0; r < l_total; r++) {
        if (region_freq[i][j][r] != 0)
          region_time[r][j] = region_time[r][j] + 1;
      }
    }
  }

  //cout<<1<<endl;

}

double calculate_consistency_cost(vector<vector<int> > &visit_time) {
  double cost = 0;

  for (int i = 0; i < n_total; i++) {

    for (int j = 0; j < v_total; j++) {
      cost = cost + visit_time[i][j] * con_weight[visit_time[i][j]];
    }
  }

  return cost;
}

double calculate_region_cost(vector<vector<int> > &region_time) {
  double cost = 0;

  for (int i = 0; i < l_total; i++) {

    for (int j = 0; j < v_total; j++) {
      cost = cost + region_time[i][j] * reg_weight[region_time[i][j]];
    }
  }

  return cost;
}

void sighandler(int sig) {
  cout << "Program was trying to abort or terminate." << sig << endl;

  if (sig == 2) {
    cout << "Printing ... Exiting... " << endl;
    exit(0);
    /*
     time_t end = time(NULL);
     time_t last = end-start;


     flag_writeall = writedataout(routebest, objbest,last, lam-1, angle, div_factor, tabu_length, delta, p);
     flag_writerecord = writerecordout(objbest, last, lam, div_factor, tabu_length, delta, p);
     */
  }
}

int main() {

  signal(SIGABRT, &sighandler);
  signal(SIGTERM, &sighandler);
  signal(SIGINT, &sighandler);

  double total_time = 0;
  time_t start = time(NULL); /*starting time*/

  bool flag_read;
  bool flag_writeall;
  bool flag_writerecord;

  flag_read = readindata();

  for (int i = 0; i <= d_total; i++) {
    con_weight[i] *= alpha;
    cout << con_weight[i] << " ";
  }
  cout << endl;

  for (int i = 0; i <= d_total; i++) {
    reg_weight[i] = beta * reg_weight[i];
    cout << reg_weight[i] << " ";
  }
  cout << endl;

  vector<int> candidate;

  /*
   * Rho records the number of times that an attribute is selected as part of a solution
   */

  vector<vector<vector<int> > > rho(n_total + 4,
      vector<vector<int> >(d_total + 4, vector<int>(v_total, 0))); /*rho matrix*/

  /*
   * Tau records data for the tabu list
   */

  vector<vector<vector<int> > > tau(n_total + 4,
      vector<vector<int> >(d_total + 4, vector<int>(v_total, 0))); /*tau matrix*/

  /*
   * Sigma records the aspiration level for each attribute
   */

  vector<vector<vector<double> > > aspiration(n_total + 4,
      vector<vector<double> >(d_total + 4,
          vector<double>(v_total + 4, DBL_MAX))); /*sigma matrix*/

  /* routes for each vehicle on each day */
  vector<vector<Path> > route(d_total + 4, vector<Path>(v_total + 4));
  vector<vector<Path> > routetemp(d_total + 4, vector<Path>(v_total + 4));
  vector<vector<Path> > routebest(d_total + 4, vector<Path>(v_total + 4));
  vector<vector<Path> > routetempbest(d_total + 4, vector<Path>(v_total + 4));

  //TODO: change the type to assignment

  /* assignment of nodes to routes on each day */
  vector<vector<int> > assign(n_total + 4, vector<int>(d_total + 4, 0));

  vector<vector<int> > assigntemp(n_total + 4, vector<int>(d_total + 4, 0));
  vector<vector<int> > assigntempbest(n_total + 4, vector<int>(d_total + 4, 0));
  vector<vector<int> > assignbest(n_total + 4, vector<int>(d_total + 4, 0));
  vector<vector<double> > angle(n_total + 4, vector<double>(4, 0));

  int i, j, r;
  int fi, fj, fr, fk;
  int freq;
  int n;
  int lam;
  list<int>::iterator it;
  int vtemp = 0; /*temp vehicle number*/
  //int fday = 0;           /*the first day that the node is served*/
  int bestv; /*record the best vehicle for each node, vehicle id = bestv+1*/
  Path addroute, delroute; /*modified routes*/
  double diff = DBL_MAX; /*difference of best objectives in consequent iterations*/
  int stop_counter = 0;
  int totalserve = 0; /*total number of nodes needed to be served on each day*/
  int tabu_length = (int) floor(tabu_weight * log10((double) n_total));

  /*variables related to objective*/
  double objective = 0; /*objective*/
  double totalcost = 0; /*total cost-w/o capacity penalty*/
  //double travelcost = 0;        /*cost of entire routing*/
  double *daycost; /*cost on each day*/
  daycost = new double[d_total];

  double objbest = DBL_MAX; /*global best objective value*/
  double objtempbest = DBL_MAX; /*local best objective value*/
  double totaltempbest = DBL_MAX; /*local best travel cost*/
  double overcap_old = 0; /* the overcap of the last solution (s) */
  double overcap_temp_best = 0; /* the overcap of the best temp solution */
  double oldroute_total = 0; /* the cost of the last solution (s) */
  double objtemp, totaltemp = 0; /*record local objective and total cost*/
  double overcap_temp = 0; /* the overcap of the temp solution */
  double consist_temp = 0; /* record local customer consistency penalty */
  double region_temp = 0; /* record local region consistency penalty */

  double travelbest = 0;
  double consistbest = 0;
  double regionbest = 0;

  //double traveltemp = 0 ;         /*record local travel cost*/
  double *daytemp; /*record temporarily cost for each day*/
  daytemp = new double[d_total];

  /*record objective value of last iteration*/
  double objlast = DBL_MAX;

  int *myints = NULL;

  /*initial assignment--Sweep algorithm*/
  angle = rankNodes(); //rank the nodes according to angles
  assign = assign_initial(init_type); //make assignment

  initialize_visit_time(visit_time);
  best_visit_time = visit_time;

  /*initialize region times to zero*/
  initialize_region_time(region_time);
  initialize_region_time(best_region_time);

  initialize_region_freq(region_freq);

  /*build initial route*/
  route = route_initial(assign);

  for (i = 0; i < n_total; i++) {
    for (r = 0; r < d_total; r++) {
      vtemp = assign[i][r];
      if (vtemp != 0) {
        tau[i][r][vtemp - 1] = tau[i][r][vtemp - 1] + tabu_length; /*tabu tenure*/
        rho[i][r][vtemp - 1] = rho[i][r][vtemp - 1] + 1; /*frequency of choosing this route*/
      }
    }
  }

  calculate_neighbor(route);
  print_neighbours();

  calculate_visit_time(route, visit_time);
  double total_consist = calculate_consistency_cost(visit_time);

  calculate_region_freq(route);

  calculate_region_time(route, region_time);
  double total_region = calculate_region_cost(region_time);

  //assert(abs(total_consist-total_region)<EP);

  /*-------------------------------------TABU Search---------------------------------------------------*/

  /*-----------------------------initialization for tabu----------------------------------------------*/
  /*calculate the initial cost*/

  totalcost = cal_cost(route, assign);
  oldroute_total = totalcost;
  overcap_old = overcap(d_total, v_total, route);
  objective = cal_obj(totalcost, route);

  objtempbest = objective;
  objtemp = objective;

  routetempbest = route;
  routetemp = route;

  assigntempbest = assign;
  assigntemp = assign;

  /*initialization of parameters*/
  bool initfeas; /*feasibility of initial solution*/
  initfeas = isfeasible(assign, cap, de_rank, d_total, v_total, n_total);

  if (initfeas) /*if initial solution is feasible*/
  {
    /*reset aspiration level*/
    for (i = 0; i < n_total; i = i + 1) /*for every node*/
    {
      for (j = 0; j < d_total; j = j + 1) /*for every day*/
      {
        vtemp = assign[i][j];
        if (vtemp != 0)
          aspiration[i][j][vtemp - 1] = objbest;
      }
    }
  }

  /*------------------------------------------Searching----------------------------------------------*/
  for (lam = 1; lam < i_total; lam = lam + 1) {
    double starts = clock();
    cout << "lambda: " << lam << endl;

    objtempbest = DBL_MAX;

    //time_t end = time(NULL);          /*ending time*/   
    //time_t last = end-start;

    //flag_writeall = writedataout(routebest, objbest,last, lam-1, angle, div_factor, tabu_length, delta, p);
    //flag_writerecord = writerecordout(objbest, last, lam, div_factor, tabu_length, delta, p);

    if (DEBUG) {
      print_matrix(best_visit_time, n_total, v_total);
      cout << calculate_consistency_cost(best_visit_time) << endl;
      for (int vi = 0; vi <= d_total; vi++) {
        cout << con_weight[vi] << " ";
      }
      cout << endl;
      cout << alpha << endl;

      print_all_routes(route);
    }

    else {
      time_t start_iter = time(NULL); /*starting time*/

      for (i = 0; i < n_total; i = i + 1) /*for every node: node id = i+1*/
      {
        if (DEBUG) {
          cout << "node: " << i + 1 << endl;
        }

        for (r = 0; r < d_total; r = r + 1) /*operate on each day*/
        {
          if (DEBUG) {
            cout << "day r: " << r + 1 << endl;
          }

          /*----------------------------------------initialization----------------------------------*/
          /*current assigned vehicle*/
          bestv = assign[i][r];

          if (assign[i][r] != 0) {
            objtemp = objective;
            totaltemp = totalcost;

            /*-------------------------------------Find best vehicle------------------------------------*/

            for (j = 0; j < v_total; j = j + 1) /*for every vehicle: vehicle id = j+1*/
            {
              if (DEBUG) {
                cout << "vehicle: " << j + 1 << endl;
              }

              /*initialization temp variables*/
              assigntemp = assign;
              totaltemp = totalcost;
              routetemp = route;
              objtemp = objective;
              vtemp = assign[i][r];

              /*for the potential solutions*/
              if ((assign[i][r] != j + 1))
              /*(assign[i][r]!=j+1):not current solution;
               (isnear(route[k], (i+1), neighbor, p, q)): route j contains q nearest neighbor of current node*/
              {

                /*calculate travel cost*/
                if (!route[r][j].nodes.empty()) /*if route is not empty*/
                {
                  addroute.nodes.clear(); /*reinitialize*/
                  addroute.cost = 0;
                  addroute.demand = 0;
                  myints = new int[route[r][j].nodes.size()];
                  n = 0;
                  for (it = route[r][j].nodes.begin();
                      it != route[r][j].nodes.end(); it++) {
                    myints[n] = *it;
                    n = n + 1;
                  }
                  addroute = geni_insert_cached(route[r][j], (i + 1), r, j, p,
                      new_d, myints, n);

                  delete[] myints;

                  /*if insertion succeeds, then calculate cost*/
                  if (addroute.nodes.size() != route[r][j].nodes.size()) {
                    myints = new int[route[r][vtemp - 1].nodes.size()];
                    n = 0;
                    for (it = route[r][vtemp - 1].nodes.begin();
                        it != route[r][vtemp - 1].nodes.end(); ++it) {

                      myints[n] = *it;
                      n = n + 1;
                    }

                    delroute.nodes.clear();
                    delroute.cost = 0;
                    delroute.demand = 0;

                    delroute = geni_remove_cached(i + 1, vtemp - 1, r,
                        route[r][vtemp - 1], p, new_d, myints, n);

                    delete[] myints;

                    totaltemp = oldroute_total + addroute.cost + delroute.cost
                        - routetemp[r][j].cost - routetemp[r][vtemp - 1].cost;

                    overcap_temp =
                        overcap_old
                            - ((routetemp[r][j].demand > cap) ?
                                (routetemp[r][j].demand - cap) : 0)
                            - ((routetemp[r][vtemp - 1].demand > cap) ?
                                (routetemp[r][vtemp - 1].demand - cap) : 0)
                            + ((addroute.demand > cap) ?
                                (addroute.demand - cap) : 0)
                            + ((delroute.demand > cap) ?
                                (delroute.demand - cap) : 0);

                    /*create new route--routetemp*/
                    routetemp[r][j] = addroute;
                    routetemp[r][vtemp - 1] = delroute;
                    routetemp[r][vtemp - 1] = check(routetemp[r][vtemp - 1]);

                    /*create new assignment--assigntemp*/
                    assigntemp[i][r] = j + 1;

                    /*calculate cost*/
                    //assert(abs(totaltemp - cal_cost(routetemp, assigntemp)) < EP); 
                    //totaltemp = cal_cost(routetemp, assigntemp); 
                    consist_temp = total_consist
                        - visit_time[i][vtemp - 1]
                            * con_weight[visit_time[i][vtemp - 1]]
                        - visit_time[i][j] * con_weight[visit_time[i][j]]
                        + (visit_time[i][vtemp - 1] - 1)
                            * con_weight[visit_time[i][vtemp - 1] - 1]
                        + (visit_time[i][j] + 1)
                            * con_weight[visit_time[i][j] + 1];

                    region_temp =
                        total_region
                            + ((region_freq[r][vtemp - 1][reg[i]] == 1) ?
                                ((region_time[reg[i]][vtemp - 1] - 1)
                                    * reg_weight[region_time[reg[i]][vtemp - 1]
                                        - 1]
                                    - region_time[reg[i]][vtemp - 1]
                                        * reg_weight[region_time[reg[i]][vtemp
                                            - 1]]) :
                                0)
                            + ((region_freq[r][j][reg[i]] == 0) ?
                                ((region_time[reg[i]][j] + 1)
                                    * reg_weight[region_time[reg[i]][j] + 1]
                                    - region_time[reg[i]][j]
                                        * reg_weight[region_time[reg[i]][j]]) :
                                0);

                    //assert(abs(consist_temp-region_temp)<EP);

                    objtemp = totaltemp + consist_temp + region_temp
                        + cap_factor * overcap_temp;
                    if (DEBUG)
                      cout << "Obj cost : " << objtemp << " "
                          << cal_obj((totaltemp + consist_temp), routetemp)
                          << endl;

                    //assert(abs(objtemp - cal_obj((totaltemp+consist_temp), routetemp)) < EP );
                    //objtemp = cal_obj((totaltemp+consist_temp), routetemp); 

                  }
                } //end of if(candidate route is not empty)

                else /*if route is empty*/
                {
                  /*create a route with only node (i+1)*/
                  addroute.nodes.clear();
                  addroute.nodes.push_back(0);
                  addroute.nodes.push_back(i + 1);
                  addroute.cost = tcost(addroute, new_d);
                  addroute.calc_demand(r, de_rank);

                  myints = new int[route[r][vtemp - 1].nodes.size()];
                  n = 0;
                  for (it = route[r][vtemp - 1].nodes.begin();
                      it != route[r][vtemp - 1].nodes.end(); it++) {
                    myints[n] = *it;
                    n = n + 1;
                  }

                  delroute.nodes.clear();
                  delroute.cost = 0;
                  delroute.demand = 0;
                  delroute = geni_remove_cached(i + 1, vtemp - 1, r,
                      route[r][vtemp - 1], p, new_d, myints, n);

                  delete[] myints;

                  totaltemp = oldroute_total + addroute.cost + delroute.cost
                      - routetemp[r][j].cost - routetemp[r][vtemp - 1].cost;

                  overcap_temp = overcap_old
                      - ((routetemp[r][j].demand > cap) ?
                          (routetemp[r][j].demand - cap) : 0)
                      - ((routetemp[r][vtemp - 1].demand > cap) ?
                          (routetemp[r][vtemp - 1].demand - cap) : 0)
                      + ((addroute.demand > cap) ? (addroute.demand - cap) : 0)
                      + ((delroute.demand > cap) ? (delroute.demand - cap) : 0);

                  /*create new route--routetemp*/
                  routetemp[r][j] = addroute;
                  routetemp[r][vtemp - 1] = delroute;
                  routetemp[r][vtemp - 1] = check(routetemp[r][vtemp - 1]);

                  /*create new assignment--assigntemp*/
                  assigntemp[i][r] = j + 1;

                  /*calculate cost*/
                  //assert(abs(totaltemp - cal_cost(routetemp, assigntemp)) < EP);
                  //totaltemp = cal_cost(routetemp, assigntemp);  //objective cost--without capacity penalty  
                  consist_temp = total_consist
                      - visit_time[i][vtemp - 1]
                          * con_weight[visit_time[i][vtemp - 1]]
                      - visit_time[i][j] * con_weight[visit_time[i][j]]
                      + (visit_time[i][vtemp - 1] - 1)
                          * con_weight[visit_time[i][vtemp - 1] - 1]
                      + (visit_time[i][j] + 1)
                          * con_weight[visit_time[i][j] + 1];

                  if (DEBUG) {
                    cout << visit_time[i][vtemp - 1] << endl;
                    cout << region_time[reg[i]][vtemp - 1] << endl;
                    cout << visit_time[i][j] << endl;
                    cout << region_time[reg[i]][j] << endl;

                    for (int i = 0; i <= d_total; i++) {
                      cout << con_weight[i] << " ";
                    }
                    cout << endl;

                    for (int i = 0; i <= d_total; i++) {
                      cout << reg_weight[i] << " ";
                    }
                    cout << endl;
                  }

                  region_temp =
                      total_region
                          + ((region_freq[r][vtemp - 1][reg[i]] == 1) ?
                              ((region_time[reg[i]][vtemp - 1] - 1)
                                  * reg_weight[region_time[reg[i]][vtemp - 1]
                                      - 1]
                                  - region_time[reg[i]][vtemp - 1]
                                      * reg_weight[region_time[reg[i]][vtemp - 1]]) :
                              0)
                          + ((region_freq[r][j][reg[i]] == 0) ?
                              ((region_time[reg[i]][j] + 1)
                                  * reg_weight[region_time[reg[i]][j] + 1]
                                  - region_time[reg[i]][j]
                                      * reg_weight[region_time[reg[i]][j]]) :
                              0);

                  //assert(abs(consist_temp-region_temp)<EP);

                  objtemp = totaltemp + consist_temp + region_temp
                      + cap_factor * overcap_temp;

                  //assert(abs(objtemp - cal_obj((totaltemp+consist_temp), routetemp)) < EP );

                  //objtemp = cal_obj((totaltemp+consist_temp), routetemp);

                } //end of if(candidate route is empty)

                /*Diversification*/
                if (objtemp >= objective) {
                  /*calculate frequency*/
                  freq = 0;
                  double start_fre = clock();
                  for (fi = 0; fi < n_total; fi = fi + 1) {
                    for (fr = 0; fr < d_total; fr = fr + 1) {
                      fj = assigntemp[fi][fr];
                      fk = assign[fi][fr];
                      if ((fj != 0) && (fj != fk))
                        freq = freq + rho[fi][fr][fj - 1];
                    }
                  }
                  frequency_time += clock() - start_fre;

                  /*add penalty cost*/
                  objtemp = objtemp
                      + div_factor * (sqrt(double(n_total * v_total * d_total)))
                          * freq
                          * (oldroute_total + total_consist + total_region)
                          / lam;

                } //no penalty if not entering this if statement

                if ((tau[i][r][j] < lam)
                    || (tau[i][r][j] >= lam
                        && isfeasible(assigntemp, cap, de_rank, d_total,
                            v_total, n_total) && aspiration[i][r][j] > totaltemp)) {
                  if ((objtemp < objtempbest) && (objtemp != 0)) {
                    objtempbest = objtemp;
                    totaltempbest = totaltemp;
                    routetempbest = routetemp;
                    assigntempbest = assigntemp;
                    overcap_temp_best = overcap_temp;
                  }
                }
              } //end of if(route j is a candidate solution)
            } //end of if (vehicle j is not tabued)
          } //end of for(every vehicle)
        } //end of if(the solution is different from current solution
      } //end of each day r

      /*-------------------------------------Update------------------------------------------------*/

      /*transition to the best solution*/

      solution_change change;
      change = find_change(assign, assigntempbest);
      if (1 || DEBUG) {
        cout << "Change " << change.moved_node << " " << change.old_route << " "
            << change.new_route << " " << change.day << endl;
        //  system("pause");
      }

      //change consistency
      total_consist = total_consist
          - visit_time[change.moved_node - 1][change.old_route]
              * con_weight[visit_time[change.moved_node - 1][change.old_route]]
          - visit_time[change.moved_node - 1][change.new_route]
              * con_weight[visit_time[change.moved_node - 1][change.new_route]]
          + (visit_time[change.moved_node - 1][change.old_route] - 1)
              * con_weight[visit_time[change.moved_node - 1][change.old_route]
                  - 1]
          + (visit_time[change.moved_node - 1][change.new_route] + 1)
              * con_weight[visit_time[change.moved_node - 1][change.new_route]
                  + 1];

      visit_time[change.moved_node - 1][change.old_route] =
          visit_time[change.moved_node - 1][change.old_route] - 1;
      visit_time[change.moved_node - 1][change.new_route] =
          visit_time[change.moved_node - 1][change.new_route] + 1;

      if (region_freq[change.day][change.old_route][reg[change.moved_node - 1]]
          == 1) {
        total_region =
            total_region
                + (region_time[reg[change.moved_node - 1]][change.old_route] - 1)
                    * reg_weight[region_time[reg[change.moved_node - 1]][change.old_route]
                        - 1]
                - region_time[reg[change.moved_node - 1]][change.old_route]
                    * reg_weight[region_time[reg[change.moved_node - 1]][change.old_route]];

        region_time[reg[change.moved_node - 1]][change.old_route] =
            region_time[reg[change.moved_node - 1]][change.old_route] - 1;
      }

      if (region_freq[change.day][change.new_route][reg[change.moved_node - 1]]
          == 0) {
        total_region =
            total_region
                + (region_time[reg[change.moved_node - 1]][change.new_route] + 1)
                    * reg_weight[region_time[reg[change.moved_node - 1]][change.new_route]
                        + 1]
                - region_time[reg[change.moved_node - 1]][change.new_route]
                    * reg_weight[region_time[reg[change.moved_node - 1]][change.new_route]];

        region_time[reg[change.moved_node - 1]][change.new_route] =
            region_time[reg[change.moved_node - 1]][change.new_route] + 1;
      }

      region_freq[change.day][change.old_route][reg[change.moved_node - 1]] =
          region_freq[change.day][change.old_route][reg[change.moved_node - 1]]
              - 1;
      region_freq[change.day][change.new_route][reg[change.moved_node - 1]] =
          region_freq[change.day][change.new_route][reg[change.moved_node - 1]]
              + 1;

      assert(abs(total_region-calculate_region_cost(region_time))<EP);

      //invalidate the cache
      for (int vi = 1; vi <= n_total; vi++) {
        // the remove cache
        delete remove_cache[vi][change.old_route][change.day];
        remove_cache[vi][change.old_route][change.day] = NULL;
        delete remove_cache[vi][change.new_route][change.day];
        remove_cache[vi][change.new_route][change.day] = NULL;

        // the insert cache
        delete insert_cache[vi][change.old_route][change.day];
        insert_cache[vi][change.old_route][change.day] = NULL;
        delete insert_cache[vi][change.new_route][change.day];
        insert_cache[vi][change.new_route][change.day] = NULL;
      }

      tau[change.moved_node - 1][change.day][change.old_route] = lam
          + tabu_length;
      rho[change.moved_node - 1][change.day][change.new_route]++;

      assign = assigntempbest;
      route = routetempbest;

      fix_neighbours(change);

      //objtempbest = cal_cost(route, assign);
      //assert(abs(totaltempbest - cal_cost(route, assign)) < EP);
      oldroute_total = totaltempbest;
      overcap_old = overcap_temp_best;
      //objtempbest = cal_cost(route, assign) + total_consist;
      objtempbest = totaltempbest + total_consist + total_region;

      cout << "Best Cost: " << objbest << endl;

      cout << "Current Cost: " << objtempbest << endl;
      cout << "Current Travel Cost: " << totaltempbest << endl;
      cout << "Current Consistency Cost: " << total_consist << endl;
      cout << "Current Regional Cost: " << total_region << endl;

      if (isfeasible(assign, cap, de_rank, d_total, v_total, n_total))/*check if the solution is feasible*/
      {

        cap_factor = cap_factor / (1 + delta);
        //isfeasible(assign, cap, de_rank, d_total, v_total, n_total);

        /*record solution*/
        if (objtempbest < objbest) {
          objbest = objtempbest;
          travelbest = totaltempbest;
          consistbest = total_consist;
          regionbest = total_region;
          routebest = routetempbest;
          assignbest = assigntempbest;
          best_visit_time = visit_time;
          best_region_time = region_time;
        }

        /*update aspiration level*/
        for (fi = 0; fi < n_total; fi = fi + 1) {
          for (fr = 0; fr < d_total; fr = fr + 1) {
            fj = assign[fi][fr];
            if (fj != 0) {
              //cout<<aspiration[1][0][0];
              if (objtempbest < aspiration[fi][fr][fj - 1]) {
                aspiration[fi][fr][fj - 1] = objtempbest;
              }
            }
          }
        }

      } //end of if (the solution is feasible)
      else {
        cap_factor = cap_factor * (1 + delta);
      } //end of if (the solution is not feasible)

      diff = abs(objbest - objlast);
      objlast = objbest;

      /*debug information*/
      //time_t end_iter = time(NULL);         /*ending time*/   
      //time_t last_iter = end_iter-start_iter;
    } //end for stopping criteria isn't met

    total_time += clock() - starts;
    cout << "Iteration Time: " << (clock() - starts) / CLOCKS_PER_SEC << endl;
    cout << "Total Time: " << total_time / CLOCKS_PER_SEC << endl;

  } //end for each lambda

  //time_t end = time(NULL);          /*ending time*/
  //time_t last = end-start;

  //calculate statistics

  total_time = total_time / CLOCKS_PER_SEC;
  con_weight = con_weight_copy;
  consistbest = calculate_consistency_cost(best_visit_time);
  regionbest = calculate_region_cost(best_region_time);
  flag_writeall = writedataout(routebest, objbest, angle, tabu_length,
      total_time, travelbest, consistbest, regionbest);
  flag_writerecord = writerecordout(objbest, tabu_length, total_time,
      travelbest, consistbest, regionbest);
  //print_all_routes(route);

  free(arr_d); /*pointer to the array of distance matrix*/
  free(new_d); /*pointer to the ordered distance matrix*/
  free(co_old); /*coordinate date input - raw*/
  free(co); /*ordered coordinate*/
  free(arr_reg); /*pointer to the array of region information*/
  free(demand); /*pointer to the array of demand matrix*/

  return 0;
}
