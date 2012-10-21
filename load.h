#include <stdio.h>
#include <stdlib.h>

/*This function reads in the data from data file
 * returns
 *   1 - cannot open input file
 *   0 - file opened successfully
 *
 * */
bool readindata() {
  FILE *file;
  int i, j;

  file = fopen(name, "r");
  //filew = fopen("debug.txt", "w");

  if (file == NULL) {
    printf("Error: can't open file.\n");
    return 1;
  } else {
    printf("File opened successfully.\n");
  }

  fscanf(file, "%d", &n_total);
  fscanf(file, "%d", &d_total);
  fscanf(file, "%d", &l_total);
  fscanf(file, "%d", &v_total);
  fscanf(file, "%d", &cap);

  arr_d = new double *[n_total + 3];

  /*read in distance matrix*/
  float temp;
  int temp_d;

  /*for (i = 0;i<=n_total; i = i+1)
   {
   //cout<<"d[i][]: ";
   fprintf(filew, "d[%d][]: ", i);
   arr_d[i] = new double[n_total+4];
   for (j = 0;j<=n_total; j = j+1)
   {
   fscanf(file, "%f", &temp);
   arr_d[i][j] = temp;
   //cout<<arr_d[i][j]<<" ";
   fprintf(filew, "%f ", arr_d[i][j]);
   }
   //cout<<endl;
   fprintf(filew, "\n");
   }*/

  /*read in coordinate matrix*/
  co_old = new double *[n_total + 4];

  for (i = 0; i <= n_total; i = i + 1) {
    //cout<<"co[i][]: ";
    //fprintf(filew, "co[%d][]: ", i);
    co_old[i] = new double[3];

    fscanf(file, "%f", &temp);
    co_old[i][0] = temp;
    //cout<<co_old[i][0]<<" ";
    //fprintf(filew, "%f ", co_old[i][0]);
    fscanf(file, "%f", &temp);
    co_old[i][1] = temp;
    //cout<<co_old[i][1]<<" ";
    //fprintf(filew, "%f ", co_old[i][1]);

    //cout<<endl;
    //fprintf(filew, "\n");
  }

  /*read in region information*/
  arr_reg = new int[n_total + 3];

  for (i = 0; i < n_total; i = i + 1) {
    fscanf(file, "%d", &temp_d);
    arr_reg[i] = temp_d;
    //cout<<arr_reg[i]<<endl;
    //fprintf(filew, "%d ", arr_reg[i]);
  }

  /*read in demand matrix*/
  demand = new double *[n_total + 3];

  for (i = 0; i < n_total; i = i + 1) {
    //cout<<"demand[i][]: ";
    //fprintf(filew, "demand[%d][]:: ", i);
    demand[i] = new double[d_total + 3];
    for (j = 0; j < d_total; j = j + 1) {
      fscanf(file, "%f", &temp);
      demand[i][j] = temp;
      //cout<<demand[i][j]<<" ";
      //	fprintf(filew, "%f ", demand[i][j]);
    }

    //fprintf(filew, "\n");
    //cout<<endl;
  }

  //fclose(filew);

  fclose(file);

  /*read in parameter*/

  file = fopen(paraname, "r");

  if (file == NULL) {
    printf("Error: can't open file.\n");
    return 1;
  } else {
    printf("File opened successfully.\n");
  }

  con_weight = new double[d_total + 1];
  con_weight_copy = new double[d_total + 1];
  cout << "con_weight" << endl;
  for (i = 0; i <= d_total; i = i + 1) {
    fscanf(file, "%f", &temp);
    con_weight[i] = temp;
    con_weight_copy[i] = temp;

    cout << con_weight[i] << " " << endl;
  }

  reg_weight = new double[d_total + 1];
  reg_weight_copy = new double[d_total + 1];

  cout << "reg_weight" << endl;
  for (i = 0; i <= d_total; i = i + 1) {
    fscanf(file, "%f", &temp);
    reg_weight[i] = temp;
    reg_weight_copy[i] = temp;

    cout << reg_weight[i] << " " << endl;
  }

  fclose(file);

  return 0;
}

/*This function write the total result into result file*/
bool writedataout(vector<vector<Path> > route, double cost,
    vector<vector<double> > angle, int tabu_length, double total_time,
    double travelbest, double consistbest, double regionbest) {
  assert(abs(cost - travelbest - consistbest - regionbest) < EP);
  FILE *file;
  int i, j;
  list<int>::iterator it;

  int node;

  file = fopen(routename, "a");

  if (file == NULL) {
    printf("Error: can't open file.\n");
    return 1;
  }

  fprintf(file, name);
  fprintf(file, "\n\n");

  int index;
  int start = 0;
  int end = 0;

  for (i = 0; i < n_total + 3; i = i + 1) {
    if (angle[i][1] == 0) {
      end = i;
    }
  }
  start = end - 2;

  for (i = 0; i < d_total; i = i + 1) {
    //fprintf(file,"On Day %d:\n", i+1);
    for (j = 0; j < v_total; j = j + 1) {
      if (!route[i][j].nodes.empty()) {
        //fprintf(file,"Vehicle %d: ", j+1);
        for (it = route[i][j].nodes.begin(); it != route[i][j].nodes.end();
            it++) {
          index = int(*it) - 1;
          if (index < 0) {
          } else {
            if (index < start) {
              node = (int) angle[index][1];
              fprintf(file, "%d ", node);
            }
            if (index >= start) {
              index = index + end - start + 1;
              node = (int) angle[index][1];
              fprintf(file, "%d ", node);
            }
          }
        }
        fprintf(file, "%d\n", -100);
      } else {
        fprintf(file, "%d\n", -100);
      }
    }
    fprintf(file, "\n");
  }

  fprintf(file, "div_factor = %f  \n", div_factor);
  fprintf(file, "tabu length = %d \n", tabu_length);
  fprintf(file, "delta =  %f \n", delta);
  fprintf(file, "alpha =  %f \n", alpha);
  fprintf(file, "beta =  %f \n", beta);
  fprintf(file, "neighborhood size = %d \n", p);
  fprintf(file, "%d iterations are performed. \n", i_total);
  fprintf(file, "The total cost is %f. \n", cost);
  fprintf(file, "The total time is %f seconds. \n", total_time);
  fprintf(file, "The travel cost is %f. \n", travelbest);
  fprintf(file, "The node consistency cost is %f. \n", consistbest);
  fprintf(file, "The region complexity cost is %f. \n\n\n\n", regionbest);

  fclose(file);

  return 0;

}

bool writerecordout(double cost, int tabu_length, double total_time,
    double travelbest, double consistbest, double regionbest) {
  assert(abs(cost - travelbest - consistbest - regionbest) < EP);
  FILE *file;
  file = fopen(recordname, "a");

  if (file == NULL) {
    printf("Error: can't open file.\n");
    return 1;
  }

  fprintf(file, name);
  fprintf(file, "\n");

  //fprintf(file, "%f  ", div_factor);
  fprintf(file, "%d  ", tabu_length);
  //fprintf(file, "%f  ", delta);
  //fprintf(file, "%d  ", p);
  fprintf(file, "%f  ", cost);
  fprintf(file, "%f  ", travelbest);
  fprintf(file, "%f  ", consistbest);
  fprintf(file, "%f  ", regionbest);
  //fprintf(file, "%d  ", lam);
  fprintf(file, "%f  \n\n\n", total_time);

  fclose(file);

  return 0;
}
/*This function output the result at current stage*/
