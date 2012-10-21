
/*This function calculates the penalty for consistency*/
double consistency (int n_total, int v_total, double alpha, int **con_r)
{
	int i, j, k;
	double cost = 0;

	for(i = 0; i<n_total; i = i+1)	/*each node*/
	{
		for(j=0; j<v_total; j = j+1)	/*each vehicle*/
		{
			k = con_r[i][j];
			cost = cost+k*alpha*con_weight[k];
		}
	}

	return cost;
}

/*This function calculates the penalty for region familiarity*/
double region (int l_total, int v_total, double beta, int **reg_s)
{
       cout<<"region--begin"<<endl;
	int i, j, k;
	double cost = 0;

	for(i = 0; i<l_total; i = i+1)	/*each region*/
	{
		for(j=0; j<v_total; j = j+1)	/*each vehicle*/
		{
			k = reg_s[i][j];
			cost = cost+k*beta*reg_weight[k];
			//cout<<"i "<<i<<" "<<"j "<<j<<endl;
			//cout<<"reg_weight "<<reg_weight[k]<<endl;
		}
	}
	//cout<<"region--end"<<endl;

	return cost;
}


/*This function calculates the penalty for over capacity*/
double overcap (int d_total, int v_total, const vector<vector<Path> >& route)
{
	double detemp = 0;
	double over = 0;		
	list<int>::iterator it;
	Path temproute;
	int vtemp;
	
	int r, k;
	for(r = 0; r<d_total; r = r+1)		/*each day*/
	{
		for(k = 0; k<v_total; k = k+1)	/*for each vehicle*/
		{
			if(!route[r][k].nodes.empty())	/*if the route is not empty*/
			{
				detemp = 0;
				temproute = route[r][k];
				for(it = temproute.nodes.begin(); it!=temproute.nodes.end(); it++)
				{	
					vtemp = *it;
					if(vtemp!=0)/*except for depot*/
					detemp = detemp+de_rank[vtemp-1][r];
				}
				if(detemp>cap) 
				{
					over = over+detemp-cap;
				}
				else 
				{
					over = over;
				}
			}
		}
	}

	//temproute.~list();

	return over;
}


/*This function calculates the total cost*/
double cal_cost(const vector<vector<Path> >& route, const vector<vector<int> >& assign)
{
	double start = clock();
	double travelcost = 0;
	double concost, regcost, totalcost = 0;
	Path temproute;
	int i, j, r, k;


	/*travel cost*/
	for (i = 0; i<d_total; i = i+1)			//for every day
	{
		//daycost[i] = 0;					//initialization for everyday cost
		for (j=0; j<v_total; j = j+1)		//for every route
		{
			if(!route[i][j].nodes.empty())	//if the route is not empty
			{
				temproute = route[i][j];
				travelcost = travelcost+tcost(temproute, new_d);
			}
		}
	} 
	

	/*int** con_r;						//consistency r matrix
	con_r = new int *[n_total+3];
	for(i = 0; i<n_total; i = i+1)
	{
		con_r[i] = new int[v_total+3];

		for(j = 0; j<v_total; j=j+1)
		{
			con_r[i][j] = 0;
		}
	}*/
	
	//cout<<"r matrix"<<endl;

	/*int** reg_s;						//consistency s matrix
	reg_s = new int *[l_total+3];
	for(i = 0; i<l_total; i = i+1)
	{
		reg_s[i] = new int[v_total+3];
		for(j = 0; j<v_total; j=j+1)
		{
			reg_s[i][j] = 0;
		}
	}*/

	//cout<<sizeof(route)<<" "<<sizeof(travelcost)<<endl;
    //cout<<"s matrix"<<endl;

	

	/*calculate r_ik_n and s_lk_n*/
	/*for(i = 0; i<n_total; i = i+1)	//for each node
	{
		for(r = 0; r<d_total; r = r+1)	// for each day
		{
			j = assign[i][r]-1;			//assigned vehicle
			if(j>=0)					//if j<0, assign = 0
			{
				k = arr_reg[i]-1;				//the region node i lies in
				con_r[i][j] = con_r[i][j]+1;
				reg_s[k][j] = reg_s[k][j]+1;
			}
		}
	}*/
	
	//cout<<"r_ik_n"<<endl;
	/*calculate penalty*/
	//concost = consistency(n_total, v_total, alpha, con_r);		/*consistency*/
	//cout<<"concost"<<endl;
	//regcost = region(l_total, v_total, beta, reg_s);			/*region familiarity*/
	//cout<<"regcost"<<endl;

	/*total cost-objective*/	

	totalcost = travelcost;//+concost+regcost;
	
	//temproute.~list();
	cal_cost_time += clock() - start;

	return totalcost;
}


/*This function calculates the total objective value*/
double cal_obj(double totalcost, const vector<vector<Path> >& oldroute)
{
	double start = clock();
	double objective = 0;

	 vector<vector<Path> > route;
	 route = oldroute;

	objective = totalcost+cap_factor*overcap(d_total, v_total, route);

	if(objective ==0)
	{
		cout<<"error!!!"<<endl;
		exit(1);
	}

	//route.~vector();
	cal_obj_time += clock() - start;
	return objective;
}

/*This function checks if the newly created route is "legal" and do corresponding corrections if it is not */
vector<vector<list<int> > >  check_all (const vector<vector<list<int> > >& route, int d_total, int v_total)
{
	vector<vector<list<int> > > route_new(d_total+3, vector<list<int> >(v_total+3, list<int>(4,-100)));
	route_new = route;

	int i, j;
	list<int>::iterator it;

	for(i = 0; i<d_total; i = i+1)	/*for every day*/
	{
		for(j = 0; j<v_total; j = j+1)	/*for every vehicle*/
		{
			if (route_new[i][j].begin() == route_new[i][j].end())
			{
				route_new[i][j].push_back(-100);

			}
			else 
			{
				it = route_new[i][j].begin();
				++it;
				if (it==route_new[i][j].end())
				{
					route_new[i][j].clear();
					route_new[i][j].push_back(-100);
				}

			}

		}
	}
	return route_new;
}


/*This function checks if one route is "legal".*/
Path  check (const Path &  route)
{
	Path route_new;
	route_new = route;

	list<int>::iterator it;

	
	if (!route_new.nodes.empty())
	{
		it = route_new.nodes.begin();
		++it;
		if (it==route_new.nodes.end())
		{
			route_new.nodes.clear();
		}

	}
	return route_new;
}

vector<vector<Path> > route_initial(const vector<vector<int> >& assign)
{
	int i, j;
	int vtemp = 0;

	vector<vector<Path > > route (d_total+3, vector<Path>(v_total+3));

	for(i = 0; i<d_total; i = i+1)		/*on each day*/
	{
		for(j = 0; j<n_total; j = j+1)	/*for each node*/
		{
			
			if(assign[j][i]!=0)	/*if the node is visited*/
			{
				vtemp = assign[j][i];
				if(route[i][vtemp-1].nodes.empty())	/*if j is the first node in the route*/
				{
					route[i][vtemp-1].nodes.clear();
					route[i][vtemp-1].nodes.push_back(0);		/*every route starts with depot*/
					route[i][vtemp-1].nodes.push_back(j+1);	/*node id is (j+1)*/
				}
				else route[i][vtemp-1].nodes.push_back(j+1);
			}
		}
	}

	for(i = 0; i<d_total; i = i+1)		/*on each day*/
	{
		/* calculate initial cost and demand*/
		for(j = 0; j<v_total; j = j+1)
		{
			route[i][j].cost = tcost(route[i][j], new_d);
			route[i][j].calc_demand(i, de_rank);
		}
	}

	return route;
}

vector<vector<int> > assign_initial(int type)
{
	
	vector<vector<int> > assign(n_total+3, vector<int>(d_total+3,0));
	
	double *cum_de;			/*pointer to the array of cumulative demand*/
	cum_de = new double [n_total+3];

	double *avg_de;
	avg_de = new double [n_total+3];

	int curr_v = 1;		/*current vehicle id*/
	double fill = cap;	/*left capacity of current vehicle*/

	int i, j;

	/*calculate the average demand for each customer*/
	//cout<<"cum_de"<<endl;
	for(i = 0; i<n_total; i = i+1)
	{
		cum_de[i]=0;
		for (j = 0; j<d_total; j = j+1)
		{
			cum_de[i] = cum_de[i]+de_rank[i][j];
		}
		//cout<<i <<" "<<cum_de[i]<<endl;
		avg_de[i] = cum_de[i]/d_total;	/*average task for each day*/
	}

	

	/*initial assignment--Sweep algorithm*/
	if(type==1)			//remain consistency
	{
		for(i = 0;i < n_total;i = i+1)	/*for each node*/
		{	
			if(curr_v<v_total)		/*before the last vehicle*/
			{
				if(fill>0)
				{
					for(j = 0;j < d_total;j = j+1)		/*on each day*/
					{
						if(de_rank[i][j]!=0) assign[i][j] = curr_v;		/*assign current vehicle to node--remain consistency*/
					}
					fill = fill-avg_de[i];
					cout<<avg_de[i]<<endl;
				}
				else
				{
					curr_v = curr_v+1;		/*start to assign to another vehicle*/
					fill = cap;
					for(j = 0;j < d_total;j = j+1)
					{
						if(de_rank[i][j]!=0) 
							assign[i][j] = curr_v;		/*assign current vehicle to node--remain consistency*/
					}
					fill = fill-avg_de[i];
				}
			}
			else			/*if reaches the last vehicle*/
			{
				for(j = 0;j<d_total;j = j+1)
				{
					if(de_rank[i][j]!=0) assign[i][j] = curr_v;		/*assign current vehicle to node--remain consistency*/
				}
			}
		}
			/*print the initial assignment*/
			//cout<<i<<" "<<assign[i][0]<<" "<<fill<<endl;
	}
	else
	{
		for(int j = 0; j<d_total; j++){
			//initialization
			fill = cap;
			curr_v = 1;

			for(int i = 0; i<n_total; i++){
				if(curr_v<v_total){
					if(fill>0){
						if(de_rank[i][j]!=0) 
							assign[i][j] = curr_v;	
						fill = fill - de_rank[i][j];
					}
					else{
						curr_v += 1;
						fill = cap;
						if(de_rank[i][j]!=0) 
							assign[i][j] = curr_v;	
						fill = fill - de_rank[i][j];
					}
				}
				else{
					if(de_rank[i][j]!=0) 
							assign[i][j] = curr_v;	
				}
			}
			/*print the initial assignment*/
			//cout<<i<<" "<<assign[i][0]<<" "<<fill<<endl;
		}

	}

	delete []cum_de;
	delete []avg_de;

	return assign;
}

/*This function determines if the current solution is not the same as the original one. 
given that the "node" is the only node that has changed its assignment*/
bool isnew(const vector<vector<int> > &assigntemp, const vector<vector<int> >& assign, int node)
{	
	int i;
	for(i = 0; i<d_total; i = i+1)
	{
		if(assigntemp[node][i]!=assign[node][i])
			return 1;
	}

	return 0;

}


double sum(const vector<vector<vector<int> > >& rho)
{
	int i, j, k;
	double total = 0;

	for(i = 0; i<n_total; i = i+1)
	{
		for (j = 0; j<d_total; j = j+1)
		{
			for(k = 0; k<v_total; k = k+1)
			{
				total = total+rho[i][j][k];
			}
		}
	}

	return total;
}

vector<vector<double> > rankNodes()
{

	vector<vector<double> > angle(n_total+3, vector<double>(4,0));

	//FILE *file;

	int i = 0;
	int j = 0;
	
	
	
	/*for (i = 0;i<n_total; i = i+1)
	{
		cout<<"demand[i][]: ";
		for (j = 0;j<d_total; j = j+1)
		{
			cout<<demand[i][j]<<" ";
		}
		cout<<endl;	
	}*/
	
	
	


	for (i = 1; i<=n_total; i = i+1)
	{
        		//cout<<i<<endl;
		if(co_old[i][0]>=0)
		{
			angle[i-1][0] = atan(co_old[i][1]/co_old[i][0]);
			//cout<<"tan "<<angle[i-1][0]<<" ";
		}
		else
		{
			angle[i-1][0] = 3.14159+atan(co_old[i][1]/co_old[i][0]);
			//cout<<"tan "<<angle[i-1][0]<<" ";
		}
		angle[i-1][1] = i;
		//cout<<endl;
	}

	//debug
	//cout<<"angle:"<<endl;

	//not debug
	sort(angle.begin(), angle.end());

    /*for (i = 1; i<=n_total+3; i = i+1)
	{
        cout<<i<<": "<<endl;
        cout<<"tan "<<angle[i-1][0]<<" ";
        cout<<"rank "<<angle[i-1][1]<<" ";
		cout<<endl;
	}*/


	
	de_rank = new double* [n_total];
	reg = new int [n_total];

	int t = 0;
	int rank=0;//, ranki, rankj = 0;

	j = 0;
	for(i = 0; i<n_total+3; i = i+1)
	{
		//cout<<"de_rank[%d]:"<<i<<endl;
		de_rank[j] = new double[d_total];
		
		
		//cout<<"angle "<<endl;
        //cout<<i<<" "<<angle[i][0]<<" ";
		rank = int(angle[i][1])-1;
		//cout<<"rank: "<<rank<<endl;
		
		cout<<"here"<<endl;
		if(rank!=-1)
        {
            for(t = 0; t<d_total; t = t+1)
    		{
    			de_rank[j][t] = demand[rank][t];
				reg[j] = arr_reg[rank];
    			cout<<de_rank[j][t]<<" ";
				cout<<reg[j];
    		}
			j = j+1;
			cout<<endl;
                    
		}
	}

	/*order coordinate*/

	new_d = new double* [n_total+4];
	co = new double* [n_total+4];


	co [0] = new double [5];
	co[0][0] = 0;
	co[0][1] = 0;
	

   int start, end; 
   int index;
   for(i= 0; i<n_total+3; i++)
   {
	   if(angle[i][1] ==0)
	   {
		   end = i;
	   }
   }

   start = end-2;
	
	
	for(i = 1; i<=n_total; i = i+1)
	{
		co [i] = new double [5];

		if(i>=start+1)
		{
			index = i+2;
		}
		else
		{
			index = i-1;
		}

		rank = int(angle[index][1]);
		
			co[i][0] = co_old[rank][0];
			co[i][1] = co_old[rank][1];
			//cout<<co[i][0]<<" "<<co[i][1]<<endl;
	
	}

	/*from coordinates*/
	
	double x;
	double y;


	for(i = 0; i<=n_total; i = i+1)
	{
		//file = fopen("test.txt", "a");
		new_d[i] = new double[n_total+4];

		for(j = 0; j<=n_total; j = j+1)
		{
			if(i==j)
			{
				//fprintf(file, "i=j=%d\n", i);
				new_d[i][j]=0;
				//cout<<new_d[i][j]<<" ";
			}
			else
			{
				x = co[i][0]-co[j][0];
				y = co[i][1]-co[j][1];

				//fprintf(file, "%f %f %f %f\n", co[i][0], co[i][1], co[j][0], co[j][1]);
				
				new_d[i][j] =sqrt(pow(x,2)+pow(y,2));
				//cout<<new_d[i][j]<<" ";
			}
		}

	}


	
	return angle;

}
