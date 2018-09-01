// File:RILP_simplex.cc
// Purpose: Implements a randomized algorithm based on the Integer Linear
// Programming formulation for
// the minimum discriminative set cover problem.
//
// Inputs: Three files (command line parameters) --- foreground.txt
//                                               --- background.txt
//                                               --- motif_mapping.txt
//         + Two real numbers:
//           KP == Need to cover 1-KP percent of the foreground
//           JP == Should cover at most JP percent of the background
//                 Note: randomized algorithm will always produce a solution
//                 that covers at least 1-KP percent of the foreground.   It may
//                 cover more than JP percent of the background.
//
// This version doesn't fail if it covers too much of the background.
// Two integers k and j
//          
//
// Added debugging information
// 
// 
#include <iostream>
#include <fstream>
#include <glpk.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <map>
#include <cassert>
#include "ReadMotifs.cc"

using namespace std;
string strip(string &x) {
  string temp;
  for (size_t i=0;i<x.size();i++) {
    if (!isspace(x[i])) {
      temp+=x[i];
    }
  }
  return temp;
}

// We are going to solve the set cover problem using 
// the MIP version of the integer linear programming
// formulation.
// Let n = number of motifs (the size of the adj_list.)
// Let m = number of sequences.
// Let m' = number of sequences in the first (foreground set = U)
// Let m-m' = number of sequences in the background set (second input file).
//
// The matrix A will have 2*m + 2 rows.
// Row 2i+1 and 2i+2 will consist of the constraints
// sum (i in Sj) u_j - v_i geq 0 and
// sum (i in Sj) uj  - K_i * vi
//
// The final two rows consist of 
// sum_{l=n+1}^{n+m'} v_i \geq m' - k
// and
// sum_{l=m'+1}^{m} v_i \leq j.
//
// Now, the matrix a is stored 
// AS A SPARSE MATRIX.
// The array ia gives the i index of each element.
// The array ja gives the j index of each element.
// The array ar give the value for each entry.
// In this case, the total number of elements is the sum of the sizes of the
// sets.
//
// Notice that we want gene i to be covered.  So, we want
// 1<= xi <= n
// Similarly, we want 0<= x_m+j <=1.
// 

vector<double> solve(vector<set<int> > &adj_list, int back_n, int fore_m, int MIP_k, int MIP_j) {
  vector<double> x;
  int n;
  int m;
  n = adj_list.size();
  m = back_n+fore_m;

  // DWJ added 07/31/17
  //cout << "***" << n+m << endl;

  vector<int> K;
  for (int i=0;i<m;i++) {
    int count = 0;
    for (int j=0;j<adj_list.size();j++) {
      if (adj_list[j].count(i) > 0) count++;
    }
    K.push_back(count);
  }

  // Calculate the size of the sparse matrix for A.
  // For each element i (0 .. m-1), there will be 
  // 2*K[i] + 2 elements.
  // Also, there will be an additional m entries in the 
  // matrix (for the two constraints.
  //
  int ASIZE=0;
  for (int i=0;i<m;i++) {
    ASIZE+=2*K[i]+2;
  }
  ASIZE+= m;

  
  glp_prob *lp;

  int *ia;
  int *ja;
  double *ar;
  double Z;



  ia = new int [ASIZE+2];
  ja = new int [ASIZE+2];
  ar = new double [ASIZE + 2];

  int index_count = 0;
  int row_count = 0;
  for (int i=0;i<m;i++) {
    // Enter the two rows for sequence i
    for (int j=0;j<adj_list.size();j++) {
      if ((adj_list[j].count(i)) > 0) {
	index_count++;
	ia[index_count] = row_count+1;
	ja[index_count] = j+1;
	ar[index_count] = 1.0;
      }
    }
    index_count++;
    ia[index_count] = row_count+1;
    ja[index_count] = n+i+1;
    ar[index_count] = -1.0;
    // Sum j=1 n (u_j containing element i) - v_i >=0.
    row_count++;
    for (int j=0;j<adj_list.size();j++) {
      if (adj_list[j].count(i) > 0) {
	index_count++;
	ia[index_count] = row_count+1;
	ja[index_count] = j+1;
	ar[index_count] = 1.0;
      }
    }
    index_count++;
    ia[index_count] = row_count+1;
    ja[index_count] = n+i+1;
    ar[index_count] = -K[i];
    row_count++;
  }
  for (int i=0;i<back_n;i++) {
    index_count++;
    ia[index_count] = row_count+1;
    ja[index_count] = n+i+1;
    ar[index_count] = 1;
  }
  row_count++;
  for (int i=back_n;i<m;i++) {
    index_count++;
    ia[index_count] = row_count+1;
    ja[index_count] = n+i+1;
    ar[index_count] = 1;
  }
  row_count++;
  // DWJ added 07/31/17
  //cout << "Number of rows = " << row_count << endl;

  //cout << "Matrix" << endl;
  //for (int i=1;i<=index_count; i++ ) {
  //  cout << ia[i] << " " << ja[i] << " " << ar[i] << endl;
  //}
  assert(index_count <ASIZE+2);
 
  lp = glp_create_prob();
  glp_set_prob_name(lp,"MINIMUM DISCRIMINATING SET COVER");
  glp_set_obj_dir(lp,GLP_MIN);

  //glp_set_int_parm(lp,LPX_K_MSGLEV,1);

  glp_add_rows(lp,row_count);
  int ri = 0;
  //cout << "Bounds" << endl;
  while (ri < row_count-2) {
    //cout << "Row (Lower) " << ri+1 << " " << 0 << endl;
    glp_set_row_bnds(lp,ri+1,GLP_LO,0,n);
    //cout << "Row (Upper) " << ri+2 << " " << 0 << endl;
    glp_set_row_bnds(lp,ri+2,GLP_UP,-n,0);
    ri++;
    ri++;
  }
  glp_set_row_bnds(lp,row_count-1,GLP_LO,back_n-MIP_k, back_n);
  glp_set_row_bnds(lp,row_count,GLP_UP,0,MIP_j);
  glp_add_cols(lp,n+m);
  for (int i=1;i<=n+m;i++) {
    glp_set_col_bnds(lp,i,GLP_DB,0.0,1.0);
    if (i<=n) {
      glp_set_obj_coef(lp,i,1.0);
    } else {
      glp_set_obj_coef(lp,i,0.0);
    }
    glp_set_col_kind(lp, i, GLP_BV);
  }
  //cout << n << " " << m << endl;
  
  glp_load_matrix(lp,index_count,ia,ja,ar);

  glp_smcp parm;
  glp_init_smcp(&parm);
  //DWJ changed 07/31/17
  // Set parm.msg_lev to GLP_MSG_ON to get debugging information.
  //parm.msg_lev = GLP_MSG_OFF;
  parm.msg_lev = GLP_MSG_OFF;
  //DWJ added 07/31/17
  parm.meth=GLP_DUALP;
  //glp_iocp parm;
  parm.tm_lim = 30*60*1000;
  

  //glp_init_iocp(&parm);
  //parm.presolve = GLP_ON;
  
  int err=glp_simplex(lp,&parm);
  if (err!=0) {
    //cout << "Error Flag = " << err << endl;
    //cout << "Likely infeasible" << endl;
  }
    //int err = glp_intopt(lp, &parm);
    //out << "Error Flag " << err << endl;
  double z = glp_get_obj_val(lp);
  //cout << "Optimum Relaxed ILP Solution = " << z << endl;
  x.resize(n);
  for (int i=1;i<=n;i++) {
    double x1 = glp_get_col_prim(lp, i);    
    //cout << x1 << endl;
    x[i-1] = x1;
  }
  
  glp_delete_prob(lp);
  delete [] ar;
  delete [] ja;
  delete [] ia;
  return x;
}

set<int> randomized(vector<double> &x, vector<set<int> > &adj_list, int back_n, int MIP_k,int MIP_j, int MIP_jp, bool &failed) {

  set<int> covered;
  set<int> cover;
  set<int> bad_covered;
  int tries = 0;
  while ((covered.size() < back_n - MIP_k) && (tries < 10)) {
    for (int i=0;i<x.size();i++) {
      double r = rand()/(RAND_MAX*1.0);
      //cout << r << " " << x[i] << endl;
      if ( r  < x[i]) {
	for (set<int>::iterator p = adj_list[i].begin();
	     p!=adj_list[i].end();++p) {
	  if ((*p)<back_n) {
	    covered.insert(*p);
	  } else {
	    bad_covered.insert(*p);
	  }
	}
	cover.insert(i);
      }
      if (covered.size()>= (back_n - MIP_k)) break; 
    }
    tries++;
  }
  //cout << "Tries = " << tries << endl;
  failed = false;
  if (tries == 10) { failed = true;}
  //cout << "Couldn't find one" << endl;}
  if (bad_covered.size() > MIP_jp) {
    // failed = true;  
  }
  return cover;
}



void get_data_first(istream &rin, 
		    map<string, int> &sequence_map,
		    vector<string> &sequences,
		    vector<string> &motifs,
		    map<string,int> &motif_map,
		    vector<set<int> > &adj_list) {

  string line;
  getline(rin,line);
  //cout << line << endl;
  istringstream in(line);
  string blank;
  getline(in,blank,',');
  //cout << "Should be blank = " << blank << endl;
  int count = 0;
  if (in.peek()==',') {
    char c;
    in.get(c);
  }
  while (!in.fail()) {
    string next;
    getline(in,next,',');
    next=strip(next);
    //cout << "*" << next << "*" <<endl;
    if (!in.fail()) {
      sequences.push_back(next);
      sequence_map[next]=count;
      //cout << "Seq " << next << "=" << count << endl;
      count++;
    }
  }
  //cout << "Here" << endl;
  count = 0;
  while (!rin.fail()) {
    string line;
    getline(rin,line);
    if (!rin.fail()) {
      istringstream in1(line);
      string m;
      getline(in1,m,',');
      m=strip(m);
      motifs.push_back(m);
      motif_map[m]=count;
      //cout << "Motif" << m << endl;
      count++;
      set<int> t;
      adj_list.push_back(t);
      int counter=0;
      while (!in1.fail()) {
	string next;
	getline(in1,next,',');
	next=strip(next);
	//cout << "*" << next << "*" <<endl;
	if (!in1.fail()) {
	  if (next=="1") {
	    adj_list[count-1].insert(counter);
	  }
	}
	counter++;
      }
    }
  }
}

// Read the second file

void get_data_second(istream &rin, 
		    map<string, int> &sequence_map,
		    vector<string> &sequences,
		    vector<string> &motifs,
		    map<string,int> &motif_map,
		    vector<set<int> > &adj_list) {

  string line;
  getline(rin,line);
  //cout << "Second line "<< endl;
  istringstream in(line);
  string blank;
  getline(in,blank,',');
  //cout << "Should be blank = " << blank << endl;
  int count = sequences.size();
  int n = count;
  if (in.peek()==',') {
    char c;
    in.get(c);
  }
  while (!in.fail()) {
    string next;
    getline(in,next,',');
    next=strip(next);
    //cout << "*" << next << "*" <<endl;
    if (!in.fail()) {
      if (sequence_map.count(next)==0) {
	sequences.push_back(next);
	sequence_map[next]=count;
	//cout << "Seq " << next << "=" << count << endl;
	count++;
      }
    }
  }
  //cout << "Here" << endl;
  count = 0;
  while (!rin.fail()) {
    string line;
    getline(rin,line);
    if (!rin.fail()) {
      istringstream in1(line);
      string m;
      getline(in1,m,',');
      m=strip(m);
      if (motif_map.count(m)>0) {
	count = motif_map[m];
      } else {
	cout << "Error???" << endl;
	motifs.push_back(m);
	motif_map[m]=count;
	cout << "Motif" << m << endl;
	count++;
      }
      set<int> t;
      
      //
      int counter=n;
      while (!in1.fail()) {
	string next;
	getline(in1,next,',');
	next=strip(next);
	//cout << "*" << next << "*" <<endl;
	if (!in1.fail()) {
	  if (next=="1") {
	    adj_list[count].insert(counter);
	  }
	}
	counter++;
      }
    }
  }
}



int main(int argc, char *argv[]) {

  ReadMotifs obj;
  if (argc!=6) {
    cout << "Invalid # of args" << argc << endl;
    exit(-1);
  }

  string f1 = argv[1];
  string f2 = argv[2];
  string f3 = argv[3];

  string KPS = argv[4];
  string JPS = argv[5];

  istringstream ikP(KPS);
  istringstream ijP(JPS);

  double KP;
  double JP;
  ikP >> KP;
  ijP >> JP;
  //ReadMotifs obj;

  obj.ReadForegroundSet(f1);
  obj.ReadBackgroundSet(f2,obj.n_f);
  obj.ReadMotifData(f3);


  

  
  //map<string, int> sequence_map;
  //vector<string> sequences;
  //vector<string> motifs;
  //map<string,int> motif_map;
  //vector<set<int> > adj_list;

  int MIP_k, MIP_j;

  MIP_k = (int) (obj.n_f*KP);
  MIP_j = (int) (obj.n_b*JP);
  
  //cin >> MIP_k;
  //cin >> MIP_j;
  
  //ifstream fin,fin1;
  //fin.open(argv[1]);

  
  //get_data_first(fin,sequence_map, sequences,motifs, motif_map,adj_list);

  //int back_n = sequences.size();
  // All numbers < back_n are in the universe.  
  // All numbers >= are in the foreground.

  //fin1.open(argv[2]);

  //get_data_second(fin1,sequence_map, sequences,motifs, motif_map,adj_list);
  //int for_n = sequences.size();
  //cout << "foreground sequence #'s = " << 0 << " to " << back_n - 1 << endl;
  //cout << "background sequence #'s = " << back_n << "to " << for_n-1 << endl;
  //cout << "number of motifs = " << motifs.size() << endl;


  vector<double> solution;
  set<int> sol_set;
  set<int> best_set;

  solution = solve(obj.adj_list, obj.n_f, obj.n_b, MIP_k, MIP_j);
  
  bool failed;
  bool first = true;
  //sol_set = randomized(solution,adj_list, back_n, MIP_k,MIP_j,failed);
  
  //best_set = sol_set;
 
  for (int i= 0;i<500;i++) {
    sol_set = randomized(solution,obj.adj_list, obj.n_f, MIP_k,MIP_j,MIP_j*2,failed);
    if (first) {
      if (!failed) {
	best_set=sol_set;
	first = false;
      }
    } else {
      if (!failed) {
	if (sol_set.size() < best_set.size()) {
	  best_set = sol_set;
	}
      }
    }
  }
  if (!first) {
    //cout << "Found a solution with " << best_set.size() << " elements " << endl;
    set<int> covered_foreground;
    set<int> covered_background;
    for (auto p=best_set.begin();p!=best_set.end();++p) {
      //cout << obj.motif_names[*p] << endl;
      for (auto p1=obj.adj_list[*p].begin();p1!=obj.adj_list[*p].end();++p1) {
	if (*p1>=obj.n_f) {
	  covered_background.insert(*p1);
	} else {
	  covered_foreground.insert(*p1);
	}
      }
    }
    //    cout << "Covers " << covered_foreground.size() << " elements of foreground " << endl;
    cout << covered_foreground.size() << "\t" << covered_background.size() << endl;
    
    //cout << "Covers " << covered_background.size() << " elements of background " << endl;
    for (auto p=best_set.begin();p!=best_set.end();++p) {
      cout << obj.motif_names[*p] << endl;
    }
    cout << "EOF" << endl;
    
  } else {
    cout << 0 << "\t" << 0 << endl;
    cout << "Failed to find a solution " << endl;
    cout << "EOF" << endl;
  }
}
