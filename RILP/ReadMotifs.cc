// File: ReadMotifs.cc
// Purpose: Class ReadMotifs is a wrapper for 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cstdlib>
#include <sstream>
using namespace std;
class ReadMotifs {
public:
  vector<set<int> > adj_list;
  map<string,int> sequences;
  
  
  vector<string> sequence_names;
  //  vector<string> background_names;


  map<string,int> motifs;
  vector<string> motif_names;
  
  
  int n; // Total sequences;
  int n_f; // Total foreground sequences
  int n_b; // Total background sequences
  int m; // Total number of motifs;
  
  string strip_whitespace(string &in) {
    string temp;
    for (int i=0;i<in.length();i++) {
      if (!isspace(in[i])) {
	  temp+=in[i];
	}
    }
    return temp;
  }
  
  void ReadForegroundSet(string &filename1) {
    fstream fin1;
    fin1.open(filename1);
    if (fin1.fail()) {
      cout << "Couldn't open " << filename1 <<endl;
      exit(-1);
    }
    string line;
    string stripped;
    int count=0;
    while (!fin1.eof()) {
      getline(fin1,line);
      stripped = strip_whitespace(line);
      if (!fin1.fail()) {
	if (sequences.count(stripped)>0) {
	  cout << "Saw foreground " << stripped << " before -- ignoring " << endl;
	} else {
	  sequences[stripped]=count;
	  sequence_names.push_back(stripped);
	  count++;
	}
      }
    }
    n_f = count;
    


  }
  void ReadBackgroundSet(string &filename2, int n_f) {
    fstream fin1;
    fin1.open(filename2);
    if (fin1.fail()) {
      cout << "Couldn't open " << filename2 <<endl;
      exit(-1);
    }
    string line;
    string stripped;
    int count=0;
    while (!fin1.eof()) {
      getline(fin1,line);
      stripped = strip_whitespace(line);
      if (!fin1.fail()) {
	if (sequences.count(stripped)>0) {
	  cout << "Saw background " << stripped << " before -- ignoring " << endl;
	} else {
	  sequences[stripped]=n_f+count;
	  sequence_names.push_back(stripped);
	  count++;
	}
      }
    }
    n_b = count;
    
    n = n_f+n_b;
    
  }
  
  void ReadMotifData(string &filename3) {
    fstream fin1;
    fin1.open(filename3);
    if (fin1.fail()) {
      cout << "Couldn't open " << filename3 <<endl;
      exit(-1);
    }

    string line;
    getline(fin1,line);
    //cout << "Header line --- =" << endl;
    //cout << line;
    //getline(fin1,line);
    //cout << line << endl;
    int motif_count =0;
    while (!fin1.eof() ) {
      getline(fin1,line);
      if (!fin1.fail()) {
	
	istringstream in(line);
	string motif;
	string sequence;
	in >> motif;
	in >> sequence;
	if (sequences.count(sequence) == 0) {
	  cout << "Fatal Error -- sequence not found " << sequence << endl;
	  exit(-1);
	} else {
	  if (motifs.count(motif) == 0) {
	    motifs[motif]=motif_count;
	    motif_names.push_back(motif);
	    motif_count++;
	  }
	  adj_list.resize(motif_names.size());
	  adj_list[motifs[motif]].insert(sequences[sequence]);
	}
      }
    }
    m=motif_count;
  }
};
