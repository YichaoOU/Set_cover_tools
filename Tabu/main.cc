//
//  main3.cpp
//  
//
//  Created by Yating Liu on 4/1/16.
//
//

// this version is used to get run data with different parameters automatically untill no best result
#include <iostream>
#include <stdlib.h>  //srand, rand
#include <time.h>    // time
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <map>
#include <unordered_map>
#include "model.h"
#include "moves.h"
#include "neighborhoods.h"
using namespace std;

double TARGET = 1;
int iterations = 1; 
int maxIteration = 1000; // Max non-improving iterations
int NUM, MOTIF, NUM_b;
double DELTA = 0.05;
double alpha = 0.5;
double beta = 0;
double penalty = 10000;
double tenure = 0.2;
//int num_f, num_b; // used to store num of seq in forground and background


struct logger : public mets::search_listener<full_neighborhood>
{
    explicit
    logger(std::ofstream& o)
    : mets::search_listener<full_neighborhood> (), iteration(0), os(o)
    {
    }
    
    
    
    void
    update(mets::abstract_search<full_neighborhood>  * as)
    {
        const my_sol& ss = dynamic_cast<const my_sol&>(as->working());
        if(as->step() == mets::abstract_search<full_neighborhood>::MOVE_MADE)
        {
            os  << iteration++ << ": " << "total foreground coverage  " << ss.total_fore_coverage() << endl << "set size  " << ss.set_size()<<endl<< "cost function  " << ss.cost_function() << endl;
            /*
             for (it = mymotif.begin(); it != mymotif.end(); ++it)
             {
             if (it->second == ss)
             {
             key = it->first;
             break;
             }
             }
             */
            //os << ss << "\n";
        }
    }
    
protected:
    int iteration;
    ofstream& os;
};

void randomizeData(vector<motif> ori_motif, vector<string> ori_name, vector<pair<string,motif> >& mymotif, vector<pair<string,motif> > tmp_motif, vector<pair<string,motif> >& mymotif_b, vector<motif>& set_m, vector<motif>& set_b) {
    /*randomize the order of mymotif*/
    mymotif.clear();
    mymotif_b.clear();
    set_m.clear();
    set_b.clear();
    int random = -1;
    while (!ori_name.empty()) {
        srand(time(NULL));
        random = rand() % ori_name.size();
        // cout << "random num = " << random << " "<< endl;
        //cout << ori_name[random] << "     ";
        set_m.push_back(ori_motif[random]);
        mymotif.push_back(make_pair(ori_name[random],ori_motif[random]));
        ori_motif.erase(ori_motif.begin()+random);
        ori_name.erase(ori_name.begin()+random);
    }
    
    /* order background motif in the same order as foreground */
    bool flag = 0;  //initialized to 0
    for(vector<pair<string, motif> >::iterator it = mymotif.begin(); it != mymotif.end(); it++)
    {
        flag = 0;  //initialized to 0
        for(vector<pair<string, motif> >::iterator it_b = tmp_motif.begin();it_b != tmp_motif.end();it_b++){
            if (it_b->first == it->first) {
                mymotif_b.push_back(make_pair(it_b->first,it_b->second));
                flag = 1;  //set to 1 if found the same motif
                set_b.push_back(it_b->second);
                break;
            }
        }
        if (!flag) {
            cout << "Error!! Not found that motif in background file\n\n";
            exit(1);
        }
        
    }
}

void printMotifs(ofstream& outfile, vector<pair<string,motif> > mymotif, vector<pair<string,motif> > mymotif_b) {
    vector<pair<string, motif> >::iterator iter_b = mymotif_b.begin();
    outfile << "Motif_Name, Foreground Coverage, Foreground Seq, Background Coverage, Background Seq, Sensitivity, Speficity, Accuracy" << endl;
    for(vector<pair<string, motif> >::iterator iter = mymotif.begin(); iter != mymotif.end(); iter++)
    {
        outfile << iter->first << "," << (iter->second).coverage() << "," << NUM << ","; //fore coverage
        outfile << (iter_b->second).coverage() << "," << NUM_b << "," //back coverage
                << (iter->second).coverage() / NUM << "," // sensitivity
                << 1 - (iter_b->second).coverage() / NUM_b << "," //speficity
                << ((iter->second).coverage() + NUM_b - (iter_b->second).coverage()) / (NUM + NUM_b) //acc
               << endl;
        iter_b++;
        
    }
    

    
}

void tabuSearch(vector<motif> set_m, vector<motif> set_b, double& best_solution, vector<pair<string,motif> >& selected_motifs, vector<pair<string,motif> >& selected_motifs_b,vector<pair<string,motif> > mymotif, vector<pair<string,motif> > mymotif_b, ofstream& out)
{
    
    my_sol model(set_m,set_b, TARGET);
    my_sol best(model);
    full_neighborhood neigh(model.size());
    logger g(out);
    mets::simple_tabu_list tabu_list(tenure * MOTIF);
    mets::best_ever_criteria aspiration_criteria;
    mets::noimprove_termination_criteria noimprove(maxIteration);
    mets::threshold_termination_criteria threshold_noimprove(&noimprove, 0);
    mets::best_ever_solution best_recorder(best);
    
    mets::tabu_search<full_neighborhood> algorithm(model,
                                                   best_recorder,
                                                   neigh,
                                                   tabu_list,
                                                   aspiration_criteria,
                                                   threshold_noimprove);
    
    algorithm.attach(g);
    algorithm.search();
    cout << "current cost: " << best_recorder.best_seen().cost_function() << endl;
    if (best_solution > cost((const my_sol&)best_recorder.best_seen())) {
       // cout << "previous cost: " << best_solution << endl;
        best_solution = cost((const my_sol&)best_recorder.best_seen());
        vector<pair<string,motif> > tmp_f = output_sol((const my_sol&)best_recorder.best_seen(), mymotif);
        selected_motifs.assign(tmp_f.begin(),tmp_f.end());
        vector<pair<string,motif> > tmp_b = output_sol((const my_sol&)best_recorder.best_seen(), mymotif_b);
        selected_motifs_b.assign(tmp_b.begin(), tmp_b.end());
       // cout << "changed cost: " << best_solution << endl;
    }
    
    
}

void argumentsParsing(int argc, char* argv[], string &input, string &output) {
    if (argc == 1) {
        cout <<"Tabu-PNPSC command line: "<< endl
        <<"-in inputFile" <<endl
        <<"-out outputFile" <<endl
        <<"-a alpha"<<endl
        <<"-b beta" <<endl
        <<"-t tenure" <<endl
        <<"-max maxIteration" <<endl
        <<"-d delta" <<endl
        <<"-p penalty"<<endl
        <<"-iter iterations" << endl;
    }
    string argument; 
    for (int count = 1; count < argc; count++) {
        cout << argv[count] << endl;
        argument = argv[count];
        if (argument == "-in") {
            input = argv[++count];
            cout << input << endl;
        }
        else if (argument  == "-out") {
            output = argv[++count];
            
        }
        else if (argument == "-a") {
            alpha = atof(argv[++count]);
            
        }
        else if (argument == "-b") {
            beta = atof(argv[++count]);
            
        }
        else if (argument == "-t") {
            tenure = atof(argv[++count]);
            
        }
        else if (argument == "-max") {
            maxIteration = atoi(argv[++count]);
            
        }
        else if (argument == "-d") {
            DELTA = atof(argv[++count]);            
        }
        else if (argument == "-p") 
            penalty = atof(argv[++count]);
        else if (argument == "-iter") 
            iterations = atoi(argv[++count]);
        else {
            cout << "invalid arguments! \n";
            exit(1);
        }               
    }
   
    if (input == "" || output == "") {
        cout<<"no input or output file"<<endl;
        exit(1);
    }
    
}

void printArguments(string input, string output) {
     cout << "Input file = "<< input << endl
    << "Output file = " << output << endl
    << "alpha = " << alpha << endl
    << "beta = " << beta << endl
    << "tenure = " << tenure << endl
    << "max non-improving iteration = " << maxIteration << endl
    << "delta = " << DELTA << endl
    << "penalty = " << penalty << endl
    << "iterations (running the tabu-PNPSC program $iterations$ times) = " << iterations << endl;   
    
}

vector<string> split(string line, string delimiter) {
    vector<string> list;
    string s = line;
    size_t pos = 0;
    string token;
    while((pos = s.find(delimiter)) != string::npos) {
        //cout << "pos = " << pos << "\n";
        token = s.substr(0, pos);
        //cout << token << "--";
        token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
        if (token != "") 
            
            list.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    list.push_back(s);
    return list;
}

int main (int argc, char* argv[]) {
    
    vector<pair<string,motif> > mymotif;
    vector<pair<string,motif> > tmp_motif;
    vector<pair<string,motif> > mymotif_b;
    vector<pair<string,motif> >::iterator iter = mymotif.begin();
    vector<pair<string,motif> >::iterator iter_b = mymotif_b.begin();
    double best_solution;
    vector<pair<string,motif> > selected_motifs;
    vector<pair<string,motif> > selected_motifs_b;
    string input = "", output = "";
    argumentsParsing(argc, argv, input, output);
    
    ifstream in;   //input matrix
    ofstream out;
   // ofstream outfile("motifcoverage.csv", ios::app);
    ofstream log("log.txt");
    in.open(input);
    out.open(output);
    std::string line;
    std::string line2;
    vector<motif> set_m;
    vector<motif> set_b;
    vector<motif> ori_motif;
    vector<string> ori_name;
    vector<string> seq_name;
    if (in.is_open()) {
        string mot_name;
        vector<string> seq;
        string seq_type;
        string s;
        getline(in, line);
        ori_name = split(line, ",");
        ori_name.pop_back();
        MOTIF = ori_name.size();
        beta = MOTIF;
        //cout << MOTIF << "\n";
        vector< vector<int>> matrix_f(MOTIF);
        vector< vector<int>> matrix_b(MOTIF);
        while (in >> s) {
            seq = split(s, ",");
            seq_name.push_back(seq.front());
            seq_type = seq.back();
            //cout << s << "\n" << seq_type << "\n" << seq.size() << "\n";
            if (stod(seq_type) > 0) {
                for (int j = 1; j <= MOTIF; j++) {
                   // cout << stoi(seq[j]) << ",";
                    matrix_f[j-1].push_back(stoi(seq[j]));
                }
            }
            else if (stod(seq_type) < 0) {
                for (int j = 1; j <= MOTIF; j++) {
                    matrix_b[j-1].push_back(stoi(seq[j]));
                }
            }
            else {
                cout << "Error: Wrong type! " << seq_type << "\n";
            }
        }
        NUM = matrix_f[0].size();
        NUM_b = matrix_b[0].size();
            
        for (int i = 0; i < MOTIF; i++) {
            motif m_fore(matrix_f[i]);
            motif m_back(matrix_b[i]);
            ori_motif.push_back(m_fore);
            mot_name = ori_name[i];
            tmp_motif.push_back(make_pair(mot_name, m_back));
        }
        in.close();
    }
     else
        cout<<"Foreground Input file cannot open\n";
    
        
    
    printArguments(input, output);
    /* output file */
    cout<<"NUM of motif: "<< MOTIF << "Num of foreground seq :"<< NUM << "Num of background seq : "
    << NUM_b << endl;
    
    cout <<"=====================================================\n";
    best_solution = MOTIF + NUM + NUM_b;  //initialize cost to max
    //iterations = MOTIF; //set the iteration to num of motifs
    for (int k = 0; k < iterations; k++) {
       // cout << "iter " << k << endl;
        /*randomize the order of motifs */
        randomizeData(ori_motif, ori_name, mymotif, tmp_motif, mymotif_b, set_m, set_b);
        /* print each motif info to outfile */
       // printMotifs(outfile, mymotif, mymotif_b);
        /*tabu search*/
        tabuSearch(set_m, set_b, best_solution, selected_motifs, selected_motifs_b ,mymotif, mymotif_b,log);
        cout << "Objective function: " << best_solution << endl;
    }
    

        string tmp;
        int j = 0;
        iter = mymotif.begin();
        iter_b = mymotif_b.begin();
        motif selected_b(NUM_b);  //used to store background motifs
        motif selected_f(NUM); //used to store foreground motifs
    vector<pair<string,motif> >::iterator it = selected_motifs.begin();
    vector<pair<string,motif> >::iterator it_b = selected_motifs_b.begin();
    for (int i = 0; i < selected_motifs.size(); i++) {
        out << it->first << endl;
        cout << it->first << endl << "Foreground Covered sequences: " << (it->second).coverage() << endl
        << "Background covered sequences: " << (it_b->second).coverage() << endl << "==========" << endl;
        selected_b.add(it_b->second);
            selected_f.add(it->second);
        it++;
        it_b++;
        }
        cout << "Best solution: " << best_solution << endl;
        cout << "Num of features: " << selected_motifs.size() << endl;
        cout << "Total foreground coverage: " << selected_f.coverage() / NUM << endl
            << "Total background coverage: " << selected_b.coverage() / NUM_b << endl;
 
    return 0;
}
