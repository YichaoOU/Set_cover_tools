//
//  model.h
//  
//
//  Created by Yating Liu on 3/9/15.

#pragma once


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <numeric>
#include <cassert>
#include <algorithm>
#include <metslib/mets.hh>
#include <string>
#include <cmath>
#include <unordered_map>

using namespace std;

extern int NUM, NUM_b, MOTIF;
extern double alpha;
extern double DELTA;
extern double penalty;
extern double beta;
//#define NUM 1041       //sequence number
//#define MOTIF 12   //motif number
typedef double gol_type;


class motif
{
public:
    int* mot;
    int seq;
    motif& add(const motif& other){
        for (int i = 0; i < seq; i++) {
            if (mot[i] == 0 && other.mot[i] ==0)
                mot[i] = 0;
            else
                mot[i] = 1;
        }
        return *this;
    }
    // num of coverage of this sequence
    double coverage() {
        double s = 0.0;
        for (int i = 0; i < seq; i++) {
            if (mot[i] == 1)
                s += 1.0;
        }
        return s;
    }
    
    motif(int NUM_seq) {
        seq = NUM_seq;
        mot = new int[seq];
        for (int i = 0; i < seq; i++) {
            mot[i] = 0;
        }
    }
    
    
    motif() {
        mot = new int[NUM];
        seq = NUM;
        for (int i = 0; i < NUM; i++) {
            mot[i] = 0;
        }
        
    }
    
    motif(vector<int>& m) {
        seq = m.size();
        mot = new int[seq];
        for (int i = 0; i < seq; i++) {
            mot[i] = m[i];
        }
    }
};



class my_sol : public mets::evaluable_solution,public motif
{
    
public:
    
    std::vector<motif> set;
    std::vector<motif> set_b;
    std::vector<bool> delta_m;
    double target_sum_m;
    double current_sum_m;

 my_sol(const std::vector<motif>& Set, const std::vector<motif>& Set_b, double sum)
    : delta_m(MOTIF, true),
    set(Set.begin(), Set.end()),
    set_b(Set_b.begin(), Set_b.end()),
    target_sum_m(sum),
    current_sum_m(MOTIF)
    { }
    
    bool accept (const feasible_solution& sol);
    double value(my_sol sol) const;
    double total_fore_coverage() const;
    double total_back_coverage() const;
    mets::gol_type cost_function() const;
    void copy_from(const mets::copyable& o);
    size_t size() const;
    size_t set_size() const;
    mets::gol_type what_if(int i, bool val) const;
    bool delta(int i) const;
    void delta(int i, bool val);
    motif element_fore(int i) const;
    
    friend ostream& operator<<(ostream& o, const motif& Motif);
    friend vector<pair<string,motif> > output_sol (const my_sol& s, vector<pair<string,motif> > mymotif);
    friend double cost(const my_sol& s);
   
};






