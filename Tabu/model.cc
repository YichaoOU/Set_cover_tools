//
//  model.cpp
//  
//
//  Created by Yating Liu on 4/11/16.
//
//

#include <stdio.h>
#include "model.h"


bool my_sol::accept (const feasible_solution& sol)
{
    if (total_fore_coverage() < 1.0)
        return false;
    else
        return true;
}

double my_sol::value(my_sol sol) const 
//calculate optimize function value, goal: minimizing the value
{
    
    return sol.set_size() + beta * (alpha * (1 - sol.total_fore_coverage()) + (1 - alpha) * sol.total_back_coverage());
}



double my_sol::total_fore_coverage() const
{
    motif* sum_fore = new motif(NUM);
    
    for(int j = 0; j < MOTIF; ++j) {
        if (delta_m[j]) {
            sum_fore->add(this->set[j]);
            
        }
    }
    double coverage_fore = sum_fore->coverage()/NUM;
    delete sum_fore;
    return coverage_fore;
}

double my_sol::total_back_coverage() const
{
    motif* sum_back = new motif(NUM_b);
    
    for(int j = 0; j < MOTIF; ++j) {
        if (delta_m[j]) {
            sum_back->add(this->set_b[j]);
            
        }
    }
    double coverage_back = sum_back->coverage()/NUM_b;
    delete sum_back;
    return coverage_back;
}


gol_type my_sol::cost_function() const
{
    double cost = current_sum_m - target_sum_m;
    return cost;
}



/// @brief This method is needed by the algorithm to record the best solution.
void my_sol::copy_from(const mets::copyable& o)
{
    const my_sol& s = dynamic_cast<const my_sol&>(o);
    delta_m = s.delta_m;
    set = s.set;
    set_b = s.set_b;
    target_sum_m = s.target_sum_m;
    current_sum_m = s.current_sum_m;
}

size_t my_sol::size() const
{ return delta_m.size(); }

// number of motifs in the set
size_t my_sol::set_size() const     //set size
{
    size_t s = 0;
    for (int i = 0; i < MOTIF; i++) {
        if (delta_m[i] == 1)
            s += 1;
    }
    return s;
}
/// @brief Evaluates the cost of a change without actually doing it.

mets::gol_type my_sol::what_if(int i, bool val) const
{
    double newcost = 0.0;
    double oricost = 0.0;
    double oricover = total_fore_coverage();
    double newcover = 0.0;
    double diff = 0;
    my_sol testsol(set, set_b, target_sum_m);
    testsol.delta_m = delta_m;
    oricost = value(testsol);
    if (delta_m[i] && !val) {
        testsol.delta_m[i] = val;
        newcost = value(testsol);
        
        diff = newcost - target_sum_m;
        newcover = testsol.total_fore_coverage();
    }
    
    else if (!delta_m[i] && val) {
        testsol.delta_m[i] = val;
        newcost = value(testsol);
        diff = newcost - target_sum_m;
        newcover = testsol.total_fore_coverage();
    }
    
    if (newcover - oricover > 0 && newcover - oricover < DELTA && val)
        return penalty * diff;
    else
        return diff;
}

bool my_sol::delta(int i) const { return delta_m[i];}

void my_sol::delta(int i, bool val)
{
    if (delta_m[i] && !val) {
        delta_m[i] = val;
        current_sum_m = set_size() + beta * (alpha * (1 - total_fore_coverage()) + (1 - alpha) * total_back_coverage());
    }
    else if (!delta_m[i] && val) {
        delta_m[i] = val;
        current_sum_m = set_size() + beta * (alpha * (1 - total_fore_coverage()) + (1 - alpha) * total_back_coverage());
    }
    
}

motif my_sol::element_fore(int i) const { return set[i]; }


ostream& operator<<(ostream& o, const motif& Motif) {
    for (int i = 0; i < Motif.seq; i++) {
        o<<Motif.mot[i]<<' ';
    }
    return o;
    
}

vector<pair<string,motif> > output_sol (const my_sol& s, vector<pair<string,motif> > mymotif)
{
    vector<pair<string,motif> >::iterator it = mymotif.begin();
    vector<pair<string,motif> > output;
    for (int ii = 0; ii < MOTIF; ++ii) {
        if (s.delta(ii)) {
            output.push_back(make_pair(it->first,it->second));
            //o << s.element_fore(ii) << endl;
        }
        it++;
    }
    
    return output;
    
}

double cost(const my_sol& s) {
    return s.value(s);
}

