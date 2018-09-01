//
//  neighborhoods.h
//  
//
//  Created by Yating Liu on 3/10/15.
//
//
#pragma once

#include <iostream>
#include "moves.h"

class full_neighborhood
{
public:
    std::vector<toggle*> moves_m;
    typedef std::vector<toggle*>::iterator iterator;
    iterator begin() { return moves_m.begin(); }
    iterator end() { return moves_m.end(); }
    
    full_neighborhood(int problem_size)
    {
        for(unsigned int ii = 0; ii != problem_size; ++ii)
            moves_m.push_back(new toggle(ii));
    }
    
    ~full_neighborhood()
    {
        // delete all moves
        for(iterator ii = begin(); ii != end(); ++ii)
            delete (*ii);
    }
    
    void refresh(mets::feasible_solution& s)
    {
        // This method can be used to adapt the neighborhood to the
        // current problem instance.  In our simple case there is no need
        // to update the neighborhood since its moves does not depend on
        // the current solution considered.
    }
    
};
