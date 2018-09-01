//
//  move.h
//  
//
//  Created by Yating Liu on 3/10/15.
//
//

#pragma once

#include <typeinfo>
#include <metslib/mets.hh>
#include "model.h"

/// @brief A move that toogles a boolean value.
///
/// A move is a way to modify one solution into a neighboring one. The
/// actual neighborhood explored depends on the set (or subset) of
/// moves that is used in the move_manager.
///

class toggle : public mets::mana_move {
    int index_m;
public:
    toggle(int i) : index_m(i) {}
    /// @brief Evaluate the cost after the move without actually
    /// performing it.
    ///
    mets::gol_type evaluate(const mets::feasible_solution& cs) const
    {
        const my_sol& model = dynamic_cast<const my_sol&>(cs);
        return model.what_if(index_m, !model.delta(index_m));
    }
    /// @brief Apply this move
    void apply(mets::feasible_solution& s) const
    {
        my_sol& model = dynamic_cast<my_sol&>(s);
        model.delta(index_m, !model.delta(index_m));
    }
    /// @brief Virtual method used by the tabu list to keep a copy of a move
    mana_move* clone() const { return new toggle(index_m);}
    
    /// @brief Virtual method used by the tabu list to keep a copy of
    /// the opposite of a move
    // mana_move* opposite_of() const { return new toggle(index_m); }
    
    /// @brief A number identifying this move as much as possible. It's
    /// used for quick testing presence of moves in the tabu list via an
    /// hash map.
    size_t hash() const { return index_m; }
    
    /// @brief Comparison of moves: used to test if a move in in a tabu list
    bool operator==(const mets::mana_move& o) const
    {
        //  We first check if the move if of our same type (different move
        // types can coexist in the same tabu list).
        try {
            const toggle& other = dynamic_cast<const toggle&>(o);
            // Then we check for equality
            return (this->index_m == other.index_m);
        } catch (std::bad_cast& e) {
            std::cerr << "bad cast?" << std::endl;
            return false;
        }
    }
    
};
