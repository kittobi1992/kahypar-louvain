/***************************************************************************
 *  Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_LIB_DATASTRUCTURE_NEIGHBORHOODHYPERGRAPH_H_
#define SRC_LIB_DATASTRUCTURE_NEIGHBORHOODHYPERGRAPH_H_

#include <vector>
#include <algorithm>
#include <sstream>

#include "lib/definitions.h"
#include "lib/datastructure/SparseSet.h"

using defs::Hypergraph;
using defs::HypernodeID;

namespace datastructure {

class NeighborhoodHypergraph {
 public:
  explicit NeighborhoodHypergraph(Hypergraph& hypergraph) :
      _hg(hypergraph),
      _neighbors(hypergraph.initialNumNodes()),
      _is_active(hypergraph.initialNumNodes(),true) {
    SparseSet<int> neighbors(hypergraph.initialNumNodes());

    for (const auto hn : hypergraph.nodes()) {
      for (const auto he : hypergraph.incidentEdges(hn)) {
        for (const auto pin : hypergraph.pins(he)) {
          neighbors.add(pin);
        }
      }

      _neighbors[hn].swap(neighbors.elements());
      std::sort(_neighbors[hn].begin(), _neighbors[hn].end());

      neighbors.clear();
    }
  }
  
  void contract(const HypernodeID x, const HypernodeID y) {
    ASSERT(x < _neighbors.size(), "HN " << x << " isn't a valid hypernode!");
    ASSERT(y < _neighbors.size(), "HN " << y << " isn't a valid hypernode!");
    int i = 0, j = 0;
    int size_x = _neighbors[x].size(), size_y = _neighbors[y].size();
    std::vector<int> new_x;
    //Merging the two neighborhoods vector of hypernode x and y
    while(i < size_x || j < size_y) {
      //If we reach the end of the list of hypernode x we can insert the remaining elements
      //from hypernode y at end of new_x
      if(i == size_x) {
	new_x.insert(new_x.end(),_neighbors[y].begin()+j,_neighbors[y].end());
	j = size_y;
	continue;
      } 
      //If we reach the end of the list of hypernode y we can insert the remaining elements
      //from hypernode x at end of new_x
      else if(j == size_y) {
	new_x.insert(new_x.end(),_neighbors[x].begin()+i,_neighbors[x].end());
	i = size_x;
	continue;
      }
      
      //Ignoring the contraction partner of y in both list
      if(_neighbors[x][i] == y) {
	i++; continue;
      } else if(_neighbors[y][j] == y) {
	j++; continue;
      }
      
      //Merging the two lists
      //Note: Special case if both elements are equal
      if(_neighbors[x][i] < _neighbors[y][j]) {
	  new_x.push_back(_neighbors[x][i++]);
      } else if(_neighbors[x][i] > _neighbors[y][j]) {
	  new_x.push_back(_neighbors[y][j++]);
      } 
      //If both elements are equal, than we have to mark this in the list of hypernode y
      //to restore the correct list if we uncontract both hypernodes.
      else {
	new_x.push_back(_neighbors[x][i++]);
	_neighbors[y][j] = -_neighbors[y][j++];
      }
    }
    
    //Validate if the values in the inactive hypernode y are correct
    ASSERT([&]() {
	for(int k = 0; k < _neighbors[y].size(); ++k) {
	  if(_neighbors[y][k] == y)
	    continue;
	  int search_pin = std::abs(_neighbors[y][k]);
	  auto find = std::lower_bound(_neighbors[x].begin(),_neighbors[x].end(),search_pin);
	  if(find != _neighbors[x].end()) {
	    //If search_pin is contained in both neighborhood lists, the value in _neighbors[y][k]
	    //should be negative.
	    if(*find == search_pin && _neighbors[y][k] > 0) {
		LOG("Neighborhood list of HN " << x << " and HN " << y << " contains both HN " << search_pin
		    << "! Therefore value in HN " << y << " neighborhood list should be negative.");
		return false;
	    } 
	    //If search_pin is only contained in the neighborhood list of y, the value in _neighbors[y][k]
	    //should be positive.
	    else if(*find != search_pin && _neighbors[y][k] < 0) {
	      LOG("Neighborhood list of HN " << x << " didn't contain HN " << search_pin
		  << "! Therefore value in HN " << y << " neighborhood list should be positive.");
	      return false;	      
	    }
	  } 
	  //If search_pin is only contained in the neighborhood list of y, the value in _neighbors[y][k]
	  //should be positive.
	  else if(_neighbors[y][k] < 0) {
	      LOG("Neighborhood list of HN " << x << " didn't contain HN " << search_pin
		  << "! Therefore value in HN " << y << " neighborhood list should be positive.");
	      return false;
	  }
	}
        return true;
      } (), "Wrong value in inactive neighborhood list of HN "<<y<<"!");
    
    //Remove hypernode y in all neighborhood lists of y's neighbors
    for(int k = 0; k < _neighbors[y].size(); ++k) {
      HypernodeID pin = static_cast<HypernodeID>(std::abs(_neighbors[y][k]));
      if(_is_active[pin] && pin != y) {
	//TODO(heuer): Naive removing variant. Another approach is to use lists instead of vectors,
	//but this makes compression very difficult.
	auto pos = std::lower_bound(_neighbors[pin].begin(),_neighbors[pin].end(),y);
	_neighbors[pin].erase(pos);
      }
    }
    
    //Apply new neighborhood list of hypernode x
    swap(_neighbors[x],new_x);
    //Disable hypernode y
    _is_active[y] = false;
    
    //Validates if hypernode y is removed from neighborhood list of neighbors of y.
    ASSERT([&]() {
	for(int k = 0; k < _neighbors[y].size(); ++k) {
	  HypernodeID pin = static_cast<HypernodeID>(std::abs(_neighbors[y][k]));
	  if(_is_active[pin]) {
	    auto find = std::lower_bound(_neighbors[pin].begin(),_neighbors[pin].end(),y);
	    if(find != _neighbors[pin].end() && *find == y) {
		LOG("Neighborhood list of HN " << pin << " contains HN " << y << "!");
		printNeighborhood(pin);
		return false;
	    }
	  }
	}
        return true;
      } (), "HN " << y << " should be removed from all neighbors!");
    
    //Validates if the neighborhood list of hypernode x contains the correct hypernodes
    ASSERT([&]() {
	std::vector<bool> match(_neighbors[x].size(),false);
	for(HyperedgeID he : _hg.incidentEdges(x)) {
	  for(HypernodeID pin : _hg.pins(he)) {
	    bool found = false;
	    for(int k = 0; k < _neighbors[x].size(); ++k) {
	      if(_neighbors[x][k] == pin) {
		match[k] = true; found = true;
		break;
	      }
	    }
	    //Hypernode pin should be contained in the neighborhood list of hypernode x.
	    if(!found) {
	      LOG("HN " << pin << " should be in the neighborhood list of HN " << x << "!");
	      return false;
	    }
	  }
	}
	for(int k = 0; k < _neighbors[x].size(); ++k) {
	  //Hypernode _neighbors[x][k] isn't in the neighborhood of hypernode x.
	  if(!match[k]) {
	    LOG("HN " << _neighbors[x][k] << " isn't in the neighborhood of HN " << x << "!");
	    return false;
	  }
	}
        return true;
      } (), "Creating new neighborhood list of HN " << x << " failed!");
    
  }
  
  
  std::string getNeighborhoodHypergraphStats() {
    std::ostringstream oss;
    oss << "\nNeighborhoodHypergraph Stats" << std::endl;
    oss << "----------------------------" << std::endl;
    oss << "#Hypernodes/#Hyperedges: " << _neighbors.size() << std::endl;
    size_t pins = 0;
    size_t min_pins = std::numeric_limits<size_t>::max();
    size_t max_pins = 0;
    for(int i = 0; i < _neighbors.size(); i++) {
      pins += _neighbors[i].size();
      min_pins = std::min(min_pins,_neighbors[i].size());
      max_pins = std::max(max_pins,_neighbors[i].size());
    }
    oss << "#Pins: " << pins << std::endl;
    oss << "Average Neighborhood size: " << static_cast<double>(pins)/_neighbors.size() << " pins" << std::endl;
    oss << "Maximum Neighborhood: " << max_pins << " pins"<< std::endl;
    oss << "Minimum Neighborhood: " << min_pins << " pins" << std::endl;
    return oss.str();
  }


  
private:
  
  void printNeighborhood(const HypernodeID u) {
    std::cout << "Neighborhood of HN " << u << std::endl;
    for(int i = 0; i < _neighbors[u].size(); i++) {
      std::cout << _neighbors[u][i] << " ";
    }
    std::cout << std::endl << "---------------------------------" << std::endl;
  }
  
  Hypergraph& _hg;
  std::vector<std::vector<int>> _neighbors;
  std::vector<bool> _is_active;

};

}  // namespace datastructure
#endif  // SRC_LIB_DATASTRUCTURE_NEIGHBORHOODHYPERGRAPH_H_
