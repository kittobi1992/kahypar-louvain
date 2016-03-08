/***************************************************************************
 *  Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_LIB_DATASTRUCTURE_NEIGHBORHOODHYPERGRAPH_H_
#define SRC_LIB_DATASTRUCTURE_NEIGHBORHOODHYPERGRAPH_H_

#include <vector>
#include <algorithm>
#include <sstream>

#include "lib/core/Mandatory.h"
#include "lib/datastructure/SparseSet.h"
#include "lib/definitions.h"

using HypernodeID = std::uint32_t;

namespace datastructure {

class NeighborhoodHypergraph {
 public:
  explicit NeighborhoodHypergraph(size_t initial_num_nodes) :
      _neighbors(initial_num_nodes),
      _is_active(initial_num_nodes,true) { }
      
  NeighborhoodHypergraph() :
      _neighbors(0),
      _is_active(0) { }
  
  
  void setNeighborhoodOfHypernode(HypernodeID hn, std::vector<int>& neighbor) {
    _neighbors[hn].clear();
    for(int i = 0; i < neighbor.size(); i++) {
      _neighbors[hn].push_back(neighbor[i]);
    }
  }
  
  void resize(size_t initial_num_nodes) {
    _neighbors.clear();
    _is_active.clear();
    _neighbors.resize(initial_num_nodes);
    _is_active.assign(initial_num_nodes,true);
  }
  
  void contract(const HypernodeID x, const HypernodeID y) {
    ASSERT(x < _neighbors.size(), "HN " << x << " isn't a valid hypernode!");
    ASSERT(y < _neighbors.size(), "HN " << y << " isn't a valid hypernode!");
    ASSERT(_is_active[x] && _is_active[y], "Both hypernodes has to be active before contraction!");
    
    int i = 0, j = 0;
    int size_x = _neighbors[x].size(), size_y = _neighbors[y].size();
    std::vector<int> new_x;
    //Merging the two neighborhoods vector of hypernode x and y
    while(i < size_x || j < size_y) {
      
      //Ignoring the contraction partner of y in both list
      if(i < size_x && _neighbors[x][i] == y) {
	i++; continue;
      } else if(j < size_y && _neighbors[y][j] == y) {
	j++; continue;
      }
      
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
    
    //Replace hypernode y in all neighborhood lists of y's neighbors with hypernode x
    for(int k = 0; k < _neighbors[y].size(); ++k) {
      HypernodeID pin = static_cast<HypernodeID>(std::abs(_neighbors[y][k]));
      if(_is_active[pin] && pin != y && pin != x) {
	//Validates if neighborhood list of pin is sorted.
	ASSERT(isSorted(pin), "Neighborhood list of HN " << pin << " should be in sorted order!");
	auto pos = std::lower_bound(_neighbors[pin].begin(),_neighbors[pin].end(),y);
	ASSERT(*pos == y, "Hypernode which should be replaced is not HN " << y << " ("<<*pos <<")!");
	*pos = x; bool remove = false;
	while((pos != _neighbors[pin].begin() && *pos <= *(pos-1)) 
	      ||  (pos+1 != _neighbors[pin].end() && *pos >= *(pos+1))) {
	    if(pos != _neighbors[pin].begin() && *pos < *(pos-1)) {
	      auto tmp = *pos; *pos = *(pos-1); *(pos-1) = tmp; pos--;
	    } else if(pos+1 != _neighbors[pin].end() && *pos > *(pos+1)) {
	      auto tmp = *pos; *pos = *(pos+1); *(pos+1) = tmp; pos++;
	    } else if(*pos == *(pos-1) || *pos == *(pos+1)) {
	      _neighbors[pin].erase(pos);
	      break;
	    }
	}
	//Validates if neighborhood list of pin is sorted.
	ASSERT(isSorted(pin), "Neighborhood list of HN " << pin << " should be in sorted order!");
      }
    }
    
    //Apply new neighborhood list of hypernode x
    swap(_neighbors[x],new_x);
    //Disable hypernode y
    _is_active[y] = false;
    
    auto pos = std::lower_bound(_neighbors[x].begin(),_neighbors[x].end(),y);
    if(pos != _neighbors[x].end() && *pos == y) {
     _neighbors[x].erase(pos);
    }
   
    //Validates if neighborhood list of x is sorted.
    ASSERT(isSorted(x), "Neighborhood list of HN " << x << " should be in sorted order!");

    
    
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
    
    
  }
  
  void uncontract(const HypernodeID x, const HypernodeID y) {
    ASSERT(x < _neighbors.size(), "HN " << x << " isn't a valid hypernode!");
    ASSERT(y < _neighbors.size(), "HN " << y << " isn't a valid hypernode!");
    ASSERT(_is_active[x] && !_is_active[y], "HN "<<x<<" has to be active and HN "<<y<<" has to be inactive to uncontract them!");
    
    int i = 0, j = 0;
    int size_x = _neighbors[x].size(), size_y = _neighbors[y].size();
    std::vector<int> new_x;
    while(i < size_x || j < size_y) {
      if(j == size_y) {
	new_x.insert(new_x.end(),_neighbors[x].begin()+i,_neighbors[x].end());
	i = size_x;
	continue;
      }
      
      int cur_y = std::abs(_neighbors[y][j]);
      if(cur_y == y) {
	j++; continue;
      }
      if(_neighbors[x][i] < cur_y) {
	new_x.push_back(_neighbors[x][i++]);
      }
      else if(_neighbors[x][i] == cur_y) {
	if(_neighbors[y][j] < 0) {
	  new_x.push_back(_neighbors[x][i++]);
	  _neighbors[y][j] = -_neighbors[y][j++];
	}
	else {
	  i++; j++;
	}
      }
    }
    
    //Apply new neighborhood list of hypernode x
    swap(_neighbors[x],new_x);
    //Activate hypernode y
    _is_active[y] = true;
    
    for(int k = 0; k < _neighbors[y].size(); ++k) {
      HypernodeID pin = _neighbors[y][k];
      if(_is_active[pin] && pin != y) {
	//TODO(heuer): Naive insert variant. Another approach is to use lists instead of vectors,
	//but this makes compression very difficult.
	auto pos = std::lower_bound(_neighbors[pin].begin(),_neighbors[pin].end(),y);
	_neighbors[pin].insert(pos,y);
      }
    }
    
    //Validates if hypernode y is inserted in each neighborhood list of neighbors of y.
    ASSERT([&]() {
	for(int k = 0; k < _neighbors[y].size(); ++k) {
	  HypernodeID pin = static_cast<HypernodeID>(std::abs(_neighbors[y][k]));
	  if(_is_active[pin]) {
	    auto find = std::lower_bound(_neighbors[pin].begin(),_neighbors[pin].end(),y);
	    if(find == _neighbors[pin].end() || *find != y) {
		LOG("Neighborhood list of HN " << pin << " didn't contain HN " << y << "!");
		printNeighborhood(pin);
		return false;
	    }
	  }
	}
        return true;
      } (), "HN " << y << " should be removed from all neighbors!"); 
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

  std::vector<std::vector<int>> _neighbors;

  
  void printNeighborhood(const HypernodeID u) {
    std::cout << "Neighborhood of HN " << u << std::endl;
    for(int i = 0; i < _neighbors[u].size(); i++) {
      std::cout << _neighbors[u][i] << " ";
    }
    std::cout << std::endl << "---------------------------------" << std::endl;
  }
  
private:
  
  bool isSorted(HypernodeID x) {
	for(int k = 1; k < _neighbors[x].size(); ++k) {
	  if(_neighbors[x][k-1] > _neighbors[x][k]) {
	    printNeighborhood(x);
	    return false;
	  }
	}
        return true;
  }
  
  
  std::vector<bool> _is_active;

};

}  // namespace datastructure
#endif  // SRC_LIB_DATASTRUCTURE_NEIGHBORHOODHYPERGRAPH_H_
