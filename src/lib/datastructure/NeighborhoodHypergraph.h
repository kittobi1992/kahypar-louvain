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
    //Precondition: Validates if neighborhood list of x and y is sorted.
    ASSERT(std::is_sorted(_neighbors[x].begin(),_neighbors[x].end()), "Neighborhood list of HN " << x << " should be in sorted order!");
    ASSERT(std::is_sorted(_neighbors[y].begin(),_neighbors[y].end()), "Neighborhood list of HN " << y << " should be in sorted order!");
    
    
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
      //Let c = _neighbors[y][j], we store -(c+1) to denote that hypernode c is contained in both list.
      //Because hypernode id 0 exists, we add +1 to the id before we negate the value, to guarantee that
      //that the new value is strict negative.
      else {
	new_x.push_back(_neighbors[x][i++]);
	_neighbors[y][j] = -(_neighbors[y][j++]+1);
      }
    }
    
    //Validate if the values in the inactive hypernode y are correct
    ASSERT([&]() {
	for(int k = 0; k < _neighbors[y].size(); ++k) {
	  if(_neighbors[y][k] == y)
	    continue;
          int search_pin = _neighbors[y][k] < 0 ? std::abs(_neighbors[y][k]+1) : _neighbors[y][k];
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
      HypernodeID pin = static_cast<HypernodeID>(_neighbors[y][k] < 0 ? std::abs(_neighbors[y][k]+1) : _neighbors[y][k]);
      if(_is_active[pin] && pin != y && pin != x) {
	//Validates if neighborhood list of pin is sorted.
	ASSERT(std::is_sorted(_neighbors[pin].begin(),_neighbors[pin].end()), "Neighborhood list of HN " << pin << " should be in sorted order!");
	auto pos = std::lower_bound(_neighbors[pin].begin(),_neighbors[pin].end(),y);
	ASSERT(*pos == y, "Hypernode which should be replaced is not HN " << y << " (Actual: "<<*pos <<")!");
	*pos = x;
	swapToSortedPosition(_neighbors[pin],pos);
	//Validates if neighborhood list of pin is sorted.
	ASSERT(std::is_sorted(_neighbors[pin].begin(),_neighbors[pin].end()), "Neighborhood list of HN " << pin << " should be in sorted order!");
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
    ASSERT(std::is_sorted(_neighbors[x].begin(),_neighbors[x].end()), "Neighborhood list of HN " << x << " should be in sorted order!");
    //Validates if all neighborhood lists of HN x neighbors are sorted.
    ASSERT([&]() {
      for(int k = 0; k < _neighbors[x].size(); ++k) {
	HypernodeID pin = _neighbors[x][k];
	if(!std::is_sorted(_neighbors[pin].begin(),_neighbors[pin].end())) {
	  LOG("Neighborhood list of HN " << pin << " is unsorted!");
	  printNeighborhood(pin);
	  return false;
	}
      }
      return true;
    }(),"A neighborhood list of HN " << x << " neighbors is unsorted!");
    //Validates if hypernode x is inserted in each neighborhood list of neighbors of x.
    ASSERT(validateIfNeighborhoodContainHypernode(x),"HN " << x << " should inserted of each neighbors list!");
    //Validates if hypernode y is removed from neighborhood list of neighbors of y.
    ASSERT(validateIfNeighborhoodContainHypernode(y,true),"HN " << y << " should be removed from all neighbors!");
    
    
  }
  
  void uncontract(const HypernodeID x, const HypernodeID y) {
    ASSERT(x < _neighbors.size(), "HN " << x << " isn't a valid hypernode!");
    ASSERT(y < _neighbors.size(), "HN " << y << " isn't a valid hypernode!");
    ASSERT(_is_active[x] && !_is_active[y], "HN "<<x<<" has to be active and HN "<<y<<" has to be inactive to uncontract them!");
    //Precondition: Validates if neighborhood list of x is sorted.
    ASSERT(std::is_sorted(_neighbors[x].begin(),_neighbors[x].end()), "Neighborhood list of HN " << x << " should be in sorted order!"); 
    
    int i = 0, j = 0;
    int size_x = _neighbors[x].size(), size_y = _neighbors[y].size();
    std::vector<int> new_x;
    bool insert_y_into_x = false;
    while(i < size_x || j < size_y) {
      //If we reach the end of the list of hypernode y we can insert the remaining elements
      //from hypernode x at end of new_x
      //Note: We only have to consider this case, because every element which is in _neighbors[y] is also
      //in _neighbors[x], therefore we have the invariant _neighbors[x][i] <= _neighbors[y][i]
      //=> end of list y is reached before end of list x.
      if(j == size_y) {
	new_x.insert(new_x.end(),_neighbors[x].begin()+i,_neighbors[x].end());
	i = size_x;
	continue;
      }
      
      //Decode value in inactive list of hypernode y
      int cur_y = _neighbors[y][j] < 0 ? std::abs(_neighbors[y][j]+1) : _neighbors[y][j];
      if(cur_y == y) {
	j++; continue;
      }
     
      if(_neighbors[x][i] < cur_y) {
	new_x.push_back(_neighbors[x][i++]);
      }
      //If the values in both list are the same, cur_y is only element of the new neighborhood list
      //of HN x, if _neighbors[y][j] < 0 (see contraction).
      else if(_neighbors[x][i] == cur_y) {
	if(_neighbors[y][j] < 0) {
	  new_x.push_back(_neighbors[x][i++]);
	  _neighbors[y][j] = cur_y;
	  HypernodeID pin = _neighbors[y][j++];
	  //We have to add hypernode y again to the neighborhood list of HN pin.
	  if(_is_active[pin]) {
	    if(pin != x) {
	      _neighbors[pin].push_back(y);
	      swapToSortedPosition(_neighbors[pin],--_neighbors[pin].end());
	      ASSERT(std::is_sorted(_neighbors[pin].begin(),_neighbors[pin].end()), "Neighborhood list of HN " << pin << " should be in sorted order!");
	    }
	    else {
	      insert_y_into_x = true;
	    }
	  }
	}
	//Hypernode cur_y is not longer elements of the neighborhood list of HN x, therefore
	//we have to replace x in cur_y neighborhood with y.
	else {
	  HypernodeID pin = cur_y;
	  if(_is_active[pin]) {
	    auto pos = std::lower_bound(_neighbors[pin].begin(),_neighbors[pin].end(),x);
	    ASSERT(*pos == x, "Hypernode which should be replaced is not HN " << x << " ("<<*pos <<")!");
	    *pos = y;
	    swapToSortedPosition(_neighbors[pin],pos);
	    ASSERT(std::is_sorted(_neighbors[pin].begin(),_neighbors[pin].end()), "Neighborhood list of HN " << pin << " should be in sorted order!");
	  }
	  i++; j++;
	}
      }
    }
    
    //This case is treated sepereate, because if we insert y in the loop before directly into
    //HN x neighborhood it can occur that list is not sorted at the end.
    if(insert_y_into_x) {
      new_x.push_back(y);
      swapToSortedPosition(new_x,--new_x.end()); 
    }
    
    //Apply new neighborhood list of hypernode x
    swap(_neighbors[x],new_x);
    //Activate hypernode y
    _is_active[y] = true;
   
    //Validates if neighborhood list of x is sorted.
    ASSERT(std::is_sorted(_neighbors[x].begin(),_neighbors[x].end()), "Neighborhood list of HN " << x << " should be in sorted order!");
    //Validates if all neighborhood lists of HN x neighbors are sorted.
    ASSERT([&]() {
      for(int k = 0; k < _neighbors[x].size(); ++k) {
	HypernodeID pin = _neighbors[x][k];
	if(!std::is_sorted(_neighbors[pin].begin(),_neighbors[pin].end())) {
	  LOG("Neighborhood list of HN " << pin << " is unsorted!");
	  printNeighborhood(pin);
	  return false;
	}
      }
      return true;
    }(),"A neighborhood list of HN " << x << " neighbors is unsorted!");
    //Validates if neighborhood list of y is sorted.
    ASSERT(std::is_sorted(_neighbors[y].begin(),_neighbors[y].end()), "Neighborhood list of HN " << y << " should be in sorted order!");
    //Validates if all neighborhood lists of HN y neighbors are sorted.
    ASSERT([&]() {
      for(int k = 0; k < _neighbors[y].size(); ++k) {
	HypernodeID pin = _neighbors[y][k];
	if(!std::is_sorted(_neighbors[pin].begin(),_neighbors[pin].end())) {
	  LOG("Neighborhood list of HN " << pin << " is unsorted!");
	  printNeighborhood(pin);
	  return false;
	}
      }
      return true;
    }(),"A neighborhood list of HN " << x << " neighbors is unsorted!");
    //Validates if hypernode x is inserted in each neighborhood list of neighbors of x.
    ASSERT(validateIfNeighborhoodContainHypernode(x),"HN " << x << " should inserted of each neighbors list!");
    //Validates if hypernode y is inserted in each neighborhood list of neighbors of y.
    ASSERT(validateIfNeighborhoodContainHypernode(y),"HN " << y << " should inserted of each neighbors list!");
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

=======
  std::vector<std::vector<int>> _neighbors;

  
  void printNeighborhood(const HypernodeID u) {
    std::cout << "Neighborhood of HN " << u << std::endl;
    for(int i = 0; i < _neighbors[u].size(); i++) {
      std::cout << _neighbors[u][i] << " ";
    }
    std::cout << std::endl << "---------------------------------" << std::endl;
  }
  
private:
  
  void swapToSortedPosition(std::vector<int>& neighbor, std::vector<int>::iterator pos) {
	while((pos != neighbor.begin() && *pos <= *(pos-1)) 
	      ||  (pos+1 != neighbor.end() && *pos >= *(pos+1))) {
	    if(pos != neighbor.begin() && *pos < *(pos-1)) {
	      auto tmp = *pos; *pos = *(pos-1); *(pos-1) = tmp; pos--;
	    } else if(pos+1 != neighbor.end() && *pos > *(pos+1)) {
	      auto tmp = *pos; *pos = *(pos+1); *(pos+1) = tmp; pos++;
	    } else if(*pos == *(pos-1) || *pos == *(pos+1)) {
	      neighbor.erase(pos);
	      break;
	    }
	}    
  }
  
  
  bool validateIfNeighborhoodContainHypernode(HypernodeID x, bool removed = false) {
 	for(int k = 0; k < _neighbors[x].size(); ++k) {
	  HypernodeID pin = _neighbors[x][k] < 0 ? std::abs(_neighbors[x][k]+1) : _neighbors[x][k];
	  if(_is_active[pin]) {
	    auto find = std::lower_bound(_neighbors[pin].begin(),_neighbors[pin].end(),x);
	    if(removed) {
	      if(find != _neighbors[pin].end() && *find == x) {
		LOG("Neighborhood list of HN " << pin << " contains HN " << x << "!");
		printNeighborhood(pin);
		return false;
	      } 
	    }
	    else {
	      if(find == _neighbors[pin].end() || *find != x) {
		LOG("Neighborhood list of HN " << pin << " didn't contain HN " << x << "!");
		printNeighborhood(pin);
		return false;
	      } 
	    }
	  }
	}
        return true;   
  }
  
  std::vector<bool> _is_active;

};

}  // namespace datastructure
#endif  // SRC_LIB_DATASTRUCTURE_NEIGHBORHOODHYPERGRAPH_H_
