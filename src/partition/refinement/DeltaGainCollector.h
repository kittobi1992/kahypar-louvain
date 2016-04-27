/***************************************************************************
 *  Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_PARTITION_REFINEMENT_DELTAGAINCOLLECTOR_H_
#define SRC_PARTITION_REFINEMENT_DELTAGAINCOLLECTOR_H_

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "lib/core/Mandatory.h"
#include "lib/definitions.h"

using defs::PartitionID;
using defs::HypernodeID;

namespace partition {
  
template <typename Value = Mandatory>
class DeltaGainCollector {

using Key = std::pair<HypernodeID,PartitionID>;
  
 private:
  struct MapElement {
    Key key;
    Value value; 
    size_t valid_entry;

    MapElement() :
      key(0,0),
      value(0),
      valid_entry(0) {}
    
    MapElement(Key k, Value val) :
      key(k),
      value(val),
      valid_entry(0) { }

    MapElement(MapElement&&) = default;
    MapElement& operator= (MapElement&&) = default;
    
  };

 public:
  explicit DeltaGainCollector(HypernodeID num_hypernodes, PartitionID k) :
    _dense(std::make_unique<Key[]>(num_hypernodes*k)),
    _sparse(std::make_unique<MapElement[]>(num_hypernodes*k)),
    _k(k),
    _size(0),
    _threshold(1) { }

  DeltaGainCollector(const DeltaGainCollector&) = delete;
  DeltaGainCollector& operator= (const DeltaGainCollector&) = delete;

  DeltaGainCollector(DeltaGainCollector&&) = default;
  DeltaGainCollector& operator= (DeltaGainCollector&&) = default;

  void swap(DeltaGainCollector& other) noexcept {
    using std::swap;
    swap(_dense, other._dense);
    swap(_sparse, other._sparse);
    swap(_size, other._size);
  }
  

  bool contains(const Key key) {
    return isValidEntry(key);
  }

  Value& operator[] (const Key key) {
    size_t hash = getHash(key);
    if (!isValidEntry(key)) {
      _dense[_size++] = key;
      _sparse[hash].value = Value();
      _sparse[hash].valid_entry = _threshold;
      return _sparse[hash].value;
    }
    return _sparse[hash].value;
  }

  Value& get(const Key key) {
    size_t hash = getHash(key);
    ASSERT(contains(key), V(getHash(key)));
    return _sparse[hash].value;
  }

  void add(const Key key, const Value value) {
    size_t hash = getHash(key);
    if (!isValidEntry(key)) {
      _dense[_size++] = key;
      _sparse[hash].value = value;
      _sparse[hash].valid_entry = _threshold;
    }
  }
  
  void addDelta(const Key key, const Value delta) {
    size_t hash = getHash(key);
    if(!isValidEntry(key)) {
      _dense[_size++] = key;
      _sparse[hash].value = delta;
      _sparse[hash].valid_entry = _threshold;
    } else {
      _sparse[hash].value += delta;
    }
  }

  size_t size() const {
    return _size;
  }

  auto begin() const {
    return _dense.get();
  }

  auto end() const {
    return _dense.get()+_size;
  }

  auto crbegin() const {
    return _dense;
  }

  auto* crend() const {
    return _dense.get() + _size;
  }

  void clear() {
    _size = 0;
    _threshold++;
  }
  
     
  bool isValidEntry(Key key) {
    return _sparse[getHash(key)].valid_entry == _threshold;
  }
   

 private:
  size_t getHash(const Key key) {
    return key.first*_k + key.second;
  }
   
  std::unique_ptr<Key[]> _dense;
  std::unique_ptr<MapElement[]> _sparse;
  PartitionID _k;
  size_t _size;
  size_t _threshold;
};

template <typename Value>
void swap(DeltaGainCollector<Value>& a,
          DeltaGainCollector<Value>& b) noexcept {
  a.swap(b);
}
}  // namespace partition
#endif  //  SRC_PARTITION_REFINEMENT_DELTAGAINCOLLECTOR_H_
