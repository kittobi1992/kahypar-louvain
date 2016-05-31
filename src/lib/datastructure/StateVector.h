/***************************************************************************
 *  Copyright (C) 2015 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/
#ifndef SRC_LIB_DATASTRUCTURE_STATEVECTOR_H_
#define SRC_LIB_DATASTRUCTURE_STATEVECTOR_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "lib/macros.h"

// based on http://upcoder.com/9/fast-resettable-flag-vector/

using std::size_t;
using std::uint16_t;

namespace datastructure {
template <typename UnderlyingType = std::int16_t,
	  size_t num_states = 1>
class StateVector {
 public:
  StateVector(const size_t size, const bool initialiser = false) :
    _v(std::make_unique<UnderlyingType[]>(size)),
    _threshold(1),
    _size(size) {
    memset(_v.get(), (initialiser ? 1 : 0), size * sizeof(UnderlyingType));
  }

  StateVector() :
    _v(nullptr),
    _threshold(0),
    _size(0) { }

  StateVector(const StateVector&) = delete;
  StateVector& operator= (const StateVector&) = delete;

  StateVector(StateVector&&) = default;
  StateVector& operator= (StateVector&&) = default;

  void swap(StateVector& other) noexcept {
    using std::swap;
    swap(_v, other._v);
    swap(_threshold, other._threshold);
  }

  __attribute__ ((always_inline)) UnderlyingType operator[] (const size_t i) const {
    UnderlyingType tmp = _v[i];
    return (tmp > _threshold) ? tmp - _threshold : 0;
  }
  
  __attribute__ ((always_inline)) bool isEntryValid(const size_t i) const {
    return _v[i] > _threshold;
  }
  
  __attribute__ ((always_inline)) void setState(const size_t i, const UnderlyingType state) const {
    ASSERT(value <= num_states, "Value is greater than number of states!");
    _v[i] = _threshold  + state;
  }


  __attribute__ ((always_inline)) void reset() {
    if (_threshold >= std::numeric_limits<UnderlyingType>::max() - num_states) {
      for (size_t i = 0; i != _size; ++i) {
        _v[i] = 0;
      }
      _threshold = 0;
    }
    else {
      _threshold += num_states;
    }
  }

  void setSize(const size_t size, const bool initialiser = false) {
    ASSERT(_v == nullptr, "Error");
    _v = std::make_unique<UnderlyingType[]>(size);
    _size = size;
    memset(_v.get(), (initialiser ? 1 : 0), size * sizeof(UnderlyingType));
  }

 private:

  std::unique_ptr<UnderlyingType[]> _v;
  UnderlyingType _threshold;
  size_t _size;
};

template <typename UnderlyingType>
void swap(StateVector<UnderlyingType>& a,
          StateVector<UnderlyingType>& b) noexcept {
  a.swap(b);
}
}  // namespace datastructure
#endif  // SRC_LIB_DATASTRUCTURE_STATEVECTOR_H_
