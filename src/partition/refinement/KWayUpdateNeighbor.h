/***************************************************************************
 *  Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_PARTITION_REFINEMENT_KWAYUPDATENEIGHBOR_H_
#define SRC_PARTITION_REFINEMENT_KWAYUPDATENEIGHBOR_H_

#include <memory>
#include <utility>
#include <vector>

#include "lib/definitions.h"
#include "lib/datastructure/SparseSet.h"
#include "lib/datastructure/KWayPriorityQueue.h"

using datastructure::SparseSet;
using datastructure::KWayPriorityQueue;

namespace partition {
 
class KWayUpdateNeighbor {
using Gain = HyperedgeWeight;
using Byte = char;
using KWayRefinementPQ = KWayPriorityQueue<HypernodeID, Gain,
                                             std::numeric_limits<Gain> >;

struct Element {
  PartitionID index;
  Gain gain;
  Element(const PartitionID i, const Gain g) :
    index(i),
    gain(g) { }
};  
  

private:
class UpdateElement {

public:
  UpdateElement(PartitionID k) : 
    _k(k), 
    _hypernodeIsConnectedToFromPart(kInvalidPart),
    _valid(0),
    _size(0) {  
     static_assert(sizeof(Gain) == sizeof(PartitionID), "Size is not correct");
     for (PartitionID i = 0; i < k; ++i) {
       new(&_size + 1 + i)PartitionID(std::numeric_limits<PartitionID>::max());
       new(reinterpret_cast<Element*>(&_size + _k + 1) + i)Element(kInvalidPart,0);
     }   
  }

  UpdateElement(const UpdateElement&) = delete;
  UpdateElement(UpdateElement&&) = delete;
  UpdateElement& operator= (const UpdateElement&) = delete;
  UpdateElement& operator= (UpdateElement&&) = delete;
  
  const PartitionID* begin()  const {
    return &_size + 1;
  }

  const PartitionID* end() const {
    ASSERT(_size <= _k, V(_size));
    return &_size + 1 + _size;
  }
  
  __attribute__ ((always_inline)) void update(PartitionID part, Gain delta) {
    ASSERT(part < _k, V(part));
    if(sparse(part).index == kInvalidPart) {
      sparse(part).index = _size;
      dense(_size++) = part;
    }
    sparse(part).gain += delta;
  }
  
  __attribute__ ((always_inline)) Gain gain(PartitionID part) {
    Gain gain = sparse(part).gain;
    sparse(part).gain = 0;
    sparse(part).index = kInvalidPart;
    return gain;
  }
  
  PartitionID isHypernodeConnectedToFromPart() {
    return _hypernodeIsConnectedToFromPart;
  }
  
  void setHypernodeIsConnectedToFromPart(bool hypernodeIsConnectedToFromPart) {
    _hypernodeIsConnectedToFromPart = hypernodeIsConnectedToFromPart;
  }
  
  __attribute__ ((always_inline)) void clear() {
    _size = 0;
    _hypernodeIsConnectedToFromPart = kInvalidPart;
  }
  
  __attribute__ ((always_inline)) bool makeValid(const size_t threshold) {
    _valid = threshold;
  }
  
  __attribute__ ((always_inline)) bool isValid(const size_t threshold) {
    return _valid >= threshold;
  }
  
private:
    // To avoid code duplication we implement non-const version in terms of const version
    PartitionID & dense(const PartitionID part) {
      return const_cast<PartitionID&>(static_cast<const UpdateElement&>(*this).dense(part));
    }

    const PartitionID & dense(const PartitionID part) const {
      ASSERT(part < _k, V(part));
      return *(&_size + 1 + part);
    }

    // To avoid code duplication we implement non-const version in terms of const version
    Element & sparse(const PartitionID part) {
      return const_cast<Element&>(static_cast<const UpdateElement&>(*this).sparse(part));
    }

    const Element & sparse(const PartitionID part) const {
      ASSERT(part < _k, V(part));
      return reinterpret_cast<const Element*>(&_size + _k + 1)[part];
    }
  
  PartitionID _k;
  PartitionID _hypernodeIsConnectedToFromPart;
  size_t _valid;
  PartitionID _size;
};

public:
  static constexpr Gain kNotCached = std::numeric_limits<Gain>::max();
  static constexpr PartitionID kInvalidPart = std::numeric_limits<PartitionID>::max();
  
  explicit KWayUpdateNeighbor(HypernodeID num_hypernodes, PartitionID k) :
    _sparse(nullptr),
    _dense(std::make_unique<HypernodeID[]>(num_hypernodes)),
    _k(k),
    _num_hypernodes(num_hypernodes),
    _threshold(1),
    _size(0) { 
      _sparse = static_cast<Byte*>(malloc(num_hypernodes * sizeOfUpdateElement()));
      for (HypernodeID hn = 0; hn < num_hypernodes; ++hn) {
	new(updateElement(hn))UpdateElement(k);
      }
    }

  KWayUpdateNeighbor(const KWayUpdateNeighbor&) = delete;
  KWayUpdateNeighbor& operator= (const KWayUpdateNeighbor&) = delete;

  KWayUpdateNeighbor(KWayUpdateNeighbor&&) = default;
  KWayUpdateNeighbor& operator= (KWayUpdateNeighbor&&) = default;
  
  ~KWayUpdateNeighbor() {
    // Since CacheElemment only contains PartitionIDs and these are PODs,
    // we do not need to call destructors of CacheElement cacheElement(i)->~CacheElement();
    static_assert(std::is_pod<Gain>::value, "Gain is not a POD");
    static_assert(std::is_pod<PartitionID>::value, "PartitionID is not a POD");
    free(_sparse);
  }
  
  HypernodeID* begin() const {
    return _dense.get();
  }

  HypernodeID* end() const {
    return _dense.get() + _size;
  }  

  __attribute__ ((always_inline)) void update(HypernodeID hn, PartitionID part, Gain delta) {
      UpdateElement* ue = updateElement(hn);
      ue->update(part,delta);
      //std::cout << "("<<hn<<","<<!ue->isValid(_threshold) << ","<<_size<<")" << std::endl;
      if(!ue->isValid(_threshold)) {
	_dense[_size++] = hn;
	ue->makeValid(_threshold);
      }
  }
  
  
  __attribute__ ((always_inline)) void updatePQ(KWayRefinementPQ& pq) {
    for(const HypernodeID* hn = begin(); hn != end(); hn++) {
     //std::cout << *hn << std::endl;
      UpdateElement* ue = updateElement(*hn);
      for(const PartitionID* part = ue->begin(); part != ue->end(); part++) {
	pq.updateKeyBy(*hn,*part,ue->gain(*part));
      }
      ue->clear();
    }
    _size = 0;
    _threshold++;
  }

private:
  
  const UpdateElement* updateElement(const HypernodeID hn) const {
    return reinterpret_cast<UpdateElement*>(_sparse + hn * sizeOfUpdateElement());
  }

  // To avoid code duplication we implement non-const version in terms of const version
  UpdateElement* updateElement(const HypernodeID hn) {
    return const_cast<UpdateElement*>(static_cast<const KWayUpdateNeighbor&>(*this).updateElement(hn));
  }
  
  size_t sizeOfUpdateElement() const {
    return sizeof(UpdateElement) + _k * sizeof(Element) + _k * sizeof(PartitionID);
  }
  
  Byte* _sparse;
  std::unique_ptr<HypernodeID[]> _dense;
  PartitionID _k;
  HypernodeID _num_hypernodes;
  size_t _threshold;
  size_t _size;
};

}  // namespace partition
#endif  //  SRC_PARTITION_REFINEMENT_KWAYUPDATENEIGHBOR_H_
