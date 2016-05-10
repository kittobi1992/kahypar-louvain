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

using defs::Hypergraph;
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
    _hypernodeIsConnectedToFromPart(0) {  
     static_assert(sizeof(Gain) == sizeof(PartitionID), "Size is not correct");
     for (PartitionID i = 0; i < k; ++i) {
       new(&_hypernodeIsConnectedToFromPart + 1 + i)Gain(kInvalidGain);
     }   
  }

  UpdateElement(const UpdateElement&) = delete;
  UpdateElement(UpdateElement&&) = delete;
  UpdateElement& operator= (const UpdateElement&) = delete;
  UpdateElement& operator= (UpdateElement&&) = delete;
  
  __attribute__ ((always_inline)) PartitionID setConnectedToFromPart(PartitionID connectedToFromPart) {
    _hypernodeIsConnectedToFromPart = connectedToFromPart;
  }
  
   __attribute__ ((always_inline)) PartitionID isConnectedToFromPart() {
    return _hypernodeIsConnectedToFromPart;
  }
  
  __attribute__ ((always_inline)) bool update(PartitionID part, Gain delta) {
    ASSERT(part < _k, V(part));
    Gain gain = sparse(part);
    sparse(part) = (gain == kInvalidPart ? delta : gain + delta);
    return gain == kInvalidPart;
  }
  
  __attribute__ ((always_inline)) Gain gain(PartitionID part) {
    Gain gain = sparse(part);
    sparse(part) = kInvalidGain;
    return gain;
  }
  
  
private:
    // To avoid code duplication we implement non-const version in terms of const version
    Gain & sparse(const PartitionID part) {
      return const_cast<Gain&>(static_cast<const UpdateElement&>(*this).sparse(part));
    }

    const Gain & sparse(const PartitionID part) const {
      ASSERT(part < _k, V(part));
      return *(&_hypernodeIsConnectedToFromPart + 1 + part);
    }


  PartitionID _k;
  PartitionID _hypernodeIsConnectedToFromPart;
};

public:
  static constexpr Gain kInvalidGain = std::numeric_limits<Gain>::max();
  static constexpr PartitionID kInvalidPart = std::numeric_limits<PartitionID>::max();
  
  explicit KWayUpdateNeighbor(HypernodeID num_hypernodes, PartitionID k) :
    _sparse(nullptr),
    _dense(std::make_unique<std::pair<HypernodeID,PartitionID>[]>(num_hypernodes)),
    _k(k),
    _num_hypernodes(num_hypernodes),
    _size(0),
    _threshold(2) { 
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
  
  std::pair<HypernodeID,PartitionID>* begin() const {
    return _dense.get();
  }

  std::pair<HypernodeID,PartitionID>* end() const {
    return _dense.get() + _size;
  }  

  __attribute__ ((always_inline)) bool hypernodeIsConnectedToPart(const Hypergraph& hg, 
								  const HypernodeID pin, 
								  const PartitionID part) {
    UpdateElement* ue = updateElement(pin);
    PartitionID isConnectedToFromPart = ue->isConnectedToFromPart();
    if(isConnectedToFromPart > _threshold-2) {
      return isConnectedToFromPart == _threshold;
    }
    for (const HyperedgeID he : hg.incidentEdges(pin)) {
      if (hg.pinCountInPart(he, part) > 0) {
	ue->setConnectedToFromPart(_threshold);
        return true;
      }
    }
    ue->setConnectedToFromPart(_threshold-1);
    return false;
  }
  
  __attribute__ ((always_inline)) void update(HypernodeID hn, PartitionID part, Gain delta) {
      UpdateElement* ue = updateElement(hn);
      if(ue->update(part,delta)) {
	_dense[_size++] = std::make_pair(hn,part);
      }
  }
  
  __attribute__ ((always_inline)) Gain deltaGain(const HypernodeID hn, const PartitionID part) {
    return updateElement(hn)->gain(part);
  }
  
  
  __attribute__ ((always_inline)) void clear() {
    _size = 0;
    _threshold += 2;
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
    return sizeof(UpdateElement) + _k * sizeof(Gain);
  }
  
  Byte* _sparse;
  std::unique_ptr<std::pair<HypernodeID,PartitionID>[]> _dense;
  PartitionID _k;
  HypernodeID _num_hypernodes;
  size_t _size;
  PartitionID _threshold;
};

}  // namespace partition
#endif  //  SRC_PARTITION_REFINEMENT_KWAYUPDATENEIGHBOR_H_
