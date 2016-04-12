/***************************************************************************
 *  Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_LIB_DATASTRUCTURE_GENERICHYPERGRAPH_H_
#define SRC_LIB_DATASTRUCTURE_GENERICHYPERGRAPH_H_

#include <algorithm>
#include <bitset>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "partition/ConfigurationEnumClasses.h"

#include "gtest/gtest_prod.h"
#include "lib/core/Empty.h"
#include "lib/core/Int2Type.h"
#include "lib/core/Mandatory.h"
#include "lib/datastructure/SparseSet.h"
#include "lib/datastructure/ConnectivitySets.h"
#include "lib/datastructure/FastResetBitVector.h"
#include "lib/datastructure/NeighborhoodHypergraph.h"
#include "lib/definitions.h"
#include "lib/macros.h"


using core::Empty;
using core::Int2Type;
using partition::RefinementAlgorithm;

struct UncontractionGainChanges;

namespace datastructure {
template <typename Iterator>
Iterator begin(std::pair<Iterator, Iterator>& x) {
  return x.first;
}

template <typename Iterator>
Iterator end(std::pair<Iterator, Iterator>& x) {
  return x.second;
}

static const bool dbg_hypergraph_uncontraction = false;
static const bool dbg_hypergraph_contraction = false;
static const bool dbg_hypergraph_restore_edge = false;

template <typename HypernodeType_ = Mandatory,
          typename HyperedgeType_ = Mandatory,
          typename HypernodeWeightType_ = Mandatory,
          typename HyperedgeWeightType_ = Mandatory,
          typename PartitionIDType_ = Mandatory,
          class HypernodeData_ = Empty,
          class HyperedgeData_ = Empty
          >
class GenericHypergraph {
  // export template parameters

 public:
  using HypernodeID = HypernodeType_;
  using HyperedgeID = HyperedgeType_;
  using PartitionID = PartitionIDType_;
  using HypernodeWeight = HypernodeWeightType_;
  using HyperedgeWeight = HyperedgeWeightType_;
  using HypernodeData = HypernodeData_;
  using HyperedgeData = HyperedgeData_;

  enum { kInvalidPartition = -1,
         kDeletedPartition = -2 };

  enum class Type : int8_t {
    Unweighted = 0,
    EdgeWeights = 1,
    NodeWeights = 10,
    EdgeAndNodeWeights = 11,
  };

 private:
  // forward delarations
  struct Memento;
  struct PartInfo;
  template <typename T, class D>
  class InternalVertex;
  template <typename T>
  class VertexIterator;
  struct HypernodeTraits;
  struct HyperedgeTraits;

  struct AdditionalHyperedgeData : public HyperedgeData {
    PartitionID connectivity = 0;
  };

  struct AdditionalHypernodeData : public HypernodeData {
    PartitionID part_id = kInvalidPartition;
    HyperedgeID num_incident_cut_hes = 0;
    std::uint32_t state = 0;
  };

  struct Dummy {
    void push_back(HypernodeID) { }
  };

  // internal
  using VertexID = unsigned int;
  using HypernodeVertex = InternalVertex<HypernodeTraits, AdditionalHypernodeData>;
  using HyperedgeVertex = InternalVertex<HyperedgeTraits, AdditionalHyperedgeData>;
  using PinHandleIterator = typename std::vector<VertexID>::iterator;

 public:
  using HyperedgeIndexVector = std::vector<size_t>;
  using HyperedgeVector = std::vector<HypernodeID>;
  using HypernodeWeightVector = std::vector<HypernodeWeight>;
  using HyperedgeWeightVector = std::vector<HyperedgeWeight>;
  using ContractionMemento = Memento;
  using IncidenceIterator = typename std::vector<VertexID>::const_iterator;
  using HypernodeIterator = VertexIterator<const std::vector<HypernodeVertex> >;
  using HyperedgeIterator = VertexIterator<const std::vector<HyperedgeVertex> >;

 private:
  static const HypernodeID kInvalidCount = std::numeric_limits<HypernodeID>::max();

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
  template <typename VertexTypeTraits, class InternalVertexData>
  class InternalVertex : public InternalVertexData {
 public:
    using WeightType = typename VertexTypeTraits::WeightType;
    using IDType = typename VertexTypeTraits::IDType;

    InternalVertex(const IDType begin, const IDType size,
                   const WeightType weight) noexcept :
      _begin(begin),
      _size(size),
      _weight(weight),
      _valid(true) { }

    InternalVertex() noexcept :
      _begin(0),
      _size(0),
      _weight(1),
      _valid(true) { }

    InternalVertex(const InternalVertex&) = default;
    InternalVertex& operator= (const InternalVertex&) = default;

    InternalVertex(InternalVertex&&) = default;
    InternalVertex& operator= (InternalVertex&&) = default;

    void disable() noexcept {
      ASSERT(!isDisabled(), "Vertex is already disabled");
      _valid = false;
    }

    bool isDisabled() const noexcept {
      return _valid == false;
    }

    void enable() noexcept {
      ASSERT(isDisabled(), "Vertex is already enabled");
      _valid = true;
    }

    IDType firstEntry() const noexcept {
      return _begin;
    }

    void setFirstEntry(IDType begin) noexcept {
      ASSERT(!isDisabled(), "Vertex is disabled");
      _begin = begin;
      _valid = true;
    }

    IDType firstInvalidEntry() const noexcept {
      return _begin + _size;
    }

    IDType size() const noexcept {
      return _size;
    }
    void setSize(IDType size) noexcept {
      ASSERT(!isDisabled(), "Vertex is disabled");
      _size = size;
    }

    void increaseSize() noexcept {
      ASSERT(!isDisabled(), "Vertex is disabled");
      ++_size;
    }

    void decreaseSize() noexcept {
      ASSERT(!isDisabled(), "Vertex is disabled");
      ASSERT(_size > 0, "Size out of bounds");
      --_size;
    }

    WeightType weight() const noexcept {
      return _weight;
    }

    void setWeight(WeightType weight) noexcept {
      ASSERT(!isDisabled(), "Vertex is disabled");
      _weight = weight;
    }

    bool operator== (const InternalVertex& rhs) const noexcept {
      return _begin == rhs._begin && _size == rhs._size && _weight == rhs._weight;
    }

    bool operator!= (const InternalVertex& rhs) const noexcept {
      return !operator== (this, rhs);
    }

    bool operator< (const InternalVertex& rhs) const noexcept {
      return _begin < rhs._begin;
    }

    bool operator> (const InternalVertex& rhs) const noexcept {
      return operator< (rhs, this);
    }

    bool operator<= (const InternalVertex& rhs) const noexcept {
      return !operator> (this, rhs);
    }

    bool operator>= (const InternalVertex& rhs) const noexcept {
      return !operator< (this, rhs);
    }

 private:
    IDType _begin;
    IDType _size;
    WeightType _weight;
    bool _valid;
  };
#pragma GCC diagnostic pop

  template <typename ContainerType>
  class VertexIterator {
    using IDType = typename ContainerType::value_type::IDType;
    using ConstPointer = typename ContainerType::const_pointer;

 public:
    VertexIterator(const VertexIterator& other) = default;
    VertexIterator& operator= (const VertexIterator& other) = default;

    VertexIterator(VertexIterator&& other) = default;
    VertexIterator& operator= (VertexIterator&& other) = default;

    VertexIterator() noexcept :
      _id(0),
      _max_id(0),
      _vertex(nullptr) { }

    VertexIterator(ConstPointer start_vertex, IDType id, IDType max_id) noexcept :
      _id(id),
      _max_id(max_id),
      _vertex(start_vertex) {
      if (_id != _max_id && _vertex->isDisabled()) {
        operator++ ();
      }
    }

    IDType operator* () const noexcept {
      return _id;
    }

    VertexIterator& operator++ () noexcept {
      ASSERT(_id < _max_id, "Hypernode iterator out of bounds");
      do {
        ++_id;
        ++_vertex;
      } while (_id < _max_id && _vertex->isDisabled());
      return *this;
    }

    VertexIterator operator++ (int) noexcept {
      VertexIterator copy = *this;
      operator++ ();
      return copy;
    }

    friend VertexIterator end<>(std::pair<VertexIterator, VertexIterator>& x);
    friend VertexIterator begin<>(std::pair<VertexIterator, VertexIterator>& x);

    VertexIterator& operator-- () noexcept {
      ASSERT(_id > 0, "Hypernode iterator out of bounds");
      do {
        --_id;
        --_vertex;
      } while (_id > 0 && _vertex->isDisabled());
      return *this;
    }

    VertexIterator operator-- (int) noexcept {
      VertexIterator copy = *this;
      operator-- ();
      return copy;
    }

    bool operator!= (const VertexIterator& rhs) noexcept {
      return _id != rhs._id;
    }

 private:
    IDType _id;
    const IDType _max_id;
    ConstPointer _vertex;
  };

  struct Memento {
    Memento(const Memento& other) = delete;
    Memento& operator= (const Memento& other) = delete;

    Memento(Memento&& other) = default;
    Memento& operator= (Memento&& other) = default;

    Memento(const HypernodeID u_, const HypernodeID u_first_entry_,
            const HypernodeID u_size_, const HypernodeID v_) noexcept :
      u(u_),
      u_first_entry(u_first_entry_),
      u_size(u_size_),
      v(v_) { }

    const HypernodeID u;
    const HypernodeID u_first_entry;
    const HypernodeID u_size;
    const HypernodeID v;
  };

  struct PartInfo {
    HypernodeWeight weight;
    HypernodeID size;

    bool operator== (const PartInfo& other) const {
      return weight == other.weight && size == other.size;
    }
  };

  struct HypernodeTraits {
    using WeightType = HypernodeWeight;
    using IDType = HypernodeID;
  };

  struct HyperedgeTraits {
    using WeightType = HyperedgeWeight;
    using IDType = HyperedgeID;
  };

 public:
  GenericHypergraph(const HypernodeID num_hypernodes,
                    const HyperedgeID num_hyperedges,
                    const HyperedgeIndexVector& index_vector,
                    const HyperedgeVector& edge_vector,
                    const PartitionID k = 2,
                    const HyperedgeWeightVector* hyperedge_weights = nullptr,
                    const HypernodeWeightVector* hypernode_weights = nullptr) noexcept :
    _num_hypernodes(num_hypernodes),
    _num_hyperedges(num_hyperedges),
    _num_pins(edge_vector.size()),
    _total_weight(0),
    _k(k),
    _type(Type::Unweighted),
    _current_num_hypernodes(_num_hypernodes),
    _current_num_hyperedges(_num_hyperedges),
    _current_num_pins(_num_pins),
    _threshold_active(1),
    _threshold_marked(2),
    _hypernodes(_num_hypernodes, HypernodeVertex(0, 0, 1)),
    _hyperedges(_num_hyperedges, HyperedgeVertex(0, 0, 1)),
    _incidence_array(2 * _num_pins, 0),
    _part_info(_k),
    _pins_in_part(_num_hyperedges * k),
    _connectivity_sets(_num_hyperedges, k),
        _hes_not_containing_u(_num_hyperedges, false),
  _n_hypergraph() {
    VertexID edge_vector_index = 0;
    for (HyperedgeID i = 0; i < _num_hyperedges; ++i) {
      hyperedge(i).setFirstEntry(edge_vector_index);
      for (VertexID pin_index = index_vector[i]; pin_index < index_vector[i + 1]; ++pin_index) {
        hyperedge(i).increaseSize();
        _incidence_array[pin_index] = edge_vector[pin_index];
        hypernode(edge_vector[pin_index]).increaseSize();
        ++edge_vector_index;
      }
      // _pins_in_part[i].resize(_k);
    }

    hypernode(0).setFirstEntry(_num_pins);
    for (HypernodeID i = 0; i < _num_hypernodes - 1; ++i) {
      hypernode(i + 1).setFirstEntry(hypernode(i).firstInvalidEntry());
      hypernode(i).setSize(0);
    }
    hypernode(num_hypernodes - 1).setSize(0);

    for (HyperedgeID i = 0; i < _num_hyperedges; ++i) {
      for (VertexID pin_index = index_vector[i]; pin_index < index_vector[i + 1]; ++pin_index) {
        const HypernodeID pin = edge_vector[pin_index];
        _incidence_array[hypernode(pin).firstInvalidEntry()] = i;
        hypernode(pin).increaseSize();
      }
    }

    bool has_hyperedge_weights = false;
    if (hyperedge_weights != nullptr && !hyperedge_weights->empty()) {
      has_hyperedge_weights = true;
      for (HyperedgeID i = 0; i < _num_hyperedges; ++i) {
        hyperedge(i).setWeight((*hyperedge_weights)[i]);
      }
    }

    bool has_hypernode_weights = false;
    if (hypernode_weights != nullptr && !hypernode_weights->empty()) {
      has_hypernode_weights = true;
      for (HypernodeID i = 0; i < _num_hypernodes; ++i) {
        hypernode(i).setWeight((*hypernode_weights)[i]);
        _total_weight += (*hypernode_weights)[i];
      }
    } else {
      _total_weight = _num_hypernodes;
    }

    if (has_hyperedge_weights && has_hypernode_weights) {
      _type = Type::EdgeAndNodeWeights;
    } else if (has_hyperedge_weights) {
      _type = Type::EdgeWeights;
    } else if (has_hypernode_weights) {
      _type = Type::NodeWeights;
    }
    
    createNeighborhoodHypergraph();
   
  }

  GenericHypergraph(const GenericHypergraph&) = delete;
  GenericHypergraph& operator= (const GenericHypergraph&) = delete;

  GenericHypergraph(GenericHypergraph&&) = default;
  GenericHypergraph& operator= (GenericHypergraph&&) = delete;

  void printHyperedgeInfo() const {
    for (HyperedgeID i = 0; i < _num_hyperedges; ++i) {
      if (!hyperedge(i).isDisabled()) {
        std::cout << "hyperedge " << i << ": begin=" << hyperedge(i).firstEntry() << " size="
        << hyperedge(i).size() << " weight=" << hyperedge(i).weight() << std::endl;
      }
    }
  }

  void printHypernodeInfo() const {
    for (HypernodeID i = 0; i < _num_hypernodes; ++i) {
      if (!hypernode(i).isDisabled()) {
        std::cout << "hypernode " << i << ": begin=" << hypernode(i).firstEntry() << " size="
        << hypernode(i).size() << " weight=" << hypernode(i).weight() << std::endl;
      }
    }
  }

  void printIncidenceArray() const {
    for (VertexID i = 0; i < _incidence_array.size(); ++i) {
      std::cout << "_incidence_array[" << i << "]=" << _incidence_array[i] << std::endl;
    }
  }

  void printHyperedges() const {
    std::cout << "Hyperedges:" << std::endl;
    for (HyperedgeID i = 0; i < _num_hyperedges; ++i) {
      if (!hyperedge(i).isDisabled()) {
        printEdgeState(i);
      }
    }
  }

  void printHypernodes() const {
    std::cout << "Hypernodes:" << std::endl;
    for (HypernodeID i = 0; i < _num_hypernodes; ++i) {
      if (!hypernode(i).isDisabled()) {
        printNodeState(i);
      }
    }
  }

  void printGraphState() const {
    printHypernodeInfo();
    printHyperedgeInfo();
    printHypernodes();
    printHyperedges();
    // printIncidenceArray();
  }

  void printEdgeState(const HyperedgeID e) const {
    if (!hyperedge(e).isDisabled()) {
      std::cout << "HE " << e << " w= " << edgeWeight(e) << "c=" << connectivity(e) << ": ";
      for (const HypernodeID pin : pins(e)) {
        std::cout << pin << " ";
      }
      for (PartitionID i = 0; i != _k; ++i) {
        std::cout << " Part[" << i << "] =" << pinCountInPart(e, i);
      }
    } else {
      std::cout << e << " -- invalid --";
    }
    std::cout << std::endl;
  }

  void printNodeState(const HypernodeID u) const {
    if (!hypernode(u).isDisabled()) {
      std::cout << "HN " << u << " (" << hypernode(u).part_id << "): ";
      for (const HyperedgeID he : incidentEdges(u)) {
        std::cout << he << " ";
      }
    } else {
      std::cout << u << " -- invalid --";
    }
    std::cout << std::endl;
  }

  std::pair<IncidenceIterator, IncidenceIterator> incidentEdges(const HypernodeID u) const noexcept {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode " << u << " is disabled");
    return std::make_pair(_incidence_array.cbegin() + hypernode(u).firstEntry(),
                          _incidence_array.cbegin() + hypernode(u).firstInvalidEntry());
  }

  std::pair<IncidenceIterator, IncidenceIterator> pins(const HyperedgeID e) const noexcept {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge " << e << " is disabled");
    return std::make_pair(_incidence_array.cbegin() + hyperedge(e).firstEntry(),
                          _incidence_array.cbegin() + hyperedge(e).firstInvalidEntry());
  }

  std::pair<HypernodeIterator, HypernodeIterator> nodes() const noexcept {
    return std::make_pair(HypernodeIterator(_hypernodes.data(), 0, _num_hypernodes),
                          HypernodeIterator((_hypernodes.data() + _num_hypernodes),
                                            _num_hypernodes, _num_hypernodes));
  }

  std::pair<HyperedgeIterator, HyperedgeIterator> edges() const noexcept {
    return std::make_pair(HyperedgeIterator(_hyperedges.data(), 0, _num_hyperedges),
                          HyperedgeIterator((_hyperedges.data() + _num_hyperedges),
                                            _num_hyperedges, _num_hyperedges));
  }

  const typename ConnectivitySets<PartitionID, HyperedgeID>::ConnectivitySet &
  connectivitySet(const HyperedgeID he) const noexcept {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
    return _connectivity_sets[he];
  }

  Memento contract(const HypernodeID u, const HypernodeID v) noexcept {
    using std::swap;
    ASSERT(!hypernode(u).isDisabled(), "Hypernode " << u << " is disabled");
    ASSERT(!hypernode(v).isDisabled(), "Hypernode " << v << " is disabled");
    ASSERT(partID(u) == partID(v), "Hypernodes " << u << " & " << v << " are in different parts: "
           << partID(u) << " & " << partID(v));

    DBG(dbg_hypergraph_contraction, "contracting (" << u << "," << v << ")");

    hypernode(u).setWeight(hypernode(u).weight() + hypernode(v).weight());
    const HypernodeID u_offset = hypernode(u).firstEntry();
    const HypernodeID u_size = hypernode(u).size();

    // The first call to connectHyperedgeToRepresentative copies the old incidence array of the
    // representative u to the end of the _incidence_array (even in the case that the old entries
    // are already located at the end, as the result of a previous contraction operation) and then
    // appends the new hyperedge. The variable is then set to false and succesive calls within the
    // same contraction operation just append at the end of the _incidence_array.
    // This behavior is necessary in order to be able to use the old entries during uncontraction.
    bool first_call = true;

    PinHandleIterator slot_of_u, last_pin_slot;
    PinHandleIterator pins_begin, pins_end;
    // Use index-based iteration because case 2 might lead to reallocation!
    for (HyperedgeID he_it = hypernode(v).firstEntry(); he_it != hypernode(v).firstInvalidEntry();
         ++he_it) {
      std::tie(pins_begin, pins_end) = pinHandles(_incidence_array[he_it]);
      ASSERT(pins_begin != pins_end, "Hyperedge " << _incidence_array[he_it] << " is empty");
      slot_of_u = last_pin_slot = pins_end - 1;
      for (PinHandleIterator pin_iter = pins_begin; pin_iter != last_pin_slot; ++pin_iter) {
        if (*pin_iter == v) {
          swap(*pin_iter, *last_pin_slot);
          --pin_iter;
        } else if (*pin_iter == u) {
          slot_of_u = pin_iter;
        }
      }

      ASSERT(*last_pin_slot == v, "v is not last entry in adjacency array!");

      if (slot_of_u != last_pin_slot) {
        // Case 1:
        // Hyperedge e contains both u and v. Thus we don't need to connect u to e and
        // can just cut off the last entry in the edge array of e that now contains v.
        hyperedge(_incidence_array[he_it]).decreaseSize();
        if (partID(v) != kInvalidPartition) {
          decreasePinCountInPart(_incidence_array[he_it], partID(v));
        }
        --_current_num_pins;
      } else {
        // Case 2:
        // Hyperedge e does not contain u. Therefore we  have to connect e to the representative u.
        // This reuses the pin slot of v in e's incidence array (i.e. last_pin_slot!)
        connectHyperedgeToRepresentative(_incidence_array[he_it], u, first_call);
      }
    }
    hypernode(v).disable();
    --_current_num_hypernodes;
    
    _n_hypergraph.contract(u,v);
    ASSERT(verifyNeighborhoodList(u),"Neighborhood list of HN " << u << " didn't match with real neighborhood!");
    ASSERT([&]() {
      for(HyperedgeID he : incidentEdges(u)) {
	for(HypernodeID pin : pins(he)) {
	  if(!verifyNeighborhoodList(pin)) {
	    return false;
	  }
	}
      }
      return true;
    } (), "Neighborhood list of HN"<<u<<" neighbors didn't match with the real neighborhood!");
    
    return Memento(u, u_offset, u_size, v);
  }

  template <typename GainChanges>
  void uncontract(const Memento& memento, GainChanges& changes, Int2Type<static_cast<int>(RefinementAlgorithm::twoway_fm)>) noexcept {
    ASSERT(!hypernode(memento.u).isDisabled(), "Hypernode " << memento.u << " is disabled");
    ASSERT(hypernode(memento.v).isDisabled(), "Hypernode " << memento.v << " is not invalid");
    ASSERT(changes.representative.size() == 1, V(changes.representative.size()));
    ASSERT(changes.contraction_partner.size() == 1, V(changes.contraction_partner.size()));

    HyperedgeWeight& changes_u = changes.representative[0];
    HyperedgeWeight& changes_v = changes.contraction_partner[0];

    DBG(dbg_hypergraph_uncontraction, "uncontracting (" << memento.u << "," << memento.v << ")");
    hypernode(memento.v).enable();
    ++_current_num_hypernodes;
    hypernode(memento.v).part_id = hypernode(memento.u).part_id;
    ++_part_info[partID(memento.u)].size;

    ASSERT(partID(memento.v) != kInvalidPartition,
           "PartitionID " << partID(memento.u) << " of representative HN " << memento.u <<
           " is INVALID - therefore wrong partition id was inferred for uncontracted HN "
           << memento.v);

    _hes_not_containing_u.resetAllBitsToFalse();
    // Assume all HEs did not contain u and we have to undo Case 2 operations.
    for (const HyperedgeID he : incidentEdges(memento.v)) {
      _hes_not_containing_u.setBit(he, true);
    }

    for (HyperedgeID i = memento.u_first_entry; i < memento.u_first_entry + memento.u_size; ++i) {
      const HyperedgeID he = _incidence_array[i];
      if (!_hes_not_containing_u[he]) {
        // These are hyperedges that are not connected to v after the uncontraction operation,
        // because they initially were only connected to u before the contraction.
        if (connectivity(he) > 1) {
          // because after uncontraction v is not connected to that HE anymore
          changes_v -= pinCountInPart(he, partID(memento.u)) == 1 ? edgeWeight(he) : 0;
        } else {
          // because after uncontraction v is not connected to that HE anymore
          ASSERT(pinCountInPart(he, partID(memento.u)) > 1, "Found Single-Node HE!");
          changes_v += edgeWeight(he);
        }
      }
      // Those HEs actually contained u and therefore will result in a Case 1 undo operation.
      _hes_not_containing_u.setBit(he, false);
    }

    if (hypernode(memento.u).size() - memento.u_size > 0) {
      // Undo case 2 opeations (i.e. Entry of pin v in HE e was reused to store connection to u):
      // Set incidence entry containing u for this HE e back to v, because this slot was used
      // to store the new edge to representative u during contraction as u was not a pin of e.
      for (const HyperedgeID he : incidentEdges(memento.u)) {
        if (_hes_not_containing_u[he]) {
          DBG(dbg_hypergraph_uncontraction, "resetting reused Pinslot of HE " << he << " from "
              << memento.u << " to " << memento.v);
          resetReusedPinSlotToOriginalValue(he, memento);

          if (connectivity(he) > 1) {
            --hypernode(memento.u).num_incident_cut_hes;    // because u is not connected to that cut HE anymore
            ++hypernode(memento.v).num_incident_cut_hes;    // because v is connected to that cut HE
            // because after uncontraction, u is not connected to that HE anymore
            changes_u -= pinCountInPart(he, partID(memento.u)) == 1 ? edgeWeight(he) : 0;
          } else {
            // because after uncontraction, u is not connected to that HE anymore
            ASSERT(pinCountInPart(he, partID(memento.u)) > 1, "Found Single-Node HE!");
            changes_u += edgeWeight(he);
          }
          // The state of this hyperedge now resembles the state before contraction.
          // Thus we don't need to process them any further.
        }
      }

      // Check if we can remove the dynamically added incidence entries for u:
      // If the old incidence entries are located before the current ones, than the
      // contraction we are currently undoing was responsible moving the incidence entries
      // of u to the end of the incidence array. Thus we can remove these entries now.
      // Otherwise these entries will still be needed by upcoming uncontract operations.
      if (memento.u_first_entry < hypernode(memento.u).firstEntry()) {
        _incidence_array.erase(_incidence_array.begin() + hypernode(memento.u).firstEntry(),
                               _incidence_array.end());
      }
    }

    restoreRepresentative(memento);

    // Undo case 1 operations (i.e. Pin v was just cut off by decreasing size of HE e):
    // Thus it is sufficient to just increase the size of the HE e to re-add the entry of v.
    for (const HyperedgeID he : incidentEdges(memento.v)) {
      if (!_hes_not_containing_u[he]) {
        DBG(dbg_hypergraph_uncontraction, "increasing size of HE " << he);
        ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
        hyperedge(he).increaseSize();
        increasePinCountInPart(he, partID(memento.v));
        ASSERT(_incidence_array[hyperedge(he).firstInvalidEntry() - 1] == memento.v,
               "Incorrect case 1 restore of HE " << he << ": "
               << _incidence_array[hyperedge(he).firstInvalidEntry() - 1] << "!=" << memento.v
               << "(while uncontracting: (" << memento.u << "," << memento.v << "))");

        if (connectivity(he) > 1) {
          ++hypernode(memento.v).num_incident_cut_hes;     // because v is connected to that cut HE
        }

        // Either the HE could have been removed from the cut before the move, or the HE
        // was not present before uncontraction, because it was a removed single-node HE
        // In both cases, there is no positive gain anymore, because the HE can't be removed
        // from the cut or a move now would result in a new cut edge.
        changes_u -= pinCountInPart(he, partID(memento.u)) == 2 ? edgeWeight(he) : 0;
        changes_v -= pinCountInPart(he, partID(memento.u)) == 2 ? edgeWeight(he) : 0;

        ++_current_num_pins;
      }
    }

    ASSERT(hypernode(memento.u).num_incident_cut_hes == numIncidentCutHEs(memento.u),
           V(memento.u) << V(hypernode(memento.u).num_incident_cut_hes) << V(numIncidentCutHEs(memento.u)));
    ASSERT(hypernode(memento.v).num_incident_cut_hes == numIncidentCutHEs(memento.v),
           V(memento.v) << V(hypernode(memento.v).num_incident_cut_hes) << V(numIncidentCutHEs(memento.v)));
  }


  template <typename GainChanges>
  void uncontract(const Memento& memento, GainChanges& changes, Int2Type<static_cast<int>(RefinementAlgorithm::kway_fm_km1)>) noexcept {
    LOG("yeah");
    exit(0);
  }

  void uncontract(const Memento& memento) noexcept {
    ASSERT(!hypernode(memento.u).isDisabled(), "Hypernode " << memento.u << " is disabled");
    ASSERT(hypernode(memento.v).isDisabled(), "Hypernode " << memento.v << " is not invalid");

    DBG(dbg_hypergraph_uncontraction, "uncontracting (" << memento.u << "," << memento.v << ")");
    hypernode(memento.v).enable();
    ++_current_num_hypernodes;
    hypernode(memento.v).part_id = hypernode(memento.u).part_id;
    ++_part_info[partID(memento.u)].size;

    ASSERT(partID(memento.v) != kInvalidPartition,
           "PartitionID " << partID(memento.u) << " of representative HN " << memento.u <<
           " is INVALID - therefore wrong partition id was inferred for uncontracted HN "
           << memento.v);

    _hes_not_containing_u.resetAllBitsToFalse();
    // Assume all HEs did not contain u and we have to undo Case 2 operations.
    for (const HyperedgeID he : incidentEdges(memento.v)) {
      _hes_not_containing_u.setBit(he, true);
    }

    for (HyperedgeID i = memento.u_first_entry; i < memento.u_first_entry + memento.u_size; ++i) {
      const HyperedgeID he = _incidence_array[i];
      // Those HEs actually contained u and therefore will result in a Case 1 undo operation.
      _hes_not_containing_u.setBit(he, false);
    }

    if (hypernode(memento.u).size() - memento.u_size > 0) {
      // Undo case 2 opeations (i.e. Entry of pin v in HE e was reused to store connection to u):
      // Set incidence entry containing u for this HE e back to v, because this slot was used
      // to store the new edge to representative u during contraction as u was not a pin of e.
      for (const HyperedgeID he : incidentEdges(memento.u)) {
        if (_hes_not_containing_u[he]) {
          DBG(dbg_hypergraph_uncontraction, "resetting reused Pinslot of HE " << he << " from "
              << memento.u << " to " << memento.v);
          resetReusedPinSlotToOriginalValue(he, memento);

          if (connectivity(he) > 1) {
            --hypernode(memento.u).num_incident_cut_hes;    // because u is not connected to that cut HE anymore
            ++hypernode(memento.v).num_incident_cut_hes;    // because v is connected to that cut HE
          }
          // The state of this hyperedge now resembles the state before contraction.
          // Thus we don't need to process them any further.
        }
      }

      // Check if we can remove the dynamically added incidence entries for u:
      // If the old incidence entries are located before the current ones, than the
      // contraction we are currently undoing was responsible moving the incidence entries
      // of u to the end of the incidence array. Thus we can remove these entries now.
      // Otherwise these entries will still be needed by upcoming uncontract operations.
      if (memento.u_first_entry < hypernode(memento.u).firstEntry()) {
        _incidence_array.erase(_incidence_array.begin() + hypernode(memento.u).firstEntry(),
                               _incidence_array.end());
      }
    }

    restoreRepresentative(memento);

    // Undo case 1 operations (i.e. Pin v was just cut off by decreasing size of HE e):
    // Thus it is sufficient to just increase the size of the HE e to re-add the entry of v.
    for (const HyperedgeID he : incidentEdges(memento.v)) {
      if (!_hes_not_containing_u[he]) {
        DBG(dbg_hypergraph_uncontraction, "increasing size of HE " << he);
        ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
        hyperedge(he).increaseSize();
        increasePinCountInPart(he, partID(memento.v));
        ASSERT(_incidence_array[hyperedge(he).firstInvalidEntry() - 1] == memento.v,
               "Incorrect case 1 restore of HE " << he << ": "
               << _incidence_array[hyperedge(he).firstInvalidEntry() - 1] << "!=" << memento.v
               << "(while uncontracting: (" << memento.u << "," << memento.v << "))");

        if (connectivity(he) > 1) {
          ++hypernode(memento.v).num_incident_cut_hes;     // because v is connected to that cut HE
        }

        ++_current_num_pins;
      }
    }

    ASSERT(hypernode(memento.u).num_incident_cut_hes == numIncidentCutHEs(memento.u),
           V(memento.u) << V(hypernode(memento.u).num_incident_cut_hes) << V(numIncidentCutHEs(memento.u)));
    ASSERT(hypernode(memento.v).num_incident_cut_hes == numIncidentCutHEs(memento.v),
           V(memento.v) << V(hypernode(memento.v).num_incident_cut_hes) << V(numIncidentCutHEs(memento.v)));
  }

  void changeNodePart(const HypernodeID hn, const PartitionID from, const PartitionID to) noexcept {
    Dummy dummy;  // smart compiler hopefully optimizes this away
    changeNodePart(hn, from, to, dummy);
  }

  // Special version for FM algorithm that identifies hypernodes that became non-border
  // hypernodes because of the move.
  template <typename Container>
  void changeNodePart(const HypernodeID hn, const PartitionID from, const PartitionID to,
                      Container& non_border_hns_to_remove) noexcept {
    ASSERT(!hypernode(hn).isDisabled(), "Hypernode " << hn << " is disabled");
    ASSERT(partID(hn) == from, "Hypernode" << hn << " is not in partition " << from);
    ASSERT(to < _k && to != kInvalidPartition, "Invalid to_part:" << to);
    ASSERT(from != to, "from part " << from << " == " << to << " part");
    updatePartInfo(hn, from, to);
    for (const HyperedgeID he : incidentEdges(hn)) {
      const bool no_pins_left_in_source_part = decreasePinCountInPart(he, from);
      const bool only_one_pin_in_to_part = increasePinCountInPart(he, to);

      if ((no_pins_left_in_source_part && !only_one_pin_in_to_part)) {
        if (pinCountInPart(he, to) == edgeSize(he)) {
          for (const HypernodeID pin : pins(he)) {
            --hypernode(pin).num_incident_cut_hes;
            if (hypernode(pin).num_incident_cut_hes == 0) {
              // ASSERT(std::find(non_border_hns_to_remove.cbegin(),
              //                  non_border_hns_to_remove.cend(), pin) ==
              //        non_border_hns_to_remove.end(),
              //        V(pin));
              non_border_hns_to_remove.push_back(pin);
            }
          }
        }
      } else if (!no_pins_left_in_source_part && only_one_pin_in_to_part) {
        if (pinCountInPart(he, from) == edgeSize(he) - 1) {
          for (const HypernodeID pin : pins(he)) {
            ++hypernode(pin).num_incident_cut_hes;
          }
        }
      }
      /**ASSERT([&]() -> bool {
         HypernodeID num_pins = 0;
         for (PartitionID i = 0; i < _k; ++i) {
         num_pins += pinCountInPart(he, i);
         }
         return num_pins == edgeSize(he);
         } (),
         "Incorrect calculation of pin counts");**/
    }
    // ASSERT([&]() {
    //    for (const HyperedgeID he : incidentEdges(hn)) {
    //    for (const HypernodeID pin : pins(he)) {
    //      if (pin == 1891) {
    //        LOG(V(hypernode(pin).num_incident_cut_hes));
    //      }

    //    if (hypernode(pin).num_incident_cut_hes != numIncidentCutHEs(pin)) {
    //    LOGVAR(pin);
    //    LOGVAR(hypernode(pin).num_incident_cut_hes);
    //    LOGVAR(numIncidentCutHEs(pin));
    //    return false;
    //    }
    //    }
    //    }
    //    return true;
    //    } (), "Inconsisten #CutHEs state");
  }

  bool isBorderNode(const HypernodeID hn) const {
    ASSERT(!hypernode(hn).isDisabled(), "Hypernode " << hn << " is disabled");
    ASSERT(hypernode(hn).num_incident_cut_hes == numIncidentCutHEs(hn), V(hn));
    ASSERT((hypernode(hn).num_incident_cut_hes > 0) == isBorderNodeInternal(hn), V(hn));
    return hypernode(hn).num_incident_cut_hes > 0;
  }


  // Used to initially set the partition ID of a HN after initial partitioning
  void setNodePart(const HypernodeID hn, const PartitionID id) noexcept {
    ASSERT(!hypernode(hn).isDisabled(), "Hypernode " << hn << " is disabled");
    ASSERT(partID(hn) == kInvalidPartition, "Hypernode" << hn << " is not unpartitioned: "
           << partID(hn));
    ASSERT(id < _k && id != kInvalidPartition, "Invalid part:" << id);
    updatePartInfo(hn, id);
    for (const HyperedgeID he : incidentEdges(hn)) {
      increasePinCountInPart(he, id);
    }
  }
  
  void createNeighborhoodHypergraph() {
    _n_hypergraph.resize(_num_hypernodes);
    
    //Creating Neighborhood Hypergraph
    SparseSet<int> neighbors(_num_hypernodes);

    for (const auto hn : nodes()) {
      for (const auto he : incidentEdges(hn)) {
        for (const auto pin : pins(he)) {
          neighbors.add(pin);
        }
      }
      std::vector<int> neighbor;
      neighbor.swap(neighbors.elements());
      std::sort(neighbor.begin(), neighbor.end());
      _n_hypergraph.setNeighborhoodOfHypernode(hn,neighbor);

      neighbors.clear();
    }
    ASSERT([&]() {
      for(HypernodeID hn : nodes()) {
	if(!verifyNeighborhoodList(hn)) {
	  return false;
	}
      }
      return true;
    } (),"Neighborhood hypergraph creation failed!");
  }

  // Deleting a hypernode amounts to removing the undirected internal edge between
  // the hypernode vertex and each of its incident hyperedge vertices as well as
  // disabling the hypernode vertex.
  // Since disabling the vertex ensures that it won't be visible/accessible from
  // the outside, we do NOT explicitely remove the directed internal edges
  // (HypernodeVertex,HyperedgeVertex). Instead we only remove the directed
  // internal edges (HyperedgeVertex,HypernodeVertex) to indicate that the hypernode
  // is no longer associated with the corresponding Hyperedge.
  // This _partial_ deletion of the internal incidence information allows us to
  // efficiently restore a removed Hypernode (currently not implemented):
  // After re-enabling the hypernode, we can directly access the information about the
  // hyperedges it was incident to (since we did not delete this information in the first
  // place): Thus it is possible to iterate over the incident hyperedges and just restore
  // the corresponding internal edge (HyperedgeVertex, HypernodeVertex) which was cut off
  // the hyperedge-vertex.
  // ATTENTION: In order for this implementation produce correct restore results, it is
  //            necessary that the restoreNode calls have to replay the removeNode calls
  //            in __reversed__ order.
  void removeNode(const HypernodeID u) noexcept  {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode is disabled!");
    for (const HyperedgeID e : incidentEdges(u)) {
      removeInternalEdge(e, u, _hyperedges);
      if (partID(u) != kInvalidPartition) {
        decreasePinCountInPart(e, partID(u));
      }
      --_current_num_pins;
    }
    hypernode(u).disable();
    --_current_num_hypernodes;
  }

  // Deleting a hyperedge amounts to removing the undirected internal edge between
  // the hypernode vertex and each of its incident hyperedge vertices as well as
  // disabling the hypernode vertex itself.
  // Since disabling the vertex ensures that it won't be visible/accessible from
  // the outside, we do NOT explicitely remove the directed internal edges
  // (HyperedgeVertex,HypernodeVertex). Instead we only remove the directed
  // internal edges (HypernodeVertex,HyperedgeVertex) to indicate that the hyperedge
  // is no longer associated with the corresponding hypernode.
  // This _partial_ deletion of the internal incidence information allows us to
  // efficiently restore a removed hyperedge (see restoreEdge(HyperedgeID he)).
  // The flag "disable_unconnected_hypernodes" can be used to differentiate between two
  // following two intentions / use-cases of removeEdge:
  // 1.) During coarsening, we want do remove the hyperedge, but leave hypernodes
  //     intact. This is used to remove any single-node hyperedges. The hypernode
  //     which was the only pin of this hyperedge has to stay in the graph, because
  //     it contains information about the graph structure (i.e. its weight represents)
  //     the number of hypernodes that have been contracted with it. Thus is achieved
  //     by setting the flag to _false_.
  // 2.) In order to avoid tedious reevaluation of ratings for really large hyperedges
  //     we want to provide an option do _really_ delete these edges from the graph before
  //     starting the actual n-level partitioning. In this case, we _do_ want unconnected
  //     hypernodes to disappear from the graph. After the partitioning is finished, we then
  //     reintegrate these edges into the graph. As we sacrificed theses edges in the beginning
  //     we are willing to pay the price that these edges now inevitably will become cut-edges.
  //     Setting the flag to _true_ removes any hypernode that is unconntected after the removal
  //     of the hyperedge.
  void removeEdge(const HyperedgeID he, const bool disable_unconnected_hypernodes) noexcept {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge is disabled!");
    for (const HypernodeID pin : pins(he)) {
      removeInternalEdge(pin, he, _hypernodes);
      disableHypernodeIfUnconnected(pin, disable_unconnected_hypernodes);
      --_current_num_pins;
    }
    hyperedge(he).disable();
    invalidatePartitionPinCounts(he);
    --_current_num_hyperedges;
  }

  // Restores the deleted Hyperedge. Since the hyperedge information is left intact on the
  // hyperedge vertex, we reuse this information to restore the information on the incident
  // hypernodes (i.e. pins). Since the removal of the internal edge (HN_Vertex,HE_Vertex)
  // was done by swapping the HyperedgeID to the end of the edgelist of the hypernode and
  // decrementing the size, it is necessary to perform the restore operations __in reverse__
  // order as the removal operations occurred!
  void restoreEdge(const HyperedgeID he) noexcept {
    ASSERT(hyperedge(he).isDisabled(), "Hyperedge is enabled!");
    enableEdge(he);
    resetPartitionPinCounts(he);
    for (const HypernodeID pin : pins(he)) {
      ASSERT(std::count(_incidence_array.begin() + hypernode(pin).firstEntry(),
                        _incidence_array.begin() + hypernode(pin).firstInvalidEntry(), he)
             == 0,
             "HN " << pin << " is already connected to HE " << he);
      DBG(dbg_hypergraph_restore_edge, "re-adding pin  " << pin << " to HE " << he);
      enableHypernodeIfPreviouslyUnconnected(pin);
      hypernode(pin).increaseSize();
      if (partID(pin) != kInvalidPartition) {
        increasePinCountInPart(he, partID(pin));
      }

      ASSERT(_incidence_array[hypernode(pin).firstInvalidEntry() - 1] == he,
             "Incorrect restore of HE " << he);
      ++_current_num_pins;
    }
  }

  // Used during uncontraction to restore parallel HEs and simultaneously set
  // the correct number of incident cut HEs for all pins.
  void restoreEdge(const HyperedgeID he, const HyperedgeID old_representative) noexcept {
    ASSERT(hyperedge(he).isDisabled(), "Hyperedge is enabled!");
    enableEdge(he);
    resetPartitionPinCounts(he);
    for (const HypernodeID pin : pins(he)) {
      ASSERT(std::count(_incidence_array.begin() + hypernode(pin).firstEntry(),
                        _incidence_array.begin() + hypernode(pin).firstInvalidEntry(), he)
             == 0,
             "HN " << pin << " is already connected to HE " << he);
      DBG(dbg_hypergraph_restore_edge, "re-adding pin  " << pin << " to HE " << he);
      enableHypernodeIfPreviouslyUnconnected(pin);
      hypernode(pin).increaseSize();
      if (partID(pin) != kInvalidPartition) {
        increasePinCountInPart(he, partID(pin));
      }

      if (connectivity(old_representative) > 1) {
        ++hypernode(pin).num_incident_cut_hes;
      }

      ASSERT(_incidence_array[hypernode(pin).firstInvalidEntry() - 1] == he,
             "Incorrect restore of HE " << he);
      ++_current_num_pins;
    }
  }

  void resetPartitioning() noexcept {
    for (HypernodeID i = 0; i < _num_hypernodes; ++i) {
      hypernode(i).part_id = kInvalidPartition;
    }
    std::fill(_part_info.begin(), _part_info.end(), PartInfo());
    std::fill(_pins_in_part.begin(), _pins_in_part.end(), 0);
    for (HyperedgeID i = 0; i < _num_hyperedges; ++i) {
      hyperedge(i).connectivity = 0;
      _connectivity_sets[i].clear();
    }
    for (HypernodeID i = 0; i < _num_hypernodes; ++i) {
      hypernode(i).num_incident_cut_hes = 0;
    }
  }

  Type type() const noexcept {
    if (isModified()) {
      return Type::EdgeAndNodeWeights;
    } else {
      return _type;
    }
  }

  std::string typeAsString() const noexcept {
    return typeToString(type());
  }

  // Accessors and mutators.
  HyperedgeID nodeDegree(const HypernodeID u) const noexcept {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode " << u << " is disabled");
    return hypernode(u).size();
  }

  HypernodeID edgeSize(const HyperedgeID e) const noexcept {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge " << e << " is disabled");
    return hyperedge(e).size();
  }

  HypernodeWeight nodeWeight(const HypernodeID u) const noexcept {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode " << u << " is disabled");
    return hypernode(u).weight();
  }

  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) noexcept {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode " << u << " is disabled");
    hypernode(u).setWeight(weight);
  }

  HyperedgeWeight edgeWeight(const HyperedgeID e) const noexcept {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge " << e << " is disabled");
    return hyperedge(e).weight();
  }

  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) noexcept {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge " << e << " is disabled");
    hyperedge(e).setWeight(weight);
  }

  PartitionID partID(const HypernodeID u) const noexcept {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode " << u << " is disabled");
    return hypernode(u).part_id;
  }

  bool nodeIsEnabled(const HypernodeID u) const noexcept {
    return !hypernode(u).isDisabled();
  }

  bool edgeIsEnabled(const HyperedgeID e) const noexcept {
    return !hyperedge(e).isDisabled();
  }

  HypernodeID initialNumNodes() const noexcept {
    return _num_hypernodes;
  }

  HyperedgeID initialNumEdges() const noexcept {
    return _num_hyperedges;
  }

  HypernodeID initialNumPins()  const noexcept {
    return _num_pins;
  }

  HypernodeID numNodes() const noexcept {
    ASSERT([&]() {
        HypernodeID count = 0;
        for (HypernodeID i = 0; i < _num_hypernodes; ++i) {
          if (!hypernode(i).isDisabled()) {
            ++count;
          }
        }
        return count;
      } () == _current_num_hypernodes,
           "Inconsistent Hypergraph State:" << "current_num_hypernodes=" << _current_num_hypernodes
           << "!= # enabled hypernodes=" <<[&]() {
        HypernodeID count = 0;
        for (HypernodeID i = 0; i < _num_hypernodes; ++i) {
          if (!hypernode(i).isDisabled()) {
            ++count;
          }
        }
        return count;
      } ());
    return _current_num_hypernodes;
  }

  HyperedgeID numEdges() const noexcept {
    ASSERT([&]() {
        HyperedgeID count = 0;
        for (HyperedgeID i = 0; i < _num_hyperedges; ++i) {
          if (!hyperedge(i).isDisabled()) {
            ++count;
          }
        }
        return count;
      } () == _current_num_hyperedges,
           "Inconsistent Hypergraph State:" << "current_num_hyperedges=" << _current_num_hyperedges
           << "!= # enabled hyperedges=" <<[&]() {
        HyperedgeID count = 0;
        for (HyperedgeID i = 0; i < _num_hyperedges; ++i) {
          if (!hyperedge(i).isDisabled()) {
            ++count;
          }
        }
        return count;
      } ());
    return _current_num_hyperedges;
  }

  HypernodeID numPins() const noexcept {
    return _current_num_pins;
  }

  PartitionID k() const noexcept {
    return _k;
  }

  bool active(const HypernodeID u) const {
    return hypernode(u).state == _threshold_active;
  }

  bool marked(const HypernodeID u) const {
    return hypernode(u).state == _threshold_marked;
  }

  void mark(const HypernodeID u) {
    ASSERT(hypernode(u).state == _threshold_active, V(u));
    hypernode(u).state = _threshold_marked;
  }

  void markRebalanced(const HypernodeID u) {
    hypernode(u).state = _threshold_marked;
  }

  void activate(const HypernodeID u) {
    ASSERT(hypernode(u).state < _threshold_active, V(u));
    hypernode(u).state = _threshold_active;
  }

  void deactivate(const HypernodeID u) {
    ASSERT(hypernode(u).state == _threshold_active, V(u));
    --hypernode(u).state;
  }

  void resetHypernodeState() {
    if (_threshold_marked == std::numeric_limits<std::uint32_t>::max()) {
      for (HypernodeID hn = 0; hn < _num_hypernodes; ++hn) {
        hypernode(hn).state = 0;
      }
      _threshold_active = -2;
      _threshold_marked = -1;
    }
    _threshold_active += 2;
    _threshold_marked += 2;
  }

  // Needs to be called after initial partitioning in order to provide
  // correct borderNode checks.
  void initializeNumCutHyperedges() {
    // std::fill(_num_incident_cut_hes.begin(), _num_incident_cut_hes.end(), 0);
    // We need to fill all fields with 0, otherwise V-cycles provide incorrect
    // results. Since hypernodes might be disabled, we bypass assertion in
    // hypernode(.) here and directly access _hypernodes.
    for (HypernodeID hn = 0; hn < _num_hypernodes; ++hn) {
      _hypernodes[hn].num_incident_cut_hes = 0;
    }
    for (const HyperedgeID he : edges()) {
      if (connectivity(he) > 1) {
        for (const HypernodeID pin : pins(he)) {
          ++hypernode(pin).num_incident_cut_hes;
        }
      }
    }
  }

  HypernodeWeight weightOfHeaviestNode() const noexcept {
    HypernodeWeight max_weight = std::numeric_limits<HypernodeWeight>::min();
    for (const HypernodeID hn : nodes()) {
      max_weight = std::max(nodeWeight(hn), max_weight);
    }
    return max_weight;
  }

  HypernodeID pinCountInPart(const HyperedgeID he, const PartitionID id) const noexcept {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
    ASSERT(id < _k && id != kInvalidPartition, "Partition ID " << id << " is out of bounds");
    ASSERT(_pins_in_part[he * _k + id] != kInvalidCount, V(he) << V(id));
    return _pins_in_part[he * _k + id];
  }

  PartitionID connectivity(const HyperedgeID he) const noexcept {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
    return hyperedge(he).connectivity;
  }

  const std::vector<PartInfo> & partInfos() const noexcept {
    return _part_info;
  }

  HypernodeWeight totalWeight() const noexcept {
    return _total_weight;
  }

  HypernodeWeight partWeight(const PartitionID id) const noexcept {
    ASSERT(id < _k && id != kInvalidPartition, "Partition ID " << id << " is out of bounds");
    return _part_info[id].weight;
  }

  HypernodeID partSize(const PartitionID id) const noexcept {
    ASSERT(id < _k && id != kInvalidPartition, "Partition ID " << id << " is out of bounds");
    return _part_info[id].size;
  }

  HypernodeData & hypernodeData(const HypernodeID hn) noexcept {
    ASSERT(!hypernode(hn).isDisabled(), "Hypernode " << hn << " is disabled");
    return hypernode(hn);
  }

  HyperedgeData & hyperedgeData(const HyperedgeID he) noexcept {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
    return hyperedge(he);
  }

 private:
  FRIEND_TEST(AHypergraph, DisconnectsHypernodeFromHyperedge);
  FRIEND_TEST(AHypergraph, RemovesHyperedges);
  FRIEND_TEST(AHypergraph, DecrementsHypernodeDegreeOfAffectedHypernodesOnHyperedgeRemoval);
  FRIEND_TEST(AnIncidenceIterator, AllowsIterationOverIncidentHyperedges);
  FRIEND_TEST(AnIncidenceIterator, AllowsIterationOverPinsOfHyperedge);
  FRIEND_TEST(AHypergraphMacro, IteratesOverAllIncidentHyperedges);
  FRIEND_TEST(AHypergraphMacro, IteratesOverAllPinsOfAHyperedge);
  FRIEND_TEST(AContractionMemento, StoresOldStateOfInvolvedHypernodes);
  FRIEND_TEST(AnUncontractionOperation, DeletesIncidenceInfoAddedDuringContraction);
  FRIEND_TEST(AHypergraph, InvalidatesPartitionPinCountsOnHyperedgeRemoval);
  FRIEND_TEST(AHypergraph, CalculatesPinCountsOfAHyperedge);
  FRIEND_TEST(APartitionedHypergraph, StoresPartitionPinCountsForHyperedges);
  FRIEND_TEST(AHypergraph, ExtractedFromAPartitionedHypergraphHasInitializedPartitionInformation);

  GenericHypergraph() noexcept :
    _num_hypernodes(0),
    _num_hyperedges(0),
    _num_pins(0),
    _total_weight(0),
    _k(2),
    _type(Type::Unweighted),
    _current_num_hypernodes(0),
    _current_num_hyperedges(0),
    _current_num_pins(0),
    _threshold_active(1),
    _threshold_marked(2),
    _hypernodes(),
    _hyperedges(),
    _incidence_array(),
    _part_info(_k),
    _pins_in_part(),
    _connectivity_sets(),
    _hes_not_containing_u() { }

  std::pair<HyperedgeWeight, HyperedgeWeight> uncontract_impl(const Memento& memento, std::true_type) noexcept {
    ASSERT(!hypernode(memento.u).isDisabled(), "Hypernode " << memento.u << " is disabled");
    ASSERT(hypernode(memento.v).isDisabled(), "Hypernode " << memento.v << " is not invalid");

    DBG(dbg_hypergraph_uncontraction, "uncontracting (" << memento.u << "," << memento.v << ")");
    hypernode(memento.v).enable();
    ++_current_num_hypernodes;
    hypernode(memento.v).part_id = hypernode(memento.u).part_id;
    ++_part_info[partID(memento.u)].size;

    ASSERT(partID(memento.v) != kInvalidPartition,
           "PartitionID " << partID(memento.u) << " of representative HN " << memento.u <<
           " is INVALID - therefore wrong partition id was inferred for uncontracted HN "
           << memento.v);

    _hes_not_containing_u.resetAllBitsToFalse();
    // Assume all HEs did not contain u and we have to undo Case 2 operations.
    for (const HyperedgeID he : incidentEdges(memento.v)) {
      _hes_not_containing_u.setBit(he, true);
    }

    std::pair<HyperedgeWeight, HyperedgeWeight> ret = { 0, 0 };
    for (HyperedgeID i = memento.u_first_entry; i < memento.u_first_entry + memento.u_size; ++i) {
      const HyperedgeID he = _incidence_array[i];
      if (!_hes_not_containing_u[he]) {
        // These are hyperedges that are not connected to v after the uncontraction operation,
        // because they initially were only connected to u before the contraction.
        if (connectivity(he) > 1) {
          // because after uncontraction v is not connected to that HE anymore
          ret.second -= pinCountInPart(he, partID(memento.u)) == 1 ? edgeWeight(he) : 0;
        } else {
          // because after uncontraction v is not connected to that HE anymore
          ASSERT(pinCountInPart(he, partID(memento.u)) > 1, "Found Single-Node HE!");
          ret.second += edgeWeight(he);
        }
      }
      // Those HEs actually contained u and therefore will result in a Case 1 undo operation.
      _hes_not_containing_u.setBit(he, false);
    }

    if (hypernode(memento.u).size() - memento.u_size > 0) {
      // Undo case 2 opeations (i.e. Entry of pin v in HE e was reused to store connection to u):
      // Set incidence entry containing u for this HE e back to v, because this slot was used
      // to store the new edge to representative u during contraction as u was not a pin of e.
      for (const HyperedgeID he : incidentEdges(memento.u)) {
        if (_hes_not_containing_u[he]) {
          DBG(dbg_hypergraph_uncontraction, "resetting reused Pinslot of HE " << he << " from "
              << memento.u << " to " << memento.v);
          resetReusedPinSlotToOriginalValue(he, memento);

          if (connectivity(he) > 1) {
            --hypernode(memento.u).num_incident_cut_hes;    // because u is not connected to that cut HE anymore
            ++hypernode(memento.v).num_incident_cut_hes;    // because v is connected to that cut HE
            // because after uncontraction, u is not connected to that HE anymore
            ret.first -= pinCountInPart(he, partID(memento.u)) == 1 ? edgeWeight(he) : 0;
          } else {
            // because after uncontraction, u is not connected to that HE anymore
            ASSERT(pinCountInPart(he, partID(memento.u)) > 1, "Found Single-Node HE!");
            ret.first += edgeWeight(he);
          }
          // The state of this hyperedge now resembles the state before contraction.
          // Thus we don't need to process them any further.
        }
      }

      // Check if we can remove the dynamically added incidence entries for u:
      // If the old incidence entries are located before the current ones, than the
      // contraction we are currently undoing was responsible moving the incidence entries
      // of u to the end of the incidence array. Thus we can remove these entries now.
      // Otherwise these entries will still be needed by upcoming uncontract operations.
      if (memento.u_first_entry < hypernode(memento.u).firstEntry()) {
        _incidence_array.erase(_incidence_array.begin() + hypernode(memento.u).firstEntry(),
                               _incidence_array.end());
      }
    }

    restoreRepresentative(memento);

    // Undo case 1 operations (i.e. Pin v was just cut off by decreasing size of HE e):
    // Thus it is sufficient to just increase the size of the HE e to re-add the entry of v.
    for (const HyperedgeID he : incidentEdges(memento.v)) {
      if (!_hes_not_containing_u[he]) {
        DBG(dbg_hypergraph_uncontraction, "increasing size of HE " << he);
        ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
        hyperedge(he).increaseSize();
        increasePinCountInPart(he, partID(memento.v));
        ASSERT(_incidence_array[hyperedge(he).firstInvalidEntry() - 1] == memento.v,
               "Incorrect case 1 restore of HE " << he << ": "
               << _incidence_array[hyperedge(he).firstInvalidEntry() - 1] << "!=" << memento.v
               << "(while uncontracting: (" << memento.u << "," << memento.v << "))");

        if (connectivity(he) > 1) {
          ++hypernode(memento.v).num_incident_cut_hes;     // because v is connected to that cut HE
        }

        // Either the HE could have been removed from the cut before the move, or the HE
        // was not present before uncontraction, because it was a removed single-node HE
        // In both cases, there is no positive gain anymore, because the HE can't be removed
        // from the cut or a move now would result in a new cut edge.
        ret.first -= pinCountInPart(he, partID(memento.u)) == 2 ? edgeWeight(he) : 0;
        ret.second -= pinCountInPart(he, partID(memento.u)) == 2 ? edgeWeight(he) : 0;

        ++_current_num_pins;
      }
    }

    ASSERT(hypernode(memento.u).num_incident_cut_hes == numIncidentCutHEs(memento.u),
           V(memento.u) << V(hypernode(memento.u).num_incident_cut_hes) << V(numIncidentCutHEs(memento.u)));
    ASSERT(hypernode(memento.v).num_incident_cut_hes == numIncidentCutHEs(memento.v),
           V(memento.v) << V(hypernode(memento.v).num_incident_cut_hes) << V(numIncidentCutHEs(memento.v)));

    
    _n_hypergraph.uncontract(memento.u,memento.v);
    ASSERT(verifyNeighborhoodList(memento.u),"Neighborhood list of HN " << memento.u << " didn't match with real neighborhood!");
    ASSERT(verifyNeighborhoodList(memento.v),"Neighborhood list of HN " << memento.v << " didn't match with real neighborhood!");
    ASSERT([&]() {
      for(HyperedgeID he : incidentEdges(memento.u)) {
	for(HypernodeID pin : pins(he)) {
	  if(!verifyNeighborhoodList(pin)) {
	    return false;
	  }
	}
      }
      return true;
    } (), "Neighborhood list of HN"<<memento.u<<" neighbors didn't match with the real neighborhood!");
    ASSERT([&]() {
      for(HyperedgeID he : incidentEdges(memento.v)) {
	for(HypernodeID pin : pins(he)) {
	  if(!verifyNeighborhoodList(pin)) {
	    return false;
	  }
	}
      }
      return true;
    } (), "Neighborhood list of HN"<<memento.v<<" neighbors didn't match with the real neighborhood!");
    return ret;
  }

  void uncontract_impl(const Memento& memento, std::false_type) noexcept {
    ASSERT(!hypernode(memento.u).isDisabled(), "Hypernode " << memento.u << " is disabled");
    ASSERT(hypernode(memento.v).isDisabled(), "Hypernode " << memento.v << " is not invalid");

    DBG(dbg_hypergraph_uncontraction, "uncontracting (" << memento.u << "," << memento.v << ")");
    hypernode(memento.v).enable();
    ++_current_num_hypernodes;
    hypernode(memento.v).part_id = hypernode(memento.u).part_id;
    ++_part_info[partID(memento.u)].size;

    ASSERT(partID(memento.v) != kInvalidPartition,
           "PartitionID " << partID(memento.u) << " of representative HN " << memento.u <<
           " is INVALID - therefore wrong partition id was inferred for uncontracted HN "
           << memento.v);

    _hes_not_containing_u.resetAllBitsToFalse();
    // Assume all HEs did not contain u and we have to undo Case 2 operations.
    for (const HyperedgeID he : incidentEdges(memento.v)) {
      _hes_not_containing_u.setBit(he, true);
    }

    for (HyperedgeID i = memento.u_first_entry; i < memento.u_first_entry + memento.u_size; ++i) {
      const HyperedgeID he = _incidence_array[i];
      // Those HEs actually contained u and therefore will result in a Case 1 undo operation.
      _hes_not_containing_u.setBit(he, false);
    }

    if (hypernode(memento.u).size() - memento.u_size > 0) {
      // Undo case 2 opeations (i.e. Entry of pin v in HE e was reused to store connection to u):
      // Set incidence entry containing u for this HE e back to v, because this slot was used
      // to store the new edge to representative u during contraction as u was not a pin of e.
      for (const HyperedgeID he : incidentEdges(memento.u)) {
        if (_hes_not_containing_u[he]) {
          DBG(dbg_hypergraph_uncontraction, "resetting reused Pinslot of HE " << he << " from "
              << memento.u << " to " << memento.v);
          resetReusedPinSlotToOriginalValue(he, memento);

          if (connectivity(he) > 1) {
            --hypernode(memento.u).num_incident_cut_hes;    // because u is not connected to that cut HE anymore
            ++hypernode(memento.v).num_incident_cut_hes;    // because v is connected to that cut HE
          }
          // The state of this hyperedge now resembles the state before contraction.
          // Thus we don't need to process them any further.
        }
      }

      // Check if we can remove the dynamically added incidence entries for u:
      // If the old incidence entries are located before the current ones, than the
      // contraction we are currently undoing was responsible moving the incidence entries
      // of u to the end of the incidence array. Thus we can remove these entries now.
      // Otherwise these entries will still be needed by upcoming uncontract operations.
      if (memento.u_first_entry < hypernode(memento.u).firstEntry()) {
        _incidence_array.erase(_incidence_array.begin() + hypernode(memento.u).firstEntry(),
                               _incidence_array.end());
      }
    }

    restoreRepresentative(memento);

    // Undo case 1 operations (i.e. Pin v was just cut off by decreasing size of HE e):
    // Thus it is sufficient to just increase the size of the HE e to re-add the entry of v.
    for (const HyperedgeID he : incidentEdges(memento.v)) {
      if (!_hes_not_containing_u[he]) {
        DBG(dbg_hypergraph_uncontraction, "increasing size of HE " << he);
        ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
        hyperedge(he).increaseSize();
        increasePinCountInPart(he, partID(memento.v));
        ASSERT(_incidence_array[hyperedge(he).firstInvalidEntry() - 1] == memento.v,
               "Incorrect case 1 restore of HE " << he << ": "
               << _incidence_array[hyperedge(he).firstInvalidEntry() - 1] << "!=" << memento.v
               << "(while uncontracting: (" << memento.u << "," << memento.v << "))");

        if (connectivity(he) > 1) {
          ++hypernode(memento.v).num_incident_cut_hes;     // because v is connected to that cut HE
        }

        ++_current_num_pins;
      }
    }

    ASSERT(hypernode(memento.u).num_incident_cut_hes == numIncidentCutHEs(memento.u),
           V(memento.u) << V(hypernode(memento.u).num_incident_cut_hes) << V(numIncidentCutHEs(memento.u)));
    ASSERT(hypernode(memento.v).num_incident_cut_hes == numIncidentCutHEs(memento.v),
           V(memento.v) << V(hypernode(memento.v).num_incident_cut_hes) << V(numIncidentCutHEs(memento.v)));
    
    _n_hypergraph.uncontract(memento.u,memento.v);
    ASSERT(verifyNeighborhoodList(memento.u),"Neighborhood list of HN " << memento.u << " didn't match with real neighborhood!");
    ASSERT(verifyNeighborhoodList(memento.v),"Neighborhood list of HN " << memento.v << " didn't match with real neighborhood!");
    ASSERT([&]() {
      for(HyperedgeID he : incidentEdges(memento.u)) {
	for(HypernodeID pin : pins(he)) {
	  if(!verifyNeighborhoodList(pin)) {
	    return false;
	  }
	}
      }
      return true;
    } (), "Neighborhood list of HN"<<memento.u<<" neighbors didn't match with the real neighborhood!");
    ASSERT([&]() {
      for(HyperedgeID he : incidentEdges(memento.v)) {
	for(HypernodeID pin : pins(he)) {
	  if(!verifyNeighborhoodList(pin)) {
	    return false;
	  }
	}
      }
      return true;
    } (), "Neighborhood list of HN"<<memento.v<<" neighbors didn't match with the real neighborhood!");
  }

  bool isBorderNodeInternal(const HypernodeID hn) const {
    ASSERT(!hypernode(hn).isDisabled(), "Hypernode " << hn << " is disabled");
    for (const HyperedgeID he : incidentEdges(hn)) {
      if (connectivity(he) > 1) {
        return true;
      }
    }
    return false;
  }


  HyperedgeID numIncidentCutHEs(const HypernodeID hn) const {
    HyperedgeID num_cut_hes = 0;
    for (const HyperedgeID he : incidentEdges(hn)) {
      if (connectivity(he) > 1) {
        ++num_cut_hes;
      }
    }
    return num_cut_hes;
  }

  void updatePartInfo(const HypernodeID u, const PartitionID id) noexcept {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode " << u << " is disabled");
    ASSERT(id < _k && id != kInvalidPartition, "Part ID" << id << " out of bounds!");
    ASSERT(hypernode(u).part_id == kInvalidPartition, "HN " << u << " is already assigned to part " << id);
    hypernode(u).part_id = id;
    _part_info[id].weight += nodeWeight(u);
    ++_part_info[id].size;
  }

  void updatePartInfo(const HypernodeID u, const PartitionID from, const PartitionID to) noexcept {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode " << u << " is disabled");
    ASSERT(from < _k && from != kInvalidPartition, "Part ID" << from << " out of bounds!");
    ASSERT(to < _k && to != kInvalidPartition, "Part ID" << to << " out of bounds!");
    ASSERT(hypernode(u).part_id == from, "HN " << u << " is not in part " << from);
    hypernode(u).part_id = to;
    _part_info[from].weight -= nodeWeight(u);
    --_part_info[from].size;
    _part_info[to].weight += nodeWeight(u);
    ++_part_info[to].size;
  }

  bool decreasePinCountInPart(const HyperedgeID he, const PartitionID id) noexcept {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
    ASSERT(pinCountInPart(he, id) > 0,
           "HE " << he << "does not have any pins in partition " << id);
    ASSERT(id < _k && id != kInvalidPartition, "Part ID" << id << " out of bounds!");
    ASSERT(_pins_in_part[he * _k + id] > 0, "invalid decrease");
    const size_t offset = he * _k + id;
    _pins_in_part[offset] -= 1;
    const bool connectivity_decreased = _pins_in_part[offset] == 0;
    if (connectivity_decreased) {
      _connectivity_sets[he].remove(id);
      hyperedge(he).connectivity -= 1;
    }
    return connectivity_decreased;
  }

  bool increasePinCountInPart(const HyperedgeID he, const PartitionID id) noexcept {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
    ASSERT(pinCountInPart(he, id) <= edgeSize(he),
           "HE " << he << ": pin_count[" << id << "]=" << pinCountInPart(he, id)
           << "edgesize=" << edgeSize(he));
    ASSERT(id < _k && id != kInvalidPartition, "Part ID" << id << " out of bounds!");
    const size_t offset = he * _k + id;
    _pins_in_part[offset] += 1;
    const bool connectivity_increased = _pins_in_part[offset] == 1;
    if (connectivity_increased) {
      hyperedge(he).connectivity += 1;
      _connectivity_sets[he].add(id);
    }
    return connectivity_increased;
  }

  void invalidatePartitionPinCounts(const HyperedgeID he) noexcept {
    ASSERT(hyperedge(he).isDisabled(),
           "Invalidation of pin counts only allowed for disabled hyperedges");
    for (PartitionID part = 0; part < _k; ++part) {
      _pins_in_part[he * _k + part] = kInvalidCount;
    }
    hyperedge(he).connectivity = 0;
    _connectivity_sets[he].clear();
  }

  void resetPartitionPinCounts(const HyperedgeID he) noexcept {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
    for (PartitionID part = 0; part < _k; ++part) {
      _pins_in_part[he * _k + part] = 0;
    }
  }

  bool isModified() const noexcept {
    return _current_num_pins != _num_pins || _current_num_hypernodes != _num_hypernodes ||
           _current_num_hyperedges != _num_hyperedges;
  }

  void enableEdge(const HyperedgeID e) noexcept {
    ASSERT(hyperedge(e).isDisabled(), "HE " << e << " is already enabled!");
    hyperedge(e).enable();
    ++_current_num_hyperedges;
  }

  void restoreRepresentative(const Memento& memento) noexcept {
    ASSERT(!hypernode(memento.u).isDisabled(), "Hypernode " << memento.u << " is disabled");
    hypernode(memento.u).setFirstEntry(memento.u_first_entry);
    hypernode(memento.u).setSize(memento.u_size);
    hypernode(memento.u).setWeight(hypernode(memento.u).weight() - hypernode(memento.v).weight());
  }


  void resetReusedPinSlotToOriginalValue(const HyperedgeID he, const Memento& memento) noexcept {
    ASSERT(!hyperedge(he).isDisabled(), "Hyperedge " << he << " is disabled");
    PinHandleIterator pin_begin, pin_end;
    std::tie(pin_begin, pin_end) = pinHandles(he);
    ASSERT(pin_begin != pin_end, "Accessed empty hyperedge");
    --pin_end;
    while (*pin_end != memento.u) {
      ASSERT(pin_end != pin_begin, "Pin " << memento.u << " not found in pinlist of HE " << he);
      --pin_end;
    }
    ASSERT(*pin_end == memento.u && std::distance(_incidence_array.begin(), pin_begin)
           <= std::distance(_incidence_array.begin(), pin_end), "Reused slot not found");
    *pin_end = memento.v;
  }

  void connectHyperedgeToRepresentative(const HyperedgeID e, const HypernodeID u,
                                        bool& first_call) noexcept {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode " << u << " is disabled");
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge " << e << " is disabled");
    ASSERT(partID(_incidence_array[hyperedge(e).firstInvalidEntry() - 1]) == partID(u),
           "Contraction target " << _incidence_array[hyperedge(e).firstInvalidEntry() - 1]
           << "& representative " << u << "are in different parts");
    // Hyperedge e does not contain u. Therefore we use the entry of v (i.e. the last entry
    // -- this is ensured by the contract method) in e's edge array to store the information
    // that u is now connected to e and add the edge (u,e) to indicate this conection also from
    // the hypernode's point of view.
    _incidence_array[hyperedge(e).firstInvalidEntry() - 1] = u;

    HypernodeVertex& nodeU = hypernode(u);
    if (first_call) {
      _incidence_array.insert(_incidence_array.cend(), _incidence_array.cbegin() + nodeU.firstEntry(),
                              _incidence_array.cbegin() + nodeU.firstInvalidEntry());
      nodeU.setFirstEntry(_incidence_array.size() - nodeU.size());
      first_call = false;
    }
    ASSERT(nodeU.firstInvalidEntry() == _incidence_array.size(), "Incidence Info of HN " << u
           << " is not at the end of the incidence array");
    _incidence_array.push_back(e);
    nodeU.increaseSize();
  }

  void disableHypernodeIfUnconnected(const HypernodeID hn,
                                     const bool disable_unconnected_hypernode) noexcept {
    if (disable_unconnected_hypernode && hypernode(hn).size() == 0) {
      hypernode(hn).disable();
      --_current_num_hypernodes;
    }
  }

  void enableHypernodeIfPreviouslyUnconnected(const HypernodeID hn) noexcept {
    if (hypernode(hn).isDisabled()) {
      hypernode(hn).enable();
      ++_current_num_hypernodes;
    }
  }

  template <typename Handle1, typename Handle2, typename Container>
  void removeInternalEdge(const Handle1 u, const Handle2 v, Container& container) noexcept {
    using std::swap;
    typename Container::reference vertex = container[u];
    ASSERT(!vertex.isDisabled(), "InternalVertex " << u << " is disabled");

    auto begin = _incidence_array.begin() + vertex.firstEntry();
    ASSERT(vertex.size() > 0, "InternalVertex is empty!");
    auto last_entry = _incidence_array.begin() + vertex.firstInvalidEntry() - 1;
    while (*begin != v) {
      ++begin;
    }
    ASSERT(begin < _incidence_array.begin() + vertex.firstInvalidEntry(), "Iterator out of bounds");
    swap(*begin, *last_entry);
    vertex.decreaseSize();
  }

  std::string typeToString(Type hypergraph_type) const noexcept {
    std::string typestring;
    switch (hypergraph_type) {
      case Type::Unweighted:
        typestring = "edgeWeights=false nodeWeights=false";
        break;
      case Type::EdgeWeights:
        typestring = "edgeWeights=true nodeWeights=false";
        break;
      case Type::NodeWeights:
        typestring = "edgeWeights=false nodeWeights=true";
        break;
      case Type::EdgeAndNodeWeights:
        typestring = "edgeWeights=true nodeWeights=true";
        break;
    }
    return typestring;
  }

  // Accessor for handles of hypernodes contained in hyperedge (aka pins)
  std::pair<PinHandleIterator, PinHandleIterator> pinHandles(const HyperedgeID he) noexcept {
    return std::make_pair(_incidence_array.begin() + hyperedge(he).firstEntry(),
                          _incidence_array.begin() + hyperedge(he).firstInvalidEntry());
  }

  // Accessor for hypernode-related information
  const HypernodeVertex & hypernode(const HypernodeID u) const noexcept {
    ASSERT(u < _num_hypernodes, "Hypernode " << u << " does not exist");
    return _hypernodes[u];
  }

  // Accessor for hyperedge-related information
  const HyperedgeVertex & hyperedge(const HyperedgeID e) const noexcept {
    ASSERT(e < _num_hyperedges, "Hyperedge " << e << " does not exist");
    return _hyperedges[e];
  }

  // To avoid code duplication we implement non-const version in terms of const version
  HypernodeVertex & hypernode(const HypernodeID u) noexcept {
    return const_cast<HypernodeVertex&>(static_cast<const GenericHypergraph&>(*this).hypernode(u));
  }

  HyperedgeVertex & hyperedge(const HyperedgeID e) noexcept {
    return const_cast<HyperedgeVertex&>(static_cast<const GenericHypergraph&>(*this).hyperedge(e));
  }
  
  bool verifyNeighborhoodList(HypernodeID x) {
	std::vector<bool> match(_n_hypergraph._neighbors[x].size(),false);
	for(HyperedgeID he : incidentEdges(x)) {
	  for(HypernodeID pin : pins(he)) {
	    bool found = false;
	    for(int k = 0; k < _n_hypergraph._neighbors[x].size(); ++k) {
	      if(_n_hypergraph._neighbors[x][k] == pin) {
		match[k] = true; found = true;
		break;
	      }
	    }
	    //Hypernode pin should be contained in the neighborhood list of hypernode x.
	    if(!found) {
	      LOG("HN " << pin << " should be in the neighborhood list of HN " << x << "!");
	      _n_hypergraph.printNeighborhood(x);
	      LOG("Real neighborhoodlist:");
	      for(HyperedgeID he : incidentEdges(x)) {
		std::cout << "HE " << he << ": ";
		for(HypernodeID pin : pins(he)) {
		  std::cout << pin << " ";
		}
		std::cout << std::endl;
	      }
	      return false;
	    }
	  }
	}
	for(int k = 0; k < _n_hypergraph._neighbors[x].size(); ++k) {
	  //Hypernode _neighbors[x][k] isn't in the neighborhood of hypernode x.
	  if(!match[k]) {
	    LOG("HN " << _n_hypergraph._neighbors[x][k] << " isn't in the neighborhood of HN " << x << "!");
	    _n_hypergraph.printNeighborhood(x);
	    return false;
	  }
	}
        return true;   
  }

  HypernodeID _num_hypernodes;
  HyperedgeID _num_hyperedges;
  HypernodeID _num_pins;
  HypernodeWeight _total_weight;
  int _k;
  Type _type;

  HypernodeID _current_num_hypernodes;
  HyperedgeID _current_num_hyperedges;
  HypernodeID _current_num_pins;

  std::uint32_t _threshold_active;
  std::uint32_t _threshold_marked;

  std::vector<HypernodeVertex> _hypernodes;
  std::vector<HyperedgeVertex> _hyperedges;
  std::vector<VertexID> _incidence_array;

  std::vector<PartInfo> _part_info;
  // for each hyperedge we store the connectivity set,
  // i.e. the parts it connects and the number of pins in that part
  std::vector<HypernodeID> _pins_in_part;
  ConnectivitySets<PartitionID, HyperedgeID> _connectivity_sets;

  // Used during uncontraction to decide how to perform the uncontraction operation.
  // Incident HEs of the representative either already contained u and v before the contraction
  // or only contained v. In the latter case, we use _hes_not_containing_u[he]=true, to
  // indicate that he have to undo a "Case 2" Operation, i.e. one, where the pin slot of
  // v was re-used during contraction.
  FastResetBitVector<> _hes_not_containing_u;
  
  NeighborhoodHypergraph _n_hypergraph;

  template <typename Hypergraph>
  friend std::pair<std::unique_ptr<Hypergraph>,
                   std::vector<typename Hypergraph::HypernodeID> > extractPartAsUnpartitionedHypergraphForBisection(const Hypergraph& hypergraph,
                                                                                                                    const typename Hypergraph::PartitionID part,
                                                                                                                    const bool split_nets = false);

  template <typename Hypergraph>
  friend bool verifyEquivalenceWithoutPartitionInfo(const Hypergraph& expected,
                                                    const Hypergraph& actual);

  template <typename Hypergraph>
  friend bool verifyEquivalenceWithPartitionInfo(const Hypergraph& expected,
                                                 const Hypergraph& actual);

  template <typename Hypergraph>
  friend std::pair<std::unique_ptr<Hypergraph>,
                   std::vector<typename Hypergraph::HypernodeID> > reindex(const Hypergraph& hypergraph);
};

template <typename Hypergraph>
bool verifyEquivalenceWithoutPartitionInfo(const Hypergraph& expected, const Hypergraph& actual) {
  ASSERT(expected._num_hypernodes == actual._num_hypernodes,
         V(expected._num_hypernodes) << V(actual._num_hypernodes));
  ASSERT(expected._num_hyperedges == actual._num_hyperedges,
         V(expected._num_hyperedges) << V(actual._num_hyperedges));
  ASSERT(expected._num_pins == actual._num_pins, V(expected._num_pins) << V(actual._num_pins));
  ASSERT(expected._total_weight == actual._total_weight,
         V(expected._total_weight) << V(actual._total_weight));
  ASSERT(expected._k == actual._k, V(expected._k) << V(actual._k));
  ASSERT(expected._type == actual._type, "Error!");
  ASSERT(expected._current_num_hypernodes == actual._current_num_hypernodes,
         V(expected._current_num_hypernodes) << V(actual._current_num_hypernodes));
  ASSERT(expected._current_num_hyperedges == actual._current_num_hyperedges,
         V(expected._current_num_hyperedges) << V(actual._current_num_hyperedges));
  ASSERT(expected._current_num_pins == actual._current_num_pins,
         V(expected._current_num_pins) << V(actual._current_num_pins));
  ASSERT(expected._hypernodes == actual._hypernodes, "Error!");
  ASSERT(expected._hyperedges == actual._hyperedges, "Error!");

  std::vector<unsigned int> expected_incidence_array(expected._incidence_array);
  std::vector<unsigned int> actual_incidence_array(actual._incidence_array);
  std::sort(expected_incidence_array.begin(), expected_incidence_array.end());
  std::sort(actual_incidence_array.begin(), actual_incidence_array.end());

  ASSERT(expected_incidence_array == actual_incidence_array,
         "expected._incidence_array != actual._incidence_array");

  return expected._num_hypernodes == actual._num_hypernodes &&
         expected._num_hyperedges == actual._num_hyperedges &&
         expected._num_pins == actual._num_pins &&
         expected._k == actual._k &&
         expected._type == actual._type &&
         expected._current_num_hypernodes == actual._current_num_hypernodes &&
         expected._current_num_hyperedges == actual._current_num_hyperedges &&
         expected._current_num_pins == actual._current_num_pins &&
         expected._hypernodes == actual._hypernodes &&
         expected._hyperedges == actual._hyperedges &&
         expected_incidence_array == actual_incidence_array;
}

template <typename Hypergraph>
bool verifyEquivalenceWithPartitionInfo(const Hypergraph& expected, const Hypergraph& actual) {
  using HyperedgeID = typename Hypergraph::HyperedgeID;
  using HypernodeID = typename Hypergraph::HypernodeID;

  ASSERT(expected._part_info == actual._part_info, "Error");
  ASSERT(expected._pins_in_part == actual._pins_in_part, "Error");

  bool connectivity_sets_valid = true;
  for (const HyperedgeID he : actual.edges()) {
    ASSERT(expected.hyperedge(he).connectivity == actual.hyperedge(he).connectivity, V(he));
    if (expected.hyperedge(he).connectivity != actual.hyperedge(he).connectivity) {
      connectivity_sets_valid = false;
      break;
    }
  }

  bool num_incident_cut_hes_valid = true;
  for (const HypernodeID hn : actual.nodes()) {
    ASSERT(expected.hypernode(hn).num_incident_cut_hes == actual.hypernode(hn).num_incident_cut_hes
           , V(hn));
    if (expected.hypernode(hn).num_incident_cut_hes != actual.hypernode(hn).num_incident_cut_hes) {
      num_incident_cut_hes_valid = false;
      break;
    }
  }

  return verifyEquivalenceWithoutPartitionInfo(expected, actual) &&
         expected._part_info == actual._part_info &&
         expected._pins_in_part == actual._pins_in_part &&
         num_incident_cut_hes_valid &&
         connectivity_sets_valid;
}


template <typename Hypergraph>
std::pair<std::unique_ptr<Hypergraph>,
          std::vector<typename Hypergraph::HypernodeID> >
reindex(const Hypergraph& hypergraph) {
  using HypernodeID = typename Hypergraph::HypernodeID;
  using HyperedgeID = typename Hypergraph::HyperedgeID;

  std::unordered_map<HypernodeID, HypernodeID> original_to_reindexed;
  std::vector<HypernodeID> reindexed_to_original;
  std::unique_ptr<Hypergraph> reindexed_hypergraph(new Hypergraph());

  reindexed_hypergraph->_k = hypergraph._k;

  HypernodeID num_hypernodes = 0;
  for (const HypernodeID hn : hypergraph.nodes()) {
    original_to_reindexed[hn] = reindexed_to_original.size();
    reindexed_to_original.push_back(hn);
    ++num_hypernodes;
  }

  reindexed_hypergraph->_hypernodes.resize(num_hypernodes);
  reindexed_hypergraph->_num_hypernodes = num_hypernodes;

  HyperedgeID num_hyperedges = 0;
  HypernodeID pin_index = 0;
  for (const HyperedgeID he : hypergraph.edges()) {
    reindexed_hypergraph->_hyperedges.emplace_back(0, 0, hypergraph.edgeWeight(he));
    ++reindexed_hypergraph->_num_hyperedges;
    reindexed_hypergraph->_hyperedges[num_hyperedges].setFirstEntry(pin_index);
    for (const HypernodeID pin : hypergraph.pins(he)) {
      reindexed_hypergraph->hyperedge(num_hyperedges).increaseSize();
      reindexed_hypergraph->_incidence_array.push_back(original_to_reindexed[pin]);
      reindexed_hypergraph->hypernode(original_to_reindexed[pin]).increaseSize();
      ++pin_index;
    }
    ++num_hyperedges;
  }

  const HypernodeID num_pins = pin_index;
  reindexed_hypergraph->_num_pins = num_pins;
  reindexed_hypergraph->_current_num_hypernodes = num_hypernodes;
  reindexed_hypergraph->_current_num_hyperedges = num_hyperedges;
  reindexed_hypergraph->_current_num_pins = num_pins;
  reindexed_hypergraph->_type = hypergraph.type();

  reindexed_hypergraph->_incidence_array.resize(hypergraph._k * num_pins);
  reindexed_hypergraph->_pins_in_part.resize(num_hyperedges * hypergraph._k);
  reindexed_hypergraph->_hes_not_containing_u.setSize(num_hyperedges);

  reindexed_hypergraph->_connectivity_sets.initialize(num_hyperedges, hypergraph._k);

  reindexed_hypergraph->hypernode(0).setFirstEntry(num_pins);
  for (HypernodeID i = 0; i < num_hypernodes - 1; ++i) {
    reindexed_hypergraph->hypernode(i + 1).setFirstEntry(
      reindexed_hypergraph->hypernode(i).firstInvalidEntry());
    reindexed_hypergraph->hypernode(i).setSize(0);
    reindexed_hypergraph->hypernode(i).setWeight(hypergraph.nodeWeight(reindexed_to_original[i]));
    reindexed_hypergraph->_total_weight += reindexed_hypergraph->hypernode(i).weight();
  }
  reindexed_hypergraph->hypernode(num_hypernodes - 1).setSize(0);
  reindexed_hypergraph->hypernode(num_hypernodes - 1).setWeight(
    hypergraph.nodeWeight(reindexed_to_original[num_hypernodes - 1]));
  reindexed_hypergraph->_total_weight +=
    reindexed_hypergraph->hypernode(num_hypernodes - 1).weight();

  for (const HyperedgeID he : reindexed_hypergraph->edges()) {
    for (const HypernodeID pin : reindexed_hypergraph->pins(he)) {
      reindexed_hypergraph->_incidence_array[
        reindexed_hypergraph->hypernode(pin).firstInvalidEntry()] = he;
      reindexed_hypergraph->hypernode(pin).increaseSize();
    }
  }

  reindexed_hypergraph->_part_info.resize(reindexed_hypergraph->_k);
  
  reindexed_hypergraph->createNeighborhoodHypergraph();

  return std::make_pair(std::move(reindexed_hypergraph), reindexed_to_original);
}

template <typename Hypergraph>
std::pair<std::unique_ptr<Hypergraph>,
          std::vector<typename Hypergraph::HypernodeID> >
extractPartAsUnpartitionedHypergraphForBisection(const Hypergraph& hypergraph,
                                                 const typename Hypergraph::PartitionID part,
                                                 const bool split_nets = false) {
  using HypernodeID = typename Hypergraph::HypernodeID;
  using HyperedgeID = typename Hypergraph::HyperedgeID;

  std::unordered_map<HypernodeID, HypernodeID> hypergraph_to_subhypergraph;
  std::vector<HypernodeID> subhypergraph_to_hypergraph;
  std::unique_ptr<Hypergraph> subhypergraph(new Hypergraph());

  HypernodeID num_hypernodes = 0;
  for (const HypernodeID hn : hypergraph.nodes()) {
    if (hypergraph.partID(hn) == part) {
      hypergraph_to_subhypergraph[hn] = subhypergraph_to_hypergraph.size();
      subhypergraph_to_hypergraph.push_back(hn);
      ++num_hypernodes;
    }
  }

  if (num_hypernodes > 0) {
    subhypergraph->_hypernodes.resize(num_hypernodes);
    subhypergraph->_num_hypernodes = num_hypernodes;

    HyperedgeID num_hyperedges = 0;
    HypernodeID pin_index = 0;
    if (split_nets) {
      // Cut-Net Splitting is used to optimize connectivity-1 metric.
      for (const HyperedgeID he : hypergraph.edges()) {
        ASSERT(hypergraph.edgeSize(he) > 1, V(he));
        if (hypergraph.connectivity(he) == 1 && *hypergraph.connectivitySet(he).begin() != part) {
          // HE is completely contained in one of the other parts
          continue;
        }
        if (hypergraph.pinCountInPart(he, part) <= 1) {
          // Single-node HEs have to be removed
          // Furthermore the hyperedge might span other parts
          // (not including the part to be extracted)
          continue;
        }
        subhypergraph->_hyperedges.emplace_back(0, 0, hypergraph.edgeWeight(he));
        ++subhypergraph->_num_hyperedges;
        subhypergraph->_hyperedges[num_hyperedges].setFirstEntry(pin_index);
        for (const HypernodeID pin : hypergraph.pins(he)) {
          if (hypergraph.partID(pin) == part) {
            subhypergraph->hyperedge(num_hyperedges).increaseSize();
            subhypergraph->_incidence_array.push_back(hypergraph_to_subhypergraph[pin]);
            subhypergraph->hypernode(hypergraph_to_subhypergraph[pin]).increaseSize();
            ++pin_index;
          }
        }
        ASSERT(subhypergraph->hyperedge(num_hyperedges).size() > 1, V(he));
        ++num_hyperedges;
      }
    } else {
      for (const HyperedgeID he : hypergraph.edges()) {
        if (hypergraph.connectivity(he) > 1) {
          continue;
        }
        if (*hypergraph.connectivitySet(he).begin() == part) {
          ASSERT(hypergraph.hyperedge(he).connectivity == 1,
                 V(he) << V(hypergraph.hyperedge(he).connectivity));
          ASSERT(hypergraph.edgeSize(he) > 1, V(he));
          subhypergraph->_hyperedges.emplace_back(0, 0, hypergraph.edgeWeight(he));
          ++subhypergraph->_num_hyperedges;
          subhypergraph->_hyperedges[num_hyperedges].setFirstEntry(pin_index);
          for (const HypernodeID pin : hypergraph.pins(he)) {
            ASSERT(hypergraph.partID(pin) == part, V(pin));
            subhypergraph->hyperedge(num_hyperedges).increaseSize();
            subhypergraph->_incidence_array.push_back(hypergraph_to_subhypergraph[pin]);
            subhypergraph->hypernode(hypergraph_to_subhypergraph[pin]).increaseSize();
            ++pin_index;
          }
          ++num_hyperedges;
        }
      }
    }


    const HypernodeID num_pins = pin_index;
    subhypergraph->_num_pins = num_pins;
    subhypergraph->_current_num_hypernodes = num_hypernodes;
    subhypergraph->_current_num_hyperedges = num_hyperedges;
    subhypergraph->_current_num_pins = num_pins;
    subhypergraph->_type = hypergraph.type();

    subhypergraph->_incidence_array.resize(2 * num_pins);
    subhypergraph->_pins_in_part.resize(num_hyperedges * 2);
    subhypergraph->_hes_not_containing_u.setSize(num_hyperedges);

    subhypergraph->_connectivity_sets.initialize(num_hyperedges, 2);

    subhypergraph->hypernode(0).setFirstEntry(num_pins);
    for (HypernodeID i = 0; i < num_hypernodes - 1; ++i) {
      subhypergraph->hypernode(i + 1).setFirstEntry(subhypergraph->hypernode(i).firstInvalidEntry());
      subhypergraph->hypernode(i).setSize(0);
      subhypergraph->hypernode(i).setWeight(hypergraph.nodeWeight(subhypergraph_to_hypergraph[i]));
      subhypergraph->_total_weight += subhypergraph->hypernode(i).weight();
    }
    subhypergraph->hypernode(num_hypernodes - 1).setSize(0);
    subhypergraph->hypernode(num_hypernodes - 1).setWeight(
      hypergraph.nodeWeight(subhypergraph_to_hypergraph[num_hypernodes - 1]));
    subhypergraph->_total_weight += subhypergraph->hypernode(num_hypernodes - 1).weight();

    for (const HyperedgeID he : subhypergraph->edges()) {
      for (const HypernodeID pin : subhypergraph->pins(he)) {
        subhypergraph->_incidence_array[subhypergraph->hypernode(pin).firstInvalidEntry()] = he;
        subhypergraph->hypernode(pin).increaseSize();
      }
    }
  }
  
  subhypergraph->createNeighborhoodHypergraph();
  
  return std::make_pair(std::move(subhypergraph),
                        subhypergraph_to_hypergraph);
}
}  // namespace datastructure
#endif  // SRC_LIB_DATASTRUCTURE_GENERICHYPERGRAPH_H_
