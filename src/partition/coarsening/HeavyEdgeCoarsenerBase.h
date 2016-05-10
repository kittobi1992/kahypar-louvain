/***************************************************************************
 *  Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_PARTITION_COARSENING_HEAVYEDGECOARSENERBASE_H_
#define SRC_PARTITION_COARSENING_HEAVYEDGECOARSENERBASE_H_

#include <algorithm>
#include <limits>
#include <stack>
#include <unordered_map>
#include <utility>
#include <vector>

#include "lib/core/Int2Type.h"
#include "lib/datastructure/heaps/NoDataBinaryMaxHeap.h"
#include "lib/definitions.h"
#include "partition/Configuration.h"
#include "partition/Metrics.h"
#include "partition/coarsening/CoarsenerBase.h"
#include "partition/coarsening/Rater.h"
#include "partition/refinement/IRefiner.h"

using core::Int2Type;
using datastructure::NoDataBinaryMaxHeap;
using defs::Hypergraph;
using defs::HypernodeID;
using defs::HyperedgeID;
using defs::HypernodeWeight;
using defs::HyperedgeWeight;

namespace partition {
struct CoarseningMemento {
  int one_pin_hes_begin;        // start of removed single pin hyperedges
  int one_pin_hes_size;         // # removed single pin hyperedges
  int parallel_hes_begin;       // start of removed parallel hyperedges
  int parallel_hes_size;        // # removed parallel hyperedges
  Hypergraph::ContractionMemento contraction_memento;
  explicit CoarseningMemento(Hypergraph::ContractionMemento&& contraction_memento_) noexcept :
    one_pin_hes_begin(0),
    one_pin_hes_size(0),
    parallel_hes_begin(0),
    parallel_hes_size(0),
    contraction_memento(std::move(contraction_memento_)) { }
};

template <class Rater = Mandatory,
          class PrioQueue = NoDataBinaryMaxHeap<HypernodeID,
                                                typename Rater::RatingType> >
class HeavyEdgeCoarsenerBase : public CoarsenerBase<CoarseningMemento>{
 protected:
  using Base = CoarsenerBase<CoarseningMemento>;
  using Base::_hg;
  using Base::_config;
  using Base::_history;
  using Base::_max_hn_weights;
  using Base::CurrentMaxNodeWeight;
  using Base::restoreSingleNodeHyperedges;
  using Base::restoreParallelHyperedges;
  using Base::performLocalSearch;
  using Base::initializeRefiner;
  using Rating = typename Rater::Rating;
  using RatingType = typename Rater::RatingType;

 public:
  HeavyEdgeCoarsenerBase(Hypergraph& hypergraph, const Configuration& config,
                         const HypernodeWeight weight_of_heaviest_node) noexcept :
    Base(hypergraph, config, weight_of_heaviest_node),
    _rater(_hg, _config),
    _pq(_hg.initialNumNodes()) { }

  ~HeavyEdgeCoarsenerBase() { }

  HeavyEdgeCoarsenerBase(const HeavyEdgeCoarsenerBase&) = delete;
  HeavyEdgeCoarsenerBase& operator= (const HeavyEdgeCoarsenerBase&) = delete;

  HeavyEdgeCoarsenerBase(HeavyEdgeCoarsenerBase&&) = delete;
  HeavyEdgeCoarsenerBase& operator= (HeavyEdgeCoarsenerBase&&) = delete;

 protected:
  FRIEND_TEST(ACoarsener, SelectsNodePairToContractBasedOnHighestRating);

  void performContraction(const HypernodeID rep_node, const HypernodeID contracted_node) noexcept {
    _history.emplace_back(_hg.contract(rep_node, contracted_node));
    if (_hg.nodeWeight(rep_node) > _max_hn_weights.back().max_weight) {
      _max_hn_weights.emplace_back(_hg.numNodes(), _hg.nodeWeight(rep_node));
    }
  }

  bool doUncoarsen(IRefiner& refiner) noexcept {
    Metrics current_metrics = { metrics::hyperedgeCut(_hg),
                                metrics::kMinus1(_hg),
                                metrics::imbalance(_hg, _config) };
    HyperedgeWeight initial_objective = std::numeric_limits<HyperedgeWeight>::min();

    switch (_config.partition.objective) {
      case Objective::cut:
        initial_objective = current_metrics.cut;
        Stats::instance().add(_config, "initialCut", initial_objective);
        LOG("initial cut =" << initial_objective);
        break;
      case Objective::km1:
        initial_objective = current_metrics.km1;
        Stats::instance().add(_config, "initialKm1", initial_objective);
        LOG("initial km1 =" << initial_objective);
        break;
      default:
        LOG("Unknown Objective");
        exit(-1);
    }

    Stats::instance().add(_config, "initialImbalance", current_metrics.imbalance);
    LOG("initial imbalance=" << current_metrics.imbalance);
    LOG("target  weights (RB): w(0)=" << _config.partition.max_part_weights[0]
        << " w(1)=" << _config.partition.max_part_weights[1]);
    LOG("initial weights (RB): w(0)=" << _hg.partWeight(0) << " w(1)=" << _hg.partWeight(1));

    initializeRefiner(refiner);
    std::vector<HypernodeID> refinement_nodes(2, 0);
    UncontractionGainChanges changes;
    changes.representative.push_back(0);
    changes.contraction_partner.push_back(0);
    while (!_history.empty()) {
      restoreParallelHyperedges();
      restoreSingleNodeHyperedges();

      DBG(dbg_coarsening_uncoarsen, "Uncontracting: (" << _history.back().contraction_memento.u << ","
          << _history.back().contraction_memento.v << ")");

      refinement_nodes.clear();
      refinement_nodes.push_back(_history.back().contraction_memento.u);
      refinement_nodes.push_back(_history.back().contraction_memento.v);

      if (_hg.numNodes() > _max_hn_weights.back().num_nodes) {
        _max_hn_weights.pop_back();
      }

      switch (_config.partition.refinement_algorithm) {
        case RefinementAlgorithm::twoway_fm:
          _hg.uncontract(_history.back().contraction_memento, changes,
                         Int2Type<static_cast<int>(RefinementAlgorithm::twoway_fm)>());
          break;
        default:
          _hg.uncontract(_history.back().contraction_memento);
      }
      performLocalSearch(refiner, refinement_nodes, current_metrics, changes);
      changes.representative[0] = 0;
      changes.contraction_partner[0] = 0;
      _history.pop_back();
    }

    // This currently cannot be guaranteed for RB-partitioning and k != 2^x, since it might be
    // possible that 2FM cannot re-adjust the part weights to be less than Lmax0 and Lmax1.
    // In order to guarantee this, 2FM would have to force rebalancing by sacrificing cut-edges.
    // ASSERT(current_imbalance <= _config.partition.epsilon,
    //        "balance_constraint is violated after uncontraction:" << metrics::imbalance(_hg, _config)
    //        << " > " << _config.partition.epsilon);
    Stats::instance().add(_config, "finalImbalance", current_metrics.imbalance);
    LOG("final weights (RB):   w(0)=" << _hg.partWeight(0) << " w(1)=" << _hg.partWeight(1));
    bool improvement_found = false;
    switch (_config.partition.objective) {
      case Objective::cut:
        Stats::instance().add(_config, "finalCut", current_metrics.cut);
        LOG("final cut: " << current_metrics.cut);
        improvement_found = initial_objective < current_metrics.cut;
        break;
      case Objective::km1:
        if (_config.partition.mode == Mode::recursive_bisection) {
          // In recursive bisection-based (initial) partitioning, km1
          // is optimized using TwoWayFM and cut-net splitting. Since
          // TwoWayFM optimizes cut, current_metrics.km1 is not updated
          // during local search (it is currently only updated/maintained
          // during k-way k-1 refinement). In order to provide correct outputs,
          // we explicitly calculated the metric after uncoarsening.
          current_metrics.km1 = metrics::kMinus1(_hg);
        }
        Stats::instance().add(_config, "finalKm1", current_metrics.km1);
        LOG("final km1: " << current_metrics.km1);
        improvement_found = initial_objective < current_metrics.km1;
        break;
      default:
        LOG("Unknown Objective");
        exit(-1);
    }

    return improvement_found;
  }

  template <typename Map>
  void rateAllHypernodes(std::vector<HypernodeID>& target, Map& sources) noexcept {
    std::vector<HypernodeID> permutation;
    createHypernodePermutation(permutation);
    for (size_t i = 0; i < permutation.size(); ++i) {
      const Rating rating = _rater.rate(permutation[i]);
      if (rating.valid) {
        _pq.push(permutation[i], rating.value);
        target[permutation[i]] = rating.target;
        sources.insert({ rating.target, permutation[i] });
      }
    }
  }

  void createHypernodePermutation(std::vector<HypernodeID>& permutation) noexcept {
    permutation.reserve(_hg.initialNumNodes());
    for (HypernodeID hn : _hg.nodes()) {
      permutation.push_back(hn);
    }
    Randomize::shuffleVector(permutation, permutation.size());
  }

  Rater _rater;
  PrioQueue _pq;
};
}  // namespace partition

#endif  // SRC_PARTITION_COARSENING_HEAVYEDGECOARSENERBASE_H_
