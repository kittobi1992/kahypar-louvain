/***************************************************************************
 *  Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_PARTITION_COARSENING_COARSENERBASE_H_
#define SRC_PARTITION_COARSENING_COARSENERBASE_H_

#include <algorithm>
#include <map>
#include <stack>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "lib/core/Mandatory.h"
#include "lib/definitions.h"
#include "lib/utils/Stats.h"
#include "partition/Configuration.h"
#include "partition/Metrics.h"
#include "partition/coarsening/HypergraphPruner.h"
#include "partition/refinement/IRefiner.h"

using defs::Hypergraph;
using defs::HypernodeID;
using defs::HyperedgeID;
using defs::HypernodeWeight;
using defs::HyperedgeWeight;
using utils::Stats;

namespace partition {
static const bool dbg_coarsening_coarsen = false;
static const bool dbg_coarsening_rating = false;
static const bool dbg_coarsening_no_valid_contraction = false;
static const bool dbg_coarsening_uncoarsen = false;
static const bool dbg_coarsening_uncoarsen_improvement = false;

template <class CoarseningMemento = Mandatory>
class CoarsenerBase {
 protected:
  struct CurrentMaxNodeWeight {
    CurrentMaxNodeWeight(const HypernodeID num_hns, const HypernodeWeight weight) :
      num_nodes(num_hns),
      max_weight(weight) { }
    HypernodeID num_nodes;
    HypernodeWeight max_weight;
  };

 public:
  CoarsenerBase(Hypergraph& hypergraph, const Configuration& config,
                const HypernodeWeight weight_of_heaviest_node) noexcept :
    _hg(hypergraph),
    _config(config),
    _history(),
    _max_hn_weights(),
    _hypergraph_pruner(_hg.initialNumNodes()) {
    _history.reserve(_hg.initialNumNodes());
    _max_hn_weights.reserve(_hg.initialNumNodes());
    _max_hn_weights.emplace_back(_hg.numNodes(), weight_of_heaviest_node);
  }

  virtual ~CoarsenerBase() { }

  CoarsenerBase(const CoarsenerBase&) = delete;
  CoarsenerBase& operator= (const CoarsenerBase&) = delete;

  CoarsenerBase(CoarsenerBase&&) = delete;
  CoarsenerBase& operator= (CoarsenerBase&&) = delete;

 protected:
  void removeSingleNodeHyperedges(const HypernodeID rep_node) noexcept {
    HyperedgeWeight removed_he_weight =
      _hypergraph_pruner.removeSingleNodeHyperedges(_hg, rep_node,
                                                    _history.back().one_pin_hes_begin,
                                                    _history.back().one_pin_hes_size);
    Stats::instance().add(_config, "removedSingleNodeHEWeight", removed_he_weight);
  }

  void removeParallelHyperedges(const HypernodeID rep_node) noexcept {
    HyperedgeID removed_parallel_hes =
      _hypergraph_pruner.removeParallelHyperedges(_hg, rep_node,
                                                  _history.back().parallel_hes_begin,
                                                  _history.back().parallel_hes_size);
    Stats::instance().add(_config, "numRemovedParalellHEs", removed_parallel_hes);
  }

  void restoreParallelHyperedges() noexcept {
    _hypergraph_pruner.restoreParallelHyperedges(_hg,
                                                 _history.back().parallel_hes_begin,
                                                 _history.back().parallel_hes_size);
  }

  void restoreSingleNodeHyperedges() noexcept {
    _hypergraph_pruner.restoreSingleNodeHyperedges(_hg,
                                                   _history.back().one_pin_hes_begin,
                                                   _history.back().one_pin_hes_size);
  }

  void initializeRefiner(IRefiner& refiner) noexcept {
#ifdef USE_BUCKET_PQ
    HyperedgeID max_degree = 0;
    for (const HypernodeID hn : _hg.nodes()) {
      max_degree = std::max(max_degree, _hg.nodeDegree(hn));
    }
    HyperedgeWeight max_he_weight = 0;
    for (const HyperedgeID he : _hg.edges()) {
      max_he_weight = std::max(max_he_weight, _hg.edgeWeight(he));
    }
    LOGVAR(max_degree);
    LOGVAR(max_he_weight);
    refiner.initialize(static_cast<HyperedgeWeight>(max_degree * max_he_weight));
#else
    refiner.initialize();
#endif
  }

  void performLocalSearch(IRefiner& refiner, std::vector<HypernodeID>& refinement_nodes,
                          Metrics& current_metrics,
                          UncontractionGainChanges& changes) noexcept {
    ASSERT(changes.representative.size() != 0, "0");
    ASSERT(changes.contraction_partner.size() != 0, "0");

    bool improvement_found = performLocalSearchIteration(refiner, refinement_nodes, changes,
                                                         current_metrics);

    UncontractionGainChanges no_changes;
    no_changes.representative.push_back(0);
    no_changes.contraction_partner.push_back(0);

    int iteration = 1;
    while ((iteration < _config.fm_local_search.num_repetitions) && improvement_found) {
      improvement_found = performLocalSearchIteration(refiner, refinement_nodes, no_changes,
                                                      current_metrics);
      ++iteration;
    }
  }

  bool performLocalSearchIteration(IRefiner& refiner,
                                   std::vector<HypernodeID>& refinement_nodes,
                                   UncontractionGainChanges& current_changes,
                                   Metrics& current_metrics) {
    const HyperedgeWeight old_cut = current_metrics.cut;
    const HyperedgeWeight old_km1 = current_metrics.km1;
    bool improvement_found = refiner.refine(refinement_nodes,
                                            { _config.partition.max_part_weights[0]
                                              + _max_hn_weights.back().max_weight,
                                              _config.partition.max_part_weights[1]
                                              + _max_hn_weights.back().max_weight },
                                            current_changes,
                                            current_metrics);
    ASSERT(_config.partition.objective != Objective::cut || current_metrics.cut <= old_cut,
           V(current_metrics.cut) << V(old_cut));
    ASSERT(_config.partition.objective != Objective::km1 ||
           current_metrics.km1 <= old_km1,
           V(current_metrics.km1) << V(old_km1));
    ASSERT(current_metrics.cut == metrics::hyperedgeCut(_hg), "Inconsistent cut values");
    DBG(dbg_coarsening_uncoarsen && (_config.partition.objective == Objective::cut),
        old_cut << "-->" << current_metrics.cut);
    DBG(dbg_coarsening_uncoarsen && (_config.partition.objective == Objective::km1),
        old_km1 << "-->" << current_metrics.km1);
    return improvement_found;
  }

  Hypergraph& _hg;
  const Configuration& _config;
  std::vector<CoarseningMemento> _history;
  std::vector<CurrentMaxNodeWeight> _max_hn_weights;
  HypergraphPruner _hypergraph_pruner;
};
}  // namespace partition

#endif  // SRC_PARTITION_COARSENING_COARSENERBASE_H_
