/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
******************************************************************************/

#pragma once

#include <algorithm>
#include <map>
#include <stack>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/partition/coarsening/coarsening_memento.h"
#include "kahypar/partition/coarsening/hypergraph_pruner.h"
#include "kahypar/partition/configuration.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/utils/stats.h"

namespace kahypar {
static const bool dbg_coarsening_coarsen = false;
static const bool dbg_coarsening_rating = false;
static const bool dbg_coarsening_no_valid_contraction = false;
static const bool dbg_coarsening_uncoarsen = false;
static const bool dbg_coarsening_uncoarsen_improvement = false;

class CoarsenerBase {
 private:
  struct CurrentMaxNodeWeight {
    CurrentMaxNodeWeight(const HypernodeID num_hns, const HypernodeWeight weight) :
      num_nodes(num_hns),
      max_weight(weight) { }
    HypernodeID num_nodes;
    HypernodeWeight max_weight;
  };

 public:
  CoarsenerBase(Hypergraph& hypergraph, const Configuration& config,
                const HypernodeWeight weight_of_heaviest_node) :
    _hg(hypergraph),
    _config(config),
    _history(),
    _max_hn_weights(),
    _hypergraph_pruner(_hg.initialNumNodes()) {
    _history.reserve(_hg.initialNumNodes());
    _max_hn_weights.reserve(_hg.initialNumNodes());
    _max_hn_weights.emplace_back(_hg.initialNumNodes(), weight_of_heaviest_node);
  }

  virtual ~CoarsenerBase() { }

  CoarsenerBase(const CoarsenerBase&) = delete;
  CoarsenerBase& operator= (const CoarsenerBase&) = delete;

  CoarsenerBase(CoarsenerBase&&) = delete;
  CoarsenerBase& operator= (CoarsenerBase&&) = delete;

 protected:
  void performContraction(const HypernodeID rep_node, const HypernodeID contracted_node) {
    _history.emplace_back(_hg.contract(rep_node, contracted_node));
    if (_hg.nodeWeight(rep_node) > _max_hn_weights.back().max_weight) {
      _max_hn_weights.emplace_back(_hg.currentNumNodes(), _hg.nodeWeight(rep_node));
    }
    removeSingleNodeHyperedges();
    removeParallelHyperedges();
  }

  void removeSingleNodeHyperedges() {
    const HyperedgeWeight removed_he_weight =
      _hypergraph_pruner.removeSingleNodeHyperedges(_hg, _history.back());
    Stats::instance().add(_config, "removedSingleNodeHEWeight", removed_he_weight);
  }

  void removeParallelHyperedges() {
    const HyperedgeID removed_parallel_hes =
      _hypergraph_pruner.removeParallelHyperedges(_hg, _history.back());
    Stats::instance().add(_config, "numRemovedParalellHEs", removed_parallel_hes);
  }

  void restoreParallelHyperedges() {
    _hypergraph_pruner.restoreParallelHyperedges(_hg, _history.back());
  }

  void restoreSingleNodeHyperedges() {
    _hypergraph_pruner.restoreSingleNodeHyperedges(_hg, _history.back());
  }

  void initializeRefiner(IRefiner& refiner) {
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
                          const UncontractionGainChanges& changes) {
    ASSERT(changes.representative.size() != 0, "0");
    ASSERT(changes.contraction_partner.size() != 0, "0");
    bool improvement_found = performLocalSearchIteration(refiner, refinement_nodes, changes,
                                                         current_metrics);
    UncontractionGainChanges no_changes;
    no_changes.representative.push_back(0);
    no_changes.contraction_partner.push_back(0);

    int iteration = 1;
    while ((iteration < _config.local_search.iterations_per_level) && improvement_found) {
      improvement_found = performLocalSearchIteration(refiner, refinement_nodes, no_changes,
                                                      current_metrics);
      ++iteration;
    }
  }

  bool performLocalSearchIteration(IRefiner& refiner,
                                   std::vector<HypernodeID>& refinement_nodes,
                                   const UncontractionGainChanges& current_changes,
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

    ASSERT(_config.partition.objective != Objective::cut ||
           (current_metrics.cut <= old_cut && current_metrics.cut == metrics::hyperedgeCut(_hg)),
           V(current_metrics.cut) << V(old_cut));
    // Km1 is optimized indirectly in recursive bisection mode via cut-net splitting and optimizing
    // cut. In this case current_metrics.km1 is not used.
    ASSERT((_config.partition.mode != Mode::direct_kway ||
            _config.partition.objective != Objective::km1) ||
           (current_metrics.km1 <= old_km1 && current_metrics.km1 == metrics::km1(_hg)),
           V(current_metrics.km1) << V(old_km1) << V(metrics::km1(_hg)));

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
}  // namespace kahypar
