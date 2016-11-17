/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <string>
#include <utility>
#include <vector>

#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/definitions.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/meta/template_parameter_to_string.h"
#include "kahypar/partition/coarsening/i_coarsener.h"
#include "kahypar/partition/coarsening/vertex_pair_coarsener_base.h"
#include "kahypar/utils/stats.h"

namespace kahypar {
template <class Rater = Mandatory>
class LazyVertexPairCoarsener final : public ICoarsener,
                                      private VertexPairCoarsenerBase<>{
 private:
  using Base = VertexPairCoarsenerBase;
  using Rating = typename Rater::Rating;

 public:
  LazyVertexPairCoarsener(Hypergraph& hypergraph, const Configuration& config,
                          const HypernodeWeight weight_of_heaviest_node) :
    Base(hypergraph, config, weight_of_heaviest_node),
    _rater(_hg, _config),
    _outdated_rating(hypergraph.initialNumNodes()),
    _target(_hg.initialNumNodes()) { }

  virtual ~LazyVertexPairCoarsener() { }

  LazyVertexPairCoarsener(const LazyVertexPairCoarsener&) = delete;
  LazyVertexPairCoarsener& operator= (const LazyVertexPairCoarsener&) = delete;

  LazyVertexPairCoarsener(LazyVertexPairCoarsener&&) = delete;
  LazyVertexPairCoarsener& operator= (LazyVertexPairCoarsener&&) = delete;

 private:
  FRIEND_TEST(ALazyUpdateCoarsener, InvalidatesAdjacentHypernodesInsteadOfReratingThem);

  void coarsenImpl(const HypernodeID limit) override final {

    int pass_nr = 0;
    
    while (_hg.currentNumNodes() > limit) {
      
      LOGVAR(pass_nr);
      LOGVAR(_hg.currentNumNodes());
      LOGVAR(_hg.currentNumEdges());
      
      _pq.clear();
      if(_config.preprocessing.use_louvain) {
        _rater.performLouvainCommunityDetection();         
      }
      rateAllHypernodes(_rater, _target);
      
      const HypernodeID num_hns_before_pass = _hg.currentNumNodes();
      
      while(!_pq.empty() && _hg.currentNumNodes() > limit) {
        const HypernodeID rep_node = _pq.top();

        if (_outdated_rating[rep_node]) {
          DBG(dbg_coarsening_coarsen,
              "Rating for HN" << rep_node << " is invalid: " << _pq.topKey() << "--->"
              << _rater.rate(rep_node).value << " target=" << _rater.rate(rep_node).target
              << " valid= " << _rater.rate(rep_node).valid);
          updatePQandContractionTarget(rep_node, _rater.rate(rep_node));
        } else {
          const HypernodeID contracted_node = _target[rep_node];

          DBG(dbg_coarsening_coarsen, "Contracting: (" << rep_node << ","
              << _target[rep_node] << ") prio: " << _pq.topKey() <<
              " deg(" << rep_node << ")=" << _hg.nodeDegree(rep_node) <<
              " deg(" << contracted_node << ")=" << _hg.nodeDegree(contracted_node));

          ASSERT(_hg.nodeWeight(rep_node) + _hg.nodeWeight(_target[rep_node])
                <= _rater.thresholdNodeWeight());
          ASSERT(_pq.topKey() == _rater.rate(rep_node).value,
                V(_pq.topKey()) << V(_rater.rate(rep_node).value));

          performContraction(rep_node, contracted_node);
          if(_config.preprocessing.use_multilevel_louvain) {
              _rater.contractHypernodes(rep_node,contracted_node);
              size_t N = _hg.initialNumNodes();
              int one_pin_hes_begin = _history.back().one_pin_hes_begin;
              int one_pin_hes_size = _history.back().one_pin_hes_size;
              for(int i = one_pin_hes_begin; i < one_pin_hes_begin+one_pin_hes_size; ++i) {
                  _rater.contractHypernodes(rep_node,N+_hypergraph_pruner.removedSingleNodeHyperedges()[i]);
              }
              
              int parallel_hes_begin = _history.back().parallel_hes_begin;
              int parallel_hes_size = _history.back().parallel_hes_size;
              for(int i = parallel_hes_begin; i < parallel_hes_begin+parallel_hes_size; ++i) {
                  _rater.contractHypernodes(N+_hypergraph_pruner.removedParallelHyperedges()[i].representative_id,
                                                N+_hypergraph_pruner.removedParallelHyperedges()[i].removed_id);
              }
          }
          ASSERT(_pq.contains(contracted_node), V(contracted_node));
          _pq.remove(contracted_node);

          // this also invalidates rep_node, however rep_node
          // will be re-rated and updated afterwards
          invalidateAffectedHypernodes(rep_node);

          updatePQandContractionTarget(rep_node, _rater.rate(rep_node));
        }
      }
      
      if (num_hns_before_pass == _hg.currentNumNodes()) {
        break;
      }
      
      ++pass_nr;
    }
  }

  
  bool uncoarsenImpl(IRefiner& refiner) override final {
    return Base::doUncoarsen(refiner);
  }

  std::string policyStringImpl() const override final {
    return std::string(" ratingFunction=" + meta::templateToString<Rater>());
  }

  void invalidateAffectedHypernodes(const HypernodeID rep_node) {
    for (const HyperedgeID he : _hg.incidentEdges(rep_node)) {
      for (const HypernodeID pin : _hg.pins(he)) {
        _outdated_rating.set(pin, true);
      }
    }
  }

  void updatePQandContractionTarget(const HypernodeID hn, const Rating& rating) {
    _outdated_rating.set(hn, false);
    if (rating.valid) {
      ASSERT(_pq.contains(hn), V(hn));
      _pq.updateKey(hn, rating.value);
      _target[hn] = rating.target;
    } else {
      // In this case, no explicit contaiment check is necessary because the
      // method is only called on rep_node, which is definetly in the PQ.
      ASSERT(_pq.contains(hn), V(hn));
      _pq.remove(hn);
      DBG(dbg_coarsening_no_valid_contraction, "Progress [" << _hg.currentNumNodes() << "/"
          << _hg.initialNumNodes() << "]:HN " << hn
          << " \t(w=" << _hg.nodeWeight(hn) << "," << " deg=" << _hg.nodeDegree(hn)
          << ") did not find valid contraction partner.");
#ifdef GATHER_STATS
      Stats::instance().add(_config, "numHNsWithoutValidContractionPartner", 1);
#endif
    }
  }

  using Base::_pq;
  using Base::_hg;
  using Base::_config;
  using Base::_history;
  Rater _rater;
  ds::FastResetFlagArray<> _outdated_rating;
  std::vector<HypernodeID> _target;
};
}              // namespace kahypar
