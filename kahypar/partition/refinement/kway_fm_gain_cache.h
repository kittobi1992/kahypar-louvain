/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <limits>
#include <memory>
#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/partition/refinement/gain_cache_element.h"

namespace kahypar {
template <typename Gain = Mandatory>
class KwayGainCache {
 private:
  static const bool debug = false;
  static const HypernodeID hn_to_debug = 2225;

  using Byte = char;
  using KFMCacheElement = CacheElement<Gain>;

 public:
  static constexpr HyperedgeWeight kNotCached = KFMCacheElement::kNotCached;
  static constexpr PartitionID kInvalidPart = KFMCacheElement::kInvalidPart;

  KwayGainCache(const HypernodeID num_hns, const PartitionID k) :
    _k(k),
    _num_hns(num_hns),
    _cache_element_size(sizeof(KFMCacheElement) +
                        _k * sizeof(typename KFMCacheElement::Element) +
                        _k * sizeof(PartitionID)),
    _cache(nullptr),
    _deltas() {
    _cache = static_cast<Byte*>(malloc(num_hns * _cache_element_size));
    for (HypernodeID hn = 0; hn < _num_hns; ++hn) {
      new(cacheElement(hn))KFMCacheElement(k);
    }
  }

  ~KwayGainCache() {
    // Since CacheElemment only contains PartitionIDs and these are PODs,
    // we do not need to call destructors of CacheElement cacheElement(i)->~CacheElement();
    static_assert(std::is_pod<Gain>::value, "Gain is not a POD");
    static_assert(std::is_pod<PartitionID>::value, "PartitionID is not a POD");
    free(_cache);
  }

  KwayGainCache(const KwayGainCache&) = delete;
  KwayGainCache& operator= (const KwayGainCache&) = delete;

  KwayGainCache(KwayGainCache&&) = default;
  KwayGainCache& operator= (KwayGainCache&&) = default;

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Gain entry(const HypernodeID hn, const PartitionID part) const {
    DBG(debug && (hn == hn_to_debug), "entry access for HN " << hn << " and part " << part);
    ASSERT(part < _k, V(part));
    return cacheElement(hn)->gain(part);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool entryExists(const HypernodeID hn,
                                                   const PartitionID part) const {
    ASSERT(part < _k, V(part));
    DBG(debug && (hn == hn_to_debug), "existence check for HN " << hn << " and part " << part
        << "=" << cacheElement(hn)->contains(part));
    return cacheElement(hn)->contains(part);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void removeEntryDueToConnectivityDecrease(const HypernodeID hn,
                                                                            const PartitionID part) {
    ASSERT(part < _k, V(part));
    _deltas.emplace_back(hn, part, cacheElement(hn)->gain(part), RollbackAction::do_add);
    DBG(debug && (hn == hn_to_debug), "removeEntryDueToConnectivityDecrease for " << hn
        << "and part " << part << " previous cache entry =" << cacheElement(hn)->gain(part)
        << " now= " << kNotCached);
    cacheElement(hn)->remove(part);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void addEntryDueToConnectivityIncrease(const HypernodeID hn,
                                                                         const PartitionID part,
                                                                         const Gain gain) {
    ASSERT(part < _k, V(part));
    ASSERT(!entryExists(hn, part), V(hn) << V(part));
    cacheElement(hn)->add(part, gain);
    DBG(debug && (hn == hn_to_debug), "addEntryDueToConnectivityIncrease for " << hn
        << "and part " << part << " new cache entry =" << cacheElement(hn)->gain(part));
    DBG(debug && (hn == hn_to_debug), "delta emplace = " << kNotCached - gain);
    _deltas.emplace_back(hn, part, kNotCached - gain, RollbackAction::do_remove);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updateFromAndToPartOfMovedHN(const HypernodeID moved_hn,
                                                                    const PartitionID from_part,
                                                                    const PartitionID to_part,
                                                                    const bool remains_connected_to_from_part) {
    if (remains_connected_to_from_part) {
      DBG(debug && (moved_hn == hn_to_debug),
          "updateFromAndToPartOfMovedHN(" << moved_hn << "," << from_part << "," << to_part << ")");
      const Gain to_part_gain = cacheElement(moved_hn)->gain(to_part);
      _deltas.emplace_back(moved_hn, from_part,
                           cacheElement(moved_hn)->gain(from_part) + to_part_gain,
                           RollbackAction::do_remove);
      cacheElement(moved_hn)->add(from_part, -to_part_gain);
    } else {
      DBG(debug && (moved_hn == hn_to_debug), "pseudoremove  for " << moved_hn
          << "and part " << from_part << " previous cache entry ="
          << cacheElement(moved_hn)->gain(from_part)
          << " now= " << kNotCached);
      ASSERT(cacheElement(moved_hn)->gain(from_part) == kNotCached, V(moved_hn) << V(from_part));
    }
    removeEntryDueToConnectivityDecrease(moved_hn, to_part);
  }


  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void clear(const HypernodeID hn) {
    DBG(debug && (hn == hn_to_debug), "clear(" << hn << ")");
    cacheElement(hn)->clear();
  }

  void initializeEntry(const HypernodeID hn, const PartitionID part, const Gain value) {
    ASSERT(part < _k, V(part));
    DBG(debug && (hn == hn_to_debug),
        "initializeEntry(" << hn << "," << part << "," << value << ")");
    cacheElement(hn)->add(part, value);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updateEntryIfItExists(const HypernodeID hn,
                                                             const PartitionID part,
                                                             const Gain delta) {
    ASSERT(part < _k, V(part));
    if (entryExists(hn, part)) {
      updateExistingEntry(hn, part, delta);
    }
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updateExistingEntry(const HypernodeID hn,
                                                           const PartitionID part,
                                                           const Gain delta) {
    ASSERT(part < _k, V(part));
    ASSERT(entryExists(hn, part), V(hn) << V(part));
    ASSERT(cacheElement(hn)->gain(part) != kNotCached, V(hn) << V(part));
    DBG(debug && (hn == hn_to_debug),
        "updateEntryAndDelta(" << hn << ", " << part << "," << delta << ")");
    cacheElement(hn)->update(part, delta);
    _deltas.emplace_back(hn, part, -delta, RollbackAction::do_nothing);
  }


  void rollbackDelta() {
    for (auto rit = _deltas.crbegin(); rit != _deltas.crend(); ++rit) {
      const HypernodeID hn = rit->hn;
      const PartitionID part = rit->part;
      const Gain delta = rit->delta;
      if (cacheElement(hn)->contains(part)) {
        DBG(debug && (hn == hn_to_debug), "rollback: " << "G[" << hn << "," << part << "]="
            << cacheElement(hn)->gain(part) << "+" << delta << "="
            << (cacheElement(hn)->gain(part) + delta));
        cacheElement(hn)->update(part, delta);
        if (rit->action == RollbackAction::do_remove) {
          ASSERT(cacheElement(hn)->gain(part) == kNotCached, V(hn));
          cacheElement(hn)->remove(part);
        }
      } else {
        DBG(debug && (hn == hn_to_debug), "rollback: SET " << "set G[" << hn << "," << part << "]="
            << delta);
        cacheElement(hn)->set(part, delta);
        if (rit->action == RollbackAction::do_add) {
          cacheElement(hn)->addToActiveParts(part);
        }
      }
    }
    _deltas.clear();
  }

  void resetDelta() {
    _deltas.clear();
  }

  const KFMCacheElement & adjacentParts(const HypernodeID hn) const {
    return *cacheElement(hn);
  }

  void clear() {
    for (HypernodeID hn = 0; hn < _num_hns; ++hn) {
      new(cacheElement(hn))KFMCacheElement(_k);
    }
  }

 private:
  const KFMCacheElement* cacheElement(const HypernodeID hn) const {
    return reinterpret_cast<KFMCacheElement*>(_cache + hn * _cache_element_size);
  }

  // To avoid code duplication we implement non-const version in terms of const version
  KFMCacheElement* cacheElement(const HypernodeID he) {
    return const_cast<KFMCacheElement*>(static_cast<const KwayGainCache&>(*this).cacheElement(he));
  }

  PartitionID _k;
  HypernodeID _num_hns;
  const size_t _cache_element_size;
  Byte* _cache;
  std::vector<RollbackElement> _deltas;
};
}  // namespace kahypar
