/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2015 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "gtest/gtest.h"

#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "tools/cnf_to_hgr_conversion.h"

using::testing::Test;

using namespace kahypar;

namespace cnfconversion {
TEST(ACnfToHgrConversionRoutine, ConvertsCNFinstancesIntoHypergraphInstances) {
  std::string cnf_filename("test_instances/SampleSAT.cnf");
  std::string hgr_filename("test_instances/SampleSAT.cnf.hgr");
  convertInstance(cnf_filename, hgr_filename);

  Hypergraph hypergraph = io::createHypergraphFromFile(hgr_filename, 2);

  ASSERT_EQ(hypergraph.initialNumNodes(), 8);
  ASSERT_EQ(hypergraph.initialNumPins(), 9);
  ASSERT_EQ(hypergraph.currentNumEdges(), 3);

  std::vector<HypernodeID> pins_he_0({ 0, 1, 2 });
  size_t i = 0;
  for (const HypernodeID pin : hypergraph.pins(0)) {
    ASSERT_EQ(pin, pins_he_0[i++]);
  }

  i = 0;
  std::vector<HypernodeID> pins_he_1({ 3, 4, 5, 2 });
  for (const HypernodeID pin : hypergraph.pins(1)) {
    ASSERT_EQ(pin, pins_he_1[i++]);
  }

  i = 0;
  std::vector<HypernodeID> pins_he_2({ 6, 7 });
  for (const HypernodeID pin : hypergraph.pins(2)) {
    ASSERT_EQ(pin, pins_he_2[i++]);
  }
}
}  // namespace cnfconversion
