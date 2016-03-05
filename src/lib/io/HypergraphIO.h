/***************************************************************************
 *  Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_LIB_IO_HYPERGRAPHIO_H_
#define SRC_LIB_IO_HYPERGRAPHIO_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "lib/definitions.h"

using defs::PartitionID;
using defs::Hypergraph;
using defs::HypernodeID;
using defs::HyperedgeID;
using defs::HyperedgeWeight;
using defs::HypernodeWeight;
using defs::HyperedgeIndexVector;
using defs::HyperedgeVector;
using defs::HyperedgeWeightVector;
using defs::HypernodeWeightVector;
using defs::HypergraphType;

namespace io {
using Mapping = std::unordered_map<HypernodeID, HypernodeID>;

static inline void readHGRHeader(std::ifstream& file, HyperedgeID& num_hyperedges,
                                 HypernodeID& num_hypernodes, HypergraphType& hypergraph_type) {
  std::string line;
  std::getline(file, line);

  // skip any comments
  while (line[0] == '%') {
    std::getline(file, line);
  }

  std::istringstream sstream(line);
  int i = 0;
  sstream >> num_hyperedges >> num_hypernodes >> i;
  hypergraph_type = static_cast<HypergraphType>(i);
}

static inline void readHypergraphFile(const std::string& filename, HypernodeID& num_hypernodes,
                                      HyperedgeID& num_hyperedges,
                                      HyperedgeIndexVector& index_vector,
                                      HyperedgeVector& edge_vector,
                                      HyperedgeWeightVector* hyperedge_weights = nullptr,
                                      HypernodeWeightVector* hypernode_weights = nullptr) {
  ASSERT(!filename.empty(), "No filename for hypergraph file specified");
  HypergraphType hypergraph_type = HypergraphType::Unweighted;
  std::ifstream file(filename);
  if (file) {
    readHGRHeader(file, num_hyperedges, num_hypernodes, hypergraph_type);
    ASSERT(hypergraph_type == HypergraphType::Unweighted ||
           hypergraph_type == HypergraphType::EdgeWeights ||
           hypergraph_type == HypergraphType::NodeWeights ||
           hypergraph_type == HypergraphType::EdgeAndNodeWeights,
           "Hypergraph in file has wrong type");

    bool has_hyperedge_weights = hypergraph_type == HypergraphType::EdgeWeights ||
                                 hypergraph_type == HypergraphType::EdgeAndNodeWeights ?
                                 true : false;
    bool has_hypernode_weights = hypergraph_type == HypergraphType::NodeWeights ||
                                 hypergraph_type == HypergraphType::EdgeAndNodeWeights ?
                                 true : false;

    index_vector.reserve(num_hyperedges +  /*sentinel*/ 1);
    index_vector.push_back(edge_vector.size());

    std::string line;
    for (HyperedgeID i = 0; i < num_hyperedges; ++i) {
      std::getline(file, line);
      std::istringstream line_stream(line);
      if (line_stream.peek() == EOF) {
        std::cerr << "Error: Hyperedge " << i << " is empty" << std::endl;
        exit(1);
      }

      if (has_hyperedge_weights) {
        HyperedgeWeight edge_weight;
        line_stream >> edge_weight;
        if (hyperedge_weights == nullptr) {
          LOG(" ****** ignoring hyperedge weights ******");
        } else {
          ASSERT(hyperedge_weights != nullptr, "Hypergraph has hyperedge weights");
          hyperedge_weights->push_back(edge_weight);
        }
      }
      HypernodeID pin;
      while (line_stream >> pin) {
        // Hypernode IDs start from 0
        --pin;
        ASSERT(pin < num_hypernodes, "Invalid hypernode ID");
        edge_vector.push_back(pin);
      }
      index_vector.push_back(edge_vector.size());
    }

    if (has_hypernode_weights) {
      if (hypernode_weights == nullptr) {
        LOG(" ****** ignoring hypernode weights ******");
      } else {
        ASSERT(hypernode_weights != nullptr, "Hypergraph has hypernode weights");
        for (HypernodeID i = 0; i < num_hypernodes; ++i) {
          std::getline(file, line);
          std::istringstream line_stream(line);
          HypernodeWeight node_weight;
          line_stream >> node_weight;
          hypernode_weights->push_back(node_weight);
        }
      }
    }
    file.close();
  } else {
    std::cerr << "Error: File not found: " << std::endl;
  }
}

static inline Hypergraph createHypergraphFromFile(const std::string& filename,
                                                  const PartitionID num_parts) {
  HypernodeID num_hypernodes;
  HyperedgeID num_hyperedges;
  HyperedgeIndexVector index_vector;
  HyperedgeVector edge_vector;
  HypernodeWeightVector hypernode_weights;
  HyperedgeWeightVector hyperedge_weights;
  readHypergraphFile(filename, num_hypernodes, num_hyperedges,
                     index_vector, edge_vector, &hyperedge_weights, &hypernode_weights);
  return std::move(Hypergraph(num_hypernodes, num_hyperedges, index_vector, edge_vector,
                              num_parts, &hyperedge_weights, &hypernode_weights));
}


static inline void writeHypernodeWeights(std::ofstream& out_stream, const Hypergraph& hypergraph) {
  for (const HypernodeID hn : hypergraph.nodes()) {
    out_stream << hypergraph.nodeWeight(hn) << std::endl;
  }
}

static inline void writeHGRHeader(std::ofstream& out_stream, const Hypergraph& hypergraph) {
  out_stream << hypergraph.numEdges() << " " << hypergraph.numNodes() << " ";
  if (hypergraph.type() != HypergraphType::Unweighted) {
    out_stream << static_cast<int>(hypergraph.type());
  }
  out_stream << std::endl;
}

static inline void writeHypergraphFile(const Hypergraph& hypergraph, const std::string& filename) {
  ASSERT(!filename.empty(), "No filename for hypergraph file specified");
  std::ofstream out_stream(filename.c_str());
  writeHGRHeader(out_stream, hypergraph);

  for (const HyperedgeID he : hypergraph.edges()) {
    if (hypergraph.type() == HypergraphType::EdgeWeights ||
        hypergraph.type() == HypergraphType::EdgeAndNodeWeights) {
      out_stream << hypergraph.edgeWeight(he) << " ";
    }
    for (const HypernodeID pin : hypergraph.pins(he)) {
      out_stream << pin + 1 << " ";
    }
    out_stream << std::endl;
  }

  if (hypergraph.type() == HypergraphType::NodeWeights ||
      hypergraph.type() == HypergraphType::EdgeAndNodeWeights) {
    writeHypernodeWeights(out_stream, hypergraph);
  }
  out_stream.close();
}

static inline void writeHypergraphForhMetisPartitioning(const Hypergraph& hypergraph,
                                                        const std::string& filename,
                                                        const Mapping& mapping) {
  ASSERT(!filename.empty(), "No filename for hMetis initial partitioning file specified");
  std::ofstream out_stream(filename.c_str());

  // coarse graphs always have edge and node weights, even if graph wasn't coarsend
  out_stream << hypergraph.numEdges() << " " << hypergraph.numNodes() << " ";
  out_stream << static_cast<int>(HypergraphType::EdgeAndNodeWeights);
  out_stream << std::endl;

  for (const HyperedgeID he : hypergraph.edges()) {
    out_stream << hypergraph.edgeWeight(he) << " ";
    for (const HypernodeID pin : hypergraph.pins(he)) {
      ASSERT(mapping.find(pin) != mapping.end(), "No mapping found for pin " << pin);
      out_stream << mapping.find(pin)->second + 1 << " ";
    }
    out_stream << std::endl;
  }

  writeHypernodeWeights(out_stream, hypergraph);
  out_stream.close();
}

static inline void writeHypergraphForPaToHPartitioning(const Hypergraph& hypergraph,
                                                       const std::string& filename,
                                                       const Mapping& mapping) {
  ASSERT(!filename.empty(), "No filename for PaToH initial partitioning file specified");
  std::ofstream out_stream(filename.c_str());
  out_stream << 1;                     // 1-based indexing
  out_stream << " " << hypergraph.numNodes() << " " << hypergraph.numEdges() << " " << hypergraph.numPins();
  out_stream << " " << 3 << std::endl;  // weighting scheme: both edge and node weights

  for (const HyperedgeID he : hypergraph.edges()) {
    out_stream << hypergraph.edgeWeight(he) << " ";
    for (const HypernodeID pin : hypergraph.pins(he)) {
      ASSERT(mapping.find(pin) != mapping.end(), "No mapping found for pin " << pin);
      out_stream << mapping.find(pin)->second + 1 << " ";
    }
    out_stream << std::endl;
  }

  for (const HypernodeID hn : hypergraph.nodes()) {
    out_stream << hypergraph.nodeWeight(hn) << " ";
  }
  out_stream << std::endl;
  out_stream.close();
}

static inline void readPartitionFile(const std::string& filename, std::vector<PartitionID>& partition) {
  ASSERT(!filename.empty(), "No filename for partition file specified");
  ASSERT(partition.empty(), "Partition vector is not empty");
  std::ifstream file(filename);
  if (file) {
    int part;
    while (file >> part) {
      partition.push_back(part);
    }
    file.close();
  } else {
    std::cerr << "Error: File not found: " << std::endl;
  }
}

static inline void writePartitionFile(const Hypergraph& hypergraph, const std::string& filename) {
  ASSERT(!filename.empty(), "No filename for partition file specified");
  std::ofstream out_stream(filename.c_str());
  for (const HypernodeID hn : hypergraph.nodes()) {
    out_stream << hypergraph.partID(hn) << std::endl;
  }
  out_stream.close();
}

static inline void writeHyperedgeVectorFile(const HyperedgeVector& edge_vector, const std::string& filename) {
  ASSERT(!filename.empty(), "No filename for partition file specified");
  std::ofstream out_stream(filename.c_str());
  out_stream << edge_vector.size() << std::endl;
  for(HypernodeID hn : edge_vector) {
      out_stream << hn << std::endl;
  }
  out_stream << std::endl;
  out_stream.close();
}

}  // namespace io

#endif  // SRC_LIB_IO_HYPERGRAPHIO_H_
