/***************************************************************************
 *  Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "lib/macros.h"
#include "lib/definitions.h"
#include "lib/io/HypergraphIO.h"

using defs::Hypergraph;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cout << "No .hgr file specified" << std::endl;
  }
  std::string hgr_filename(argv[1]);
  std::string graphml_filename(hgr_filename + ".graphml");

  Hypergraph hypergraph(io::createHypergraphFromFile(hgr_filename,2));

  std::ofstream out_stream(graphml_filename.c_str());
  out_stream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  out_stream << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"";
  out_stream << " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"";
  out_stream << " xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns";
  out_stream << " http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">"
            << std::endl;
  out_stream << "<graph id=\"G\" edgedefault=\"undirected\">" << std::endl;
  for (const defs::HypernodeID hn : hypergraph.nodes()){
    out_stream << "<node id=\"n" << hn << "\"/>"   << std::endl;
  }
  for (const defs::HyperedgeID he : hypergraph.edges()) {
    out_stream << "<hyperedge>"  << std::endl;
    for (const defs::HypernodeID pin : hypergraph.pins(he)) {
      out_stream << "<endpoint node=\"n" <<  pin << "\"/>" << std::endl;
    }
    out_stream << "</hyperedge>"  << std::endl;
  }

  out_stream << "</graph>" << std::endl;
  out_stream << "</graphml>" << std::endl;
  out_stream.close();
  std::cout << " ... done!" << std::endl;
  return 0;
}
