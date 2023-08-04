#ifndef PEEL_H_
#define PEEL_H_

#include "common.h"

// class Peel is to peel D-truss on a graph
class Peel final {
 private:
  // adjacency array entry type
  typedef struct final {
    uint32_t vid;
    uint32_t eid;
  } ArrayEntry;
  // data members
  uint32_t n_;  // the # of vertices
  uint32_t m_;  // the # of edges

  //query parameters
  uint32_t m_k_c_;
  uint32_t m_k_f_;

  // the adjacency array representation
  vector<vector<ArrayEntry>> adj_in;
  vector<vector<ArrayEntry>> adj_out;

  vector<uint32_t> m_Sup_c;
  vector<uint32_t> m_Sup_f;
  // the set of edges
  vector<pair<uint32_t, uint32_t>> edges_;
  vector<uint32_t> nodes_;

  static uint32_t h_index(vector<::uint32_t> &neigh_cycle_support);
  void findNeib(vector<ArrayEntry> &vAdj1, vector<ArrayEntry> &vAdj2, vector<pair<uint32_t, uint32_t>> & tris);


 public:
  Peel(vector<pair<uint32_t, uint32_t> > & edges, bool flag, uint32_t iK_c, uint32_t iK_f);
  Peel(const Peel&) = delete;
  Peel& operator=(const Peel&) = delete;

  void start(int iK_c, int iK_f);
  pair<double,double> getDtQuality(vector<pair<uint32_t, uint32_t> > & vEdges);

  // the set of edges in D-truss
  vector<pair<uint32_t, uint32_t> > m_vResE;
};

#endif
