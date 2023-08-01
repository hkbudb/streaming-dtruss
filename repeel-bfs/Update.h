#ifndef PEEL_H_
#define PEEL_H_

#include "common.h"

// class Update is to update D-truss on a streaming graph
class Update final {
 private:
  // adjacency array entry type
  typedef struct stAdjInfo {
    uint32_t pid;
    uint32_t eid;
  } AdjEntry;
  // edge array entry type
  typedef struct stEdgeInfo {
    uint32_t cnt;
    bool bInDT;
    bool bUsed;
  } EdgeEntry;
  // data members
  uint32_t m_ui_n;  // the maximum # of vertices
  uint32_t m_ui_m;  // the maximum # of edges
  uint32_t m_ui_kc;  // cycle triangles
  uint32_t m_ui_kf;  // flow triangles

  // the adjacency array representation, size: max pid, capacity: n
  vector<vector<AdjEntry> > m_vv_adj_in;
  vector<vector<AdjEntry> > m_vv_adj_out;

  // the set of nodes
  vector<uint32_t> m_v_rePid;   // size: n
  vector<bool> m_v_bPid;        // size: n
  vector<uint32_t> m_v_nodes;   // size: max pid, capacity: m
  vector<uint32_t> m_v_avaPid;  // size: dynamic, capacity: m

  // the set of edges
  vector<pair<uint32_t, uint32_t> > m_vp_edges; // size: max eid, capacity: m
  vector<uint32_t> m_v_avaEid;  // size: dynamic, capacity: m
  vector<EdgeEntry> m_v_EInfo;  // size: max eid, capacity: m


  void findNeib(vector<AdjEntry> &vAdj1, vector<AdjEntry> &vAdj2, vector<pair<uint32_t, uint32_t> > & vTrdP);
  bool findEid(uint32_t x, uint32_t y, uint32_t *piEid);
  uint32_t renamePid(uint32_t x);
  bool addEInfo(uint32_t x, uint32_t y, uint32_t *piEid);

  uint32_t h_index(vector<::uint32_t> &neigh_cycle_support);
  uint32_t insertOptDT(vector<uint32_t> &vInsE);
  uint32_t insertDT(vector<uint32_t> &vInsE);
  uint32_t rmDT(vector<uint32_t> &vRmE);

 public:
  Update(vector<pair<uint32_t, uint32_t> > & vEdges, uint32_t uiMaxN, uint32_t uiMaxM, uint32_t uiKc, uint32_t uiKf);
  Update(const Update&) = delete;
  Update& operator=(const Update&) = delete;

  uint32_t m_uiDTCnt;
  uint32_t m_uiEdgeCnt;


  void peeling(vector<uint32_t> &vEdges, vector<uint32_t> &vFixE, vector<uint32_t> &vResE);
  void peeling_with_constrcuted(vector<uint32_t> &vEdges, vector<uint32_t> &vFixE, vector<uint32_t> &vResE,
               vector<vector<AdjEntry>> &adj_in, vector<vector<AdjEntry>> &adj_out, vector<bool> &vbTarget,
                                vector<bool> &vbFix);
  void add(vector<pair<uint32_t, uint32_t> > & vEdges);
  void remove(vector<pair<uint32_t, uint32_t> > & vEdges);
  void save(vector<pair<uint32_t, uint32_t> > & vEdges);
  pair<double,double> getDtQuality(vector<pair<uint32_t, uint32_t> > & vEdges);

  // the set of edges in D-truss
  //vector<pair<uint32_t, uint32_t> > m_vResE;
};

#endif
