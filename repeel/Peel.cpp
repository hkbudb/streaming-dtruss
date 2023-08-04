
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <utility>

#include "Peel.h"
/**
 calculate the h-index a vector
**/
uint32_t Peel::h_index(vector<::uint32_t> &neigh_cycle_support) {
    int n = neigh_cycle_support.size();
    vector <int> bucket(n + 1);
    for(int i = 0; i < n; i++){
        int x = neigh_cycle_support[i];
        if(x >= n){
            bucket[n]++;
        } else {
            bucket[x]++;
        }
    }
    int cnt = 0;
    for(int i = n; i >= 0; i--){
        cnt += bucket[i];
        if(cnt >= i)return i;
    } return -1;
}
/**
    find neighbor
**/
void Peel::findNeib(vector<ArrayEntry> &vAdj1, vector<ArrayEntry> &vAdj2, vector<pair<uint32_t, uint32_t>> & tris)
{
  size_t p1 = 0, p2 = 0;
  while (p1 < vAdj1.size() && p2 < vAdj2.size()) {
    if (vAdj1[p1].vid == vAdj2[p2].vid) {
      tris.push_back({vAdj1[p1].eid, vAdj2[p2].eid});
      ++p1; ++p2;
    } else if (vAdj1[p1].vid < vAdj2[p2].vid) {
      ++p1;
    } else {
      ++p2;
    }
  }
}
/**
    init
**/
Peel::Peel(vector<pair<uint32_t, uint32_t> > & vEdges, bool flag, uint32_t iK_c, uint32_t iK_f)
{
    m_k_c_ = iK_c;
    m_k_f_ = iK_f;
    //do not use that h-index pruning
    if(!flag){
        edges_ = vEdges;
        sort(edges_.begin(), edges_.end());
        auto it = unique(edges_.begin(), edges_.end());
        edges_.erase(it, edges_.end());

        m_ = edges_.size();
        nodes_.reserve(2 * m_);
        map<uint32_t, vector<pair<uint32_t, uint32_t> > > mpRec;
        map<uint32_t, vector<pair<uint32_t, uint32_t> > >::iterator itmpRec;
        vector<pair<uint32_t, uint32_t> >::iterator itvE;
        for (uint32_t eid = 0; eid < m_; ++eid)
        {
            const uint32_t x = edges_[eid].first;
            const uint32_t y = edges_[eid].second;

            //printf("PEEL get (%d, %d)\n", x, y);

            mpRec[x].push_back(pair<uint32_t, uint32_t>(eid, 1));
            mpRec[y].push_back(pair<uint32_t, uint32_t>(eid, 2));
        }
        int iRePid = 0;
        for (itmpRec = mpRec.begin(); itmpRec != mpRec.end(); ++itmpRec, ++iRePid)
        {
            nodes_.push_back(itmpRec->first);
            for (itvE = itmpRec->second.begin(); itvE != itmpRec->second.end(); ++itvE)
            {
                int iEid = itvE->first;
                if (1 == itvE->second)
                {
                    edges_[iEid].first = iRePid;
                }
                else if (2 == itvE->second)
                {
                    edges_[iEid].second = iRePid;
                }
            }
        }
        n_ = nodes_.size();
    }
    //with h-index pruning, calculate the cycle_h_index / flow_h_index to prune first
    else{
        vector<pair<uint32_t, uint32_t>> tmp_edges_;
        vector<uint32_t> tmp_nodes_;
        tmp_edges_ = vEdges;
        sort(tmp_edges_.begin(), tmp_edges_.end());
        auto tmp_it = unique(tmp_edges_.begin(), tmp_edges_.end());
        tmp_edges_.erase(tmp_it, tmp_edges_.end());

        uint32_t tmp_n_ = 0;
        uint32_t tmp_m_ = 0;
        tmp_m_ = tmp_edges_.size();
        tmp_nodes_.reserve(2 * tmp_m_);
        map<uint32_t, vector<pair<uint32_t, uint32_t> > > tmp_mpRec;
        map<uint32_t, vector<pair<uint32_t, uint32_t> > >::iterator tmp_itmpRec;
        vector<pair<uint32_t, uint32_t> >::iterator tmp_itvE;
        for (uint32_t eid = 0; eid < tmp_m_; ++eid)
        {
            const uint32_t x = tmp_edges_[eid].first;
            const uint32_t y = tmp_edges_[eid].second;

            //printf("PEEL get (%d, %d)\n", x, y);

            tmp_mpRec[x].push_back(pair<uint32_t, uint32_t>(eid, 1));
            tmp_mpRec[y].push_back(pair<uint32_t, uint32_t>(eid, 2));
        }
        int tmp_iRePid = 0;
        for (tmp_itmpRec = tmp_mpRec.begin(); tmp_itmpRec != tmp_mpRec.end(); ++tmp_itmpRec, ++tmp_iRePid)
        {
            tmp_nodes_.push_back(tmp_itmpRec->first);
            for (tmp_itvE = tmp_itmpRec->second.begin(); tmp_itvE != tmp_itmpRec->second.end(); ++tmp_itvE)
            {
                int iEid = tmp_itvE->first;
                if (1 == tmp_itvE->second)
                {
                    tmp_edges_[iEid].first = tmp_iRePid;
                }
                else if (2 == tmp_itvE->second)
                {
                    tmp_edges_[iEid].second = tmp_iRePid;
                }
            }
        }
        tmp_n_ = tmp_nodes_.size();

        //pruning based on h-index in the pruned edges
        vector<pair<uint32_t, uint32_t> > pruned_targe_edges_support_vec;
        /*count unique vertex number*/

        /*constrct the adj_lists first*/

        vector<vector<ArrayEntry> > tmp_adj_in(tmp_n_);
        vector<vector<ArrayEntry> > tmp_adj_out(tmp_n_);
        for (uint32_t eid = 0; eid < tmp_m_; ++eid) {
            const uint32_t v1 = tmp_edges_[eid].first;
            const uint32_t v2 = tmp_edges_[eid].second;
            tmp_adj_out[v1].push_back({v2, eid});
            tmp_adj_in[v2].push_back({v1, eid});
        }

        for (uint32_t vid = 0; vid < tmp_n_; ++vid) {
            std::sort(tmp_adj_in[vid].begin(), tmp_adj_in[vid].end(),
                      [](const ArrayEntry& ae1, const ArrayEntry& ae2) {
                          return ae1.vid < ae2.vid;
                      });
            std::sort(tmp_adj_out[vid].begin(), tmp_adj_out[vid].end(),
                      [](const ArrayEntry& ae1, const ArrayEntry& ae2) {
                          return ae1.vid < ae2.vid;
                      });
        }
        /*calculate the edge support and cycle/flow neighbors in the inducted graph*/
        vector<::uint32_t> sub_cycle_support(tmp_m_, 0);
        vector<::uint32_t> sub_flow_support(tmp_m_, 0);
        vector<vector<pair<::uint32_t,::uint32_t>>> eid2_cycle_neighbor(tmp_m_), eid2_flow_neighbor(tmp_m_);
        vector<pair<uint32_t, uint32_t> > vCE;
        vector<pair<uint32_t, uint32_t> > vFE;
        pruned_targe_edges_support_vec.clear();
        for (uint32_t eid = 0; eid < tmp_m_; ++eid) {
            ::uint32_t x = tmp_edges_[eid].first;
            ::uint32_t y = tmp_edges_[eid].second;
            vCE.clear();
            findNeib(tmp_adj_in[x], tmp_adj_out[y], vCE);
            vFE.clear();
            findNeib(tmp_adj_in[x], tmp_adj_in[y], vFE);
            findNeib(tmp_adj_out[x], tmp_adj_out[y], vFE);
            findNeib(tmp_adj_out[x], tmp_adj_in[y], vFE);
            sub_cycle_support[eid] = vCE.size();
            sub_flow_support[eid] = vFE.size();
            eid2_cycle_neighbor[eid] = vCE;
            eid2_flow_neighbor[eid] = vFE;
        }
        /*calculte the h-index*/
        for(unsigned int eid = 0; eid < tmp_m_; eid++){
            ::uint32_t cycle_neighbor_h_index = 0;
            vector<::uint32_t> neighbor_edge_support(eid2_cycle_neighbor[eid].size(),0);
            if(m_k_c_ > 1){
                //tranverse the cycle/flow neighbors of eid to get its cycle support vector
                for(::uint32_t j = 0; j < eid2_cycle_neighbor[eid].size(); j++){
                    pair<::uint32_t,::uint32_t> edges = eid2_cycle_neighbor[eid][j];
                    neighbor_edge_support[j] = std::min(sub_cycle_support[edges.first], sub_cycle_support[edges.second]);
                }
                cycle_neighbor_h_index = h_index(neighbor_edge_support);
            }
            else{
                cycle_neighbor_h_index = sub_cycle_support[eid];
            }

            ::uint32_t flow_neighbor_h_index = 0;
            if(m_k_f_ > 1){
                neighbor_edge_support.clear();
                neighbor_edge_support.resize(eid2_flow_neighbor[eid].size(),0);
                for(::uint32_t j = 0; j < eid2_flow_neighbor[eid].size(); j++){
                    pair<::uint32_t,::uint32_t> edges = eid2_flow_neighbor[eid][j];
                    neighbor_edge_support[j] = std::min(sub_flow_support[edges.first], sub_flow_support[edges.second]);
                }
                flow_neighbor_h_index = h_index(neighbor_edge_support);
            }
            else{
                flow_neighbor_h_index = sub_flow_support[eid];
            }
            ASSERT_MSG(cycle_neighbor_h_index <= sub_cycle_support[eid] && flow_neighbor_h_index <= sub_flow_support[eid], "h-index error");

            if((cycle_neighbor_h_index >= m_k_c_ && flow_neighbor_h_index >= m_k_f_)){
                pruned_targe_edges_support_vec.push_back(tmp_edges_[eid]);
            }
        }
        //do the initialization
        //cout << "pruned_targe_edges_support_vec.size() = " << pruned_targe_edges_support_vec.size()
        //     << "vEsize: " << vEdges.size() << endl;
        ASSERT_MSG(pruned_targe_edges_support_vec.size() <= vEdges.size(), "h-index pruning size error");
        edges_ = pruned_targe_edges_support_vec;
        sort(edges_.begin(), edges_.end());
        auto it = unique(edges_.begin(), edges_.end());
        edges_.erase(it, edges_.end());

        m_ = edges_.size();
        nodes_.reserve(2 * m_);
        map<uint32_t, vector<pair<uint32_t, uint32_t> > > mpRec;
        map<uint32_t, vector<pair<uint32_t, uint32_t> > >::iterator itmpRec;
        vector<pair<uint32_t, uint32_t> >::iterator itvE;
        for (uint32_t eid = 0; eid < m_; ++eid)
        {
            const uint32_t x = edges_[eid].first;
            const uint32_t y = edges_[eid].second;

            //printf("PEEL get (%d, %d)\n", x, y);

            mpRec[x].push_back(pair<uint32_t, uint32_t>(eid, 1));
            mpRec[y].push_back(pair<uint32_t, uint32_t>(eid, 2));
        }
        int iRePid = 0;
        for (itmpRec = mpRec.begin(); itmpRec != mpRec.end(); ++itmpRec, ++iRePid)
        {
            nodes_.push_back(itmpRec->first);
            for (itvE = itmpRec->second.begin(); itvE != itmpRec->second.end(); ++itvE)
            {
                int iEid = itvE->first;
                if (1 == itvE->second)
                {
                    edges_[iEid].first = iRePid;
                }
                else if (2 == itvE->second)
                {
                    edges_[iEid].second = iRePid;
                }
            }
        }
        n_ = nodes_.size();
    }


    //printf("PEEL init n: %d, m: %d\n", n_, m_);
}
/**
    calculate the CMSin and CMSout of finnal D-truss
**/
pair<double,double> Peel::getDtQuality(vector<pair<uint32_t, uint32_t>> & vEdges)
{
    
    double CMSin = 0, CMSout =  0;
    set<uint32_t> VertexSet;
    map<uint32_t,set<uint32_t>> OutGraph, InGraph;
    
    for(const auto i : vEdges){  //obtain all the vertices in the community
        uint32_t x = i.first;
        uint32_t y = i.second;

        VertexSet.insert(x);
        VertexSet.insert(y);



        if(OutGraph[x].empty()){
            set<uint32_t> tmp;
            tmp.insert(y);
            OutGraph[x] = tmp;
        }
        else{
            OutGraph[x].insert(y);
        }

        if(InGraph[y].empty()){
            set<uint32_t> tmp;
            tmp.insert(x);
            InGraph[y] = tmp;
        }
        else{
            InGraph[y].insert(x);
        }
    }
   

    int VertexNum = VertexSet.size();

    for(auto i = VertexSet.cbegin(); i != VertexSet.cend(); i++){
        for(auto j = VertexSet.cbegin();j != VertexSet.cend();j++){
            double InterIn = 0, UnionIn = 0, InterOut = 0, UnionOut = 0;

            vector<uint32_t> tmpResult1;
            set_intersection(InGraph[*i].begin(),InGraph[*i].end(),InGraph[*j].begin(),InGraph[*j].end(),
                                inserter(tmpResult1, tmpResult1.begin()));
            InterIn = (double) tmpResult1.size();
            
           
            vector<uint32_t> tmpResult2;
            set_union(InGraph[*i].begin(),InGraph[*i].end(),InGraph[*j].begin(),InGraph[*j].end(),
                        std::back_inserter(tmpResult2));
            UnionIn = (double) tmpResult2.size();

            
            vector<uint32_t> tmpResult3;
            set_intersection(OutGraph[*i].begin(),OutGraph[*i].end(),OutGraph[*j].begin(),OutGraph[*j].end(),
                            inserter(tmpResult3, tmpResult3.begin()));
            InterOut = (double) tmpResult3.size();
  

            vector<uint32_t> tmpResult4;
            set_union(OutGraph[*i].begin(),OutGraph[*i].end(),OutGraph[*j].begin(),OutGraph[*j].end(),
                        std::back_inserter(tmpResult4));
            UnionOut = (double) tmpResult4.size();
         
        
            if(UnionOut != 0){
                CMSout += (double) InterOut / UnionOut;
            }
            if(UnionIn != 0){
                CMSin += (double) InterIn / UnionIn;
            }
            
        }
        
    }
    CMSin = CMSin / VertexNum / VertexNum;
    CMSout = CMSout / VertexNum / VertexNum;
    return {CMSin,CMSout};
}
/**
    get D-truss by peeling
**/
void Peel::start(int iK_c, int iK_f)
{
    ASSERT_MSG(iK_c == m_k_c_, "iK_c != m_k_c_");
    ASSERT_MSG(iK_f == m_k_f_, "iK_f != m_k_f_");
  const auto beg = std::chrono::steady_clock::now();
  // initialize adjacency arrays
  adj_in.resize(n_);
  adj_out.resize(n_);
  for (uint32_t eid = 0; eid < m_; ++eid) {
    const uint32_t v1 = edges_[eid].first;
    const uint32_t v2 = edges_[eid].second;
    adj_out[v1].push_back({v2, eid});
    adj_in[v2].push_back({v1, eid});
  }
  for (uint32_t vid = 0; vid < n_; ++vid) {
    std::sort(adj_in[vid].begin(), adj_in[vid].end(),
              [](const ArrayEntry& ae1, const ArrayEntry& ae2) {
                return ae1.vid < ae2.vid;
              });
    std::sort(adj_out[vid].begin(), adj_out[vid].end(),
              [](const ArrayEntry& ae1, const ArrayEntry& ae2) {
                return ae1.vid < ae2.vid;
              });
  }
  // peeling for D-truss
  // 1. compute the support of each edge by triangle listing
  // 1.1. define a total order over the vertices
  const auto pred = [this](const uint32_t v1, const uint32_t v2) {
    const size_t deg1 = adj_in[v1].size() + adj_out[v1].size();
    const size_t deg2 = adj_in[v2].size() + adj_out[v2].size();
    if (deg1 != deg2) return deg1 > deg2;
    else return v1 > v2;
  };
  // 1.2. sort the vertices in non-ascending order of degree
  vector<uint32_t> verts(n_);
  std::iota(verts.begin(), verts.end(), 0);
  std::sort(verts.begin(), verts.end(), pred);
  // 1.3. call the "forward" algorithm to list cycle triangles
  const auto begCnt = std::chrono::steady_clock::now();
  m_Sup_c.resize(m_, 0);
  vector<vector<ArrayEntry>> A_in(n_);
  vector<vector<ArrayEntry>> A_out(n_);
  for (const uint32_t v : verts) {
    for (const auto ae : adj_in[v]) {
      const uint32_t u = ae.vid;
      const uint32_t e = ae.eid; 
      if (pred(u, v)) continue;
      size_t pv = 0, pu = 0;
      while (pv < A_out[v].size() && pu < A_in[u].size()) {
        if (A_out[v][pv].vid == A_in[u][pu].vid) {
          ++m_Sup_c[e]; ++m_Sup_c[A_out[v][pv].eid]; ++m_Sup_c[A_in[u][pu].eid];
          ++pv; ++pu;
        } else if (pred(A_out[v][pv].vid, A_in[u][pu].vid)) {
          ++pv;
        } else {
          ++pu;
        }
      }
      A_out[u].push_back({v, e});
    }
    for (const auto ae : adj_out[v]) {
      const uint32_t u = ae.vid;
      const uint32_t e = ae.eid;
      if (pred(u, v)) continue;
      size_t pv = 0, pu = 0;
      while (pv < A_in[v].size() && pu < A_out[u].size()) {
        if (A_in[v][pv].vid == A_out[u][pu].vid) {
          ++m_Sup_c[e]; ++m_Sup_c[A_in[v][pv].eid]; ++m_Sup_c[A_out[u][pu].eid];
          ++pv; ++pu;
        } else if (pred(A_in[v][pv].vid, A_out[u][pu].vid)) {
          ++pv;
        } else {
          ++pu;
        }
      }
      A_in[u].push_back({v, e});
    }
  }
  decltype(A_in)().swap(A_in);
  decltype(A_out)().swap(A_out);
  decltype(verts)().swap(verts);

  // 1.4. list flow triangles
  m_Sup_f.resize(m_, 0);
  for (uint32_t eid = 0; eid < m_; ++eid) {
    const uint32_t v1 = edges_[eid].first;
    const uint32_t v2 = edges_[eid].second;
    vector<pair<uint32_t, uint32_t>> tris;
    findNeib(adj_out[v1], adj_in[v2], tris);
    for (const auto tri : tris)
    {
        const uint32_t e1 = tri.first;
        const uint32_t e2 = tri.second;
        ++m_Sup_f[eid];++m_Sup_f[e1];++m_Sup_f[e2];
    }
  }

  const auto endCnt = std::chrono::steady_clock::now();
  const auto difCnt = endCnt - begCnt;
  /*printf("Counting costs \x1b[1;31m%f\x1b[0m ms.\n",
         std::chrono::duration<double, std::milli>(difCnt).count());*/

    // 2. peeling
    vector<bool> removed(m_, false);
    vector<bool> vWaitRmFlag(m_, false);
    uint32_t uiRmCnt = 0;
    vector<uint32_t> vWait(m_);
    vector<uint32_t> vWaitRmCache(m_);
    // init cache
    vWaitRmCache.clear();
    for (uint32_t eid = 0; eid < m_; ++eid)
    {
      /*printf("(%d, %d) c_sup: %d f_sup: %d bool: %d %d\n", edges_[eid].first, edges_[eid].second,
             m_Sup_c[eid], m_Sup_f[eid],
             m_Sup_c[eid] < iK_c, m_Sup_f[eid] < iK_f);*/
      if (m_Sup_c[eid] < iK_c)
      {
        vWaitRmFlag[eid] = true;
        vWaitRmCache.push_back(eid);
      }
      else if (m_Sup_f[eid] < iK_f)
      {
        vWaitRmFlag[eid] = true;
        vWaitRmCache.push_back(eid);
      }
    }
    while (uiRmCnt < m_)
    {
        vWait.clear();
        if (vWaitRmCache.empty())
        {
            /* all remained are D-truss */
            break;
        }
        else
        {
            vWait.swap(vWaitRmCache);
        }

        for (auto eid : vWait)
        {
            removed[eid] = true;
            ++uiRmCnt;
            //printf("start find triangles\n");
            vector<pair<uint32_t, uint32_t>> tris;
          const uint32_t v1 = edges_[eid].first;
          const uint32_t v2 = edges_[eid].second;
            // find triangles containing the edge with ID eid, flow type 1: v1 -> v2
            findNeib(adj_out[v1], adj_in[v2], tris);
            // find triangles containing the edge with ID eid, flow type 2: v1 -> w
            findNeib(adj_out[v1], adj_out[v2], tris);
            // find triangles containing the edge with ID eid, flow type 3: w -> v1
            findNeib(adj_in[v1], adj_in[v2], tris);
            for (const auto tri : tris)
            {
                const uint32_t e1 = tri.first;
                const uint32_t e2 = tri.second;
                if (removed[e1] || removed[e2]) continue;
                if (!vWaitRmFlag[e1])
                {
                    --m_Sup_f[e1];
                    if (m_Sup_f[e1] < iK_f)
                  {
                    vWaitRmFlag[e1] = true;
                    vWaitRmCache.push_back(e1);
                  }
                }
                if (!vWaitRmFlag[e2])
                {
                    --m_Sup_f[e2];
                    if (m_Sup_f[e2] < iK_f)
                  {
                    vWaitRmFlag[e2] = true;
                    vWaitRmCache.push_back(e2);
                  }
                }
            }
            // find triangles containing the edge with ID eid, cycle
            tris.clear();
            findNeib(adj_in[v1], adj_out[v2], tris);
            for (const auto tri : tris)
            {
                const uint32_t e1 = tri.first;
                const uint32_t e2 = tri.second;
                if (removed[e1] || removed[e2]) continue;
                if (!vWaitRmFlag[e1])
                {
                    --m_Sup_c[e1];
                    if (m_Sup_c[e1] < iK_c)
                  {
                    vWaitRmFlag[e1] = true;
                    vWaitRmCache.push_back(e1);
                  }
                }
                if (!vWaitRmFlag[e2])
                {
                    --m_Sup_c[e2];
                    if (m_Sup_c[e2] < iK_c)
                  {
                    vWaitRmFlag[e2] = true;
                    vWaitRmCache.push_back(e2);
                  }
                }
            }
        }

    }

  // clear adj_
  decltype(adj_in)().swap(adj_in);
  decltype(adj_out)().swap(adj_out);

  // keep result
  for (uint32_t eid = 0; eid < m_; ++eid)
  {
      if (!removed[eid])
      {
        uint32_t x = edges_[eid].first;
        uint32_t y = edges_[eid].second;
        m_vResE.push_back(pair<uint32_t, uint32_t>(nodes_[x], nodes_[y]));
      }
  }

  const auto end = std::chrono::steady_clock::now();
  const auto dif = end - beg;
  /*printf("Peeling costs \x1b[1;31m%f\x1b[0m ms.\n",
         std::chrono::duration<double, std::milli>(dif).count());*/
    //printf("DEBUG index done last: (%d, %d)\n", edges_[m_ - 1].first, edges_[m_ - 1].second);
}
