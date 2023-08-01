
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <utility>

#include "Update.h"

::uint32_t total_pruned_size = 0, total_unpruned_size = 0;
chrono::duration<double> peeling_time;
chrono::duration<double> g_tUptIns;
chrono::duration<double> g_tUptRm;
int SearchSpace = 0;
/**
    find neighbor
**/
void Update::findNeib(vector<AdjEntry> &vAdj1, vector<AdjEntry> &vAdj2, vector<pair<uint32_t, uint32_t> > & vTrdE)
{
  size_t p1 = 0, p2 = 0;
  while (p1 < vAdj1.size() && p2 < vAdj2.size()) {
    if (vAdj1[p1].pid == vAdj2[p2].pid) {
      vTrdE.push_back(pair<uint32_t, uint32_t>(vAdj1[p1].eid, vAdj2[p2].eid));
      ++p1; ++p2;
    } else if (vAdj1[p1].pid < vAdj2[p2].pid) {
      ++p1;
    } else {
      ++p2;
    }
  }
}
/**
 calculate the h-index a vector
**/
uint32_t Update::h_index(vector<::uint32_t> &neigh_cycle_support) {
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
    init
**/
Update::Update(vector<pair<uint32_t, uint32_t> > & vEdges, uint32_t uiMaxN, uint32_t uiMaxM, uint32_t uiKc, uint32_t uiKf)
{
    m_ui_n = uiMaxN;
    m_ui_m = uiMaxM;
    m_ui_kc = uiKc;
    m_ui_kf = uiKf;

    m_vv_adj_in.reserve(m_ui_n);
    m_vv_adj_out.reserve(m_ui_n);
    m_v_rePid.resize(m_ui_n);
    m_v_bPid.resize(m_ui_n);
    m_v_nodes.reserve(m_ui_n);
    m_v_avaPid.reserve(m_ui_n);

    m_vp_edges.reserve(m_ui_m);
    m_v_avaEid.reserve(m_ui_m);
    m_v_EInfo.reserve(m_ui_m);

    vector<uint32_t> vTargetE;
    vTargetE.reserve(vEdges.size());
    uint32_t eid = 0;
    for (auto atE : vEdges)
    {
        bool bRes = addEInfo(atE.first, atE.second, &eid);
        /*printf("NEW get (%d, %d) rename (%d, %d) bool: %d\n",
               atE.first, atE.second,
               m_vp_edges[eid].first, m_vp_edges[eid].second,
               bRes);*/
        if (bRes)
        {
            vTargetE.push_back(eid);
        }
        /*else
        {
            printf("NEW in it (%d, %d)\n",
                   m_v_nodes[m_vp_edges[eid].first], m_v_nodes[m_vp_edges[eid].second]);
        }*/
    }

    vector<uint32_t> vFixE;
    vector<uint32_t> vResE;

    peeling(vTargetE, vFixE, vResE);  //tmp is not used
    for (auto atE : vResE)
    {
        m_v_EInfo[atE].bInDT = true;
        //printf("D-truss eid: %d\n", atE);
    }
    m_uiDTCnt = vResE.size();
    printf("Init D-truss size: %d\n", m_uiDTCnt);

}
/**
    binary search eid
    if it is a new edge, return the position to be inserted
**/
bool Update::findEid(uint32_t x, uint32_t y, uint32_t *piEid)
{
    if (m_vv_adj_out[x].empty())
    {
        return false;
    }
    if (m_vv_adj_in[y].empty())
    {
        return false;
    }
    /* the minimum degree */
    vector<AdjEntry> *pv = &(m_vv_adj_out[x]);
    uint32_t v = y;
    if (m_vv_adj_out[x].size() > m_vv_adj_in[y].size())
    {
        pv = &(m_vv_adj_in[y]);
        v = x;
    }

    vector<AdjEntry>::iterator it = lower_bound(pv->begin(), pv->end(), v,
              [](const AdjEntry& e, uint32_t tpPid) {
                return e.pid < tpPid;
              });
    if ((it != pv->end()) && (it->pid == v))
    {
        *piEid = it->eid;
        /*printf("(%d, %d) v: %d x size: %d y size: %d pv size: %d e: (%d, %d) eid: %d\n", x, y, v,
               m_vv_adj_out[x].size(),
               m_vv_adj_in[y].size(), pv->size(),
               m_vp_edges[it->eid].first, m_vp_edges[it->eid].second, it->eid);
        for (auto atE : m_vv_adj_in[y])
        {
            printf("(%d, %d) ", atE.pid, atE.eid);
        }
        printf("\n");*/
        return true;
    }
    return false;
}

/**
    rename pid
**/
uint32_t Update::renamePid(uint32_t x)
{
    uint32_t u = 0;
    if (m_v_bPid[x])
    {
        u = m_v_rePid[x];
    }
    else if (!m_v_avaPid.empty())
    {
        u = m_v_avaPid.back();
        m_v_avaPid.pop_back();
        m_v_bPid[x] = true;
        m_v_rePid[x] = u;
        m_v_nodes[u] = x;
    }
    else
    {
        u = m_v_nodes.size();
        m_v_nodes.push_back(x);
        m_v_bPid[x] = true;
        m_v_rePid[x] = u;
        m_vv_adj_in.resize(u + 1);
        m_vv_adj_out.resize(u + 1);
    }
    return u;
}
/**
    add one edge info
    return eid
**/
bool Update::addEInfo(uint32_t x, uint32_t y, uint32_t *puiEid)
{
    /* rename */
    uint32_t u = renamePid(x);
    uint32_t v = renamePid(y);
    if (findEid(u, v, puiEid))
    {
        /* not new */
        ++m_v_EInfo[*puiEid].cnt;
        return false;
    }
    /* new */
    if (!m_v_avaEid.empty())
    {
        *puiEid = m_v_avaEid.back();
        m_v_avaEid.pop_back();
        m_vp_edges[*puiEid] = pair<uint32_t, uint32_t>(u, v);
    }
    else
    {
        *puiEid = m_vp_edges.size();
        m_vp_edges.push_back(pair<uint32_t, uint32_t>(u, v));
        EdgeEntry stEInfo = {0};
        m_v_EInfo.push_back(stEInfo);
    }
    m_v_EInfo[*puiEid].cnt = 1;
    m_v_EInfo[*puiEid].bInDT = false;
    m_v_EInfo[*puiEid].bUsed = true;

    /*printf("EDGE_INFO (%d, %d) edge (%d, %d) original: (%d, %d) eid: %d\n",
           u, v, m_vp_edges[*puiEid].first, m_vp_edges[*puiEid].second,
           m_v_nodes[u], m_v_nodes[v], *puiEid);*/

    ASSERT(u == m_vp_edges[*puiEid].first);
    ASSERT(v == m_vp_edges[*puiEid].second);

    /* adj info */
    vector<AdjEntry>::iterator it = upper_bound(m_vv_adj_out[u].begin(), m_vv_adj_out[u].end(), v,
              [](uint32_t tpPid, const AdjEntry& e) {
                return tpPid < e.pid;
              });
    m_vv_adj_out[u].insert(it, (AdjEntry){v, *puiEid});
    it = upper_bound(m_vv_adj_in[v].begin(),m_vv_adj_in[v].end(), u,
              [](uint32_t tpPid, const AdjEntry& e) {
                return tpPid < e.pid;
              });
    m_vv_adj_in[v].insert(it, (AdjEntry){u, *puiEid});

    return true;
}
/**
  add one edge info, with roughly calculated skyline trussness pruning technique
**/
uint32_t Update::insertOptDT(vector<uint32_t> &vInsE)
{
    list<uint32_t> lsQ;
    vector<pair<uint32_t, uint32_t> > vCE;
    vector<pair<uint32_t, uint32_t> > vFE;

    for (uint32_t eid : vInsE)
    {
        uint32_t x = m_vp_edges[eid].first;
        uint32_t y = m_vp_edges[eid].second;
        /* cycle support */
        vCE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_out[y], vCE);
        uint32_t uiCSup = vCE.size();
        //printf("BFS queue get CSup: %d\n", uiCSup);
        if (uiCSup < m_ui_kc)
        {
            /* cannot in D-truss */
            continue;
        }
        /* flow support */
        vFE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_in[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_out[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_in[y], vFE);
        uint32_t uiFSup = vFE.size();
        //printf("BFS queue get FSup: %d\n", uiFSup);
        if (uiFSup < m_ui_kf)
        {
            /* cannot in D-truss */
            continue;
        }

        lsQ.push_back(eid);
    }

    if (lsQ.empty())
    {
        /* nothing to do */
        return 0;
    }

    uint32_t m_ = m_vp_edges.size();


//    uint32_t iteration_num = 1;
//    vector<vector<pair<uint32_t, uint32_t> > > eid2_cycle_neighbor(m_);
//    vector<::uint32_t> cycle_support(m_, 0);
//    vector<bool> cycle_support_calculated(m_, false);
//
//    vector<vector<pair<uint32_t, uint32_t> > > eid2_flow_neighbor(m_);
//    vector<::uint32_t> flow_support(m_, 0);
//    vector<bool> flow_support_calculated(m_, false);

    vector<uint32_t> vTargetE;
    vTargetE.reserve(m_);
    vector<uint32_t> vDetE;
    vDetE.reserve(m_);
    vector<bool> vbInQ(m_, false);
    vector<bool> vbInPool(m_, false);
    for (uint32_t eid : lsQ)
    {
        vbInQ[eid] = true;
    }

    while (!lsQ.empty())
    {
        uint32_t uiCurEid = lsQ.front();
        lsQ.pop_front();
        uint32_t x = m_vp_edges[uiCurEid].first;
        uint32_t y = m_vp_edges[uiCurEid].second;

        //printf("BFS queue get original (%d, %d)\n", m_v_nodes[x], m_v_nodes[y]);
        //printf("BFS queue get (%d, %d) in it: %d\n", x, y, m_v_EInfo[uiCurEid].bInDT);

        if (m_v_EInfo[uiCurEid].bInDT)
        {
            vDetE.push_back(uiCurEid);
            vbInPool[uiCurEid] = true;
            continue;
        }
        /* check support */
        /* cycle support */
        vCE.clear();
        uint32_t uiCSup = 0;
//        if(!cycle_support_calculated[uiCurEid]){
        findNeib(m_vv_adj_in[x], m_vv_adj_out[y], vCE);
        uiCSup = vCE.size();
//            cycle_support[uiCurEid] = uiCSup;
//            eid2_cycle_neighbor[uiCurEid] = vCE;
//            cycle_support_calculated[uiCurEid] = true;
//        }
//        else{
//            uiCSup = cycle_support[uiCurEid];
//            vCE = eid2_cycle_neighbor[uiCurEid];
//        }



        //printf("BFS queue get CSup: %d\n", uiCSup);
        if (uiCSup < m_ui_kc)
        {
            /* cannot in D-truss */
            continue;
        }

        /* flow support */
        vFE.clear();
        uint32_t uiFSup = 0;
//        if(!flow_support_calculated[uiCurEid]){
        findNeib(m_vv_adj_in[x], m_vv_adj_in[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_out[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_in[y], vFE);
        uiFSup = vFE.size();
//            flow_support[uiCurEid] = uiFSup;
//            eid2_flow_neighbor[uiCurEid] = vFE;
//            flow_support_calculated[uiCurEid] = true;
//        }
//        else{
//            uiFSup = flow_support[uiCurEid];
//            vFE = eid2_flow_neighbor[uiCurEid];
//        }


        //printf("BFS queue get FSup: %d\n", uiFSup);
        if (uiFSup < m_ui_kf)
        {
            /* cannot in D-truss */
            continue;
        }

//        //traverse the neighbors to get its cycle-h-index / flow-h-index to see if larger than kc/kf
//        vector<::uint32_t> neighbor_edge_support(cycle_support[uiCurEid],0);
//        #pragma omp parallel for num_threads(1)
//        for(::uint32_t j = 0; j < vCE.size(); j++){
//            //handling cycle support neighbors
//            ::uint32_t neighbor_support_1 = 0, neighbor_support_2 = 0;
//            if(cycle_support_calculated[vCE[j].first]) {
//                neighbor_support_1 = cycle_support[vCE[j].first];
//            }
//            else{
//                vector<pair<uint32_t, uint32_t> > tmp_vCE;
//                tmp_vCE.reserve(m_);
//                findNeib(m_vv_adj_in[m_vp_edges[vCE[j].first].first], m_vv_adj_out[m_vp_edges[vCE[j].first].second], tmp_vCE);
//                neighbor_support_1 = tmp_vCE.size();
//                eid2_cycle_neighbor[vCE[j].first] = tmp_vCE;
//                cycle_support[vCE[j].first] = neighbor_support_1;
//                cycle_support_calculated[vCE[j].first] = true;
//            }
//
//            if(cycle_support_calculated[vCE[j].second]) {
//                neighbor_support_2 = cycle_support[vCE[j].second];
//            }
//            else{
//                vector<pair<uint32_t, uint32_t> > tmp_vCE;
//                tmp_vCE.reserve(m_);
//                findNeib(m_vv_adj_in[m_vp_edges[vCE[j].second].first], m_vv_adj_out[m_vp_edges[vCE[j].second].second], tmp_vCE);
//                neighbor_support_2 = tmp_vCE.size();
//                eid2_cycle_neighbor[vCE[j].second] = tmp_vCE;
//                cycle_support[vCE[j].second] = neighbor_support_2;
//                cycle_support_calculated[vCE[j].second] = true;
//            }
//            neighbor_edge_support[j] = std::min(neighbor_support_1, neighbor_support_2);
//        }
//        ::uint32_t cycle_neighbor_h_index = h_index(neighbor_edge_support);
//
////        neighbor_edge_support.clear();
////        neighbor_edge_support.resize(flow_support[uiCurEid],0);
////        #pragma omp parallel for num_threads(1)
////        for(::uint32_t j = 0; j < vFE.size(); j++){
////            //handling flow support neighbors
////            ::uint32_t neighbor_support_1 = 0, neighbor_support_2 = 0;
////            if(flow_support_calculated[vFE[j].first]) {
////                neighbor_support_1 = flow_support[vFE[j].first];
////            }
////            else{
////                vector<pair<uint32_t, uint32_t> > tmp_vFE;
////                tmp_vFE.reserve(m_);
////                findNeib(m_vv_adj_in[m_vp_edges[vFE[j].first].first], m_vv_adj_in[m_vp_edges[vFE[j].first].second], tmp_vFE);
////                findNeib(m_vv_adj_out[m_vp_edges[vFE[j].first].first], m_vv_adj_out[m_vp_edges[vFE[j].first].second], tmp_vFE);
////                findNeib(m_vv_adj_out[m_vp_edges[vFE[j].first].first], m_vv_adj_in[m_vp_edges[vFE[j].first].second], tmp_vFE);
////                neighbor_support_1 = tmp_vFE.size();
////                eid2_flow_neighbor[vFE[j].first] = tmp_vFE;
////                flow_support[vFE[j].first] = neighbor_support_1;
////                flow_support_calculated[vFE[j].first] = true;
////            }
////
////            if(flow_support_calculated[vFE[j].second]) {
////                neighbor_support_2 = flow_support[vFE[j].second];
////            }
////            else{
////                vector<pair<uint32_t, uint32_t> > tmp_vFE;
////                tmp_vFE.reserve(m_);
////                findNeib(m_vv_adj_in[m_vp_edges[vFE[j].second].first], m_vv_adj_in[m_vp_edges[vFE[j].second].second], tmp_vFE);
////                findNeib(m_vv_adj_out[m_vp_edges[vFE[j].second].first], m_vv_adj_out[m_vp_edges[vFE[j].second].second], tmp_vFE);
////                findNeib(m_vv_adj_out[m_vp_edges[vFE[j].second].first], m_vv_adj_in[m_vp_edges[vFE[j].second].second], tmp_vFE);
////                neighbor_support_2 = tmp_vFE.size();
////                eid2_flow_neighbor[vFE[j].second] = tmp_vFE;
////                flow_support[vFE[j].second] = neighbor_support_2;
////                flow_support_calculated[vFE[j].second] = true;
////            }
////            neighbor_edge_support[j] = std::min(neighbor_support_1, neighbor_support_2);
////        }
////        ::uint32_t flow_neighbor_h_index = h_index(neighbor_edge_support);
//
//        if(cycle_neighbor_h_index < m_ui_kc /*|| flow_neighbor_h_index< m_ui_kf*/){
//            continue;
//        }

        /* may be in the D-truss */
        vTargetE.push_back(uiCurEid);
        vbInPool[uiCurEid] = true;

        for (auto atTpE : vCE)
        {
            for (auto atTpEid : {atTpE.first, atTpE.second})
            {
                if (!vbInQ[atTpEid])
                {
                    lsQ.push_back(atTpEid);
                    vbInQ[atTpEid] = true;
                }
            }
        }
        for (auto atTpE : vFE)
        {
            for (auto atTpEid : {atTpE.first, atTpE.second})
            {
                if (!vbInQ[atTpEid])
                {
                    lsQ.push_back(atTpEid);
                    vbInQ[atTpEid] = true;
                }
            }
        }
    }

    vector<bool> pruned_targe_edges(m_,false);
    vector<::uint32_t> pruned_targe_edges_support_vec;

    //pruning based on h-index
//    for(uint32_t iter = 0; iter < iteration_num; iter++){
//        //traverse edges in the vTargetE
//        #pragma omp parallel for num_threads(1)
//        for(uint32_t i = 0; i < vTargetE.size(); i++){
//            ::uint32_t eid = vTargetE[i];
//            vector<::uint32_t> neighbor_edge_support(eid2_cycle_neighbor[eid].size(),0);
//
//            //tranverse the cycle/flow neighbors of eid to get its cycle support vector
//            for(::uint32_t j = 0; j < eid2_cycle_neighbor[eid].size(); j++){
//                pair<::uint32_t,::uint32_t> edges = eid2_cycle_neighbor[eid][j];
//                neighbor_edge_support[j] = std::min(cycle_support[edges.first], cycle_support[edges.second]);
//            }
//            ::uint32_t cycle_neighbor_h_index = h_index(neighbor_edge_support);
//
//            neighbor_edge_support.clear();
//            neighbor_edge_support.resize(eid2_flow_neighbor[eid].size(),0);
//            for(::uint32_t j = 0; j < eid2_flow_neighbor[eid].size(); j++){
//                pair<::uint32_t,::uint32_t> edges = eid2_flow_neighbor[eid][j];
//                neighbor_edge_support[j] = std::min(flow_support[edges.first], flow_support[edges.second]);
//            }
//            ::uint32_t flow_neighbor_h_index = h_index(neighbor_edge_support);
//
//            if(cycle_neighbor_h_index >= m_ui_kc && flow_neighbor_h_index >= m_ui_kf){
//                pruned_targe_edges[eid] = true;
//            }
//        }
//    }
//    for(::uint32_t eid = 0; eid < m_; eid++){
//        if(pruned_targe_edges[eid]){
//            pruned_targe_edges_support_vec.push_back(eid);
//        }
//    }


    //pruning based on h-index in the pruned edges
    /*constrct the adj_lists first*/
    uint32_t n_ = m_v_nodes.size();
    vector<vector<AdjEntry> > adj_in(n_);
    vector<vector<AdjEntry> > adj_out(n_);
    vector<bool> vbTarget(m_, false);
    vector<bool> vbFix(m_, false);
    for (uint32_t eid : vTargetE) {
        const uint32_t v1 = m_vp_edges[eid].first;
        const uint32_t v2 = m_vp_edges[eid].second;
        vbTarget[eid] = true;
        adj_out[v1].push_back({v2, eid});
        adj_in[v2].push_back({v1, eid});
    }
    for (uint32_t eid : vDetE) {
        const uint32_t v1 = m_vp_edges[eid].first;
        const uint32_t v2 = m_vp_edges[eid].second;
        vbFix[eid] = true;
        adj_out[v1].push_back({v2, eid});
        adj_in[v2].push_back({v1, eid});
    }
    for (uint32_t vid = 0; vid < n_; ++vid) {
        std::sort(adj_in[vid].begin(), adj_in[vid].end(),
                  [](const AdjEntry& ae1, const AdjEntry& ae2) {
                      return ae1.pid < ae2.pid;
                  });
        std::sort(adj_out[vid].begin(), adj_out[vid].end(),
                  [](const AdjEntry& ae1, const AdjEntry& ae2) {
                      return ae1.pid < ae2.pid;
                  });
    }
    /*calculate the edge support and cycle/flow neighbors in the inducted graph*/
    vector<::uint32_t> sub_cycle_support(m_, 0);
    vector<::uint32_t> sub_flow_support(m_, 0);
    vector<vector<pair<::uint32_t,::uint32_t>>> eid2_cycle_neighbor(m_), eid2_flow_neighbor(m_);
    pruned_targe_edges_support_vec.clear();
    for(unsigned int eid : vTargetE){
        ::uint32_t x = m_vp_edges[eid].first;
        ::uint32_t y = m_vp_edges[eid].second;
        vCE.clear();
        findNeib(adj_in[x], adj_out[y], vCE);
        vFE.clear();
        findNeib(adj_in[x], adj_in[y], vFE);
        findNeib(adj_out[x], adj_out[y], vFE);
        findNeib(adj_out[x], adj_in[y], vFE);
        sub_cycle_support[eid] = vCE.size();
        sub_flow_support[eid] = vFE.size();
        eid2_cycle_neighbor[eid] = vCE;
        eid2_flow_neighbor[eid] = vFE;
        // if( vCE.size() >= m_ui_kc &&  vFE.size() >= m_ui_kf){
        //     pruned_targe_edges_support_vec.push_back(eid);
        // }
    }
    for(unsigned int eid : vDetE){
        ::uint32_t x = m_vp_edges[eid].first;
        ::uint32_t y = m_vp_edges[eid].second;
        vCE.clear();
        findNeib(adj_in[x], adj_out[y], vCE);
        vFE.clear();
        findNeib(adj_in[x], adj_in[y], vFE);
        findNeib(adj_out[x], adj_out[y], vFE);
        findNeib(adj_out[x], adj_in[y], vFE);
        sub_cycle_support[eid] = vCE.size();
        sub_flow_support[eid] = vFE.size();
        eid2_cycle_neighbor[eid] = vCE;
        eid2_flow_neighbor[eid] = vFE;
        // if( vCE.size() >= m_ui_kc &&  vFE.size() >= m_ui_kf){
        //     pruned_targe_edges_support_vec.push_back(eid);
        // }
    }
    /*calculte the h-index*/
    #pragma omp parallel for num_threads(64)
    for(unsigned int i = 0; i < vTargetE.size(); i++){
       uint32_t eid = vTargetE[i];
       vector<::uint32_t> neighbor_edge_support(eid2_cycle_neighbor[eid].size(),0);

       //tranverse the cycle/flow neighbors of eid to get its cycle support vector
       for(::uint32_t j = 0; j < eid2_cycle_neighbor[eid].size(); j++){
           pair<::uint32_t,::uint32_t> edges = eid2_cycle_neighbor[eid][j];
           neighbor_edge_support[j] = std::min(sub_cycle_support[edges.first], sub_cycle_support[edges.second]);
       }
       ::uint32_t cycle_neighbor_h_index = h_index(neighbor_edge_support);

       neighbor_edge_support.clear();
       neighbor_edge_support.resize(eid2_flow_neighbor[eid].size(),0);
       for(::uint32_t j = 0; j < eid2_flow_neighbor[eid].size(); j++){
           pair<::uint32_t,::uint32_t> edges = eid2_flow_neighbor[eid][j];
           neighbor_edge_support[j] = std::min(sub_flow_support[edges.first], sub_flow_support[edges.second]);
       }

       ::uint32_t flow_neighbor_h_index = h_index(neighbor_edge_support);
       ASSERT_MSG(cycle_neighbor_h_index <= sub_cycle_support[eid] && flow_neighbor_h_index <= sub_flow_support[eid], "h-index error");
    //    if((sub_cycle_support[eid] >= m_ui_kc && sub_flow_support[eid] >= m_ui_kf) /*&& (cycle_neighbor_h_index < m_ui_kc || flow_neighbor_h_index < m_ui_kf)*/){
    //        cout << "eid: " << eid << " cycle_neighbor_h_index: " << cycle_neighbor_h_index << " cycle support " << sub_cycle_support[eid]
    //             << " flow_neighbor_h_index: " << flow_neighbor_h_index << " flow support " << sub_flow_support[eid]  << endl;
    //    }

       if((cycle_neighbor_h_index >= m_ui_kc && flow_neighbor_h_index >= m_ui_kf)){
           pruned_targe_edges[eid] = true;
           //pruned_targe_edges_support_vec.push_back(eid);
       }

   }

   #pragma omp parallel for num_threads(64)
   for(uint32_t i = 0; i < vDetE.size(); i++){
       ::uint32_t eid = vDetE[i];
       vector<::uint32_t> neighbor_edge_support(eid2_cycle_neighbor[eid].size(),0);

       //tranverse the cycle/flow neighbors of eid to get its cycle support vector
       for(::uint32_t j = 0; j < eid2_cycle_neighbor[eid].size(); j++){
           pair<::uint32_t,::uint32_t> edges = eid2_cycle_neighbor[eid][j];
           neighbor_edge_support[j] = std::min(sub_cycle_support[edges.first], sub_cycle_support[edges.second]);
       }
       ::uint32_t cycle_neighbor_h_index = h_index(neighbor_edge_support);

       neighbor_edge_support.clear();
       neighbor_edge_support.resize(eid2_flow_neighbor[eid].size(),0);
       for(::uint32_t j = 0; j < eid2_flow_neighbor[eid].size(); j++){
           pair<::uint32_t,::uint32_t> edges = eid2_flow_neighbor[eid][j];
           neighbor_edge_support[j] = std::min(sub_flow_support[edges.first], sub_flow_support[edges.second]);
       }
       ::uint32_t flow_neighbor_h_index = h_index(neighbor_edge_support);

       if(cycle_neighbor_h_index >= m_ui_kc && flow_neighbor_h_index >= m_ui_kf){
           pruned_targe_edges[eid] = true;
           //pruned_targe_edges_support_vec.push_back(eid);
       }
   }

   for(uint32_t eid = 0; eid < pruned_targe_edges.size(); eid++){
        if(pruned_targe_edges[eid]){
            pruned_targe_edges_support_vec.push_back(eid);
        }
   }

    SearchSpace += vTargetE.size() + vDetE.size();

    vector<uint32_t> vResE;
    vResE.reserve(m_);
    //peeling_with_constrcuted(/*vTargetE*/pruned_targe_edges_support_vec, vDetE, vResE, adj_in, adj_out, vbTarget, vbFix);
    peeling(/*vTargetE*/pruned_targe_edges_support_vec, vDetE, vResE);
    total_pruned_size += pruned_targe_edges_support_vec.size();
    total_unpruned_size+=vTargetE.size();
    total_unpruned_size+=vDetE.size();
    ASSERT_MSG((double)pruned_targe_edges_support_vec.size() / (vTargetE.size()+vDetE.size()) <= 1, "ratio should be less than 1 " << (double)pruned_targe_edges_support_vec.size() / vTargetE.size());
    //cout << "ratio: " << (double)pruned_targe_edges_support_vec.size() / (vTargetE.size()+vDetE.size()) << " vTargetE.size + vDetE.size(): " << vTargetE.size() + vDetE.size()
    //    << " pruned_target_size: " << pruned_targe_edges_support_vec.size() << " vResE.size: " << vResE.size() << endl;
    //cout << "ratio: " << (double)total_pruned_size/total_unpruned_size << " total_pruned_size: " << total_pruned_size << " vResE.size: " << vResE.size() << endl;
    for (auto atE : vResE)
    {
        m_v_EInfo[atE].bInDT = true;
    }


    //pruning based on support

//    uint32_t n_ = m_v_nodes.size();
//    vector<vector<AdjEntry> > adj_in(n_);
//    vector<vector<AdjEntry> > adj_out(n_);
//    vector<bool> vbTarget(m_, false);
//    vector<bool> vbFix(m_, false);
//    for (uint32_t eid : vTargetE) {
//        const uint32_t v1 = m_vp_edges[eid].first;
//        const uint32_t v2 = m_vp_edges[eid].second;
//        vbTarget[eid] = true;
//        adj_out[v1].push_back({v2, eid});
//        adj_in[v2].push_back({v1, eid});
//    }
//    for (uint32_t eid : vDetE) {
//        const uint32_t v1 = m_vp_edges[eid].first;
//        const uint32_t v2 = m_vp_edges[eid].second;
//        vbFix[eid] = true;
//        adj_out[v1].push_back({v2, eid});
//        adj_in[v2].push_back({v1, eid});
//    }
//    for (uint32_t vid = 0; vid < n_; ++vid) {
//        std::sort(adj_in[vid].begin(), adj_in[vid].end(),
//                  [](const AdjEntry& ae1, const AdjEntry& ae2) {
//                      return ae1.pid < ae2.pid;
//                  });
//        std::sort(adj_out[vid].begin(), adj_out[vid].end(),
//                  [](const AdjEntry& ae1, const AdjEntry& ae2) {
//                      return ae1.pid < ae2.pid;
//                  });
//    }
//
//    //vector<::uint32_t> sub_cycle_support(m_, 0);
//    //vector<::uint32_t> sub_flow_support(m_, 0);
//    //vector<vector<pair<::uint32_t,::uint32_t>>> eid2_cycle_neighbor(m_), eid2_flow_neighbor(m_);
//    /*check the support in vTargetE*/
//    pruned_targe_edges_support_vec.clear();
//    for(unsigned int eid : vTargetE){
//        ::uint32_t x = m_vp_edges[eid].first;
//        ::uint32_t y = m_vp_edges[eid].second;
//        vCE.clear();
//        findNeib(adj_in[x], adj_out[y], vCE);
//        vFE.clear();
//        findNeib(adj_in[x], adj_in[y], vFE);
//        findNeib(adj_out[x], adj_out[y], vFE);
//        findNeib(adj_out[x], adj_in[y], vFE);
//        //sub_cycle_support[eid] = vCE.size();
//        //sub_flow_support[eid] = vFE.size();
//        //eid2_cycle_neighbor[eid] = vCE;
//        //eid2_flow_neighbor[eid] = vFE;
//        if( vCE.size() >= m_ui_kc &&  vFE.size() >= m_ui_kf){
//            //pruned_targe_edges[eid] = true;
//            pruned_targe_edges_support_vec.push_back(eid);
//        }
//    }
//    for(unsigned int eid : vDetE){
//        ::uint32_t x = m_vp_edges[eid].first;
//        ::uint32_t y = m_vp_edges[eid].second;
//        vCE.clear();
//        findNeib(adj_in[x], adj_out[y], vCE);
//        vFE.clear();
//        findNeib(adj_in[x], adj_in[y], vFE);
//        findNeib(adj_out[x], adj_out[y], vFE);
//        findNeib(adj_out[x], adj_in[y], vFE);
//        if( vCE.size() >= m_ui_kc &&  vFE.size() >= m_ui_kf){
//            pruned_targe_edges_support_vec.push_back(eid);
//        }
//    }
//
//    //SearchSpace += vTargetE.size();
//    vector<uint32_t> vResE;
//    vResE.reserve(m_);
//    //peeling_with_constrcuted(/*vTargetE*/pruned_targe_edges_support_vec, vDetE, vResE, adj_in, adj_out,vbTarget, vbFix);
//    peeling(pruned_targe_edges_support_vec, vDetE, vResE);
//    total_pruned_size+=pruned_targe_edges_support_vec.size();
//    total_unpruned_size+=vTargetE.size();
//    total_unpruned_size+=vDetE.size();
//    ASSERT_MSG((double)pruned_targe_edges_support_vec.size() / (vTargetE.size()+vDetE.size()) <= 1, "ratio should be less than 1 " << (double)pruned_targe_edges_support_vec.size() / vTargetE.size());
////    cout << "ratio: " << (double)pruned_targe_edges_support_vec.size() / (vTargetE.size()+vDetE.size()) << " vTargetE.size + vDetE.size(): " << vTargetE.size() + vDetE.size()
////         << " pruned_target_size: " << pruned_targe_edges_support_vec.size() << " vResE.size: " << vResE.size() << endl;
//    cout << "ratio: " << (double)total_pruned_size/total_unpruned_size << " total_pruned_size: " << total_pruned_size << " vResE.size: " << vResE.size() << endl;
//    for (auto atE : vResE)
//    {
//        m_v_EInfo[atE].bInDT = true;
//    }

    //printf("Insert single finish D-truss size: %d\n", vResE.size());

    return vResE.size();
}

/**
    add one edge info
**/
uint32_t Update::insertDT(vector<uint32_t> &vInsE)
{
    list<uint32_t> lsQ;
    vector<pair<uint32_t, uint32_t> > vCE;
    vector<pair<uint32_t, uint32_t> > vFE;

    for (uint32_t eid : vInsE)
    {
        uint32_t x = m_vp_edges[eid].first;
        uint32_t y = m_vp_edges[eid].second;
        /* cycle support */
        vCE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_out[y], vCE);
        uint32_t uiCSup = vCE.size();
        //printf("BFS queue get CSup: %d\n", uiCSup);
        if (uiCSup < m_ui_kc)
        {
            /* cannot in D-truss */
            continue;
        }
        /* flow support */
        vFE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_in[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_out[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_in[y], vFE);
        uint32_t uiFSup = vFE.size();
        //printf("BFS queue get FSup: %d\n", uiFSup);
        if (uiFSup < m_ui_kf)
        {
            /* cannot in D-truss */
            continue;
        }
        lsQ.push_back(eid);
    }

    if (lsQ.empty())
    {
        /* nothing to do */
        return 0;
    }

    uint32_t m_ = m_vp_edges.size();
    vector<uint32_t> vTargetE;
    vTargetE.reserve(m_);
    vector<uint32_t> vDetE;
    vDetE.reserve(m_);
    vector<bool> vbInQ(m_, false);
    vector<bool> vbInPool(m_, false);
    for (uint32_t eid : lsQ)
    {
        vbInQ[eid] = true;
    }

    while (!lsQ.empty())
    {
        uint32_t uiCurEid = lsQ.front();
        lsQ.pop_front();
        uint32_t x = m_vp_edges[uiCurEid].first;
        uint32_t y = m_vp_edges[uiCurEid].second;

        //printf("BFS queue get original (%d, %d)\n", m_v_nodes[x], m_v_nodes[y]);
        //printf("BFS queue get (%d, %d) in it: %d\n", x, y, m_v_EInfo[uiCurEid].bInDT);

        if (m_v_EInfo[uiCurEid].bInDT)
        {
            vDetE.push_back(uiCurEid);
            vbInPool[uiCurEid] = true;
            continue;
        }
        /* check support */
        /* cycle support */
        vCE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_out[y], vCE);
        uint32_t uiCSup = vCE.size();
        //printf("BFS queue get CSup: %d\n", uiCSup);
        if (uiCSup < m_ui_kc)
        {
            /* cannot in D-truss */
            continue;
        }

        /* flow support */
        vFE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_in[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_out[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_in[y], vFE);
        uint32_t uiFSup = vFE.size();
        //printf("BFS queue get FSup: %d\n", uiFSup);
        if (uiFSup < m_ui_kf)
        {
            /* cannot in D-truss */
            continue;
        }

        /* may be in the D-truss */
        vTargetE.push_back(uiCurEid);
        vbInPool[uiCurEid] = true;

        for (auto atTpE : vCE)
        {
            for (auto atTpEid : {atTpE.first, atTpE.second})
            {
                if (!vbInQ[atTpEid])
                {
                    lsQ.push_back(atTpEid);
                    vbInQ[atTpEid] = true;
                }
            }
        }
        for (auto atTpE : vFE)
        {
            for (auto atTpEid : {atTpE.first, atTpE.second})
            {
                if (!vbInQ[atTpEid])
                {
                    lsQ.push_back(atTpEid);
                    vbInQ[atTpEid] = true;
                }
            }
        }    
    }

    SearchSpace += vTargetE.size();
    vector<uint32_t> vResE;
    vResE.reserve(m_);
    peeling(vTargetE, vDetE, vResE);
    for (auto atE : vResE)
    {
        m_v_EInfo[atE].bInDT = true;
    }
    //printf("Insert single finish D-truss size: %d\n", vResE.size());

    return vResE.size();
}
/**
    add edges
**/
void Update::add(vector<pair<uint32_t, uint32_t> > & vEdges)
{
    vector<uint32_t> vInsE;
    vInsE.reserve(vEdges.size());
    for (auto atE : vEdges)
    {
        uint32_t eid = 0;
        bool bRes = addEInfo(atE.first, atE.second, &eid);
        if (bRes)
        {
            vInsE.push_back(eid);
        }
    }
    const auto beg = std::chrono::steady_clock::now();
    m_uiDTCnt += insertDT/*insertOptDT*/(vInsE);
    const auto end = std::chrono::steady_clock::now();
    g_tUptIns += end - beg;
    
     m_uiEdgeCnt = m_vp_edges.size();
    /* DEBUG */
    /*int iDTCnt = 0;
    for (uint32_t eid = 0; eid < m_vp_edges.size(); ++eid)
    {
        if (!m_v_EInfo[eid].bUsed)
        {
            continue;
        }
        if (m_v_EInfo[eid].bInDT)
        {
            ++iDTCnt;
        }
    }
    ASSERT(m_uiDTCnt == iDTCnt);*/
    //printf("Insert finish D-truss size: %d\n", iDTCnt);
}
/**
    save community
**/
void Update::save(vector<pair<uint32_t, uint32_t> > & vEdges)
{
    //printf("Total search space: %d\n", SearchSpace);
    //printf("DEBUG save total size: %d\n", m_vp_edges.size());
    for (uint32_t eid = 0; eid < m_vp_edges.size(); ++eid)
    {
        if (!m_v_EInfo[eid].bUsed)
        {
            continue;
        }
        if (m_v_EInfo[eid].bInDT)
        {
            uint32_t x = m_vp_edges[eid].first;
            uint32_t y = m_vp_edges[eid].second;
            vEdges.push_back(pair<uint32_t, uint32_t>(m_v_nodes[x], m_v_nodes[y]));
        }
    }
}
/**
    calculate the CMSin and CMSout of finnal D-truss
**/
pair<double,double> Update::getDtQuality(vector<pair<uint32_t, uint32_t>> & vEdges)
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
    remove signle edge in D-truss
**/
uint32_t Update::rmDT(vector<uint32_t> &vRmE)
{
    list<uint32_t> lsQ;
    for (uint32_t eid : vRmE)
    {
        if (!m_v_EInfo[eid].bInDT)
        {
            /* not in D-truss, remove directly */
            continue;
        }
        lsQ.push_back(eid);
    }

    if (lsQ.empty())
    {
        /* nothing to do */
        return 0;
    }

    uint32_t uiRmCnt = 0;

    vector<pair<uint32_t, uint32_t> > vCE;
    vector<pair<uint32_t, uint32_t> > vFE;
    vector<pair<uint32_t, uint32_t> > vTpCE;
    vector<pair<uint32_t, uint32_t> > vTpFE;
    uint32_t m_ = m_vp_edges.size();
    vector<bool> vbInQ(m_, false);
    vector<bool> vbInPool(m_, false);

    for (uint32_t eid : lsQ)
    {
        vbInQ[eid] = true;
    }

    /* BFS all edges connected and in D-truss */
    while (!lsQ.empty())
    {
        uint32_t uiCurEid = lsQ.front();
        lsQ.pop_front();
        uint32_t x = m_vp_edges[uiCurEid].first;
        uint32_t y = m_vp_edges[uiCurEid].second;

        m_v_EInfo[uiCurEid].bInDT = false;
        ++uiRmCnt;

        //printf("BFS queue get original (%d, %d)\n", m_v_nodes[x], m_v_nodes[y]);
        //printf("BFS queue get (%d, %d) in it: %d\n", x, y, m_v_EInfo[uiCurEid].bInDT);

        /* check support */
        /* cycle support */
        vCE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_out[y], vCE);

        for (auto atTpE : vCE)
        {
            if ((!m_v_EInfo[atTpE.first].bInDT) || (!m_v_EInfo[atTpE.second].bInDT))
            {
                /* not in D-truss */
                continue;
            }
            for (auto atTpEid : {atTpE.first, atTpE.second})
            {
                if (vbInQ[atTpEid])
                {
                    continue;
                }
                uint32_t u = m_vp_edges[atTpEid].first;
                uint32_t v = m_vp_edges[atTpEid].second;
                vTpCE.clear();
                findNeib(m_vv_adj_in[u], m_vv_adj_out[v], vTpCE);
                uint32_t uiCSup = 0;
                for (auto atTpTpE : vTpCE)
                {
                    if ((m_v_EInfo[atTpTpE.first].bInDT) && (m_v_EInfo[atTpTpE.second].bInDT))
                    {
                        ++uiCSup;
                    }
                }
                //printf("BFS queue get CSup: %d\n", uiCSup);
                if (uiCSup < m_ui_kc)
                {
                    /* cannot in D-truss */
                    lsQ.push_back(atTpEid);
                    vbInQ[atTpEid] = true;
                }
            }
        }

        /* flow support */
        vFE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_in[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_out[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_in[y], vFE);

        for (auto atTpE : vFE)
        {
            if ((!m_v_EInfo[atTpE.first].bInDT) || (!m_v_EInfo[atTpE.second].bInDT))
            {
                /* not in D-truss */
                continue;
            }
            for (auto atTpEid : {atTpE.first, atTpE.second})
            {
                if (vbInQ[atTpEid])
                {
                    continue;
                }
                uint32_t u = m_vp_edges[atTpEid].first;
                uint32_t v = m_vp_edges[atTpEid].second;

                vTpFE.clear();
                findNeib(m_vv_adj_in[u], m_vv_adj_in[v], vTpFE);
                findNeib(m_vv_adj_out[u], m_vv_adj_out[v], vTpFE);
                findNeib(m_vv_adj_out[u], m_vv_adj_in[v], vTpFE);
                uint32_t uiFSup = 0;
                for (auto atTpTpE : vTpFE)
                {
                    if ((m_v_EInfo[atTpTpE.first].bInDT) && (m_v_EInfo[atTpTpE.second].bInDT))
                    {
                        ++uiFSup;
                    }
                }
                //printf("BFS queue get FSup: %d\n", uiFSup);
                if (uiFSup < m_ui_kf)
                {
                    /* cannot in D-truss */
                    lsQ.push_back(atTpEid);
                    vbInQ[atTpEid] = true;
                }
            }
        }
    }

    return uiRmCnt;
}

/**
    remove edges
**/
void Update::remove(vector<pair<uint32_t, uint32_t> > & vEdges)
{
    uint32_t eid = 0;

    vector<uint32_t> vRmE;
    vRmE.reserve(vEdges.size());

    for (auto atE : vEdges)
    {
        /* rename */
        uint32_t u = renamePid(atE.first);
        uint32_t v = renamePid(atE.second);
        if (findEid(u, v, &eid))
        {
            /* found */
            --m_v_EInfo[eid].cnt;
            if (1 > m_v_EInfo[eid].cnt)
            {
                /* remove edge in D-truss */
                vRmE.push_back(eid);
            }
        }
    }

    const auto beg = std::chrono::steady_clock::now();
    uint32_t uiRmCnt = rmDT(vRmE);
    const auto end = std::chrono::steady_clock::now();
    g_tUptRm += end - beg;

    m_uiDTCnt -= uiRmCnt;
    for (auto atEid : vRmE)
    {
        uint32_t u =  m_vp_edges[atEid].first;
        uint32_t v =  m_vp_edges[atEid].second;

        /* remove edge information */
        m_v_avaEid.push_back(atEid);
        m_vp_edges[atEid] = pair<uint32_t, uint32_t>(0, 0);
        m_v_EInfo[atEid].bUsed = false;

        /* remove adj information */
        vector<AdjEntry>::iterator it = lower_bound(m_vv_adj_out[u].begin(), m_vv_adj_out[u].end(), v,
                  [](const AdjEntry& e, uint32_t tpPid) {
                    return e.pid < tpPid;
                  });
        if ((m_vv_adj_out[u].end() != it) && (it->pid == v))
        {
            m_vv_adj_out[u].erase(it);
        }
        else
        {
            ASSERT(0);
        }
        it = lower_bound(m_vv_adj_in[v].begin(), m_vv_adj_in[v].end(), u,
                  [](const AdjEntry& e, uint32_t tpPid) {
                    return e.pid < tpPid;
                  });
        if ((m_vv_adj_in[v].end() != it) && (it->pid == u))
        {
            m_vv_adj_in[v].erase(it);
        }
        else
        {
            ASSERT(0);
        }

        /* remove nodes */
        for (uint32_t uiP : {u, v})
        {
            if (1 > m_vv_adj_out[uiP].size() + m_vv_adj_in[uiP].size())
            {
                /* should be removed */
                uint32_t uiInitPid = m_v_nodes[uiP];
                m_v_avaPid.push_back(uiP);
                m_v_nodes[uiP] = 0;
                m_v_bPid[uiInitPid] = false;
                m_v_rePid[uiInitPid] = 0;
            }
        }
    }

    m_uiEdgeCnt = m_vp_edges.size();
    //printf("RM %d edges removed from D-truss\n", uiRmCnt);
}
/**
    get D-truss by peeling
**/
void Update::peeling(vector<uint32_t> &vEdges, vector<uint32_t> &vFixE, vector<uint32_t> &vResE)
{
    const auto beg = std::chrono::steady_clock::now();
    uint32_t n_ = m_v_nodes.size();
    uint32_t m_ = m_vp_edges.size();

    vector<vector<AdjEntry> > adj_in(n_);
    vector<vector<AdjEntry> > adj_out(n_);
    vector<bool> vbTarget(m_, false);
    vector<bool> vbFix(m_, false);

    for (uint32_t eid : vEdges) {
        const uint32_t v1 = m_vp_edges[eid].first;
        const uint32_t v2 = m_vp_edges[eid].second;
        vbTarget[eid] = true;
        adj_out[v1].push_back({v2, eid});
        adj_in[v2].push_back({v1, eid});
    }


    for (uint32_t eid : vFixE) {
        const uint32_t v1 = m_vp_edges[eid].first;
        const uint32_t v2 = m_vp_edges[eid].second;
        vbFix[eid] = true;
        adj_out[v1].push_back({v2, eid});
        adj_in[v2].push_back({v1, eid});
    }

    for (uint32_t vid = 0; vid < n_; ++vid) {
        std::sort(adj_in[vid].begin(), adj_in[vid].end(),
                  [](const AdjEntry& ae1, const AdjEntry& ae2) {
                      return ae1.pid < ae2.pid;
                  });
        std::sort(adj_out[vid].begin(), adj_out[vid].end(),
                  [](const AdjEntry& ae1, const AdjEntry& ae2) {
                      return ae1.pid < ae2.pid;
                  });
    }


    // peeling for D-truss
    // 1. compute the support of each edge by triangle listing
    // 1.1. define a total order over the vertices
    const auto pred = [this, &adj_in, &adj_out](const uint32_t v1, const uint32_t v2) {
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
    vector<uint32_t> v_Sup_c(m_, 0);
    vector<vector<AdjEntry> > A_in(n_);
    vector<vector<AdjEntry> > A_out(n_);


    for (const uint32_t v : verts)
    {
        //printf("PEEL d(%d) = %d\n", v, adj_in[v].size() + adj_out[v].size());
        if (adj_in[v].size() + adj_out[v].size() < 2)
        {
            break;
        }
        for (const auto ae : adj_in[v])
        {
            const uint32_t u = ae.pid;
            const uint32_t e = ae.eid;
            if (pred(u, v)) continue;
            size_t pv = 0, pu = 0;
            while (pv < A_out[v].size() && pu < A_in[u].size())
            {
                if (A_out[v][pv].pid == A_in[u][pu].pid)
                {
                    ++v_Sup_c[e]; ++v_Sup_c[A_out[v][pv].eid]; ++v_Sup_c[A_in[u][pu].eid];
                    ++pv; ++pu;
                }
                else if (pred(A_out[v][pv].pid, A_in[u][pu].pid))
                {
                    ++pv;
                }
                else
                {
                    ++pu;
                }
            }
            A_out[u].push_back({v, e});
        }
        for (const auto ae : adj_out[v])
        {
            const uint32_t u = ae.pid;
            const uint32_t e = ae.eid;
            if (pred(u, v)) continue;
            size_t pv = 0, pu = 0;
            while (pv < A_in[v].size() && pu < A_out[u].size())
            {
                /*printf("PEEL counting start %d, %d, %d, %d\n",
                       v, u, A_in[v][pv].pid, A_out[u][pu].pid);*/
                if (A_in[v][pv].pid == A_out[u][pu].pid)
                {
                    /*printf("PEEL counting get triangle (%d, %d, %d, %d) original (%d, %d, %d)\n",
                           v, u, A_in[v][pv].pid, A_out[u][pu].pid,
                            m_v_nodes[v], m_v_nodes[u], m_v_nodes[A_in[v][pv].pid]);*/

                    ++v_Sup_c[e]; ++v_Sup_c[A_in[v][pv].eid]; ++v_Sup_c[A_out[u][pu].eid];
                    ++pv; ++pu;
                }
                else if (pred(A_in[v][pv].pid, A_out[u][pu].pid))
                {
                    ++pv;
                }
                else
                {
                    ++pu;
                }
            }
            A_in[u].push_back({v, e});

            /*printf("PEEL counting get (%d, %d) original (%d, %d), sup: %d\n",
                v, u,
                m_v_nodes[v], m_v_nodes[u],
                v_Sup_c[e]);*/
        }
    }
    decltype(A_in)().swap(A_in);
    decltype(A_out)().swap(A_out);
    decltype(verts)().swap(verts);

  // 1.4. list flow triangles
  vector<uint32_t> v_Sup_f(m_, 0);

  for (uint32_t eid = 0; eid < m_; ++eid) {
    if ((!vbTarget[eid]) && (!vbFix[eid]))
    {
        continue;
    }
    const uint32_t v1 = m_vp_edges[eid].first;
    const uint32_t v2 = m_vp_edges[eid].second;
    vector<pair<uint32_t, uint32_t> > tris;
    findNeib(adj_out[v1], adj_in[v2], tris);
    for (const auto tri : tris)
    {
        const uint32_t e1 = tri.first;
        const uint32_t e2 = tri.second;
        ++v_Sup_f[eid];++v_Sup_f[e1];++v_Sup_f[e2];
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
    uint32_t uiENum = 0;
    vector<uint32_t> vWait(m_);
    vector<uint32_t> vWaitRmCache(m_);
    // init cache
    vWaitRmCache.clear();
    for (uint32_t eid : vEdges)
    {
        ++uiENum;

        const uint32_t v1 = m_vp_edges[eid].first;
        const uint32_t v2 = m_vp_edges[eid].second;
        /*printf("PEEL init get (%d, %d) original (%d, %d) reject: (%d, %d), sup: %d %d\n",
            v1, v2,
            m_v_nodes[v1], m_v_nodes[v2],
            m_v_rePid[m_v_nodes[v1]], m_v_rePid[m_v_nodes[v2]],
            v_Sup_c[eid], v_Sup_f[eid]);*/
        if (v_Sup_c[eid] < m_ui_kc)
        {
            vWaitRmFlag[eid] = true;
            vWaitRmCache.push_back(eid);
        }
        else if (v_Sup_f[eid] < m_ui_kf)
        {
            vWaitRmFlag[eid] = true;
            vWaitRmCache.push_back(eid);
        }
    }
    while (uiRmCnt < uiENum)
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
            vector<pair<uint32_t, uint32_t>> tris;
            const uint32_t v1 = m_vp_edges[eid].first;
            const uint32_t v2 = m_vp_edges[eid].second;

            /*printf("PEEL current (%d, %d) sup: %d %d\n",
                m_v_nodes[v1], m_v_nodes[v2],
                v_Sup_c[eid], v_Sup_f[eid]);*/

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
                if ((!vWaitRmFlag[e1]) && (vbTarget[e1]))
                {
                    --v_Sup_f[e1];
                    if (v_Sup_f[e1] < m_ui_kf)
                    {
                        vWaitRmFlag[e1] = true;
                        vWaitRmCache.push_back(e1);
                    }
                }
                if ((!vWaitRmFlag[e2]) && (vbTarget[e2]))
                {
                    --v_Sup_f[e2];
                    if (v_Sup_f[e2] < m_ui_kf)
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
                if ((!vWaitRmFlag[e1]) && (vbTarget[e1]))
                {
                    --v_Sup_c[e1];
                    /*printf("PEEL decrease (%d, %d) sup: %d due to (%d, %d)\n",
                        m_v_nodes[m_vp_edges[e1].first], m_v_nodes[m_vp_edges[e1].second],
                        v_Sup_c[e1],
                        m_v_nodes[v1], m_v_nodes[v2]);*/
                    if (v_Sup_c[e1] < m_ui_kc)
                    {
                        vWaitRmFlag[e1] = true;
                        vWaitRmCache.push_back(e1);
                    }
                }
                if ((!vWaitRmFlag[e2]) && (vbTarget[e2]))
                {
                    --v_Sup_c[e2];
                    /*printf("PEEL decrease (%d, %d) sup: %d due to (%d, %d)\n",
                        m_v_nodes[m_vp_edges[e2].first], m_v_nodes[m_vp_edges[e2].second],
                        v_Sup_c[e2],
                        m_v_nodes[v1], m_v_nodes[v2]);*/
                    if (v_Sup_c[e2] < m_ui_kc)
                    {
                        vWaitRmFlag[e2] = true;
                        vWaitRmCache.push_back(e2);
                    }
                }
            }
        }
    }

  // keep result
  for (uint32_t eid : vEdges)
  {
      if (!removed[eid])
      {
          /*const uint32_t v1 = m_vp_edges[eid].first;
          const uint32_t v2 = m_vp_edges[eid].second;
          printf("PEEL in it %d (%d, %d)\n", eid,
                   m_v_nodes[v1], m_v_nodes[v2]);*/

          vResE.push_back(eid);
      }
  }

  const auto end = std::chrono::steady_clock::now();
  const auto dif = end - beg;
    peeling_time += end - beg;
  /*printf("Peeling costs \x1b[1;31m%f\x1b[0m ms.\n",
         std::chrono::duration<double, std::milli>(dif).count());*/
    //printf("DEBUG index done last: (%d, %d)\n", m_vp_edges[m_ - 1].first, m_vp_edges[m_ - 1].second);
}
/**
    get D-truss by peeling with constructed adjacency list
**/
void Update::peeling_with_constrcuted(vector<uint32_t> &vEdges, vector<uint32_t> &vFixE, vector<uint32_t> &vResE,
                                      vector<vector<Update::AdjEntry>> &adj_in,
                                      vector<vector<Update::AdjEntry>> &adj_out,
                                      vector<bool> &vbTarget,
                                      vector<bool> &vbFix) {
    const auto beg = std::chrono::steady_clock::now();
    uint32_t n_ = m_v_nodes.size();
    uint32_t m_ = m_vp_edges.size();

    // peeling for D-truss
    // 1. compute the support of each edge by triangle listing
    // 1.1. define a total order over the vertices
    const auto pred = [this, &adj_in, &adj_out](const uint32_t v1, const uint32_t v2) {
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
    vector<uint32_t> v_Sup_c(m_, 0);
    vector<vector<AdjEntry> > A_in(n_);
    vector<vector<AdjEntry> > A_out(n_);


    for (const uint32_t v : verts)
    {
        //printf("PEEL d(%d) = %d\n", v, adj_in[v].size() + adj_out[v].size());
        if (adj_in[v].size() + adj_out[v].size() < 2)
        {
            break;
        }
        for (const auto ae : adj_in[v])
        {
            const uint32_t u = ae.pid;
            const uint32_t e = ae.eid;
            if (pred(u, v)) continue;
            size_t pv = 0, pu = 0;
            while (pv < A_out[v].size() && pu < A_in[u].size())
            {
                if (A_out[v][pv].pid == A_in[u][pu].pid)
                {
                    ++v_Sup_c[e]; ++v_Sup_c[A_out[v][pv].eid]; ++v_Sup_c[A_in[u][pu].eid];
                    ++pv; ++pu;
                }
                else if (pred(A_out[v][pv].pid, A_in[u][pu].pid))
                {
                    ++pv;
                }
                else
                {
                    ++pu;
                }
            }
            A_out[u].push_back({v, e});
        }
        for (const auto ae : adj_out[v])
        {
            const uint32_t u = ae.pid;
            const uint32_t e = ae.eid;
            if (pred(u, v)) continue;
            size_t pv = 0, pu = 0;
            while (pv < A_in[v].size() && pu < A_out[u].size())
            {
                /*printf("PEEL counting start %d, %d, %d, %d\n",
                       v, u, A_in[v][pv].pid, A_out[u][pu].pid);*/
                if (A_in[v][pv].pid == A_out[u][pu].pid)
                {
                    /*printf("PEEL counting get triangle (%d, %d, %d, %d) original (%d, %d, %d)\n",
                           v, u, A_in[v][pv].pid, A_out[u][pu].pid,
                            m_v_nodes[v], m_v_nodes[u], m_v_nodes[A_in[v][pv].pid]);*/

                    ++v_Sup_c[e]; ++v_Sup_c[A_in[v][pv].eid]; ++v_Sup_c[A_out[u][pu].eid];
                    ++pv; ++pu;
                }
                else if (pred(A_in[v][pv].pid, A_out[u][pu].pid))
                {
                    ++pv;
                }
                else
                {
                    ++pu;
                }
            }
            A_in[u].push_back({v, e});

            /*printf("PEEL counting get (%d, %d) original (%d, %d), sup: %d\n",
                v, u,
                m_v_nodes[v], m_v_nodes[u],
                v_Sup_c[e]);*/
        }
    }
    decltype(A_in)().swap(A_in);
    decltype(A_out)().swap(A_out);
    decltype(verts)().swap(verts);

    // 1.4. list flow triangles
    vector<uint32_t> v_Sup_f(m_, 0);

    for (uint32_t eid = 0; eid < m_; ++eid) {
        if ((!vbTarget[eid]) && (!vbFix[eid]))
        {
            continue;
        }
        const uint32_t v1 = m_vp_edges[eid].first;
        const uint32_t v2 = m_vp_edges[eid].second;
        vector<pair<uint32_t, uint32_t> > tris;
        findNeib(adj_out[v1], adj_in[v2], tris);
        for (const auto tri : tris)
        {
            const uint32_t e1 = tri.first;
            const uint32_t e2 = tri.second;
            ++v_Sup_f[eid];++v_Sup_f[e1];++v_Sup_f[e2];
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
    uint32_t uiENum = 0;
    vector<uint32_t> vWait(m_);
    vector<uint32_t> vWaitRmCache(m_);
    // init cache
    vWaitRmCache.clear();
    for (uint32_t eid : vEdges)
    {
        ++uiENum;

        const uint32_t v1 = m_vp_edges[eid].first;
        const uint32_t v2 = m_vp_edges[eid].second;
        /*printf("PEEL init get (%d, %d) original (%d, %d) reject: (%d, %d), sup: %d %d\n",
            v1, v2,
            m_v_nodes[v1], m_v_nodes[v2],
            m_v_rePid[m_v_nodes[v1]], m_v_rePid[m_v_nodes[v2]],
            v_Sup_c[eid], v_Sup_f[eid]);*/
        if (v_Sup_c[eid] < m_ui_kc)
        {
            vWaitRmFlag[eid] = true;
            vWaitRmCache.push_back(eid);
        }
        else if (v_Sup_f[eid] < m_ui_kf)
        {
            vWaitRmFlag[eid] = true;
            vWaitRmCache.push_back(eid);
        }
    }
    while (uiRmCnt < uiENum)
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
            vector<pair<uint32_t, uint32_t>> tris;
            const uint32_t v1 = m_vp_edges[eid].first;
            const uint32_t v2 = m_vp_edges[eid].second;

            /*printf("PEEL current (%d, %d) sup: %d %d\n",
                m_v_nodes[v1], m_v_nodes[v2],
                v_Sup_c[eid], v_Sup_f[eid]);*/

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
                if ((!vWaitRmFlag[e1]) && (vbTarget[e1]))
                {
                    --v_Sup_f[e1];
                    if (v_Sup_f[e1] < m_ui_kf)
                    {
                        vWaitRmFlag[e1] = true;
                        vWaitRmCache.push_back(e1);
                    }
                }
                if ((!vWaitRmFlag[e2]) && (vbTarget[e2]))
                {
                    --v_Sup_f[e2];
                    if (v_Sup_f[e2] < m_ui_kf)
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
                if ((!vWaitRmFlag[e1]) && (vbTarget[e1]))
                {
                    --v_Sup_c[e1];
                    /*printf("PEEL decrease (%d, %d) sup: %d due to (%d, %d)\n",
                        m_v_nodes[m_vp_edges[e1].first], m_v_nodes[m_vp_edges[e1].second],
                        v_Sup_c[e1],
                        m_v_nodes[v1], m_v_nodes[v2]);*/
                    if (v_Sup_c[e1] < m_ui_kc)
                    {
                        vWaitRmFlag[e1] = true;
                        vWaitRmCache.push_back(e1);
                    }
                }
                if ((!vWaitRmFlag[e2]) && (vbTarget[e2]))
                {
                    --v_Sup_c[e2];
                    /*printf("PEEL decrease (%d, %d) sup: %d due to (%d, %d)\n",
                        m_v_nodes[m_vp_edges[e2].first], m_v_nodes[m_vp_edges[e2].second],
                        v_Sup_c[e2],
                        m_v_nodes[v1], m_v_nodes[v2]);*/
                    if (v_Sup_c[e2] < m_ui_kc)
                    {
                        vWaitRmFlag[e2] = true;
                        vWaitRmCache.push_back(e2);
                    }
                }
            }
        }
    }

    // keep result
    for (uint32_t eid : vEdges)
    {
        if (!removed[eid])
        {
            /*const uint32_t v1 = m_vp_edges[eid].first;
            const uint32_t v2 = m_vp_edges[eid].second;
            printf("PEEL in it %d (%d, %d)\n", eid,
                     m_v_nodes[v1], m_v_nodes[v2]);*/

            vResE.push_back(eid);
        }
    }

    const auto end = std::chrono::steady_clock::now();
    const auto dif = end - beg;
    peeling_time += end - beg;
    /*printf("Peeling costs \x1b[1;31m%f\x1b[0m ms.\n",
           std::chrono::duration<double, std::milli>(dif).count());*/
    //printf("DEBUG index done last: (%d, %d)\n", m_vp_edges[m_ - 1].first, m_vp_edges[m_ - 1].second);
}