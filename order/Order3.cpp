//
// Created by xkliao on 2023/6/9.
//
//
// Created by xkliao on 2023/6/7.
//
#include "Order3.h"

#include <cstdlib>
#include <fstream>
#include <numeric>
#include <utility>


//::uint32_t total_pruned_size = 0, total_unpruned_size = 0;
//chrono::duration<double> peeling_time;
chrono::duration<double> g_tUptIns_order_3;
chrono::duration<double> g_tUptRm_order_3;
//int SearchSpace = 0;
//long long int MAX_LABEL_2 = 0xfffffffffffffff;  //the maximum tag of each edge
//long int INIT_TAG_GAP_2 = 10; //32 bit integer
/**
    find neighbor
**/
void Order3::findNeib(vector<AdjEntry> &vAdj1, vector<AdjEntry> &vAdj2, vector<pair<uint32_t, uint32_t> > & vTrdE)
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
    remove a given item from a vector
 **/
void Order3::removeItem(vector<uint32_t> &v, uint32_t item) {
    for (auto it = v.begin(); it != v.end(); ++it) {
        if (*it == item) {
            v.erase(it);
            break;
        }
    }
    ASSERT_MSG(true, "removeItem do not found");
}
/**
 compress the layers when finished edge insertion/deletion, including removing empty layers and compress layers to the L_{mu}
**/
void Order3::compressUpdatedLayers() {
    vector<pair<uint32_t, uint32_t> > vCE;
    vector<pair<uint32_t, uint32_t> > vFE;

    //remove empty layers
    for(uint32_t layer_id = 0; layer_id <= m_mu; ++layer_id){
        if(m_vv_ordered_layers[layer_id].edges.empty()){
            m_vv_ordered_layers.erase(m_vv_ordered_layers.begin() + layer_id);
            --m_mu;
            for(uint32_t sub_layer_id = layer_id; sub_layer_id <= m_mu; ++sub_layer_id){
                m_vv_ordered_layers[sub_layer_id].layer_number--;
                for(uint32_t sub_eid : m_vv_ordered_layers[sub_layer_id].edges){
                    if(m_vv_ordered_list[sub_eid].layer_number > 0)m_vv_ordered_list[sub_eid].layer_number--;
                    ASSERT_MSG(m_vv_ordered_list[sub_eid].layer_number < 50000, " sub_eid " << sub_eid << " " << m_vv_ordered_list[sub_eid].layer_number);
                }
                ASSERT(m_vv_ordered_layers[sub_layer_id].layer_number >= 0);
            }
        }
    }
    /*compress qualified layers to L_{mu}*/
    if(m_mu > 0){
        for(uint32_t layer_id = m_mu - 1; layer_id >= 0; --layer_id){
            bool be_compressed = true;
            if(!(m_vv_ordered_layers[layer_id].edges.empty())){
                for(uint32_t cur_eid : m_vv_ordered_layers[layer_id].edges){
                    uint32_t u = m_vp_edges[cur_eid].first;
                    uint32_t v = m_vp_edges[cur_eid].second;
                    vector<uint32_t> sub_support_cnt_c(m_mu + 1 , 0);
                    vector<uint32_t> sub_support_cnt_f(m_mu + 1, 0);
                    vCE.clear();
                    findNeibAndLayer(m_vv_adj_in[u], m_vv_adj_out[v], vCE, sub_support_cnt_c);
                    vFE.clear();
                    findNeibAndLayer(m_vv_adj_in[u], m_vv_adj_in[v], vFE, sub_support_cnt_f);
                    findNeibAndLayer(m_vv_adj_out[u], m_vv_adj_out[v], vFE, sub_support_cnt_f);
                    findNeibAndLayer(m_vv_adj_out[u], m_vv_adj_in[v], vFE, sub_support_cnt_f);
                    if(!((sub_support_cnt_c[m_mu -1] >= m_ui_kc && sub_support_cnt_f[m_mu -1] >= m_ui_kf) &&
                         (sub_support_cnt_c[m_mu ] < m_ui_kc || sub_support_cnt_f[m_mu ] < m_ui_kf))){
                        be_compressed = false;
                        break;
                    }
                }
            }

            if(be_compressed){
                for(uint32_t sub_eid : m_vv_ordered_layers[m_mu].edges){
                    if(m_vv_ordered_list[sub_eid].layer_number > 0)m_vv_ordered_list[sub_eid].layer_number--;
                    ASSERT_MSG(m_vv_ordered_list[sub_eid].layer_number < 50000, " sub_eid " << sub_eid << " " << m_vv_ordered_list[sub_eid].layer_number);
                }
                m_vv_ordered_layers[layer_id].edges.insert(m_vv_ordered_layers[layer_id].edges.end(), m_vv_ordered_layers[m_mu].edges.begin(), m_vv_ordered_layers[m_mu].edges.end());
                m_vv_ordered_layers.pop_back();
                --m_mu;
                m_vv_ordered_layers[m_mu].layer_number = m_mu;
                if(layer_id == 0) break;
            }
            else{
                break;
            }
        }
    }

}
/**
  when new a layer is formed right next to L_{mu}, update m_mu, update ordered list and layers
**/
void Order3::insertNewLayer(uint32_t eid) {
    //orderInsert(m_vv_ordered_layers[m_mu - 1].edges.back(), eid);
    m_vv_ordered_list[eid].layer_number = m_mu;
    m_vv_ordered_layers.back().layer_number++;

    vector<uint32_t> tmp;
    tmp.push_back(eid);
    m_vv_ordered_layers.insert(m_vv_ordered_layers.end() -1, LayerEntry{tmp,m_mu/*, (int)eid*/});
    ++m_mu;

    for(uint32_t sub_eid : m_vv_ordered_layers[m_mu].edges){
        m_vv_ordered_list[sub_eid].layer_number++;
    }
}
/**
 find the layer where edges should be inserted, do the update accordingly
**/
void Order3::findLayerToInsert(vector<uint32_t> &support_cnt_c, vector<uint32_t> &support_cnt_f, uint32_t eid) {
    uint32_t layer_id = 0;
    if(support_cnt_c[0] < m_ui_kc || support_cnt_f[0] < m_ui_kf){
        layer_id = -1;
    }
    else if(support_cnt_c[m_mu] >= m_ui_kc && support_cnt_f[m_mu] >= m_ui_kf){
        layer_id = m_mu - 1;
    }
    else{
        for(int i = 0; i < m_mu ; ++i){
            if((support_cnt_c[i] >= m_ui_kc && support_cnt_f[i] >= m_ui_kf) &&
               (support_cnt_c[i + 1] < m_ui_kc || support_cnt_f[i + 1] < m_ui_kf)){
                layer_id = i;
                break;
            }
        }
    }

    /* insert the edge into the new layer*/
    //orderInsert(m_vv_ordered_layers[layer_id + 1].edges[0], eid);
    m_vv_ordered_layers[layer_id + 1].edges.push_back(eid);
    m_vv_ordered_list[eid].layer_number = layer_id + 1;
}
/**
insert the affected neighboring edges into the priority queue for further checking
**/
void Order3::insertNeighborsToQueue(vector<bool> &be_in_pq, PQ &m_PQ, vector<vector<pair<uint32_t, uint32_t>>> &neighbors_of_different_layer,
                                    uint32_t lower_bound, uint32_t upper_bound) {
//    if(lower_bound >= upper_bound){
//        return;
//    }
    ASSERT_MSG(lower_bound <= upper_bound, "lower_bound " << lower_bound << " upper_bound " << upper_bound << " " << m_mu << " " << neighbors_of_different_layer.size());
    // lower_bound = 0;
    // upper_bound = neighbors_of_different_layer.size() - 1;
    for(uint32_t i = lower_bound; i <= upper_bound; ++i){
        for(auto neighbor : neighbors_of_different_layer[i]){
            //ASSERT_MSG(neighbor.first < be_in_pq.size() && neighbor.second < be_in_pq.size(), "neighbor.first " << neighbor.first << " neighbor.second " << neighbor.second << " be_in_pq.size() " << be_in_pq.size());
            ASSERT_MSG(neighbor.first < m_vv_ordered_list.size() && neighbor.second < m_vv_ordered_list.size(), "neighbor.first " << neighbor.first << " neighbor.second " << neighbor.second << " m_vv_ordered_list.size() " << m_vv_ordered_list.size());
            if(!be_in_pq[neighbor.first]){
                m_PQ.push(m_vv_ordered_list[neighbor.first]);
                be_in_pq[neighbor.first] = true;
            }
            if(!be_in_pq[neighbor.second]){
                m_PQ.push(m_vv_ordered_list[neighbor.second]);
                be_in_pq[neighbor.second] = true;
            }
        }
    }
}

/**
 calculate the h-index a vector
**/
uint32_t Order3::h_index(vector<::uint32_t> &neigh_cycle_support) {
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
Order3::Order3(vector<pair<uint32_t, uint32_t> > & vEdges, uint32_t uiMaxN, uint32_t uiMaxM, uint32_t uiKc, uint32_t uiKf)
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
        if (bRes)
        {
            vTargetE.push_back(eid);
        }
    }

    m_vv_ordered_list.resize(m_vp_edges.size());
    m_vv_ordered_layers.reserve(m_vp_edges.size());

    vector<uint32_t> vFixE;
    vector<uint32_t> vResE;

    // m_vv_ordered_list and m_vv_ordered_layers are initialized
    peeling(vTargetE, vFixE, vResE);

    printf("m_vv_ordered_layers.size: %d, m_vv_ordered_list.size: %d \n", m_vv_ordered_layers.size(), m_vv_ordered_list.size());
    printf("L[mu] size: %d\n", m_vv_ordered_layers[m_vv_ordered_layers.size() - 1].edges.size());

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
bool Order3::findEid(uint32_t x, uint32_t y, uint32_t *piEid)
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
uint32_t Order3::renamePid(uint32_t x)
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
bool Order3::addEInfo(uint32_t x, uint32_t y, uint32_t *puiEid)
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
        m_vv_ordered_list[*puiEid].eid = *puiEid;
        m_vv_ordered_list[*puiEid].layer_number = 0;
    }
    else
    {
        *puiEid = m_vp_edges.size();
        m_vp_edges.push_back(pair<uint32_t, uint32_t>(u, v));
        OrderEdgeEntry stOrderE = {*puiEid,0};
        m_vv_ordered_list.push_back(stOrderE);
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
    add one edge info
**/
uint32_t Order3::insertDT(vector<uint32_t> &vInsE)
{
    uint32_t m_ = m_vp_edges.size();

    //for process edge one by one
    vector<bool> be_in_graph(m_, true);
    for (uint32_t eid : vInsE)
    {
        be_in_graph[eid] = false;
    }
    //process edges one by one
    for (uint32_t eid : vInsE)
    {
        be_in_graph[eid] = true;
        /*to check which layer should we put the new edge in*/
        vector<pair<uint32_t, uint32_t> > vCE;
        vector<pair<uint32_t, uint32_t> > vFE;
        vector<pair<uint32_t,uint32_t>> v_sub_CE;
        vector<pair<uint32_t,uint32_t>> v_sub_FE;
        PQ m_PQ;
        vector<bool> be_in_pq(m_, false);
        be_in_pq[eid] = true;
        uint32_t x = m_vp_edges[eid].first;
        uint32_t y = m_vp_edges[eid].second;
        v_sub_CE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_out[y], v_sub_CE);
        v_sub_FE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_in[y], v_sub_FE);
        findNeib(m_vv_adj_out[x], m_vv_adj_out[y], v_sub_FE);
        findNeib(m_vv_adj_out[x], m_vv_adj_in[y], v_sub_FE);

        /*get its' support cnt in with different layers*/
        vector<uint32_t> support_cnt_c(m_mu + 1,0); // for example, support_cnt_c[1] is the cycle_support in G \ L_{1}
        vector<uint32_t> support_cnt_f(m_mu + 1,0);
        vector<vector<pair< uint32_t, uint32_t>>> neighbors_of_different_layer(m_mu + 1);
        /*cycle triangle*/
        for (auto atTpE : v_sub_CE)
        {
            if (!(be_in_graph[atTpE.first]) || !(be_in_graph[atTpE.second]))
            {
                /* cannot form triangles */
                continue;
            }
            for(uint32_t layer_id = 0; layer_id <= COMMON_MIN(m_vv_ordered_list[atTpE.first].layer_number, m_vv_ordered_list[atTpE.second].layer_number) ; ++layer_id) {
                support_cnt_c[layer_id] += 1;
                neighbors_of_different_layer[layer_id].push_back(atTpE);
            }
            if(m_vv_ordered_list[atTpE.first].layer_number < m_mu){
                //if(!be_in_pq[atTpE.first]) {
                    m_PQ.push(m_vv_ordered_list[atTpE.first]);
                    be_in_pq[atTpE.first] = true;
                //}
            }
            if(m_vv_ordered_list[atTpE.second].layer_number < m_mu){
                //if(!be_in_pq[atTpE.second]) {
                    m_PQ.push(m_vv_ordered_list[atTpE.second]);
                    be_in_pq[atTpE.second] = true;
                //}
            }
        }
        /*flow triangle*/
        for (auto atTpE : v_sub_FE)
        {
            if (!(be_in_graph[atTpE.first]) || !(be_in_graph[atTpE.second]))
            {
                /* cannot form triangles */
                continue;
            }
            for(uint32_t layer_id = 0; layer_id <= COMMON_MIN(m_vv_ordered_list[atTpE.first].layer_number, m_vv_ordered_list[atTpE.second].layer_number) ; ++layer_id) {
                support_cnt_f[layer_id] += 1;
                neighbors_of_different_layer[layer_id].push_back(atTpE);
            }
            if(m_vv_ordered_list[atTpE.first].layer_number < m_mu){
                //if(!be_in_pq[atTpE.first]) {
                    m_PQ.push(m_vv_ordered_list[atTpE.first]);
                    be_in_pq[atTpE.first] = true;
                //}
            }
            if(m_vv_ordered_list[atTpE.second].layer_number < m_mu){
                //if(!be_in_pq[atTpE.second]) {
                    m_PQ.push(m_vv_ordered_list[atTpE.second]);
                    be_in_pq[atTpE.second] = true;
               //}
            }
        }
        /*new layer next to L_{mu} is formed*/
        m_vv_ordered_list[eid].eid = eid;
        // if((support_cnt_c[m_mu -1] >= m_ui_kc && support_cnt_f[m_mu -1] >= m_ui_kf) &&
        //    (support_cnt_c[m_mu] < m_ui_kc || support_cnt_f[m_mu] < m_ui_kf) ){
        //     insertNewLayer(eid);
        // }
        // else{
        //     findLayerToInsert(support_cnt_c,support_cnt_f, eid);
        // }
        findLayerToInsert(support_cnt_c,support_cnt_f, eid);
        while (!m_PQ.empty()){
            uint32_t cur_eid =  m_PQ.top().eid;
            m_PQ.pop();
            be_in_pq[cur_eid] = true;
            uint32_t original_layer = m_vv_ordered_list[cur_eid].layer_number;
            /*get the neighbor of cur_eid*/
            x = m_vp_edges[cur_eid].first;
            y = m_vp_edges[cur_eid].second;
            vCE.clear();
            findNeib(m_vv_adj_in[x], m_vv_adj_out[y], vCE);
            vFE.clear();
            findNeib(m_vv_adj_in[x], m_vv_adj_in[y], vFE);
            findNeib(m_vv_adj_out[x], m_vv_adj_out[y], vFE);
            findNeib(m_vv_adj_out[x], m_vv_adj_in[y], vFE);
            /*get its' support cnt in with different layers*/
            vector<uint32_t> sub_support_cnt_c(m_mu + 1,0); // for example, support_cnt_c[1] is the cycle_support in G \ L_{1}
            vector<uint32_t> sub_support_cnt_f(m_mu + 1,0);
            vector<vector<pair<uint32_t, uint32_t>>> sub_neighbors_of_different_layer(m_mu + 1);
            /*cycle triangle*/
            for (auto atTpE : vCE)
            {
                if (!(be_in_graph[atTpE.first]) || !(be_in_graph[atTpE.second]))
                {
                    /* cannot form triangles */
                    continue;
                }

                for(uint32_t layer_id = 0; layer_id <= COMMON_MIN(m_vv_ordered_list[atTpE.first].layer_number, m_vv_ordered_list[atTpE.second].layer_number) ; ++layer_id) {
                    sub_support_cnt_c[layer_id] += 1;
                    sub_neighbors_of_different_layer[layer_id].push_back(atTpE);
                }
            }
            /*flow triangle*/
            for (auto atTpE : vFE)
            {
                if (!(be_in_graph[atTpE.first]) || !(be_in_graph[atTpE.second]))
                {
                    /* cannot form triangles */
                    continue;
                }
                for(uint32_t layer_id = 0; layer_id <= COMMON_MIN(m_vv_ordered_list[atTpE.first].layer_number, m_vv_ordered_list[atTpE.second].layer_number) ; ++layer_id) {
                    sub_support_cnt_f[layer_id] += 1;
                    sub_neighbors_of_different_layer[layer_id].push_back(atTpE);
                }
            }
            /*new layer next to L_{mu} is formed*/
//            if((sub_support_cnt_c[m_mu - 1] >= m_ui_kc && sub_support_cnt_f[m_mu - 1] >= m_ui_kf) &&
//               (sub_support_cnt_c[m_mu] < m_ui_kc || sub_support_cnt_f[m_mu] < m_ui_kf) ){
//                cout << "insert test 7" << endl;
//                //orderDelete(cur_eid);
//                removeItem(m_vv_ordered_layers[original_layer].edges,cur_eid);
//                insertNewLayer(cur_eid);
//                cout << "insert test 7.1" << endl;
//                /*insert neighbors of cur_eid into the Queue*/
//                uint32_t layer_lower_bound = COMMON_MIN(original_layer + 1, m_mu - 1);
//                uint32_t layer_upper_bound = m_mu - 1;
//                cout << "insert test 7.2" << endl;
//                if(original_layer == m_vv_ordered_list[cur_eid].layer_number) continue;
//                insertNeighborsToQueue(be_in_pq, m_PQ, sub_neighbors_of_different_layer,
//                                       layer_lower_bound, layer_upper_bound);
//                cout << "insert test 7.3" << endl;
//            }
//            else{

            /*find the layer that the edge should be inserted*/
            //orderDelete(cur_eid);
            removeItem(m_vv_ordered_layers[original_layer].edges,cur_eid);
            findLayerToInsert(sub_support_cnt_c,sub_support_cnt_f, cur_eid);
            if(original_layer == m_vv_ordered_list[cur_eid].layer_number) continue;
            /*insert neighbors of cur_eid into the Queue*/
            uint32_t layer_lower_bound = COMMON_MIN(original_layer + 1, m_mu);
            uint32_t layer_upper_bound = m_mu;
            insertNeighborsToQueue(be_in_pq, m_PQ, sub_neighbors_of_different_layer,
                                   layer_lower_bound, layer_upper_bound);

            //}
        }

        //compress the updated layers
        /*delete empty layers*/
        //compressUpdatedLayers();

    }
    //return the DT size
    vector<uint32_t> vResE;
    vResE.reserve(m_);
    for(auto eid : m_vv_ordered_layers.back().edges){
        vResE.push_back(eid);
        m_v_EInfo[eid].bInDT = true;
    }
    return vResE.size();
}
/**
    add edges
**/
void Order3::add(vector<pair<uint32_t, uint32_t> > & vEdges)
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
    m_uiDTCnt += insertDT/*insert_hindex_pruned_DT*/(vInsE);
    const auto end = std::chrono::steady_clock::now();
    g_tUptIns_order_3 += end - beg;

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
void Order3::save(vector<pair<uint32_t, uint32_t> > & vEdges)
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
pair<double,double> Order3::getDtQuality(vector<pair<uint32_t, uint32_t>> & vEdges)
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
uint32_t Order3::rm_bfs_DT(vector<uint32_t> &vRmE)
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
    remove single edge in D-truss
**/
uint32_t Order3::rmDT(vector<uint32_t> &vRmE)
{
    //order based method
    uint32_t old_community_size = m_vv_ordered_layers.back().edges.size();
    uint32_t m_ = m_vp_edges.size();
    /*for correctly count the support*/
    vector<bool> be_deleted(m_, false);
    /*process edge by edge*/
    for(uint32_t eid : vRmE){
        m_v_EInfo[eid].bInDT = false;
        be_deleted[eid] = true;
        m_vv_ordered_list[eid].layer_number = 0;
        m_vv_ordered_list[eid].eid = 0;
        vector<pair<uint32_t, uint32_t> > vCE;
        vector<pair<uint32_t, uint32_t> > vFE;
        vector<pair<uint32_t, uint32_t> > v_sub_CE;
        vector<pair<uint32_t, uint32_t> > v_sub_FE;
        PQ m_PQ;
        vector<bool> be_in_pq(m_, false);
        be_in_pq[eid] = true;
        /*find neighbors*/
        uint32_t x = m_vp_edges[eid].first;
        uint32_t y = m_vp_edges[eid].second;
        vCE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_out[y], vCE);
        for(auto neighbors : vCE){
            //if(!be_in_pq[neighbors.first] && !be_deleted[neighbors.first] /*&& m_vv_ordered_list[neighbors.first].layer_number > 0*/){
                m_PQ.push(m_vv_ordered_list[neighbors.first]);
                be_in_pq[neighbors.first] = true;
            //}
            //if(!be_in_pq[neighbors.second]  && !be_deleted[neighbors.second]/*&& m_vv_ordered_list[neighbors.second].layer_number > 0*/ ){
                m_PQ.push(m_vv_ordered_list[neighbors.second]);
                be_in_pq[neighbors.second] = true;
            //}
        }
        vFE.clear();
        findNeib(m_vv_adj_in[x], m_vv_adj_in[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_out[y], vFE);
        findNeib(m_vv_adj_out[x], m_vv_adj_in[y], vFE);
        for(auto neighbors : vFE){
            //if(!be_in_pq[neighbors.first]  && !be_deleted[neighbors.first] /*&& m_vv_ordered_list[neighbors.first].layer_number > 0*/){
                m_PQ.push(m_vv_ordered_list[neighbors.first]);
                be_in_pq[neighbors.first] = true;
            //}
            //if(!be_in_pq[neighbors.second]  && !be_deleted[neighbors.second] /*&& m_vv_ordered_list[neighbors.second].layer_number > 0*/){
                m_PQ.push(m_vv_ordered_list[neighbors.second]);
                be_in_pq[neighbors.second] = true;
            //}
        }
        while (!m_PQ.empty()){
            uint32_t cur_eid =  m_PQ.top().eid;
            m_PQ.pop();
            be_in_pq[cur_eid] = true;
            uint32_t original_layer = m_vv_ordered_list[cur_eid].layer_number;
            /*get the neighbor of cur_eid*/
            x = m_vp_edges[cur_eid].first;
            y = m_vp_edges[cur_eid].second;
            v_sub_CE.clear();
            findNeib(m_vv_adj_in[x], m_vv_adj_out[y], v_sub_CE);
            v_sub_FE.clear();
            findNeib(m_vv_adj_in[x], m_vv_adj_in[y], v_sub_FE);
            findNeib(m_vv_adj_out[x], m_vv_adj_out[y], v_sub_FE);
            findNeib(m_vv_adj_out[x], m_vv_adj_in[y], v_sub_FE);

            /*get its' support cnt in with different layers*/
            vector<uint32_t> support_cnt_c(m_mu + 1,0); // for example, support_cnt_c[1] is the cycle_support in L_{1} U L_{...} U L_[mu]
            vector<uint32_t> support_cnt_f(m_mu + 1,0);
            vector<vector<pair<uint32_t, uint32_t>>> neighbors_of_different_layer(m_mu + 1);
            /*cycle triangle*/
            for (auto atTpE : v_sub_CE)
            {
                if ((be_deleted[atTpE.first]) || (be_deleted[atTpE.second]))
                {
                    /* cannot form triangles */
                    continue;
                }

                for(uint32_t layer_id = 0; layer_id <= COMMON_MIN(m_vv_ordered_list[atTpE.first].layer_number, m_vv_ordered_list[atTpE.second].layer_number) ; ++layer_id) {
                    support_cnt_c[layer_id] += 1;
                    neighbors_of_different_layer[layer_id].push_back(atTpE);
                }
            }
            /*flow triangle*/
            for (auto atTpE : v_sub_FE)
            {
                if ((be_deleted[atTpE.first]) || (be_deleted[atTpE.second]))
                {
                    /* cannot form triangles */
                    continue;
                }
                for(uint32_t layer_id = 0; layer_id <= COMMON_MIN(m_vv_ordered_list[atTpE.first].layer_number, m_vv_ordered_list[atTpE.second].layer_number) ; ++layer_id) {
                    support_cnt_f[layer_id] += 1;
                    neighbors_of_different_layer[layer_id].push_back(atTpE);
                }
            }

            if(support_cnt_c.back() >= m_ui_kc && support_cnt_f.back() >= m_ui_kf){
                continue;
            }

            /*the cur_eid is not in the L_{mu}, otherise it will not affect its neighbors*/
            else{
                removeItem(m_vv_ordered_layers[original_layer].edges,cur_eid);
                uint32_t layer_lower_bound = 0 , layer_upper_bound = 0;

                /*new layers is formed*/
//                if(support_cnt_c[m_mu - 1] >= m_ui_kc && support_cnt_f[m_mu - 1] >= m_ui_kf){
//                    cout << "test 7.2" << endl;
//                    insertNewLayer(cur_eid);
//                    //layer_lower_bound = COMMON_MIN(m_vv_ordered_list[cur_eid].layer_number, 0);
//                    layer_lower_bound = COMMON_MAX(m_vv_ordered_list[cur_eid].layer_number, 1);
//                    layer_upper_bound = COMMON_MIN(original_layer + 1 , m_mu - 1);
//                }
//                else{    //find layer to insert

//                uint32_t test = 0;
//                bool flag = verify(original_layer,cur_eid,test);
//                if(!flag){
//                    //ASSERT_MSG(test <= original_layer , "verify wrong " << original_layer << " " << cur_eid << " " << test << " " << m_mu);
//                }
                findLayerToInsert(support_cnt_c, support_cnt_f, cur_eid);
                // layer_lower_bound = COMMON_MIN(m_vv_ordered_list[cur_eid].layer_number, 0);
                layer_lower_bound = COMMON_MAX(m_vv_ordered_list[cur_eid].layer_number, 1);
                layer_upper_bound = COMMON_MIN(original_layer + 1 , m_mu) ;
//                }

                if(original_layer <= m_vv_ordered_list[cur_eid].layer_number) continue;
                ASSERT_MSG(m_vv_ordered_list[cur_eid].layer_number <= original_layer, "layer info wrong " << " eid " << cur_eid << " " << m_vv_ordered_list[cur_eid].layer_number << " " << original_layer);

                insertNeighborsToQueue(be_in_pq, m_PQ, neighbors_of_different_layer,
                                       layer_lower_bound, layer_upper_bound);
            }

        }


        //compress the updated layers
        //compressUpdatedLayers();
    }


    /*the new L_{mu} is not D-truss*/
    uint32_t uiRmCnt = 0;
    if(m_vv_ordered_layers.back().edges.size() < old_community_size){
        uiRmCnt =  old_community_size - m_vv_ordered_layers.back().edges.size();
    }
    else{
        uiRmCnt =  old_community_size ;
    }

    return uiRmCnt;
}

/**
    remove edges
**/
void Order3::remove(vector<pair<uint32_t, uint32_t> > & vEdges)
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
    g_tUptRm_order_3 += end - beg;

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

        /*remove order information*/
        m_vv_ordered_list[atEid].layer_number = 0;
        m_vv_ordered_list[atEid].eid = 0;
        //m_vv_ordered_list[atEid].prev = m_vv_ordered_list[atEid].next = -1;
        //m_vv_ordered_list[atEid].layer_number = 0;

    }

    m_uiEdgeCnt = m_vp_edges.size();
    //printf("RM %d edges removed from D-truss\n", uiRmCnt);
}


/**
verify the corectness of the laryer_number
**/
bool Order3::verify(uint32_t peel_layer_number, uint32_t eid,  uint32_t &layer_number_test){
    uint32_t x = 0, y = 0;
    vector<pair<uint32_t,uint32_t>> vCE, vFE;

    x = m_vp_edges[eid].first;
    y = m_vp_edges[eid].second;
    vCE.clear();
    findNeib(m_vv_adj_in[x], m_vv_adj_out[y], vCE);
    vFE.clear();
    findNeib(m_vv_adj_in[x], m_vv_adj_in[y], vFE);
    findNeib(m_vv_adj_out[x], m_vv_adj_out[y], vFE);
    findNeib(m_vv_adj_out[x], m_vv_adj_in[y], vFE);

    /*get its' support cnt in with different layers*/
    vector<uint32_t> support_cnt_c(m_mu + 1,0); // for example, support_cnt_c[1] is the cycle_support in G \ L_{1}
    vector<uint32_t> support_cnt_f(m_mu + 1,0);
    vector<vector<pair<uint32_t, uint32_t>>> neighbors_of_different_layer(m_mu + 1);
    /*cycle triangle*/
    for (auto atTpE : vCE)
    {
        for(uint32_t layer_id = 0; layer_id <= COMMON_MIN(m_vv_ordered_list[atTpE.first].layer_number, m_vv_ordered_list[atTpE.second].layer_number) ; ++layer_id) {
            support_cnt_c[layer_id] += 1;
            neighbors_of_different_layer[layer_id].push_back(atTpE);
        }
    }
    /*flow triangle*/
    for (auto atTpE : vFE)
    {
        for(uint32_t layer_id = 0; layer_id <= COMMON_MIN(m_vv_ordered_list[atTpE.first].layer_number, m_vv_ordered_list[atTpE.second].layer_number) ; ++layer_id) {
            support_cnt_f[layer_id] += 1;
            neighbors_of_different_layer[layer_id].push_back(atTpE);
        }
    }

    /*find the layer number*/
    uint32_t layer_id = 0;
    if(support_cnt_c[0] < m_ui_kc || support_cnt_f[0] < m_ui_kf){
        layer_id = -1;
    }
    else if(support_cnt_c[m_mu] >= m_ui_kc && support_cnt_f[m_mu] >= m_ui_kf){
        layer_id = m_mu - 1;
    }
    else{
        for(int i = 0; i < m_mu; i++){
            if((support_cnt_c[i] >= m_ui_kc && support_cnt_f[i] >= m_ui_kf) &&
               (support_cnt_c[i + 1] < m_ui_kc || support_cnt_f[i + 1] < m_ui_kf)){
                layer_id = i;
                break;
            }
        }
    }
    uint32_t layer_number = 0;
    layer_number = layer_id + 1;
    layer_number_test = layer_number;
    if(peel_layer_number != layer_number) {
        return false;
    }
    return true;
}

/**
   get D-truss by peeling
**/
void Order3::peeling(vector<uint32_t> &vEdges, vector<uint32_t> &vFixE, vector<uint32_t> &vResE)
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
                if (A_in[v][pv].pid == A_out[u][pu].pid)
                {
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
    vector<uint32_t> edge_orders;  //record the order of edge deletions
    uint32_t layer_number_cnt = 0;
    vector<bool> removed(m_, false);
    vector<bool> vWaitRmFlag(m_, false);
    uint32_t uiRmCnt = 0;
    uint32_t uiENum = 0;
    vector<uint32_t> vWait(m_);
    vector<uint32_t> vWaitRmCache(m_);
    // init cache
    vWaitRmCache.clear();
    vector<uint32_t> edges_in_this_layer;
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

            edges_in_this_layer.push_back(eid);

            m_vv_ordered_list[eid].layer_number = layer_number_cnt;
            edge_orders.push_back(eid);

        }
        else if (v_Sup_f[eid] < m_ui_kf)
        {
            vWaitRmFlag[eid] = true;
            vWaitRmCache.push_back(eid);

            edges_in_this_layer.push_back(eid);

            m_vv_ordered_list[eid].layer_number = layer_number_cnt;
            edge_orders.push_back(eid);

        }
    }
    if(!edges_in_this_layer.empty()){
        stLayerInfo layer = {edges_in_this_layer, (int)layer_number_cnt/*, (int)edges_in_this_layer.front()*/};
        m_vv_ordered_layers.push_back(layer);
    }
    edges_in_this_layer.clear();
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
        ++layer_number_cnt;
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

                        edges_in_this_layer.push_back(e1);
                        m_vv_ordered_list[e1].layer_number = layer_number_cnt;
                        edge_orders.push_back(e1);
                    }
                }
                if ((!vWaitRmFlag[e2]) && (vbTarget[e2]))
                {
                    --v_Sup_f[e2];
                    if (v_Sup_f[e2] < m_ui_kf)
                    {
                        vWaitRmFlag[e2] = true;
                        vWaitRmCache.push_back(e2);

                        edges_in_this_layer.push_back(e1);
                        m_vv_ordered_list[e2].layer_number = layer_number_cnt;
                        edge_orders.push_back(e2);
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

                        edges_in_this_layer.push_back(e1);
                        m_vv_ordered_list[e1].layer_number = layer_number_cnt;
                        edge_orders.push_back(e1);
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

                        edges_in_this_layer.push_back(e1);
                        m_vv_ordered_list[e2].layer_number = layer_number_cnt;
                        edge_orders.push_back(e2);
                    }
                }
            }
        }
        if(!edges_in_this_layer.empty()){
            stLayerInfo layer = {edges_in_this_layer, (int)layer_number_cnt/*, (int)edges_in_this_layer.front()*/};
            m_vv_ordered_layers.push_back(layer);
        }
        edges_in_this_layer.clear();
    }

    // keep result
    for (uint32_t eid : vEdges)
    {
        if (!removed[eid])
        {
            vResE.push_back(eid);
            m_mu = layer_number_cnt;
            edges_in_this_layer.push_back(eid);
            m_vv_ordered_list[eid].layer_number = layer_number_cnt;
            edge_orders.push_back(eid);
        }
    }
    if(!edges_in_this_layer.empty()){
        stLayerInfo layer = {edges_in_this_layer, (int)layer_number_cnt/*, (int)edges_in_this_layer.front()*/};
        m_vv_ordered_layers.push_back(layer);
    }
    edges_in_this_layer.clear();


    //assign labels to edges, do the link among edges / layers
    ASSERT_MSG(edge_orders.size() == m_vp_edges.size(), "some edges are not considered in the peeling process " << edge_orders.size() << " " << m_vp_edges.size());
    /*label the edges, link them*/
    //label_t label = 0;
    for(uint32_t i = 1; i < edge_orders.size()  - 1; ++i){
        uint32_t eid = edge_orders[i];
        //label += INIT_TAG_GAP_2;
        ASSERT(m_vv_ordered_list[eid].layer_number >= 0);
        m_vv_ordered_list[eid].eid = eid;
        //m_vv_ordered_list[eid].label = label;
        //m_vv_ordered_list[eid].prev = edge_orders[i - 1] ;
        //m_vv_ordered_list[eid].next = edge_orders[i + 1] ;
        //cout << "eid: " << eid << " " << m_v_nodes[m_vp_edges[eid].first] << " " << m_v_nodes[m_vp_edges[eid].second] << " " << m_vv_ordered_list[i].label << endl;
    }
    m_vv_ordered_list[edge_orders[0]].eid = edge_orders[0];
    //m_vv_ordered_list[edge_orders[0]].prev = -1;
    //m_vv_ordered_list[edge_orders[0]].next = edge_orders[1];
    //m_vv_ordered_list[edge_orders[0]].label = 0;

    m_vv_ordered_list[edge_orders.back()].eid = edge_orders.back();
    //m_vv_ordered_list[edge_orders.back()].prev = edge_orders[edge_orders.size() - 2];
    //m_vv_ordered_list[edge_orders.back()].next = -1;
    //m_vv_ordered_list[edge_orders.back()].label = MAX_LABEL;

    //printf the orderd list
//    for(uint32_t i = 0; i < m_vv_ordered_list.size(); ++i){
//        cout << "eid: " << m_vv_ordered_list[i].eid << " " << m_v_nodes[m_vp_edges[m_vv_ordered_list[i].eid].first] << " " << m_v_nodes[m_vp_edges[m_vv_ordered_list[i].eid].second] << " " << m_vv_ordered_list[i].label <<
//             " prev: " <<  m_vv_ordered_list[i].prev << " next: " << m_vv_ordered_list[i].next  << endl;
//    }

    //verify the layer_number
    // for(uint32_t i = 0; i < m_vv_ordered_list.size();i++){
    //     ASSERT(i == m_vv_ordered_list[i].eid);
    //     uint32_t layer_number_test = 0;
    //     bool flag = verify(m_vv_ordered_list[i].layer_number,i, layer_number_test);
    //     ASSERT_MSG(flag, "verify failed, eid: " << i << " " << m_vv_ordered_list[i].layer_number << " " << layer_number_test);
    //     //ASSERT_MSG(i!=1780, "eid: " << i << " " << m_vv_ordered_list[i].layer_number << " " << layer_number_test);
    // }

    const auto end = std::chrono::steady_clock::now();
    const auto dif = end - beg;
    //peeling_time += end - beg;
    /*printf("Peeling costs \x1b[1;31m%f\x1b[0m ms.\n",
           std::chrono::duration<double, std::milli>(dif).count());*/
    //printf("DEBUG index done last: (%d, %d)\n", m_vp_edges[m_ - 1].first, m_vp_edges[m_ - 1].second);
}
/**
    find neighbor, at the same time find in which layer the edge should be moved to
**/
void Order3::findNeibAndLayer(vector<Order3::AdjEntry> &vAdj1, vector<Order3::AdjEntry> &vAdj2, vector<pair<uint32_t, uint32_t>> &vTrdE,
                              vector<uint32_t> &support_cnt) {
    size_t p1 = 0, p2 = 0;
    while (p1 < vAdj1.size() && p2 < vAdj2.size()) {
        if (vAdj1[p1].pid == vAdj2[p2].pid) {
            // record the support count in the sub-graph starting from each layer
            for(uint32_t layer_id = 0; layer_id <= COMMON_MIN(m_vv_ordered_list[vAdj1[p1].eid].layer_number, m_vv_ordered_list[vAdj2[p2].eid].layer_number) ; ++layer_id) {
                support_cnt[layer_id] += 1;
            }
            vTrdE.push_back(pair<uint32_t, uint32_t>(vAdj1[p1].eid, vAdj2[p2].eid));
            ++p1; ++p2;
        } else if (vAdj1[p1].pid < vAdj2[p2].pid) {
            ++p1;
        } else {
            ++p2;
        }
    }
}
