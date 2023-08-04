
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <utility>
#include <queue>
    

#include "Update.h"

chrono::duration<double> g_tUptIns;
chrono::duration<double> g_tUptRm;
int SearchSpace = 0;
/**
 Find the cycle/flow neighbr & assign value for lifetime vector
**/
void Update::findNeibAndAssignLifeT(vector<AdjEntry> &vAdj1, vector<AdjEntry> &vAdj2, vector<pair<uint32_t, uint32_t> > & vTrdE,
const bool &flag ,const int &iStart, const int &iStep, const uint32_t &nEid)
{
    size_t p1 = 0, p2 = 0;
    int vSize = m_v_EInfo[nEid].lifeTime.size();
    while (p1 < vAdj1.size() && p2 < vAdj2.size()) { //cycle trianle
        if (vAdj1[p1].pid == vAdj2[p2].pid) {
        //if(m_v_EInfo[vAdj1[p1].eid].t < m_v_EInfo[nEid].t  && m_v_EInfo[vAdj2[p2].eid].t < m_v_EInfo[nEid].t)
        //{
        vTrdE.push_back(pair<uint32_t, uint32_t>(vAdj1[p1].eid, vAdj2[p2].eid));
        //}
        uint32_t smallerTime = COMMON_MIN(m_v_EInfo[vAdj1[p1].eid].t, m_v_EInfo[vAdj2[p2].eid].t);
        if(!flag){
            /*update lifetime cycle support for the new edge*/
            
            for(int i = 0; i < COMMON_RUP((smallerTime - iStart), iStep); i++)
            {       
                ++m_v_EInfo[nEid].lifeTime[i].first;
            }
            /*update lifetime cycle support for the old edge*/
           
            for(int i = vSize - 1; i >= vSize - COMMON_RUP((m_v_EInfo[vAdj1[p1].eid].t - iStart), iStep); i--) 
            {       
                ++m_v_EInfo[vAdj1[p1].eid].lifeTime[i].first;
            }
           
            for(int i = vSize - 1; i >= vSize - COMMON_RUP((m_v_EInfo[vAdj2[p2].eid].t - iStart), iStep); i--) 
            {       
                ++m_v_EInfo[vAdj2[p2].eid].lifeTime[i].first;
            }
        }
        else{
            /*update lifetime flow support for the new edge*/
            for(int i = 0; i < COMMON_RUP((smallerTime - iStart), iStep); i++) 
            {       
                ++m_v_EInfo[nEid].lifeTime[i].second;
            }
            /*update lifetime flow support for the old edge*/
            for(int i = vSize - 1; i >= vSize - COMMON_RUP((m_v_EInfo[vAdj1[p1].eid].t - iStart), iStep); i--) 
            {       
                ++m_v_EInfo[vAdj1[p1].eid].lifeTime[i].second;
            }
            for(int i = vSize - 1; i >= vSize - COMMON_RUP((m_v_EInfo[vAdj2[p2].eid].t - iStart), iStep); i--) 
            {       
                ++m_v_EInfo[vAdj2[p2].eid].lifeTime[i].second;
            }
        }
        ++p1; ++p2;
        } else if (vAdj1[p1].pid < vAdj2[p2].pid) {
        ++p1;
        } else {
        ++p2;
        }
    }
}
/**
   assign value for lifetime vector
   nEid: the new edge eid
   flag = 0: cycle triangle, = 1: flow triangle
**/
void Update::assignLifeT(vector<AdjEntry> &vAdj1, vector<AdjEntry> &vAdj2, const bool &flag ,const int &iStart, const int &iStep, const uint32_t &nEid)
{
    size_t p1 = 0, p2 = 0;
    int vSize = m_v_EInfo[nEid].lifeTime.size();
    while (p1 < vAdj1.size() && p2 < vAdj2.size()) { //cycle trianle
        if (vAdj1[p1].pid == vAdj2[p2].pid) {
        //vTrdE.push_back(pair<uint32_t, uint32_t>(vAdj1[p1].eid, vAdj2[p2].eid));
        uint32_t smallerTime = COMMON_MIN(m_v_EInfo[vAdj1[p1].eid].t, m_v_EInfo[vAdj2[p2].eid].t);
        if(!flag){
            /*update lifetime cycle support for the new edge*/
            for(int i = 0; i < COMMON_RUP((smallerTime - iStart), iStep); i++) 
            {       
                ++m_v_EInfo[nEid].lifeTime[i].first;
            }
            /*update lifetime cycle support for the old edge*/
            for(int i = vSize - 1; i >= vSize - COMMON_RUP((m_v_EInfo[vAdj1[p1].eid].t - iStart), iStep); i--) 
            {       
                ++m_v_EInfo[vAdj1[p1].eid].lifeTime[i].first;
            }
            for(int i = vSize - 1; i >= vSize - COMMON_RUP((m_v_EInfo[vAdj2[p2].eid].t - iStart), iStep); i--) 
            {       
                ++m_v_EInfo[vAdj2[p2].eid].lifeTime[i].first;
            }
        }
        else{
            /*update lifetime flow support for the new edge*/
            for(int i = 0; i < COMMON_RUP((smallerTime - iStart), iStep); i++) 
            {       
                ++m_v_EInfo[nEid].lifeTime[i].second;
            }
            /*update lifetime flow support for the old edge*/
            for(int i = vSize - 1; i >= vSize - COMMON_RUP((m_v_EInfo[vAdj1[p1].eid].t - iStart), iStep); i--) 
            {       
                ++m_v_EInfo[vAdj1[p1].eid].lifeTime[i].second;
            }
            for(int i = vSize - 1; i >= vSize - COMMON_RUP((m_v_EInfo[vAdj2[p2].eid].t - iStart), iStep); i--) 
            {       
                ++m_v_EInfo[vAdj2[p2].eid].lifeTime[i].second;
            }
        }
        ++p1; ++p2;
        } else if (vAdj1[p1].pid < vAdj2[p2].pid) {
        ++p1;
        } else {
        ++p2;
        }
    }
}
/**
    find nos-neighbor
**/
void Update::findNosNeib(vector<AdjEntry> &vAdj1, vector<AdjEntry> &vAdj2, vector<pair<uint32_t, uint32_t> > & vTrdE, u_int32_t &edgeT)
{
  size_t p1 = 0, p2 = 0;
  while (p1 < vAdj1.size() && p2 < vAdj2.size()) {
    if (vAdj1[p1].pid == vAdj2[p2].pid) {
      if(m_v_EInfo[vAdj1[p1].eid].t < edgeT && m_v_EInfo[vAdj2[p2].eid].t < edgeT)
      {
      vTrdE.push_back(pair<uint32_t, uint32_t>(vAdj1[p1].eid, vAdj2[p2].eid));
      }
      ++p1; ++p2;
    } else if (vAdj1[p1].pid < vAdj2[p2].pid) {
      ++p1;
    } else {
      ++p2;
    }
  }
}
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
    init
**/
Update::Update(vector<ST_EDGE_ENTRY> & vEdges, uint32_t uiMaxN, uint32_t uiMaxM, uint32_t uiKc, uint32_t uiKf, int &iWd, int &iStep,int &iStart)
{
    m_ui_n = uiMaxN;
    m_ui_m = uiMaxM;
    m_ui_kc = uiKc;
    m_ui_kf = uiKf;
    m_ui_wd = iWd;
    m_ui_step = iStep;

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
        bool bRes = addEInfo(atE.x, atE.y, atE.t, &eid);
        /*printf("NEW get (%d, %d) rename (%d, %d) bool: %d\n",
               atE.first, atE.second,
               m_vp_edges[eid].first, m_vp_edges[eid].second,
               bRes);*/
        if (bRes)
        {
            vTargetE.push_back(eid);
            /*initialize lifetime vector for initial window*/
        }
        /*else
        {
            printf("NEW in it (%d, %d)\n",
                   m_v_nodes[m_vp_edges[eid].first], m_v_nodes[m_vp_edges[eid].second]);
        }*/
    }
    /*assign value for edges' lifetime vector*/
   
    for(uint32_t eid : vTargetE)
    {
        assignLifeT(m_vv_adj_in[m_vp_edges[eid].first], m_vv_adj_out[m_vp_edges[eid].second], false, iStart , iStep, eid);

        assignLifeT(m_vv_adj_in[m_vp_edges[eid].first], m_vv_adj_in[m_vp_edges[eid].second], true, iStart , iStep, eid);
        assignLifeT(m_vv_adj_out[m_vp_edges[eid].first], m_vv_adj_out[m_vp_edges[eid].second], true, iStart , iStep, eid);
        assignLifeT(m_vv_adj_out[m_vp_edges[eid].first], m_vv_adj_in[m_vp_edges[eid].second], true, iStart , iStep, eid);
    }

    vector<uint32_t> vFixE;
    vector<uint32_t> vResE;
    peeling(vTargetE, vFixE, vResE);
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
bool Update::addEInfo(uint32_t x, uint32_t y, uint32_t t, uint32_t *puiEid)
{
    /* rename */
    uint32_t u = renamePid(x);
    uint32_t v = renamePid(y);
    if (findEid(u, v, puiEid))
    {
        /* not new */
        ++m_v_EInfo[*puiEid].cnt;
        m_v_EInfo[*puiEid].t = t;
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
    m_v_EInfo[*puiEid].t = t;
    m_v_EInfo[*puiEid].bInDT = false;
    m_v_EInfo[*puiEid].bUsed = true;
    vector<pair<uint32_t,uint32_t>> tmpLifeTime(COMMON_RUP(m_ui_wd,m_ui_step));
    m_v_EInfo[*puiEid].lifeTime = tmpLifeTime;

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
uint32_t Update::insertDT(vector<uint32_t> &vInsE, int &iStart)
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
        //findNosNeib(m_vv_adj_in[x], m_vv_adj_out[y], vCE, m_v_EInfo[eid].t);
        findNeibAndAssignLifeT(m_vv_adj_in[x], m_vv_adj_out[y], vCE, false, iStart, m_ui_step, eid);
        uint32_t uiCSup = vCE.size();
        //printf("BFS queue get CSup: %d\n", uiCSup);
        if (uiCSup < m_ui_kc)
        {
            /* cannot in D-truss */
            continue;
        }
        /* flow support */
        vFE.clear();
        //findNosNeib(m_vv_adj_in[x], m_vv_adj_in[y], vFE, m_v_EInfo[eid].t);
        //findNosNeib(m_vv_adj_out[x], m_vv_adj_out[y], vFE, m_v_EInfo[eid].t);
        //findNosNeib(m_vv_adj_out[x], m_vv_adj_in[y], vFE, m_v_EInfo[eid].t);
        findNeibAndAssignLifeT(m_vv_adj_in[x], m_vv_adj_in[y], vFE, true, iStart, m_ui_step, eid);
        findNeibAndAssignLifeT(m_vv_adj_out[x], m_vv_adj_out[y], vFE, true, iStart, m_ui_step, eid);
        findNeibAndAssignLifeT(m_vv_adj_out[x], m_vv_adj_in[y], vFE, true, iStart, m_ui_step, eid);
        uint32_t uiFSup = vFE.size();
        //printf("BFS queue get FSup: %d\n", uiFSup);
        if (uiFSup < m_ui_kf)
        {
            /* cannot in D-truss */
            continue;
        }
        /*update the lifetime vector for new edges and its neighbor edges, can be paralleled?*/
        /*can be intergereted with nos-support check since they all involve cycle/flow neighbor accessing*/
        // assignLifeT(m_vv_adj_in[m_vp_edges[eid].first], m_vv_adj_out[m_vp_edges[eid].second], false, iStart, iStep, eid);

        // assignLifeT(m_vv_adj_in[m_vp_edges[eid].first], m_vv_adj_in[m_vp_edges[eid].second], true, iStart, iStep, eid);
        // assignLifeT(m_vv_adj_out[m_vp_edges[eid].first], m_vv_adj_out[m_vp_edges[eid].second], true, iStart, iStep, eid);
        // assignLifeT(m_vv_adj_out[m_vp_edges[eid].first], m_vv_adj_in[m_vp_edges[eid].second], true, iStart, iStep, eid);
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


    for(uint32_t uiCurEid = 0; uiCurEid < m_ ; uiCurEid++)
    {
        if(m_v_EInfo[uiCurEid].bUsed)
        {
            if (m_v_EInfo[uiCurEid].bInDT)
            {
                vDetE.push_back(uiCurEid);
                vbInPool[uiCurEid] = true;
                continue;
            }
            /*check support*/

            int index = m_v_EInfo[uiCurEid].lifeTime.size() - COMMON_RUP((m_v_EInfo[uiCurEid].t - iStart),m_ui_step);
            if(index < 0 || index > m_v_EInfo[uiCurEid].lifeTime.size())cout << index << " " << m_v_EInfo[uiCurEid].t << " " << iStart <<  " " << m_v_EInfo[uiCurEid].lifeTime.size() << " " << COMMON_RUP((m_v_EInfo[uiCurEid].t - iStart),m_ui_step) << endl;
            if((m_v_EInfo[uiCurEid].lifeTime[index].first < m_ui_kc)||(m_v_EInfo[uiCurEid].lifeTime[index].second < m_ui_kf))
            {
                continue;
            }
            vTargetE.push_back(uiCurEid);
            vbInPool[uiCurEid] = true;
        }
    }
     

    // while (!lsQ.empty())
    // {
    //     uint32_t uiCurEid = lsQ.front();
    //     lsQ.pop_front();
    //     uint32_t x = m_vp_edges[uiCurEid].first;
    //     uint32_t y = m_vp_edges[uiCurEid].second;

    //     //printf("BFS queue get original (%d, %d)\n", m_v_nodes[x], m_v_nodes[y]);
    //     //printf("BFS queue get (%d, %d) in it: %d\n", x, y, m_v_EInfo[uiCurEid].bInDT);

    //     if (m_v_EInfo[uiCurEid].bInDT)
    //     {
    //         vDetE.push_back(uiCurEid);
    //         vbInPool[uiCurEid] = true;
    //         continue;
    //     }
    //     /* check support */
    //     int index = m_v_EInfo[uiCurEid].lifeTime.size() - COMMON_RUP((m_v_EInfo[uiCurEid].t - iStart),m_ui_step);
    //     if((m_v_EInfo[uiCurEid].lifeTime[index].first < m_ui_kc)||(m_v_EInfo[uiCurEid].lifeTime[index].second < m_ui_kf))
    //     {
    //         continue;
    //     }
    //     /* may be in the D-truss */
    //     vTargetE.push_back(uiCurEid);
    //     vbInPool[uiCurEid] = true;

    //     for (auto atTpE : vCE)
    //     {
    //         for (auto atTpEid : {atTpE.first, atTpE.second})
    //         {
    //             if (!vbInQ[atTpEid])
    //             {
    //                 lsQ.push_back(atTpEid);
    //                 vbInQ[atTpEid] = true;
    //             }
    //         }
    //     }
    //     for (auto atTpE : vFE)
    //     {
    //         for (auto atTpEid : {atTpE.first, atTpE.second})
    //         {
    //             if (!vbInQ[atTpEid])
    //             {
    //                 lsQ.push_back(atTpEid);
    //                 vbInQ[atTpEid] = true;
    //             }
    //         }
    //     }
    // }
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
void Update::add(vector<ST_EDGE_ENTRY> & vEdges, int &iStart)
{
    vector<uint32_t> vInsE;
    vInsE.reserve(vEdges.size());
   
    for (auto atE : vEdges)
    {
        uint32_t eid = 0;
        bool bRes = addEInfo(atE.x, atE.y, atE.t, &eid);
        if (bRes)
        {
            vInsE.push_back(eid);
        }
    }
   
    const auto beg = std::chrono::steady_clock::now();
    int realStartTime = iStart - m_ui_wd;
    m_uiDTCnt += insertDT(vInsE, /*realStartTime*/iStart );
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
    add edges
**/
void Update::save(vector<pair<uint32_t, uint32_t> > & vEdges)
{
     printf("Total search space: %d\n", SearchSpace);
    printf("DEBUG save total size: %d\n", m_vp_edges.size());
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
void Update::remove(vector<ST_EDGE_ENTRY> & vEdges)
{
    uint32_t eid = 0;

    vector<uint32_t> vRmE;
    vRmE.reserve(vEdges.size());

    for (auto atE : vEdges)
    {
        /* rename */
        uint32_t u = renamePid(atE.x);
        uint32_t v = renamePid(atE.y);
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

    // const auto beg = std::chrono::steady_clock::now();
    // uint32_t uiRmCnt = rmDT(vRmE);
    // const auto end = std::chrono::steady_clock::now();
    // g_tUptRm += end - beg;

    m_uiDTCnt -= /*uiRmCnt*/vRmE.size();
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

    // initialize adjacency arrays
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
  /*printf("Peeling costs \x1b[1;31m%f\x1b[0m ms.\n",
         std::chrono::duration<double, std::milli>(dif).count());*/
    //printf("DEBUG index done last: (%d, %d)\n", m_vp_edges[m_ - 1].first, m_vp_edges[m_ - 1].second);
}
