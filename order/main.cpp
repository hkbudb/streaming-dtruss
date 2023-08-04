/*
 * This is the BFS-based algorithm.
 *
 * 1.For inserted edges, we use BFS to collect the adjacent edges with (cupp >= kc and fupp >= kf) the
 *   new edges, then we peel the collected result and merge it with the original Dtruss
 *
 * 2. For deleted edges, we use a verify algorithm to delete the disqualified edges in the original
 *    community
 */

#include "common.h"
#include "Stream.h"
#include "Update.h"
#include "Order.h"
#include "Order2.h"
#include "Order3.h"

extern std::chrono::duration<double> g_tUptIns_order_3;
extern std::chrono::duration<double> g_tUptRm_order_3;
extern std::chrono::duration<double> g_tUptIns_order;
extern std::chrono::duration<double> g_tUptRm_order;
extern std::chrono::duration<double> g_tUptIns;
extern std::chrono::duration<double> g_tUptRm;
extern std::chrono::duration<double> peeling_time;
extern ::uint32_t total_pruned_size, total_unpruned_size;
extern long long int MAX_LABEL;
extern long int INIT_TAG_GAP;

int main(int argc, char** argv) {
    // read the graph and peel D-truss
    //ASSERT_MSG(5 < argc, "wrong parameters");
    int iWd = 9000;
    int iStep = 900;
    int iKC = 1;
    int iKF = 1;
    string Q;
    char *filepath = "/Users/xkliao/CLionProjects/StreamingNew/data/sx.txt";
    sscanf(argv[2], "%d", &iWd);
    sscanf(argv[3], "%d", &iStep);
    sscanf(argv[4], "%d", &iKC);
    sscanf(argv[5], "%d", &iKF);
    sscanf(argv[6], "%s", &Q);

    Stream oGraph;
    oGraph.readFromFile(argv[1]/*filepath*/);

    printf("Graph info: [%d, %d], n: %d m: %d\n", oGraph.m_uiTS, oGraph.m_uiTE,
           oGraph.m_uiN, oGraph.m_uiM);

    printf("Parameters: %d, %d, %d, %d\n", iWd, iStep,
           iKC, iKF);


    const auto beg = std::chrono::steady_clock::now();
    vector<pair<uint32_t, uint32_t> > vWEdges;
    /*int iEndT = oGraph.m_uiTS;
    int iStartT = iEndT - iWd + 1;
    iStartT = iStartT > 0? iStartT: 0;*/

    int iEndT = oGraph.m_uiTS + iWd - 1;

    int iStartT = oGraph.m_uiTS;
    iStartT = iStartT > 0? iStartT: 0;
    int iBatch = iStep;
    
    int CurrentTime = 0;
    oGraph.getEdges(iStartT, iEndT, vWEdges);
    //printf("get edges: %d from [%d, %d]\n", vWEdges.size(), iStartT, iEndT);
    //Update oUpdateG(vWEdges, oGraph.m_uiN, oGraph.m_uiM, iKC, iKF);
    Order3 oOrderG(vWEdges, oGraph.m_uiN, oGraph.m_uiM, iKC, iKF);

    vector<pair<uint32_t, uint32_t>> vResE;
    oOrderG.save(vResE);
    std::sort(vResE.begin(), vResE.end());
    oGraph.writeToFile("../out-order-initial.txt", vResE);
    vResE.clear();

    vResE.clear();
    //oUpdateG.save(vResE);
    //std::sort(vResE.begin(), vResE.end());
    oGraph.writeToFile("../out-peeling-initial.txt", vResE);
    vResE.clear();

    while (iEndT < oGraph.m_uiTE)
    {
        int ProcssedEdgeNum = 0;
        const auto preTime1 = g_tUptRm_order_3, preTime2 = g_tUptIns_order_3;  //for throughput calculation
        iBatch = COMMON_MIN(iBatch, oGraph.m_uiTE - iEndT);
        int iRmEnd = iStartT + iBatch - 1;
        vWEdges.clear();
        oGraph.getEdges(iStartT, iRmEnd, vWEdges);
        //oUpdateG.remove(vWEdges);
        oOrderG.remove(vWEdges);
        iStartT += iBatch;
        ProcssedEdgeNum += /*vWEdges.size()*/oOrderG.m_uiEdgeCnt;


        int iInsEnd = iEndT + iBatch;
        vWEdges.clear();
        oGraph.getEdges(iEndT + 1, iInsEnd, vWEdges);
        //printf("get edges: %d from [%d, %d]\n", vWEdges.size(), iEndT + 1, iInsEnd);
        //oUpdateG.add(vWEdges);
        oOrderG.add(vWEdges);
        iEndT = iInsEnd;
        ProcssedEdgeNum+= /*vWEdges.size()*/oOrderG.m_uiEdgeCnt;
        const auto endTime1 = g_tUptRm_order_3, endTime2 = g_tUptIns_order_3;  //for throughput calculation
        //if( CurrentTime % (((int) oGraph.m_uiTE / iBatch) / 20 + 1) == 0) 
        //{
              //cout << CurrentTime*iBatch << " " << (long) 1000*ProcssedEdgeNum / (std::chrono::duration<double, std::milli>(endTime1 - preTime1).count() + std::chrono::duration<double, std::milli>(endTime2 - preTime2).count()) << endl;
        //}
        CurrentTime++;

        //break;
    }
    const auto end = std::chrono::steady_clock::now();
    const auto dif = end - beg;


    //vertification
    //vector<pair<uint32_t, uint32_t>> vResE;
    vResE.clear();
    //oUpdateG.save(vResE);
    std::sort(vResE.begin(), vResE.end());
    oGraph.writeToFile("../out-bfs.txt", vResE);
    vResE.clear();

    vResE.clear();
    oOrderG.save(vResE);
    std::sort(vResE.begin(), vResE.end());
    oGraph.writeToFile("../out-order.txt", vResE);
    vResE.clear();

    vWEdges.clear();
    oGraph.getEdges(iStartT, iEndT, vWEdges);
    Update oTpPeel(vWEdges, oGraph.m_uiN, oGraph.m_uiM, iKC, iKF);
    oTpPeel.save(vResE);
    std::sort(vResE.begin(), vResE.end());
    oGraph.writeToFile("../out-repeel.txt", vResE);


    //print the community quality of the final community
    //pair<double,double> DtQuality = oUpdateG.getDtQuality(vResE);
    //printf("Save size: %d CMSin: %f CMSout: %f\n", vResE.size() , DtQuality.first, DtQuality.second);

    //print running time
    printf("Query costs \x1b[1;31m%f\x1b[0m ms.\n",
           std::chrono::duration<double, std::milli>(dif).count());

    printf("peeling_time: \x1b[1;31m%f\x1b[0m ms.\n",
           std::chrono::duration<double, std::milli>(peeling_time).count());

    printf("g_tUptIns: \x1b[1;31m%f\x1b[0m ms.\n",
           std::chrono::duration<double, std::milli>(g_tUptIns).count());

    printf("g_tUptRm: \x1b[1;31m%f\x1b[0m ms.\n",
           std::chrono::duration<double, std::milli>(g_tUptRm).count());

    printf("g_tUptIns+g_tUptRm: %.3f ms\n",
           std::chrono::duration<double, std::milli>(g_tUptIns+g_tUptRm).count());

    printf("order_g_tUptIns: \x1b[1;31m%f\x1b[0m ms.\n",
           std::chrono::duration<double, std::milli>(g_tUptIns_order_3).count());

    printf("order_g_tUptRm: \x1b[1;31m%f\x1b[0m ms.\n",
           std::chrono::duration<double, std::milli>(g_tUptRm_order_3).count());

    printf("order_g_tUptIns+order_g_tUptRm: %.3f ms\n",
           std::chrono::duration<double, std::milli>(g_tUptIns_order_3+g_tUptRm_order_3).count());

    cout << "pruning efficiency: " << 100 * ((double)(1.0 -  1.0 * total_pruned_size/total_unpruned_size))
        << "% " << " insert original search size： " << total_unpruned_size << endl;



    if (7 == argc)
    {
        vector<pair<uint32_t, uint32_t> > vResE;
        printf("Save %s\n", argv[6]/*"out.txt"*/);
        vWEdges.clear();
        oGraph.getEdges(iStartT, iEndT, vWEdges);
        printf("get edges: %d from [%d, %d]\n", vWEdges.size(), iStartT, iEndT);
        Update oTpPeel(vWEdges, oGraph.m_uiN, oGraph.m_uiM, iKC, iKF);
        oTpPeel.save(vResE);
        printf("Save size: %d\n", vResE.size());
        oGraph.writeToFile(argv[6]/*"out.txt"*/, vResE);
    }
}


