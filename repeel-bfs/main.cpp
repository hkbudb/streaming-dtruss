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

extern std::chrono::duration<double> g_tUptIns;
extern std::chrono::duration<double> g_tUptRm;


int main(int argc, char** argv) {
    // read the graph and peel D-truss
    //ASSERT_MSG(5 < argc, "wrong parameters");
    int iWd = 6000;
    int iStep = 600;
    int iKC = 1;
    int iKF = 1;
    char *filepath = "/Users/xkliao/CLionProjects/StreamingNew/data/sx.txt";
    sscanf(argv[2], "%d", &iWd);
    sscanf(argv[3], "%d", &iStep);
    sscanf(argv[4], "%d", &iKC);
    sscanf(argv[5], "%d", &iKF);

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
    Update oUpdateG(vWEdges, oGraph.m_uiN, oGraph.m_uiM, iKC, iKF);
    while (iEndT < oGraph.m_uiTE)
    {
        int ProcssedEdgeNum = 0;
        const auto preTime1 = g_tUptRm, preTime2 = g_tUptIns;  //for throughput calculation
        iBatch = COMMON_MIN(iBatch, oGraph.m_uiTE - iEndT);
        int iRmEnd = iStartT + iBatch - 1;
        vWEdges.clear();
        oGraph.getEdges(iStartT, iRmEnd, vWEdges);
        oUpdateG.remove(vWEdges);
        iStartT += iBatch;
        ProcssedEdgeNum += vWEdges.size();

        int iInsEnd = iEndT + iBatch;
        vWEdges.clear();
        oGraph.getEdges(iEndT + 1, iInsEnd, vWEdges);
        //printf("get edges: %d from [%d, %d]\n", vWEdges.size(), iEndT + 1, iInsEnd);
        oUpdateG.add(vWEdges);
        iEndT = iInsEnd;
        ProcssedEdgeNum+= vWEdges.size();
        const auto endTime1 = g_tUptRm, endTime2 = g_tUptIns;  //for throughput calculation
        if( CurrentTime % (((int) oGraph.m_uiTE / iBatch) / 20 + 1) == 0) 
        {
              //cout << CurrentTime*iBatch << " " << (long) 1000*ProcssedEdgeNum / (std::chrono::duration<double, std::milli>(endTime1 - preTime1).count() + std::chrono::duration<double, std::milli>(endTime2 - preTime2).count()) << endl;
        }
        CurrentTime++;

        //break;
    }
    const auto end = std::chrono::steady_clock::now();
    const auto dif = end - beg;


    //vertification
    vector<pair<uint32_t, uint32_t>> vResE;
    oUpdateG.save(vResE);
    std::sort(vResE.begin(), vResE.end());
    oGraph.writeToFile("../out-bfs.txt", vResE);
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

    printf("g_tUptIns: \x1b[1;31m%f\x1b[0m ms.\n",
           std::chrono::duration<double, std::milli>(g_tUptIns).count());

    printf("g_tUptRm: \x1b[1;31m%f\x1b[0m ms.\n",
           std::chrono::duration<double, std::milli>(g_tUptRm).count());

    printf("g_tUptIns+g_tUptRm: %.3f ms\n",
           std::chrono::duration<double, std::milli>(g_tUptIns+g_tUptRm).count());


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

