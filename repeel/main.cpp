
/**
*  This is the most basic algorithm baseline, whenver there comes a single/batch of edges, repeel the whole thing
*/

//void process_mem_usage(double& vm_usage, double& resident_set)
//{
//    using std::ios_base;
//    using std::ifstream;
//    using std::string;
//
//    vm_usage     = 0.0;
//    resident_set = 0.0;
//
//    // 'file' stat seems to give the most reliable results
//    //
//    ifstream stat_stream("/proc/self/stat",ios_base::in);
//
//    // dummy vars for leading entries in stat that we don't care about
//    //
//    string pid, comm, state, ppid, pgrp, session, tty_nr;
//    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
//    string utime, stime, cutime, cstime, priority, nice;
//    string O, itrealvalue, starttime;
//
//    // the two fields we want
//    //
//    unsigned long vsize;
//    long rss;
//
//    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
//                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
//                >> utime >> stime >> cutime >> cstime >> priority >> nice
//                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
//
//    stat_stream.close();
//
//    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
//    vm_usage     = vsize / 1024.0;
//    resident_set = rss * page_size_kb;
//}




#include "common.h"
#include "Stream.h"
#include "Peel.h"

int main(int argc, char** argv) {
    // read the graph and peel D-truss
    //ASSERT_MSG(4 < argc, "wrong parameters");
    bool h_index_optimization_on = false;
    int iWd = 11000;
    int iStep = 1100;
    int iKC = 1;
    int iKF = 0;
    char *filepath = "/Users/xkliao/CLionProjects/StreamingNew/data/ubuntu.txt";
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
    int iEndT = oGraph.m_uiTS + iWd - 1;
    int iStartT = oGraph.m_uiTS;
    iStartT = iStartT > 0? iStartT: 0;
    int iBatch = iStep;
    oGraph.getEdges(iStartT, iEndT, vWEdges);

    Peel oPeel(vWEdges, h_index_optimization_on, iKC, iKF);
    oPeel.start(iKC, iKF);
    oGraph.writeToFile("../initial-community.txt", oPeel.m_vResE);
    
    int CurrentTime = 0;
    vector<long> Throughputs;
    while (iEndT <= oGraph.m_uiTE)
    {
        iBatch = COMMON_MIN(iBatch, oGraph.m_uiTE - iEndT);
        iStartT += iBatch;
        iEndT += iBatch;
        //int iStartT = iEndT - iWd + 1;
        //iStartT = iStartT > 0? iStartT: 0;
        vWEdges.clear();
        const auto preTime1 = std::chrono::steady_clock::now();
        oGraph.getEdges(iStartT, iEndT, vWEdges);
        Peel mPeel(vWEdges, h_index_optimization_on, iKC, iKF);
        mPeel.start(iKC, iKF);
        const auto endTime1 = std::chrono::steady_clock::now();


        //for case study
//        string path = "out-" + to_string(iEndT) + ".txt";
//        char pathArr[path.length()+1];
//        strcpy(pathArr, path.c_str());
//        oGraph.writeToFile(pathArr, oPeel.m_vResE);
//        cout << path << endl;

//        if( CurrentTime % (((int) oGraph.m_uiTE / iBatch) / 20 + 1) == 0) {
//            cout << CurrentTime*iBatch << " " << (long) 1000*vWEdges.size() / (std::chrono::duration<double, std::milli>(endTime1 - preTime1).count()) << endl;
//            Throughputs.push_back( (long) 1000*vWEdges.size() / (std::chrono::duration<double, std::milli>(endTime1 - preTime1).count()));
//        }
        ++iEndT;
        ++CurrentTime;
    }
    const auto end = std::chrono::steady_clock::now();
    const auto dif = end - beg;

    //print community quality
    //pair<double,double> DtQuality = oPeel.getDtQuality(oPeel.m_vResE);
    //printf("Save size: %d CMSin: %f CMSout: %f\n", oPeel.m_vResE.size() , DtQuality.first, DtQuality.second);

    printf("Query costs \x1b[1;31m%f\x1b[0m ms.\n",
           std::chrono::duration<double, std::milli>(dif).count()*2);



    //print the thoughput
//    long avgToughput = 0, ThroughputsCnt = 0;
//    for(long thp : Throughputs)
//    {
//        ThroughputsCnt++;
//        avgToughput += thp;
//    }
//    avgToughput = (long) avgToughput / ThroughputsCnt;
//    cout << "throughput: " <<  avgToughput / 3  << endl;

    //printf("Get window edges: %d\n", vWEdges.size());
    //printf("Get D-truss edges: %d\n", oPeel.m_vResE.size());
    // if (6 == argc)
    // {
    //     printf("Save %s\n", argv[6]);
    //     printf("save size:%d\n",oPeel.m_vResE.size());
    //     oGraph.writeToFile(argv[6], oPeel.m_vResE);
    // }
}
