#ifndef STREAM_H_
#define STREAM_H_

#include "common.h"

// class Peel is to peel D-truss on a graph
class Stream final {
 private:
    // order
    static bool cmp(const ST_EDGE_ENTRY &e1, const ST_EDGE_ENTRY &e2){
        if (e1.t != e2.t)
        {
            return e1.t < e2.t;
        }
        else if (e1.x != e2.x)
        {
            return e1.x < e2.x;
        }
        else
        {
            return e1.y < e2.y;
        }
    }
    // the set of edges <timestamps, <edges> >
    //map<int, vector<pair<uint32_t, uint32_t> > > m_mpE;
    vector<ST_EDGE_ENTRY> m_vE;

 public:
    uint32_t m_uiN;  // the # of vertices
    uint32_t m_uiM;  // the # of edges
    uint32_t m_uiTS;  // start timestamps
    uint32_t m_uiTE;  // end timestamps

    Stream();

    void getEdges(uint32_t uiTS, uint32_t uiTE, vector<ST_EDGE_ENTRY> &vDesEdges);
    void readFromFile(char* pcFile);
    void writeToFile(char* pcFile, vector<pair<uint32_t, uint32_t> > &vEdges) const;
};

#endif
