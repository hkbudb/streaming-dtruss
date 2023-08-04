#ifndef DECOMPOSITION_H_
#define DECOMPOSITION_H_

#include "common.h"

// Do intiali truss decomposition to get the 
class Decomposition final {
 private:
    typedef struct {
        int u,v;
    } TEdge;

    bool operator<(TEdge e1, TEdge e2) {
        return e1.u<e2.u || (e1.u==e2.u && e1.v<e2.v);
    }

    typedef map<int,int> MII;
    typedef vector<int> VI;
    typedef MII::iterator EdgeIter;

    const int maxClass = 1000;

    ifstream fin;
    ofstream fout;
    string infile, outfile;

    int n, m;
    VI mapto;
    VI deg, bin;
    vector<TEdge> binEdge;
    vector<VI> A;
    vector<MII> adj, pos;

    int cntClass[maxClass];



    inline bool compVertex(int i, int j);
    inline void printClass(int u, int v, int cls);
    inline void updateSupport(int u, int v, int delta);
    inline void removeEdge(int u, int v);
    


 public:
    uint32_t m_uiN;  // the # of vertices
    uint32_t m_uiM;  // the # of edges
    uint32_t m_uiTS;  // start timestamps
    uint32_t m_uiTE;  // end timestamps


    

     //void getEdges(uint32_t uiTS, uint32_t uiTE, vector<pair<uint32_t, uint32_t> > &vDesEdges);
     //void readFromFile(char* pcFile);
     //void writeToFile(char* pcFile, vector<pair<uint32_t, uint32_t> > &vEdges) const;
};

#endif
