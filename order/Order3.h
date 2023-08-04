//
// Created by xkliao on 2023/6/9.
//

#ifndef ORDER_ORDER3_H
#define ORDER_ORDER3_H
#include "common.h"
#include "queue"


// class order is to update D-truss on a streaming graph with order-based manner
class Order3 final {
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
    typedef struct stOrderEdgeInfo{
        uint32_t eid;   //the eid of current edge, and also the
        uint32_t layer_number;   //which layer this edges belongs to
        //label_t label;      //label in the whole ordered list, in the range [0, m_ui_m^2]
        //int prev;  // the index of the edge before this edge in the ordered list
        //int next;  // the index of the  edge after this edge in the ordered list
    } OrderEdgeEntry;
    // ordered layer entry type
    typedef struct stLayerInfo {
        std::vector<uint32_t> edges;
        int layer_number;
        //int head;  //first edge in the layer
        //uint32_t tail;  //last edge in the layer
        //struct stLayerInfo *prev = nullptr;  // the previous layer
        //struct stLayerInfo *next = nullptr;  //the next layer (preceding this layer)
    } LayerEntry;
    // compare function for the priority queue， min first
    struct CompareOrder {
        bool operator()(const OrderEdgeEntry &a, const OrderEdgeEntry &b) {
            if(a.layer_number < b.layer_number) return true;
            else return false;
            //return a.label < b.label;
        }
    };

    //priority queue, sorted by the order of edges in the ordered list (min-first),
    // for the processing of insertion / deletion
    typedef std::priority_queue<OrderEdgeEntry, std::vector<OrderEdgeEntry>, CompareOrder> PQ;
    PQ pq_;

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

    //ordered list
    vector<LayerEntry> m_vv_ordered_layers;  //the ordered list in the form of layers
    vector<OrderEdgeEntry> m_vv_ordered_list;  //the ordered list of edges
    int m_mu = -1;  //the layer number of the last layer in the ordered layers, -1 means no D-truss retrived

    //order related functions
    void findNeibAndLayer(vector<AdjEntry> &vAdj1, vector<AdjEntry> &vAdj2, vector<pair<uint32_t, uint32_t> > & vTrdP,
                          vector<uint32_t> &layer_support_cnt);
    void removeItem(vector<uint32_t> &vector, uint32_t item);  //delete item from vector
    void compressUpdatedLayers();  //delete empty layers and compress layers to L_{mu}
    void insertNewLayer(uint32_t eid);  //when new layer is created, update list/layers
    void findLayerToInsert(vector<uint32_t> &support_cnt_c, vector<uint32_t> &support_cnt_f, uint32_t eid);
    void insertNeighborsToQueue(vector<bool> &be_in_pq, PQ &m_PQ, vector<vector<pair<uint32_t,uint32_t> >> &neighbors_of_different_layer,
                                uint32_t lower_bound, uint32_t upper_bound);

    bool verify(uint32_t peel_layer_number, uint32_t eid, uint32_t &layer_number_test);

    void findNeib(vector<AdjEntry> &vAdj1, vector<AdjEntry> &vAdj2, vector<pair<uint32_t, uint32_t> > & vTrdP);
    bool findEid(uint32_t x, uint32_t y, uint32_t *piEid);
    uint32_t renamePid(uint32_t x);
    bool addEInfo(uint32_t x, uint32_t y, uint32_t *piEid);


    uint32_t h_index(vector<::uint32_t> &neigh_cycle_support);
    uint32_t insertDT(vector<uint32_t> &vInsE);
    uint32_t rmDT(vector<uint32_t> &vRmE);
    uint32_t rm_bfs_DT(vector<uint32_t> &vRmE);

public:
    Order3(vector<pair<uint32_t, uint32_t> > & vEdges, uint32_t uiMaxN, uint32_t uiMaxM,
           uint32_t uiKc, uint32_t uiKf);
    Order3(const Order3&) = delete;
    Order3& operator=(const Order3&) = delete;

    uint32_t m_uiDTCnt;
    uint32_t m_uiEdgeCnt;

    void peeling(vector<uint32_t> &vEdges, vector<uint32_t> &vFixE, vector<uint32_t> &vResE);
    void add(vector<pair<uint32_t, uint32_t> > & vEdges);
    void remove(vector<pair<uint32_t, uint32_t> > & vEdges);
    void save(vector<pair<uint32_t, uint32_t> > & vEdges);
    pair<double,double> getDtQuality(vector<pair<uint32_t, uint32_t> > & vEdges);

    //priority queue for the processing of insertion /deletion
    PQ PRIORITY_Q(size_t size) {
        std::vector<OrderEdgeEntry> container;
        container.reserve(size);
        CompareOrder compare;
        pq_ = PQ(compare, std::move(container));
    }
    PQ PRIORITY_Q(){}
    inline void push(OrderEdgeEntry d) { pq_.push(d);}
    inline OrderEdgeEntry top() {return pq_.top();}
    inline void pop() {pq_.pop();}
    inline bool empty() {return pq_.empty(); }
    inline void clear() { while (!pq_.empty()){pq_.pop();} }
    inline uint32_t size() {return pq_.size();}

    // the set of edges in D-truss
    //vector<pair<uint32_t, uint32_t> > m_vResE;
};

#endif //ORDER_ORDER3_H
