//
// Created by Christos on 9/25/2025.
//

#ifndef CONVEX_HULL_3D_BIPARTITEGRAPH_H
#define CONVEX_HULL_3D_BIPARTITEGRAPH_H

#endif //CONVEX_HULL_3D_BIPARTITEGRAPH_H
#include "iostream"
#include "vector"
#include "unordered_map"
#include "Face.h"

struct NodeLeft;
struct NodeRight;

struct NodeLeft{
    Face face;
    std::vector<NodeRight*> connections;

    explicit NodeLeft(Face& f) : face(f){};

    bool has_no_edges() const {
        return connections.empty();
    }
};

struct NodeRight{
    int point_id;
    std::vector<NodeLeft*> connections;

    explicit NodeRight(int p) : point_id(p){};

    bool has_no_edges() const {
        return connections.empty();
    }


};

class BipartiteGraph{

    std::vector<NodeLeft> left_nodes;
    std::vector<NodeRight> right_nodes;

public:

    std::vector<NodeLeft>& get_left_nodes(){
        return left_nodes;
    }

    std::vector<NodeRight>& get_right_nodes(){
        return right_nodes;
    }

    void insert_left(Face& face){
        left_nodes.emplace_back((face));
    }

    void insert_right(int point_id){
        right_nodes.emplace_back((point_id));
    }


    void connect_nodes(NodeLeft* node_left, NodeRight* node_right){
        node_left->connections.push_back(node_right);
        node_right->connections.push_back(node_left);
    }

    bool has_no_edges(){
        for (auto& node : left_nodes){
            if (!node.connections.empty()) return false;
        }
        return true;
    }








};
