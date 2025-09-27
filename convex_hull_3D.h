//
// Created by Christos on 9/22/2025.
//
#pragma once
#ifndef CONVEX_HULL_3D_CONVEX_HULL_3D_H
#define CONVEX_HULL_3D_CONVEX_HULL_3D_H

#endif //CONVEX_HULL_3D_CONVEX_HULL_3D_H

#include <random>
#include <algorithm>
#include <set>
#include <unordered_set>
#include "Face.h"
#include "iostream"
#include "vector"
#include "map"
#include "array"
#include "TetrahedronCreator.h"
#include "BipartiteGraph.h"
#include "CustomUnorderedSet.h"

template<size_t N>
class convex_hull_3D {

    TetrahedronCreator<N> tetrahedronCreator;


public:

    std::vector<Face> create_convex_hull_3D() {

        std::array<double, 3> limits = {10, 10, 10};

        std::array<Point_3D, N> points = generate_points(limits, N);

        Tetrahedron tetrahedron = tetrahedronCreator.create_tetrahedron(points);

        std::vector<Face> tetrahedron_faces = tetrahedron.faces;

        BipartiteGraph bGraph;

        for (int i = 0; i < points.size(); ++i) bGraph.insert_right(i);

        for (auto &f: tetrahedron_faces) bGraph.insert_left(f);

        update_conflict_graph(bGraph, tetrahedron.centroid);

        return tetrahedron_faces;

    }

    void update_conflict_graph(BipartiteGraph &bGraph, Point_3D &centroid, std::array<Point_3D, N> &points) {

        if (bGraph.has_no_edges()) return;

        for (auto &left_node: bGraph.get_left_nodes()) {
            Face &face = left_node.face;
            tetrahedronCreator.calculate_normal(face, centroid);

            for (auto &right_node: bGraph.get_right_nodes()) {
                Point_3D &candidate_point = points[right_node.point_id];
                Point_3D A = points[face[0]];
                Point_3D n = face.normal_vec;

                Point_3D AP{candidate_point.x - A.x, candidate_point.y - A.y, candidate_point.z - A.z};
                double dist = n.x * AP.x + n.y * AP.y + n.z * AP.z;
                if (dist > 0) bGraph.connect_nodes(&left_node, &right_node);
            }

        }

        if (bGraph.has_no_edges()) return;

        std::vector<Face> new_faces = find_new_faces(bGraph);
        update_conflict_graph(bGraph, centroid, points);
    }

    std::vector<Face> find_new_faces(BipartiteGraph &bGraph) {

        CustomUnorderedSet<Edge> edge_set;

        for (auto &right_node: bGraph.get_right_nodes()) {

            if (right_node.has_no_edges()) continue;

        }

    }


    std::array<Point_3D, N> generate_points(const std::array<double, 3> &lims, int num_points) {
        std::array<Point_3D, N> points;
        points.reserve(num_points);

        std::random_device rd;
        std::mt19937 gen(rd());

        std::uniform_real_distribution<double> dist_x(0.0, lims[0]);
        std::uniform_real_distribution<double> dist_y(0.0, lims[1]);
        std::uniform_real_distribution<double> dist_z(0.0, lims[2]);

        for (int i = 0; i < num_points; ++i) {
            points[i] = {dist_x(gen), dist_y(gen), dist_z(gen)};
        }

        return points;
    }




//    static std::vector<Point_3D> find_outside_points(const Face &face, const std::set<Point_3D> &points) {
//
//        std::vector<Point_3D> outside_points;
//
//        // Compute normal
//        Point_3D A = face.vertices[0];
//        Point_3D n = face.normal_vec;
//
//        // Find outside points
//        for (const auto &P: points) {
//            if (P == face.vertices[0] || P == face.vertices[1] || P == face.vertices[2]) continue;
//            Point_3D AP{P.x - A.x, P.y - A.y, P.z - A.z};
//            double dist = n.x * AP.x + n.y * AP.y + n.z * AP.z;
//            if (dist > 0) outside_points.push_back(P);
//        }
//
//        return outside_points;
//    }



//    static std::vector<Face> find_new_faces(BipartiteGraph& bGraph) {
//
//
//        for (auto &left_node: bGraph.get_left_nodes()) {
//
//            if (left_node.has_no_edges()) continue;
//
//            Face& face = left_node.face;
//            double max_distance = 0;
//            NodeRight* furthest_point_node;
//            for (auto right_node : left_node.edges){
//                Point_3D& out_point = right_node->out_point;
//
//                double distance_from_face = get_distance_from_face(out_point, face);
//                if (distance_from_face > max_distance) {
//                    max_distance = distance_from_face;
//                    furthest_point_node = right_node;
//                }
//            }
//
//        }
//
//
//
//
//    }


};
