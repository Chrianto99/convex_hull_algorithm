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
#include "ConflictGraph.h"

template<size_t N>
class convex_hull_3D {

    TetrahedronCreator<N> tetrahedronCreator;

public:

    std::vector<Face> create_convex_hull_3D() {

        std::array<double, 3> limits = {10, 10, 10};

        std::array<Point_3D, N> points = generate_points(limits, N);

        Tetrahedron tetrahedron = tetrahedronCreator.create_tetrahedron(points);

        std::vector<Face> tetrahedron_faces = tetrahedron.faces;

        ConflictGraph conflictGraph;


        for (auto point: points) conflictGraph.insert_point(point);

        for (auto &f: tetrahedron_faces) conflictGraph.insert_face(f);

        update_conflict_graph(conflictGraph, tetrahedron.centroid);

        return tetrahedron_faces;

    }

    void update_conflict_graph(ConflictGraph &conflictGraph, Point_3D &centroid) {

        if (conflictGraph.get_points_map().empty()) return;

        // Calculate normals
        for (auto &[fid ,face] : conflictGraph.get_face_map()) {
            tetrahedronCreator.calculate_normal(face, centroid);

            for (auto &[pid ,point] : conflictGraph.get_points_map()) {

                Point_3D A =  conflictGraph.get_points_map()[face[0]];
                Point_3D n = face.normal_vec;

                Point_3D AP{point.x - A.x, point.y - A.y, point.z - A.z};
                double dist = n.x * AP.x + n.y * AP.y + n.z * AP.z;
                if (dist > 0) {
                    face.outside_pids.insert(pid);
                    point.visible_fids.insert(fid);
                }
            }

        }

        std::vector<Face> new_faces = find_new_faces(conflictGraph);
        update_conflict_graph(conflictGraph, centroid);
    }

    std::vector<Face> find_new_faces(ConflictGraph &conflictGraph) {

        std::unordered_map<int, Face> face_map = conflictGraph.get_face_map();
        std::unordered_map<int, Point_3D> point_map = conflictGraph.get_points_map();

        //Iterate all faces in conflict graph
        for (auto &[fid, face] : face_map) {

            // Find the furthest outside point
            int furthest_point_id = conflictGraph.get_point_furthest_from_face(face);

            Point_3D& furthest_point = conflictGraph.get_points_map().at(furthest_point_id);

            // Find the horizon edges
            std::unordered_set<Edge> horizon_edges;

            for (auto fid : furthest_point.visible_fids){

                Face& visible_face = face_map.at(fid);

                for (auto& edge : visible_face.edges) {
                    if (horizon_edges.contains(edge)) horizon_edges.erase(edge);
                    else horizon_edges.insert(edge);
                }
            }

            // Connect new faces
            for (auto edge : horizon_edges){
                face_map.emplace((edge[0], edge[1], furthest_point_id));
            }
            // Erase previous face
            face_map.erase(fid);
            point_map.erase(furthest_point_id);

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



};
