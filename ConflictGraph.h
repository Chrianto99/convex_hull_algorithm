//
// Created by Christos on 10/11/2025.
//

#ifndef CONVEX_HULL_3D_CONFLICTGRAPH_H
#define CONVEX_HULL_3D_CONFLICTGRAPH_H

#endif //CONVEX_HULL_3D_CONFLICTGRAPH_H

#include "unordered_map"
#include "unordered_set"
#include "Face.h"

struct ConflictGraph{
    int pid = 0, fid = 0;
    std::unordered_map<int, Point_3D> points_map;
    std::unordered_map<int, Face> faces_map;

    ConflictGraph() = default;

    std::unordered_map<int, Point_3D> &get_points_map(){
        return points_map;
    }

    std::unordered_map<int, Face> &get_face_map(){
        return faces_map;
    }

    void insert_point(Point_3D &point) {
        points_map.insert({pid++, point});

    }

    void insert_face(Face &face) {
        faces_map.insert({fid++, face});

    }

    int get_point_furthest_from_face(Face& face){
        int furthest_point_id;
        Point_3D A = points_map[face[0]];
        double max_distance = -1;

        for (auto &pid: face.outside_pids){
                // Vector AP
                Point_3D P = points_map[pid];
                Point_3D AP{P.x - A.x, P.y - A.y, P.z - A.z};

                // Signed distance = dot(AP, normal)
                double distance = std::abs(AP.x * face.normal_vec.x + AP.y * face.normal_vec.y + AP.z * face.normal_vec.z);

                if (distance > max_distance) {
                    max_distance = distance;
                    furthest_point_id = pid;
                }

        }

        return furthest_point_id;

    }

};