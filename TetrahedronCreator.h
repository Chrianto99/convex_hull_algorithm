//
// Created by Christos on 9/22/2025.
//

#pragma once
#ifndef CONVEX_HULL_3D_TETRAHEDRON_H
#define CONVEX_HULL_3D_TETRAHEDRON_H

#endif //CONVEX_HULL_3D_TETRAHEDRON_H

#include "iostream"
#include "array"
#include "vector"
#include "Face.h"
#include "algorithm"

template<size_t N>
class TetrahedronCreator {

    std::array<Point_3D,N> points;

public:

     Tetrahedron create_tetrahedron() {

        std::array<int, 6> candidate_points = find_candidate_points();

        Edge longest_edge = get_longest_edge(candidate_points);

        int pid_furthest_from_line = get_pid_furthest_from_line(candidate_points, longest_edge);

        Face face_formed = Face(longest_edge[0], longest_edge[1], pid_furthest_from_line);

        int pid_furthest_from_face = get_pid_furthest_from_face(candidate_points, face_formed);

        std::vector<Face> tetrahedron_faces;

        tetrahedron_faces.push_back(face_formed);
        tetrahedron_faces.emplace_back(face_formed[0], face_formed[1], pid_furthest_from_face);
        tetrahedron_faces.emplace_back(face_formed[0], face_formed[2], pid_furthest_from_face);
        tetrahedron_faces.emplace_back(face_formed[1], face_formed[2], pid_furthest_from_face);

        std::array<int, 4> tetrahedron_points = {face_formed[0], face_formed[1], face_formed[2],
                                                 pid_furthest_from_face};

        Point_3D centroid = get_centroid(tetrahedron_points);

        Tetrahedron tetra = Tetrahedron(tetrahedron_faces, tetrahedron_points, centroid);

        return tetra;


    }

    std::vector<int> find_candidate_points() {
        int x_min = std::numeric_limits<int>::max();
        int y_min = std::numeric_limits<int>::max();
        int z_min = std::numeric_limits<int>::max();

        int x_max = std::numeric_limits<int>::lowest();
        int y_max = std::numeric_limits<int>::lowest();
        int z_max = std::numeric_limits<int>::lowest();

        for (const auto& p : points) {
            if (p.x < x_min) x_min = p.x;
            if (p.y < y_min) y_min = p.y;
            if (p.z < z_min) z_min = p.z;

            if (p.x > x_max) x_max = p.x;
            if (p.y > y_max) y_max = p.y;
            if (p.z > z_max) z_max = p.z;
        }

        return {x_min, y_min, z_min, x_max, y_max, z_max};
    }

    Edge get_longest_edge(std::vector<int> &candidate_points) {

        Edge longest_edge;

        double maxDistance = 0;

        for (auto pid1: candidate_points) {
            for (auto pid2: candidate_points) {

                if (pid1 == pid2) continue;

                double euclidean_distance = Helper::calculate_euclidean_distance(points[pid1], points[pid2]);
                if (euclidean_distance > maxDistance) {
                    maxDistance = euclidean_distance;
                    longest_edge = Edge(pid1, pid2);
                }
            }
        }
        candidate_points.erase(std::remove_if(candidate_points.begin(), candidate_points.end(),
                                              [&](const int pid) { return (pid == points[longest_edge[0]] || pid == points[longest_edge[1]]); }),
                               candidate_points.end());

        return longest_edge;
    }



    Point_3D get_centroid(std::array<int, 4> tetrahedron_vertices) {

        Point_3D centroid{0.0, 0.0, 0.0};

        for (const auto &pid : tetrahedron_vertices) {
            Point_3D p = points[pid];
            centroid.x += p.x;
            centroid.y += p.y;
            centroid.z += p.z;
        }

        centroid.x /= 4.0;
        centroid.y /= 4.0;
        centroid.z /= 4.0;

        return centroid;

    }

    int get_pid_furthest_from_line(std::vector<int> &candidate_points, const Edge line) {

        Point_3D A = points[line[0]];
        Point_3D B = points[line[1]];

        double maxDistance = -1.0;
        int furthest_point_index;

        // Direction vector of the line
        Point_3D AB{B.x - A.x, B.y - A.y, B.z - A.z};
        int index = 0;
        for (const auto pid: candidate_points) {

            // Vector AP
            Point_3D P = points[pid];
            Point_3D AP{P.x - A.x, P.y - A.y, P.z - A.z};

            // Cross product AB × AP
            Point_3D cross{
                    AB.y * AP.z - AB.z * AP.y,
                    AB.z * AP.x - AB.x * AP.z,
                    AB.x * AP.y - AB.y * AP.x
            };

            double crossNorm = std::sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
            double abNorm = std::sqrt(AB.x * AB.x + AB.y * AB.y + AB.z * AB.z);
            double distance = crossNorm / abNorm;

            if (distance > maxDistance) {
                maxDistance = distance;
                furthest_point_index = pid;
            }
        }

        candidate_points.erase(std::remove_if(candidate_points.begin(), candidate_points.end(),
                                              [&](const int p) { return (p == furthest_point_index); }),
                               candidate_points.end());

        return furthest_point_index;
    }

    int get_pid_furthest_from_face(std::vector<int> &candidate_points, const Face &face) {
        // Step 1: compute face normal using cross product

        Point_3D A = points[face[0]];
        // Step 2: find point with max distance to plane
        double maxDistance = -1.0;
        int furthest_point_index;

        for (const auto &pid: candidate_points) {
            // Vector AP
            Point_3D P = points[pid];
            Point_3D AP{P.x - A.x, P.y - A.y, P.z - A.z};

            // Signed distance = dot(AP, normal)
            double distance = std::abs(AP.x * face.normal_vec.x + AP.y * face.normal_vec.y + AP.z * face.normal_vec.z);

            if (distance > maxDistance) {
                maxDistance = distance;
                furthest_point_index = pid;
            }
        }

        return furthest_point_index;
    }

    Point_3D calculate_normal(Face &face, const Point_3D &centroid) {
        const Point_3D &p0 = points[face[0]];
        const Point_3D &p1 = points[face[1]];
        const Point_3D &p2 = points[face[2]];

        // Edge vectors
        double ux = p1.x - p0.x;
        double uy = p1.y - p0.y;
        double uz = p1.z - p0.z;

        double vx = p2.x - p0.x;
        double vy = p2.y - p0.y;
        double vz = p2.z - p0.z;

        // Cross product u × v
        Point_3D normal{
                uy * vz - uz * vy,
                uz * vx - ux * vz,
                ux * vy - uy * vx
        };

        Point_3D centroidVec{centroid.x - p0.x, centroid.y - p0.y, centroid.z - p0.z};
        double dotProd = normal.x * centroidVec.x + normal.y * centroidVec.y + normal.z * centroidVec.z;

        double length = std::sqrt(normal.x * normal.x +
                                  normal.y * normal.y +
                                  normal.z * normal.z);
        if (length > 1e-12) { // avoid division by zero
            normal.x /= length;
            normal.y /= length;
            normal.z /= length;
        }

        if (dotProd > 0) face.set_normal_vec(normal.x * (-1), normal.y * (-1), normal.z * (-1));

        else face.set_normal_vec(normal);

        return normal;
    }



};
