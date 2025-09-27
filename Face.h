//
// Created by Christos on 9/22/2025.
//
#pragma once
#ifndef CONVEX_HULL_3D_FACE_H
#define CONVEX_HULL_3D_FACE_H

#endif //CONVEX_HULL_3D_FACE_H

#include "iostream"
#include "array"
#include "cmath"

struct Point_3D {
    double x, y, z;

    Point_3D() = default;

    Point_3D(double x, double y, double z) : x(x), y(y), z(z) {}

    bool operator<(const Point_3D &other) const {
        for (int i = 0; i < 3; ++i) {
            if (x != other.x) return x < other.x;
            if (y != other.y) return y < other.y;
            if (z != other.z) return z < other.z;
        }
        return false;
    }

    bool operator==(Point_3D p2) const  {

        if (x == p2.x && y == p2.y && z == p2.z) return true;

        return false;

    }
};

namespace std {
    template<>
    struct hash<Point_3D> {
        size_t operator()(const Point_3D &p) const {
            size_t hx = hash<double>{}(p.x);
            size_t hy = hash<double>{}(p.y);
            size_t hz = hash<double>{}(p.z);

            // Combine hashes
            size_t seed = hx;
            seed ^= hy + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= hz + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };
}




struct Edge {

public:
    int vertex_indices[2];


    Edge() = default;

    Edge(const int &p1, const int &p2)
            : vertex_indices{p1, p2} {}



    bool operator<(const Edge &other) const {
        for (int i = 0; i < 2; ++i) {
            if (vertex_indices[i] != other.vertex_indices[i]) return vertex_indices[i] < other.vertex_indices[i];
            if (vertex_indices[i] != other.vertex_indices[i]) return vertex_indices[i] < other.vertex_indices[i];
            if (vertex_indices[i] != other.vertex_indices[i]) return vertex_indices[i] < other.vertex_indices[i];
        }
        return false;
    }

    bool operator==(const Edge &other) const {
        if (vertex_indices[0] == other.vertex_indices[0] && vertex_indices[1] == other.vertex_indices[1]) return true;
        return false;
    }

    int operator[](int i) const {
        return vertex_indices[i];
    }



};

namespace std {
    template<>
    struct hash<Edge> {
        size_t operator()(const Edge &f) const noexcept {
            // Copy vertex indices
            std::array<int,3> verts = {f.vertex_indices[0], f.vertex_indices[1]};

            // Sort the indices so that order doesn’t matter
            std::sort(verts.begin(), verts.end());

            // Combine hashes of the three indices
            size_t seed = 0;
            for (int v : verts) {
                seed ^= std::hash<int>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}

struct Face {

public:
    int vertex_indices[3];
    Point_3D normal_vec;


    Face(const int &p1, const int &p2, const int &p3)
            : vertex_indices{p1, p2, p3} {}

    void set_normal_vec(double x, double y, double z) {
        normal_vec = Point_3D(x, y, z);

    }

    void set_normal_vec(Point_3D normal) {
        normal_vec = normal;
    }


    bool operator<(const Face &other) const {
        for (int i = 0; i < 3; ++i) {
            if (vertex_indices[i] != other.vertex_indices[i]) return vertex_indices[i] < other.vertex_indices[i];
            if (vertex_indices[i] != other.vertex_indices[i]) return vertex_indices[i] < other.vertex_indices[i];
            if (vertex_indices[i] != other.vertex_indices[i]) return vertex_indices[i] < other.vertex_indices[i];
        }
        return false;
    }

    bool operator==(const Face &other) const {
        if (vertex_indices[0] == other.vertex_indices[0] && vertex_indices[1] == other.vertex_indices[1] && vertex_indices[2] == other.vertex_indices[2]) return true;
        return false;
    }

    int operator[](int i) const {
        return vertex_indices[i];
    }



};

namespace std {
    template<>
    struct hash<Face> {
        size_t operator()(const Face &f) const noexcept {
            // Copy vertex indices
            std::array<int,3> verts = {f.vertex_indices[0], f.vertex_indices[1], f.vertex_indices[2]};

            // Sort the indices so that order doesn’t matter
            std::sort(verts.begin(), verts.end());

            // Combine hashes of the three indices
            size_t seed = 0;
            for (int v : verts) {
                seed ^= std::hash<int>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}

struct Tetrahedron {
public:
    std::vector<Face> faces;
    std::array<int, 4> vertices;
    Point_3D centroid;

//    Tetrahedron() = default;

    Tetrahedron(const std::vector<Face> &f, const std::array<int, 4> &v, const Point_3D &c)
            : faces(f), vertices(v), centroid(c) {}
};

class Helper {
public:

    static double calculate_euclidean_distance(const Point_3D &p1, const Point_3D &p2) {
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        double dz = p1.z - p2.z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }






};






