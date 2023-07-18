#pragma once

#ifndef _USER_DEFINED_STRUCTS_H_
#define _USER_DEFINED_STRUCTS_H_
#include <iostream>
#include <unordered_map>

enum class NodeType {

    Halo,
    Solid,
    Fluid
};


struct Point {

    double x;
    double y;
    double z;
    NodeType Type;
};

struct Line {
    Point start;
    Point end;
};

struct Triangle {
    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;
    double nx, ny, nz;
};


struct Tetrahedron {
    Point vertices[4];
};

struct BoundingBox {
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double min_z;
    double max_z;
};

struct NodeData {
    int global_x;
    int global_y;
    int global_z;
    double velX;
    double velY;
    double velZ;
    double velMag;
};

struct PairHash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

#endif // !_USER_DEFINED_STRUCTS_H_