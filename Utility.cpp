#pragma once
#include "Utility.h"

Utility::Utility() {


}

Utility::~Utility() {


}

bool Utility::CompareNodeData(const NodeData& a, const NodeData& b) {

    if (a.global_x != b.global_x)
        return a.global_x < b.global_x;
    if (a.global_y != b.global_y)
        return a.global_y < b.global_y;
    return a.global_z < b.global_z;
}


std::vector<std::tuple<int, int, int>> Utility::removeDuplicates(const std::vector<std::tuple<int, int, int>>& inputVector) {
    std::vector<std::tuple<int, int, int>> resultVector;
    std::set<std::tuple<int, int, int>> uniqueSet(inputVector.begin(), inputVector.end());

    std::copy(uniqueSet.begin(), uniqueSet.end(), std::back_inserter(resultVector));

    return resultVector;
}


std::map<std::tuple<int, int, int>, int> Utility::findRepeatedTuples(std::vector<std::tuple<int, int, int>>& tuples) {
    std::map<std::tuple<int, int, int>, int> repeatedTuples;

    for (const auto& t : tuples) {
        if (repeatedTuples.count(t)) {
            repeatedTuples[t]++;
        }
        else {
            repeatedTuples[t] = 1;
        }
    }

    for (const auto& rt : repeatedTuples) {
        std::cout << "Tuple (" << std::get<0>(rt.first) << ", " << std::get<1>(rt.first) << ", " << std::get<2>(rt.first) << ") repeated " << rt.second << " times" << std::endl;
    }

    return repeatedTuples;
}

//bool Geometry::isPointInsideTetrahedra(const C3t3& c3t3, Point& point) {
//    Point_3 queryPoint(point.x, point.y, point.z);
//
//    typedef C3t3::Triangulation Triangulation_3;
//    typedef typename C3t3::Triangulation::Finite_edges_iterator Edge_iterator;
//    typedef typename C3t3::Triangulation::Geom_traits::Tetrahedron_3 Tetrahedron_3;
//    typedef typename C3t3::Triangulation::Geom_traits::Oriented_side Oriented_side;
//
//    for (C3t3::Cell_iterator it = c3t3.cells_begin(); it != c3t3.cells_end(); ++it) {
//
//        //for (int i = 0; i < 4; ++i) {
//        //    if (it->vertex(i)->point() == queryPoint)
//        //        return true;
//        //}
//
//        //for (Edge_iterator edge_it = c3t3.triangulation().finite_edges_begin(); edge_it != c3t3.triangulation().finite_edges_end(); ++edge_it) {
//        //    const typename C3t3::Triangulation::Edge& edge = *edge_it;
//        //    const typename C3t3::Triangulation::Vertex_handle& v1 = edge.first->vertex(edge.second);
//        //    const typename C3t3::Triangulation::Vertex_handle& v2 = edge.first->vertex(edge.third);
//
//        //    if (c3t3.triangulation().collinear(v1->point(), v2->point(), CGAL::Weighted_point_3<K>(queryPoint)))
//        //        return true;
//
//        //}
//
//        for (C3t3::Facet_iterator facet_it = c3t3.facets_begin(); facet_it != c3t3.facets_end(); ++facet_it) {
//            const C3t3::Cell_handle& cell = facet_it->first;
//            const int index = facet_it->second;
//
//            if (c3t3.triangulation().triangle(cell, index).has_on(queryPoint))
//                return true;
//        }
//
//        Tetrahedron_3 tetrahedron = c3t3.triangulation().tetrahedron(it);
//        bool inside = tetrahedron.has_on_bounded_side(queryPoint);
//        if (inside)
//            return true;
//    }
//
//    return false;
//}
