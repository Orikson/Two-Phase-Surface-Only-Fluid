#include "fluid.h"

double TPFluid::partitionVolume(const vector<size_t>& faces) const {
    double volume = 0;
    vec3 xref(0, 0, 0);
    for (size_t i = 0; i < faces.size(); i ++) {
        LosTopos::Vec3st t = mesh().m_tris[faces[i]];
        vec3 x0 = pos(t[0]);
        vec3 x1 = pos(t[1]);
        vec3 x2 = pos(t[2]);
        double v = (x0 - xref).cross(x1 - xref).dot(x2 - xref) / 6;

        LosTopos::Vec2i l = mesh().get_triangle_label(faces[i]);
        volume += v * (l[0] > l[1] ? 1 : -1);
    }

    return volume;
}

double TPFluid::vertIntegralMeanCurvature(size_t v) const {
    double H = 0;
    for (size_t i = 0; i < mesh().m_vertex_to_edge_map[v].size(); i ++) {
        size_t e = mesh().m_vertex_to_edge_map[v][i];
        assert(mesh().m_edge_to_triangle_map[e].size() == 2);
        size_t f0 = mesh().m_edge_to_triangle_map[e][0];
        size_t f1 = mesh().m_edge_to_triangle_map[e][1];
        bool o0 = (mesh().oriented(mesh().m_edges[e][0], mesh().m_edges[e][1], mesh().m_tris[f0]) == mesh().get_triangle_label(f0)[0] > mesh().get_triangle_label(f0)[1]);
        bool o1 = (mesh().oriented(mesh().m_edges[e][0], mesh().m_edges[e][1], mesh().m_tris[f1]) == mesh().get_triangle_label(f1)[0] > mesh().get_triangle_label(f1)[1]);
        assert(o0 != o1);
        if (o1)
            std::swap(f0, f1),
            std::swap(o0, o1);
        assert(o0);
        assert(!o1);
        vec3 n0 = faceOutwardNormal(f0);
        vec3 n1 = faceOutwardNormal(f1);

        double angle = angleAroundAxis(n0, n1, edgeTangent(e));

        H += angle * edgeLength(e);
    }

    return H / 2;
}

LosTopos::Vec3st TPFluid::getShuffledTriangle(const LosTopos::Vec3st& t, size_t vertex_to_be_front) const {
    assert(t[0] == vertex_to_be_front || t[1] == vertex_to_be_front || t[2] == vertex_to_be_front);

    LosTopos::Vec3st nt = t;
    while (nt[0] != vertex_to_be_front) {
        size_t tmp = nt[0];
        nt[0] = nt[1];
        nt[1] = nt[2];
        nt[2] = tmp;
    }
    assert(nt[0] == vertex_to_be_front);

    return nt;
}

mat3 TPFluid::getVertexPositions(const LosTopos::Vec3st& t) const {
    mat3 x;
    x.col(0) = pos(t[0]);
    x.col(1) = pos(t[1]);
    x.col(2) = pos(t[2]);
    return x;
}

mat3 TPFluid::getVertexVelocities(const LosTopos::Vec3st& t) const {
    mat3 v;
    v.col(0) = vel(t[0]);
    v.col(1) = vel(t[1]);
    v.col(2) = vel(t[2]);
    return v;
}

double TPFluid::vertInteriorSolidAngle(size_t v) const {
    double sa = 0;
    for (size_t i = 0; i < mesh().m_vertex_to_edge_map[v].size(); i++) {
        size_t ei = mesh().m_vertex_to_edge_map[v][i];
        assert(mesh().m_edge_to_triangle_map[ei].size() == 2);
        LosTopos::Vec2st e = mesh().m_edges[ei];
        vec3 et = edgeTangent(ei);

        size_t f0 = mesh().m_edge_to_triangle_map[ei][0];
        size_t f1 = mesh().m_edge_to_triangle_map[ei][1];
        bool o0 = (mesh().oriented(e[0], e[1], mesh().m_tris[f0]) == (mesh().get_triangle_label(f0)[0] > mesh().get_triangle_label(f0)[1]));
        bool o1 = (mesh().oriented(e[0], e[1], mesh().m_tris[f1]) == (mesh().get_triangle_label(f1)[0] > mesh().get_triangle_label(f1)[1]));
        assert(o0 != o1);
        if (o1)
            std::swap(f0, f1), std::swap(o0, o1);
        assert(o0 && !o1);

        vec3 n0 = faceOutwardNormal(f0);
        vec3 n1 = faceOutwardNormal(f1);
        double dihedral_angle = M_PI + (n0.cross(n1).dot(et) > 0 ? -1 : 1) * acos(std::min(1.0, std::max(-1.0, n0.dot(n1))));

        sa += dihedral_angle;
    }

    sa -= (mesh().m_vertex_to_edge_map[v].size() - 2) * M_PI;

    return sa;
}

