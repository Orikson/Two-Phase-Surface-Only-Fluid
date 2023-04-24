#include "fluid.h"

// Callbacks from BPS3D.cpp, see license under LosTopos

bool TPFluid::generate_collapsed_position(LosTopos::SurfTrack& st, size_t v0, size_t v1, LosTopos::Vec3d& pos)
{
    if (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1))
    {
        if (vertexIsOnTripleJunction(v0) && vertexIsOnTripleJunction(v1))
        {
            pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
        }
        else if (vertexIsOnTripleJunction(v0))
        {
            pos = st.pm_positions[v0];
        }
        else if (vertexIsOnTripleJunction(v1))
        {
            pos = st.pm_positions[v1];
        }
        else
        {
            pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
        }
    }
    else if (st.vertex_is_any_solid(v0))
    {
        pos = st.pm_positions[v0];
    }
    else if (st.vertex_is_any_solid(v1))
    {
        pos = st.pm_positions[v1];
    }
    else
    {
        pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
    }

    return true;
}

bool TPFluid::generate_split_position(LosTopos::SurfTrack& st, size_t v0, size_t v1, LosTopos::Vec3d& pos)
{
    callback_log() << "solid callback: generate split position: " << v0 << " " << v1 << " " << (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1)) << std::endl;
    pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
    return true;
}

bool TPFluid::generate_snapped_position(LosTopos::SurfTrack& st, size_t v0, size_t v1, LosTopos::Vec3d& pos)
{
    bool v0solid = st.vertex_is_any_solid(v0);
    bool v1solid = st.vertex_is_any_solid(v1);

    // scan the neighborhood right now to prune out the case of a singular solid vertex with no incident solid faces. This can happen as T1 pulls apart a solid vertex which is a singular vertex shared by a free surface and a solid surface (e.g. the middle of the thin sheet in floor splash; see 1450021611).
    bool v0_solid_face = false;
    for (size_t j = 0; j < surf->m_mesh.m_vertex_to_triangle_map[v0].size(); j++)
        if (faceIsSolid(surf->m_mesh.m_vertex_to_triangle_map[v0][j]))
            v0_solid_face = true;
    if (!v0_solid_face)
        v0solid = false;   // a solid vertex incident to no solid faces is not a true solid vertex.

    bool v1_solid_face = false;
    for (size_t j = 0; j < surf->m_mesh.m_vertex_to_triangle_map[v1].size(); j++)
        if (faceIsSolid(surf->m_mesh.m_vertex_to_triangle_map[v1][j]))
            v1_solid_face = true;
    if (!v1_solid_face)
        v1solid = false;   // a solid vertex incident to no solid faces is not a true solid vertex.

    if (v0solid && v1solid)
    {
        if (vertexIsOnTripleJunction(v0) && vertexIsOnTripleJunction(v1))
        {
            pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
        }
        else if (vertexIsOnTripleJunction(v0))
        {
            pos = st.pm_positions[v0];
        }
        else if (vertexIsOnTripleJunction(v1))
        {
            pos = st.pm_positions[v1];
        }
        else
        {
            pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
        }
    }
    else if (v0solid)
    {
        pos = st.pm_positions[v0];
    }
    else if (v1solid)
    {
        pos = st.pm_positions[v1];
    }
    else
    {
        pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
    }

    return true;
}

LosTopos::Vec3c TPFluid::generate_collapsed_solid_label(LosTopos::SurfTrack& st, size_t v0, size_t v1, const LosTopos::Vec3c& label0, const LosTopos::Vec3c& label1)
{
    return LosTopos::Vec3c(label0[0] || label1[0], label0[1] || label1[1], label0[2] || label1[2]);
}

LosTopos::Vec3c TPFluid::generate_split_solid_label(LosTopos::SurfTrack& st, size_t v0, size_t v1, const LosTopos::Vec3c& label0, const LosTopos::Vec3c& label1)
{
    return LosTopos::Vec3c(label0[0] && label1[0], label0[1] && label1[1], label0[2] && label1[2]);
}

LosTopos::Vec3c TPFluid::generate_snapped_solid_label(LosTopos::SurfTrack& st, size_t v0, size_t v1, const LosTopos::Vec3c& label0, const LosTopos::Vec3c& label1)
{
    return LosTopos::Vec3c(label0[0] || label1[0], label0[1] || label1[1], label0[2] || label1[2]);
}

bool TPFluid::generate_edge_popped_positions(LosTopos::SurfTrack& st, size_t oldv, const LosTopos::Vec2i& cut, LosTopos::Vec3d& pos_upper, LosTopos::Vec3d& pos_lower)
{
    return false;
}

bool TPFluid::generate_vertex_popped_positions(LosTopos::SurfTrack& st, size_t oldv, int A, int B, LosTopos::Vec3d& pos_a, LosTopos::Vec3d& pos_b)
{
    double floorz = options.floorz;
    if (pos_a[2] < floorz) pos_a[2] = floorz;
    if (pos_b[2] < floorz) pos_b[2] = floorz;

    return true;
}

bool TPFluid::solid_edge_is_feature(const LosTopos::SurfTrack& st, size_t e)
{
    return edgeIsOnTripleJunction(e);
}

LosTopos::Vec3d TPFluid::sampleVelocity(LosTopos::Vec3d& pos)
{
    return LosTopos::Vec3d(0, 0, 0);
}

bool TPFluid::sampleDirectionalDivergence(const LosTopos::Vec3d& pos, const LosTopos::Vec3d& dir, double& output)
{
    return false;
}

struct CollapseTempData
{
    size_t v0;
    size_t v1;

    vec3 old_x0;
    vec3 old_x1;

    vec3 old_u0;
    vec3 old_u1;
};

void TPFluid::pre_collapse(const LosTopos::SurfTrack& st, size_t e, void** data)
{
    CollapseTempData* td = new CollapseTempData;
    td->v0 = st.m_mesh.m_edges[e][0];
    td->v1 = st.m_mesh.m_edges[e][1];

    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);

    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);

    *data = (void*)td;
    callback_log() << "pre collapse: " << e << ": " << td->v0 << " " << td->v1 << std::endl;
}

void TPFluid::post_collapse(const LosTopos::SurfTrack& st, size_t e, size_t merged_vertex, void* data)
{
    CollapseTempData* td = (CollapseTempData*)data;
    callback_log() << "post collapse: " << e << ": " << td->v0 << " " << td->v1 << " => " << merged_vertex << std::endl;
    assert((st.m_mesh.vertex_is_deleted(td->v0) && merged_vertex == td->v1) || (st.m_mesh.vertex_is_deleted(td->v1) && merged_vertex == td->v0));

    vec3 merged_x = vc(st.pm_positions[merged_vertex]);
    double s = (merged_x - td->old_x0).dot(td->old_x1 - td->old_x0) / (td->old_x1 - td->old_x0).squaredNorm();
    if (s > 1) s = 1;
    if (s < 0) s = 0;
    vec3 new_u = td->old_u0 * (1 - s) + td->old_u1 * s;

    vel(merged_vertex) = new_u;
}

struct SplitTempData
{
    size_t v0;
    size_t v1;

    vec3 old_x0;
    vec3 old_x1;

    vec3 old_u0;
    vec3 old_u1;
};

void TPFluid::pre_split(const LosTopos::SurfTrack& st, size_t e, void** data)
{
    SplitTempData* td = new SplitTempData;
    td->v0 = st.m_mesh.m_edges[e][0];
    td->v1 = st.m_mesh.m_edges[e][1];

    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);

    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);

    *data = (void*)td;
    callback_log() << "pre split: " << e << ": " << td->v0 << " " << td->v1 << std::endl;
}

void TPFluid::post_split(const LosTopos::SurfTrack& st, size_t e, size_t new_vertex, void* data)
{
    SplitTempData* td = (SplitTempData*)data;
    callback_log() << "post split: " << e << ": " << td->v0 << " " << td->v1 << " => " << new_vertex << std::endl;

    vec3 midpoint_x = vc(st.pm_positions[new_vertex]);
    double s = (midpoint_x - td->old_x0).dot(td->old_x1 - td->old_x0) / (td->old_x1 - td->old_x0).squaredNorm();
    if (s > 1) s = 1;
    if (s < 0) s = 0;
    vec3 new_u = td->old_u0 * (1 - s) + td->old_u1 * s;

    vel(new_vertex) = new_u;
}

void TPFluid::pre_flip(const LosTopos::SurfTrack& st, size_t e, void** data)
{

}

void TPFluid::post_flip(const LosTopos::SurfTrack& st, size_t e, void* data)
{

}

struct T1TempData
{

};

void TPFluid::pre_t1(const LosTopos::SurfTrack& st, size_t v, void** data)
{
    callback_log() << "pre t1: " << v << std::endl;
}

void TPFluid::post_t1(const LosTopos::SurfTrack& st, size_t v, size_t a, size_t b, void* data)
{
    callback_log() << "post t1: " << v << " => " << a << " " << b << std::endl;

    vel(a) = vel(v);
    vel(b) = vel(v);
}

struct FaceSplitTempData
{
    size_t v0;
    size_t v1;
    size_t v2;

    vec3 old_x0;
    vec3 old_x1;
    vec3 old_x2;

    vec3 old_u0;
    vec3 old_u1;
    vec3 old_u2;
};

void TPFluid::pre_facesplit(const LosTopos::SurfTrack& st, size_t f, void** data)
{
    FaceSplitTempData* td = new FaceSplitTempData;

    td->v0 = st.m_mesh.get_triangle(f)[0];
    td->v1 = st.m_mesh.get_triangle(f)[1];
    td->v2 = st.m_mesh.get_triangle(f)[2];

    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);
    td->old_x2 = vc(st.pm_positions[td->v2]);

    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);
    td->old_u2 = vel(td->v2);

    *data = (void*)td;
    callback_log() << "pre facesplit: " << f << ": " << td->v0 << " " << td->v1 << " " << td->v2 << std::endl;
}

void TPFluid::post_facesplit(const LosTopos::SurfTrack& st, size_t f, size_t new_vertex, void* data)
{
    FaceSplitTempData* td = (FaceSplitTempData*)data;
    callback_log() << "post facesplit: " << f << " => " << new_vertex << std::endl;

    vec3 new_x = vc(st.pm_positions[new_vertex]);
    vec3 c = vec3::Zero();
    vec3 n = (td->old_x1 - td->old_x0).cross(td->old_x2 - td->old_x0);
    double nsq = n.squaredNorm();
    c[0] = 1 - (new_x - td->old_x0).dot(n.cross(td->old_x1 - td->old_x2)) / nsq;
    c[1] = 1 - (new_x - td->old_x1).dot(n.cross(td->old_x2 - td->old_x0)) / nsq;
    if (c[0] > 1)        c[0] = 1;
    if (c[0] < 0)        c[0] = 0;
    if (c[1] > 1 - c[0]) c[1] = 1 - c[0];
    if (c[1] < 0)        c[1] = 0;
    c[2] = 1 - c[0] - c[1];

    vel(new_vertex) = td->old_u0 * c[0] + td->old_u1 * c[1] + td->old_u2 * c[2];
}

struct SnapTempData
{
    size_t v0;
    size_t v1;

    vec3 old_x0;
    vec3 old_x1;

    vec3 old_u0;
    vec3 old_u1;
};

void TPFluid::pre_snap(const LosTopos::SurfTrack& st, size_t v0, size_t v1, void** data)
{
    SnapTempData* td = new SnapTempData;
    td->v0 = v0;
    td->v1 = v1;

    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);

    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);

    *data = (void*)td;
    callback_log() << "pre snap: " << td->v0 << " " << td->v1 << std::endl;
}

void TPFluid::post_snap(const LosTopos::SurfTrack& st, size_t v_kept, size_t v_deleted, void* data)
{
    SnapTempData* td = (SnapTempData*)data;
    callback_log() << "post snap: " << td->v0 << " " << td->v1 << " => " << v_kept << std::endl;
    assert((td->v0 == v_kept && td->v1 == v_deleted) || (td->v1 == v_kept && td->v0 == v_deleted));
    assert(v_kept != v_deleted);
    //    assert(st.m_mesh.vertex_is_deleted(v_deleted));
    //    assert(!st.m_mesh.vertex_is_deleted(v_kept));

    vec3 merged_x = vc(st.pm_positions[v_kept]);
    double s = (merged_x - td->old_x0).dot(td->old_x1 - td->old_x0) / (td->old_x1 - td->old_x0).squaredNorm();
    if (s > 1) s = 1;
    if (s < 0) s = 0;
    vec3 new_u = td->old_u0 * (1 - s) + td->old_u1 * s;

    vel(v_kept) = new_u;
}

void TPFluid::pre_smoothing(const LosTopos::SurfTrack& st, void** data)
{

}

void TPFluid::post_smoothing(const LosTopos::SurfTrack& st, void* data)
{

}
