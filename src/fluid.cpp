#include "fluid.h"

TPFluid::TPFluid(LosTopos::SurfTrack* surface, vector<vec3>& velocity, struct FluidOptions& options) : surf(surface), options(options), simTime(0.0) {
    vels = new LosTopos::NonDestructiveTriMesh::VertexData<vec3>(&(surf->m_mesh));
    for (size_t i = 0; i < surf->m_mesh.nv(); i++)
        (*vels)[i] = velocity[i];

    surf->m_solid_vertices_callback = this;
    surf->m_mesheventcallback = this;
}

TPFluid::~TPFluid() {
    delete surf;
}

void TPFluid::render() {
    std::cout << "Rendering..." << std::endl;
    vector<RenderVertex> vertices;
    vector<vec3ui> indices;

    // Get vertex positions
    vector<vec3> xs(surf->m_mesh.nv());
    for (size_t i = 0; i < surf->m_mesh.nv(); i++)
        xs[i] = vc(surf->pm_positions[i]);

    // Compute normals at each vertex
    vector<vec3> vn(surf->m_mesh.nv(), vec3(0, 0, 0));
    for (size_t i = 0; i < surf->m_mesh.nt(); i++) {
        LosTopos::Vec3st t = surf->m_mesh.get_triangle(i);
        vec3 x0 = xs[t[0]];
        vec3 x1 = xs[t[1]];
        vec3 x2 = xs[t[2]];

        vec3 nt = (x1 - x0).cross(x2 - x0);
        if (surf->m_mesh.get_triangle_label(i)[0] < surf->m_mesh.get_triangle_label(i)[1])
            nt = -nt;

        vn[t[0]] += nt;
        vn[t[1]] += nt;
        vn[t[2]] += nt;
    }

    // Generate render vertices
    for (size_t i = 0; i < surf->m_mesh.nv(); i++) {
        vn[i].normalize();
        RenderVertex rv;
        rv.position = vec3f(xs[i][0], xs[i][1], xs[i][2]);
        rv.normal = vec3f(vn[i][0], vn[i][1], vn[i][2]);
        vertices.push_back(rv);
    }

    // Get triangle indices
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < surf->m_mesh.nt(); i++) {
        LosTopos::Vec3st t = surf->m_mesh.get_triangle(i);
        indices.push_back(vec3ui(t[0], t[1], t[2]));
        LosTopos::Vec3d x0 = vc(xs[t[0]]);
        LosTopos::Vec3d x1 = vc(xs[t[1]]);
        LosTopos::Vec3d x2 = vc(xs[t[2]]);
        LosTopos::Vec3d n0 = vc(vn[t[0]]);
        LosTopos::Vec3d n1 = vc(vn[t[1]]);
        LosTopos::Vec3d n2 = vc(vn[t[2]]);
        glNormal3d(n0[0], n0[1], n0[2]);    glVertex3d(x0[0], x0[1], x0[2]);
        glNormal3d(n1[0], n1[1], n1[2]);    glVertex3d(x1[0], x1[1], x1[2]);
        glNormal3d(n2[0], n2[1], n2[2]);    glVertex3d(x2[0], x2[1], x2[2]);
    }
    glEnd();

    /*VertexArrayObject vao;
    vao.init(vertices, indices);
    vao.draw();
    vao.remove();*/
}

void TPFluid::step() {
    double actual_dt = advect(options.dt);
    
    std::cout << "Improving mesh..." << std::endl;
    improveMesh();

    std::cout << "Partitioning mesh..." << std::endl;
    Partitioning partitioning;
    fullPartition(partitioning);

    // mesh no longer changes, so we can have a more useful, static data structure
    vecX v = vecX::Zero(nv() * 3);
    for (size_t i = 0; i < nv(); i++)
        v.segment<3>(i * 3) = (*vels)[i];

    std::cout << "Helmholtz Decomposition..." << std::endl;
    helmholtzDecomposition(actual_dt, v, partitioning);
    
    double smoothing_coef = options.smoothing_coefficient * actual_dt - 0.5;
    std::cout << "Smoothing velocity..." << std::endl;
    smoothVelocity(actual_dt, v, partitioning, smoothing_coef);

    // concaveSmoothing(...);

    std::cout << "Add gravity..." << std::endl;
    addGravity(actual_dt, v);

    std::cout << "Pressure solve..." << std::endl;
    solvePressure(actual_dt, v, partitioning);

    std::cout << "Final velocity..." << std::endl;
    updateFinalVelocity(actual_dt, v);

    simTime += actual_dt;
}

void TPFluid::renderBBox() {
    vec3 org = options.bbox_center - options.bbox_dim / 2;
    vec3 fin = options.bbox_center + options.bbox_dim / 2;
    if (!options.bbox && options.floorOnly)
        org[1] = options.floorz;

#define OG org[0], org[1], org[2]
#define PX fin[0], org[1], org[2]
#define PY org[0], fin[1], org[2]
#define PZ org[0], org[1], fin[2]
#define PXPY fin[0], fin[1], org[2]
#define PXPZ fin[0], org[1], fin[2]
#define PYPZ org[0], fin[1], fin[2]
#define PXPYPZ fin[0], fin[1], fin[2]

    glLineWidth(1);
    glColor3f(0.5, 0.5, 0.5);
    glBegin(GL_LINES);
    if (options.bbox) {
        glVertex3f(OG);
        glVertex3f(PZ);
        glVertex3f(PX);
        glVertex3f(PXPZ);
        glVertex3f(PY);
        glVertex3f(PYPZ);
        glVertex3f(PXPY);
        glVertex3f(PXPYPZ);

        glVertex3f(OG);
        glVertex3f(PX);
        glVertex3f(PX);
        glVertex3f(PXPY);
        glVertex3f(PXPY);
        glVertex3f(PY);
        glVertex3f(PY);
        glVertex3f(OG);

        glVertex3f(PZ);
        glVertex3f(PXPZ);
        glVertex3f(PXPZ);
        glVertex3f(PXPYPZ);
        glVertex3f(PXPYPZ);
        glVertex3f(PYPZ);
        glVertex3f(PYPZ);
        glVertex3f(PZ);
    }
    else if (options.floorOnly) {
        glVertex3f(OG);
        glVertex3f(PX);
        glVertex3f(PX);
        glVertex3f(PXPZ);
        glVertex3f(PXPZ);
        glVertex3f(PZ);
        glVertex3f(PZ);
        glVertex3f(OG);
    }
    glEnd();
}