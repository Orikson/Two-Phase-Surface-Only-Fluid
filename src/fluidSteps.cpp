#include "fluid.h"

//**********************************************************
// Mesh advection
//**********************************************************
double TPFluid::advect(double dt) {
    double actual_dt = 0;

    //**********************************************************
    // Remove Solid Influx
    //**********************************************************
    double influx = currentInflux();
    vecX nsolid = vecX::Zero(nv() * 3);
    if (influx != 0) {
        for (size_t i = 0; i < nf(); i++) {
            if (faceIsSolid(i)) {
                vec3 n = faceOutwardNormal(i);
                LosTopos::Vec3st t = mesh().m_tris[i];
                nsolid.segment<3>(t[0] * 3) = vec3(nsolid.segment<3>(t[0] * 3) + n).normalized();
                nsolid.segment<3>(t[1] * 3) = vec3(nsolid.segment<3>(t[1] * 3) + n).normalized();
                nsolid.segment<3>(t[2] * 3) = vec3(nsolid.segment<3>(t[2] * 3) + n).normalized();
            }
        }

        for (size_t i = 0; i < nv(); i++) {
            if (vertexIsSolid(i)) {
                vec3 n = nsolid.segment<3>(i * 3);
                assert(n.squaredNorm() != 0);
                vel(i) = vel(i) - vel(i).dot(n) * n; // now the velocity of solid vertices represents the solid velocity, not the liquid velocity
            }
        }
    }

    //**********************************************************
    // Mesh advection
    //**********************************************************
    for (size_t i = 0; i < nv(); i++)
        surf->pm_newpositions[i] = surf->pm_positions[i] + vc(vel(i)) * dt;

    surf->rebuild_continuous_broad_phase();
    surf->integrate(dt, actual_dt); // advection
    if (actual_dt != dt)
        std::cout << "Warning: SurfTrack::integrate() failed to step the full length of the time step!" << std::endl;

    // test for floor contact
    bool solid_integrate_needed = false;
    double floorz = options.floorz;
    std::vector<bool> floor_contact(nv(), false);
    for (size_t i = 0; i < nv(); i++)
        if (surf->pm_newpositions[i][2] <= floorz) {
            if (surf->pm_newpositions[i][2] < floorz)
                solid_integrate_needed = true;  // if a vertex has come strictly below the floor, it'll be moved onto the floor and thus another SurfTrack()::integrate() call will be necessary
            surf->pm_newpositions[i][2] = floorz;
            surf->m_masses[i] = LosTopos::Vec3d(1, 1, std::numeric_limits<double>::infinity());
        }

    // unlabel any single-vertex-solid-contact as not solid, because the pressure solve assumes any solid contact region will have at least one face.
    for (size_t i = 0; i < nv(); i++) {
        bool solid_face = false;
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
            if (faceIsSolid(mesh().m_vertex_to_triangle_map[i][j]))
                solid_face = true;
        if (!solid_face)
            surf->m_masses[i] = LosTopos::Vec3d(1, 1, 1);   // a solid vertex incident to no solid faces is not a true solid vertex.
    }

    if (solid_integrate_needed) {
        surf->rebuild_continuous_broad_phase();
        surf->integrate(dt, actual_dt);
    }

    for (size_t i = 0; i < nv(); i++)
        surf->pm_velocities[i] = vc(vel(i));

    //**********************************************************
    // Replace solid influx
    //**********************************************************
    if (influx != 0) {
        for (size_t i = 0; i < nv(); i++) {
            if (vertexIsSolid(i)) {
                vec3 n = nsolid.segment<3>(i * 3);
                assert(n.squaredNorm() != 0);
                vel(i) = vel(i) - vel(i).dot(n) * n - influx * n; // now the velocity of solid vertices represents the liquid velocity, not the solid velocity
            }
        }
    }

    return actual_dt;
}

//**********************************************************
// Mesh improvement
//**********************************************************
void TPFluid::improveMesh() {
    for (size_t i = 0; i < nv(); i++)
        surf->pm_velocities[i] = vc(vel(i));
    for (int i = 0; i < options.remesh_iter; i++) {
        // recompute the target edge lengths (these will be maintained incrementally when performing the remeshing operations, but incremental updates may be too conservative so they are recomputed from scratch here)
        surf->compute_all_vertex_target_edge_lengths();
        surf->topology_changes();
        surf->compute_all_vertex_target_edge_lengths();
        surf->improve_mesh();
    }

    std::vector<size_t> dummy;
    surf->defrag_mesh_from_scratch(dummy);
    
    // confirm floor contact (there may be cases where a vertex is exactly on the floor but can't be labeled as solid by the floor detection code above integrate() because it'd be a branch from the triple junction curve (i.e. it is incident to no solid face) and the remeshing steps above may have made it okay to label it as solid. if we don't label it as solid in such cases, the pressure solve will have trouble because that the triple junction will have a vertex whose air and solid normals are identical (because the air faces are really already on the floor; they just haven't been labeled as solid yet), which makes the triple junciton tangent NaN and the BEM solve badly conditioned. see the crash in 1449952789)
    double floorz = options.floorz;
    for (size_t i = 0; i < nv(); i++)
        if (surf->pm_positions[i][2] == floorz)
            surf->m_masses[i] = LosTopos::Vec3d(1, 1, std::numeric_limits<double>::infinity());

    // unlabel any single-vertex-solid-contact as not solid, because the pressure solve assumes any solid contact region will have at least one face.
    for (size_t i = 0; i < nv(); i++) {
        bool solid_face = false;
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
            if (faceIsSolid(mesh().m_vertex_to_triangle_map[i][j]))
                solid_face = true;
        if (!solid_face)
            surf->m_masses[i] = LosTopos::Vec3d(1, 1, 1);   // a solid vertex incident to no solid faces is not a true solid vertex.
    }

    m_pre_stepping_geometry = vecX::Zero(nv() * 3);
    for (size_t i = 0; i < nv(); i++)
        m_pre_stepping_geometry.segment<3>(i * 3) = pos(i);
}

//**********************************************************
// Mesh partitioning
//**********************************************************
void TPFluid::partitionMesh(Partitioning& partitioning) {
    std::vector<std::vector<size_t>>& p2v = partitioning.p2v;
    std::vector<std::vector<size_t>>& p2f = partitioning.p2f;
    std::vector<int>& v2p = partitioning.v2p;
    std::vector<int>& f2p = partitioning.f2p;
    std::vector<size_t>& flattened_partition_vertices = partitioning.flattened_partition_vertices;
    std::vector<size_t>& indices_in_partitions = partitioning.indices_in_partitions;

    p2v.clear();
    p2f.clear();
    v2p.assign(nv(), -1);
    f2p.assign(nf(), -1);

    // build the partition of vertices first
    while (true) {
        size_t seed_vertex = nv();
        for (size_t i = 0; i < nv(); i++)
            if (v2p[i] < 0) {
                seed_vertex = i;
                break;
            }
        if (seed_vertex == nv())    // can't find a vertex that doesn't belong to any partition. we're done.
            break;

        int new_partition_id = (int)p2v.size();
        p2v.push_back(std::vector<size_t>());

        std::vector<size_t> openset;
        openset.push_back(seed_vertex);
        while (openset.size() > 0) {   // DPS to traverse this paritition, finding all its vertices
            size_t current = openset.back();
            openset.pop_back();

            if (v2p[current] >= 0) {
                assert(v2p[current] == new_partition_id);
                continue;
            }

            v2p[current] = new_partition_id;
            p2v.back().push_back(current);

            for (size_t i = 0; i < mesh().m_vertex_to_edge_map[current].size(); i++) {
                LosTopos::Vec2st e = mesh().m_edges[mesh().m_vertex_to_edge_map[current][i]];
                size_t vother = (e[0] == current ? e[1] : e[0]);

                openset.push_back(vother);
            }
        }
    }

    // sanity check on the vertex counts
    assert(p2v.size() > 0);
    size_t partition_vertex_sum = 0;
    for (size_t i = 0; i < p2v.size(); i++)
        partition_vertex_sum += p2v[i].size();
    assert(partition_vertex_sum == nv());

    // build the partition of faces from the partition of vertices
    p2f.resize(p2v.size());
    for (size_t i = 0; i < nf(); i++) {
        LosTopos::Vec3st t = mesh().m_tris[i];
        int p = v2p[t[0]];
        assert(p == v2p[t[1]]);
        assert(p == v2p[t[2]]);
        assert(p >= 0);
        f2p[i] = p;
        p2f[p].push_back(i);
    }

    // sanity check on the face counts
    assert(p2f.size() > 0);
    size_t partition_face_sum = 0;
    for (size_t i = 0; i < p2f.size(); i++)
        partition_face_sum += p2f[i].size();
    assert(partition_face_sum == nf());

    for (size_t i = 0; i < nv(); i++)
        assert(v2p[i] >= 0);
    for (size_t i = 0; i < nf(); i++)
        assert(f2p[i] >= 0);

    // build the flattened array of vertices which helps with parallelization
    flattened_partition_vertices.clear();
    flattened_partition_vertices.reserve(nv());
    indices_in_partitions.clear();
    indices_in_partitions.reserve(nv());
    for (size_t pi = 0; pi < p2v.size(); pi++) {
        flattened_partition_vertices.insert(flattened_partition_vertices.end(), p2v[pi].begin(), p2v[pi].end());
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
            indices_in_partitions.push_back(pvi);
    }
    assert(flattened_partition_vertices.size() == nv());
    assert(indices_in_partitions.size() == nv());
}

//**********************************************************
// Mesh partitioning with bubble removal
//**********************************************************
void TPFluid::fullPartition(Partitioning& partitioning) {
    partitionMesh(partitioning);

    // find the partition with the largest positive volume: it may be used below
    size_t largest_partition = 0;
    double largest_volume = 0;
    for (size_t pi = 0; pi < partitioning.p2v.size(); pi++) {
        double volume = partitionVolume(partitioning.p2v[pi]);
        if (volume > largest_volume)
            largest_volume = volume,
            largest_partition = pi;
    }

    bool deletion_happened = false; 
    bool merging_happened = false;
    for (size_t pi = 0; pi < partitioning.p2f.size(); pi++) {
        double volume = partitionVolume(partitioning.p2f[pi]);

        bool outofview = true;
        double viewrange = std::abs(options.view_range);
        for (size_t i = 0; i < partitioning.p2v[pi].size(); i++)
            if (pos(partitioning.p2v[pi][i]).norm() < viewrange)
                outofview = false;

        if ((options.ibr && volume < 0) || outofview) {
            // IBR
            // remove any closed interior bubbles (isolated closed surfaces with negative liquid volume). They are usually by-products of merging and are created at close to zero volumes, but due to our air boundary condition, our solver doesn't have the ability to handle them in a meaningful way.
            // also remove any droplet that's far away from the center of the scene
            for (size_t i = 0; i < partitioning.p2f[pi].size(); i++)
                surfTrack()->remove_triangle(partitioning.p2f[pi][i]);
            for (size_t i = 0; i < partitioning.p2v[pi].size(); i++)
                surfTrack()->remove_vertex(partitioning.p2v[pi][i]);
            deletion_happened = true;

            // clear this partition
            for (size_t i = 0; i < partitioning.p2v[pi].size(); i++)
                partitioning.v2p[partitioning.p2v[pi][i]] = -1;
            for (size_t i = 0; i < partitioning.p2f[pi].size(); i++)
                partitioning.f2p[partitioning.p2f[pi][i]] = -1;
            partitioning.p2v[pi].clear();
            partitioning.p2f[pi].clear();

        } else if (volume < 0) {
            // this is the case of an internal bubble, but it hasn't been deleted in the above branch because IBR is turned off.
            // although we don't have a theory for this, we should do it by treating the entire connected liquid body together, i.e. doing HD on the inward facing (internal
            //  bubble) partitions and the outward facing partition together. this is because the outward facing partition's interior velocity field is not zero, and it affects
            //  the internal bubble partitions. at least this will get rid of the sudden change of behavior upon closing of a concavity which forms a bubble.
            // TODO: ideally we should be able to do a point-mesh inside/outside test to figure out which outward facing partition contains this inward facing partition.
            // For now we don't have time to do this. Let's just say the outward facing partition containing this bubble is the one with the largest positive volume.
            merging_happened = true;
            assert(pi != largest_partition);    // a negative volume partition can't be the largest partition.

            partitioning.p2v[largest_partition].insert(partitioning.p2v[largest_partition].end(), partitioning.p2v[pi].begin(), partitioning.p2v[pi].end());
            partitioning.p2f[largest_partition].insert(partitioning.p2f[largest_partition].end(), partitioning.p2f[pi].begin(), partitioning.p2f[pi].end());
            for (size_t i = 0; i < partitioning.p2v[pi].size(); i++)
                partitioning.v2p[partitioning.p2v[pi][i]] = largest_partition;
            for (size_t i = 0; i < partitioning.p2f[pi].size(); i++)
                partitioning.f2p[partitioning.p2f[pi][i]] = largest_partition;
            partitioning.p2v[pi].clear();
            partitioning.p2f[pi].clear();
        }
    }

    if (deletion_happened) {
        // if there has been deletion, the mesh needs to be defragged; as a result, the partitioning's indices will be invalidated and need to be mapped/rebuilt
        vector<size_t> defrag_vmap;
        defrag_vmap.reserve(nv());
        for (size_t i = 0; i < nv(); i++)
            defrag_vmap.push_back(i);
        surf->defrag_mesh_from_scratch(defrag_vmap);

        // rebuild p2v, p2f, v2p, f2p
        // map the vertex indices
        for (size_t pi = 0; pi < partitioning.p2v.size(); pi++)
            for (size_t i = 0; i < partitioning.p2v[pi].size(); i++)
            {
                partitioning.p2v[pi][i] = defrag_vmap[partitioning.p2v[pi][i]];
                assert(partitioning.p2v[pi][i] >= 0);   // if a vertex is deleted, its partition should have already been emptied, so this shouldn't happen
            }

        // rebuild the inverse maps
        std::vector<int> f2p(nf(), -1);
        std::vector<int> v2p(nf(), -1);
        for (size_t pi = 0; pi < partitioning.p2v.size(); pi++)
            for (size_t i = 0; i < partitioning.p2v[pi].size(); i++) {
                v2p[partitioning.p2v[pi][i]] = pi;
                for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[partitioning.p2v[pi][i]].size(); j++)
                    f2p[mesh().m_vertex_to_triangle_map[partitioning.p2v[pi][i]][j]] = pi;
            }

        // rebuild the face indices
        for (size_t pi = 0; pi < partitioning.p2f.size(); pi++)
            partitioning.p2f[pi].clear();
        for (size_t i = 0; i < nf(); i++) {
            assert(f2p[i] >= 0);    // if a face is deleted, its partition should have alrady been emptied, so this shouldn't happen
            partitioning.p2f[f2p[i]].push_back(i);
        }

        partitioning.v2p = v2p;
        partitioning.f2p = f2p;
    }

    if (deletion_happened || merging_happened) {
        // merging doesn't require rebuilding p2v, p2f, v2p, f2p like deletion does, but both require rebuilding the flattened array (flattened_partition_vertices) and its inverse map (indices_in_partitions)
        partitioning.flattened_partition_vertices.clear();
        partitioning.flattened_partition_vertices.reserve(nv());
        partitioning.indices_in_partitions.clear();
        partitioning.indices_in_partitions.reserve(nv());
        for (size_t pi = 0; pi < partitioning.p2v.size(); pi++) {
            partitioning.flattened_partition_vertices.insert(partitioning.flattened_partition_vertices.end(), partitioning.p2v[pi].begin(), partitioning.p2v[pi].end());
            for (size_t pvi = 0; pvi < partitioning.p2v[pi].size(); pvi++)
                partitioning.indices_in_partitions.push_back(pvi);
        }
        assert(partitioning.flattened_partition_vertices.size() == nv());
        assert(partitioning.indices_in_partitions.size() == nv());
    }
}

//**********************************************************
// Helmholtz Decomposition to project back onto a divergence free velocity field
//**********************************************************
void TPFluid::helmholtzDecomposition(double dt, vecX& v, Partitioning& partitioning) {
    const vector<vector<size_t>>& p2v = partitioning.p2v;
    const vector<vector<size_t>>& p2f = partitioning.p2f;
    const vector<int>& v2p = partitioning.v2p;
    const vector<int>& f2p = partitioning.f2p;
    const vector<size_t>& flattened_partition_vertices = partitioning.flattened_partition_vertices;
    const vector<size_t>& indices_in_partitions = partitioning.indices_in_partitions;

    // remove the global (rigid translation) component, improving accuracy by only reconstructing a smaller component
    vector<vec3> vglobal(p2v.size(), vec3::Zero());
    for (size_t pi = 0; pi < p2v.size(); pi++) {
        double areasum = 0;
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++) {
            size_t i = p2v[pi][pvi];
            vglobal[pi] += v.segment<3>(i * 3) * vertexArea(i);
            areasum += vertexArea(i);
        }
        vglobal[pi] /= areasum;
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++) {
            size_t i = p2v[pi][pvi];
            v.segment<3>(i * 3) -= vglobal[pi];
        }
    }

    for (size_t i = 0; i < nv(); i++)
        vel(i) = v.segment<3>(i * 3);

    // cache geometry some data
    vector<double> face_areas(nf());
    vector<vec3> face_outward_normals(nf());
    vector<vec3> face_centers(nf());
    vector<vec3> face_center_velocities(nf());
    for (size_t i = 0; i < nf(); i++)
    {
        LosTopos::Vec3st t = mesh().m_tris[i];

        face_areas[i] = faceArea(i);
        face_outward_normals[i] = faceOutwardNormal(i);
        face_centers[i] = (pos(t[0]) + pos(t[1]) + pos(t[2])) / 3;
        face_center_velocities[i] = (v.segment<3>(t[0] * 3) + v.segment<3>(t[1] * 3) + v.segment<3>(t[2] * 3)) / 3;
    }

    vec3 y;
    double jacobian_j;
    const double oneover4pi = 1 / (4 * M_PI);

    // compute the gradient field component
    vector<vec3> dPhidn(nf(), vec3::Zero()); // one integral for each face as a neighbor of each vertex

    const std::vector<vec2> quadrature_line = BoundaryIntegral::quadrature_line();
    const std::vector<vec3> quadrature_square = BoundaryIntegral::quadrature_square();
    const std::vector<vec3> quadrature_triangle = BoundaryIntegral::quadrature_triangle();
    const vec3 qt0 = quadrature_triangle[0];

#pragma omp parallel default (none) shared (quadrature_square, dPhidn, face_outward_normals, face_centers, face_areas, face_center_velocities, p2v, p2f, v2p, f2p, flattened_partition_vertices, indices_in_partitions, std::cout) private (y, jacobian_j)
    {
#pragma omp master
        {
            //        std::cout << "HD dPhidn loop: spawning " << omp_get_num_threads() << " threads." << std::endl;
        }

#pragma omp for schedule (dynamic, 10)
        for (long fpi = 0; fpi < nv(); fpi++)   // degree of freedom on vertex i
        {
            size_t i = flattened_partition_vertices[fpi];
            size_t pvi = indices_in_partitions[fpi];
            size_t pi = v2p[i];

            for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++) {
                size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];

                LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
                mat3 xis = getVertexPositions(fi);
                vec3 n_i = face_outward_normals[face_i];

                double I = 0;
                for (size_t qi = 0; qi < quadrature_square.size(); qi++) {
                    vec2 qik = quadrature_square[qi].segment<2>(0);   // coordinates in the square ref domain
                    double qiw = quadrature_square[qi].z();

                    vec3 x;
                    double jacobian_i;
                    vec3 c_i;
                    BoundaryIntegral::duffyTransform(xis, 0, qik, x, jacobian_i, c_i);

                    double theta_i = c_i[0];    // the pyramid function for vertex i

                    double I_dSLP = 0;
                    for (size_t pfj = 0; pfj < p2f[pi].size(); pfj++) { // inner integral: SLP normal derivative
                        size_t j = p2f[pi][pfj];

                        if (j == face_i)
                            continue;   // the inner integral is the normal derivative of SLP, which receives no contribution from the face where the evaluation point x is located in

                        vec3 n_j = face_outward_normals[j];

                        y = face_centers[j];
                        vec3 dx = x - y;
                        double dxn = dx.norm();
                        double dGdx = dx.dot(n_i) / (dxn * dxn * dxn) * oneover4pi;

                        I_dSLP += face_areas[j] * face_center_velocities[j].dot(n_j) * dGdx;
                    }

                    I += qiw * jacobian_i * theta_i * I_dSLP;
                }

                LosTopos::Vec2st dummy;
                dPhidn[face_i][mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy)] = I;
            }
        }
        //#pragma omp for

    }
    //#pragma omp parallel

    // compute the curl field component
    vector<vec3> curl_A(nf(), vec3::Zero());  // one integral for each face as a neighbor of each vertex

#pragma omp parallel default (none) shared (quadrature_square, curl_A, face_outward_normals, face_centers, face_areas, face_center_velocities, p2v, p2f, v2p, f2p, flattened_partition_vertices, indices_in_partitions, std::cout) private (y, jacobian_j)
    {
#pragma omp master
        {
            //        std::cout << "HD curl_A loop: spawning " << omp_get_num_threads() << " threads." << std::endl;
        }

#pragma omp for schedule (dynamic, 10)
        // TODO: merge this loop with the one for dPhidn
        for (long fpi = 0; fpi < nv(); fpi++) {  // degree of freedom on vertex i
            size_t i = flattened_partition_vertices[fpi];
            size_t pvi = indices_in_partitions[fpi];
            size_t pi = v2p[i];

            // integrate over the faces
            for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++) {
                size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];

                LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
                mat3 xis = getVertexPositions(fi);
                mat3 vis = getVertexVelocities(fi);
                vec3 n_i = faceOutwardNormal(face_i);
                double area_fi = faceArea(face_i);

                double I = 0;
                for (size_t qi = 0; qi < quadrature_square.size(); qi++) {
                    vec2 qik = quadrature_square[qi].segment<2>(0);   // coordinates in the square ref domain
                    double qiw = quadrature_square[qi].z();

                    vec3 x;
                    double jacobian_i;
                    vec3 c_i;
                    BoundaryIntegral::duffyTransform(xis, 0, qik, x, jacobian_i, c_i);

                    double theta_i = c_i[0];
                    vec3 grad_theta_i = (xis.col(1) - xis.col(0)).cross(xis.col(2) - xis.col(0)).normalized().cross(xis.col(2) - xis.col(1)) / (area_fi * 2);

                    vec3 v_x = vis * c_i;  // find the velocity at x by interpolation

                    vec3 I_SLP = vec3::Zero();
                    for (size_t pfj = 0; pfj < p2f[pi].size(); pfj++) {  // inner integral: SLP
                        size_t j = p2f[pi][pfj];

                        if (j == face_i) {
                            // for the coincident case (j == face_i), subdivide the face at position x to create three faces, each integrated by quadrature individually
                            // otherwise, treat the entire face as a whole
                            for (int k = 0; k < 3; k++) {
                                LosTopos::Vec3st fj = mesh().m_tris[j];
                                mat3 xjs = getVertexPositions(fj);
                                mat3 vjs = getVertexVelocities(fj);
                                vec3 n_j = faceOutwardNormal(j);

                                xjs.col(k) = x;     // override one of the vertices because this triangle is one of the three subdivision triangles of face j
                                vjs.col(k) = v_x;

                                for (size_t qj = 0; qj < quadrature_square.size(); qj++)
                                {
                                    vec2 qjk = quadrature_square[qj].segment<2>(0);   // coordinates in the square ref domain
                                    double qjw = quadrature_square[qj].z();

                                    vec3 y;
                                    double jacobian_j;
                                    vec3 c_j;
                                    BoundaryIntegral::duffyTransform(xjs, k, qjk, y, jacobian_j, c_j);   // use duffy transform only for the coincident case

                                    vec3 v_y = vjs * c_j;  // find the velocity at y by interpolation

                                    I_SLP += qjw * jacobian_j * n_j.cross(v_y) * BoundaryIntegral::G(x, y);
                                }
                            }
                        } else {
                            vec3 n_j = face_outward_normals[j];

                            vec3 y = face_centers[j];
                            vec3 dx = x - y;
                            double dxn = dx.norm();
                            double G = -oneover4pi / dxn;

                            I_SLP += face_areas[j] * n_j.cross(face_center_velocities[j]) * G;
                        }
                    }

                    I += -qiw * jacobian_i * n_i.dot(grad_theta_i.cross(I_SLP));
                }

                LosTopos::Vec2st dummy;
                curl_A[face_i][mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy)] = I;
            }

        }
        //#pragma omp for

    }
    //#pragma omp parallel

    // the jump term (i.e. half the SLP charge in the scalar potential term)
    vector<vec3> charge(nf(), vec3::Zero());  // one integral for each face as a neighbor of each vertex

    // TODO: merge this loop with the ones for dPhidn and curl_A
    for (size_t i = 0; i < nv(); i++) {
        for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++) {
            size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];

            LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
            mat3 xis = getVertexPositions(fi);
            mat3 vis = getVertexVelocities(fi);
            vec3 n_i = faceOutwardNormal(face_i);
            double area_fi = faceArea(face_i);

            double I = -(vis.col(0) / 2 + vis.col(1) / 4 + vis.col(2) / 4).dot(n_i) * (area_fi / 3) / 2;   // jump is half the charge

            LosTopos::Vec2st dummy;
            charge[face_i][mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy)] = I;
        }
    }
    
    // composition
    vecX v_dPhidn = vecX::Zero(nv() * 3);
    vecX v_charge = vecX::Zero(nv() * 3);
    vecX v_curl_A = vecX::Zero(nv() * 3);
    vecX v_n = vecX::Zero(nv() * 3);
    vecX v_t = vecX::Zero(nv() * 3);
    vecX newv = v;
    for (size_t i = 0; i < nv(); i++) {
        double area_vi = vertexArea(i);

        vec3 n_basis = vec3::Zero();
        for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++) {
            size_t fi = mesh().m_vertex_to_triangle_map[i][ii];
            vec3 fn = faceOutwardNormal(fi);
            double fa = faceArea(fi);
            n_basis += fn * fa / 3; // n_basis is the volume gradient, or equivalently, the area-weighted average of face normals
        }
        // note that n_basis is not normalized here, retaining its magnitude as volume gradient

        double n_component = 0;             // the normal component (scalar) integrated over the vertex neighborhood (weighted by theta_i), interpreted as volume change rate
        vec3 v_integral = vec3::Zero();   // old velocity (full vector) integrated over the vertex neighborhood (weighted by theta_i)
        double nc_dPhidn = 0;
        double nc_charge = 0;
        double nc_curl_A = 0;
        for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++) {
            size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];

            LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
            mat3 xis = getVertexPositions(fi);
            mat3 vis = getVertexVelocities(fi);
            vec3 n_i = faceOutwardNormal(face_i);
            double area_fi = faceArea(face_i);

            LosTopos::Vec2st dummy;
            size_t vi = mesh().index_in_triangle(mesh().m_tris[face_i], i, dummy);

            //            n_component += (curl_A[face_i][vi] - (dPhidn[face_i][vi] + charge[face_i][vi]));   // the normal component (integral), which is also the volume change rate
            nc_dPhidn += -dPhidn[face_i][vi];
            nc_charge += -charge[face_i][vi];
            nc_curl_A += curl_A[face_i][vi];

            v_integral += (vis.col(0) / 2 + vis.col(1) / 4 + vis.col(2) / 4) * (area_fi / 3);   // the integral of old velocity times theta_i
        }

        n_component = nc_dPhidn + nc_charge + nc_curl_A;

        mat3 P = mat3::Identity() - n_basis * n_basis.transpose() / n_basis.squaredNorm();    // tangent projector

        vec3 vt = P * v_integral / area_vi;
        vec3 vn = n_component * n_basis / n_basis.squaredNorm();

        v_dPhidn.segment<3>(i * 3) = nc_dPhidn * n_basis / n_basis.squaredNorm();
        v_charge.segment<3>(i * 3) = nc_charge * n_basis / n_basis.squaredNorm();
        v_curl_A.segment<3>(i * 3) = nc_curl_A * n_basis / n_basis.squaredNorm();
        v_n.segment<3>(i * 3) = vn;
        v_t.segment<3>(i * 3) = vt;

        newv.segment<3>(i * 3) = vt + vn;
    }

    v = newv;

    // added back the previously removed global (rigid translation) component
    for (size_t pi = 0; pi < p2v.size(); pi++) {
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++) {
            size_t i = p2v[pi][pvi];
            v.segment<3>(i * 3) += vglobal[pi];
        }
    }

    for (size_t i = 0; i < nv(); i++)
        vel(i) = v.segment<3>(i * 3);
}

//**********************************************************
// Velocity smoothing (sharpening)
//**********************************************************
void TPFluid::smoothVelocity(double dt, vecX& v, Partitioning& partitioning, double coef) {
    // velocity smoothing (implicit)
    // u_new = u_old + smoothing_coef * dt * Laplacian u_new
    // the equation is constructed as if u_new and u_old are scalar fields, and solved three times, each with a different coordinate component as a different rhs
    // note that the Laplacian discretization is consistent with the smoothing introduced by HD, so the effect should be exactly canceling if smoothing_coef * dt = -0.5 here.
    // note also that although the incoming arguments include a partitioning, the implementation below ignores it and does a global solve including all connected components of the mesh.
    Eigen::SparseMatrix<double> smoothing_A(nv(), nv());
    Eigen::Matrix<double, Eigen::Dynamic, 3> smoothing_rhs(nv(), 3);   // the three columns are for x, y and z components respectively
    smoothing_rhs.setZero();
    assert(smoothing_rhs.rows() == nv());
    assert(smoothing_rhs.cols() == 3);

    std::vector<Eigen::Triplet<double> > smoothing_A_triplets;
    for (size_t i = 0; i < nv(); i++) {
        double ai = vertexArea(i) * 3;
        for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++) {
            size_t e = mesh().m_vertex_to_edge_map[i][j];
            size_t vother = edgeOtherVertex(e, i);

            assert(mesh().m_edge_to_triangle_map[e].size() == 2);
            double aj = 0;
            for (size_t k = 0; k < mesh().m_edge_to_triangle_map[e].size(); k++)
                aj += faceArea(mesh().m_edge_to_triangle_map[e][k]);
            double w = aj / (ai * 2);

            smoothing_A_triplets.push_back(Eigen::Triplet<double>(i, vother, -w * coef));
        }

        smoothing_A_triplets.push_back(Eigen::Triplet<double>(i, i, 1 + coef));

        smoothing_rhs.row(i) = v.segment<3>(i * 3).transpose();
    }

    smoothing_A.setFromTriplets(smoothing_A_triplets.begin(), smoothing_A_triplets.end());

    vecX newv = vecX::Zero(nv() * 3);
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > smoothing_solver(smoothing_A);
    vecX smoothing_solution = smoothing_solver.solve(smoothing_rhs.col(0));     // the x component of u_new
    for (size_t i = 0; i < nv(); i++) newv[i * 3 + 0] = smoothing_solution[i];
    smoothing_solution = smoothing_solver.solve(smoothing_rhs.col(1));          // the y component of u_new
    for (size_t i = 0; i < nv(); i++) newv[i * 3 + 1] = smoothing_solution[i];
    smoothing_solution = smoothing_solver.solve(smoothing_rhs.col(2));          // the y component of u_new
    for (size_t i = 0; i < nv(); i++) newv[i * 3 + 2] = smoothing_solution[i];

    v = newv;
}

//**********************************************************
// Add gravity
//**********************************************************
void TPFluid::addGravity(double dt, vecX& v) {
    for (size_t i = 0; i < nv(); i++)
        v.segment<3>(i * 3) += options.gravity * dt;
}

//**********************************************************
// Pressure solve
//**********************************************************
void TPFluid::solvePressure(double dt, vecX& v, Partitioning& partitioning) {
    bool TPCF = options.tpcf;
    bool TDMC = options.tdmc;

    const vector<vector<size_t> >& p2v = partitioning.p2v;
    const vector<vector<size_t> >& p2f = partitioning.p2f;
    const vector<int>& v2p = partitioning.v2p;
    const vector<int>& f2p = partitioning.f2p;
    const vector<size_t>& flattened_partition_vertices = partitioning.flattened_partition_vertices;
    const vector<size_t>& indices_in_partitions = partitioning.indices_in_partitions;

    for (size_t i = 0; i < nv(); i++)
        vel(i) = v.segment<3>(i * 3);

    size_t N = nv();

    // compute the boundary conditions
    // note that both p and dpdn here are time integrated (multiplied by dt)
    vecX BC_p = vecX::Zero(N);       // pressure at vertex (set 0 for solid interior vertices, which will not be used anyway)
    vecX BC_dpdn = vecX::Zero(N);    // dpdn at vertex (the solid side of the vertex for triple junction vertices, and 0 for air interior vertices which will not be used anyway)
    vecX BC_dpdna = vecX::Zero(N);   // dpdn at triple junction vertex on the air side, in the amount just enough to keep the triple junction stationary (dpdna = v_tangent.dot(na) * rho = ((I - ns ns^T) v).dot(na) * rho)
    double influx = currentInflux();
    for (size_t i = 0; i < nv(); i++) {
        vec3 n_i = vertOutwardNormal(i);

        if (!vertexIsSolid(i)) {
            // vertex i is on the interior of the liquid-air interface
            //  BC_p[i] = the pressure jump due to surface tension; BC_dpdn[i] = 0 (unused)
            BC_p[i] = vertIntegralMeanCurvature(i) * options.sigma / vertexArea(i) * dt;     // multiplied by dt because the pressure is time integrated
            BC_dpdn[i] = 0;
            BC_dpdna[i] = 0;
            assert(BC_p[i] == BC_p[i]);
        } else if (vertexIsSolid(i) && !vertexIsOnTripleJunction(i)) {
            // vertex i is on the interior of the liquid-solid interface
            //  BC_p[i] = 0 (unused); BC_dpdn[i] = the solid boundary condition
            BC_p[i] = 0;
            BC_dpdn[i] = (v.segment<3>(i * 3).dot(n_i) + influx) * options.rho;
            BC_dpdna[i] = 0;
        } else {
            // vertex i is on triple junction
            //  BC_p[i] = the pressure jump due to the triple junction surface tension; BC_dpdn[i] = the solid side boundary condition
            // the pressure jump comes from two sources: the curvature along the triple junction tangent direction, and the balance of the
            //  three-phase surface tension forces in the plane orthogonal to the triple junction (which can be seen as the curvature along
            //  the direction orthogonal to the triple junction).
            vec3 n_solid = vec3::Zero();
            vec3 n_air = vec3::Zero();
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++) {
                size_t fj = mesh().m_vertex_to_triangle_map[i][j];
                if (faceIsSolid(fj))  n_solid += faceInteriorAngle(fj, i) * faceOutwardNormal(fj);
                else                    n_air += faceInteriorAngle(fj, i) * faceOutwardNormal(fj);
            }
            n_solid.normalize();
            n_air.normalize();

            vec3 triple_junction_tangent = n_solid.cross(n_air).normalized();
            vec3 t_solid_outward = triple_junction_tangent.cross(n_solid).normalized();
            vec3 t_air_outward = -triple_junction_tangent.cross(n_air).normalized();

            // compute the pressure jump in the direction orthogonal to the triple junction (due to combined three-face surface tension forces)
            vec3 surface_tension_force_sl = -t_solid_outward * options.sigma * options.sigma_sl;
            vec3 surface_tension_force_sa = t_solid_outward * options.sigma * options.sigma_sa;
            vec3 surface_tension_force_la = -t_air_outward * options.sigma;

            vec3 surface_tension_force_combined = surface_tension_force_sl + surface_tension_force_sa + surface_tension_force_la;

            double pressure_jump_normal_to_triple_junction = surface_tension_force_combined.dot(-t_solid_outward) / tripleJunctionVirtualWidth();  // project the combined force in the plane orthogonal to the triple junction in the direction of the solid tangent, because the component normal to the solid surface will be canceled by the solid normal force

            // compute the pressure jump in the direction tangent to the triple junction (due to the curvature of the triple junction curve)
            double triple_junction_solid_angle = 0;
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
                if (faceIsSolid(mesh().m_vertex_to_triangle_map[i][j]))
                    triple_junction_solid_angle += faceInteriorAngle(mesh().m_vertex_to_triangle_map[i][j], i);
            double triple_junction_curve_integral_curvature = M_PI - triple_junction_solid_angle;

            std::vector<size_t> triple_junction_edges;
            for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++)
                if (edgeIsOnTripleJunction(mesh().m_vertex_to_edge_map[i][j]))
                    triple_junction_edges.push_back(mesh().m_vertex_to_edge_map[i][j]);
            //            assert(triple_junction_edges.size() == 2);

            double pressure_jump_tangent_to_triple_junction = triple_junction_curve_integral_curvature * options.sigma / ((edgeLength(triple_junction_edges[0]) + edgeLength(triple_junction_edges[1])) / 2);

            if (triple_junction_edges.size() != 2)  // this can only happen transiently when a droplet is beginning to come into contact with a solid surface but hasn't fully formed a contiguous contact region yet
            {
                pressure_jump_tangent_to_triple_junction = 0;   // don't do anything in this case because the system is lacking DOF and there's no right thing to do. hopefully the inertial motion will resolve this situation very soon.
                pressure_jump_normal_to_triple_junction = 0;
            }
            assert(pressure_jump_normal_to_triple_junction == pressure_jump_normal_to_triple_junction);
            assert(pressure_jump_tangent_to_triple_junction == pressure_jump_tangent_to_triple_junction);

            // combined the two pressure jumps to find the pressure Dirichlet BC
            BC_p[i] = (pressure_jump_normal_to_triple_junction + pressure_jump_tangent_to_triple_junction) * dt;
            BC_dpdn[i] = (v.segment<3>(i * 3).dot(n_solid) + influx) * options.rho;
            BC_dpdna[i] = ((mat3::Identity() - n_solid * n_solid.transpose()) * v.segment<3>(i * 3)).dot(n_air) * options.rho;
        }

    }
    
    assert(BC_p == BC_p);
    assert(BC_dpdn == BC_dpdn);
    assert(BC_dpdna == BC_dpdna);

    // collocation equation:
    //  omega_j p_j = \int p dGdny - \int dpdn G
    //  omega_j p_j = \sum p_i \int theta_i dGdny - \sum dpdn_i \int theta_i G
    //  omega_j p_j = collocA_{ij} p_i - collocB_{ij} dpdn_i
    //  omega_j p_j = colloc_A^T p - collocB^T dpdn
    std::vector<matX> A(p2v.size());
    std::vector<vecX> rhs(p2v.size());
    for (size_t pi = 0; pi < p2v.size(); pi++)
    {
        size_t Npi = p2v[pi].size();
        A[pi] = matX::Zero(Npi, Npi);
        rhs[pi] = vecX::Zero(Npi);
    }

    // allocate one VecXd per vertex, so that the OMP loop below doesn't need to access the same variable in different loop iterations which causes a race. these per-vertex rhs are summed into one VecXd (per partition) after the OMP loop.
    std::vector<std::vector<vecX>> rhs_per_vertex(p2v.size());
    for (size_t pi = 0; pi < p2v.size(); pi++) {
        size_t Npi = p2v[pi].size();
        rhs_per_vertex[pi].assign(Npi, vecX::Zero(Npi));
    }

    // temporary data (code level optimization for the loop below)
    const double oneover4pi = 1 / (4 * M_PI);

    const std::vector<vec3> quadrature_square = BoundaryIntegral::quadrature_square();

    std::vector<bool> face_is_solids(nf(), false);
    std::vector<vec3> face_outward_normals(nf(), vec3::Zero());
    std::vector<double> face_areas(nf(), 0);
    for (size_t i = 0; i < nf(); i++)
    {
        face_is_solids[i] = faceIsSolid(i);
        face_outward_normals[i] = faceOutwardNormal(i);
        face_areas[i] = faceArea(i);
    }

    std::vector<vec3> vert_solid_normals(nv(), vec3::Zero());
    std::vector<vec3> vert_air_normals(nv(), vec3::Zero());
    std::vector<double> vert_solid_areas(nv(), 0);
    std::vector<double> vert_air_areas(nv(), 0);
    for (size_t i = 0; i < nv(); i++)
    {
        for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
        {
            size_t face_i = mesh().m_vertex_to_triangle_map[i][j];
            vec3 n = face_outward_normals[face_i];
            double a = face_areas[face_i];
            if (face_is_solids[face_i])
                vert_solid_normals[i] += a * n,
                vert_solid_areas[i] += a;   // don't divide by 3 here, because we're looking at the support of the vertex basis function theta_i
            else
                vert_air_normals[i] += a * n,
                vert_air_areas[i] += a;     // don't divide by 3 here, because we're looking at the support of the vertex basis function theta_i
        }
        if (vert_solid_areas[i] > 0)
            vert_solid_normals[i].normalize();
        if (vert_air_areas[i] > 0)
            vert_air_normals[i].normalize();
    }

#pragma omp parallel default (none) shared (quadrature_square, oneover4pi, A, rhs, rhs_per_vertex, BC_p, BC_dpdn, BC_dpdna, TPCF, p2v, p2f, v2p, f2p, flattened_partition_vertices, indices_in_partitions, vert_solid_normals, vert_air_normals, vert_solid_areas, vert_air_areas, face_areas, face_outward_normals, face_is_solids, std::cout)
    {
#pragma omp master
        {
            //        std::cout << "BEM matrix assembly: spawning " << omp_get_num_threads() << " threads." << std::endl;
        }

#pragma omp for schedule (dynamic, 10)

        for (long fpi = 0; fpi < nv(); fpi++)   // dof on vertex i, where the unknown is either p or dpdn depending on where vertex i is sitting
        {
            size_t i = flattened_partition_vertices[fpi];
            size_t pvi = indices_in_partitions[fpi];
            size_t pi = v2p[i];
            size_t Npi = p2v[pi].size();

            double omega = vertInteriorSolidAngle(i) / (4 * M_PI);

            double mel = 0;
            for (size_t j = 0; j < mesh().m_vertex_to_edge_map[i].size(); j++)
                mel += edgeLength(mesh().m_vertex_to_edge_map[i][j]);
            mel /= mesh().m_vertex_to_edge_map[i].size();

            // one row of the original collocA_solid, collocA_air, collocB_solid and collocB_air matrices, corresponding to DOF i
            // component j in one of these vectors corresponds to collocation point j
            vecX collocA_solid_rowi = vecX::Zero(Npi);
            vecX collocA_air_rowi = vecX::Zero(Npi);
            vecX collocB_solid_rowi = vecX::Zero(Npi);
            vecX collocB_air_rowi = vecX::Zero(Npi);

            for (size_t pvj = 0; pvj < Npi; pvj++) {  // the collocation point is x_j
                size_t j = p2v[pi][pvj];

                vec3 x = pos(j);

                if ((x - pos(i)).norm() > 5 * mel) {
                    // for far-field vertices, use one sample point at the vertex to represent the whole integral domain
                    vec3 dx = pos(i) - x;
                    double dxn = dx.norm();
                    vec3 dGdy = oneover4pi * dx / (dxn * dxn * dxn);
                    double G = -oneover4pi / dxn;

                    if (vert_solid_areas[i] > 0)
                        collocA_solid_rowi[pvj] = vert_solid_areas[i] * dGdy.dot(vert_solid_normals[i]) / 3,  // the factor 1/3 is the integral of the basis function theta_i over an incident face of unit area
                        collocB_solid_rowi[pvj] = vert_solid_areas[i] * G / 3;
                    if (vert_air_areas[i] > 0)
                        collocA_air_rowi[pvj] = vert_air_areas[i] * dGdy.dot(vert_air_normals[i]) / 3,    // the factor 1/3 is the integral of the basis function theta_i over an incident face of unit area
                        collocB_air_rowi[pvj] = vert_air_areas[i] * G / 3;
                } else {
                    // for near-field vertices, use Gaussian quadrature (with Duffy transform)
                    for (size_t ii = 0; ii < mesh().m_vertex_to_triangle_map[i].size(); ii++) {
                        size_t face_i = mesh().m_vertex_to_triangle_map[i][ii];

                        LosTopos::Vec3st fi = getShuffledTriangle(mesh().m_tris[face_i], i);
                        vec3 n_i = face_outward_normals[face_i];
                        double a_i = face_areas[face_i];

                        vec3 x0 = pos(fi[0]);
                        vec3 x1 = pos(fi[1]);
                        vec3 x2 = pos(fi[2]);

                        double Iij = 0;
                        for (size_t qi = 0; qi < quadrature_square.size(); qi++) {
                            vec2 qik = quadrature_square[qi].segment<2>(0);   // coordinates in the square ref domain
                            double qiw = quadrature_square[qi].z();

                            vec3 y;
                            double jacobian_i;
                            vec3 c_i;
                            //                        BoundaryIntegral::duffyTransform(xis, 0, qik, y, jacobian_i, c_i);

                            vec2 t((qik.x() + 1) / 2, (qik.x() + 1) / 2 * (qik.y() + 1) / 2);

                            jacobian_i = t.x() * a_i * 0.5;
                            y = x0 * (1 - t.x()) + x1 * (t.x() - t.y()) + x2 * t.y();

                            double theta_i = 1 - t.x();     // the pyramid function for vertex i

                            vec3 dx = y - x;
                            double dxn = dx.norm();
                            double dGdy = oneover4pi * dx.dot(n_i) / (dxn * dxn * dxn);
                            double G = -oneover4pi / dxn;

                            double Ia = qiw * jacobian_i * theta_i * dGdy;
                            double Ib = qiw * jacobian_i * theta_i * G;

                            if (face_is_solids[face_i])
                                collocA_solid_rowi[pvj] += Ia,
                                collocB_solid_rowi[pvj] += Ib;
                            else
                                collocA_air_rowi[pvj] += Ia,
                                collocB_air_rowi[pvj] += Ib;
                        }
                    }
                }
            }

            if (!vertexIsSolid(i)) {
                // vertex i is in the interior of the liquid-air interface
                //  dpdn is the unknown; p is known
                A[pi].col(pvi) = -(collocB_solid_rowi + collocB_air_rowi);
                rhs_per_vertex[pi][pvi] = -(collocA_solid_rowi + collocA_air_rowi) * BC_p[i];
                rhs[pi][pvi] += omega * BC_p[i];

            } else if (vertexIsSolid(i) && !vertexIsOnTripleJunction(i)) {
                // vertex i is in the interior of the liquid-solid interface
                //  p is the unknown; dpdn is known
                A[pi].col(pvi) = (collocA_solid_rowi + collocA_air_rowi);
                rhs_per_vertex[pi][pvi] = (collocB_solid_rowi + collocB_air_rowi) * BC_dpdn[i];
                A[pi](pvi, pvi) += -omega;

            } else {
                // vertex i is on triple junction
                if (TPCF) {
                    // TPCF: the triple junction has a positional constraint, which makes p unknown but dpdn on the air side is known (whatever's necessary to cancel the triple junction's lateral motion)
                    A[pi].col(pvi) = (collocA_solid_rowi + collocA_air_rowi);
                    rhs_per_vertex[pi][pvi] = collocB_solid_rowi * BC_dpdn[i] + collocB_air_rowi * BC_dpdna[i];
                    A[pi](pvi, pvi) += -omega;
                } else {
                    //  dpdn on the air side is the unknown; dpdn on the solid side and p are both known
                    A[pi].col(pvi) = -collocB_air_rowi;
                    rhs_per_vertex[pi][pvi] = -(collocA_solid_rowi + collocA_air_rowi) * BC_p[i] + collocB_solid_rowi * BC_dpdn[i];
                    rhs[pi][pvi] += omega * BC_p[i];
                }
            }
        }
        //#pragma omp for

    }
    //#pragma omp parallel
    
    for (size_t pi = 0; pi < p2v.size(); pi++)
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
            rhs[pi] += rhs_per_vertex[pi][pvi];
    rhs_per_vertex.clear();

    // solve for the unknowns and interpret the result
    vector<vecX> solution(p2v.size());
    for (size_t pi = 0; pi < p2v.size(); pi++) {
        //        solution[pi] = A[pi].partialPivLu().solve(rhs[pi]);
        Eigen::BiCGSTAB<matX> s(A[pi]);
        solution[pi] = s.solve(rhs[pi]);
    }
    assert(solution[0] == solution[0]);

    vecX p = BC_p;
    vecX dpdn = BC_dpdn;
    for (size_t pi = 0; pi < p2v.size(); pi++) {
        for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++) {
            size_t i = p2v[pi][pvi];   // dof on vertex i, where the unknown is either p or dpdn depending on where vertex i is sitting

            if (!vertexIsSolid(i)) {
                // vertex i is in the interior of the liquid-air interface
                //  dpdn is the unknown; p is known
                dpdn[i] = solution[pi][pvi];
            } else if (vertexIsSolid(i) && !vertexIsOnTripleJunction(i)) {
                // vertex i is in the interior of the liquid-solid interface
                //  p is the unknown; dpdn is known
                p[i] = solution[pi][pvi]; 
            } else {
                // vertex i is on triple junction
                if (TPCF) {
                    // TPCF: for triple junction the dpdn on the air side is in BC_dpdna, not in the solution; the solution for this vertex is going to be a pressure, as it is no longer known
                    dpdn[i] = BC_dpdna[i];
                    p[i] = solution[pi][pvi];
                } else {
                    //  dpdn on the air side is the unknown; dpdn on the solid side and p are both known
                    dpdn[i] = solution[pi][pvi];  // this is for the air side; the solid side is stored in BC_dpdn[i]
                }
            }
        }
    }
    assert(p == p);
    assert(dpdn == dpdn);

    // compute the velocity update
    vecX dv = vecX::Zero(nv() * 3);
    for (size_t i = 0; i < nv(); i++) {
        if (!vertexIsOnTripleJunction(i)) {
            // for smooth regions (away from triple junctions) the averaged normal should work (at least in the limit of mesh refinement)
            vec3 n_i = vertOutwardNormal(i);

            // the normal pressure derivative comes from either the Neumann BC or the solve
            vec3 dpdn_i = dpdn[i] * n_i;

            // the tangential pressure derivatives are computed by finite differencing the p field (averaging the piecewise-constant gradient on each incident face)
            vec3 dpdt_i = vec3::Zero();
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++) {
                LosTopos::Vec3st t = mesh().m_tris[mesh().m_vertex_to_triangle_map[i][j]];
                vec3 x0 = pos(t[0]);
                vec3 x1 = pos(t[1]);
                vec3 x2 = pos(t[2]);
                vec3 n = (x1 - x0).cross(x2 - x0).normalized();
                vec3 gradp_integral = (p[t[0]] * n.cross(x2 - x1) + p[t[1]] * n.cross(x0 - x2) + p[t[2]] * n.cross(x1 - x0)) / 2;  // grad p times face area

                dpdt_i += gradp_integral / 3;   // a third of the face integral goes to vertex i
            }
            dpdt_i /= vertexArea(i);

            dv.segment<3>(i * 3) = -(dpdn_i + dpdt_i) / options.rho; 
        }
        else
        {
            // the triple junction is not smooth
            // use both BC_dpdn[i] (which corresponds to the solid side) and the dpdn value in solution[i] (which corresponds to the air side) to
            //  determine a gradient in the plane normal to the triple junction.
            // note that the tangential derivative of p also contains information about this gradient. not taking into account this information is dangerous
            //  because n_solid and n_air could be near colinear (e.g. in the case of superhydrophobia surfaces) and the solve using those two alone can be
            //  ill conditioned. solving a normal equation that uses both dpdn on both sides as well as dpdt instead is a better approach.
            vec3 n_solid = vec3::Zero();
            vec3 n_air = vec3::Zero();
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++)
            {
                size_t fj = mesh().m_vertex_to_triangle_map[i][j];
                if (faceIsSolid(fj))  n_solid += faceInteriorAngle(fj, i) * faceOutwardNormal(fj);
                else                    n_air += faceInteriorAngle(fj, i) * faceOutwardNormal(fj);
            }
            n_solid.normalize();
            n_air.normalize();

            vec3 triple_junction_tangent = n_solid.cross(n_air).normalized();

            vec3 t_solid = n_solid.cross(triple_junction_tangent).normalized();
            vec3 t_air = n_air.cross(triple_junction_tangent).normalized();

            // compute the tangential pressure derivatives by finite differencing the p field, for the air and solid sides separately
            vec3 dpdt_i_solid = vec3::Zero();
            vec3 dpdt_i_air = vec3::Zero();
            double vert_area_solid = 0;
            double vert_area_air = 0;
            for (size_t j = 0; j < mesh().m_vertex_to_triangle_map[i].size(); j++) {
                size_t fj = mesh().m_vertex_to_triangle_map[i][j];
                LosTopos::Vec3st t = mesh().m_tris[fj];
                vec3 x0 = pos(t[0]);
                vec3 x1 = pos(t[1]);
                vec3 x2 = pos(t[2]);
                vec3 n = (x1 - x0).cross(x2 - x0).normalized();
                vec3 gradp_integral = (p[t[0]] * n.cross(x2 - x1) + p[t[1]] * n.cross(x0 - x2) + p[t[2]] * n.cross(x1 - x0)) / 2;  // grad p times face area
                double fa = faceArea(fj);

                if (faceIsSolid(fj)) {
                    dpdt_i_solid += gradp_integral / 3;   // a third of the face integral goes to vertex i
                    vert_area_solid += fa;
                } else {
                    dpdt_i_air += gradp_integral / 3;   // a third of the face integral goes to vertex i
                    vert_area_air += fa;
                }
            }
            double dpdt_triple_junction = (dpdt_i_solid + dpdt_i_air).dot(triple_junction_tangent) / vertexArea(i);
            dpdt_i_solid /= vert_area_solid;
            dpdt_i_air /= vert_area_air;


            // solve a normal equation that includes both dpdn and dpdt information
            matX gradp_proj = matX::Zero(5, 3);
            vecX gradp_proj_rhs = vecX::Zero(5);
            gradp_proj.row(0) = n_solid.transpose();        gradp_proj_rhs[0] = BC_dpdn[i];                 // dpdn on the solid side
            gradp_proj.row(1) = n_air.transpose();          gradp_proj_rhs[1] = dpdn[i];                    // dpdn on the air side
            gradp_proj.row(2) = t_solid.transpose();        gradp_proj_rhs[2] = dpdt_i_solid.dot(t_solid);  // dpdt on the solid side
            gradp_proj.row(3) = t_air.transpose();          gradp_proj_rhs[3] = dpdt_i_air.dot(t_air);      // dpdn on the air side
            gradp_proj.row(4) = triple_junction_tangent;    gradp_proj_rhs[4] = dpdt_triple_junction;       // dpdt along the triple junction

            vec3 gradp = (gradp_proj.transpose() * gradp_proj).partialPivLu().solve(gradp_proj.transpose() * gradp_proj_rhs);

            // project the pressure gradient to remove the solid normal component, because that component should be the prescribed BC_dpdn[i]
            gradp = (mat3::Identity() - n_solid * n_solid.transpose()) * gradp + n_solid * BC_dpdn[i];

            dv.segment<3>(i * 3) = -gradp / options.rho;
        }
    }
    assert(dv == dv);

    // TPCF: set the tangential velocity of all solid vertices to the negative of current velocity. this is equivalent to a solid tangential force applied explicitly.
    if (TPCF) {
        for (size_t i = 0; i < nf(); i++) {
            if (faceIsSolid(i)) {
                vec3 nsolid = faceOutwardNormal(i);
                LosTopos::Vec3st t = mesh().m_tris[i];
                for (int j = 0; j < 3; j++) {
                    mat3 P = mat3::Identity() - nsolid * nsolid.transpose();
                    dv.segment<3>(t[j] * 3) = -P * v.segment<3>(t[j] * 3) + dv.segment<3>(t[j] * 3).dot(nsolid) * nsolid;    // set tangential velocity to cancel the original velocity's tangential component, without changing the normal component
                }
            }
        }
    }

    // TDMCV: manual (pseudo-)momentum conservation for tiny free droplets where the numerics can't hope to be accurate.
    if (TDMC) {
        for (size_t pi = 0; pi < p2v.size(); pi++) {
            double meanedgelen = 0;
            int pine = 0;
            for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++) {
                for (size_t j = 0; j < mesh().m_vertex_to_edge_map[p2v[pi][pvi]].size(); j++)
                    meanedgelen += edgeLength(mesh().m_vertex_to_edge_map[p2v[pi][pvi]][j]);
                pine += mesh().m_vertex_to_edge_map[p2v[pi][pvi]].size();
            }
            meanedgelen /= pine;

            double pivolume = partitionVolume(p2f[pi]);
            if (pivolume > pow(meanedgelen * 2, 3) * M_PI * 4 / 3)  // this is not a tiny droplet
                continue;

            std::cout << "TDMCV: processing tiny droplet (" << p2v[pi].size() << " vertices, " << p2f[pi].size() << " faces)." << std::endl;

            bool solid = false;
            for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
                if (vertexIsSolid(p2v[pi][pvi]))
                    solid = true;
            if (solid)      // this is not a free droplet
                continue;

            // remove the global rigid translation mode from dv for this droplet
            vec3 dvmean = vec3(0, 0, 0);
            for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
                dvmean += dv.segment<3>(p2v[pi][pvi] * 3);
            dvmean /= p2v[pi].size();
            for (size_t pvi = 0; pvi < p2v[pi].size(); pvi++)
                dv.segment<3>(p2v[pi][pvi] * 3) -= dvmean;
        }
    }
    assert(dv == dv);


    // update the velocity field
    vecX newv = v + dv;
    assert(v == v);
    assert(newv == newv);

    v = newv;

    for (size_t i = 0; i < nv(); i++)
        vel(i) = v.segment<3>(i * 3);
}



//**********************************************************
// Solve final velocity
//**********************************************************
void TPFluid::updateFinalVelocity(double dt, vecX& v) {
    for (size_t i = 0; i < nv(); i++)
        (*vels)[i] = v.segment<3>(i * 3);
}
