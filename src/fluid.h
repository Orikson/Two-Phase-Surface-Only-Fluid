#ifndef FLUID_H
#define FLUID_H

#include <vector>
using std::vector;

#include <collisionpipeline.h>
#include <surftrack.h>
#include <dynamicsurface.h>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <GL/glew.h>

#include <eigenHeaders.h>

#include <fluidOptions.h>
#include <vao.h>
#include <integrals/BoundaryIntegral.h>

// Surface partition
struct Partitioning {
	vector<vector<size_t>> p2v;		// a list of vertices for each partition
	vector<vector<size_t>> p2f;		// a list of faces for each partition
	vector<int> v2p;				// a partition index for each vertex
	vector<int> f2p;				// a partition index for each face

	vector<size_t> flattened_partition_vertices;	// flatten out all partitions of vertices into one linear array
	vector<size_t> indices_in_partitions;			// for each vertex in flattened_partition_vertices, record its index inside its respective partition
};

class TPFluid : public LosTopos::SurfTrack::SolidVerticesCallback, public LosTopos::T1Transition::VelocityFieldCallback, public LosTopos::SurfTrack::MeshEventCallback {
	public:
		TPFluid() {}
		TPFluid(LosTopos::SurfTrack* surface, vector<vec3>& vels, struct FluidOptions& options);
		~TPFluid();

		void render();
		void step();

	private:
		LosTopos::SurfTrack* surf;
		LosTopos::NonDestructiveTriMesh::VertexData<vec3>* vels;
		vector<vec3> vn_tmp;
		vecX m_pre_stepping_geometry;
		struct FluidOptions options;

		int frame = 0;
		double simTime;

		// Primary timestep operations (in order)
		double advect(double dt);						// returns the final dt accepted by LosTopos
		void improveMesh();								// ensures that the mesh is properly defragged
		void fullPartition(Partitioning& partitioning);	// partition mesh into a set of (naive) contractible domains, then remove internal bubbles
		void helmholtzDecomposition(double dt, vecX& v, Partitioning& partitioning);
		void smoothVelocity(double dt, vecX& v, Partitioning& partitioning, double coef);
		void concaveSmoothing(double dt, vecX& v, vecX& oldv);		// might not be needed? (CSOC is false)
		void addGravity(double dt, vecX& v);
		void solvePressure(double dt, vecX& v, Partitioning& partitioning);
		void updateFinalVelocity(double dt, vecX& v);

		// Useful operations
		void computeNormals();
		void partitionMesh(Partitioning& partitioning);
		double partitionVolume(const vector<size_t>& faces) const;
		double vertIntegralMeanCurvature(size_t v) const;
		double currentInflux() const { return options.influx + options.influx_inc * simTime; }

		const LosTopos::SurfTrack* surfTrack() const { return surf; }
			  LosTopos::SurfTrack* surfTrack()       { return surf; }
		const LosTopos::NonDestructiveTriMesh& mesh() const { return surf->m_mesh; }
			  LosTopos::NonDestructiveTriMesh& mesh()       { return surf->m_mesh; }
		
		size_t nv() const { return mesh().nv(); }
		size_t nf() const { return mesh().nt(); }

		const vec3& vel(size_t i) const { return (*vels)[i]; }
			  vec3& vel(size_t i)		{ return (*vels)[i]; }
		vec3 pos(size_t i) const { return vc(surf->pm_positions[i]); }

		size_t edgeOtherVertex(size_t e, size_t v) const { LosTopos::Vec2st edge = surf->m_mesh.m_edges[e]; return edge[0] == v ? edge[1] : edge[0]; }
		
		double faceArea(size_t i) const { return surf->get_triangle_area(i); }
		double vertexArea(size_t i) const { double a = 0; for (auto f : surf->m_mesh.m_vertex_to_triangle_map[i]) a += faceArea(f); return a / 3; }
		double edgeArea(size_t i) const { double a = 0; for (auto f : surf->m_mesh.m_edge_to_triangle_map[i]) a += faceArea(f); return a / 3; }
		double edgeLength(size_t e) const { return surf->get_edge_length(e); }
		double angleAroundAxis(const vec3& v0, const vec3& v1, const vec3& a) const { double asq = a.squaredNorm(); assert(asq != 0); vec3 u = v0 - v0.dot(a) * a / asq; assert(u.squaredNorm() != 0); u.normalize(); vec3 v = a.cross(u).normalized(); return atan2(v1.dot(v), v1.dot(u)); }
		double faceInteriorAngle(size_t f, size_t v) const { assert(surf->m_mesh.triangle_contains_vertex(surf->m_mesh.m_tris[f], v)); LosTopos::Vec3st t = getShuffledTriangle(surf->m_mesh.m_tris[f], v); double as = (pos(t[1]) - pos(t[0])).squaredNorm(); double bs = (pos(t[2]) - pos(t[0])).squaredNorm(); double cs = (pos(t[2]) - pos(t[1])).squaredNorm(); double a = sqrt(as); double b = sqrt(bs); double cosa = (as + bs - cs) / (2 * a * b); if (cosa < -1) cosa = -1; if (cosa > 1)  cosa = 1; return acos(cosa); }
		double vertInteriorSolidAngle(size_t v) const;
		double tripleJunctionVirtualWidth() const { return options.mean_edge_length * 0.25; }

		vec3 faceOutwardNormal(size_t f) const { LosTopos::Vec3st t = surf->m_mesh.m_tris[f]; vec3 n = (pos(t[1]) - pos(t[0])).cross(pos(t[2]) - pos(t[0])).normalized(); LosTopos::Vec2i l = surf->m_mesh.get_triangle_label(f); if (l[0] < l[1]) return -n; else return n; }
		vec3 vertOutwardNormal(size_t v) const { vec3 n = vec3::Zero(); for (size_t i = 0; i < surf->m_mesh.m_vertex_to_triangle_map[v].size(); i ++) n += faceArea(surf->m_mesh.m_vertex_to_triangle_map[v][i]) * faceOutwardNormal(surf->m_mesh.m_vertex_to_triangle_map[v][i]); return n.normalized(); }
		vec3 edgeTangent(size_t e)  const { return (pos(surf->m_mesh.m_edges[e][1]) - pos(surf->m_mesh.m_edges[e][0])).normalized(); }
		
		bool vertexIsSolid(size_t i) const { return surf->vertex_is_any_solid(i); }
		bool faceIsSolid(size_t i) const { LosTopos::Vec3st t = surf->m_mesh.m_tris[i]; return vertexIsSolid(t[0]) && vertexIsSolid(t[1]) && vertexIsSolid(t[2]); }
		bool vertexIsOnTripleJunction(size_t v) const { if (!vertexIsSolid(v)) return false; for (size_t i = 0; i < surf->m_mesh.m_vertex_to_triangle_map[v].size(); i++) if (!faceIsSolid(surf->m_mesh.m_vertex_to_triangle_map[v][i])) return true; return false; }
		bool edgeIsOnTripleJunction(size_t e) const { bool solid = false; bool air = false; for (size_t i = 0; i < surf->m_mesh.m_edge_to_triangle_map[e].size(); i++) if (faceIsSolid(surf->m_mesh.m_edge_to_triangle_map[e][i])) solid = true; else air = true; return (solid && air); }

		LosTopos::Vec3st getShuffledTriangle(const LosTopos::Vec3st& t, size_t vertex_to_be_front) const;
		mat3 getVertexPositions(const LosTopos::Vec3st& t) const;
		mat3 getVertexVelocities(const LosTopos::Vec3st& t) const;

		// Render functions
		void renderBBox();

	public:
		// SurfTrack::SolidVerticesCallback method
		bool            generate_collapsed_position(LosTopos::SurfTrack& st, size_t v0, size_t v1, LosTopos::Vec3d& pos);
		bool            generate_split_position(LosTopos::SurfTrack& st, size_t v0, size_t v1, LosTopos::Vec3d& pos);
		bool            generate_snapped_position(LosTopos::SurfTrack& st, size_t v0, size_t v1, LosTopos::Vec3d& pos);
		LosTopos::Vec3c generate_collapsed_solid_label(LosTopos::SurfTrack& st, size_t v0, size_t v1, const LosTopos::Vec3c& label0, const LosTopos::Vec3c& label1);
		LosTopos::Vec3c generate_split_solid_label(LosTopos::SurfTrack& st, size_t v0, size_t v1, const LosTopos::Vec3c& label0, const LosTopos::Vec3c& label1);
		LosTopos::Vec3c generate_snapped_solid_label(LosTopos::SurfTrack& st, size_t v0, size_t v1, const LosTopos::Vec3c& label0, const LosTopos::Vec3c& label1);
		bool            generate_edge_popped_positions(LosTopos::SurfTrack& st, size_t oldv, const LosTopos::Vec2i& cut, LosTopos::Vec3d& pos_upper, LosTopos::Vec3d& pos_lower);
		bool            generate_vertex_popped_positions(LosTopos::SurfTrack& st, size_t oldv, int A, int B, LosTopos::Vec3d& pos_a, LosTopos::Vec3d& pos_b);
		bool            solid_edge_is_feature(const LosTopos::SurfTrack& st, size_t e);

		// T1Transition::VelocityFieldCallback methods
		LosTopos::Vec3d sampleVelocity(LosTopos::Vec3d& pos);
		bool sampleDirectionalDivergence(const LosTopos::Vec3d& pos, const LosTopos::Vec3d& dir, double& output);

		// SurfTrack::MeshEventCallback
		void pre_collapse(const LosTopos::SurfTrack& st, size_t e, void** data);
		void post_collapse(const LosTopos::SurfTrack& st, size_t e, size_t merged_vertex, void* data);

		void pre_split(const LosTopos::SurfTrack& st, size_t e, void** data);
		void post_split(const LosTopos::SurfTrack& st, size_t e, size_t new_vertex, void* data);

		void pre_flip(const LosTopos::SurfTrack& st, size_t e, void** data);
		void post_flip(const LosTopos::SurfTrack& st, size_t e, void* data);

		void pre_t1(const LosTopos::SurfTrack& st, size_t v, void** data);
		void post_t1(const LosTopos::SurfTrack& st, size_t v, size_t a, size_t b, void* data);

		void pre_facesplit(const LosTopos::SurfTrack& st, size_t f, void** data);
		void post_facesplit(const LosTopos::SurfTrack& st, size_t f, size_t new_vertex, void* data);

		void pre_snap(const LosTopos::SurfTrack& st, size_t v0, size_t v1, void** data);
		void post_snap(const LosTopos::SurfTrack& st, size_t v_kept, size_t v_deleted, void* data);

		void pre_smoothing(const LosTopos::SurfTrack& st, void** data);
		void post_smoothing(const LosTopos::SurfTrack& st, void* data);

		std::ostream& log() { return std::cout; }
		std::ostream& callback_log() { return std::cout; }
};

#endif