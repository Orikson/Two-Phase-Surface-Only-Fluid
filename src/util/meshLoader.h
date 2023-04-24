#ifndef MESH_LOADER_H
#define MESH_LOADER_H

#include <vector>
using std::vector;
#include <string>
using std::string;
#include <algorithm>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include <surftrack.h>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <GL/glew.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include <GL/glut.h>

#include <eigenHeaders.h>

inline LosTopos::Vec3d vc(const aiVector3D& v) {
    return LosTopos::Vec3d(v.x, v.y, v.z);
}

// Mesh loading interface using Assimp
// Only loads vertices and face indices (not textures)
class MeshLoader {
    public:
        // Load a mesh as a solid or as a fluid
        bool loadModel(const string& path, vec3 translate, float scale, vector<LosTopos::Vec3d>& vertices, vector<LosTopos::Vec3st>& faces, vector<LosTopos::Vec2i>& labels, vector<char>& solids, const vec2i& label = vec2i(1, 0), bool isSolid = false) {
            Assimp::Importer importer;
            const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_GenSmoothNormals | aiProcess_CalcTangentSpace | aiProcess_JoinIdenticalVertices);

            if (!scene) {
                cout << "Failed to open mesh at " << path << std::endl;
                return false;
            }

            return init(scene, translate, scale, vertices, faces, labels, solids, label, isSolid);
        }

    private:
        // Initialize all sub meshes of a loaded model
        bool init(const aiScene* scene, vec3 translate, float scale, vector<LosTopos::Vec3d>& vertices, vector<LosTopos::Vec3st>& faces, vector<LosTopos::Vec2i>& labels, vector<char>& solids, const vec2i& label = vec2i(1, 0), bool isSolid = false) {
            for (unsigned int i = 0; i < scene->mNumMeshes; i ++) {
                const aiMesh* mesh = scene->mMeshes[i];

                for (size_t j = 0; j < mesh->mNumVertices; j ++) solids.push_back(isSolid);
                for (size_t j = 0; j < mesh->mNumFaces; j ++) labels.push_back(LosTopos::Vec2i(label[0], label[1]));

                initMesh(mesh, vertices, faces, translate, scale);
            }
            return true;
        }

        // Initialize a single sub mesh of a model
        bool initMesh(const aiMesh* mesh, vector<LosTopos::Vec3d>& vertices, vector<LosTopos::Vec3st>& faces, vec3 translate, float scale) {
            size_t lastIndex = vertices.size();
            
            for (unsigned int i = 0; i < mesh->mNumVertices; i ++) {
                LosTopos::Vec3d vert = vc(mesh->mVertices[i]) * scale + vc(translate);
                vertices.push_back(vert);
            }

            for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
                const aiFace& face = mesh->mFaces[i];
                assert(face.mNumIndices == 3);
                faces.push_back(LosTopos::Vec3st(face.mIndices[0] + lastIndex, face.mIndices[1] + lastIndex, face.mIndices[2] + lastIndex));
            }

            return true;
        }
};

#endif