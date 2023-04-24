#ifndef VAO_H
#define VAO_H

#include "shader.h"

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <GL/glew.h>

#include <iostream>
#include <functional>
#include <utility>
#include <fstream>
#include <sstream>
#include <string>
using std::string;

#include <vector>
using std::vector;

#include <eigenHeaders.h>

// For interleaving different vertex attributes
struct RenderVertex {
    vec3f position;
    vec3f normal;
};

/**
 * @brief Reduced overhead for (dynamic?) vertex array objects
 */
class VertexArrayObject {
    public:
        VertexArrayObject() {
            VAO = 0;
            VBO = 0;
            EBO = 0;
        }
        ~VertexArrayObject() {
            remove();
        }

        void remove() {
            glDeleteVertexArrays(1, &VAO);
            glDeleteBuffers(1, &VBO);
            glDeleteBuffers(1, &EBO);

            VAO = 0; VBO = 0; EBO = 0;
        }

        void init(vector<RenderVertex>& vertices, vector<vec3ui>& indices) {
            vsize = vertices.size(); verts = vsize;
            isize = indices.size() * 3; elems = isize * 3;

            glGenVertexArrays(1, &VAO);
            glBindVertexArray(VAO);

            glGenBuffers(1, &VBO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER,
                vertices.size() * sizeof(RenderVertex),     // size of vertices buffer
                &vertices[0],                               // pointer to first element
                GL_STATIC_DRAW);

            glGenBuffers(1, &EBO);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                indices.size() * sizeof(vec3ui),            // size of indices buffer
                &indices[0],                                // pointer to first element
                GL_STATIC_DRAW);

            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);

            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
        }

        void bind() {
            glBindVertexArray(VAO);
        }

        void draw() {
            bind();
            glDrawElements(GL_TRIANGLES, elems, GL_UNSIGNED_INT, (void*)0);
        }

        void update(vector<float>& vertices, vector<unsigned int>& indices) {
            bind();

            if (vertices.size() > vsize) {
                vsize = vertices.size();
                verts = vertices.size();

                glDeleteBuffers(1, &VBO);
                
                glGenBuffers(1, &VBO);
                glBindBuffer(GL_ARRAY_BUFFER, VBO);
                glBufferData(GL_ARRAY_BUFFER,
                    vertices.size() * sizeof(float),            // size of vertices buffer
                    &(vertices)[0],                             // pointer to first element
                    GL_STATIC_DRAW);
            }
            else {
                verts = vertices.size();
                glBindBuffer(GL_ARRAY_BUFFER, VBO);
                glBufferSubData(GL_ARRAY_BUFFER, 0, vertices.size() * sizeof(float), &vertices[0]);
            }

            if (indices.size() > isize) {
                isize = indices.size();
                elems = indices.size();

                glDeleteBuffers(1, &EBO);

                glGenBuffers(1, &EBO);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                    indices.size() * sizeof(unsigned int),      // size of indices buffer
                    &(indices)[0],                              // pointer to first element
                    GL_STATIC_DRAW);
            }
            else {
                elems = indices.size();
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
                glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, indices.size() * sizeof(unsigned int), &indices[0]);
            }
        }

    private:
        size_t vsize, isize;
        size_t verts, elems;
        GLuint VAO, VBO, EBO;
};

#endif
