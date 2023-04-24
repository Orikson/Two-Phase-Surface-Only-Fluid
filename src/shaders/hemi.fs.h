#ifndef HEMI_FS_H
#define HEMI_FS_H

#include <string>

const string hemi_fs = R"(
#version 430 core
out vec4 FragColor;

in vec3 Normal;  
in vec3 FragPos;

// hemispherical lighting
vec4 hemi(vec3 N, vec3 objColor, vec3 lightColor) {
    float up = N.z * 0.5 + 0.5;
    vec3 amb = vec3(0.1) + up * lightColor;

    return vec4(amb * objColor, 1.);
}

void main() {
    FragColor = hemi(Normal, vec3(0.6), vec3(1.));
}
)";

#endif