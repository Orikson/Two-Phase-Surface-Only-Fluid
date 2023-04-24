set (Directories
 .
 util
 integrals
 shaders
)

append_files (Headers "h" ${Directories})
append_files (Sources "cpp" ${Directories})

find_package (OpenGL REQUIRED)
include_directories (${OPENGL_INCLUDE_DIR})

find_package (GLUT REQUIRED glut)
include_directories (${GLUT_INCLUDE_DIR})

find_package (glm CONFIG REQUIRED)

find_package (GLEW REQUIRED)
include_directories (${GLEW_INCLUDE_DIRS})
link_libraries (${GLEW_LIBRARIES})

add_executable (TPFluid ${Headers} ${Sources})
target_link_libraries (TPFluid LosTopos glm::glm glfw ${DEFAULT_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY} ${GLEW_LIBRARIES})
