#define SAVE_IMAGES

#include <cassert>
#include <iostream>
using std::cout;
#include <string>
using std::string;

#include <shader.h>
#include <camera.h>
#include <fluid.h>
#include <shaders/shaders.h>
#include <scene.h>
#include <fluidOptions.h>

#include <subdivisionscheme.h>

#include <GLFW/glfw3.h>
#include <stb_image_write.h>

// Simulation state
struct SimState {
	string title;
	int w, h;
	GLFWwindow* window;
	
	int frame;
	int mode;
	bool m_ldown, m_rdown;
	double m_lx, m_ly, m_rx, m_ry;
	
	bool step, run;

	Shader hemi;
	ArcCamera camera;
	TPFluid* fluid;
} simState;

// Render scene
void render() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glm::mat4 projection = glm::perspective(glm::radians(simState.camera.zoom), 1.0f, 0.1f, 100.0f);
	glm::mat4 view = simState.camera.getViewMatrix();
	glm::mat4 model = glm::mat4(1.0f);

	if (simState.mode == 0) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	
	simState.hemi.use();
	simState.hemi.setMat4("projection", projection);
	simState.hemi.setMat4("view", view);
	simState.hemi.setMat4("model", model);

	simState.fluid->render();

	glFlush();
}

// Error callback
void error(int error, const char* description) {
	cout << "GLFW Error #" << error << ": " << description << std::endl;
}

// Key pressed callback
void keyPressed(GLFWwindow* window, int k, int sc, int action, int mods) {
	if (action == GLFW_PRESS) {
		if (k == GLFW_KEY_Q || k == GLFW_KEY_ESCAPE) {
			glfwSetWindowShouldClose(window, 1);
		} else if (k == GLFW_KEY_S) {
			simState.step = true;
		} else if (k == GLFW_KEY_R) {
			simState.run = !simState.run;
		} else if (k >= GLFW_KEY_0 && k <= GLFW_KEY_9) {
			simState.mode = k - GLFW_KEY_0;
		}
	}
}

// Mouse button callback
void mouseClicked(GLFWwindow* window, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		simState.m_ldown = true;
		glfwGetCursorPos(window, &simState.m_lx, &simState.m_ly);
	} else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
		simState.m_ldown = false;
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
		simState.m_rdown = true;
		glfwGetCursorPos(window, &simState.m_rx, &simState.m_ry);
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE) {
		simState.m_rdown = false;
    }
}

// Mouse motion callback
void mouseMoved(GLFWwindow* window, double x, double y) {
	if (simState.m_ldown) {
		simState.camera.updateMouse(x - simState.m_lx, y - simState.m_ly);
		simState.m_lx = x;
		simState.m_ly = y;
	}
}

// Scroll motion callback
void scrollWheel(GLFWwindow* window, double x, double y) {
	simState.camera.updateScroll(y);
}

// Save frame to PNG
void saveToPNG(string fname, GLFWwindow* window) {
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	
	GLsizei nrChannels = 3;
	GLsizei stride = nrChannels * width;
	stride += (stride % 4) ? (4 - stride % 4) : 0;
	GLsizei bufferSize = stride * height;
	std::vector<char> buffer(bufferSize);
	
	glPixelStorei(GL_PACK_ALIGNMENT, 4);
	glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());
	
	stbi_flip_vertically_on_write(true);
	stbi_write_png(fname.c_str(), width, height, nrChannels, buffer.data(), stride);
}


// Setup a GLFW window
GLFWwindow* createWindow(int w, int h, string title) {
	if (!glfwInit())
		assert(!"GLFW failed to initialize");
	
	glfwSetErrorCallback(error);

	GLFWwindow* window = glfwCreateWindow(w, h, title.c_str(), NULL, NULL);
	if (!window) {
		glfwTerminate();
		assert(!"GLFW failed to create a window");
	}

	glfwSetWindowPos(window, 100, 100);

	glfwSetKeyCallback(window, keyPressed);
	glfwSetMouseButtonCallback(window, mouseClicked);
	glfwSetCursorPosCallback(window, mouseMoved);
	glfwSetScrollCallback(window, scrollWheel);
	
	glfwMakeContextCurrent(window);
	
	if (glewInit() != GLEW_OK)
		assert(!"GLEW failed to initialize");
	
	return window;
}

// Initialize simulation
void initSim() {
	simState.title = "TPFluid";
	simState.w = 1500;
	simState.h = 1500;
	simState.window = createWindow(simState.w, simState.h, simState.title);
	
	simState.mode = 0;
	simState.camera = ArcCamera(glm::vec3(0, 0, 0), 4, 0, 45);
	simState.camera.mouseSens = 0.2;
	simState.camera.minZ = 1;
	simState.camera.maxZ = 20;
	simState.hemi = Shader(vertex_vs, hemi_fs);

	simState.m_ldown = false;
	simState.m_rdown = false;

	simState.step = false;
	simState.run = false;

	simState.frame = 1;

	glClearColor(1, 1, 1, 1);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_DEPTH_TEST);
}

// Cleanup simulation
void cleanupSim() {
	glfwTerminate();
}

// Render loop
int main(int argc, char* argv[]) {
	initSim();

	vector<LosTopos::Vec3d> vertices;
	vector<LosTopos::Vec3st> faces;
	vector<LosTopos::Vec2i> face_labels;
	vector<vec3> vels;
	vector<char> solids;
	struct FluidOptions options;
	LosTopos::SurfTrack* tmp = s_collision(vertices, faces, face_labels, vels, solids, options);
	simState.fluid = new TPFluid(tmp, vels, options);

	while (!glfwWindowShouldClose(simState.window)) {
		render();

		if (simState.step || simState.run) {
			simState.fluid->step();
			simState.step = false;
		}

		glfwSwapBuffers(simState.window);

		glfwPollEvents();
#ifdef SAVE_IMAGES
		if (simState.run)
			saveToPNG("frame" + std::to_string(simState.frame) + ".png", simState.window);
#endif
		simState.frame += 1;
	}

	cleanupSim();
}