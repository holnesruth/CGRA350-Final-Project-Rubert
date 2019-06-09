
#pragma once

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include "opengl.hpp"
#include "cgra/cgra_mesh.hpp"
#include "skeleton_model.hpp"


// Basic model that holds the shader, mesh and transform for drawing.
// Can be copied and modified for adding in extra information for drawing
// including textures for texture mapping etc.
struct basic_model {
	GLuint shader = 0;
	cgra::gl_mesh mesh;
	glm::vec3 color{0.7};
	glm::mat4 modelTransform{1.0};
	GLuint gradient;
	GLuint cubeMap;
	glm::vec2 tParams;
	glm::vec2 leParams;

	void draw(const glm::mat4 &view, const glm::mat4 proj);
	void updateParams(float minThickness, float maxThickness, float intensity, float transparency);
};


// Main application class
//
class Application {
private:
	// window
	glm::vec2 m_windowsize;
	GLFWwindow *m_window;

	bool m_view = false, m_shade = false, m_sim = false;

	// oribital camera
	float m_pitch = 0;
	float m_yaw = 0;
	float m_distance = 5;

	// last input
	bool m_leftMouseDown = false;
	glm::vec2 m_mousePosition;

	// drawing flags
	bool m_show_axis = false;
	bool m_show_grid = false;
	bool m_showWireframe = false;

	int m_min = 100;
	int m_max = 700;
	float m_intensity = 0.5;
	float m_opacity = 0.2;

	char* m_map_options[10] = { "Colosseum", "Creek", "LancellottiChapel",
								 "Lycksele", "MountainPath", "NissiBeach",
								 "PereaBeach", "SaintPetersBasilica", "Skansen",
								 "Tantolunden" };
	int m_map = 8;

	// geometry
	basic_model m_model;


public:
	// setup
	Application(GLFWwindow *);

	// disable copy constructors (for safety)
	Application(const Application&) = delete;
	Application& operator=(const Application&) = delete;

	// rendering callbacks (every frame)
	void render();
	void renderGUI();

	void showViewOptions();
	void showShaderOptions();
	void showSoftBodyOptions();

	// input callbacks
	void cursorPosCallback(double xpos, double ypos);
	void mouseButtonCallback(int button, int action, int mods);
	void scrollCallback(double xoffset, double yoffset);
	void keyCallback(int key, int scancode, int action, int mods);
	void charCallback(unsigned int c);

	unsigned int loadCubemap(std::vector<std::string> cubeFaces);
	void setUpCubeMap(char* name);
};