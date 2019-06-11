
#pragma once

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include "opengl.hpp"
#include "cgra/cgra_mesh.hpp"
#include "skeleton_model.hpp"
#include "Softbody.h"


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

    void draw(const glm::mat4 &view, const glm::mat4 proj, bool drawAsSphere);
    void updateParams(float minThickness, float maxThickness, float intensity, float transparency);
};


// Main application class
//
class Application {
private:
	// window
	glm::vec2 m_windowsize;
	GLFWwindow *m_window;

	bool m_view = true, m_shade = true, m_sim = true;

	// oribital camera
	float m_pitch = .86;
	float m_yaw = -.86;
	float m_distance = 20;

	// last input
	bool m_leftMouseDown = false;
	glm::vec2 m_mousePosition;

    // drawing flags
	bool m_show_axis = false;
	bool m_show_grid = true;
	bool m_showWireframe = false;

	// geometry
	basic_model m_model;

    /** =========================== Robert Parameters ===================== */
	// timer
	double m_lastMillis;

	// soft bodies
	std::vector<Softbody> m_softbodies;

	// utility meshes
    cgra::gl_mesh m_bbox_mesh;
    cgra::gl_mesh m_ground_plane_mesh;
    cgra::gl_mesh m_hiRez_ball;


    // soft body geometry
	float m_ball_radius = 4;
	cgra::mesh_builder m_mesh;

    // simulation parameters
    float m_gravity = 1;
    float m_mass = 1.0;
    float m_kd = 6;
    float m_ks = 6;
    float m_pressure = 100;

    // mouse mode
    bool m_place_softbodies = false;

    /** =========================== Ruth Parameters ===================== */
    // Bubbles
    int m_min = 100;
    int m_max = 700;
    float m_intensity = 0.5;
    float m_opacity = 0.2;

    char* m_map_options[10] = { "Colosseum", "Creek", "LancellottiChapel",
                                "Lycksele", "MountainPath", "NissiBeach",
                                "PereaBeach", "SaintPetersBasilica", "Skansen",
                                "Tantolunden" };
    int m_map = 8;

    /** =========================== Shared Parameters ===================== */

    GLuint m_shader_default;
    GLuint m_shader_bubble;

	enum m_mode { Shader, Simulation, FullDemo };
	char* m_mode_options[3] = { "Shader", "Simulation", "Full Demo" };
	int m_current_mode = Simulation;

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
	void showModeChanger();

	// input callbacks
	void cursorPosCallback(double xpos, double ypos);
	void mouseButtonCallback(int button, int action, int mods);
	void scrollCallback(double xoffset, double yoffset);
	void keyCallback(int key, int scancode, int action, int mods);
	void charCallback(unsigned int c);

	/** =========================== Robert Functions ===================== */

    static void cleanMesh(cgra::mesh_builder &mesh);
    void addNewSoftbody(glm::mat4 initialTransform, bool printVerts);
    void createBBox();

    /** =========================== Ruth Functions ===================== */
    // Cube mapping
    unsigned int loadCubemap(std::vector<std::string> cubeFaces);
    void setUpCubeMap(char* name);

    void createGroundplane();

    void drawModel(glm::mat4 &view, glm::mat4 &proj);
};