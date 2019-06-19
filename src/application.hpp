#pragma once

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include "opengl.hpp"
#include "cgra/cgra_mesh.hpp"
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
	GLuint noise;
    GLuint cubeMap;

    glm::vec2 thicknessParams;
    glm::vec2 lightParams;
	glm::vec2 flowSpeeds = glm::vec2(0.6, 1.6);

	float time = 100.0;
	int flowOctaves = 7;

    void draw(const glm::mat4 &view, const glm::mat4 proj, bool drawAsSphere);
    void updateParams(float minThickness, float maxThickness, float intensity, float transparency);
};


// Main application class
//
class Application {
public:
	float m_distance = 20;
private:
	// window
	glm::vec2 m_windowsize;
	GLFWwindow *m_window;

	// gui window booleans for show/hide
	bool m_view = true, m_shade = true, m_sim = true;

	// oribital camera
	float m_pitch = .86;
	float m_yaw = -.86;

	// last input
	bool m_leftMouseDown = false;
	glm::vec2 m_mousePosition;

    // drawing flags
	bool m_show_axis = false;
	bool m_show_grid = false;
	bool m_showWireframe = false;

	// geometry
	basic_model m_model;

    /** =========================== Robert's Parameters ===================== */

	// timer
	double m_lastMillis;

	// soft bodies
	std::vector<Softbody> m_softbodies;

	// utility meshes
    cgra::gl_mesh m_bbox_mesh;
    cgra::gl_mesh m_ground_plane_mesh;

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

    /** =========================== Ruth's Parameters ===================== */

    // Thickness
    int m_min = 100;
    int m_max = 700;

	// Light
    float m_intensity = 0.5;

	// Surface appearance
    float m_opacity = 0.2;
	float m_speed = 2.0;

	// Cube maps
    char* m_map_options[10] = { "Colosseum", "Creek", "LancellottiChapel",
                                "Lycksele", "MountainPath", "NissiBeach",
                                "PereaBeach", "SaintPetersBasilica", "Skansen",
                                "Tantolunden" };
    int m_map = 8;

	// Flow noise
	bool m_flow = false;

    /** =========================== Shared Parameters ===================== */

	// Shaders
    GLuint m_shader_default;
    GLuint m_shader_bubble;

	// Modes
	enum m_mode { Shader, Simulation, FullDemo };
	char* m_mode_options[3] = { "Shader", "Simulation", "Full Demo" };
	int m_current_mode = 0;

public:
	// setup
	Application(GLFWwindow *);

	// disable copy constructors (for safety)
	Application(const Application&) = delete;
	Application& operator=(const Application&) = delete;

	// rendering callbacks (every frame)
	void render();
	void renderGUI();

	// GUI options
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

	/** =========================== Robert's Functions ===================== */

    static void cleanMesh(cgra::mesh_builder &mesh);

    void addNewSoftbody(glm::mat4 initialTransform, bool printVerts);

    void createBBox();

    void createGroundplane();

    /** =========================== Ruth's Functions ===================== */

	// Draws a model in two passes using openGL face culling
	void drawModel(glm::mat4& view, glm::mat4& proj);

	// Updates the flow timer
	void updateFlow();

    // Load a cubemap from a vector of image filenames
    unsigned int loadCubemap(std::vector<std::string> cubeFaces);

	// Construct the cubemap and update the shader variable
    void setUpCubeMap(char* mapName);

	// Calculates the position of the camera
	glm::vec3 cameraPos();

	// Sorts the softbody vector by their distance from the camera
	void sortSoftBodies();
};