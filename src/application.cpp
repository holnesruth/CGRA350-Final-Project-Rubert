
// std
#include <iostream>
#include <string>
#include <chrono>
#include <algorithm>

// glm
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/random.hpp>

// project
#include "application.hpp"
#include "cgra/cgra_geometry.hpp"
#include "cgra/cgra_gui.hpp"
#include "cgra/cgra_image.hpp"
#include "cgra/cgra_shader.hpp"
#include "cgra/cgra_wavefront.hpp"


using namespace std;
using namespace cgra;
using namespace glm;


void basic_model::draw(const glm::mat4& view, const glm::mat4 proj, bool drawAsSphere) {
	mat4 modelview = view * modelTransform;

	// load shader and variables
	glUseProgram(shader);

	// matrices
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelMatrix"), 1, false, value_ptr(modelTransform));

	// colour
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));
	glUniform2fv(glGetUniformLocation(shader, "uThickness"), 1, value_ptr(tParams));
	glUniform2fv(glGetUniformLocation(shader, "uLightEffects"), 1, value_ptr(leParams));

	// flow noise
	glUniform1f(glGetUniformLocation(shader, "uTime"), time);
	glUniform2fv(glGetUniformLocation(shader, "uResolution"), 1, value_ptr(vec2(8.0, 8.0)));
	glUniform1i(glGetUniformLocation(shader, "uFlow"), flow);

	// textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, gradient);
	glUniform1i(glGetUniformLocation(shader, "uTexture"), 0);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, noise);
	glUniform1i(glGetUniformLocation(shader, "uNoise"), 1);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_CUBE_MAP, cubeMap);
	glUniform1i(glGetUniformLocation(shader, "uCubeMap"), 2);

	// decide which mesh to draw
	if (!drawAsSphere) {
		mesh.draw(); // draw the soft body mesh
	} else {
		drawSphere(); // draw a standard sphere
	}
}

// Update the models's parameters for shading
void basic_model::updateParams(float minThickness, float maxThickness, float intensity, float opacity) {
	tParams.x = minThickness;
	tParams.y = maxThickness;
	leParams.x = intensity;
	leParams.y = opacity;
}

// Load the application
Application::Application(GLFWwindow* window) : m_window(window) {

	// Create the shaders
	shader_builder sb;

	// Bubble shader
	sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//simple_vert.glsl"));
	sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//bubbles.glsl"));
	m_shader_bubble = sb.build();

	// Colour shader for simulation
	sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_vert.glsl"));
	sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_frag.glsl"));
	m_shader_default = sb.build();

	// ====================== Soft body ==========================
    m_mesh = load_wavefront_data(CGRA_SRCDIR + std::string("/res//assets//ball_480.obj"));
    cleanMesh(m_mesh);
    m_softbodies.emplace_back();

    mat4 initialPositionTransform = translate(mat4(1.0f), vec3(0, 2, 0)) * scale(mat4(1.0f), vec3(m_ball_radius));
    m_softbodies.at(0).initializeMesh(m_mesh, initialPositionTransform);

	for( auto &softbody : m_softbodies){
	    m_model.mesh = softbody.constructMesh(m_showWireframe);
	}

	createBBox();
	createGroundplane();
	m_lastMillis = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();

	m_hiRez_ball = load_wavefront_data(CGRA_SRCDIR + std::string("/res//assets//ball_1080.obj")).build();

    // ====================== Shaders ==========================
    m_model.shader = m_shader_default;
	m_model.color = vec3(1, 1, 1);

	// set up the thickness texture
	rgba_image img = rgba_image(CGRA_SRCDIR + std::string("//res//textures//grad.jpg"));
	m_model.gradient = img.uploadTexture();

	// Set up the noise texture
	img = rgba_image(CGRA_SRCDIR + std::string("//res//textures//noise.png"));
	m_model.noise = img.uploadTexture();

	// Create the reflections with a cube map
	setUpCubeMap("Skansen");
}


void Application::render() {

	// retrieve the window height
	int width, height;
	glfwGetFramebufferSize(m_window, &width, &height);

	m_windowsize = vec2(width, height); // update window size
	glViewport(0, 0, width, height); // set the viewport to draw to the entire window

	// clear the back-buffer
	glClearColor(0.3f, 0.3f, 0.4f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// enable flags for normal/forward rendering
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	// Enable tgransparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// projection matrix
	mat4 proj = perspective(1.f, float(width) / height, 0.1f, 1000.f);

	// view matrix
	mat4 view = translate(mat4(1), vec3(3, 0, -m_distance))
		* rotate(mat4(1), m_pitch, vec3(1, 0, 0))
		* rotate(mat4(1), m_yaw, vec3(0, 1, 0));

	// helpful draw options
	glPolygonMode(GL_FRONT_AND_BACK, (m_showWireframe) ? GL_LINE : GL_FILL);

	// update model parameters for shading
	m_model.updateParams(m_min, m_max, m_intensity, m_opacity);


    if (m_current_mode == Simulation || m_showWireframe) {
        m_model.shader = m_shader_default;
    }
    else {
        m_model.shader = m_shader_bubble;
    }


    if (m_show_grid) drawGrid(view, proj);
	if (m_show_axis) drawAxis(view, proj);

    double millis = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();



    // if in shader mode, dont bother with soft bodies
    if (m_current_mode == Shader) {
        m_model.mesh = m_hiRez_ball;
        m_model.modelTransform = scale(rotate(mat4(1.0f), radians(90.0f), vec3(1, 0, 0)), vec3(m_ball_radius));
        m_show_grid = false;
        drawModel(view, proj);

        // update flow noise parameter
        if (millis - m_lastMillis > DT * 1000) {
            m_model.time += m_speed / 100;
            if (m_model.time > 150) {
                m_model.time = 100; // reset
            }
            m_lastMillis = chrono::duration_cast<chrono::milliseconds>(
                    chrono::system_clock::now().time_since_epoch()).count();
        }
        return;
    }

    /** ============ Simulation Step ====================== */
    if (millis - m_lastMillis > DT * 1000) {

        m_model.time += m_speed / 100;
        if (m_model.time > 150) {
            m_model.time = 100; // reset
        }

        for (auto &softbody : m_softbodies)
            softbody.updateCentroid();

        for (auto &softbody : m_softbodies) {
            softbody.AccumulateForces();
            softbody.IntegrateForces(m_current_mode != FullDemo, m_softbodies, m_ball_radius);
        }
        m_lastMillis = chrono::duration_cast<chrono::milliseconds>(
                chrono::system_clock::now().time_since_epoch()).count();

    }

    /** ============ Draw Softbodies ====================== */
    for (auto &softbody : m_softbodies) {

        m_model.mesh = softbody.constructMesh(m_showWireframe);

        // if in wireframe draw spheres on all points
        if (m_showWireframe) {
            for (auto &point : softbody.m_points) {
                m_model.modelTransform = translate(mat4(1.0f), point.pos) * scale(mat4(1.0f), vec3(0.1f));
                m_model.draw(view, proj, true);
            }
        }

        m_model.modelTransform = mat4(1.0);

        drawModel(view, proj);

    }

    m_model.shader = m_shader_default;

    // draw bounding box
    m_model.mesh = m_bbox_mesh;
    drawModel(view, proj);

    // if in sim mode, draw the ground plane
    if (m_current_mode == Simulation) {
        m_model.mesh = m_ground_plane_mesh;
    }
    drawModel(view, proj);


}


/**  ========================================== Robert Functions ==========================================*/


void Application::drawModel(mat4 &view, mat4 &proj){
	glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);

    m_model.leParams.y *= m_intensity;
    m_model.draw(view, proj, m_current_mode == Shader);

    glCullFace(GL_BACK);

    m_model.leParams.y = m_opacity;
    m_model.draw(view, proj, m_current_mode == Shader);
}

void Application::cleanMesh(mesh_builder &mesh){
    /** ============== Clean up the mesh data ============== */
    // for every vert in the mesh
    for (int i = mesh.vertices.size() - 1; i >= 0; i--) {
        vec3 vertPos = mesh.vertices.at(i).pos;

        // go through the rest of the vertices
        for (int j = i - 1; j >= 0; --j) {
            // if a duplicate is found
            if (vertPos == mesh.vertices.at(j).pos) {
                // delete the duplicate
                mesh.vertices.erase(next(begin(mesh.vertices), i));

                // go through the list of verticies
                for (int z = i; z < mesh.indices.size(); z++) {
                    // if they were pointing at the deleted index
                    if (mesh.indices.at(z) == i) {
                        // change them to point at the found duplicate index
                        mesh.indices.at(z) = j;
                        // if they were pointing at a vertex after the deleted index
                    }
                    else if (mesh.indices.at(z) > i) {
                        // subtract one to make sure they still point at the same vertex as before
                        // (because we deleted one and the indices all shift over by one)
                        mesh.indices.at(z) = mesh.indices.at(z) - 1;
                    }
                }
                break;
            }
        }
    }
}

void Application::addNewSoftbody(glm::mat4 initialTransform, bool printVerts) {
    m_softbodies.emplace_back();
    m_softbodies.back().initializeMesh(m_mesh, initialTransform);
    if (m_softbodies.size() > 20){
        m_softbodies.erase(m_softbodies.begin());
    }
    m_softbodies.back().m_gravity = m_gravity;
    m_softbodies.back().m_mass = m_mass;
    m_softbodies.back().m_kd = m_kd;
    m_softbodies.back().m_ks = m_ks;
    m_softbodies.back().m_pressure = m_pressure;

    if (!printVerts) return;

    int totalPoints = 0;
    int totalSprings = 0;
    for (auto &softbody : m_softbodies) {
        totalPoints += softbody.m_points.size();
        totalSprings += softbody.m_springs.size();
    }

    cout << "New soft body added for a total of:" << endl;
    cout << "\t * " << totalPoints << " Simulating points"  << endl;
    cout << "\t * " << totalSprings << " Simulating springs"  << endl;

}

void Application::createBBox() {
    mesh_builder mb(GL_LINES);

    mb.push_vertex(mesh_vertex({vec3(bbox.x, bbox.y, bbox.z),
                                normalize(vec3(bbox.x, bbox.y, bbox.z)),
                                vec2(0)}));

    mb.push_vertex(mesh_vertex({vec3(-bbox.x, bbox.y, bbox.z),
                                normalize(vec3(-bbox.x, bbox.y, bbox.z)),
                                vec2(0)}));

    mb.push_vertex(mesh_vertex({vec3(bbox.x, -bbox.y, bbox.z),
                                normalize(vec3(bbox.x, -bbox.y, bbox.z)),
                                vec2(0)}));

    mb.push_vertex(mesh_vertex({vec3(-bbox.x, -bbox.y, bbox.z),
                                normalize(vec3(-bbox.x, -bbox.y, bbox.z)),
                                vec2(0)}));

    mb.push_vertex(mesh_vertex({vec3(bbox.x, bbox.y, -bbox.z),
                                normalize(vec3(bbox.x, bbox.y, -bbox.z)),
                                vec2(0)}));

    mb.push_vertex(mesh_vertex({vec3(-bbox.x, bbox.y, -bbox.z),
                                normalize(vec3(-bbox.x, bbox.y, -bbox.z)),
                                vec2(0)}));

    mb.push_vertex(mesh_vertex({vec3(bbox.x, -bbox.y, -bbox.z),
                                normalize(vec3(bbox.x, -bbox.y, -bbox.z)),
                                vec2(0)}));

    mb.push_vertex(mesh_vertex({vec3(-bbox.x, -bbox.y, -bbox.z),
                                normalize(vec3(-bbox.x, -bbox.y, -bbox.z)),
                                vec2(0)}));

    mb.push_indices({0, 1,
                     0, 2,
                     1, 3,
                     2, 3,
                     4, 5,
                     4, 6,
                     5, 7,
                     6, 7,
                     0, 4,
                     1, 5,
                     2, 6,
                     3, 7});


    m_bbox_mesh = mb.build();
}

void Application::createGroundplane() {

    mesh_builder mb;
    mb.push_vertex(mesh_vertex({vec3(200, 0, 200),
                                vec3(0, 1, 0),
                                vec2(0)}));

    mb.push_vertex(mesh_vertex({vec3(-200, 0, 200),
                                vec3(0, 1, 0),
                                vec2(0)}));

    mb.push_vertex(mesh_vertex({vec3(-200, 0, -200),
                                vec3(0, 1, 0),
                                vec2(0)}));

    mb.push_vertex(mesh_vertex({vec3(200, 0, -200),
                                vec3(0, 1, 0),
                                vec2(0)}));
    mb.push_indices({0, 2, 1,
                     0, 3, 2});

    m_ground_plane_mesh = mb.build();
}


/**  ========================================== End of Robert Functions ==========================================*/


void Application::renderGUI() {
	if (ImGui::BeginMainMenuBar()) {
		if (ImGui::BeginMenu("Options")) {
			if (ImGui::MenuItem("Toggle View Menu")) { m_view = !m_view; }

			ImGui::Separator();

			if (ImGui::MenuItem("Toggle Shader Menu")) { m_shade = !m_shade; }
			if (ImGui::MenuItem("Toggle Simulation Menu")) { m_sim = !m_sim; }

			ImGui::Separator();

			if (ImGui::MenuItem("Show all")) {
				m_view = true; m_shade = true; m_sim = true;
			}
			if (ImGui::MenuItem("Hide all")) {
				m_view = false; m_shade = false; m_sim = false;
			}
			ImGui::EndMenu();
		}
		ImGui::EndMainMenuBar();
	}

	if (m_view) { showViewOptions(); }
	if (m_shade) { showShaderOptions(); }
	if (m_sim) { showSoftBodyOptions(); }
	showModeChanger();
}


void Application::showViewOptions() {
	// setup window
	ImGui::SetNextWindowPos(ImVec2(5, 25), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ImVec2(300, 165), ImGuiSetCond_Once);
	ImGui::Begin("View Options", 0);

	// display current camera parameters
	ImGui::Text("Application %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::SliderFloat("Pitch", &m_pitch, -pi<float>() / 2, pi<float>() / 2, "%.2f");
	ImGui::SliderFloat("Yaw", &m_yaw, -pi<float>(), pi<float>(), "%.2f");
	ImGui::SliderFloat("Distance", &m_distance, 0, 100, "%.2f", 2.0f);

	// helpful drawing options
	ImGui::Checkbox("Show axis", &m_show_axis);
	ImGui::SameLine();
	ImGui::Checkbox("Show grid", &m_show_grid);
	ImGui::Checkbox("Wireframe", &m_showWireframe);
	ImGui::SameLine();
	if (ImGui::Button("Screenshot")) rgba_image::screenshot(true);

	// finish creating window
	ImGui::End();
}


void Application::showShaderOptions() {
	// setup window
	ImGui::SetNextWindowPos(ImVec2(5, 195), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ImVec2(300, 170), ImGuiSetCond_Once);
	ImGui::Begin("Shader Options", 0);

	ImGui::PushItemWidth(-120);

	ImGui::SliderInt("Min Thickness", &m_min, 10, m_max);
	ImGui::SliderInt("Max Thickness", &m_max, m_min, 2000);

	ImGui::SliderFloat("Light Intensity", &m_intensity, 0.0f, 1.0f, "%.2f");
	ImGui::SliderFloat("Opacity", &m_opacity, 0.0f, 1.0f, "%.2f");

	ImGui::Checkbox("Flow", &m_model.flow);
	ImGui::SameLine();
	ImGui::SliderFloat("Speed", &m_speed, 0.1f, 10.0f, "%.2f");

	if (ImGui::Combo("Cube Map", &m_map, m_map_options, 10)) {
		setUpCubeMap(m_map_options[m_map]);
	}

	// finish creating window
	ImGui::End();
}


void Application::showSoftBodyOptions() {
	// setup window
	ImGui::SetNextWindowPos(ImVec2(5, 370), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ImVec2(300, 189), ImGuiSetCond_Once);
	ImGui::Begin("Simulation Options", 0);

	ImGui::PushItemWidth(-120);
    if (ImGui::SliderFloat("Gravity", &m_gravity, 0, 2)){
        for (auto &softbody : m_softbodies) {
            softbody.m_gravity = m_gravity;
        }
    }
    if (ImGui::SliderFloat("Mass", &m_mass, 0.8, 1.5)){
        for (auto &softbody : m_softbodies) {
            softbody.m_mass = m_mass;
        }
    }
    if (ImGui::SliderFloat("Damping factor", &m_kd, 0, 6)){
        for (auto &softbody : m_softbodies) {
            softbody.m_kd = m_kd;
        }
    }
    if (ImGui::SliderFloat("Spring factor", &m_ks, 0, 9)){
        for (auto &softbody : m_softbodies) {
            softbody.m_ks = m_ks;
        }
    }
	if (ImGui::SliderFloat("Pressure factor", &m_pressure, 0, 400)){
        for (auto &softbody : m_softbodies) {
            softbody.m_pressure = m_pressure;
        }
    }
    ImGui::Checkbox("Place Soft Bodies", &m_place_softbodies);

	if (ImGui::Button("Restart Simulation")) {
        for (auto &softbody : m_softbodies)
            softbody.resetSimulation();
	}

	// finish creating window
	ImGui::End();
}

void Application::showModeChanger() {
	ImGui::SetNextWindowPos(ImVec2(m_windowsize.x - 230, 25), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ImVec2(225, 50), ImGuiSetCond_Once);
	ImGui::Begin("Modes", 0);

    ImGui::PushItemWidth(-60);
	if (ImGui::Combo("Mode", &m_current_mode, m_mode_options, 3)){
        for (auto &softbody : m_softbodies)
            softbody.resetSimulation();
	    if (m_current_mode == Simulation){
             m_gravity = 1;
             m_mass = 1.0;
             m_kd = 6;
             m_ks = 8;
             m_pressure = 300;

             // clear all somebodies
             m_softbodies.erase(m_softbodies.begin(), m_softbodies.end());

             m_softbodies.emplace_back();
             mat4 initialPositionTransform = translate(mat4(1.0f), vec3(0, 2, 0)) * scale(mat4(1.0f), vec3(m_ball_radius));
             m_softbodies.at(0).initializeMesh(m_mesh, initialPositionTransform);

	    } else if (m_current_mode == FullDemo) {

            m_gravity = 0.025;
            m_mass = 1.0;
            m_kd = 6;
            m_ks = 8;
            m_pressure = 300;

            vec3 min = vec3(-bbox.x * 0.75, -bbox.y * 0.75, -bbox.z * 0.75);
            vec3 max = vec3(bbox.x * 0.75, bbox.y * 0.75, bbox.z * 0.75);

            for (int i = 0; i < 21; ++i) {
                vec3 pos = glm::linearRand(min, max);
                mat4 initialPositionTransform = translate(mat4(1.0f), pos) * scale(mat4(1.0f), vec3(m_ball_radius));
                addNewSoftbody(initialPositionTransform, i == 20);
            }
	    }
	}

	ImGui::End();
}


void Application::cursorPosCallback(double xpos, double ypos) {
	if (m_leftMouseDown) {
		vec2 whsize = m_windowsize / 2.0f;

		// clamp the pitch to [-pi/2, pi/2]
		m_pitch += float(acos(glm::clamp((m_mousePosition.y - whsize.y) / whsize.y, -1.0f, 1.0f))
			- acos(glm::clamp((float(ypos) - whsize.y) / whsize.y, -1.0f, 1.0f)));
		m_pitch = float(glm::clamp(m_pitch, -pi<float>() / 2, pi<float>() / 2));

		// wrap the yaw to [-pi, pi]
		m_yaw += float(acos(glm::clamp((m_mousePosition.x - whsize.x) / whsize.x, -1.0f, 1.0f))
			- acos(glm::clamp((float(xpos) - whsize.x) / whsize.x, -1.0f, 1.0f)));
		if (m_yaw > pi<float>()) m_yaw -= float(2 * pi<float>());
		else if (m_yaw < -pi<float>()) m_yaw += float(2 * pi<float>());
	}

	// updated mouse position
	m_mousePosition = vec2(xpos, ypos);
}


void Application::mouseButtonCallback(int button, int action, int mods) {
	(void)mods; // currently un-used

	// capture is left-mouse down
	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		m_leftMouseDown = (action == GLFW_PRESS); // only other option is GLFW_RELEASE
		static double lastDown = 999999;

		if (action == GLFW_PRESS) {
			lastDown = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();
		}
		if (action == GLFW_RELEASE && chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count() - lastDown < 500) {
            // if not in simulation/final mode, skip this
		    if (m_current_mode == Shader) return;

			// https://stackoverflow.com/a/30005258
			/** ================================================================================================ */
			// Normalized Device Coordinates
			float mouseX = m_mousePosition.x / (m_windowsize.x * 0.5f) - 1.0f;
			float mouseY = m_mousePosition.y / (m_windowsize.y * 0.5f) - 1.0f;

			mat4 proj = perspective(1.f, float(m_windowsize.x) / m_windowsize.y, 0.3f, 1000.f);
			mat4 view = translate(mat4(1), vec3(3, 0, -m_distance))
				* rotate(mat4(1), m_pitch, vec3(1, 0, 0))
				* rotate(mat4(1), m_yaw, vec3(0, 1, 0));

			// matrix to get from screen coords to world coords
			glm::mat4 invVP = glm::inverse(proj * view);
			glm::vec4 screenPos = glm::vec4(mouseX, -mouseY, 1.0f, 1.0f);
			glm::vec4 worldPos = invVP * screenPos;

			glm::vec3 direction = glm::normalize(glm::vec3(worldPos));

			/** ================================================================================================ */

			vec4 viewport = vec4(0, 0, m_windowsize.x, m_windowsize.y);
			vec3 mousePos = vec3(m_mousePosition.x, viewport.w - m_mousePosition.y, 0.01);

			vec3 cameraPos = unProject(mousePos, view, proj, viewport);

			vec3 rayDestination = cameraPos + direction * glm::max((m_distance - bbox.x/2), 20.0f);

			if (m_place_softbodies){
                mat4 initialPositionTransform = translate(mat4(1.0f), rayDestination) * scale(mat4(1.0f), vec3(m_ball_radius));
                addNewSoftbody(initialPositionTransform, true);
            } else {
                for (auto &softbody : m_softbodies) {
                    softbody.applyClick(cameraPos, rayDestination, direction, m_ball_radius);
                }
			}
		}
	}

}


void Application::scrollCallback(double xoffset, double yoffset) {
	(void)xoffset; // currently un-used
	m_distance *= pow(1.1f, -yoffset);
}


void Application::keyCallback(int key, int scancode, int action, int mods) {
	(void)key, (void)scancode, (void)action, (void)mods; // currently un-used
}


void Application::charCallback(unsigned int c) {
	(void)c; // currently un-used
}

// loads a cubemap texture from 6 individual texture faces
// Tutorial source := www.learnopengl.com
// -------------------------------------------------------
unsigned int Application::loadCubemap(vector<std::string> cubeFaces) {
	unsigned int id;
	glGenTextures(1, &id);
	glBindTexture(GL_TEXTURE_CUBE_MAP, id);

	int width, height, channels;
	for (unsigned int i = 0; i < cubeFaces.size(); i++) {
		unsigned char* imageData = stbi_load(cubeFaces[i].c_str(), &width, &height, &channels, 0);
		if (imageData) {
			glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, imageData);
			stbi_image_free(imageData);
		}
		else {
			std::cout << "Cubemap texture failed to load at path: " << cubeFaces[i] << std::endl;
			stbi_image_free(imageData);
		}
	}
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

	return id;
}

void Application::setUpCubeMap(char* name) {
	// set up the cubemap
	std::string cubeMapPath = CGRA_SRCDIR + std::string("//res//textures//cube_maps//") + std::string(name);
	vector<std::string> faces{
		(cubeMapPath + std::string("//posx.jpg")),
		(cubeMapPath + std::string("//negx.jpg")),
		(cubeMapPath + std::string("//negy.jpg")),
		(cubeMapPath + std::string("//posy.jpg")),
		(cubeMapPath + std::string("//posz.jpg")),
		(cubeMapPath + std::string("//negz.jpg")),
	};
	m_model.cubeMap = loadCubemap(faces);
}

