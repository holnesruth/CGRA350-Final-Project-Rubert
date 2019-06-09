
// std
#include <iostream>
#include <string>
#include <chrono>

// glm
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

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


void basic_model::draw(const glm::mat4& view, const glm::mat4 proj) {
	mat4 modelview = view * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelMatrix"), 1, false, value_ptr(modelTransform));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));
	glUniform2fv(glGetUniformLocation(shader, "uThickness"), 1, value_ptr(tParams));
	glUniform2fv(glGetUniformLocation(shader, "uLightEffects"), 1, value_ptr(leParams));
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, gradient);
	glUniform1i(glGetUniformLocation(shader, "uTexture"), 0);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_CUBE_MAP, cubeMap);
	glUniform1i(glGetUniformLocation(shader, "uCubeMap"), 1);

	drawSphere();
}


void basic_model::updateParams(float minThickness, float maxThickness, float intensity, float opacity) {
	tParams.x = minThickness;
	tParams.y = maxThickness;
	leParams.x = intensity;
	leParams.y = opacity;
}


Application::Application(GLFWwindow* window) : m_window(window) {

	shader_builder sb;
	sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//simple_vert.glsl"));
	sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//bubbles.glsl"));
	GLuint shader = sb.build();
	m_model.shader = shader;

	// set up the thickness texture
	rgba_image img = rgba_image(CGRA_SRCDIR + std::string("//res//textures//grad.jpg"));
	m_model.gradient = img.uploadTexture();

	setUpCubeMap("Skansen");
	m_model.modelTransform = rotate(mat4(1.0), radians(90.0f), vec3(1, 0, 0));
}


void Application::render() {

	// retrieve the window hieght
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

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// projection matrix
	mat4 proj = perspective(1.f, float(width) / height, 0.1f, 1000.f);

	// view matrix
	mat4 view = translate(mat4(1), vec3(0, 0, -m_distance))
		* rotate(mat4(1), m_pitch, vec3(1, 0, 0))
		* rotate(mat4(1), m_yaw, vec3(0, 1, 0));


	// helpful draw options
	if (m_show_grid) drawGrid(view, proj);
	if (m_show_axis) drawAxis(view, proj);
	glPolygonMode(GL_FRONT_AND_BACK, (m_showWireframe) ? GL_LINE : GL_FILL);

	m_model.updateParams(m_min, m_max, m_intensity, m_opacity);

	// draw the model
	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);

	m_model.leParams.y *= m_intensity;
	m_model.draw(view, proj);

	glCullFace(GL_BACK);

	m_model.leParams.y = m_opacity;
	m_model.draw(view, proj);
}


void Application::renderGUI() {
	if (ImGui::BeginMainMenuBar()) {
		if (ImGui::BeginMenu("Options")) {
			if (ImGui::MenuItem("View")) { m_view = !m_view; }

			ImGui::Separator();

			if (ImGui::MenuItem("Shader")) { m_shade = !m_shade; }
			if (ImGui::MenuItem("Simulation")) { m_sim = !m_sim; }

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
	ImGui::SetNextWindowSize(ImVec2(300, 150), ImGuiSetCond_Once);
	ImGui::Begin("Shader Options", 0);

	ImGui::PushItemWidth(-120);
	ImGui::SliderInt("Min Thickness", &m_min, 10, m_max);
	ImGui::SliderInt("Max Thickness", &m_max, m_min, 2000);
	ImGui::SliderFloat("Light Intensity", &m_intensity, 0, 1, "%.2f");
	ImGui::SliderFloat("Opacity", &m_opacity, 0, 1, "%.2f");

	if (ImGui::Combo("Cube Map", &m_map, m_map_options, 10)) {
		setUpCubeMap(m_map_options[m_map]);
	}

	// finish creating window
	ImGui::End();
}


void Application::showSoftBodyOptions() {
	// setup window
	ImGui::SetNextWindowPos(ImVec2(5, 350), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ImVec2(300, 100), ImGuiSetCond_Once);
	ImGui::Begin("Simulation Options", 0);

	// finish creating window
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
	if (button == GLFW_MOUSE_BUTTON_LEFT)
		m_leftMouseDown = (action == GLFW_PRESS); // only other option is GLFW_RELEASE
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
		} else {
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