
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

	// Create soft body mesh
	initializeMesh();
	
	// Setg up model parameters for simulation
	m_model.mesh = constructMesh(m_points, m_springs);
	m_model.color = vec3(1, 1, 1);

	// Get the time
	m_lastMillis = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();

	// Set up the thickness texture
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
	
	// get the current time
	double millis = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();


	// if it's time for another simulation step
	if (millis - m_lastMillis > DT * 1000) {
		AccumulateForces();

		IntegrateForces();

		m_model.mesh = constructMesh(m_points, m_springs);

		m_lastMillis = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count();

		// update flow noise parameter
		m_model.time += m_speed/100;
		if (m_model.time > 150) {
			m_model.time = 100; // reset
		}
	}

	// draw spheres on all points
	if (m_current_mode != Shader && m_showWireframe) {
		GLuint originalShader = m_model.shader;
		if (originalShader != m_shader_default) {
			m_model.shader = m_shader_default;
		}

		for (auto &point : m_points) {
			m_model.modelTransform = translate(mat4(1.0f), point.pos) * scale(mat4(1.0f), vec3(0.1f));
			m_model.draw(view, proj, true);
		}

		m_model.shader = originalShader;
	}

	if (m_current_mode == Simulation) {
		m_model.shader = m_shader_default;
	}
	else {
		m_model.shader = m_shader_bubble;
	}

	// draw the mesh
	if (m_current_mode == Shader) {
		m_model.modelTransform = scale(rotate(mat4(1.0f), radians(90.0f), vec3(1, 0, 0)), vec3(m_ball_radius)); 
	} else {
		m_model.modelTransform = mat4(1.0);
	}

	if (m_show_grid) drawGrid(view, proj);
	if (m_show_axis) drawAxis(view, proj);

	// draw the model
	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);

	m_model.leParams.y *= m_intensity;
	m_model.draw(view, proj, (m_current_mode == Shader));

	glCullFace(GL_BACK);

	m_model.leParams.y = m_opacity;
	m_model.draw(view, proj, (m_current_mode == Shader));
}


/**  ========================================== Robert's Functions ==========================================*/

void Application::AccumulateForces() {

	/** ============== Accumulate gravity force ============================ */
	for (auto& m_point : m_points) {
		m_point.force = vec3(0, -m_gravity * m_mass, 0);
	}

	/** ============== Accumulate spring force ============================ */
	for (auto& spring : m_springs) {
		Point& p1 = m_points.at(spring.index1);
		Point& p2 = m_points.at(spring.index2);

		//        if (p1.mass == -1) continue;
		//        if (p2.mass == -1) continue;

		//            float dist = distance(p1.pos, p2.pos);
		//            if (dist == 0) continue;
		//
		//            vec3 forc = (dist -  spring.length) * m_ks + (p1.vel - p2.vel) * ((p1.pos - p2.pos) / dist) * m_kd;
		//
		//            p1.force -= f;
		//            p2.force += f;

				// ====================================================
		float x1 = p1.pos.x;	float y1 = p1.pos.y; float z1 = p1.pos.z;
		float x2 = p2.pos.x;	float y2 = p2.pos.y; float z2 = p2.pos.z;

		float r12d = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));

		if (r12d == 0) continue;

		float vx12 = p1.vel.x - p2.vel.x;
		float vy12 = p1.vel.y - p2.vel.y;
		float vz12 = p1.vel.z - p2.vel.z;

		float f = (r12d - spring.length) * m_ks + (vx12 * (x1 - x2) + vy12 * (y1 - y2) + vz12 * (z1 - z2)) * m_kd / r12d;


		float Fx = ((x1 - x2) / r12d) * f;
		float Fy = ((y1 - y2) / r12d) * f;
		float Fz = ((z1 - z2) / r12d) * f;

		p1.force -= vec3(Fx, Fy, Fz);
		p2.force += vec3(Fx, Fy, Fz);

	}

	/** ============== Calculate Volume ============================ */
//    Robust method for calculating volume from Frank Krueger at https://stackoverflow.com/a/1568551
	float volume = 0;
	// for each triangle
	for (int i = 0; i < m_springs.size(); i += 3) {
		vec3 p1 = m_points.at(m_springs.at(i).index1).pos;
		vec3 p2 = m_points.at(m_springs.at(i + 1).index1).pos;
		vec3 p3 = m_points.at(m_springs.at(i + 2).index1).pos;

		volume += dot(p1, cross(p2, p3)) / 6.0f;

	}


	/** ============== Accumulate pressure force ============================ */
	for (int i = 0; i < m_springs.size(); i += 3) {
		Point p1 = m_points.at(m_springs.at(i).index1);
		Point p2 = m_points.at(m_springs.at(i + 1).index1);
		Point p3 = m_points.at(m_springs.at(i + 2).index1);

		vec3 p1_p2 = p2.pos - p1.pos;
		vec3 p1_p3 = p3.pos - p1.pos;

		float faceSize = length(cross(p1_p2, p1_p3)) / 2;

		float pressureValue = faceSize * m_pressure * (1.0f / volume);

		for (int j = 0; j < 3; ++j) {
			Spring s = m_springs.at(i + j);
			m_points.at(s.index1).force += pressureValue * s.normal;
			m_points.at(s.index2).force += pressureValue * s.normal;
		}
	}
}

void Application::IntegrateForces() {
	// integrate forces
	for (auto& pnt : m_points) {
		pnt.vel += (pnt.force / m_mass) * DT;
		//        float lastY = pnt.pos.y;
		pnt.pos += pnt.vel;

		if (pnt.pos.y <= 0) {
			pnt.pos.y = 0;
			// if sliding across the ground
			pnt.vel.y *= -0.5;

		}
	}
}

void Application::initializeMesh() {
	mesh_builder mesh = load_wavefront_data(CGRA_SRCDIR + std::string("/res//assets//ball_270.obj"));

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

	/** ============== Construct points and springs from clean mesh data ============== */

	mat4 initialPositionTransform = translate(mat4(1.0f), vec3(0, 2, 0)) * scale(mat4(1.0f), vec3(m_ball_radius));

	for (auto& meshVertex : mesh.vertices) {
		vec4 initialPos = vec4(meshVertex.pos.x, meshVertex.pos.y, meshVertex.pos.z, 1);
		vec4 newPos = initialPositionTransform * initialPos;
		m_points.push_back(Point({ vec3(newPos), vec3(0.0f), vec3(0.0f), 0.0f }));
		m_restPos.push_back(Point({ vec3(newPos), vec3(0.0f), vec3(0.0f), 0.0f }));
	}

	for (int i = 0; i < mesh.indices.size(); i += 3) {
		int index1 = mesh.indices.at(i);
		int index2 = mesh.indices.at(i + 1);
		int index3 = mesh.indices.at(i + 2);

		mesh_vertex p1 = mesh.vertices.at(index1);
		mesh_vertex p2 = mesh.vertices.at(index2);
		mesh_vertex p3 = mesh.vertices.at(index3);

		vec3 pos1 = m_points.at(index1).pos;
		vec3 pos2 = m_points.at(index2).pos;
		vec3 pos3 = m_points.at(index3).pos;


		m_springs.push_back(Spring({ index1, index2, distance(pos1, pos2), (p1.norm + p2.norm) / 2.0f }));
		m_springs.push_back(Spring({ index2, index3, distance(pos2, pos3), (p2.norm + p3.norm) / 2.0f }));
		m_springs.push_back(Spring({ index3, index1, distance(pos3, pos1), (p3.norm + p1.norm) / 2.0f }));
	}


	// erase all spring duplicates
	vector<Spring>::iterator end = m_springs.end();
	for (vector<Spring>::iterator it = m_springs.begin(); it != end; ++it) {
		end = std::remove(it + 1, end, *it);
	}

	m_springs.erase(end, m_springs.end());

	m_mesh = mesh;
}

cgra::gl_mesh Application::constructMesh(std::vector<Point>& points, std::vector<Spring>& springs) {
	// TODO
	if (m_showWireframe) {
		mesh_builder mb(GL_LINES);
		vec3 normal = vec3(0, 0, 1);
		for (auto& point : points) {
			mb.push_vertex(mesh_vertex{
					point.pos,
					normal,
					vec2(0, 0)
				});
		}

		for (auto& spring : springs) {
			mb.push_index(spring.index1);
			mb.push_index(spring.index2);
		}

		return mb.build();
	}
	else {
		for (int i = 0; i < m_points.size(); i++) {
			m_mesh.vertices.at(i).pos = m_points.at(i).pos;
		}
		return m_mesh.build();
	}


}


void Application::resetSimulation() {
	// TODO
//    vec3 min = vec3(-0.05, 0.1, -0.05);
//    vec3 max = vec3(0.05, 0.4, 0.05);

	for (int i = 0; i < m_points.size(); i++) {
		//        p.vel = glm::linearRand(min, max);
		m_points.at(i) = m_restPos.at(i);
	}


}


double computeDistance(vec3 A, vec3 B, vec3 C) {
	vec3 d = (C - B) / distance(C, B);
	vec3 v = A - B;
	float t = dot(v, d);
	vec3 P = B + (t * d);
	return distance(A, P);
}

/**  ========================================== End of Robert's Functions ==========================================*/

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
	ImGui::SetNextWindowSize(ImVec2(300, 170), ImGuiSetCond_Once);
	ImGui::Begin("Simulation Options", 0);

	ImGui::PushItemWidth(-120);
	ImGui::SliderFloat("Gravity", &m_gravity, 0, 3);
	ImGui::SliderFloat("Mass", &m_mass, 0.5, 1.5);
	ImGui::SliderFloat("Damping factor", &m_kd, 0, 10);
	ImGui::SliderFloat("Spring factor", &m_ks, 0, 10);
	ImGui::SliderFloat("Pressure factor", &m_pressure, 0, 200);

	if (ImGui::Button("Restart Simulation")) {
		resetSimulation();
	}

	// finish creating window
	ImGui::End();
}

void Application::showModeChanger() {
	ImGui::SetNextWindowPos(ImVec2(5, 545), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ImVec2(300, 50), ImGuiSetCond_Once);
	ImGui::Begin("Modes", 0);

	ImGui::PushItemWidth(-120);
	ImGui::Combo("Mode", &m_current_mode, m_mode_options, 3);

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
		if (action == GLFW_RELEASE && chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count() - lastDown < 100) {


			// https://stackoverflow.com/a/30005258
			/** ================================================================================================ */
			// Normalized Device Coordinates
			float mouseX = m_mousePosition.x * 3 / (m_windowsize.x * 0.5f) - 1.0f;
			float mouseY = m_mousePosition.y * 2 / (m_windowsize.y * 0.5f) - 1.0f;

			mat4 proj = perspective(1.f, float(m_windowsize.x) / m_windowsize.x, 0.3f, 1000.f);
			mat4 view = translate(mat4(1), vec3(3, 0, -m_distance))
				* rotate(mat4(1), m_pitch, vec3(1, 0, 0))
				* rotate(mat4(1), m_yaw, vec3(0, 1, 0));

			// matrix to get from screen coords to world coords
			glm::mat4 invVP = glm::inverse(proj * view);
			glm::vec4 screenPos = glm::vec4(mouseX, -mouseY, 1.0f, 1.0f);
			glm::vec4 worldPos = invVP * screenPos;

			glm::vec3 direction = glm::normalize(glm::vec3(worldPos));

			/** ================================================================================================ */
//            vec3 cameraPos2 = inverse(view) * vec4(0, 0, 0, 1);

			vec4 viewport = vec4(0, 0, m_windowsize.x, m_windowsize.y);
			vec3 mousePos = vec3(m_mousePosition.x * 3, viewport.w - m_mousePosition.y * 2, 0.01);

			vec3 cameraPos = unProject(mousePos, view, proj, viewport);

			vec3 objPos = m_points.at(0).pos;
			vec3 C = cameraPos + direction * 30.0f;

			uint size = m_mesh.vertices.size() - 1;
			m_mesh.indices.push_back(size + 1);
			m_mesh.indices.push_back(size + 2);
			m_mesh.indices.push_back(size + 3);

			mesh_vertex meshVertex1; meshVertex1.pos = cameraPos; m_mesh.vertices.push_back(meshVertex1);
			mesh_vertex meshVertex2; meshVertex2.pos = C; m_mesh.vertices.push_back(meshVertex2);
			mesh_vertex meshVertex3; meshVertex3.pos = C + vec3(0, -0.5, 0); m_mesh.vertices.push_back(meshVertex3);



			/*float dist = computeDistance(objPos, cameraPos, C);
			cout << "dist: " << dist << endl;

			for (auto& point : m_points) {

			}


			int x = 0;*/
			//            std::cout << "objPos: x:" << objPos.x << ", y:" << objPos.y << ", z:" << objPos.z << endl;
			//            std::cout << "cameraPos: x:" << objPos.x << ", y:" << objPos.y << ", z:" << objPos.z << endl;


			//
			//            vec3 mousePos = vec3(m_mousePosition.x, viewport.w - m_mousePosition.y, 0);
			//            cout << "mouse z: " << mousePos.z << endl;
			//
			//            std::cout << "projectedPos: x:" << newPos.x<< ", y:" << newPos.y << ", z:" << newPos.z << endl << endl;

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