//
// Created by Robert Wilkins on 2019-06-11.
//

#ifndef CGRA_PROJECT_BASE_SOFTBODY_H
#define CGRA_PROJECT_BASE_SOFTBODY_H

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include "opengl.hpp"
#include "cgra/cgra_mesh.hpp"
#include "skeleton_model.hpp"

// Soft body points
struct Point {
    glm::vec3 pos, vel, force;
    float mass;
    glm::vec3 norm;
};

// Soft body springs
struct Spring {
    int index1, index2; // indices for the two points on either end of m_spring
    float length; // rest length
    glm::vec3 normal; // normal vector

    inline bool operator==(const Spring& other)
    {
        return index1 == other.index1 && index2 == other.index2;
    }
};


// simulation constants
static const float DT = 0.01666666667; // time difference
static const float KS = 1755.0f; // spring constant

class Softbody {
public:

    // simulation params
    float m_gravity = 1;
    float m_mass = 1.0;
    float m_kd = 6;
    float m_ks = 6;
    float m_pressure = 100;

    // geometry
    cgra::mesh_builder m_mesh;
    std::vector<Point> m_points;
    std::vector<Spring> m_springs;
    std::vector<Point> m_restPos;


    void applyClick(glm::vec3 cameraPos, glm::vec3 rayDestination, glm::vec3 direction, float ball_radius);


    void AccumulateForces();

    void IntegrateForces(bool useGroundPlane);

    void initializeMesh(cgra::mesh_builder mesh, glm::mat4 initialTransform);

    cgra::gl_mesh constructMesh(bool showWireframe);

    void resetSimulation();

    static double computeDistance(glm::vec3 A, glm::vec3 B, glm::vec3 C);

    static float myMap(float value, float start1, float stop1, float start2, float stop2);
};


#endif //CGRA_PROJECT_BASE_SOFTBODY_H
