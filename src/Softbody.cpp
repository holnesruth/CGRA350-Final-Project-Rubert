//
// Created by Robert Wilkins on 2019-06-11.
//

#include "Softbody.h"

// glm
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/random.hpp>
#include <algorithm>


using namespace std;
using namespace cgra;
using namespace glm;


double Softbody::computeDistance(vec3 A, vec3 B, vec3 C) {

    vec3 d = (C - B) / distance(C, B);
    vec3 v = A - B;
    float t = dot(v, d);
    vec3 P = B + (t * d);

    return distance(A, P);
}


float Softbody::myMap(float value, float start1, float stop1, float start2, float stop2) {

    float val = start2 + (stop2 - start2) * ((value - start1) / (stop1 - start1));
    if (val > stop2)
        return stop2;
    else
        return val;
}


void Softbody::applyClick(glm::vec3 cameraPos, glm::vec3 rayDestination, glm::vec3 direction, float ball_radius) {
    
    for (auto& point : m_points) {
        float dist = computeDistance(point.pos, cameraPos, rayDestination);
        float vel =  1.3 - myMap(dist, 0, ball_radius / 1.5, 0, 1.3);
        point.vel += direction * vel;
    }
}


void Softbody::AccumulateForces() {

    /** ============== Accumulate gravity force ============================ */
    for (auto& m_point : m_points) {
        m_point.force = vec3(0, -m_gravity * m_mass, 0);
    }

    /** ============== Accumulate spring force ============================ */
    for (auto& spring : m_springs) {
        Point& p1 = m_points.at(spring.index1);
        Point& p2 = m_points.at(spring.index2);

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
	// Robust method for calculating volume from Frank Krueger at https://stackoverflow.com/a/1568551
    double volume = 0;
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


void Softbody::IntegrateForces(bool useGroundPlane, std::vector<Softbody> &softbodies, float ball_radius) {

    // integrate forces
    for (auto& pnt : m_points) {
        pnt.vel += (pnt.force / m_mass) * DT;
        //        float lastY = pnt.pos.y;
        pnt.pos += pnt.vel;

        if (useGroundPlane) {
            if (pnt.pos.y <= 0) {
                pnt.pos.y = 0;
                pnt.vel.y *= -0.5;
                pnt.vel.x *= 0.99;
                pnt.vel.z *= 0.99;
            }
        }
        if (abs(pnt.pos.x) > bbox.x){
            pnt.pos.x = pnt.pos.x / abs(pnt.pos.x) * bbox.x;
            pnt.vel.x *= -0.9;
        } else if (abs(pnt.pos.z) > bbox.z){
            pnt.pos.z = pnt.pos.z / abs(pnt.pos.z) * bbox.z;
            pnt.vel.z *= -0.9;
        }


        // drag force
        pnt.vel *= 0.995;

        // collision detection and handling
        for (auto &otherSoftbody : softbodies) {
            if (&otherSoftbody != this){
                vec3 otherCentroid = otherSoftbody.m_centroid;
                float dist = distance(pnt.pos, otherCentroid);
                if (dist <= ball_radius + 1) {
                    float amountToMove = abs(dist - ball_radius);
                    vec3 directionToMove = normalize(pnt.pos - otherCentroid);
                    pnt.pos += directionToMove * amountToMove;
                    pnt.vel = directionToMove * length(pnt.vel) * 0.8f;

                }
            }
        }
    }

    // if using bbox
    if (!useGroundPlane) {
        if (abs(m_centroid.y) > bbox.y){
            float delta = normalize(m_centroid.y) * bbox.y * -2;
            for (auto &pnt : m_points){
                pnt.pos.y += delta;
            }
        }
    }
}


void Softbody::updateCentroid(){

    vec3 total = vec3(0);
    for (auto &point : m_points) {
        total+=point.pos;
    }
    m_centroid = total / (float) m_points.size();
}


/* Construct points and springs from clean mesh data */
void Softbody::initializeMesh(mesh_builder mesh, mat4 initialTransform) {

    for (auto& meshVertex : mesh.vertices) {
        vec4 initialPos = vec4(meshVertex.pos.x, meshVertex.pos.y, meshVertex.pos.z, 1);
        vec4 newPos = initialTransform * initialPos;
        m_points.push_back(Point({ vec3(newPos), vec3(0.0f), vec3(0.0f), 0.0f, meshVertex.norm}));
        m_restPos.push_back(Point({ vec3(newPos), vec3(0.0f), vec3(0.0f), 0.0f, meshVertex.norm}));
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


cgra::gl_mesh Softbody::constructMesh(bool showWireframe) {

    if (showWireframe) {
        mesh_builder mb(GL_LINES);
        vec3 normal = vec3(0, 0, 1);
        for (auto& point : m_points) {
            mb.push_vertex(mesh_vertex{
                    point.pos,
                    normal,
                    vec2(0, 0)
            });
        }

        for (auto& spring : m_springs) {
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


void Softbody::resetSimulation() {
    for (int i = 0; i < m_points.size(); i++) {
        m_points.at(i) = m_restPos.at(i);
    }
}
