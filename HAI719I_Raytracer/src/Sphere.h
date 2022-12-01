#ifndef Sphere_H
#define Sphere_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySphereIntersection
{
    bool intersectionExists;
    float t;
    float theta, phi;
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal;
    // Adding reflexion pour recursive
    Vec3 bounce_direction;
};

static Vec3 SphericalCoordinatesToEuclidean(Vec3 ThetaPhiR)
{
    return ThetaPhiR[2] * Vec3(cos(ThetaPhiR[0]) * cos(ThetaPhiR[1]), sin(ThetaPhiR[0]) * cos(ThetaPhiR[1]), sin(ThetaPhiR[1]));
}
static Vec3 SphericalCoordinatesToEuclidean(float theta, float phi)
{
    return Vec3(cos(theta) * cos(phi), sin(theta) * cos(phi), sin(phi));
}

static Vec3 EuclideanCoordinatesToSpherical(Vec3 xyz)
{
    float R = xyz.length();
    float phi = asin(xyz[2] / R);
    float theta = atan2(xyz[1], xyz[0]);
    return Vec3(theta, phi, R);
}

class Sphere : public Mesh
{
public:
    Vec3 m_center;
    float m_radius;

    Sphere() : Mesh() {}
    Sphere(Vec3 c, float r) : Mesh(), m_center(c), m_radius(r) {}

    void build_arrays()
    {
        unsigned int nTheta = 20, nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi);
        normalsArray.resize(3 * nTheta * nPhi);
        uvs_array.resize(2 * nTheta * nPhi);
        for (unsigned int thetaIt = 0; thetaIt < nTheta; ++thetaIt)
        {
            float u = (float)(thetaIt) / (float)(nTheta - 1);
            float theta = u * 2 * M_PI;
            for (unsigned int phiIt = 0; phiIt < nPhi; ++phiIt)
            {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float)(phiIt) / (float)(nPhi - 1);
                float phi = -M_PI / 2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean(theta, phi);
                positions_array[3 * vertexIndex + 0] = m_center[0] + m_radius * xyz[0];
                positions_array[3 * vertexIndex + 1] = m_center[1] + m_radius * xyz[1];
                positions_array[3 * vertexIndex + 2] = m_center[2] + m_radius * xyz[2];
                normalsArray[3 * vertexIndex + 0] = xyz[0];
                normalsArray[3 * vertexIndex + 1] = xyz[1];
                normalsArray[3 * vertexIndex + 2] = xyz[2];
                uvs_array[2 * vertexIndex + 0] = u;
                uvs_array[2 * vertexIndex + 1] = v;
            }
        }
        triangles_array.clear();
        for (unsigned int thetaIt = 0; thetaIt < nTheta - 1; ++thetaIt)
        {
            for (unsigned int phiIt = 0; phiIt < nPhi - 1; ++phiIt)
            {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt + 1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt + 1) * nTheta;
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuV);
            }
        }
    }

    RaySphereIntersection intersect(const Ray &ray) const
    {
        RaySphereIntersection intersection;
        Vec3 d = ray.direction();
        Vec3 o = ray.origin();
        Vec3 c = m_center;
        double r = m_radius;

        // Equation d'intersection rayon-sphère
        //  t2 d.d + 2t d.(o-c) + ||o-c||² – r² = 0
        double A, B, C;
        A = d.dot(d, d);
        B = 2 * (d.dot(d, o - c));
        C = ((o - c).length() * (o - c).length()) - r * r;

        // discriminant : b*b - 4*a*c
        double Discriminant;
        Discriminant = B * B - (4 * A * C);

        // Si on a des racines
        double racine1, racine2;
        if (Discriminant >= 0)
        {

            racine2 = fmax(-B + sqrt(Discriminant), -B - sqrt(Discriminant));
            racine1 = fmin(-B + sqrt(Discriminant), -B - sqrt(Discriminant));

            racine1 /= 2 * A;
            racine2 /= 2 * A;

            intersection.intersectionExists = true;
            intersection.t = racine1;
            intersection.intersection = o + (d * racine1);
            intersection.secondintersection = o + (d * racine2);
            intersection.normal = intersection.intersection - c;
            intersection.normal.normalize();
            intersection.bounce_direction = d - 2 * intersection.normal * (Vec3::dot(d,intersection.normal)); // Trouvé sur internet pour avoir le rayon reflechi
            intersection.bounce_direction.normalize();
        }
        else
        {
            intersection.intersectionExists = false;
        }
        return intersection;
    }
};
#endif
