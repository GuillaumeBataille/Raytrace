#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

double EPSILON_up = 0.00001;
double EPSILON_down = 0.0001;

struct RaySquareIntersection
{
    bool intersectionExists;
    float t;
    float u, v;
    Vec3 intersection;
    Vec3 normal;
    // Adding reflexion pour recursive
    Vec3 bounce_direction;
};

double triangleArea(Vec3 A, Vec3 B, Vec3 C)
{
    return (C[0] * B[1] - B[0] * C[1]) - (C[0] * A[1] - A[0] * C[1]) + (B[0] * A[1] - A[0] * B[1]);
}

class Square : public Mesh
{
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const &bottomLeft, Vec3 const &rightVector, Vec3 const &upVector, float width = 1., float height = 1.,
           float uMin = 0.f, float uMax = 1.f, float vMin = 0.f, float vMax = 1.f) : Mesh()
    {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad(Vec3 const &bottomLeft, Vec3 const &rightVector, Vec3 const &upVector, float width = 1., float height = 1.,
                 float uMin = 0.f, float uMax = 1.f, float vMin = 0.f, float vMax = 1.f)
    {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector, upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector * width;
        m_up_vector = m_up_vector * height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;
        vertices[0].u = uMin;
        vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;
        vertices[1].u = uMax;
        vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;
        vertices[2].u = uMax;
        vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;
        vertices[3].u = uMin;
        vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;
    }

    Vec3 bottomLeft() const { return this->vertices[0].position; }
    Vec3 bottomRight() const { return this->vertices[1].position; }
    Vec3 upRight() const { return this->vertices[2].position; }
    Vec3 upLeft() const { return this->vertices[3].position; }
    Vec3 normalvec() const { return Vec3::cross((bottomRight() - bottomLeft()), (upLeft() - bottomLeft())); }
    float getumin() const { return this->vertices[0].u; }
    float getvmin() const { return this->vertices[0].v; }
    RaySquareIntersection intersect(const Ray &ray) const
    {
        RaySquareIntersection intersection;

        Vec3 d = ray.direction();
        Vec3 o = ray.origin();
        Vec3 n = normalvec();
        n.normalize();

        Vec3 a = bottomLeft();   // D
        Vec3 a1 = bottomRight(); // C
        Vec3 a2 = upRight();     // B
        Vec3 a3 = upLeft();      // A

        Vec3 inter;

        double D = a.dot(a, n);

        /* if (Vec3::dot(n, d) <= EPSILON_up && Vec3::dot(n, d) >= -EPSILON_down)
        {
            intersection.intersectionExists = false;
            return intersection;
        }*/

        double t = (D - Vec3::dot(o, n)) / Vec3::dot(d, n); // Fonction de t
        double orientation = Vec3::dot(d,  n);
        if (t > 0 && orientation <= 0) // On est dans le plan et on est face au plan (pour ne pas render le dos des squares inutilement)
        //if (t>0)
        {
            inter = o + (t * d); // Coordonnée du point d'intersection
            // Changement de repère : on se mets dans le repère du square
            // ABx = l'axe X du carré
            // ABy = l'axe Y du carré
            Vec3 ABx = a1 - a;
            Vec3 ABy = a3 - a;
            // Le vecteur entre le 0,0,0 du carré et le point intersection
            Vec3 AP = inter - a;
            // Les normes de projection du point intersection sur l'axe X et Y
            double APx = Vec3::dot(ABx, AP) / ABx.norm();
            double APy = Vec3::dot(ABy, AP) / ABy.norm();
            double ABx_norm = ABx.norm();
            double ABy_norm = ABy.norm();
            // Si la norme Proj x < norme ABx alors on est dans le carré en X (cf en y)
            if (APx <= ABx.norm() && APx > 0 && APy <= ABy.norm() && APy > 0)
            {
                intersection.intersectionExists = true;
                intersection.t = t;
                intersection.normal = n;
                intersection.intersection = inter;
                intersection.u = APx;
                intersection.v = APy;
                intersection.bounce_direction = d - 2 * intersection.normal *( Vec3::dot(d,intersection.normal));
            }
            else
            {
                intersection.intersectionExists = false;
            }
        }
        return intersection;
    }
};
#endif // SQUARE_H
