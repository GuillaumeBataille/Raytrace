#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"



struct RayTriangleIntersection{
    bool intersectionExists;
    float t;
    float w0,w1,w2;
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
    // Adding reflexion pour recursive
    Vec3 bounce_direction;
};

class Triangle {
private:
    Vec3 m_c[3] , m_normal;
    float area;
public:
    Triangle() {}
    Triangle( Vec3 const & c0 , Vec3 const & c1 , Vec3 const & c2 ) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }
    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross( m_c[1] - m_c[0] , m_c[2] - m_c[0] );
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }
    void setC0( Vec3 const & c0 ) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1( Vec3 const & c1 ) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2( Vec3 const & c2 ) { m_c[2] = c2; } // remember to update the area and normal afterwards!
    Vec3 const & normal() const { return m_normal; }
    Vec3 projectOnSupportPlane( Vec3 const & p ) const {
        Vec3 result;
        //TODO completer
        return result;
    }
    float squareDistanceToSupportPlane( Vec3 const & p ) const {
        float result;
        //TODO completer
        return result;
    }
    float distanceToSupportPlane( Vec3 const & p ) const { return sqrt( squareDistanceToSupportPlane(p) ); }
    bool isParallelTo( Line const & L ) const {
        bool result;
        //TODO completer
        return result;
    }
    Vec3 getIntersectionPointWithSupportPlane( Line const & L ) const {
        // you should check first that the line is not parallel to the plane!
        Vec3 result;
        //TODO completer
        return result;
    }


    Vec3 gets0() const { return this->m_c[0]; } // Sommet 0 courant
    Vec3 gets1() const { return this->m_c[1]; } // Sommet 1 courant
    Vec3 gets2() const { return this->m_c[2]; } // Sommet 2 courant
    Vec3 getnorm() const { return this->m_normal; } // Normale courante
    float getarea() const { return this->area; } // area courante


    void computeBarycentricCoordinates( Vec3 const & p , float & u0 , float & u1 , float & u2 ) const {
        //TODO Complete
    }

    RayTriangleIntersection getIntersection( Ray const & ray ) const {
        
        RayTriangleIntersection intersection;
        Vec3 d = ray.direction();
        Vec3 o = ray.origin();
        Vec3 n = getnorm();
        float area = getarea();
        n.normalize();

        Vec3 s0 = gets0();
        Vec3 s1 = gets1();
        Vec3 s2 = gets2();

        float u,v;

        Vec3 s0s1 = s1 - s0;
        Vec3 s0s2 = s2 - s0;
        Vec3 N = Vec3::cross(s0s1,s0s2);
        float denom = Vec3::dot(N,N); 
        
        Vec3 inter;

        double D = Vec3::dot(s0, n);

        // 1) check that the ray is not parallel to the triangle:

         if ((Vec3::dot(n, d) <= 0.01) && (Vec3::dot(n, d) >= -0.01))
        {
            intersection.intersectionExists = false;
            return intersection;
        }

        // 2) check that the triangle is "in front of" the ray:
        double t = (D - Vec3::dot(o, n)) / Vec3::dot(d, n); // Fonction de t
        double orientation = Vec3::dot(d,  n);
        if (t > 0 && orientation <= 0) // On est dans le plan et on est face au plan (pour ne pas render le dos des triangles inutilement)
        //if(t > 0)
        {
        // 3) check that the intersection point is inside the triangle:
        // CONVENTION: compute u,v such that p = w0*c0 + w1*c1 + w2*c2, check that 0 <= w0,w1,w2 <= 1

            inter = o + (t * d); // Coordonnée du point d'intersection
            Vec3 C;
            // edge 0 
            Vec3 edge0 = s1 - s0;
            Vec3 vp0 = inter - s0;
            C = Vec3::cross(edge0,vp0);
            if (Vec3::dot(N,C) < 0){ // Sort de l'aire
            intersection.intersectionExists = false;
            return intersection;
            }

            // edge 1
            Vec3 edge1 = s2 - s1;
            Vec3 vp1 = inter - s1;
            C = Vec3::cross(edge1,vp1);
            u = Vec3::dot(N,C);
            if (u < 0){ // Sort de l'aire
            intersection.intersectionExists = false;
            return intersection;
            }

            // edge 2
            Vec3 edge2 = s0 - s2;
            Vec3 vp2 = inter - s2;
            C = Vec3::cross(edge2,vp2);
            v = Vec3::dot(N,C);
            if ( v < 0){ // Sort de l'aire
            intersection.intersectionExists = false;
            return intersection;
            }
            // ON est bien interect 
            u /= denom;
            v /= denom;
            //std::cout<<" u :"<< u <<"  v : "<<v<<"  DENOM : " <<denom <<std::endl;

                intersection.intersectionExists = true;
                intersection.t = t;
                intersection.w0 = u;
                intersection.w1 = v;
                intersection.w2 = 1 - u - v;
                intersection.normal = n;
                intersection.intersection = inter;
                intersection.bounce_direction = d - 2 * intersection.normal *( Vec3::dot(d,intersection.normal));
                intersection.tIndex;
        }

        return intersection;
    }
};
#endif
