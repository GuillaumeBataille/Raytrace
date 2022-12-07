#include <vector>
#include <string>
#include "Vec3.h"
#include "Triangle.h"

class Box{
    //Attributs
    public: 
    std::vector<Vec3> vertices;
    std::vector<Vec3> triangles;
    Vec3 BB_min;
    Vec3 BB_max;
    //Constructor : 
    Box(Vec3 BB_min, Vec3 BB_max){
        this->BB_min = BB_min;
        this->BB_max = BB_max;
        buildbox();
    }
    Box(){}

       /// @brief 
       void buildbox() {
        vertices.resize(8);
        triangles.resize(12);
        
        vertices[0]=BB_min;
        vertices[1]=Vec3(BB_min[0]  ,BB_min[1]  ,BB_max[2]  );
        vertices[2]=Vec3(BB_min[0]  ,BB_max[1]  ,BB_min[2]  );
        vertices[3]=Vec3(BB_min[0]  ,BB_max[1]  ,BB_max[2]  );
        vertices[4]=Vec3(BB_max[0]  ,BB_min[1]  ,BB_min[2]  );
        vertices[5]=Vec3(BB_max[0]  ,BB_min[1]  ,BB_max[2]  );
        vertices[6]=Vec3(BB_max[0]  ,BB_max[1]  ,BB_min[2]  );
        vertices[7]=BB_max;

        triangles[0]=Vec3(0,1,2 );// Good
        triangles[1]=Vec3(3,2,1 );

        triangles[2]=Vec3(0,4,1 );// Good
        triangles[3]=Vec3(5,1,4 );

        triangles[4]=Vec3(0,4,2 );// Good
        triangles[5]=Vec3(6,2,4 );

        triangles[6]=Vec3(7,5,6 );// Good
        triangles[7]=Vec3(4,6,5 );

        triangles[8]=Vec3(7,3,6 );// Good
        triangles[9]=Vec3(2,6,3 );

        triangles[10]=Vec3(7,5,3 );// Good
        triangles[11]=Vec3(3,1,5 );
    }

      RayTriangleIntersection intersect( Ray const & ray) {
        RayTriangleIntersection result;
        result.intersectionExists = false;
        for (Vec3 t : triangles) // pour chaque triangle de ma box
        {
            Vec3 s0 = vertices[t[0]];
            Vec3 s1 = vertices[t[1]];
            Vec3 s2 = vertices[t[2]];
            Triangle current_triangle = Triangle(s0,s1,s2);
            RayTriangleIntersection current_intersect = current_triangle.getIntersection(ray);

            if (current_intersect.intersectionExists){
                result = current_intersect;
            }
        }
        return result;
      }
};