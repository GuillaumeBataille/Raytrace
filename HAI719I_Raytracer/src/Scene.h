#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"

#include <GL/glut.h>

#define ECHANTILLONAGE_LIGHT 20

enum LightType
{
    LightType_Spherical,
    LightType_Quad
};

struct Light
{
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}
};

struct RaySceneIntersection
{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    Vec3 normal;           // ajout normal
    Vec3 position;         // ajout position
    Vec3 bounce_direction; // Ajout bounce direction
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RaySceneIntersection() : intersectionExists(false), t(FLT_MAX) {}
};

class Scene
{
    std::vector<Mesh> meshes;
    std::vector<Sphere> spheres;
    std::vector<Square> squares;
    std::vector<Light> lights;

public:
    Scene()
    {
    }

    void draw()
    {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for (unsigned int It = 0; It < meshes.size(); ++It)
        {
            Mesh const &mesh = meshes[It];
            mesh.draw();
        }
        for (unsigned int It = 0; It < spheres.size(); ++It)
        {
            Sphere const &sphere = spheres[It];
            sphere.draw();
        }
        for (unsigned int It = 0; It < squares.size(); ++It)
        {
            Square const &square = squares[It];
            square.draw();
        }
    }

    Material getRayMaterial(RaySceneIntersection RSI)
    {
        switch (RSI.typeOfIntersectedObject)
        {
        case 0:
            return spheres[RSI.objectIndex].material;
            break;
        case 1:
            return squares[RSI.objectIndex].material;
            break;
        default:
            return Material();
            break;
        }
    }

    RaySceneIntersection computeIntersection(Ray const &ray, double znear)
    {
        RaySceneIntersection result;     // L'interesection la plus proche
        RaySphereIntersection raysphere; // Definit mon IntersectionSphere
        RaySquareIntersection raysquare; // Definit mon IntersectionSquare

        for (unsigned int i = 0; i < spheres.size(); i++)
        {                                          // Pour toute les sphères rencontrées;
            raysphere = spheres[i].intersect(ray); // On teste l'intersection
            if (raysphere.intersectionExists)
            { // Si l'intersection a bien lieu
                if (raysphere.t < result.t && raysphere.t != 0 && raysphere.t > znear)
                { // On regarde si on a déja eu une intersection, si c'est pas le cas

                    result.intersectionExists = true;                     // On le mets a true
                    result.objectIndex = i;                               // On recup l'id de l'élément touché
                    result.typeOfIntersectedObject = 0;                   // On set a 0 le type d'objet rencontré (avec 0 pour sphère et 1 pour square )
                    result.t = raysphere.t;                               // On recup son t, distance entre la caméra et le point d'intersection.
                    result.raySphereIntersection = raysphere;             // On récupère l'intersection
                    result.normal = raysphere.normal;                     // On recupère la normale
                    result.position = raysphere.intersection;             // On recupère la position de l'intersection
                    result.bounce_direction = raysphere.bounce_direction; // On recupère la direction reflechie
                }
            }
        }

        for (unsigned int j = 0; j < squares.size(); j++)
        {
            raysquare = squares[j].intersect(ray); // On teste l'intersection
            if (raysquare.intersectionExists)
            { // Si l'intersection a bien lieu
                if (raysquare.t < result.t && raysquare.t != 0 && raysquare.t > znear)
                { // On regarde si on a déja eu une intersection, si c'est pas le cas

                    result.intersectionExists = true;                     // On le mets a true
                    result.objectIndex = j;                               // On recup l'id de l'élément touché
                    result.typeOfIntersectedObject = 1;                   // On set a 1 le type d'objet rencontré (avec 0 pour sphère et 1 pour square )
                    result.t = raysquare.t;                               // On recup son t, distance entre la caméra et le point d'intersection.
                    result.raySquareIntersection = raysquare;             // On récupère l'intersection
                    result.normal = raysquare.normal;                     // On récupère la normale
                    result.position = raysquare.intersection;             // On recupère la position de l'intersection
                    result.bounce_direction = raysquare.bounce_direction; // On recupère la direction reflechie
                }
            }
        }
        return result;
    }

    void phong(RaySceneIntersection RSI, Vec3 &color, Material mat, Ray &ray)
    {
        for (unsigned long int i = 0; i < lights.size(); i++)
        {
            Vec3 L = lights[i].pos - RSI.position; // Le vecteur depuis le point d'intersection vers la lumière courante
            L.normalize();
            //DIFFUSE
            double lightvalue = Vec3::dot(L, RSI.normal); // On le recup
            if (lightvalue <= 0)
                continue;
            Vec3 diffuse = Vec3::compProduct(lights[i].material, mat.diffuse_material) * lightvalue;
            //SPECULAR
            Vec3 R = 2 * (Vec3::dot(RSI.normal, L) * RSI.normal) - L;
            Vec3 V = ray.origin() - RSI.position;
            R.normalize();
            V.normalize();
            Vec3 specular = Vec3::compProduct(lights[i].material, mat.specular_material * pow(Vec3::dot(R, V), mat.shininess));
            color = diffuse + specular;
        }
        color += mat.ambient_material;
    }

    void shadow(RaySceneIntersection RSI, Vec3 &color, Material mat)
    {
        for (unsigned long int i = 0; i < lights.size(); i++)
        {
            Vec3 L = lights[i].pos - RSI.position; // Le vecteur depuis le point d'intersection vers la lumière courante
            double range = L.length();             // La distance entre la position et la lumière
            L.normalize();
            double shadowvalue = 0;           // Valeur de l'ombre courante (O si obstacle et 1 sinon)
            Ray light = Ray(RSI.position, L); // origine la position de l'intersect et direction la light

            // OMBRES DURES - cas avec type of light spherical / ponctuel
            if (lights[i].type == LightType_Spherical)
            {
                RaySceneIntersection inter2 = computeIntersection(light, 0.0001);
                if (inter2.t >= range)
                {
                    shadowvalue = 1; // Le 0.001 sur le znear me sert a ne pas considérer comme obstacle un vertice voisin ou moi même
                    // Pas besoin de shadowvalue car ici l'ombre est nette, soit 0 soit 1;
                }
            }
            else
            {
                // OMBRES DOUCES
                // création du square de lumière a partir du mesh ajouté dans l'init de la scène de cornell
                Vec3 bottom_left = lights[i].quad.vertices[0].position;
                Vec3 right = lights[i].quad.vertices[1].position - bottom_left;
                Vec3 up = lights[i].quad.vertices[3].position - bottom_left;

                for (unsigned long int j = 0; j < ECHANTILLONAGE_LIGHT; j++)
                {
                    Vec3 randx = (double(rand()) / RAND_MAX * right);
                    Vec3 randy = (double(rand()) / RAND_MAX * up);
                    Vec3 randpoint = randx + randy + bottom_left;

                    Vec3 vec_vers_randpoint = randpoint - RSI.position;
                    double rangebis = vec_vers_randpoint.norm();
                    vec_vers_randpoint.normalize();
                    light = Ray(RSI.position, vec_vers_randpoint);
                    RaySceneIntersection inter3 = computeIntersection(light, 0.0001);
                    if (!(inter3.intersectionExists && inter3.t < range && inter3.t > 0.0001))
                    {
                        shadowvalue++;
                    }
                }
                shadowvalue /= ECHANTILLONAGE_LIGHT;  // Recupère le pourcentage de light qui est passé en divisant par le nombre d'échantillons
                color *= shadowvalue / lights.size(); // On pondère l'intensité lumineuse (ambiante + diffuse + speculaire) par l'ombrage
            }
        }
    }

    //Fonctionne pas..
    Vec3 refract(const Vec3 &I, const Vec3 &N, const float &ior)
    {
        float NdotI = Vec3::dot(N, I);
        float cosi;
        NdotI <= -1 ? cosi = -1 : (NdotI >= 1 ? cosi = 1 : cosi = NdotI);
        float etai = 1, etat = ior;
        Vec3 n = N;
        if (cosi < 0)
        {
            cosi = -cosi;
        }
        else
        {
            std::swap(etai, etat);
            n = N * -1;
        }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);

        if (k >= 0)
            return (eta * I + (eta * cosi - sqrtf(k)) * n);
        else
            return Vec3(0, 0, 0);
    }

    Vec3 rayTraceRecursive(Ray ray, int NRemainingBounces, double znear)
    {

        // TODO appeler la fonction recursive
        Vec3 color;
        int nb_current_bounce = 0;
        RaySceneIntersection RSI = computeIntersection(ray, znear);
        if (!RSI.intersectionExists) // Si y'a pas d'intersection
            return Vec3(0, 0, 0);

        Material mat = getRayMaterial(RSI); // On recupère le material courant
        color = mat.ambient_material;       // Récup couleur ambiant
        phong(RSI, color, mat, ray);        // On appelle phong pour y ajouter la couleur de l'intersection (diffuse + specular) modulé par l'ombre
        shadow(RSI, color, mat);            // On appelle shadow qui va compute les ompbres dans color

        // Recursive part :
        //      Reflection sur mirroir
        if (mat.type == Material_Mirror && NRemainingBounces > 0)
        { //On teste si on est avec mirroir ET si il reste des rebonds
            // Parametrage du nouveau point de depart du rayon
            Vec3 new_origin = RSI.position;            //La position de l'intersection devient l'origine du rayon
            Vec3 new_direction = RSI.bounce_direction; // La direction bounce (reflexion via old direction et normale) devient la direction du rayon
            Ray Ray_bounce = Ray(new_origin, new_direction);
            color += rayTraceRecursive(Ray_bounce, NRemainingBounces - 1, 0.001) / NRemainingBounces; // On retire un rayon en décrémentant le nbr de rayon a tirer
        }
        //      Refraction dans glass
        /* if (mat.type == Material_Glass && NRemainingBounces > 0)
        { //On teste si on est avec mirroir ET si il reste des rebonds
            // Parametrage du nouveau point de depart du rayon
            //std::cout << "Hello" << std::endl;
            Vec3 N = RSI.normal;
            Vec3 V = ray.direction();
            double n = mat.transparency;
            double NdotV = Vec3::dot(N, V);
            double W = (1 - n * n * (1 - (NdotV * NdotV)));
            if (W < 0)
            {
                return color;
            }
            Vec3 T = (n * NdotV - sqrt(1 - n * n * (1 - (NdotV * NdotV)))) * N - n * V;
            // T = refract(V,N,n);
            Vec3 new_origin_legerement_decal = RSI.position + T * 0.01; //La position de l'intersection devient l'origine du rayon
            //Ray Ray_Refract = Ray(RSI.position, T);
            //color += rayTraceRecursive(Ray_Refract,NRemainingBounces-1, 0.001); // On retire un rayon en décrémentant le nbr de rayon a tirer
        }*/
        //std::cout<<"square size : "<<squares.size()<<std::endl;
        //TODO RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        // znear =0.001 c'est pour pas qu'il se cogne sur lui même ou sur son voisin

        return color;
    }

    Vec3 rayTrace(Ray const &rayStart, double znear)
    {
        //  std::cout << "call a ray"<<std::endl;
        // TODO appeler la fonction recursive
        Vec3 color;
        RaySceneIntersection RSI = computeIntersection(rayStart, znear);
        //      std::cout <<"Sphere t : "<< RSI.raySphereIntersection.t <<"   Square t :" <<RSI.raySquareIntersection.t<<std::endl;;
        if (RSI.intersectionExists)
        { // Si il y a bien une intersection avec un objet de la scène
            if (RSI.typeOfIntersectedObject == 0)
            {                                                               // On regarde si on est avec une sphère (0)
                color = spheres[RSI.objectIndex].material.diffuse_material; // Récup couleur
            }
            else if (RSI.typeOfIntersectedObject == 1)
            {                                                               // On regarde si on est avec une square (0)
                color = squares[RSI.objectIndex].material.diffuse_material; // Récup couleur
            }
        }

        return color;
    }

    void setup_single_sphere()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0., 0., 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1, 0, 0);
            s.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s.material.shininess = 20;

            meshes.resize(meshes.size() + 1);
            meshes[meshes.size() - 1] = s;
        }
    }

    void setup_single_square()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        {
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0., 0., 1.);
            s.material.specular_material = Vec3(0.8, 0.8, 0.8);
            s.material.shininess = 20;
        }
    }

    void setup_cornell_box()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(0.0, 1.5, 0.0);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
            // Adding quad mesh
            light.type = LightType_Quad;
            //light.type = LightType_Spherical; 
            Square s;
            s.setQuad(Vec3(-1., 0, 0.0), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(-1 * light.pos);
            //s.translate(Vec3(0.,0.5,0.0));
            //s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            //s.build_arrays();
            light.quad = s;
        }

        /*{
    squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., 0, 0.0), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            
            s.translate(-1 * x);
            //s.translate(Vec3(0.,0.5,0.0));
            //s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();

}*/
        { // Back Wall

            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 0., 0.7);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
        }

        { // Left Wall

            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
        }

        { // Right Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.0, 1.0, 0.0);
            s.material.specular_material = Vec3(0.0, 1.0, 0.0);
            s.material.shininess = 16;
        }

        { // Floor
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(0.3, 0.43, 0.5);
            s.material.specular_material = Vec3(1., 1.0, 1.0);
        }

        { // Ceiling
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1, 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.1, 1, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
        }

        { // Front Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
        }

        { // GLASS Sphere

            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }

        { // MIRRORED Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(0.4, 0.1, 0.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }
    }
};

#endif
