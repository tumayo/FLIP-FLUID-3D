#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include <stdlib.h>     /* srand, rand */
#include <time.h>    

#include "array3_utils.h"
#include "fluidsim.h"


using namespace std;

int grid_resolution = 20;
float timestep = 0.01f;
int frame = 0;

float grid_width = 1;

FluidSim sim;

float sphere_phi(const Vec3f& position, const Vec3f& centre, float radius) {
    return (dist(position, centre) - radius);
}

float box_phi(const Vec3f& position, const Vec3f& centre, Vec3f& b) {
    Vec3f dif = position - centre;
    Vec3f d = Vec3f(abs(dif[0]), abs(dif[1]), abs(dif[2])) - b;
    return dist(Vec3f(max(d[0], 0.0f), max(d[1], 0.0f), max(d[2], 0.0f)), Vec3f(0, 0, 0)) + min(max(d[0], max(d[1], d[2])), 0.0f);
}

Vec3f c0(0.5f, 0.5f, 0.5f); // centre
float rad0 = 0.4f; // radius
Vec3f b0(0.4f, 0.4f, 0.4f); // box size


float boundary_phi(const Vec3f& position) {
    //Sphere
    //return -sphere_phi(position, c0, rad0);
    //Box
    return -box_phi(position, c0, b0);
}


float liquid_phi(const Vec3f& position) {
    return sphere_phi(position, Vec3f(0.55f, 0.55f, 0.4f), 0.23f);
}

void export_particles(string path, int frame, const std::vector<Vec3f>& particles, float radius);
// Tumay -- Export particles in .ply format to import Houdini Geometry
void export_particles_houdini(string path, int frame, const std::vector<Vec3f>& particles, float radius);

//Main testing code
//-------------
int main(int argc, char** argv)
{
    string outpath("");

    printf("Initializing data\n");
    sim.initialize(grid_width, grid_resolution, grid_resolution, grid_resolution);

    printf("Initializing boundary\n");
    sim.set_boundary(boundary_phi);

    //printf("Initializing liquid\n");
    //sim.set_liquid(liquid_phi);

    printf("Exporting initial data\n");
    export_particles(outpath, 0, sim.particles, sim.particle_radius);

    for (frame = 1; frame < 1000; ++frame) {
        printf("--------------------\nFrame %d\n", frame);

        //Simulate
        printf("Simulating fluid\n");
        //sim.advance(timestep);
        sim.flip_adv_advance(timestep); // FLIP advection
   
        //Initialize random seed
        srand(time(NULL));

        //Emit particles randomly within the region (7-13)x(7-13)
        for (int k = 7; k < 13; k++) for (int i = 7; i < 13; i++)  {
            int j = 3;
            float rand_i = (((float)rand() / RAND_MAX) * 6 + 7);
            float rand_k = (((float)rand() / RAND_MAX) * 6 + 7);
            Vec3f pos( rand_i * sim.dx, j * sim.dx, rand_k * sim.dx);
            sim.density(roundf(rand_i), j, roundf(rand_k)) = 1;
            sim.temperature(roundf(rand_i), j, roundf(rand_k)) = 300;
            sim.particles.push_back(pos);
        }

        printf("Exporting particle data\n");
        export_particles(outpath, frame, sim.particles, sim.particle_radius);
    }

    return 0;
}

void export_particles(string path, int frame, const std::vector<Vec3f>& particles, float radius) {
    //Write the output

    std::stringstream strout;
    strout << path << "particles_" << frame << ".txt";
    string filepath = strout.str();

    ofstream outfile(filepath.c_str());
    //write vertex count and particle radius
    outfile << particles.size() << " " << radius << std::endl;
    //write vertices
    for (unsigned int i = 0; i < particles.size(); ++i)
        outfile << particles[i][0] << " " << particles[i][1] << " " << particles[i][2] << std::endl;

    std::cout << "Writing to: " << filepath << std::endl;
    outfile.close();
}

//Export particles in .ply format to import Houdini Geometry
void export_particles_houdini(string path, int frame, const std::vector<Vec3f>& particles, float radius) {
    //Write the output
    std::stringstream strout;
    strout << path << "particles_" << frame << "SLICE" << ".ply";
    string filepath = strout.str();
    ofstream outfile(filepath.c_str());

    //Header
    outfile << "ply" << endl;
    outfile << "format ascii 1.0" << endl;
    outfile << "element vertex " << particles.size() << endl;
    outfile << "property float x" << endl;
    outfile << "property float y" << endl;
    outfile << "property float z" << endl;
    outfile << "end_header" << endl;

    //Write vertices
    for (unsigned int i = 0; i < particles.size(); ++i)
        outfile << particles[i][0] << " " << particles[i][1] << " " << particles[i][2] << std::endl;

    std::cout << "Writing to: " << filepath << std::endl;
    outfile.close();
}


