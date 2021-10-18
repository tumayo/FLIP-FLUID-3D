#ifndef PARTICLES_H
#define PARTICLES_H

#include <iostream>
#include "array3.h"
#include "vec.h"
#include "array3_utils.h"

using namespace std;

//FLIP particles
class Particles {

public:

	int num_particles; // number of particles
	Array3f x, y, u, v; // positions and velocities
	Array3f gamma; // vorticity of each particle
	Array3f sum;   // transfer stuff

	std::vector<Vec3f> particles_pos;
	std::vector<Vec3f> particles_vel;

	void initialize(float width, int ni, int nj, int nk) {
		num_particles = ni * nj * nk; // each grid has one particle
		float dx = width / ni; // dx = width / nj

		int rando = 21753;
		int part_per_cell = 6;
		for (int i = 0; i < ni; ++i) for (int j = 0; j < nj; ++j) for (int k = 0; k < nj; ++k) {
			for (int p = 0; p < part_per_cell; ++p) {
				float r1 = randhashf(++rando) - 0.5;
				float r2 = randhashf(++rando) - 0.5;
				float r3 = randhashf(++rando) - 0.5;
				Vec3f pos = Vec3f((i + 0.5f) * dx, (j + 0.5f) * dx, (k + 0.5f) * dx) + Vec3f(r1, r2, r3) * dx;
				particles_pos.push_back(pos);
				particles_vel.push_back(Vec3f(0.0f, 0.0f, 0.0f)); // Initially particles with 0 velocity
			}
		}
	}

	void transfer_to_grid(float width, int ni, int nj, int nk, Array3f& grid_u, Array3f& grid_v, Array3f& grid_w) {
		int i, j, k;
		float fx, fy, fz;
		float dx = width / (float)ni;

		sum.resize(ni + 1, nj, nk);
		sum.set_zero();

		grid_u.set_zero();
		for (int p = 0; p < particles_pos.size(); ++p) {
			get_barycentric(particles_pos[p][0] / dx, i, fx, 0, ni + 1);
			get_barycentric(particles_pos[p][1] / dx - 0.5f, j, fy, 0, nj);
			get_barycentric(particles_pos[p][2] / dx - 0.5f, k, fz, 0, nk);
			accumulate(grid_u, particles_vel[p][0], i, j, k, fx, fy, fz);
		}

		for (i = 0; i < ni + 1; ++i) for (j = 0; j < nj; ++j) for (k = 0; k < nk; ++k) {
			if (sum(i, j, k) != 0)
				grid_u(i, j, k) /= sum(i, j, k);
			else
				grid_u(i, j, k) = 0;
		}

		sum.resize(ni, nj + 1, nk);
		sum.set_zero();
		grid_v.set_zero();
		for (int p = 0; p < particles_pos.size(); ++p) {
			get_barycentric(particles_pos[p][0] / dx - 0.5f, i, fx, 0, ni);
			get_barycentric(particles_pos[p][1] / dx, j, fy, 0, nj + 1);
			get_barycentric(particles_pos[p][2] / dx - 0.5f, k, fz, 0, nk);
			accumulate(grid_v, particles_vel[p][1], i, j, k, fx, fy, fz);

		}
		for (i = 0; i < ni; ++i) for (j = 0; j < nj + 1; ++j) for (k = 0; k < nk; ++k) {
				if (sum(i, j, k) != 0)
					grid_v(i, j, k) /= sum(i, j, k);
				else
					grid_v(i, j, k) = 0;
		}

		sum.resize(ni, nj, nk + 1);
		sum.set_zero();
		grid_w.set_zero();
		for (int p = 0; p < particles_pos.size(); ++p) {
			get_barycentric(particles_pos[p][0] / dx - 0.5f, i, fx, 0, ni);
			get_barycentric(particles_pos[p][1] / dx - 0.5f, j, fy, 0, nj);
			get_barycentric(particles_pos[p][2] / dx, k, fz, 0, nk + 1);
			accumulate(grid_w, particles_vel[p][2], i, j, k, fx, fy, fz);

		}
		for (i = 0; i < ni; ++i) for (j = 0; j < nj; ++j) for (k = 0; k < nk + 1; ++k) {
			if (sum(i, j, k) != 0)
				grid_w(i, j, k) /= sum(i, j, k);
			else
				grid_w(i, j, k) = 0;
		}
	}


	void move_particles_in_grid(float dt, Array3f grid_u, Array3f grid_v, Array3f grid_w, float xmax, float ymax, float zmax, float dx, int ni, int nj, int nk) {
		Vec3f midx, gu;
		for (int p = 0; p < particles_pos.size(); ++p) {

			// first stage of Runge-Kutta 2 (do a half Euler step)
			float x_pos = particles_pos[p][0] / dx;
			float y_pos = particles_pos[p][1] / dx - 0.5;
			float z_pos = particles_pos[p][2] / dx - 0.5;
			gu[0] = interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_u);

			x_pos = particles_pos[p][0] / dx - 0.5;
			y_pos = particles_pos[p][1] / dx;
			z_pos = particles_pos[p][2] / dx - 0.5;
			gu[1] = interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_v);

			x_pos = particles_pos[p][0] / dx - 0.5;
			y_pos = particles_pos[p][1] / dx - 0.5;
			z_pos = particles_pos[p][2] / dx;
			gu[2] = interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_w);

			midx[0] = particles_pos[p][0] + 0.5 * dt * gu[0];
			midx[1] = particles_pos[p][1] + 0.5 * dt * gu[1];
			midx[2] = particles_pos[p][2] + 0.5 * dt * gu[2];

			midx[0] = clamp(midx[0], 0.0f, xmax);
			midx[1] = clamp(midx[1], 0.0f, ymax);
			midx[2] = clamp(midx[2], 0.0f, zmax);

			// second stage of Runge-Kutta 2
			x_pos = midx[0] / dx;
			y_pos = midx[1] / dx - 0.5;
			z_pos = midx[2] / dx - 0.5;
			gu[0] = interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_u);

			x_pos = midx[0] / dx - 0.5;
			y_pos = midx[1] / dx;
			z_pos = midx[2] / dx - 0.5;
			gu[1] = interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_v);

			x_pos = midx[0] / dx - 0.5;
			y_pos = midx[1] / dx - 0.5;
			z_pos = midx[2] / dx;
			gu[2] = interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_w);

			particles_pos[p][0] += dt * gu[0];
			particles_pos[p][1] += dt * gu[1];
			particles_pos[p][2] += dt * gu[2];

			particles_pos[p][0] = clamp(particles_pos[p][0], 0.0f, xmax);
			particles_pos[p][1] = clamp(particles_pos[p][1], 0.0f, ymax);
			particles_pos[p][2] = clamp(particles_pos[p][2], 0.0f, zmax);

		}
	}


	void update_from_grid(float dx, float ni, float nj, float nk, const Array3f& grid_u, const Array3f& grid_v, const Array3f& grid_w,
		const Array3f& grid_du, const Array3f& grid_dv, const Array3f& grid_dw) {

		for (int p = 0; p < particles_pos.size(); ++p) {

			float x_pos = particles_pos[p][0] / dx;
			float y_pos = particles_pos[p][1] / dx - 0.5;
			float z_pos = particles_pos[p][2] / dx - 0.5;
			particles_vel[p][0] += interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_du); // FLIP

			particles_vel[p][0] = 0.99 * particles_vel[p][0] + 0.01 * interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_u); //PIC/FLIP mix
			//particles_vel[p][0] = interpolate_value(Vec2f(x_pos, y_pos), grid_u); //PIC

			x_pos = particles_pos[p][0] / dx - 0.5;
			y_pos = particles_pos[p][1] / dx;
			z_pos = particles_pos[p][2] / dx - 0.5;
			particles_vel[p][1] += interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_dv);

			particles_vel[p][1] = 0.99 * particles_vel[p][1] + 0.01 * interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_v); //PIC/FLIP mix
			//particles_vel[p][1] = interpolate_value(Vec2f(x_pos, y_pos), grid_v); //PIC

			x_pos = particles_pos[p][0] / dx - 0.5;
			y_pos = particles_pos[p][1] / dx - 0.5;
			z_pos = particles_pos[p][2] / dx;
			particles_vel[p][2] += interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_dw);

			particles_vel[p][2] = 0.99 * particles_vel[p][2] + 0.01 * interpolate_value(Vec3f(x_pos, y_pos, z_pos), grid_w); //PIC/FLIP mix
		}
	}

private:
	void accumulate(Array3f& accum, float q, int i, int j, int k, float fx, float fy, float fz) {

		float weight;
		weight = (1 - fx) * (1 - fy) * (1 - fz);
		accum(i, j, k) += weight * q;
		sum(i, j, k) += weight;

		weight = fx * (1 - fy) * (1 - fz);
		accum(i + 1, j, k) += weight * q;
		sum(i + 1, j, k) += weight;

		weight = (1 - fx) * fy * (1 - fz);
		accum(i, j + 1, k) += weight * q;
		sum(i, j + 1, k) += weight;

		weight = (1 - fx) * (1 - fy) * fz;
		accum(i, j, k + 1) += weight * q;
		sum(i, j, k + 1) += weight;

		weight = fx * fy * fz;
		accum(i + 1, j + 1, k + 1) += weight * q;
		sum(i + 1, j + 1, k + 1) += weight;

	}

};

#endif
