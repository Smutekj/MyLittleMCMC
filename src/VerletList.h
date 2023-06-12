//
// Created by smutekj on 26.04.23.
//

#ifndef GR_GCMC_CLION_VERLETLIST_H
#define GR_GCMC_CLION_VERLETLIST_H

class Grid;
class Systemx;
class Molecule;
class Particle;
template <typename Type>
    class CellVector;
#include <vector>
#include "vectypes.h"


struct PotentialData{
    real dr2;
    real A;
    real B;
};


class VerletList{
public:
    typedef unsigned long size_t;
    std::vector<std::vector<size_t>> particle_to_list;
    std::vector<std::vector<real>> particle_to_LJold;
    real r_cutoffsq;
    real r_buffersq;

    VerletList(size_t n_atoms, real r_cutoff, real r_buffer) :
            r_cutoffsq(r_cutoff*r_cutoff), r_buffersq(r_buffer*r_buffer)
    {
        particle_to_list.resize(n_atoms);
        particle_to_LJold.resize(n_atoms);
    }

    [[nodiscard]] const std::vector<size_t>& getNeighbours(const size_t particle_index) const{
        return particle_to_list[particle_index];
    }

    void fillPotentialData(const Systemx& system, const Molecule& molecule, const std::vector<Particle> &particles,
                                       std::vector<PotentialData>& pot_data,
                            gmx::RVec box_length) const;
    void fillPotentialData(const Systemx& system, const Molecule& molecule,
                           std::vector<PotentialData>& pot_data,
                            gmx::RVec box_length) const;
    void constructVerletListByMolecule(const Systemx& system,
                                       real box_length);
    void constructVerletListByCell(const Systemx& system,
                                      real box_length);
};


class VerletListSpecial{
public:
    typedef unsigned long size_t;
    std::vector<std::vector<size_t>> particle_to_list;
    std::vector<std::vector<gmx::RVec>> particle_to_r;

    std::vector<real> particle_to_LJold;
    real r_cutoffsq;
    real r_buffersq;

    VerletListSpecial(size_t n_atoms, real r_cutoff, real r_buffer) :
            r_cutoffsq(r_cutoff*r_cutoff), r_buffersq(r_buffer*r_buffer)
    {
        particle_to_list.resize(n_atoms);
        particle_to_LJold.resize(n_atoms);
    }

    [[nodiscard]] const std::vector<size_t>& getNeighbours(const size_t particle_index) const{
        return particle_to_list[particle_index];
    }


    void fillVerletListSpecial(const bool* __restrict__ is_within_buffer, const size_t* __restrict__  neighbouring_particle_indices,
                               gmx::RVec* r_coords,
                               size_t n_particles_in_cell, size_t particle_index,
                                                  size_t& last_index, size_t& first_index);


};



#endif //GR_GCMC_CLION_VERLETLIST_H


/*
for (size_t i = 0; i < n_changed; ++i) {
rx_coords_changed[0][i] = rx_coords[0][out[i]];
rx_coords_changed[1][i] = rx_coords[1][out[i]];
rx_coords_changed[2][i] = rx_coords[2][out[i]];
}

for (size_t i = 0; i < n_changed; ++i) {

size_t n_cycles = n_changed/8;
size_t n_rest = n_cycles*8 - (n_changed-1);

__m256 x_coords_o_s = _mm256_set1_ps(rx_coords[0][i]);
__m256 y_coords_o_s = _mm256_set1_ps(rx_coords[1][i]);
__m256 z_coords_o_s = _mm256_set1_ps(rx_coords[2][i]);

for (size_t pivot = 0; pivot < n_cycles; ++pivot) {
__m256 x_coords_s = _mm256_load_ps(&(rx_coords_changed[0][pivot*8]));
__m256 y_coords_s = _mm256_load_ps(&(rx_coords_changed[1][pivot*8]));
__m256 z_coords_s = _mm256_load_ps(&(rx_coords_changed[2][pivot*8]));
__m256 dxs = simd_dist_sq_1d_PBC(x_coords_o_s, x_coords_s, box_lengths);
__m256 dys = simd_dist_sq_1d_PBC(y_coords_o_s, y_coords_s, box_lengths);
__m256 dzs = simd_dist_sq_1d_PBC(z_coords_o_s, z_coords_s, box_lengths);
__m256 wtf = simd_dist_sq_3d(dxs, dys, dzs);

_mm256_storeu_ps(&(dr2s[pivot*8]), wtf);
//                _mm256_store_ps(&(dr2s[pivot*8]), wtf);
}
for (size_t rest = 8*n_cycles; rest < 8*n_cycles + n_rest; ++rest) {
const rvec pica = {rx_coords_changed[0][rest], rx_coords_changed[1][rest], rx_coords_changed[2][rest]};
dr2s[rest] = dist_sq_pbc(r_coords[i], pica, half_box_size, box_size);
}*/
