//
// Created by smutekj on 26.04.23.
//

#include "VerletList.h"
#include "Grid.h"
#include "Systemx.h"
#include "NeighbourSearcherImplementation.h"

/*

std::vector<real> calcDr2s(const size_t n_changed,
                           real rx_coords_new[3][700], const gmx::RVec ri,
                           const size_t n_atoms, const real box_length) {

    __m256 box_lengths = _mm256_set1_ps(box_length);

    std::vector<real> dr2s(n_changed);


    size_t n_cycles = n_changed/8;
    size_t n_rest = (n_changed-1) - n_cycles*8;

    __m256 x_coords_o_s = _mm256_set1_ps(ri[0]);
    __m256 y_coords_o_s = _mm256_set1_ps(ri[1]);
    __m256 z_coords_o_s = _mm256_set1_ps(ri[2]);

    for (size_t pivot = 0; pivot < n_cycles; ++pivot) {
        __m256 x_coords_s = _mm256_load_ps(&(rx_coords_new[0][pivot*8]));
        __m256 y_coords_s = _mm256_load_ps(&(rx_coords_new[1][pivot*8]));
        __m256 z_coords_s = _mm256_load_ps(&(rx_coords_new[2][pivot*8]));
        __m256 dxs = simd_dist_sq_1d_PBC(x_coords_o_s, x_coords_s, box_lengths);
        __m256 dys = simd_dist_sq_1d_PBC(y_coords_o_s, y_coords_s, box_lengths);
        __m256 dzs = simd_dist_sq_1d_PBC(z_coords_o_s, z_coords_s, box_lengths);
        __m256 wtf = simd_dist_sq_3d(dxs, dys, dzs);

        _mm256_storeu_ps(&(dr2s[pivot*8]), wtf);
    }
    for (size_t rest = 8*n_cycles; rest < 8*n_cycles + n_rest; ++rest) {
        const rvec pica = {rx_coords_new[0][rest], rx_coords_new[1][rest], rx_coords_new[2][rest]};
        dr2s[rest] = dist_sq_pbc(ri, pica, box_length/2.0f, box_length);
    }

    return dr2s;
}
*/


/*

void filterFarAwayParticles(
        const real* __restrict__ x_coords, const real* __restrict__ y_coords, const real* __restrict__ z_coords,
        const size_t* __restrict__ gathered_indices, std::vector<real>& dr2s, const size_t particle_index,
        const real x_new, const real y_new, const  real z_new,
        const size_t n_particle_indices,
        const real box_length, const real half_box_length, const real cutoff_sq) {

    size_t first_index = 0;
    size_t last_index = n_particle_indices - 1;
#pragma clang loop vectorize(enable)
    for (int i = 0; i < n_particle_indices; ++i) {

        const auto dx = distancePBC(x_coords[gathered_indices[i]], x_new,
                                    half_box_length, box_length) + box_length * (gathered_indices[i] == particle_index);
        const auto dy = distancePBC(y_coords[gathered_indices[i]], y_new,
                                    half_box_length, box_length);
        const auto dz = distancePBC(z_coords[gathered_indices[i]], z_new,
                                    half_box_length, box_length);
        const auto dr2 = dx * dx + dy * dy + dz * dz;
        size_t index_to_insert_to = first_index * (dr2 < cutoff_sq) + last_index*(dr2 >= cutoff_sq);
        dr2s[index_to_insert_to] = dr2;
        first_index += (dr2 < cutoff_sq) ;
        last_index -= dr2 >= cutoff_sq;
    }

    dr2s.resize(first_index);
}
*/

void fillVerletList(const bool* __restrict__ is_within_buffer, const size_t* __restrict__  neighbouring_particle_indices,
                                const size_t n_particles_in_cell, std::vector<size_t>& particle_to_list,
                                size_t& last_index, size_t& first_index){

//#pragma loop distribute(enable) // Vectorizer complains about "unsafe dependent memory operations in loop"
    for (int j = 0; j < n_particles_in_cell; ++j) {
        const size_t index_to_insert_to = last_index*(!is_within_buffer[j]) + first_index*(is_within_buffer[j]);
        particle_to_list[index_to_insert_to] = neighbouring_particle_indices[j];
        first_index += is_within_buffer[j];
        last_index -= !is_within_buffer[j];
    }
//    assert(first_index <= last_index);
}

void fillVerletListSIMD(const bool* __restrict__ is_within_buffer, const size_t* __restrict__  neighbouring_particle_indices,
                    const size_t n_particles_in_cell, size_t* __restrict__ particle_to_list, real* dr2s,
                    size_t& last_index, size_t& first_index){

    for (int j = 0; j < n_particles_in_cell; ++j) {
        const size_t index_to_insert_to = last_index*(!is_within_buffer[j]) + first_index*(is_within_buffer[j]);

        first_index += is_within_buffer[j];
        last_index -= !is_within_buffer[j];
        particle_to_list[index_to_insert_to] = neighbouring_particle_indices[j];
    }
}

void VerletListSpecial::fillVerletListSpecial(const bool* __restrict__ is_within_buffer, const size_t* __restrict__  neighbouring_particle_indices,
                                              gmx::RVec* __restrict__ r_coords,
                                                const size_t n_particles_in_cell, const size_t particle_index,
                                                size_t& last_index, size_t& first_index){

//#pragma loop distribute(enable) // Vectorizer complains about "unsafe dependent memory operations in loop"
    for (int j = 0; j < n_particles_in_cell; ++j) {
        const size_t index_to_insert_to = last_index*(!is_within_buffer[j]) + first_index*(is_within_buffer[j]);

        first_index += is_within_buffer[j];
        last_index -= !is_within_buffer[j];
        particle_to_list[particle_index][index_to_insert_to] = neighbouring_particle_indices[j];
        particle_to_r[particle_index][index_to_insert_to] = r_coords[neighbouring_particle_indices[j]];
    }
}

void VerletList::constructVerletListByMolecule(const Systemx& system,
                                              const real box_length) {
    const auto half_box_length = box_length / 2.0;
//    const auto a = system.
    const auto& active_molecule_indices = system.getAllActiveMolecules();
    const auto& cell_vectors = system.cell_vectors2_;

    for (auto acitve_mol_index : active_molecule_indices) {
        const auto molecule = system.world.all_molecules.at(acitve_mol_index);
        for (auto world_index = molecule.first_; world_index <= molecule.last_; ++world_index) {

            const auto cell_index = system.world_to_cell[world_index].first;
            const auto index_in_cell = system.world_to_cell[world_index].second;
            const auto r = cell_vectors[cell_index].r[index_in_cell].r;

            const auto &neighbouring_cells = system.grid->cell2nearest_neighbours.at(cell_index);
            size_t total_neighbouring_particle_indices(0);
            particle_to_list[world_index].clear();

            for (const auto cell_index_n: neighbouring_cells) {
                total_neighbouring_particle_indices += system.cell_vectors2_.at(cell_index_n).atoms_in_cell;
            }
            bool is_within_buffer[total_neighbouring_particle_indices];

            particle_to_list[world_index].resize(total_neighbouring_particle_indices);
            size_t last_index = total_neighbouring_particle_indices - 1;
            size_t first_index = 0;

            for (const auto cell_index_n: neighbouring_cells) {
                const auto r_cell_n =system.cell_vectors2_.at(cell_index_n).r ;
                const auto n_particles_in_cell = system.cell_vectors2_.at(cell_index_n).atoms_in_cell;

#pragma clang loop vectorize(enable)
                for (int j = 0; j < n_particles_in_cell; ++j) {
                    const auto neighbouring_index = system.cell_to_world[cell_index_n][j];
                    const auto dx = distancePBC(r_cell_n[j].r[XX], r[XX],
                                                half_box_length, box_length);
                    const auto dy = distancePBC(r_cell_n[j].r[YY], r[YY],
                                                half_box_length, box_length);
                    const auto dz = distancePBC(r_cell_n[j].r[ZZ], r[ZZ],
                                                half_box_length, box_length);
                    const auto dr2 = dx * dx + dy * dy + dz * dz;
                    is_within_buffer[j] = dr2 < r_buffersq and world_index != neighbouring_index;

                }

                fillVerletList(is_within_buffer, system.cell_to_world[cell_index_n].data(),
                               n_particles_in_cell, particle_to_list[world_index],
                               last_index, first_index);
            }
            particle_to_list[world_index].resize(first_index);
        }
    }
}


void VerletList::constructVerletListByCell(const Systemx& system, const real box_length) {
    const auto half_box_length = (box_length / 2.0f);

    const auto& cell_list = system.cell_vectors2_;
    /// Iterate over cells of grid;
    for( int cell_index = 0; cell_index < system.cell_to_world.size(); ++cell_index) {

        const size_t n_particles_in_cell = cell_list[cell_index].atoms_in_cell;
        const auto& r_cell = cell_list[cell_index].r;

        const auto &neighbouring_cells = system.grid->cell2nearest_neighbours[cell_index];
        size_t total_neighbouring_particle_indices(0);
        for (const auto cell_index_n: neighbouring_cells) {
            total_neighbouring_particle_indices += cell_list[cell_index_n].atoms_in_cell;
        }

        for (size_t j = 0; j < n_particles_in_cell; ++j) {
            const auto iw =  system.cell_to_world.at(cell_index).at(j);
            assert(iw < system.world_to_cell.size());
            const auto r = r_cell[j].r;

            particle_to_list[iw].resize(total_neighbouring_particle_indices);

            size_t last_index  = total_neighbouring_particle_indices-1;
            size_t first_index = 0;
            size_t next_index = 0;

            for (const auto cell_index_n: neighbouring_cells) {

                const auto& neighbouring_particle_indices = system.cell_to_world[cell_index_n];
                const auto& r_cell_n = cell_list[cell_index_n].r;
                const auto n_particles_in_cell_n = cell_list[cell_index_n].atoms_in_cell;
                bool is_within_buffer[n_particles_in_cell_n];

                assert(n_particles_in_cell_n == neighbouring_particle_indices.size());

#pragma clang loop vectorize(enable)
                for (int jn = 0; jn < n_particles_in_cell_n; ++jn) {

                    const auto dx = distancePBC(r_cell_n[jn].r[XX], r[XX],
                                                half_box_length, box_length);
                    const auto dy = distancePBC(r_cell_n[jn].r[YY], r[YY],
                                                half_box_length, box_length);
                    const auto dz = distancePBC(r_cell_n[jn].r[ZZ], r[ZZ],
                                                half_box_length, box_length);
                    const auto dr2 = dx * dx + dy * dy + dz * dz;
                    is_within_buffer[jn] = dr2 < r_buffersq and iw != neighbouring_particle_indices[jn];
                    next_index += (int)is_within_buffer[jn];

                }
                assert(last_index < particle_to_list.at(iw).size());

                fillVerletList(is_within_buffer, neighbouring_particle_indices.data(),
                               n_particles_in_cell_n, particle_to_list.at(iw),
                               last_index, first_index);

            }
//            std::cout << next_index << "\n";
            assert(first_index - last_index == 1);
            particle_to_list[iw].resize(first_index);
        }
    }
}



void VerletList::fillPotentialData(const Systemx& system, const Molecule& molecule, const std::vector<Particle> &particles,
                                   std::vector<PotentialData>& pot_data,
                                   const gmx::RVec box_length) const {
    const auto half_box_length = box_length / 2.0;
    const auto& cell_vectors = system.cell_vectors2_;
    std::vector<PotentialData> pot_data_tmp(0);

    for (auto world_index = molecule.first_; world_index <= molecule.last_; ++world_index) {

        const auto r_coord = particles.at(world_index - molecule.first_).r;
        const auto atom_type = particles.at(world_index - molecule.first_).atom_type;
        const auto cell_index = system.grid->coordToCell(r_coord);

        const auto neighbouring_cells = system.grid->cell2nearest_neighbours.at(cell_index);
        size_t total_neighbouring_particle_indices = 0;

        for (const auto cell_index_n: neighbouring_cells) {
            total_neighbouring_particle_indices += system.cell_vectors2_.at(cell_index_n).atoms_in_cell;
        }
        bool is_within_buffer[total_neighbouring_particle_indices];
        pot_data_tmp.resize( total_neighbouring_particle_indices);
        size_t next_index = pot_data.size();
//        pot_data.resize( pot_data.size() + total_neighbouring_particle_indices);


        for (const auto cell_index_n: neighbouring_cells) {
            const auto r_cell_n =system.cell_vectors2_.at(cell_index_n).r ;
            const auto n_particles_in_cell = cell_vectors.at(cell_index_n).atoms_in_cell;

#pragma clang loop vectorize(enable)
            for (int j = 0; j < n_particles_in_cell; ++j) {
                const auto neighbouring_index = system.cell_to_world[cell_index_n][j];
                const auto dx = distancePBC(r_cell_n[j].r[XX], r_coord[XX],
                                            half_box_length[XX], box_length[XX]);
                const auto dy = distancePBC(r_cell_n[j].r[YY], r_coord[YY],
                                            half_box_length[YY], box_length[YY]);
                const auto dz = distancePBC(r_cell_n[j].r[ZZ], r_coord[ZZ],
                                            half_box_length[ZZ], box_length[ZZ]);
                const auto dr2 = dx * dx + dy * dy + dz * dz;
                bool is_within_molecule = neighbouring_index >= molecule.first_ and neighbouring_index <= molecule.last_;
                is_within_buffer[j] = dr2 < r_buffersq and !is_within_molecule;
                pot_data_tmp[j].dr2 = dr2;
            }

            for (int j = 0; j < n_particles_in_cell; ++j) {
                if(is_within_buffer[j]){
                    pot_data.push_back({
                       pot_data_tmp[j].dr2,
                       system.lj_table.pair_coeffs[atom_type][r_cell_n[j].atom_type].A ,
                       system.lj_table.pair_coeffs[atom_type][r_cell_n[j].atom_type].B
                    });
                    next_index++;
                }
            }
        }
//        pot_data.resize(next_index);
    }
}

void VerletList::fillPotentialData(const Systemx& system, const Molecule& molecule,
                                   std::vector<PotentialData>& pot_data,
                                   const gmx::RVec box_length) const {
    const auto half_box_length = box_length / 2.0;
    const auto& cell_vectors = system.cell_vectors2_;

    std::vector<PotentialData> pot_data_tmp(0);

    for (auto world_index = molecule.first_; world_index <= molecule.last_; ++world_index) {

        const auto r_coord = system.world.getParticle(world_index).r;
        const auto atom_type =  system.world.getParticle(world_index).atom_type;
        const auto cell_index = system.grid->coordToCell(r_coord);

        const auto &neighbouring_cells = system.grid->cell2nearest_neighbours.at(cell_index);
        size_t total_neighbouring_particle_indices = 0;

        for (const auto cell_index_n: neighbouring_cells) {
            total_neighbouring_particle_indices += system.cell_vectors2_.at(cell_index_n).atoms_in_cell;
        }
        bool is_within_buffer[total_neighbouring_particle_indices];
        pot_data_tmp.resize( total_neighbouring_particle_indices);
        size_t next_index = pot_data.size();
//        pot_data.resize( pot_data.size() + total_neighbouring_particle_indices);

        for (const auto cell_index_n: neighbouring_cells) {
            const auto r_cell_n =cell_vectors.at(cell_index_n).r ;
            const auto n_particles_in_cell = cell_vectors.at(cell_index_n).atoms_in_cell;

#pragma clang loop vectorize(enable)
            for (int j = 0; j < n_particles_in_cell; ++j) {
                const auto neighbouring_index = system.cell_to_world[cell_index_n][j];
                const auto dx = distancePBC(r_cell_n[j].r[XX], r_coord[XX],
                                            half_box_length[XX], box_length[XX]);
                const auto dy = distancePBC(r_cell_n[j].r[YY], r_coord[YY],
                                            half_box_length[YY], box_length[YY]);
                const auto dz = distancePBC(r_cell_n[j].r[ZZ], r_coord[ZZ],
                                            half_box_length[ZZ], box_length[ZZ]);
                const auto dr2 = dx * dx + dy * dy + dz * dz;
                bool is_within_molecule = neighbouring_index >= molecule.first_ and neighbouring_index <= molecule.last_;
                is_within_buffer[j] = dr2 < r_buffersq and !is_within_molecule;
                pot_data_tmp[j].dr2 = dr2;
            }

            for (int j = 0; j < n_particles_in_cell; ++j) {
                if(is_within_buffer[j]){
                    pot_data.push_back({   pot_data_tmp[j].dr2,
                                           -system.lj_table.pair_coeffs[atom_type][r_cell_n[j].atom_type].A ,
                                           -system.lj_table.pair_coeffs[atom_type][r_cell_n[j].atom_type].B
                                       });
                    next_index++;
                }
            }
        }
//        pot_data.resize(next_index);
    }
}



