#include "Energy.h"
#include "Grid.h"
#include "VerletList.h"
#include "Systemx.h"


Energy::Energy(size_t n_atoms) {

    for (int i = 0; i < n_atoms; ++i){
        all_indices.push_back(i);
    }
}

real Energy::calcTotalEnergyBasic(const Systemx& system,
                                 const gmx::RVec box_length,
                                 real cutoff_sq)const{
    real energy = 0.0;
    const auto half_box_length = box_length/2.;

    const auto& active_indices = system.world.getAllActiveParticleIndices();
    const auto n_active = active_indices.size();

    for ( int i = 0; i < n_active; ++i) {
        for(  int j = i+1; j < n_active; ++j) {

            if(system.isInSameMolecule(i,j)){
                continue;
            }

            const auto atom_type_i = system.world.getParticle(active_indices.at(i)).atom_type;
            const auto& r_i = system.world.getParticle(active_indices.at(i)).r;
            const auto atom_type_j = system.world.getParticle(active_indices.at(j)).atom_type;
            const auto& r_j = system.world.getParticle(active_indices.at(j)).r;

            const auto A = system.lj_table.get_coeff_A(atom_type_i, atom_type_j);
            const auto B = system.lj_table.get_coeff_B(atom_type_i, atom_type_j);

            const auto dx = distancePBC(r_i[XX], r_j[XX],
                                        half_box_length[XX], box_length[XX]);
            const auto dy = distancePBC(r_i[YY], r_j[YY],
                                        half_box_length[YY], box_length[YY]);
            const auto dz = distancePBC(r_i[ZZ], r_j[ZZ],
                                        half_box_length[ZZ], box_length[ZZ]);
            const auto dr2 = dx*dx + dy*dy + dz*dz;
            const auto next_energy = LJ(dr2, A, B, cutoff_sq);
            energy += next_energy;
        }
    }
    return energy;
}

real Energy::calcVirial(const Systemx &system, const gmx::RVec box_length,
                        real cutoff_sq) const{
    real virial = 0.0;
    const auto half_box_length = box_length/2.;


    const auto& active_indices = system.world.getAllActiveParticleIndices();

    for (const auto  particle_index : active_indices) {
        const auto atom_type_i = system.world.getParticle(particle_index).atom_type;
        const auto& r_particle = system.world.getParticle(particle_index).r;

        const auto cell_index_n = system.world_to_cell[particle_index].first;
        const auto& neighbouring_cells = system.grid->cell2nearest_neighbours[cell_index_n];
        const auto& cell_list = system.cell_vectors2_;

        for( auto cell_index_m : neighbouring_cells){
            const auto n_particles_in_cell_m =cell_list[cell_index_m].atoms_in_cell;
            const auto& particles_cell_m = cell_list[cell_index_m].r;

            for( int i = 0; i<n_particles_in_cell_m; ++i ){

                if(system.isInSameMolecule(system.cell_to_world[cell_index_m][i], particle_index)){
                    continue;
                }

                const auto atom_type_j = particles_cell_m[i].atom_type;
                const auto A = system.lj_table.get_coeff_A(atom_type_i, atom_type_j);
                const auto B = system.lj_table.get_coeff_B(atom_type_i, atom_type_j);

                const auto dx = distancePBC(particles_cell_m[i].r[XX], r_particle[XX],
                                            half_box_length[XX], box_length[XX]) ;
                const auto dy = distancePBC(particles_cell_m[i].r[YY], r_particle[YY],
                                            half_box_length[YY], box_length[YY]);
                const auto dz = distancePBC(particles_cell_m[i].r[ZZ], r_particle[ZZ],
                                            half_box_length[ZZ], box_length[ZZ]);
                const gmx::RVec delta_r = {dx, dy, dz};
                assert(delta_r.norm() != 0.);
                auto next_vir = gmx::dot(LJ_force(delta_r, A, B, cutoff_sq), delta_r );

                virial += next_vir;
            }
        }
    }

    return virial/6.;
}

real Energy::calcVirialBasic(const Systemx &system, const gmx::RVec box_length,
                                 real cutoff_sq) {
    real virial = 0.0;
    const auto half_box_length = box_length/2.;

    const auto& active_indices = system.world.getAllActiveParticleIndices();
    const auto n_active = active_indices.size();

    for ( int i = 0; i < n_active; ++i) {
        for(  int j = i+1; j < n_active; ++j) {

            if(system.isInSameMolecule(i,j)){
                continue;
            }

            const auto atom_type_i = system.world.getParticle(active_indices.at(i)).atom_type;
            const auto& r_i = system.world.getParticle(active_indices.at(i)).r;
            const auto atom_type_j = system.world.getParticle(active_indices.at(j)).atom_type;
            const auto& r_j = system.world.getParticle(active_indices.at(j)).r;

            const auto A = system.lj_table.get_coeff_A(atom_type_i, atom_type_j);
            const auto B = system.lj_table.get_coeff_B(atom_type_i, atom_type_j);



            const gmx::RVec delta_r = deltarPBC(r_i, r_j,
                                                half_box_length, box_length);
            assert(delta_r.norm() != 0.);
            auto next_vir = gmx::dot(LJ_force(delta_r, A, B, cutoff_sq), delta_r );

            virial += next_vir;
        }
    }


    return virial/3.;
}


real Energy::calMoleculeMoveEnergyChange(  const Systemx& system,
                                             Molecule molecule,
                                             const std::vector<Particle>& particles,
                                             const gmx::RVec box_length,
                                             real cutoff_sq) {

    const auto half_box_length = box_length/2.;
    real energy_change_L6 = 0.0;
    real energy_change_L12 = 0.0;

    const auto& cell_list = system.cell_vectors2_;

    for(auto world_index = molecule.first_; world_index<=molecule.last_; world_index++) {

        const auto r_coords_new = particles.at(world_index - molecule.first_).r;
        const auto &cell_index_inserted = system.grid->coordToCell(r_coords_new);
        const auto &new_neighbouring_cells = system.grid->cell2nearest_neighbours[cell_index_inserted];

        const auto atom_type = particles[world_index - molecule.first_].atom_type;

        for (const auto cell_index_n: new_neighbouring_cells) {

            const auto &particles_cell_n = cell_list[cell_index_n].r;
            const auto n_particles_in_cell_n = cell_list[cell_index_n].atoms_in_cell;

#pragma clang loop vectorize(enable)
            for (int jn = 0; jn < n_particles_in_cell_n; ++jn) {
                const auto atom_type_neighbour = particles_cell_n[jn].atom_type;
                const auto A = system.lj_table.get_coeff_A(atom_type, atom_type_neighbour);
                const auto B = system.lj_table.get_coeff_B(atom_type, atom_type_neighbour);
                const auto LJ12_0 = system.lj_table.get_LJ12_0(atom_type, atom_type_neighbour);
                const auto LJ6_0 = system.lj_table.get_LJ6_0(atom_type, atom_type_neighbour);

                const auto neighbour_world_index = system.cell_to_world[cell_index_n][jn];
                bool is_in_same_molecule = neighbour_world_index >= molecule.first_ and  neighbour_world_index <= molecule.last_;

                const auto dx = distancePBC(particles_cell_n[jn].r[XX], r_coords_new[XX],
                                half_box_length[XX], box_length[XX]) + box_length[XX]*(is_in_same_molecule);
                const auto dy = distancePBC(particles_cell_n[jn].r[YY], r_coords_new[YY],
                                            half_box_length[YY], box_length[YY]);
                const auto dz = distancePBC(particles_cell_n[jn].r[ZZ], r_coords_new[ZZ],
                                            half_box_length[ZZ], box_length[ZZ]);
                const auto dr2 = dx * dx + dy * dy + dz * dz;
                const auto dr6 = dr2*dr2*dr2;
                const auto dr12 = dr6*dr6;

                energy_change_L6 += -B/dr6  * (dr2<=cutoff_sq) - LJ6_0;
                energy_change_L12 += A/dr12 * (dr2<=cutoff_sq) - LJ12_0; //LJ(dr2, A, B, cutoff_sq);
            }
        }

//        const auto r_coords_old = system.world.getParticle(world_index).r;
        const auto &cell_index_old = system.world_to_cell[world_index].first;
        const auto r_coords_old = system.cell_vectors2_[cell_index_old].r[system.world_to_cell[world_index].second].r;
        if(cell_index_old != system.grid->coordToCell(r_coords_old)){
            std::cout << cell_index_old << " " << system.grid->cellCoord(cell_index_old)[0] << " " << system.grid->cellCoord(cell_index_old)[1] << " " << system.grid->cellCoord(cell_index_old)[2]<<  '\n';
        }
        assert(cell_index_old == system.grid->coordToCell(r_coords_old));
        const auto &old_neighbouring_cells = system.grid->cell2nearest_neighbours[cell_index_old];
        for (const auto cell_index_n: old_neighbouring_cells) {

            const auto &particles_cell_n = cell_list[cell_index_n].r;
            const auto n_particles_in_cell_n = cell_list[cell_index_n].atoms_in_cell;

#pragma clang loop vectorize(enable)
            for (int jn = 0; jn < n_particles_in_cell_n; ++jn) {
                const auto atom_type_neighbour = particles_cell_n[jn].atom_type;
                const auto A = system.lj_table.get_coeff_A(atom_type, atom_type_neighbour);
                const auto B = system.lj_table.get_coeff_B(atom_type, atom_type_neighbour);

                const auto LJ12_0 = system.lj_table.get_LJ12_0(atom_type, atom_type_neighbour);
                const auto LJ6_0 = system.lj_table.get_LJ6_0(atom_type, atom_type_neighbour);

                const auto neighbour_world_index = system.cell_to_world[cell_index_n][jn];
                bool is_in_same_molecule = neighbour_world_index >= molecule.first_ and  neighbour_world_index <= molecule.last_;

                const auto dx = distancePBC(particles_cell_n[jn].r[XX], r_coords_old[XX],
                                            half_box_length[XX], box_length[XX]) + box_length[XX]*(is_in_same_molecule);
                const auto dy = distancePBC(particles_cell_n[jn].r[YY], r_coords_old[YY],
                                            half_box_length[YY], box_length[YY]);
                const auto dz = distancePBC(particles_cell_n[jn].r[ZZ], r_coords_old[ZZ],
                                            half_box_length[ZZ], box_length[ZZ]);
                const auto dr2 = dx * dx + dy * dy + dz * dz;
                const auto dr6 = dr2*dr2*dr2;
                const auto dr12 = dr6*dr6;

                energy_change_L6 -= -B/dr6  * (dr2<=cutoff_sq) - LJ6_0;
                energy_change_L12 -= A/dr12 * (dr2<=cutoff_sq) - LJ12_0; //LJ(dr2, A, B, cutoff_sq);
            }
        }
    }
    total_lj_energy.L6 += energy_change_L6;
    total_lj_energy.L12 += energy_change_L12;
    return energy_change_L12 + energy_change_L6 ;
}

real Energy::calcTotalEnergy(  const Systemx& system,
                                           const gmx::RVec box_length,
                                           real cutoff_sq) {

    const auto half_box_length = box_length/2.;
    real new_energy = 0.0;
    total_lj_energy.L6  = 0;
    total_lj_energy.L12 = 0;

    const auto& cell_list = system.cell_vectors2_;
    const auto& active_molecule_indices = system.world.active_molecules;

    for(const auto mol_index : active_molecule_indices) {
        auto molecule = system.world.all_molecules.at(mol_index);
        for(auto world_index = molecule.first_; world_index<=molecule.last_; ++world_index) {

            const auto particle = system.world.getParticle(world_index);
            const auto r_coords = particle.r;
            const auto atom_type = particle.atom_type;

            const auto &cell_index_inserted = system.grid->coordToCell(r_coords);
            const auto &new_neighbouring_cells = system.grid->cell2nearest_neighbours[cell_index_inserted];

            for (const auto cell_index_n: new_neighbouring_cells) {

                const auto &particles_cell_n = cell_list[cell_index_n].r;
                const auto n_particles_in_cell_n = cell_list[cell_index_n].atoms_in_cell;

#pragma clang loop vectorize(enable)
                for (int jn = 0; jn < n_particles_in_cell_n; ++jn) {

                    const auto neighbour_world_index = system.cell_to_world[cell_index_n][jn];
                    bool is_in_same_molecule =
                            neighbour_world_index >= molecule.first_ and neighbour_world_index <= molecule.last_;

                    const auto atom_type_neighbour = particles_cell_n[jn].atom_type;
                    const auto A = system.lj_table.get_coeff_A(atom_type, atom_type_neighbour);
                    const auto B = system.lj_table.get_coeff_B(atom_type, atom_type_neighbour);

                    const auto LJ12_0 = system.lj_table.get_LJ12_0(atom_type, atom_type_neighbour);
                    const auto LJ6_0 = system.lj_table.get_LJ6_0(atom_type, atom_type_neighbour);

                    const auto dx = distancePBC(particles_cell_n[jn].r[XX], r_coords[XX],
                                                half_box_length[XX], box_length[XX]) + box_length[XX] * (is_in_same_molecule);
                    const auto dy = distancePBC(particles_cell_n[jn].r[YY], r_coords[YY],
                                                half_box_length[YY], box_length[YY]);
                    const auto dz = distancePBC(particles_cell_n[jn].r[ZZ], r_coords[ZZ],
                                                half_box_length[ZZ], box_length[ZZ]);
                    const auto dr2 = dx * dx + dy * dy + dz * dz;

                    const auto idr2 = 1./dr2;
                    const auto idr6 = idr2*idr2*idr2;
//                    new_energy += (dr2<=cutoff_sq) * LJ2(idr2);
                    total_lj_energy.L6 += -B*idr6/2.0 * (dr2<=cutoff_sq) - LJ6_0;
                    total_lj_energy.L12 +=  A*idr6*idr6/2.0 * (dr2<=cutoff_sq) - LJ12_0;
                    new_energy += LJ(dr2, A, B, cutoff_sq);
                }
            }
        }
    }
    return new_energy/2.0;
}
/*

real Energy::calcScalingEnergyChange(  const Systemx& system,
                               const real box_length, const real scale,
                               real cutoff_sq) const{

    real half_box_length = box_length/2.;
    real new_energy = 0.0;

    const auto& cell_list = system.cell_vectors2_;
    const auto& active_world_indices = system.world.getAllActiveParticleIndices();

    for(const auto world_index : active_world_indices) {

        const auto particle = system.world.getParticle(world_index);
        const auto r_coords = particle.r;
        const auto &cell_index_inserted = system.grid->coordToCell(r_coords);
        const auto &new_neighbouring_cells = system.grid->cell2nearest_neighbours[cell_index_inserted];

        const auto atom_type = particle.atom_type;

        for (const auto cell_index_n: new_neighbouring_cells) {

            const auto &particles_cell_n = cell_list[cell_index_n].r;
            const auto n_particles_in_cell_n = cell_list[cell_index_n].atoms_in_cell;

#pragma clang loop vectorize(enable)
            for (int jn = 0; jn < n_particles_in_cell_n; ++jn) {
                if(system.isInSameMolecule(system.cell_to_world[cell_index_n][jn], world_index)){
                    continue;
                }
                const auto atom_type_neighbour = particles_cell_n[jn].atom_type;
                const auto A = system.lj_table.get_coeff_A(atom_type, atom_type_neighbour);
                const auto B = system.lj_table.get_coeff_B(atom_type, atom_type_neighbour);

                const auto dx = distancePBC(particles_cell_n[jn].r[XX], r_coords[XX],
                                            half_box_length, box_length) + box_length*(system.cell_to_world[cell_index_n][jn] == world_index);
                const auto dy = distancePBC(particles_cell_n[jn].r[YY], r_coords[YY],
                                            half_box_length, box_length);
                const auto dz = distancePBC(particles_cell_n[jn].r[ZZ], r_coords[ZZ],
                                            half_box_length, box_length);
                const auto dr2 = dx * dx + dy * dy + dz * dz;

                new_energy += LJ(dr2, A, B, cutoff_sq);
            }
        }
    }
    return new_energy/2.0;
}
*/



LJComponents Energy::calcMoleculeEnergy(  const Systemx& system,
                                     Molecule molecule,
                                     const gmx::RVec box_length,
                                     real cutoff_sq) const{

    const auto half_box_length = box_length/2.;
    LJComponents lj_energy_components = {0.0, 0.0};

    const auto& cell_list = system.cell_vectors2_;

    assert(system.world.isActiveMolecule(molecule.world_index_));

    for(auto world_index = molecule.first_; world_index<=molecule.last_; world_index++) {

        const auto r_coords_new = system.world.getParticle(world_index).r;
        const auto inserted_atom_type = system.world.getParticle(world_index).atom_type;

        const auto &cell_index_inserted = system.grid->coordToCell(r_coords_new);
        const auto &neighbouring_cells = system.grid->cell2nearest_neighbours[cell_index_inserted];

        for (const auto cell_index_n: neighbouring_cells) {

            const auto &r_cell_n = cell_list[cell_index_n].r;
            const auto n_particles_in_cell_n = cell_list[cell_index_n].atoms_in_cell;

#pragma clang loop vectorize(enable)
            for (int jn = 0; jn < n_particles_in_cell_n; ++jn) {
                const auto atom_type_neighbour = r_cell_n[jn].atom_type;
                const auto A = system.lj_table.get_coeff_A(inserted_atom_type, atom_type_neighbour);
                const auto B = system.lj_table.get_coeff_B(inserted_atom_type, atom_type_neighbour);

                const auto world_index_neighbour = system.cell_to_world[cell_index_n][jn];
                bool is_in_molecule = molecule.first_ <= world_index_neighbour and molecule.last_ >= world_index_neighbour ;
                if(!is_in_molecule) {
                    const auto dx = distancePBC(r_cell_n[jn].r[XX], r_coords_new[XX],
                                                half_box_length[XX], box_length[XX]);
                    const auto dy = distancePBC(r_cell_n[jn].r[YY], r_coords_new[YY],
                                                half_box_length[YY], box_length[YY]);
                    const auto dz = distancePBC(r_cell_n[jn].r[ZZ], r_coords_new[ZZ],
                                                half_box_length[ZZ], box_length[ZZ]);
                    const auto dr2 = dx * dx + dy * dy + dz * dz;
                    const auto dr6 = dr2*dr2*dr2;
                    const auto dr12 = dr6*dr6;

                    lj_energy_components.L6 += -B/dr6 * (dr2<=cutoff_sq);
                    lj_energy_components.L12 += A/dr12 * (dr2<=cutoff_sq);
                    if(!std::isfinite(lj_energy_components.L6)){
                        std::cout << "cecky\n";
                        auto p = system.world.getParticle(world_index_neighbour);
                    }
                }
            }
        }
    }
    return lj_energy_components;
}

real Energy::calcParticlesEnergy(  const Systemx& system,
                                  const std::vector<Particle>& particles,
                                  const gmx::RVec box_length,
                                  real cutoff_sq) const{

    const auto half_box_length = box_length/2.;
    real new_energy = 0.0;

    const auto& cell_list = system.cell_vectors2_;

    for(const auto& particle : particles ) {

        const auto r_coords_new =particle.r;
        const auto &cell_index_inserted = system.grid->coordToCell(r_coords_new);
        const auto &neighbouring_cells = system.grid->cell2nearest_neighbours[cell_index_inserted];

        const auto inserted_atom_type = particle.atom_type;

        for (const auto cell_index_n: neighbouring_cells) {

            const auto &r_cell_n = cell_list[cell_index_n].r;
            const auto n_particles_in_cell_n = cell_list[cell_index_n].atoms_in_cell;

#pragma clang loop vectorize(enable)
            for (int jn = 0; jn < n_particles_in_cell_n; ++jn) {
                const auto atom_type_neighbour = r_cell_n[jn].atom_type;
                const auto A = system.lj_table.get_coeff_A(inserted_atom_type, atom_type_neighbour);
                const auto B = system.lj_table.get_coeff_B(inserted_atom_type, atom_type_neighbour);

                const auto world_index_neighbour = system.cell_to_world[cell_index_n][jn];

                const auto dx = distancePBC(r_cell_n[jn].r[XX], r_coords_new[XX],
                                            half_box_length[XX], box_length[XX]);
                const auto dy = distancePBC(r_cell_n[jn].r[YY], r_coords_new[YY],
                                            half_box_length[YY], box_length[YY]);
                const auto dz = distancePBC(r_cell_n[jn].r[ZZ], r_coords_new[ZZ],
                                            half_box_length[ZZ], box_length[ZZ]);
                const auto dr2 = dx * dx + dy * dy + dz * dz;

                new_energy += LJ(dr2, A, B, cutoff_sq);
            }
        }
    }
    return new_energy;
}
real Energy::calcNewMoleculeEnergy(  const Systemx& system,
                                           Molecule new_molecule,
                                           const std::vector<Particle>& new_particles,
                                           const gmx::RVec box_length,
                                           real cutoff_sq){

    const auto half_box_length = box_length/2.;
    real new_energy_L6 = 0.0;
    real new_energy_L12 = 0.0;

    const auto& cell_list = system.cell_vectors2_;

    for(auto world_index = new_molecule.first_; world_index<=new_molecule.last_; world_index++) {

        const auto r_coords_new = new_particles[world_index-new_molecule.first_].r;
        const auto &cell_index_inserted = system.grid->coordToCell(r_coords_new);
        const auto &neighbouring_cells = system.grid->cell2nearest_neighbours[cell_index_inserted];

        const auto inserted_atom_type = new_particles[world_index-new_molecule.first_].atom_type;

        for (const auto cell_index_n: neighbouring_cells) {

            const auto &r_cell_n = cell_list[cell_index_n].r;
            const auto n_particles_in_cell_n = cell_list[cell_index_n].atoms_in_cell;

#pragma clang loop vectorize(enable)
            for (int jn = 0; jn < n_particles_in_cell_n; ++jn) {
                const auto atom_type_neighbour = r_cell_n[jn].atom_type;
                const auto A = system.lj_table.get_coeff_A(inserted_atom_type, atom_type_neighbour);
                const auto B = system.lj_table.get_coeff_B(inserted_atom_type, atom_type_neighbour);

                const auto LJ12_0 = system.lj_table.get_LJ12_0(inserted_atom_type, atom_type_neighbour);
                const auto LJ6_0 = system.lj_table.get_LJ6_0(inserted_atom_type, atom_type_neighbour);

                const auto dx = distancePBC(r_cell_n[jn].r[XX], r_coords_new[XX],
                                            half_box_length[XX], box_length[XX]);
                const auto dy = distancePBC(r_cell_n[jn].r[YY], r_coords_new[YY],
                                            half_box_length[YY], box_length[YY]);
                const auto dz = distancePBC(r_cell_n[jn].r[ZZ], r_coords_new[ZZ],
                                            half_box_length[ZZ], box_length[ZZ]);
                const auto dr2 = dx * dx + dy * dy + dz * dz;
                const auto dr6 = dr2*dr2*dr2;
                const auto dr12 = dr6*dr6;

                new_energy_L12 += A/dr12 * (dr2<=cutoff_sq) - LJ12_0;
                new_energy_L6  +=-B/dr6  * (dr2<=cutoff_sq) - LJ6_0;
//                new_energy += LJ(dr2, A, B, cutoff_sq);
            }
        }
    }
    total_lj_energy.L12 += new_energy_L12;
    total_lj_energy.L6 += new_energy_L6;

    return new_energy_L12 + new_energy_L6;
}


real evaluatePotentials(const std::vector<PotentialData>& potential_data, const real cutoff_sq ){
    real sum = 0.0;

    for(const auto& datum : potential_data){
        sum += LJ(datum.dr2, datum.A, datum.B, cutoff_sq);
    }
    return sum;
}


real energyCorrection(const Systemx& system, const Energy& energizer,
                      const real volume,
                      const real cutoff){
    real energy_correction = 0.;
    size_t n_atom_types = system.atom_type_to_lj_params.size();
    for(int atom_type_i = 0; atom_type_i < n_atom_types; ++atom_type_i){
        for(int atom_type_j = 0; atom_type_j < n_atom_types; ++atom_type_j) {

            const auto sig = system.lj_table.get_sig(atom_type_i, atom_type_j);
            const auto eps = system.lj_table.get_eps(atom_type_i, atom_type_j);

            const real n_atoms_i = system.getActiveAtomsOf(atom_type_i).size();
            const real n_atoms_j = system.getActiveAtomsOf(atom_type_j).size();

            const real rhoi_red = (n_atoms_i/volume) * std::pow(sig, 3);
            const real rhoj_red = (n_atoms_j/volume) * std::pow(sig, 3);

            energy_correction += eps * 8. / 3. * M_PI * n_atoms_i * rhoj_red *
                                 (1. / 3. * std::pow(sig / cutoff, 9) -
                                  std::pow(sig / cutoff, 3));
        }
    }
    return energy_correction ;
}

real pressureCorrection2(const Systemx& system, const Energy& energizer,
                        const gmx::RVec box_size, const real cutoff) {

    const auto volume = box_size[XX] * box_size[YY] * box_size[ZZ];

    const size_t n_atom_types = system.atom_type_to_lj_params.size();
    real pressure_correction = 0.0;
    auto total_lj_energy_comp_new = energizer.total_lj_energy_comp;
    for (int atom_type_i = 0; atom_type_i < n_atom_types; ++atom_type_i) {
        for (int atom_type_j = 0; atom_type_j < n_atom_types; ++atom_type_j) {
            const auto sig = system.lj_table.get_sig(atom_type_i, atom_type_j);
            const auto eps = system.lj_table.get_eps(atom_type_i, atom_type_j);

            const real n_atoms_i = system.getActiveAtomsOf(atom_type_i).size();
            const real n_atoms_j = system.getActiveAtomsOf(atom_type_j).size();

            const real rho_red_i = ((real) n_atoms_i / volume) * std::pow(sig, 3);
            const real rho_red_j = ((real) n_atoms_j / volume) * std::pow(sig, 3);
            pressure_correction +=
                    16. / 3. * M_PI * rho_red_i * rho_red_j * std::pow(sig, -3) * eps *
                    (2. / 3. * std::pow(cutoff / sig, -9) - std::pow(cutoff / sig, -3)) * 1.0e30 / 6.022e23 *
                    1e3 * 1e-5;
        }
    }
    return pressure_correction;
}

