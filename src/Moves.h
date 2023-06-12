#include "Grid.h"
#include "core.h"
#include <random>

#ifndef GR_GCMC_CLION_MOVES_H
#define GR_GCMC_CLION_MOVES_H

static int n_total = 0;
static int n_accepted = 0;


struct MoveData{
    size_t n_accepted;
    size_t n_total;
    const size_t n_attempts;
};

struct MovesSamplingData{
    MoveData translation;
    MoveData exchange;
    MoveData volume;
};


real energyCorrectionChange(const Systemx& system, const real V_old, const real V_new, const real cutoff){
    real energy_correction_old = 0.;
    real energy_correction_new = 0.;
    size_t n_atom_types = system.atom_type_to_lj_params.size();
    for(int i = 0; i < n_atom_types; ++i){
        for(int j = 0; j < n_atom_types; ++j){
            const auto n_atoms_i = system.getActiveAtomsOf(i).size();
            const auto n_atoms_j = system.getActiveAtomsOf(j).size();

            const auto sig = system.lj_table.get_sig(i, j);
            const auto eps = system.lj_table.get_eps(i, j);

            const real rho_red_old = ((real) n_atoms_j / V_old) * std::pow(sig, 3);
            const real rho_red_new = ((real) n_atoms_j / V_new) * std::pow(sig, 3);

            energy_correction_old += eps * 8. / 3. * M_PI * (real) n_atoms_i * rho_red_old *
                                     (1. / 3. * std::pow(sig / cutoff, 9) - std::pow(sig / cutoff, 3));
            energy_correction_new += eps * 8. / 3. * M_PI * (real) n_atoms_i * rho_red_new *
                                     (1. / 3. * std::pow(sig / cutoff, 9) - std::pow(sig / cutoff, 3));

        }
    }
    return energy_correction_new - energy_correction_old;
}


real energyCorrectionChangeScaledCutoff(const Systemx& system, const real V_old, const real V_new, const real cutoff){
    real energy_correction_old = 0.;
    real energy_correction_new = 0.;
    real scale = std::pow(V_new/V_old, 1./3.);
    size_t n_atom_types = system.atom_type_to_lj_params.size();
    for(int i = 0; i < n_atom_types; ++i){
        for(int j = 0; j < n_atom_types; ++j){
            const auto n_atoms_i = system.getActiveAtomsOf(i).size();
            const auto n_atoms_j = system.getActiveAtomsOf(j).size();

            const auto sig = system.lj_table.get_sig(i, j);
            const auto eps = system.lj_table.get_eps(i, j);

            const real rho_red_old = ((real) n_atoms_j / V_old) * std::pow(sig, 3);
            const real rho_red_new = ((real) n_atoms_j / V_new) * std::pow(sig, 3);

            energy_correction_old += eps * 8. / 3. * M_PI * (real) n_atoms_i * rho_red_old *
                                     (1. / 3. * std::pow(sig / cutoff, 9) - std::pow(sig / cutoff, 3));
            energy_correction_new += eps * 8. / 3. * M_PI * (real) n_atoms_i * rho_red_new *
                                     (1. / 3. * std::pow(sig / (cutoff * scale), 9) -
                                      std::pow(sig / (cutoff * scale), 3));

        }
    }
    return energy_correction_new - energy_correction_old;
}



void deletionMove(int& step,
                  VerletList* verlet_list,
                  Systemx& system, Energy& energizer,
                  size_t& n_atoms,
                  gmx::RVec box_size,
                  const real cutoff, const real beta, MoveData& move_statistics,
                  std::mt19937 & e1, std::uniform_real_distribution<real>& uniform_dist) {


    const size_t n_atom_types = system.atom_type_to_lj_params.size();
    const auto random_mol_type = system.gc_molecule_types.at(rand()%(system.gc_molecule_types.size()));
    const auto& molecules_of_removed_type = system.getActiveMolecules(random_mol_type);
    const size_t n_mols_of_random_type = molecules_of_removed_type.size();
    assert(n_mols_of_random_type > 0);
    Molecule random_molecule = molecules_of_removed_type.at(rand()%(n_mols_of_random_type));



    real energy_correction_change = 0.0;
    const auto volume =  box_size[XX]*box_size[YY]*box_size[ZZ];
/*
    for(int world_index = random_molecule.first_; world_index<=random_molecule.last_; world_index++) {
        const auto random_type = system.world.all_particles.at(world_index).atom_type;
        for (int i = 0; i < n_atom_types; ++i) {

            const auto sig = system.lj_table.get_sig(i, random_type);
            const auto eps = system.lj_table.get_eps(i, random_type);
            const auto factor =
                    eps * 8. / 3. * M_PI * std::pow(sig, 3) / volume * (1. / 3. * std::pow(sig / cutoff, 9) -
                                                                        std::pow(sig / cutoff, 3));

            const real n_atoms_i = system.getActiveAtomsOf(i).size();// - 1*(i==random_type);

            energy_correction_change -= (2 * n_atoms_i - 1 * (i == random_type)) * factor;
        }
    }
*/

    auto new_particles = system.world.getCopyOfMoleculeCoords(random_molecule.world_index_);

    const auto old_total_lj_energy = energizer.total_lj_energy;
    const auto removed_particle_energy_lj_components = energizer.calcMoleculeEnergy(system, random_molecule,
                                                                                     box_size, cutoff*cutoff);
    const auto removed_particle_energy = removed_particle_energy_lj_components.L12 +  removed_particle_energy_lj_components.L6;

    const real bias_factor = (real)n_mols_of_random_type /volume;
    const real chem_pot = system.chemical_potentials.at(random_mol_type);
    auto mc_factor = std::exp(-(chem_pot-removed_particle_energy + energy_correction_change) * beta) * bias_factor;
    const auto roll = uniform_dist(e1);



    if (roll < min(mc_factor, 1.0)) {

        system.removeMolecule(random_molecule, random_mol_type);

        n_atoms--;
        n_accepted++;
        move_statistics.n_accepted++;

/*        verlet_list->particle_to_list[random_particle] = verlet_list->particle_to_list[n_atoms];
        verlet_list->particle_to_list.pop_back();*/
        energizer.sum_of_energy_changed += -removed_particle_energy + energy_correction_change;
        energizer.total_lj_energy.L6 = old_total_lj_energy.L6 - removed_particle_energy_lj_components.L6;
        energizer.total_lj_energy.L12 = old_total_lj_energy.L12 - removed_particle_energy_lj_components.L12;
    } else {
        energizer.total_lj_energy = old_total_lj_energy;
    }
    n_total++;
    move_statistics.n_total++;
    step++;
}

void insertionMove(int& step,  VerletList* verlet_list,
                   Systemx& system, Energy& energizer,
                   size_t& n_atoms,
                   gmx::RVec box_size,
                   const real cutoff, const real beta, MoveData& move_statistics,
                   std::mt19937 & e1, std::uniform_real_distribution<real>& uniform_dist) {

    const auto random_mol_type = system.gc_molecule_types.at(rand()%(system.gc_molecule_types.size()));
    assert(system.world.inactive_molecules.size() >0);
    const auto random_molecule = system.world.randomInactiveMoleculeOf(random_mol_type);
    assert(system.world.isInactiveMolecule(random_molecule.world_index_));

    auto new_particles = system.molecule_type_to_structure.at(random_mol_type);
    generateRandomPosition(new_particles, box_size, uniform_dist, e1);

    const auto old_lj_component_energy = energizer.total_lj_energy;
    const auto new_molecule_energy = energizer.calcNewMoleculeEnergy(system, random_molecule, new_particles,
                                                                              box_size, cutoff*cutoff);

    real energy_correction = 0.0;
    const real volume = box_size[XX]*box_size[YY]*box_size[ZZ];
    const size_t n_atom_types = Systemx::n_atom_types;
/*    for(int world_index = random_molecule.first_; world_index<=random_molecule.last_; world_index++){
        const auto random_atom_type = system.world.all_particles.at(world_index).atom_type;
        for(int i = 0; i < n_atom_types; ++i){

            const auto sig = system.lj_table.get_sig(i, random_atom_type);
            const auto eps = system.lj_table.get_eps(i, random_atom_type);

            const auto factor = eps * 8. / 3. * M_PI * std::pow(sig, 3) / volume * (1. / 3. * std::pow(sig / cutoff, 9) -
                                                                                    std::pow(sig / cutoff, 3));

            const real n_atoms_i = system.getActiveAtomsOf(i).size();// - 1*(i==random_type);

            energy_correction += (2* n_atoms_i  - 1*(i==random_atom_type) )* factor  ;
        }
    }*/



    const real n_molecules_random_type = static_cast<real>(system.getActiveMolecules(random_mol_type).size());
    const real bias_factor = volume/(real)(n_molecules_random_type+1);
    const auto chem_pot = system.chemical_potentials.at(random_mol_type);
    const auto mc_factor = std::exp((chem_pot - (new_molecule_energy + energy_correction)) * beta) * bias_factor;
    const auto roll = uniform_dist(e1);

    if (roll < min(mc_factor, 1.0)) {

        energizer.sum_of_energy_changed += new_molecule_energy + energy_correction;
        n_atoms++;
        n_accepted++;
        move_statistics.n_accepted++;
        system.insertMolecule(new_particles, random_molecule, random_mol_type);

    } else {
        energizer.total_lj_energy = old_lj_component_energy;
    }
    n_total++;
    move_statistics.n_total++;
    step++;

}


void volumeMove(int& step,
                Systemx& system,
                Energy& energizer,
                GridSearch* grid_searcher,
                gmx::RVec& box_size,
                const real& cutoff, const real beta, const real P, MoveData& move_statistics,
                std::mt19937 & e1, std::uniform_real_distribution<real>& uniform_dist) {

    const real delta_V_max = 500;

    const auto delta_V = (uniform_dist(e1) - 0.5) * delta_V_max;
    const auto V_old =  box_size[XX]*box_size[YY]*box_size[ZZ];
    const auto V_new = V_old + delta_V;
    const auto scale = std::pow(V_new/ V_old, 1./3.);

    const auto en_before =  energizer.total_energy;
    system.scaleAllCoordinates(scale);
    const auto new_cutoff = cutoff;// * scale;
    const auto new_box_size = box_size*scale;
    const auto cell_size = system.grid->cell_size_;
    system.grid->cell_size_ *= scale;
    system.grid->box_size_ = new_box_size;
    const auto en_after =  energizer.calcTotalEnergy(system, new_box_size, new_cutoff*new_cutoff);


    real energy_correction_before  = energyCorrection(system, energizer, V_old, cutoff);
    real energy_change = en_after - (en_before - energy_correction_before);
    real energy_correction_after  = energyCorrection(system, energizer, V_new, cutoff);
    real energy_correction_change  = energy_correction_after - energy_correction_before;

    auto n_particles = system.world.active_molecules.size();
    const auto bias_factor = std::exp((real) n_particles * std::log(V_new / V_old));
    const auto mc_factor = std::exp(-( energy_change + energy_correction_change + P * (V_new - V_old) ) * beta) * bias_factor;

    const auto roll = uniform_dist(e1);
    if (roll < min(mc_factor, 1.0) and std::isfinite(mc_factor)) {

        energizer.sum_of_energy_changed += energy_change + energy_correction_change;
        energizer.total_energy += energy_change + energy_correction_change;

        box_size = new_box_size;
//        cutoff = new_cutoff;
        system.cutoff = new_cutoff;

        grid_searcher->updateNearestNeighbours(new_box_size, cutoff);

        n_accepted++;
        move_statistics.n_accepted++;
    } else {
        system.grid->cell_size_ = cell_size;
        system.grid->box_size_ = box_size;
//        system.cutoff = cutoff;
        system.scaleAllCoordinates2(scale);
    }
    n_total++;
    move_statistics.n_total++;
    step++;
}


void volumeMoveScaledCutoff(int& step,
                            Systemx& system,
                            Energy& energizer,
                            gmx::RVec& box_size, real& cutoff,
                            const real beta, const real P, MoveData& move_statistics,
                            std::mt19937 & e1, std::uniform_real_distribution<real>& uniform_dist) {

    const real delta_V_max = 500;

    const auto delta_V = (uniform_dist(e1) - 0.5) * delta_V_max;
    const auto V_old =  box_size[XX]*box_size[YY]*box_size[ZZ];
    const auto V_new = V_old + delta_V;
    const auto scale = std::pow(V_new/ V_old, 1./3.);

    real cutoff_old = system.cutoff;
    const auto en_before =  energizer.total_energy;
    system.scaleAllCoordinates(scale);
    const auto cutoff_new = cutoff_old* scale;
    const auto new_box_size = box_size*scale;
    const auto cell_size = system.grid->cell_size_;
    system.grid->cell_size_ *= scale;
    system.grid->box_size_ = new_box_size;

    const auto LJ_6_after =  energizer.total_lj_energy.L6*std::pow(scale, -6);
    const auto LJ_12_after =  energizer.total_lj_energy.L12*std::pow(scale, -12);
    const auto en_after  = LJ_12_after + LJ_6_after + energyCorrection(system, energizer, V_new, cutoff_new);

    real energy_change = en_after - en_before;

    auto n_particles = system.world.active_molecules.size();
    const auto bias_factor = std::exp((real) n_particles * std::log(V_new / V_old));
    const auto mc_factor = std::exp(-( energy_change + P * (V_new - V_old) ) * beta) * bias_factor;

    const auto roll = uniform_dist(e1);
    if (roll < min(mc_factor, 1.0) and std::isfinite(mc_factor)) {
        energizer.sum_of_energy_changed += energy_change;
        energizer.total_energy += energy_change;
        energizer.total_lj_energy.L12 = LJ_12_after;
        energizer.total_lj_energy.L6 = LJ_6_after;


        box_size = new_box_size;
        system.cutoff = cutoff_new;
        cutoff = system.cutoff;

        n_accepted++;
        move_statistics.n_accepted++;
    } else {
        system.grid->cell_size_ = cell_size;
        system.grid->box_size_ = box_size;
//        system.cutoff = cutoff;
        system.scaleAllCoordinates2(scale);

    }
    n_total++;
    move_statistics.n_total++;
    step++;
}


void translationMove2(int& step,
                       VerletList* verlet_list,
                       Systemx& system, Energy& energizer,
                       gmx::RVec box_size,
                       const real dx0, const real cutoff, MoveData& move_statistics,
                       std::mt19937 & e1, std::uniform_real_distribution<real>& uniform_dist) {
    real bias_factor = (1.0);

    step++;
    const real beta = system.beta;

    const gmx::RVec rand_displacement = {(uniform_dist(e1) - 0.5) * dx0,
                                         (uniform_dist(e1) - 0.5) * dx0,
                                         (uniform_dist(e1) - 0.5) * dx0};

    const auto random_molecule_type = rand()%system.n_molecule_types;
    const auto& active_molecules = system.getActiveMolecules(random_molecule_type);
    const auto n_active_molecules = active_molecules.size();
    assert(n_active_molecules>0);
    const auto random_index_of_active_molecule = rand()%n_active_molecules;
    const auto random_molecule = active_molecules.at(random_index_of_active_molecule);

    auto new_coordinates = system.world.getCopyOfMoleculeCoords(random_molecule.world_index_);

    std::vector<PotentialData> a(0);
    translateMolecule(new_coordinates, rand_displacement, box_size);

    const auto old_total_lj_energies = energizer.total_lj_energy_comp;
/*    real energy_change =
            energizer.calMoleculeMoveEnergyChange(system, random_molecule,
                                                  new_coordinates,
                                                  box_size, cutoff * cutoff);*/

    verlet_list->fillPotentialData(system, random_molecule, new_coordinates, a, box_size);
    verlet_list->fillPotentialData(system, random_molecule, a, box_size);

    auto energy_change =  evaluatePotentials(a, cutoff*cutoff);
//    std::cout << energy_change << " :aaa: " << energy_change2 << "\n";


    const auto mc_factor = std::exp(-(energy_change) * beta) * bias_factor;
    auto roll = uniform_dist(e1);

    if (roll < min(mc_factor, 1.0) and std::isfinite(mc_factor)) {

        system.updateMolecule(random_molecule, new_coordinates);
        energizer.sum_of_energy_changed += energy_change;
        energizer.total_energy += energy_change;
        n_accepted++;
        move_statistics.n_accepted++;

    } else {
    }
    move_statistics.n_total++;
    n_total++;
}

void translationMove(int& step,
                       Systemx& system, Energy& energizer,
                       gmx::RVec box_size,
                       const real dx0, const real cutoff, MoveData& move_statistics,
                       std::mt19937 & e1, std::uniform_real_distribution<real>& uniform_dist) {
    real bias_factor = (1.0);

    step++;
    const real beta = system.beta;

    const gmx::RVec rand_displacement = {(uniform_dist(e1) - 0.5) * dx0,
                                        (uniform_dist(e1) - 0.5) * dx0,
                                        (uniform_dist(e1) - 0.5) * dx0};

    const auto random_molecule_type = rand()%system.n_molecule_types;
    const auto& active_molecules = system.getActiveMolecules(random_molecule_type);
    const auto n_active_molecules = active_molecules.size();

    assert(n_active_molecules>0);
    const auto random_index_of_active_molecule = rand()%n_active_molecules;
    const auto random_molecule = active_molecules.at(random_index_of_active_molecule);

    auto new_coordinates = system.world.getCopyOfMoleculeCoords(random_molecule.world_index_);
    translateMolecule(new_coordinates, rand_displacement, box_size);

    const auto old_total_lj_energy = energizer.total_lj_energy;
    const real energy_change =
            energizer.calMoleculeMoveEnergyChange(system, random_molecule,
                                                new_coordinates,
                                                box_size, cutoff * cutoff);

    assert(std::isfinite(energy_change));
    if(!std::isfinite(energy_change)){
        std::cout << "your energy is infinite and you probably fucked up!" << "\n";
    }


    const auto mc_factor = std::exp(-(energy_change) * beta) * bias_factor;
    auto roll = uniform_dist(e1);

    if (roll < min(mc_factor, 1.0) and std::isfinite(mc_factor)) {

        system.updateMolecule(random_molecule, new_coordinates);
        energizer.sum_of_energy_changed += energy_change;
        energizer.total_energy += energy_change;
        n_accepted++;
        move_statistics.n_accepted++;

    } else {
        energizer.total_lj_energy = old_total_lj_energy;
    }
    move_statistics.n_total++;
    n_total++;
}

real specialMoveDisplacements(Grid& cell_list, const size_t cell_index,
                              const real&  x, const real& y, const real& z,
                              real& x_new, real& y_new, real& z_new,
                              const real dx0, const real dy0, const real dz0,
                              int particle_index, real box_size,
                              std::default_random_engine& e1, std::uniform_real_distribution<real>& rand){

    const auto cell_size = cell_list.cell_size;

    const auto i_cell_x = cell_list.cellCoordX(cell_index);
    const auto i_cell_y = cell_list.cellCoordY(cell_index);
    const auto i_cell_z = cell_list.cellCoordZ(cell_index);


    const real lx = (i_cell_x+1)*cell_size - x + (i_cell_x == cell_list.nx-1)*(box_size - (i_cell_x+1)*cell_size);
    real rx = x - i_cell_x*cell_size ;
    const real ly = (i_cell_y+1)*cell_size - y + (i_cell_y == cell_list.ny-1)*(box_size - (i_cell_y+1)*cell_size);
    real ry = y - i_cell_y*cell_size ;
    const real lz = (i_cell_z+1)*cell_size - z + (i_cell_z == cell_list.nz-1)*(box_size - (i_cell_z+1)*cell_size);
    real rz = z - i_cell_z*cell_size ;

    real dx_r = dx0 + (rx<=dx0)*(rx - dx0 )  + GMX_FLOAT_EPS;
    real dx_l = dx0 + (lx<=dx0)*(lx - dx0 );
    real dy_r = dy0 + (ry<=dy0)*(ry - dy0 ) +  GMX_FLOAT_EPS;
    real dy_l = dy0 + (ly<=dy0)*(ly - dy0 );
    real dz_r = dz0 + (rz<=dz0)*(rz - dz0 )  + GMX_FLOAT_EPS;
    real dz_l = dz0 + (lz<=dz0)*(lz - dz0 ) ;

    /// Calculate volume of possible jump destinations not leaving current cell
    const auto Vij = (dx_r + dx_l)*(dy_r + dy_l)*(dz_r + dz_l);

    /// Make the actual move.
    std::vector<real> dr(3);
    dr[0] = rand(e1)*(dx_r + dx_l) - dx_r;
    dr[1] = rand(e1)*(dy_r + dy_l) - dy_r;
    dr[2] = rand(e1)*(dz_r + dz_l) - dz_r;


    x_new = (x + dr[0]);
    x_new += - (box_size)*(x_new>=box_size) + (box_size)*(x_new<0);
    y_new = (y + dr[1]);
    y_new += - (box_size)*(y_new>=box_size) + (box_size)*(y_new<0);
    z_new = (z + dr[2]);
    z_new += - (box_size)*(z_new>=box_size) + (box_size)*(z_new<0);

/*    const auto new_cell_index = cell_list.coordToCell(x_new,y_new,z_new);
    if( cell_list.coordToCell(x_new,y_new,z_new) != cell_index ){

        const auto i_cell_new_x = cell_list.cellCoordX(new_cell_index);
        const auto i_cell_new_y = cell_list.cellCoordY(new_cell_index);
        const auto i_cell_new_z = cell_list.cellCoordZ(new_cell_index);
        std::cout<< x_new << " " << y_new << " " << z_new << "\n";
        std::cout<< i_cell_x << " " << i_cell_y<< " " << i_cell_z << "\n";
        std::cout<< i_cell_new_x << " " << i_cell_new_y<< " " << i_cell_new_z << "\n";

        std::cout << "index of problem is: " << particle_index << "\n";
    }*/
//    assert( cell_list.coordToCell(x_new,y_new,z_new) == cell_index );

    rx = x_new - i_cell_x*cell_size ;
    ry = y_new - i_cell_y*cell_size ;
    rz = z_new - i_cell_z*cell_size ;

    dx_r = dx0 + (rx<dx0)*(rx - dx0);
    dx_l = dx0 + (lx<dx0)*(lx - dx0);
    dy_r = dy0 + (ry<dy0)*(ry - dy0);
    dy_l = dy0 + (ly<dy0)*(ly - dy0);
    dz_r = dz0 + (rz<dz0)*(rz - dz0);
    dz_l = dz0 + (lz<dz0)*(lz - dz0);
    /// Calculate volume of possible jump destinations from move endpoint not leaving the current cell
    const auto Vji = (dx_r + dx_l)*(dy_r + dy_l)*(dz_r + dz_l);

    return  Vji/Vij;
}






#endif //GR_GCMC_CLION_MOVES_H
