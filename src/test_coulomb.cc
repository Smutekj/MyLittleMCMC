#include <gtest/gtest.h>
#include "Coulomb.h"
#include "vectypes.h"
#include <random>
#include "core.h"
#include "Grid.h"


// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
// Expect equality.

    real box_size = 50;
    real density = 0.055;
    size_t n_atoms = std::floor(density * box_size*box_size*box_size);
    n_atoms -= n_atoms%2;
    real r_cutoff = 10.;
    real k_cutoff = 1.0;
    real beta = 5./box_size;
    const real epsilon_rel = 10.;

    std::vector<gmx::RVec> r_coords(n_atoms);
    std::vector<real> charges(n_atoms);
    std::vector<real> As(n_atoms);
    std::vector<real> Bs(n_atoms);

    real cell_size = r_cutoff;

    std::unique_ptr<Grid> cell_list = std::make_unique<Grid>(cell_size, box_size);

    std::random_device r;
    std::default_random_engine e1(10);
    std::uniform_real_distribution<real> uniform_dist(0., 1.0);

    std::vector<unsigned long> all_indices(n_atoms);

    real total_charge = 0;

    int n_atoms_tmp = 1;
    As[0] = 4.0 * 1.0 * (std::pow(3.3, 12)) ;
    r_coords[0] = {0, 0 , 0};
    charges[0] = -1.0;
    while(  n_atoms_tmp < n_atoms) {
        auto x = box_size * uniform_dist(e1);
        auto y = box_size * uniform_dist(e1);
        auto z = box_size * uniform_dist(e1);
        auto r_coords_new = gmx::RVec(x, y, z);

        bool far_enough = true;
        for (int i = 0; i < n_atoms_tmp; ++i) {
            far_enough = far_enough and dist_sq_pbc(r_coords[i], r_coords_new, box_size/2., box_size) > 3.*3.;
        }
        if(far_enough){
            r_coords[n_atoms_tmp] = r_coords_new;
            charges[n_atoms_tmp] = (2 * (n_atoms_tmp % 2) - 1);
            As[n_atoms_tmp] = 4.0 * 1.0 * (std::pow(3.3, 12)) ;
            Bs[n_atoms_tmp] = 4.0 * 1.0 * (std::pow(3.3, 6)) ;
            total_charge += charges[n_atoms_tmp];
            n_atoms_tmp++;

        }
    }

    std::cout << "\n total charge is: " << total_charge << "\n";

    std::vector<CellVector<gmx::RVec>> cell_vectors(cell_list->nx*cell_list->ny*cell_list->nz);
    std::vector<CellVector<gmx::RVec>> cell_vectors2(cell_list->nx*cell_list->ny*cell_list->nz);
    Systemx system(r_coords.data(), n_atoms, cell_list.get(), cell_vectors.data());
    Systemx system2(r_coords.data(), n_atoms, cell_list.get(), cell_vectors2.data());

    auto lj_old = Electrostatics::BasicEnergyLj(&system, As.data(), Bs.data(), box_size);

    /// Direct
    auto total_energy_basic_old = Electrostatics::BasicEnergy(&system, charges.data(), box_size, epsilon_rel);
    auto r_new = system.world_to_rcoord[2] + gmx::RVec({0.2, 0.2, 0.2});
    auto old_cell = system.world_to_cell[2].first;
    auto new_cell = system.grid->coordToCell(r_new[XX], r_new[YY], r_new[ZZ]);
    system.updateParticle(2, r_new, cell_vectors[old_cell], cell_vectors[new_cell]);
    auto total_energy_basic_new = Electrostatics::BasicEnergy(&system, charges.data(), box_size, epsilon_rel);

    auto lj_new = Electrostatics::BasicEnergyLj(&system, As.data(), Bs.data(), box_size);

    /// Ewald
    auto total_energy_ewald_old = Electrostatics::EwaldEnergy(cell_vectors2 , &system2, charges.data(),
                                                              beta, box_size, r_cutoff, k_cutoff, epsilon_rel);

    EXPECT_TRUE(system2.world_to_rcoord[2] == system.world_to_rcoord[2] - gmx::RVec({0.2, 0.2, 0.2}));
    r_new = system2.world_to_rcoord[2] + gmx::RVec({0.2, 0.2, 0.2});
    old_cell = system2.world_to_cell[2].first;
    new_cell = system2.grid->coordToCell(r_new[XX], r_new[YY], r_new[ZZ]);
    system2.updateParticle(2, r_new, cell_vectors2[old_cell], cell_vectors2[new_cell]);
    auto total_energy_ewald_new = Electrostatics::EwaldEnergy(cell_vectors2 , &system2, charges.data(),
                                                              beta, box_size, r_cutoff, k_cutoff, epsilon_rel);
    real delta_lj = lj_new - lj_old ;
    std::cout << total_energy_basic_old + lj_old<< " " <<  total_energy_ewald_old + lj_old << "\n";
    std::cout << total_energy_basic_new - total_energy_basic_old + lj_new - lj_old
              << " ewald: " << total_energy_ewald_new - total_energy_ewald_old + lj_new - lj_old << "\n";
    std::cout << "rel. error = " << (1 - (total_energy_basic_new - total_energy_basic_old + delta_lj)/(total_energy_ewald_new - total_energy_ewald_old + delta_lj))*100 << " %\n";
}
