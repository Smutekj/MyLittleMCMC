//
// Created by smutekj on 23.04.23.
//

#ifndef GR_GCMC_CLION_ENERGY_H
#define GR_GCMC_CLION_ENERGY_H


#include <vector>
#include <memory>
#include "vectypes.h"


#include <unordered_set>

constexpr real dr = 0.2;
constexpr real r_max = 25;

class Grid;
class VerletList;
class PotentialData;
class VerletListSpecial;
template <typename Type>
    class CellVector;
class Systemx;


struct LJComponents{
    real L6;
    real L12;
};


class Particle;
class Molecule;


class Energy {

public:
    real total_energy = 0.0;

    std::vector<real> particle_energies;
    std::vector<size_t> all_indices;
    std::vector<real> LJ_table;
    real sum_of_energy_changed = 0.0;

    LJComponents total_lj_energy;

    std::vector<std::vector<LJComponents>> total_lj_energy_comp;

public:

        explicit Energy(size_t n_atoms);


    real inline LJTable(const real r){
        return LJ_table[std::floor(r/dr)];
    }
    real calcVirialBasic(const Systemx &system, gmx::RVec box_length,
                                 real cutoff_sq);
    real calcVirial(const Systemx& system, gmx::RVec box_length,  real cutoff_sq) const;

    [[nodiscard]] real calcTotalEnergyBasic(const Systemx& system,
                                            gmx::RVec box_length,
                                            real cutoff_sq)const;

    [[nodiscard]] real calcNewMoleculeEnergy(  const Systemx& system,
                                         Molecule new_molecule,
                                         const std::vector<Particle>& new_particles,
                                         gmx::RVec box_length,
                                         real cutoff_sq);

    [[nodiscard]] real calcParticlesEnergy(  const Systemx& system,
                                       const std::vector<Particle>& particles,
                                             gmx::RVec box_length,
                                       real cutoff_sq) const;

    LJComponents calcMoleculeEnergy(  const Systemx& system,
                                      Molecule molecule,
                                      gmx::RVec box_length,
                                      real cutoff_sq) const;
    real calcTotalEnergy(  const Systemx& system,
                                   gmx::RVec box_length,
                                   real cutoff_sq) ;

    real calMoleculeMoveEnergyChange(  const Systemx& system,
                                               Molecule molecule,
                                               const std::vector<Particle>& particles,
                                               gmx::RVec box_length,
                                               real cutoff_sq);

    real calcNewParticleEnergyVLNoGather(  const std::vector<size_t>& neighbouring_indices, size_t changed_index,
                                               const gmx::RVec* r_coords,
                                               const gmx::RVec& r_coords_new,
                                               real box_length,
                                               real A, real B, real cutoff_sq);
    real calcNewParticleEnergyVL(VerletList& verlet_list, size_t changed_index,
                                   const gmx::RVec* r_coords,
                                   const gmx::RVec& r_coords_new,
                                   real box_length,
                                   real A, real B, real cutoff_sq);

    };


constexpr real siga = 2;
constexpr real sig6 = siga*siga*siga*siga*siga*siga;
constexpr real A0 = 4*1*(sig6*sig6);
constexpr real B0 = 4*1*(sig6);

real inline LJ2(const real r_m2){
    const real r_m6 = r_m2*r_m2*r_m2;
    return r_m6*(A0*r_m6 - B0);
}

real inline LJ(const real r2, const real A, const real B, const real cutoff_sq){
    const real r6 = 1/(r2*r2*r2);
    return (A*r6*r6 - B*r6)*(r2 <= cutoff_sq);
}


real evaluatePotentials(const std::vector<PotentialData>& potential_data, real cutoff_sq );



real inline LJ6(const real r6, const real B){
    return (B*r6);
}

real inline LJ12(const real r6, const real A){
    return (A*r6*r6);
}

gmx::RVec inline LJ_force(const gmx::RVec delta_r,
                          const real A, const real B, const real cutoff_sq){

    const auto r = gmx::norm(delta_r);
    const auto r2 = gmx::norm2(delta_r);
    const real r7 = 1/(r*r2*r2*r2);
    const real r13 = r7*r7*r;
    const real force_norm = -(-12.0*A*r13 + 6.0*B*r7)*(r2 <= cutoff_sq);

    return gmx::unitVector(delta_r) * force_norm;
}


real energyCorrection(const Systemx& system, const Energy& energizer,
                      const real volume,
                      const real cutoff);

real pressureCorrection(const Systemx& system, const Energy& energizer,
                        const gmx::RVec box_size, const real cutoff);



#endif //GR_GCMC_CLION_ENERGY_H



