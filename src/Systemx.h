//
// Created by smutekj on 14.05.23.
//
#pragma once
#include <vector>
#include <unordered_set>
#include <fstream>
#include "vectypes.h"
#include <random>
#include <nlohmann/json.hpp>




constexpr int N_MOLECULE_TYPES = 2;
constexpr int N_ATOM_TYPES = 2;
#define MAX_ATOMS_IN_CELL 100
#define DIM 3


template<typename Type>
class CellVector{

public:
    Type r[MAX_ATOMS_IN_CELL];
    size_t atoms_in_cell = 0;

    CellVector(){}

    void insertParticle(const Type& x){
        atoms_in_cell++;
        r[atoms_in_cell-1] = x;
        // assert(MAX_ATOMS_IN_CELL*DIM > 3*atoms_in_cell
    }

    void insertParticles(const gmx::RVec* x, size_t* indices_in_cell, const size_t n_inserted){

        for(size_t i = 0; i<n_inserted; ++i){
            r[atoms_in_cell] = x[i];
            indices_in_cell[i] = atoms_in_cell;
            atoms_in_cell++;
        }
    }

    void removeParticle(size_t index_in_cell){
        assert(atoms_in_cell>0);
        r[index_in_cell] = r[atoms_in_cell-1];
        atoms_in_cell--;
    }

    void removeParticles(const size_t n_removed, size_t* indices_in_cell){
//        assert(atoms_in_cell>0);


        for(int i = 0; i<n_removed; ++i) {
            r[indices_in_cell[i]] = r[atoms_in_cell-1];
            atoms_in_cell--;
            indices_in_cell[i]  = atoms_in_cell; // tell that to whoever called the function
        }
    }
};

template<typename Type>
class CellVectorD{

public:
    std::vector<Type> r;
    size_t atoms_in_cell = 0;

    CellVectorD(){}

    void insertParticle(const Type& x){
        atoms_in_cell++;
        r.push_back(x);
        // assert(MAX_ATOMS_IN_CELL*DIM > 3*atoms_in_cell
    }

    void insertParticles(const gmx::RVec* x, size_t* indices_in_cell, const size_t n_inserted){

        for(size_t i = 0; i<n_inserted; ++i){
            r[atoms_in_cell] = x[i];
            indices_in_cell[i] = atoms_in_cell;
            atoms_in_cell++;
        }
    }

    void removeParticle(size_t index_in_cell){
        assert(atoms_in_cell>0);
        r[index_in_cell] = r[atoms_in_cell-1];
        r.pop_back();
        atoms_in_cell--;
    }

    void removeParticles(const size_t n_removed, size_t* indices_in_cell){
//        assert(atoms_in_cell>0);


        for(int i = 0; i<n_removed; ++i) {
            r[indices_in_cell[i]] = r[atoms_in_cell-1];
            atoms_in_cell--;
            indices_in_cell[i]  = atoms_in_cell; // tell that to whoever called the function
        }
    }
};



struct LJInteractionParams{
    real A;
    real B;
    real cutoff;
};


template <class DataType>
class SetToWorldMap {

    struct SetIndex {
        size_t set_index = -1;
        size_t index_in_set = -1;
    };

    std::vector<SetIndex> world_to;
    std::vector<std::vector<DataType>> set_index_to_data;

public:
    SetToWorldMap(size_t max_world, size_t max_set_indices) : world_to(max_world), set_index_to_data(max_set_indices){}
    SetToWorldMap() = default;

    void add(size_t set_index, DataType datum){
        set_index_to_data.at(set_index).push_back(datum);
        world_to.at(datum) = {set_index, set_index_to_data.at(set_index).size()-1} ;
    }

    void remove(size_t world_index){
        const auto removed = world_to.at(world_index);
        const auto removed_set_index = removed.set_index;
        const auto removed_index_in_set = removed.index_in_set;
        assert(set_index_to_data.at(removed_set_index).at(removed_index_in_set).world_index  == world_index);
        const auto index_of_last_in_removed_set = set_index_to_data.at(removed_set_index).size()-1;
        const auto world_of_last_in_set = set_index_to_data.at(removed_set_index).at(index_of_last_in_removed_set).world_index;

        set_index_to_data.at(removed_set_index).at(removed_index_in_set) = *(set_index_to_data.at(removed_set_index).end()-1);
        set_index_to_data.at(removed_set_index).pop_back();

        world_to.at(world_of_last_in_set) = removed;
    }

    [[nodiscard]] const std::vector<DataType>& at(size_t set_index) const{
        return set_index_to_data.at(set_index);
    }

    [[nodiscard]] const DataType& getDatum(size_t world_index) const{
        const auto set_index = world_to.at(world_index).set_index;
        const auto index_in_set = world_to.at(world_index).index_in_set;
        return set_index_to_data.at(set_index).at(index_in_set);
    }
};


struct Molecule{
    size_t  first_;
    size_t last_;
    size_t world_index_;
    Molecule(size_t first, size_t last, size_t world_index) :
            first_(first), last_(last), world_index_(world_index)
    {}
    Molecule() = default;
};

class SetToMoleculeMap {

    struct SetIndex {
        size_t set_index;
        size_t index_in_set;
    };

    std::vector<SetIndex> world_to;
    std::vector<std::vector<Molecule>> set_index_to_data;

public:
    SetToMoleculeMap(size_t max_world, size_t max_set_indices) : world_to(max_world), set_index_to_data(max_set_indices){
    }
    SetToMoleculeMap() = default;

    void add(size_t set_index, Molecule mol){
        set_index_to_data.at(set_index).push_back(mol);
        world_to.at(mol.world_index_) = {set_index, set_index_to_data.at(set_index).size()-1} ;
    }

    void remove(Molecule mol){
        const auto removed = world_to.at(mol.world_index_);
        const auto removed_set_index = removed.set_index;
        const auto removed_index_in_set = removed.index_in_set;
        assert(set_index_to_data.at(removed_set_index).at(removed_index_in_set).world_index_  == mol.world_index_);
        const auto index_of_last_in_removed_set = set_index_to_data.at(removed_set_index).size()-1;
        const auto world_of_last_in_set = set_index_to_data.at(removed_set_index).at(index_of_last_in_removed_set).world_index_;

        set_index_to_data.at(removed_set_index).at(removed_index_in_set) = *(set_index_to_data.at(removed_set_index).end()-1);
        set_index_to_data.at(removed_set_index).pop_back();

        world_to.at(world_of_last_in_set) = removed;
    }

    [[nodiscard]] const std::vector<Molecule>& getMoleculesOf(size_t mol_type_index) const{
        return set_index_to_data.at(mol_type_index);
    }

    [[nodiscard]] const Molecule& getMoleculeAt(size_t world_index) const{
        const auto set_index = world_to.at(world_index).set_index;
        const auto index_in_set = world_to.at(world_index).index_in_set;
        return set_index_to_data.at(set_index).at(index_in_set);
    }
};

template class CellVector<gmx::RVec>;

class Grid;
class SearchGrid;

struct CellIndex{
    size_t of_cell;
    size_t in_cell;
};

struct LJInteractionTable{

    LJInteractionParams  pair_coeffs[N_ATOM_TYPES][N_ATOM_TYPES];
    real  LJ6_0s[N_ATOM_TYPES][N_ATOM_TYPES];
    real  LJ12_0s[N_ATOM_TYPES][N_ATOM_TYPES];

    size_t n_atom_types = N_ATOM_TYPES;

    LJInteractionTable() = default;

    explicit LJInteractionTable(std::vector<LJInteractionParams>& atom_type_to_lj_params)
        : n_atom_types(atom_type_to_lj_params.size())
        {
//            pair_coeffs.resize(n_atom_types*n_atom_types);
            size_t indexi = 0;
            for(auto lj_params_i : atom_type_to_lj_params){
                size_t indexj=0;
                for(auto lj_params_j : atom_type_to_lj_params){
                    const auto sig_i = std::pow(lj_params_i.A/lj_params_i.B, 1./6.);
                    const auto sig_j = std::pow(lj_params_j.A/lj_params_j.B, 1./6.);
                    const auto eps_i = lj_params_i.B/std::pow(sig_i, 6.) / 4.;
                    const auto eps_j = lj_params_j.B/std::pow(sig_j, 6.) / 4.;
                    const auto sig_ij = (sig_i + sig_j)/2.;
                    const auto eps_ij = std::sqrt(eps_i*eps_j);
                    pair_coeffs[indexi][indexj].A = 4*eps_ij*std::pow(sig_ij, 12.);
                    pair_coeffs[indexi][indexj].B = 4*eps_ij*std::pow(sig_ij, 6.);
                    indexj++;
                    const auto r_cutoff_0 = lj_params_i.cutoff;
                    LJ12_0s[indexi][indexj] = 0*4.*eps_ij*std::pow(sig_ij/r_cutoff_0, 12.);
                    LJ6_0s[indexi][indexj] = -0*4.*eps_ij*std::pow(sig_ij/r_cutoff_0, 6.);
                }
                indexi++;
            }
    }

    [[nodiscard]] const real& get_coeff_A(size_t atom_type_i, size_t atom_type_j) const{
        return pair_coeffs[atom_type_i][atom_type_j].A;
    }
    [[nodiscard]] const real& get_coeff_B(size_t atom_type_i, size_t atom_type_j) const{
        return pair_coeffs[atom_type_i][atom_type_j].B;
    }
    [[nodiscard]] const real& get_LJ12_0(size_t atom_type_i, size_t atom_type_j) const{
        return LJ12_0s[atom_type_i][atom_type_j];
    }
    [[nodiscard]] const real& get_LJ6_0(size_t atom_type_i, size_t atom_type_j) const{
        return LJ6_0s[atom_type_i][atom_type_j];
    }
    [[nodiscard]] real get_sig(size_t atom_type_i, size_t atom_type_j) const{
        const auto A = pair_coeffs[atom_type_i][atom_type_j].A;
        const auto B = pair_coeffs[atom_type_i][atom_type_j].B;
        return std::pow(A/B, 1./6.);
    }
    [[nodiscard]] real get_eps(size_t atom_type_i, size_t atom_type_j) const{
        const auto A = pair_coeffs[atom_type_i][atom_type_j].A;
        const auto B = pair_coeffs[atom_type_i][atom_type_j].B;
        return B*B/(4*A);
    }
};

struct AtomIndex{
    size_t atom_type;
    size_t index_in_atom;
};


class Particle;

struct ParticleData{

    Particle* particle = nullptr;
    size_t atom_type = -1;
    size_t molecule_world_index = 0;
    real charge = 0;
    real radius = 3;
    bool is_active=false;
};

struct MoleculeTypeData{

    size_t molecule_type;
    std::vector<gmx::RVec> structure;
};


#include <iostream>
class World {

public:
    std::vector<ParticleData> all_particles; //! active and inactive particles
    std::vector<int> to_active_particles;

    std::vector<int> to_inactive_particles;
    std::vector<size_t> active_particles;

    std::vector<size_t> inactive_particles;
    std::vector<Molecule> all_molecules;
    std::vector<int> to_active_molecules; //! from mol-world to indices of active molecules in mol-world

    std::vector<int> to_inactive_molecules;
    std::vector<size_t> active_molecules; //! points to all_molecules

    std::vector<size_t> inactive_molecules;
    //! mol/atom_type specific maps
    SetToMoleculeMap molecule_type_to_active_molecules;
    SetToWorldMap<size_t> atom_type_to_active_atom_indices;

    //! mol/atom_type specific maps
    std::vector<std::vector<size_t>> molecule_type_to_inactive_molecules;

    World() = default;

    void init_world(const size_t n_max_atoms,
                    const size_t n_max_molecules,  const std::vector<size_t>& n_max_of_mol_type) {

        molecule_type_to_active_molecules = SetToMoleculeMap(n_max_molecules, N_MOLECULE_TYPES);
        molecule_type_to_inactive_molecules.resize(N_MOLECULE_TYPES);
        atom_type_to_active_atom_indices = SetToWorldMap<size_t>(n_max_atoms, N_ATOM_TYPES);

  /*      all_particles.resize(n_max_atoms);
        all_molecules.resize(n_max_molecules);*/


        std::vector<size_t> all_inactive_indices(n_max_atoms);
        std::iota(all_inactive_indices.begin(), all_inactive_indices.end(), 0);
        inactive_particles.insert(inactive_particles.end(),
                                    all_inactive_indices.begin(), all_inactive_indices.end());
        to_inactive_particles.insert(to_inactive_particles.end(),
                                  all_inactive_indices.begin(), all_inactive_indices.end());
        to_active_particles.resize(n_max_atoms, -1);


        all_inactive_indices.resize(n_max_molecules);
        std::iota(all_inactive_indices.begin(), all_inactive_indices.end(), 0);
        inactive_molecules.insert(inactive_molecules.end(),
                                  all_inactive_indices.begin(), all_inactive_indices.end());
        to_inactive_molecules.insert(to_inactive_molecules.end(),
                                  all_inactive_indices.begin(), all_inactive_indices.end());
        to_active_molecules.resize(n_max_molecules, -1);

        size_t mol_type = 0;
        size_t n_max_last = 0;

        for(auto n_max : n_max_of_mol_type){
            auto& inactive_mol_indices = molecule_type_to_inactive_molecules.at(mol_type);

            all_inactive_indices.resize(n_max);
            std::iota(all_inactive_indices.begin(), all_inactive_indices.end(), n_max_last);
            inactive_mol_indices.insert(inactive_mol_indices.end(),
                                     all_inactive_indices.begin(), all_inactive_indices.end());
//            std::make_heap(inactive_mol_indices.begin(), inactive_mol_indices.end());

                mol_type++;
            n_max_last+=n_max;
        }

    }


    [[nodiscard]] const gmx::RVec& coordsOfAtom(size_t world_index) const;

    [[nodiscard]] const Particle& getParticle(size_t world_index) const{
        return *all_particles.at(world_index).particle;
    }

    [[nodiscard]] const Particle* getParticlep(size_t world_index) const{
        return all_particles.at(world_index).particle;
    }

    void activateMolecule(size_t molecule_type, Molecule molecule){
        assert(isInactiveMolecule(molecule.world_index_));
        for(auto world_index = molecule.first_; world_index<=molecule.last_; ++world_index){
            activateParticle(world_index);
        }
        const auto activated_index_in_inactive = to_inactive_molecules.at(molecule.world_index_);
        const auto last_in_inactive = *(inactive_molecules.end()-1);
        assert(inactive_molecules.at(activated_index_in_inactive) == molecule.world_index_);


        to_inactive_molecules.at(last_in_inactive) = activated_index_in_inactive;
        inactive_molecules.at(activated_index_in_inactive) = last_in_inactive;
        inactive_molecules.pop_back();
        active_molecules.push_back(molecule.world_index_);

        to_inactive_molecules.at(molecule.world_index_) = -1;
        to_active_molecules.at(molecule.world_index_) = static_cast<int>(active_molecules.size())-1;

        molecule_type_to_active_molecules.add(molecule_type, molecule);

        auto x = std::find(molecule_type_to_inactive_molecules[molecule_type].begin(), molecule_type_to_inactive_molecules[molecule_type].end(), molecule.world_index_);
        molecule_type_to_inactive_molecules[molecule_type].erase(x);
    }

    void deactivateMolecule(size_t molecule_type, Molecule molecule){
        for(auto world_index = molecule.first_; world_index<=molecule.last_; ++world_index){
            deactivateParticle(world_index);
        }
        assert(isActiveMolecule(molecule.world_index_));
        const auto deactivated_index_in_active = to_active_molecules.at(molecule.world_index_);
        const auto last_in_active = *(active_molecules.end()-1);
        assert(active_molecules.at(deactivated_index_in_active) == molecule.world_index_);

        to_active_molecules.at(last_in_active) = deactivated_index_in_active;
        active_molecules.at(deactivated_index_in_active) = last_in_active;
        active_molecules.pop_back();
        inactive_molecules.push_back(molecule.world_index_);

        to_active_molecules.at(molecule.world_index_) = -1;
        to_inactive_molecules.at(molecule.world_index_) = static_cast<int>(inactive_molecules.size())-1;

        molecule_type_to_active_molecules.remove(molecule);
        if(molecule.world_index_ ==6091){
            std::cout <<"pica\n";
        }
        molecule_type_to_inactive_molecules[molecule_type].push_back(molecule.world_index_);
    }

    void activateParticle(size_t world_index);
    void deactivateParticle(size_t world_index);

    Molecule randomActiveMolecule(size_t random_int) const{
        return all_molecules.at(active_molecules.at(random_int));
    }

    [[nodiscard]] const std::vector<size_t>& getAllActiveParticleIndices() const{
        return active_particles;
    }

    Molecule randomInactiveMoleculeOf(size_t mol_type_index) const{
        return all_molecules.at(molecule_type_to_inactive_molecules.at(mol_type_index).at(0));
    }

    bool isActiveParticle(size_t world_index) const{
        return to_active_particles.at(world_index) != -1;
    }

    bool isActiveMolecule(size_t mol_world_index) const{
        return to_active_molecules.at(mol_world_index) != -1;
    }

    bool isInactiveMolecule(size_t mol_world_index) const{
        return !isActiveMolecule(mol_world_index);
    }

    [[nodiscard]] std::vector<Particle> getCopyOfMoleculeCoords(size_t mol_world_index) const;
    [[nodiscard]] const std::vector<Molecule>& getActiveMoleculesOf(size_t mol_type)const;
    [[nodiscard]] const std::vector<size_t>& getActiveAtomsOf(size_t atom_type)const;
};

class Systemx{

public:

    World world;

    std::map<size_t, std::string> atom_type_to_name;
    std::map<std::string, size_t> name_to_atom_type;

    std::map<size_t, std::string> molecule_type_to_name;
    std::map<std::string, size_t> name_to_molecule_type;

    std::vector<std::vector<Particle>> molecule_type_to_structure;


    std::vector<std::pair<size_t, size_t>> world_to_cell ;
    std::vector<std::vector<size_t>> cell_to_world ;

    std::vector<LJInteractionParams> atom_type_to_lj_params;
    std::vector<LJInteractionParams> world_to_lj_params;

    std::vector<size_t> gc_molecule_types;

    static constexpr size_t n_molecule_types = N_MOLECULE_TYPES;
    static constexpr size_t n_atom_types = N_ATOM_TYPES;

    std::vector<real> chemical_potentials;
    real beta;
    real p;
    real cutoff;

    LJInteractionTable lj_table;

    std::vector<CellVector<Particle>> cell_vectors2_;
    SearchGrid* grid;

    Systemx(SearchGrid* grid_, const nlohmann::json& input);

    void updateMolecule(Molecule updated_molecule, std::vector<Particle>& new_particles);
    void insertMolecule(const std::vector<Particle>& new_particles, Molecule new_molecule, size_t new_molecule_type);
    void insertInactiveMolecule(const std::vector<Particle>& new_particles, Molecule new_molecule, size_t new_molecule_type);
    void removeMolecule(Molecule, size_t);

    [[nodiscard]] const std::vector<size_t>& getAllActiveAtomIndices() const;
    [[nodiscard]] const std::vector<size_t>& getAllActiveMolecules() const;
    [[nodiscard]] const std::vector<Molecule> & getActiveMolecules(size_t mol_type) const;
    [[nodiscard]] const std::vector<size_t>& getActiveAtomsOf(size_t atom_type) const;

    [[nodiscard]] bool isInSameMolecule(size_t world_index_i, size_t world_index_j) const{
        auto b = world.all_particles.at(world_index_i).molecule_world_index
                == world.all_particles.at(world_index_j).molecule_world_index;
        return b;
    }

    void removeFromCell(size_t removed_world_index);

    void scaleAllCoordinates(real scale);
    void scaleAllCoordinates2(real scale);
};



gmx::RVec calcCom(std::vector<Particle>& particles);

void translateMolecule(std::vector<Particle>& particles, gmx::RVec trans, gmx::RVec box_size);

gmx::RVec generateRandomPosition(std::vector<Particle>& particles, gmx::RVec box_size,
                                 std::uniform_real_distribution<real>& uniform_dist, std::mt19937& e1);




