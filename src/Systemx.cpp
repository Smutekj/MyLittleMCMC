//
// Created by smutekj on 14.05.23.
//
#include "Systemx.h"
#include "Grid.h"
#include <random>
#include "core.h"


std::vector<std::string> separateLine(std::string& line, char delimiter) {
    size_t pos;
    size_t start_pos = line.find_first_not_of(' ', 0);
    std::vector<std::string> separated_line;

    while ((pos = line.find(delimiter, start_pos)) != std::string::npos) {
        separated_line.push_back(line.substr(start_pos, pos-start_pos));
        start_pos = line.find_first_not_of(delimiter, pos);
    }
    separated_line.push_back(line.substr(start_pos, line.size()));
    return separated_line;
};


std::vector<Particle> readGroStructure(std::string gro_file_name, Systemx& system){

    std::ifstream struct_gro_file(gro_file_name);

    std::cout << "Reading molecule structure in file: " << gro_file_name << "\n";

    size_t n_atoms_in_molecule = 0;
    std::string line;
    getline(struct_gro_file ,line);
    getline(struct_gro_file, line);
    n_atoms_in_molecule = std::stoi(line);
    std::vector<Particle> particles(n_atoms_in_molecule);
    for(int atom_i = 0; atom_i < n_atoms_in_molecule; ++atom_i){

        getline(struct_gro_file, line);
        const auto seperated_line = separateLine(line, ' ');
        const gmx::RVec r = {std::stof(seperated_line[3]), std::stof(seperated_line[4]), std::stof(seperated_line[5])};
        particles.at(atom_i).r = r;

        const auto atom_name = std::string(seperated_line[1]);

        if(system.name_to_atom_type.count(atom_name) == 0){
            throw std::runtime_error("atom name : " + atom_name  + " not in name_to_atom_type!");
        }
        particles.at(atom_i).atom_type = system.name_to_atom_type.at(atom_name);

    }
    return particles;
}

void rotateMolecule(std::vector<Particle>& particles, gmx::RVec trans, Systemx& system){

}

gmx::RVec calcCom(Molecule molecule, Systemx& system){
    gmx::RVec com = {0,0,0};
    for(size_t world_index = molecule.first_; world_index<=molecule.last_; ++world_index){
        const auto& particle = system.world.getParticle(world_index);
        real n = static_cast<real>(world_index - molecule.first_);
        com = (n*com + particle.r)/(n+1);
    }
    return com;
}


void readJsonAndInsertMoleculesIntoSystem(Systemx& system, const nlohmann::json& input) {

    const auto& topology = input.at("topology");

    //! Initialize atomTypes and interactions
    size_t atom_type_index = 0;
    std::vector<real> atom_type_radii;
    for (const auto &atom_type: topology.at("atomtypes")) {
        const auto &atom_type_name = atom_type.value("name", "");

        system.atom_type_to_name[atom_type_index] = atom_type_name;
        system.name_to_atom_type[atom_type_name] = atom_type_index;

        LJInteractionParams lj_params;
        const auto sig = static_cast<real>(atom_type.value("sig", 3.0));
        const auto eps = static_cast<real>(atom_type.value("eps", 1.0));
        lj_params.A = 4.0 * eps * std::pow(sig, 12.);
        lj_params.B = 4.0 * eps * std::pow(sig, 6.0);
        lj_params.cutoff = atom_type.value("cutoff", 10.0);
        system.atom_type_to_lj_params.push_back(lj_params);
        std::cout << "A =" << lj_params.A << "\n";
        std::cout << "B =" << lj_params.B << "\n";

        atom_type_radii.push_back(atom_type.value("radius", 3.0));
        atom_type_index++;
    }

    size_t mol_type_index = 0;
    size_t max_molecules = 0;
    size_t max_atoms = 0;
    std::vector<size_t> max_ns;
    //! Initialize MoleculeTypes
    for (const auto &molecule_type: topology.at("molecules")) {
        const auto molecule_type_name = molecule_type.value("name", "");
        const auto max_molecule_count = molecule_type.value("max_count", 1);

        assert(!molecule_type_name.empty());

        system.molecule_type_to_name[mol_type_index] = molecule_type_name;
        system.name_to_molecule_type[molecule_type_name] = mol_type_index;

        const auto gro_file = molecule_type.value("structure", "");
        system.molecule_type_to_structure.push_back(readGroStructure(gro_file, system));
        const auto n_atoms_in_molecule = (system.molecule_type_to_structure.end()-1)->size();

        auto first = max_atoms;
        for(int i = 0; i<max_molecule_count; ++i){
            Molecule new_molecule = {first ,first + n_atoms_in_molecule-1 ,i + max_molecules};
            system.world.all_molecules.push_back(new_molecule);
            first += n_atoms_in_molecule;
            for(int atom_i = 0; atom_i < n_atoms_in_molecule; atom_i++){
                ParticleData pd;
                auto& atom_type = system.molecule_type_to_structure.at(mol_type_index).at(atom_i).atom_type;
                pd.radius = atom_type_radii.at(atom_type);
                pd.atom_type = atom_type;
                pd.molecule_world_index = i + max_molecules;
                system.world.all_particles.push_back(pd);
            }
        }

        size_t n_atoms_in_max_molecule_count = max_molecule_count*n_atoms_in_molecule;
        max_ns.push_back(n_atoms_in_max_molecule_count);

        max_molecules += max_molecule_count;
        max_atoms += n_atoms_in_max_molecule_count;

        if(molecule_type.contains("mu")){
            system.gc_molecule_types.push_back(mol_type_index);
            system.chemical_potentials.push_back(molecule_type.value("mu", 0));
        } else{
            system.chemical_potentials.push_back(GMX_REAL_MAX);
        }
        mol_type_index++;
    }

    system.world.init_world(max_atoms, max_molecules, max_ns);
    system.world_to_cell.resize(max_atoms);
    system.cell_to_world.resize( system.grid->n_cells_[XX] * system.grid->n_cells_[YY] * system.grid->n_cells_[ZZ]);


    size_t n_molecules = 0;
    const auto *grid = system.grid;
    size_t n_molecules_tmp = 0;
    size_t at = 0;
    size_t first = 0; size_t last;


    std::uniform_real_distribution<real> uni_dist(0,1);
    std::random_device r;
    std::mt19937 e1(0*r());

    const auto box_size = grid->box_size_;

    //! Insert Molecules themselves
    size_t mol_world_index = 0;
    size_t world_index = 0;
    for (const auto &molecule: input.at("insertmolecules")) {
        const auto &molecule_type_name = molecule.value("name", "");
        const auto molecule_type = system.name_to_molecule_type.at(molecule_type_name);
        auto& mol_structure = system.molecule_type_to_structure.at(molecule_type);
        const auto n_atoms_in_molecule = mol_structure.size();

        const auto active_count = molecule.value("count", 50);
        const auto inactive_count = molecule.value("max_count", 1) - active_count;
        n_molecules += active_count;

        while (mol_world_index < n_molecules) {

            Molecule new_molecule = {world_index, world_index+n_atoms_in_molecule-1, mol_world_index};

            const auto new_com = generateRandomPosition(mol_structure, box_size, uni_dist, e1);
            auto cell_index = system.grid->coordToCell(new_com);
            const auto& nearest_cells = system.grid->cell2nearest_neighbours.at(cell_index);

            bool far_enough = true;
            for (auto cell_index_n  : nearest_cells) {
                const auto& n_particles_in_cell = system.cell_vectors2_.at(cell_index_n).atoms_in_cell;
                for(int jn = 0; jn < n_particles_in_cell; ++jn){
                    const auto world_index_n = system.cell_to_world.at(cell_index_n).at(jn);
                    auto particle_data = system.world.all_particles.at(world_index_n);
                    auto radius_sq = particle_data.radius*particle_data.radius;
                    far_enough = far_enough and
                                 dist_sq_pbc(particle_data.particle->r, new_com, box_size / 2., box_size) >  radius_sq;
                }
            }
            if (far_enough) {
                for(auto& particle : mol_structure){
                    const auto lj_params = system.atom_type_to_lj_params.at(particle.atom_type);
                    system.world_to_lj_params.push_back(lj_params);
                    world_index++;
                }
                system.insertMolecule(mol_structure, new_molecule, molecule_type);
                system.world.all_molecules.at(mol_world_index) = new_molecule;
                std::cout << "inserted: " << mol_world_index << " molecules of type: "<< molecule_type_name << "\n";
                mol_world_index++;
            }
        }
        //! We have to update the indices for all the inactive particles that were not inserted!
        n_molecules += inactive_count;
        mol_world_index += inactive_count;
        world_index += inactive_count*n_atoms_in_molecule;
    }
    system.lj_table = LJInteractionTable(system.atom_type_to_lj_params);
}



Systemx::Systemx(SearchGrid* grid_, const nlohmann::json& input) :
        cell_vectors2_( grid_->n_cells_[XX] * grid_->n_cells_[YY] * grid_->n_cells_[ZZ])
{
    grid = grid_;
    const auto n_cells = grid->n_cells_[XX] * grid->n_cells_[YY] * grid->n_cells_[ZZ];
    cell_to_world.resize(n_cells);

    readJsonAndInsertMoleculesIntoSystem(*this, input);

    for(auto molecules_of_type : world.all_molecules){

    }
}


void Systemx::insertMolecule(const std::vector<Particle>& new_particles, Molecule new_molecule, size_t new_molecule_type){
    size_t i = 0;
    for(const auto particle : new_particles){
        const auto new_cell_index = grid->coordToCell(particle.r);
        auto& cell_vector = cell_vectors2_.at(new_cell_index);
        cell_vector.insertParticle(particle);
        const auto world_index = new_molecule.first_ + i;
        assert(cell_vector.atoms_in_cell > 0);

        cell_to_world.at(new_cell_index).push_back(world_index);
        world_to_cell.at(world_index) = {new_cell_index, cell_vector.atoms_in_cell-1};

        world.all_particles.at(world_index).particle = &(cell_vector.r[cell_vector.atoms_in_cell-1]);
        world.all_particles.at(world_index).atom_type = particle.atom_type;

        i++;
    }

    world.activateMolecule(new_molecule_type, new_molecule);
    assert(world.molecule_type_to_active_molecules.getMoleculeAt(new_molecule.world_index_).world_index_ == new_molecule.world_index_);
}

void Systemx::insertInactiveMolecule(const std::vector<Particle>& new_particles, Molecule new_molecule, size_t new_molecule_type){
    size_t i = 0;
    for(const auto particle : new_particles){
        const auto world_index = new_molecule.first_ + i;
        world.all_particles.at(world_index).particle->atom_type = particle.atom_type;
        i++;
    }

    world.activateMolecule(new_molecule_type, new_molecule);
}


void Systemx::removeMolecule(Molecule removed_molecule, size_t removed_molecule_type){
    const size_t n_atoms_in_molecule = removed_molecule.last_ - removed_molecule.first_ + 1;
    for(size_t world_index = removed_molecule.first_; world_index <= removed_molecule.last_; ++world_index){
        removeFromCell(world_index);
    }

    world.deactivateMolecule(removed_molecule_type, removed_molecule);
}

void Systemx::updateMolecule(Molecule updated_molecule, std::vector<Particle>& new_particles){
    for(size_t world_index = updated_molecule.first_; world_index<=updated_molecule.last_; ++world_index){
        const auto particle = new_particles.at(world_index-updated_molecule.first_);
        *world.all_particles.at(world_index).particle = particle;

        const auto new_cell_index = grid->coordToCell(particle.r);
        const auto old_cell_index = world_to_cell.at(world_index).first;
        const auto old_index_in_cell = world_to_cell.at(world_index).second;
        if(new_cell_index != old_cell_index){
            //! remove particle from old_cell
            const auto world_of_last_in_cell = *(cell_to_world.at(old_cell_index).end()-1);
            cell_vectors2_.at(old_cell_index).removeParticle(old_index_in_cell);
            cell_to_world.at(old_cell_index).at(old_index_in_cell) = world_of_last_in_cell;
            cell_to_world.at(old_cell_index).pop_back();
            world_to_cell.at(world_of_last_in_cell) = {old_cell_index, old_index_in_cell};
            world.all_particles.at(world_of_last_in_cell).particle = world.all_particles.at(world_index).particle;

            //! add particle into new_cell
            cell_vectors2_.at(new_cell_index).insertParticle(particle);
            cell_to_world.at(new_cell_index).push_back(world_index);
            world_to_cell.at(world_index) = {new_cell_index, cell_vectors2_.at(new_cell_index).atoms_in_cell-1};
            world.all_particles.at(world_index).particle = &(cell_vectors2_.at(new_cell_index).r[0]) + cell_vectors2_.at(new_cell_index).atoms_in_cell-1;
        }
    }

}

void Systemx::removeFromCell(size_t removed_world_index){


    const auto cell_of_removed_index = world_to_cell.at(removed_world_index).first;
    const auto index_in_cell_of_removed_index = world_to_cell.at(removed_world_index).second;

    auto& cell_vector = cell_vectors2_.at(cell_of_removed_index);
    const auto last_index_in_cell = cell_vector.atoms_in_cell - 1;

    assert(removed_world_index == cell_to_world.at(cell_of_removed_index).at(index_in_cell_of_removed_index));
//    assert(world.getParticlep(removed_world_index) == &(cell_vector.r[cell_vector.atoms_in_cell-1]));


    // Remove the particle from its CellVector
    cell_vector.removeParticle(index_in_cell_of_removed_index);

    // This moved the last particle in cell_vector at last_index_in_cell to index_in_cell_of_removed_index
    auto moved_index_in_world = cell_to_world.at(cell_of_removed_index).at(last_index_in_cell)  ;

    world.all_particles.at(moved_index_in_world) = world.all_particles.at(removed_world_index);

    // tell this to world_to_cell
    if(moved_index_in_world != removed_world_index){
        world_to_cell.at(moved_index_in_world) = {cell_of_removed_index, index_in_cell_of_removed_index};
    }
    // tell the same thing to cell_to_world
    cell_to_world.at(cell_of_removed_index).at(index_in_cell_of_removed_index) = moved_index_in_world;

    // finally remove the particle by shortening the cell_to_world
    cell_to_world.at(cell_of_removed_index).pop_back();
}

void World::activateParticle(size_t world_index){

    const auto activated_index_in_inactive = to_inactive_particles.at(world_index);
    const auto last_in_inactive = *(inactive_particles.end()-1);
    assert(activated_index_in_inactive != -1 );

    to_inactive_particles.at(last_in_inactive) = activated_index_in_inactive;
    inactive_particles.at(activated_index_in_inactive) =  last_in_inactive;
    inactive_particles.pop_back();
    active_particles.push_back(world_index);

    to_inactive_particles.at(world_index) = -1;
    to_active_particles.at(world_index) = static_cast<int>(active_particles.size())-1;

    const auto atom_type = all_particles.at(world_index).atom_type;
    atom_type_to_active_atom_indices.add(atom_type, world_index);
}

void World::deactivateParticle(size_t world_index){
    const auto deactivated_index_in_active = to_active_particles.at(world_index);
    const auto last_in_active = *(active_particles.end()-1);

    to_active_particles.at(last_in_active) = deactivated_index_in_active;
    active_particles.at(deactivated_index_in_active) =last_in_active;
    active_particles.pop_back();
    inactive_particles.push_back(world_index);

    to_active_particles.at(world_index) = -1;
    to_inactive_particles.at(world_index) = static_cast<int>(inactive_particles.size())-1;

    const auto atom_type = all_particles.at(world_index).atom_type;
//    atom_type_to_active_atom_indices.remove(atom_type);
}

[[nodiscard]] std::vector<Particle> World::getCopyOfMoleculeCoords(size_t mol_world_index) const
{
    const auto first =  all_molecules.at(mol_world_index).first_;
    const auto last =  all_molecules.at(mol_world_index).last_;
    const auto n_atoms_in_molecule = (last - first + 1);
    std::vector<Particle> particles_of_molecule(n_atoms_in_molecule);
    if(isActiveMolecule(mol_world_index)) {
        for (auto world_index = first; world_index <= last; ++world_index) {
            particles_of_molecule.at(world_index - first) = *all_particles.at(world_index).particle;
        }
    }

    return particles_of_molecule;
}

[[nodiscard]] const std::vector<size_t>& Systemx::getAllActiveAtomIndices() const{
    return world.getAllActiveParticleIndices();
}

[[nodiscard]] const std::vector<size_t> & Systemx::getAllActiveMolecules() const {
    return world.active_molecules;
}

[[nodiscard]] const std::vector<Molecule> & Systemx::getActiveMolecules(size_t mol_type) const {
    return world.getActiveMoleculesOf(mol_type);
}

[[nodiscard]] const std::vector<size_t>& Systemx::getActiveAtomsOf(size_t atom_type) const{
    return world.getActiveAtomsOf(atom_type);
}

[[nodiscard]] const std::vector<Molecule>& World::getActiveMoleculesOf(size_t mol_type)const{
    return molecule_type_to_active_molecules.getMoleculesOf(mol_type);
}
[[nodiscard]] const std::vector<size_t>& World::getActiveAtomsOf(size_t atom_type)const{
    return atom_type_to_active_atom_indices.at(atom_type);
}


[[nodiscard]] const gmx::RVec& World::coordsOfAtom(size_t world_index) const{
    return all_particles.at(world_index).particle->r;
}

void Systemx::scaleAllCoordinates(real scale){

#pragma clang loop vectorize(enable)
    for(auto& cell_vector : cell_vectors2_){
        for(int i = 0; i<cell_vector.atoms_in_cell; ++i){
            cell_vector.r[i].r *= scale;
        }
    }
}

void Systemx::scaleAllCoordinates2(real scale){

#pragma clang loop vectorize(enable)
    for(auto& cell_vector : cell_vectors2_){
        for(int i = 0; i<cell_vector.atoms_in_cell; ++i){
            cell_vector.r[i].r /= scale;
        }
    }
}



gmx::RVec calcCom(std::vector<Particle>& particles){
    gmx::RVec com = {0,0,0};
    real n = 0;
    for(auto& particle : particles){
        com = (n*com + particle.r)/(n+1);
        n++;
    }
    return com;
}

void translateMolecule(std::vector<Particle>& particles, gmx::RVec trans, gmx::RVec box_size){
    for(auto& particle : particles){
        particle.r += trans;
        particle.r[XX] += -(box_size[XX]) * (particle.r[XX] >= box_size[XX]) + (box_size[XX]) * (particle.r[XX] < 0);
        particle.r[YY] += -(box_size[YY]) * (particle.r[YY] >= box_size[YY]) + (box_size[YY]) * (particle.r[YY] < 0);
        particle.r[ZZ] += -(box_size[ZZ]) * (particle.r[ZZ] >= box_size[ZZ]) + (box_size[ZZ]) * (particle.r[ZZ] < 0);
    }
}

gmx::RVec generateRandomPosition(std::vector<Particle>& particles, gmx::RVec box_size,
                                 std::uniform_real_distribution<real>& uniform_dist, std::mt19937& e1){
    const auto new_com = gmx::RVec(box_size[XX] * uniform_dist(e1),
                                   box_size[YY] * uniform_dist(e1),
                                   box_size[ZZ] * uniform_dist(e1));

    const auto mol_com = calcCom(particles);
    translateMolecule(particles, new_com - mol_com, box_size);
    return new_com;
}


