#pragma once
#include <sstream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <nlohmann/json.hpp>

class System {
    typedef size_t ParticleType;
    typedef size_t MoleculeType;

    std::vector<std::string> molecule_names;
    std::vector<std::string> atom_names;
    std::vector<size_t> num_of_molecules;
    std::unordered_map<std::string, std::vector<size_t>> resname2atominds;
    std::unordered_map<std::string, std::vector<size_t>> atomname2atominds;
    std::unordered_map<std::string, std::vector<std::string>> resname2atomnames;

public:
    size_t n_atoms;

    const std::vector<size_t>& indicesOfMolecule(std::string& molname) const;
    const std::vector<size_t>& indicesOfAtom(std::string& atomname) const;

        explicit System(const std::string& file_name, const nlohmann::json& input);
    explicit System(const nlohmann::json& input);
};

