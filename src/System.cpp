#include "System.h"
#include <iostream>

const std::vector<size_t>& System::indicesOfMolecule(std::string& molname) const{
    return resname2atominds.at(molname);
}

const std::vector<size_t>& System::indicesOfAtom(std::string& atomname) const{
    return atomname2atominds.at(atomname);
}

System::System(const std::string& file_name, const nlohmann::json& molecules_j) {

    n_atoms = 0;
    for(const auto& molecule : molecules_j){
        const auto& molecule_name = molecule.value("name", "");
        const auto n_molecules = molecule.value("count", 0);

        for(int molecule_index = 0; molecule_index<n_molecules; ++molecule_index) {
            for (const auto &atom_name: molecule.at("atoms")) {
                const auto name = atom_name.get<std::string>();
                resname2atomnames[molecule_name].push_back(name);
                resname2atominds[molecule_name].push_back(n_atoms);
                atomname2atominds[name].push_back(n_atoms);
                n_atoms++;
            }
        }
    }
    return;
    const auto extension_pos = file_name.find_last_of('.');
    if(file_name.substr(extension_pos) == ".gro"){
        std::ifstream file;
        std::string line;

        file.open(file_name);
        if (!file.is_open()) {
            throw std::runtime_error(file_name + " not found or inaccessible!");
        }

        std::getline(file, line); // this line is a comment with TITLE
        std::getline(file, line); // this line contains number of atoms
        n_atoms = std::stoi(line);

        std::string coordinate, residue_name;
        size_t i = 0;
        while (i<n_atoms) {
            std::getline(file, line);
            std::stringstream ss(line);

            std::string res;
            ss >> res;
            auto res_num_pos = res.find_first_not_of("0123456789");
            auto residue_number = std::stoi(res.substr(0, res_num_pos));
            residue_name = res.substr(res_num_pos);

            atom_names = resname2atomnames.at(residue_name);
            std::cout << residue_name << "\n";
            if(residue_name =="C3"){
                std::cout << "ola" << "\n";
            }
            for(const auto& atom_name : resname2atomnames.at(residue_name)){

                std::string atom_name_in_gro;
                ss >> atom_name_in_gro;
                std::cout << atom_name_in_gro << "\n";

                if(atom_name_in_gro.find(atom_name) == std::string::npos){
                    throw std::runtime_error("error reading .gro file, supplied atom_name in json: " + atom_name +
                                             " does not correspond to atom name in .gro file: " + atom_name_in_gro);
                }

                atomname2atominds[atom_name].push_back(i);
                resname2atominds[residue_name].push_back(i);
                i++;
                if(atom_name!=atom_names.back()){
                    std::getline(file, line);
                    ss.str(line);
                    ss >> res;
                }
            }
        }
        file.close();
    } else {
        throw std::runtime_error(file_name + " does not have a .gro extension!");
    }
}
System::System(const nlohmann::json& input) : System(input.value("filename", ""), input["molecules"])
{}
