#include <chrono>
#include <memory>
#include <algorithm>

#include "core.h"
#include "Readers.h"
#include "TrajectoryAnalysis.h"
#include "System.h"
#include <stdexcept>

#include "Energy.h"
#include "vectypes.h"
#include "Grid.h"
#include "NeighbourSearcherImplementation.h"
#include "VerletList.h"
#include "Systemx.h"

#include "Moves.h"

template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    auto size = static_cast<size_t>( size_s );
    std::unique_ptr<char[]> buf( new char[ size ] );
    std::snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

using json = nlohmann::json;


class Random{

    std::mt19937 gen;
    std::uniform_real_distribution<real> real_dist_;
    std::uniform_int_distribution<u_int32_t> int_dist_;

public:
    Random(size_t seed = std::rand()) : gen(seed), real_dist_(0, 1), int_dist_(0) {}

    auto rand(const real a = 0.0, const real b = 1.0){
        return a + (b-a)*real_dist_(gen);
    }
    auto randInt(){
        return int_dist_(gen);
    }

    void setNmax(size_t nmax){
        int_dist_ = std::uniform_int_distribution<u_int32_t>(nmax);
    }

};

std::unique_ptr<TrajectoryAnalysis> createAnalysis(const std::string& name, const json& entry,
                                                   std::shared_ptr<System> system){
    if (name== "rdf") {
        return std::make_unique<AnalysisRdf>(entry, system);
    }
    else if (name== "flucts") {
        return std::make_unique<NumberFluctsAnalysisOld>(entry, system);
    }
    else {
        throw std::runtime_error("invalid kezyword: " + json::string_t(entry) + " , it is not among supported analyses");
    }
}

std::unique_ptr<TrajectoryAnalysis> createAnalysis(const std::string& name, const json& entry,
                                                   const Systemx& system, Energy* energizer){
    if (name== "rdf") {
        return std::make_unique<AnalysisRdf>(entry, system);
    }
    else if (name== "widom") {
        return std::make_unique<AnalysisWidom>(entry, system, energizer);
    }
    else if (name== "density") {
        return std::make_unique<NumberFluctsAnalysis>(entry, system);
    }
    else if (name== "pressure") {
        return std::make_unique<Pressure>(entry, system, energizer);
    }
    else {
        throw std::runtime_error("invalid kezyword: " + json::string_t(entry) + " , it is not among supported analyses");
    }
}


std::tuple<const std::string&, const json&> jsonSingleItem(const json& j) {
    if (j.is_object() && j.size() == 1) {
        const auto it = j.cbegin();
        return {it.key(), it.value()};
    }
    throw std::runtime_error("invalid data: single key expected");
}

real relativeEnergyDrift(real current_total_energy, real sum_of_energy_changes, const real initial_energy) {
    real du = current_total_energy - initial_energy;
    if (std::isfinite(du)) {
        if (std::fabs(du) <= GMX_REAL_EPS) {
            return 0.0;
        }
        return (current_total_energy - (initial_energy + sum_of_energy_changes)) /
               (std::fabs(initial_energy) > GMX_REAL_EPS ? initial_energy : current_total_energy);
    }
    return std::numeric_limits<double>::quiet_NaN();
}

void test_simulation(size_t n_atoms, MovesSamplingData& move_sampling,
                         real cutoff,
                         gmx::RVec box_size, real T, real P,
                         long int n_steps, const long int n_equilibration_steps, const long int sampling_period,
                         const std::string& state_name,
                         real dx0, const nlohmann::json& input){

    const real R = 8.3145 / 1000. ; // kJ/K/mol
    const real beta = 1./(R*T);

    size_t n_trans = move_sampling.translation.n_attempts;
    size_t n_exc = move_sampling.exchange.n_attempts;
    size_t n_volume = move_sampling.volume.n_attempts;

    real r_buffer(cutoff);

    n_total =0;
    n_accepted = 0;

    std::shared_ptr<SearchGrid> grid = std::make_unique<SearchGrid>(box_size, gmx::RVec{cutoff, cutoff, cutoff});
    std::unique_ptr<GridSearch> grid_searcher = std::make_unique<GridSearch> ( box_size, gmx::RVec{cutoff, cutoff, cutoff}, cutoff, grid);
    std::unique_ptr<VerletList> verlet_list = std::make_unique<VerletList>(n_atoms, cutoff, r_buffer);

    Systemx system(grid.get(),  input.at("system"));
    system.beta = beta;
    system.cutoff = cutoff;

    std::random_device dev;
    std::mt19937 e1(0*dev());
    std::uniform_real_distribution<real> uniform_dist(0., 1.0);
    std::uniform_int_distribution<std::mt19937::result_type> uniform_int(1, n_exc+n_trans+n_volume);


    XTCWriter writer(state_name + ".xtc");
    TrajectoryFrame frame(n_atoms);

    real avg_pressure = 0.;
    real avg_v = 0.;
    real n_samples = 0.;


    Energy energizer(n_atoms);
    const auto n_atom_types = system.atom_type_to_lj_params.size();
    energizer.total_lj_energy_comp.resize(n_atom_types);
    for(int i = 0; i<n_atom_types; ++i){
        energizer.total_lj_energy_comp[i].resize(n_atom_types);
    }

    const real initial_energy = energizer.calcTotalEnergy(system, box_size, cutoff*cutoff) +
                                 energyCorrection(system, energizer, box_size[XX]*box_size[YY]*box_size[ZZ], cutoff);
    energizer.total_energy = initial_energy;

    //! Construct Analysis objects
    std::vector<std::unique_ptr<TrajectoryAnalysis>> analyzers;
    for( const auto& entry : input.at("analysis") ){
        const auto& [name, json_parameters] = jsonSingleItem(entry);
        analyzers.push_back(createAnalysis(name, json_parameters, system, &energizer));
    }

    auto start = std::chrono::high_resolution_clock::now();
    int step=0;

    //! MAIN MCMC CYCLE HERE
    while(step<n_steps) {
        const auto roll = uniform_dist(e1);
        const size_t rand_move = uniform_int(e1)%(n_trans + n_exc + n_volume);

        if(rand_move < n_trans) {

            translationMove(step,
//                             verlet_list.get(),
                              system, energizer,
                              box_size,
                              dx0, cutoff, move_sampling.translation,
                              e1, uniform_dist);
        } else if (rand_move < n_trans + n_volume){
            volumeMoveScaledCutoff(step,
                       system,
                       energizer,
//                       grid_searcher.get(),
                       box_size, cutoff,
                       beta, P, move_sampling.volume,
                       e1, uniform_dist);
        }
        else {
            if (roll < 0.5) {
                insertionMove(step, verlet_list.get(),
                              system, energizer,
                              n_atoms,  box_size,
                              cutoff, beta, move_sampling.exchange,
                              e1, uniform_dist);
            } else{
                    deletionMove(step, verlet_list.get(),
                                 system, energizer,
                                 n_atoms,  box_size,
                                 cutoff, beta, move_sampling.exchange,
                                 e1, uniform_dist);
            }
        }

        //!**** OUTPUT ****//
        if(step%(sampling_period) == 5){
            std::cout << "\n the simulation is in: " << 100*(step/(real)n_steps)<< " %\n";
            std::cout << "acceptance prob.: " << 100*(n_accepted/((real)n_total))<< " %\n";

            const auto remaining_fraction_of_time = (real)(n_steps - step-1)/(real)(step+1);
            auto calc_time = std::chrono::high_resolution_clock::now()- start;
            auto approx_remaining_time_m = std::chrono::duration_cast<std::chrono::minutes>(calc_time * remaining_fraction_of_time);
            auto approx_remaining_time_s = std::chrono::duration_cast<std::chrono::seconds>(calc_time * remaining_fraction_of_time);

            std::cout << "approx time left:  "  << approx_remaining_time_m.count() << "min "
                                                << approx_remaining_time_s.count() - 60*approx_remaining_time_m.count() << "s\n\n";
        }

        if(step % (sampling_period) == 0 and step >=  n_equilibration_steps) {
            n_samples = n_samples + 1.0;

            const real volume = box_size[XX] * box_size[YY] * box_size[ZZ];
            avg_v = (avg_v*(n_samples - 1) + volume)/n_samples;

            const auto total_energy_check2 = energizer.calcTotalEnergy(system, box_size, cutoff*cutoff);
            const auto a = energyCorrection(system, energizer, volume, cutoff);
            std::cout << "changed by : " << energizer.sum_of_energy_changed << "\n";
            std::cout << "\nENERGIES: " <<  energizer.sum_of_energy_changed + initial_energy  << "  " << total_energy_check2  + a<<"\n\n";
            std::cout << "By components: " << energizer.total_lj_energy.L12 + energizer.total_lj_energy.L6 << "\n";

            printf( "relative energy drift is: %4.15g\n" , relativeEnergyDrift(total_energy_check2+a, energizer.sum_of_energy_changed, initial_energy));
            printf("box_size = %.10f \n", box_size[XX]);
            printf("volume = %.10f \n", volume);

           /* const auto virial = energizer.calcVirial(system, box_size, cutoff*cutoff);
            real pressure_correction = pressureCorrection(system, energizer, box_size, cutoff);
            const auto pressure = (virial /volume + ((real)n_atoms * R * T) / volume)*1.0e30/6.022e23*1e3*1e-5 + pressure_correction;
            avg_pressure = (avg_pressure*(n_samples-1) + pressure)/n_samples;
            std::cout << "avg_pressure = " << avg_pressure << "\n";
            std::cout << "correction = " << pressure_correction << "\n";
            std::cout << "virial = " << virial /volume * 1.0e30/6.022e23*1e3*1e-5 <<"\n";*/

            const auto& all_active_atom_indices = system.getAllActiveAtomIndices();
            frame.reset(all_active_atom_indices.size());
            frame.number = step;
            frame.bounding_box[0] = box_size[XX];
            frame.bounding_box[1] = box_size[YY];
            frame.bounding_box[2] = box_size[ZZ];
            for(const auto active_index : all_active_atom_indices){
                frame.trajectory.x[active_index] = system.world.getParticle(active_index).r[XX];
                frame.trajectory.y[active_index] = system.world.getParticle(active_index).r[YY];
                frame.trajectory.z[active_index] = system.world.getParticle(active_index).r[ZZ];
            }
//            writer.writeNext(frame);

            //update analysis
            system.grid->box_size_ = box_size;
            system.cutoff = cutoff;
            for( auto& analyzer : analyzers ){
                analyzer->analyse(frame);
                if(auto* anal_widom = dynamic_cast<AnalysisWidom*>(analyzer.get()); anal_widom != nullptr  ){
                    const std::string& name =  system.molecule_type_to_name.at(anal_widom->inserted_molecule_type);
                    printf("widom  %s = %.10f \n", name.c_str(), anal_widom->widomEnergy());
                }
            }
        }
    }

    //! tell me how long the simulation took
    auto calc_time = std::chrono::high_resolution_clock::now()- start;
    auto simulation_time_m = std::chrono::duration_cast<std::chrono::minutes>(calc_time).count();
    auto simulation_time_s = std::chrono::duration_cast<std::chrono::seconds>(calc_time).count();
    std::cout << "\nThe simulation_took:  "  << simulation_time_m << "min "
              << simulation_time_s - 60*simulation_time_m << "s\n\n";

    // write analysis
    for( auto& analyzer : analyzers ){
        analyzer->writeOutput();
    }

    /// Write .gro File //
    std::ofstream gro_file;
    std::string gro_format = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f";
    gro_file.open("last_configuration" + state_name + ".gro");
    gro_file << "Writing this to a file.\n";
    gro_file << n_atoms <<"\n";
    const auto& active_particle_indices = system.getAllActiveAtomIndices();
    for( auto i : active_particle_indices){
        const auto particle = system.world.getParticle(i);
        const auto atom_type_name = system.atom_type_to_name.at(particle.atom_type);
        auto gro_line = string_format(gro_format, 1, "MOL", atom_type_name.c_str(), i+1,
                                      system.world.getParticle(i).r[XX],
                                      system.world.getParticle(i).r[YY],
                                      system.world.getParticle(i).r[ZZ]);
        gro_file << gro_line << "\n";
    }
    gro_file << box_size[XX] << " " << box_size[YY] << " " << box_size[ZZ];
    gro_file.close();
    ////////////////////
}

int main(int argc, char** argv) {

    std::string input_file_name = argv[2];
    std::ifstream input_file(input_file_name);
    json input;

    if(input_file.is_open()) {
        input = json::parse(input_file);
    } else{
        throw std::runtime_error("json input file: " + input_file_name + " does not exist or is inaccessible");
    }

    size_t n_atoms =     input.at("natoms");
    gmx::RVec box_size = {input.at("boxsizex"), input.at("boxsizey"), input.at("boxsizez")};
    real T = input.at("temp");
    real cutoff = input.at("cutoff");
    real dx0 = cutoff/5.;
    std::string state_name = input.at("statename");

    MoveData trans_data({0, 0, input.at("ntrans")});
    MoveData exchange_data({0, 0, input.at("nexc")});
    MoveData volume_data({0, 0, input.at("nvolume")});

    MovesSamplingData move_sampling({trans_data, exchange_data, volume_data});
    const long int n_steps = input.at("nsteps") ;
    const long int n_equil = input.at("nequil") ;
    const long int sampling_period = input.at("nsample") ;
    std::cout << "natoms = " << n_atoms << "\n";

    real P = 0.0 ;
    if(move_sampling.volume.n_attempts > 0){
        P = static_cast<real>(input.at("pressure")) * 0.00006022 ; // in code pressure is in reciprocal kj/mol/A^3 !!
    }

    test_simulation(n_atoms, move_sampling,
                    cutoff, box_size, T, P,
                    n_steps, n_equil, sampling_period,
                    state_name, dx0, input);


    return 1;

    std::shared_ptr<System> system;

    if(input.count("system") != 0) {
        system = std::make_shared<System>(input.at("system"));
    } else{
        throw std::runtime_error("input error! system keyword not found!");
    }

/*
    std::string input_file_name = argv[2];
    std::ifstream input_file(input_file_name);
    json input;
*/

    if(input_file.is_open()) {
        input = json::parse(input_file);
    } else{
        throw std::runtime_error("json input file: " + input_file_name + " does not exist or is inaccessible");
    }

//    std::shared_ptr<System> system;

    if(input.count("system") != 0) {
            system = std::make_shared<System>(input.at("system"));
    } else{
        throw std::runtime_error("input error! system keyword not found!");
    }

    std::vector<std::unique_ptr<TrajectoryAnalysis>> analyses;

    if(input.count("analysis") != 0) {
        if(!input.at("analysis").is_array()){
            throw std::runtime_error("json array expected");
        }
        for (const auto& j: input.at("analysis")) {
           const auto& [name, json_parameters] = jsonSingleItem(j);
           analyses.push_back(createAnalysis(name, json_parameters, system));
        }
    } else{
        throw std::runtime_error("input error! analysis keyword not found!");
    }

    for(auto& analysis : analyses){
        analysis->analyseAll();
    }

    return 0;
}


//! Jet fuel can't melt steel beams


