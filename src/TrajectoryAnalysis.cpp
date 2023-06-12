#include "core.h"
#include "Readers.h"
#include "TrajectoryAnalysis.h"
#include <memory>
#include <chrono>
#include "System.h"
#include "Systemx.h"
#include "Grid.h"
#include "Energy.h"


TrajectoryAnalysis::TrajectoryAnalysis(const std::string& trajectory_name, const Systemx* systemx) :
systemx(systemx)
{
    if(!trajectory_name.empty()){
        reader = createReader(trajectory_name, systemx->getAllActiveAtomIndices().size());
    }
}

std::vector<real> AnalysisRdf::calculateDistances(const size_t reference_particle_index,
                                                            const real* __restrict__ x, const real* __restrict__ y,
                                                            const real* __restrict__ z,
                                                            const real box_length, const size_t n_selected_atoms,
                                                            const size_t* __restrict__ selection_inds) const{
    const auto x_central = x[reference_particle_index];
    const auto y_central = y[reference_particle_index];
    const auto z_central = z[reference_particle_index];

    const auto half_box_length = box_length / 2.0f;

    std::vector<real> distances(n_selected_atoms);

    for (size_t i = 0; i<n_selected_atoms; ++i) {
        const auto dx = distancePBC(x_central, x[selection_inds[i]], half_box_length, box_length);
        const auto dy = distancePBC(y_central, y[selection_inds[i]], half_box_length, box_length);
        const auto dz = distancePBC(z_central, z[selection_inds[i]], half_box_length, box_length);
        distances[i] = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    return distances;
}


void NumberFluctsAnalysisOld::analyse(const TrajectoryFrame& frame){

}

void NumberFluctsAnalysisOld::analyseAll(){
    TrajectoryFrame frame(system->n_atoms);

    while(reader->read(frame)){
        if(frame.timestamp >= first_step) {
            n_samples++;
            updateNParticlesInSphere(frame, frame.n_atoms);
        }
    }
    writeOutput();
}



void NumberFluctsAnalysisOld::updateNParticlesInSphere(const TrajectoryFrame& frame, const size_t n_atoms){

    const auto& selection_inds = system->indicesOfAtom(selection1);

    for( auto j = 0; j<test_box_lengths.size(); ++j) {
        size_t n_particles1_in_region = 0;
        for (const auto ind: selection_inds) {

            std::vector<real> particle_coords = {frame.trajectory.x[ind], frame.trajectory.y[ind],
                                                   frame.trajectory.z[ind]};
            std::vector<real> ion_coords = {frame.trajectory.x[0], frame.trajectory.y[0], frame.trajectory.z[0]};

            n_particles1_in_region += is_in_region_sphere(particle_coords, ion_coords,
                                                          frame.bounding_box[0], probe_dr);

        }
        avg_n1_particles_in_boxes[j] = (avg_n1_particles_in_boxes[j] * (n_samples - 1)
                                          + n_particles1_in_region) / n_samples;
        avg_n1sq_particles_in_boxes[j] = (avg_n1sq_particles_in_boxes[j] * (n_samples - 1)
                                          + n_particles1_in_region * n_particles1_in_region) / n_samples;
    }
}

//void NumberFluctsAnalysisOld::updateNParticlesInSphere(const TrajectoryFrame& frame, const size_t n_atoms){
//
//    const auto& selection_inds = system->indicesOfAtom(selection1);
//
//    for( auto j = 0; j<test_box_lengths.size(); ++j) {
//        size_t n_particles1_in_region = 0;
//        for (const auto ind: selection_inds) {
//
//            std::vector<real> particle_coords = {frame.trajectory.x[ind], frame.trajectory.y[ind],
//                                                   frame.trajectory.z[ind]};
//            std::vector<real> ion_coords = {frame.trajectory.x[0], frame.trajectory.y[0], frame.trajectory.z[0]};
//
//            n_particles1_in_region += is_in_region_sphere(particle_coords, ion_coords,
//                                                          frame.bounding_box[0], probe_dr);
//
//        }
//        avg_n1_particles_in_boxes[j] = (avg_n1_particles_in_boxes[j] * (n_samples - 1)
//                                        + n_particles1_in_region) / n_samples;
//        avg_n1sq_particles_in_boxes[j] = (avg_n1sq_particles_in_boxes[j] * (n_samples - 1)
//                                          + n_particles1_in_region * n_particles1_in_region) / n_samples;
//    }
//}

void NumberFluctsAnalysisOld::updateNParticlesInBoxes(const TrajectoryFrame& frame, const size_t n_atoms){

    const auto& selection_inds1 = system->indicesOfAtom(selection1);
    const auto& selection_inds2 = system->indicesOfAtom(selection2);

    for( auto j = 0; j<test_box_lengths.size(); ++j){
        assert(test_box_llcs[j][0] + test_box_lengths[j] < frame.bounding_box[0]);
        size_t n_particles1_in_region = 0;
        size_t n_particles2_in_region = 0;
        for(const auto ind : selection_inds1){
            n_particles1_in_region += is_in_region({frame.trajectory.x[ind], frame.trajectory.y[ind], frame.trajectory.z[ind]},
                                                  test_box_lengths[j], test_box_llcs[j]); }
        for(const auto ind : selection_inds2){
            n_particles2_in_region += is_in_region({frame.trajectory.x[ind], frame.trajectory.y[ind], frame.trajectory.z[ind]},
                                                  test_box_lengths[j], test_box_llcs[j]); }

        avg_n1_particles_in_boxes[j] = (avg_n1_particles_in_boxes[j] * (n_samples-1)
                                       + n_particles1_in_region)/n_samples;
        avg_n2_particles_in_boxes[j] = (avg_n2_particles_in_boxes[j] * (n_samples-1)
                                        + n_particles2_in_region)/n_samples;
        avg_n1n2_particles_in_boxes[j] = (avg_n1n2_particles_in_boxes[j] * (n_samples-1)
                                       + n_particles1_in_region*n_particles2_in_region)/n_samples;
        avg_n1sq_particles_in_boxes[j] = (avg_n1sq_particles_in_boxes[j] * (n_samples-1)
                                          + n_particles1_in_region*n_particles1_in_region)/n_samples;
        avg_n2sq_particles_in_boxes[j] = (avg_n2sq_particles_in_boxes[j] * (n_samples-1)
                                          + n_particles2_in_region*n_particles2_in_region)/n_samples;
    }
}

NumberFluctsAnalysisOld::NumberFluctsAnalysisOld(const size_t dt, size_t first_step, const size_t n_boxes,
                                                 const std::vector<real> probe_dr,
                                                 const real min_box_length, const real max_box_length,
                                                 const std::string& traj_name, std::string output_name,
                                                 const std::string& selection1, const std::string& selection2,
                                                 const real sphere_radius, std::shared_ptr<System> system)
        :  avg_n1_particles_in_boxes(n_boxes), avg_n2_particles_in_boxes(n_boxes), avg_n1sq_particles_in_boxes(n_boxes),
           avg_n2sq_particles_in_boxes(n_boxes), avg_n1n2_particles_in_boxes(n_boxes),
           test_box_lengths(n_boxes), test_box_llcs(n_boxes, {0., 0., 0.}),
           sphere_radius(sphere_radius),
           probe_dr(probe_dr), first_step(first_step),
           selection1(selection1), selection2(selection2), output_filename(output_name) {

    traj_data.dt = dt;
    this->system = system;

    for( auto i = 0; i<n_boxes; ++i){
        test_box_lengths[i] = min_box_length + ((max_box_length - min_box_length)*i)/(n_boxes);
    }

    reader = createReader(traj_name, system->n_atoms);
}

NumberFluctsAnalysisOld::NumberFluctsAnalysisOld(const nlohmann::json& input,
                                                 std::shared_ptr<System> system)
        :
        NumberFluctsAnalysisOld(input.value("dt", 1),
                                input.value("first_step", 0),
                                input.value("n_boxes", 1),
                                {input.value("probe_dx", 0.0f), input.value("probe_dy", 0.0f), input.value("probe_dz", 0.0f)},
                                input.value("min_box_length", 1.0),
                                input.value("max_box_length", 1.0),
                                input.value("trajectory", "traj.crdbox"),
                                input.value("output", "fluct.out"),
                                input.value("sel1", ""),
                                input.value("sel2", ""),
                                input.value("probe_radius", 0.3),
                                std::move(system)) {}


NumberFluctsAnalysis::NumberFluctsAnalysis(const std::string& traj_name, std::string out_name,
                                           const std::string& selection, const Systemx& systemx)
        :  output_filename(out_name), selection(selection), TrajectoryAnalysis(traj_name, &systemx)
{}

NumberFluctsAnalysis::NumberFluctsAnalysis(const nlohmann::json& input, const Systemx& systemx)
    : NumberFluctsAnalysis(input.value("trajfile", ""), input.value("output", "n-flucts.dat"),
                           input.value("selection", ""), systemx)
{
    if(input.contains("volume")){
        const auto& volume_input = input.at("volume");
        const auto shape_name = volume_input.value("shape", "");
        if( shape_name == "sphere"){

            gmx::RVec r_center = {volume_input.at("center")[0], volume_input.at("center")[1], volume_input.at("center")[2]};
            real radius = volume_input.value("radius", 0);
            probe_volume = std::make_unique<Sphere>(r_center, radius);

        } else if( shape_name == "box"){

            gmx::RVec r_center = {volume_input.at("center")[0], volume_input.at("center")[1], volume_input.at("center")[2]};
            gmx::RVec l_sides = {volume_input.at("sides")[0], volume_input.at("sides")[1], volume_input.at("sides")[2]};
            probe_volume = std::make_unique<Box>(r_center, l_sides);
        }
    }
    else {
        probe_volume = std::make_unique<Volume>();
    }
}

void NumberFluctsAnalysis::analyse(const TrajectoryFrame &frame)
{
    const auto selection_type_index = systemx->name_to_atom_type.at(selection);
    const auto& selected_particle_indices =systemx->getActiveMolecules(selection_type_index);
    const auto n_selected_particles = selected_particle_indices.size();

    int n = 0;
    for(const auto molecule : selected_particle_indices){
        for(size_t world_index = molecule.first_; world_index<=molecule.last_; world_index++){
            n += probe_volume->isIn(systemx->world.getParticle(world_index).r);
        }
    }

    const auto box_size = systemx->grid->box_size_;
    const auto volume = box_size[XX]*box_size[YY]*box_size[ZZ];
    n_in_volume.add(n);
    volume_size.add(volume);
}

void NumberFluctsAnalysis::analyseAll()
{
}

void NumberFluctsAnalysis::writeOutput(){
    std::ofstream out_file(output_filename);

    const auto avg_V = volume_size.calcRunningAverage();
    const auto& n_samples = n_in_volume.getData();
    const auto& n_run_avg = n_in_volume.getData();

    const auto densities = n_in_volume / volume_size;
    const auto densities_avg = densities.calcRunningAverage();

    const auto& volume_samples = volume_size.getData();
    const auto& density_samples = densities.getData();

    for( int i = 0; i < volume_samples.size(); ++i ){
        out_file << i << "\t" <<  n_samples.at(i)       <<
                         "\t" <<  volume_samples.at(i)  <<
                         "\t" <<  density_samples.at(i) <<
                         "\t" <<  n_run_avg.at(i)       <<
                         "\t" <<  avg_V.at(i)           <<
                         "\t" <<  densities_avg.at(i)   <<  std::endl;
    }
    out_file.close();
}


AnalysisRdf::AnalysisRdf(const size_t dt, const real dr, const real rmax,
                         const std::string& traj_name, std::string output_name,
                         const std::string& reference, const std::string& selection, std::shared_ptr<System> system)
                            : dr(dr), hist(std::ceil(rmax/dr), 0), output_filename(output_name), reference(reference), selection(selection)
                            {

    traj_data.dt = dt;
    this->system = system;
    rdf.resize(std::ceil(rmax/dr));
    reader = createReader(traj_name, system->n_atoms);
}


AnalysisRdf::AnalysisRdf(const size_t dt, const real dr, const real rmax,
                         const std::string& traj_name, std::string output_name,
                         const std::string& reference, const std::string& selection, const Systemx& system)
        :   TrajectoryAnalysis(traj_name, &system),
            dr(dr), hist(std::ceil(rmax/dr), 0), output_filename(output_name),
            reference(reference), selection(selection) {

    traj_data.dt = dt;
    rdf.resize(std::ceil(rmax/dr));
}

AnalysisRdf::AnalysisRdf(const real dr, const real rmax,
                         std::string output_name,
                         std::string  reference, std::string  selection)
                         : dr(dr), hist(std::ceil(rmax/dr), 0), output_filename(std::move(output_name)),
                            reference(std::move(reference)), selection(std::move(selection)) {

    traj_data.dt = 1;
    this->system = system;
    rdf.resize(std::ceil(rmax/dr));
}

AnalysisRdf::AnalysisRdf(const nlohmann::json& input,
                         std::shared_ptr<System> system)
                         :
                         AnalysisRdf(input.value("dt", 1), input.value("dr", 0.02),
                                    input.value("rmax", 30),
                                input.value("trajectory", "traj.crdbox"),
                                input.value("output", "rdf.out"),
                                input.value("ref", ""),
                                input.value("sel", ""), system) {}


AnalysisRdf::AnalysisRdf(const nlohmann::json& input,
                         const Systemx& systemx)
        :
        AnalysisRdf(input.value("dt", 1), input.value("dr", 0.02),
                    input.value("rmax", 30),
                    input.value("trajectory", ""),
                    input.value("output", "rdf.dat"),
                    input.value("ref", ""),
                    input.value("sel", ""), systemx) {}

void AnalysisRdf::updateRdf(std::vector<real>& histogram,
                                 const real density, const real half_box_length, const size_t n_atoms){
    auto n_bins = histogram.size();
    real vol_sphere_prev = 0.0;
    for( int i = 0 ; i < n_bins; ++i){
        const auto r = (i+0.5)*dr;
        auto vol_sphere = sphereInBox(half_box_length, r);
        auto bin_volume = vol_sphere - vol_sphere_prev;
        vol_sphere_prev = vol_sphere;
        rdf[i] = (rdf[i]*(real)(n_samples-1) + histogram[i]/bin_volume/density/(real)n_atoms)/(real)n_samples;
    }
}

void AnalysisRdf::updateHistogramGRO(std::vector<u_int64_t>& histogram,  const std::vector<real>& values,
                                      const size_t n_bins) const{
    const auto n_distances = values.size();

    for( size_t i = 0; i < n_distances; ++i){
        const size_t bin_index = (size_t)std::floor((values[i] - dr/2.0)/dr) + 1;
        if(bin_index < n_bins ) {
            histogram[bin_index] += 1;
        }
    }
}
void AnalysisRdf::updateHistogramGRO(std::vector<real>& histogram,  const std::vector<real>& values,
                                     const size_t n_bins) const{
    const auto n_distances = values.size();

    for( size_t i = 0; i < n_distances; ++i){
        const size_t bin_index = (size_t)std::floor((values[i] - dr/2.0)/dr) + 1;
        if(bin_index < n_bins ) {
            histogram[bin_index] += 1;
        }
    }
}

void AnalysisRdf::writeOutput() {

    std::ofstream rdf_file(output_filename);
    real r = 0.0;
    for( const auto value : rdf ){
        rdf_file << r << "\t" <<  value << std::endl;
        r += dr;
    }
    rdf_file.close();
}

void AnalysisRdf::writeOutput(std::string output_name) {

    std::ofstream rdf_file(output_name);
    real r = 0.0;
    for( const auto value : rdf ){
        rdf_file << r << "\t" <<  value << std::endl;
        r += dr;
    }
    rdf_file.close();
}

std::vector<real> AnalysisRdf::integrateRDF(const std::vector<real>& gr, const real half_box_length) const{
    real prev_sphere_volume = 0;
    auto Gr = gr;
    for( size_t i = 0; i<gr.size(); ++i){
        real r = i * dr +  0.5*dr;
        auto sphere_volume =4./3.*r*r*r*M_PI; //sphereInBox(half_box_length, r);
        auto bin_volume = sphere_volume - prev_sphere_volume;
        Gr[i] = Gr[i-1] + (gr[i] - 1.0)*bin_volume;
        prev_sphere_volume = sphere_volume;
    }
    return Gr;
}

std::vector<real> AnalysisRdf::correctgrVdV(real V, real avg_density) {

    auto Gr = integrateRDF(rdf, 0.);
    auto result = Gr;

    int is_same_species = (selection == reference);

    for( size_t i = 0; i< rdf.size(); ++i){
        real r = i * dr +  0.5*dr;
        auto f = 1. - 4./3.*M_PI*r*r*r / V ;
        rdf[i] = rdf[i] * f / (f - Gr[i] / V - is_same_species * 1.0/(V*avg_density));
    }
    return result;
}

std::vector<real> AnalysisRdf::correctgrKruger(real half_box_length) const{

    std::vector<real> Gr(rdf.size(), 0.0);
    for( size_t j = 1; j<rdf.size(); ++j) {
        real R = j * dr + 0.5 * dr;
        for (size_t i = 1; i < j; ++i) {
            real r = i * dr + 0.5 * dr;
            real w = 4. * M_PI * r*r * (1 -  r*r*r/(R*R*R));
            Gr[j] +=  (rdf[i] - 1.0) * w * dr;
        }
    }
    return Gr;
}

void AnalysisRdf::analyseAll() {

    TrajectoryFrame frame(system->n_atoms);

    auto start = std::chrono::high_resolution_clock::now();

    while(reader->read(frame)){
        if(frame.number % traj_data.dt == 0 ){
            analyse(frame);
        }
        auto volume = std::pow(frame.bounding_box[0], 3);
    }

    writeOutput(output_filename);
    const auto& selection_inds = system->indicesOfAtom(selection);
    const auto dens_ref = (real)selection_inds.size() / (real)system->n_atoms * avg_density;
    rdf = correctgrVdV((real)system->n_atoms/avg_density, dens_ref);
//    rdf = correctgrVdV((real)system->n_atoms/avg_density, dens_ref, selection_inds.size());

    auto calc_time = std::chrono::high_resolution_clock::now() - start;
    std::cout << std::chrono::duration_cast<std::chrono::seconds>(calc_time).count() << std::endl;

   writeOutput(output_filename + "-VdVcorr");
}

void AnalysisRdf::analyse(const TrajectoryFrame &frame) {
    n_samples++;

    reference_inds = systemx->getActiveAtomsOf(systemx->name_to_atom_type.at(reference));
    selection_inds = systemx->getActiveAtomsOf(systemx->name_to_atom_type.at(selection));

    auto n_bins = rdf.size();
    std::vector<real> histogram(n_bins, 0.);
    auto box_length = frame.bounding_box[0];
    auto volume = box_length*box_length*box_length;
    avg_density = (avg_density * ((real)n_samples-1) + (real)frame.n_atoms/volume)/((real)n_samples);

    auto& reference_inds_used = reference_inds;
    auto& selection_inds_used = selection_inds;
    if(system.use_count() != 0){
        reference_inds_used = system->indicesOfAtom(reference);
        selection_inds_used = system->indicesOfAtom(selection);
    }

    for( const auto i  : reference_inds_used){
        auto distances = calculateDistances(i, frame.trajectory.x.data(), frame.trajectory.y.data(),
                                             frame.trajectory.z.data(), box_length,
                                             selection_inds_used.size(), selection_inds_used.data());
        updateHistogramGRO(histogram, distances, n_bins);
        if(reference == selection){
            histogram[0] = 0;}
    }
    real number_of_distances = 0;
    size_t n_atoms = reference_inds_used.size();
    for(auto x : histogram){ number_of_distances += x;}
    std::cout<< "rdf histogram has =  " << static_cast<size_t>(std::floor(number_of_distances)) << " total distances ... there should be: " <<
        n_atoms*n_atoms - n_atoms<< " distances in total\n";

    real density = real(n_atoms) / std::pow(box_length, 3);
    updateRdf(histogram, density, box_length/2.0f, selection_inds_used.size());
}

AnalysisRdf::AnalysisRdf(AnalysisRdf& rdf) {

}

void AnalysisWidom::analyseAll() {

    TrajectoryFrame frame(systemx->world.getAllActiveParticleIndices().size());

    auto start = std::chrono::high_resolution_clock::now();

    while (reader->read(frame)) {
        if (frame.number % traj_data.dt == 0) {
            n_samples++;
            analyse(frame);
        }
    }

    writeOutput(output_filename);
}

AnalysisWidom::AnalysisWidom(const nlohmann::json& input, const Systemx& systemx, Energy* energizer)
    :
        TrajectoryAnalysis(input.value("trajfile", ""), &systemx),
        n_insertions(input.value("n_insertions", 1000)),
        energizer_(energizer),
        output_filename(input.value("output", "widom.dat"))
{
    inserted_molecule_type = this->systemx->name_to_molecule_type.at(input.value("molecule", ""));
    inserted_molecule_coordinates = this->systemx->molecule_type_to_structure.at(inserted_molecule_type);
}


void AnalysisWidom::analyse(const TrajectoryFrame &frame) {

    std::random_device r;
    std::mt19937 e1(r());
    std::uniform_real_distribution<real> uniform_dist(0,1);

    real cutoff = systemx->cutoff;
    real beta = systemx->beta;

    real widom_energy = 0.0;
    real energy_correction = 0.0;
    gmx::RVec box_size = systemx->grid->box_size_;
    real volume = box_size[XX]*box_size[YY]*box_size[ZZ];
    size_t n_atom_types = systemx->n_atom_types;

    //! correction
    for(auto inserted_particle : inserted_molecule_coordinates) {
        const auto inserted_atom_type = inserted_particle.atom_type;
        for (int atom_type = 0; atom_type < n_atom_types; ++atom_type) {

            const auto sig = systemx->lj_table.get_sig(inserted_atom_type, atom_type);
            const auto eps = systemx->lj_table.get_eps(inserted_atom_type, atom_type);

            const real n_atoms_j = systemx->getActiveAtomsOf(atom_type).size();
            const real rhoj_red = (n_atoms_j / volume) * std::pow(sig, 3);


            energy_correction += eps * 8. / 3. * M_PI * rhoj_red *
                                 (1. / 3. * std::pow(sig / cutoff, 9) -
                                  std::pow(sig / cutoff, 3));
        }
    }

    //! insertion
    for (int i_insert = 0; i_insert < n_insertions; ++i_insert) {


        gmx::RVec r_insert(uniform_dist(e1) * box_size[XX],
                           uniform_dist(e1) * box_size[YY],
                           uniform_dist(e1) * box_size[ZZ]);


        generateRandomPosition(inserted_molecule_coordinates, box_size, uniform_dist, e1);

        const auto next_widom_energy =
                energizer_->calcParticlesEnergy(*systemx,
                                                  inserted_molecule_coordinates,
                                                  box_size, cutoff * cutoff);

        widom_energy = (widom_energy*i_insert + std::exp(-beta*next_widom_energy))/(i_insert+1);
    }
    exp_mbU.add(widom_energy * std::exp(-beta*energy_correction));
    V.add(volume);
}


void AnalysisWidom::writeOutput(){
    writeOutput(output_filename);
}

void AnalysisWidom::writeOutput(std::string out_file_name){
    std::ofstream out_file(out_file_name);
    real beta = systemx->beta;
    const auto& exp_mbUs = exp_mbU.getData();
    const auto Vexp_mbU = exp_mbU * V;
    const auto avg_Vexp_mBus = Vexp_mbU.calcRunningAverage();
    const auto avg_V = V.calcRunningAverage();
    for( int i = 0; i < exp_mbUs.size(); ++i ){
        out_file << i << "\t" <<  -1./beta * std::log(avg_Vexp_mBus.at(i) / avg_V.at(i)) <<
                 "\t" <<  exp_mbUs.at(i)    <<
                 "\t" <<  V.getData().at(i) <<   std::endl;
    }
    out_file.close();
}

real AnalysisWidom::widomEnergy() const{
    real beta = systemx->beta;
    const auto Vexp_mbU = exp_mbU * V;
    const auto avg_Vexp_mbU = Vexp_mbU.calcAverage();
    const auto avg_V = V.calcAverage();

    real n_i = systemx->getActiveMolecules(inserted_molecule_type).size();
    real rho_i = n_i/avg_V;
    const auto widom_ideal = 1./beta*std::log(rho_i);

    return -1./beta * std::log(avg_Vexp_mbU / avg_V) + 0*widom_ideal;
}






//// BACKUP

//void AnalysisGyration::analyseAll() {
//
 //    TrajectoryFrame frame(system->n_atoms);
//
//    auto start = std::chrono::high_resolution_clock::now();
//
//    while(reader->read(frame)){
//        if(frame.number % traj_data.dt == 0 ){
//            n_samples++;
//            analyse(frame);
//        }
//        auto volume = std::pow(frame.bounding_box[0], 3);
//        avg_density = avg_density * (n_samples-1)/n_samples + frame.n_atoms/volume/n_samples;
//    }
//
//    writeOutput(output_filename);
//    const auto& selection_inds = system->indicesOfAtom(selection);
//    const auto dens_ref = (real)selection_inds.size() / (real)system->n_atoms * avg_density;
//
//    auto calc_time = std::chrono::high_resolution_clock::now() - start;
//    std::cout << std::chrono::duration_cast<std::chrono::seconds>(calc_time).count() << std::endl;
//
//    writeOutput(output_filename + "-asdf");
//}
//
//void AnalysisGyration::analyse(const TrajectoryFrame &frame) {
//
//
//    auto box_length = frame.bounding_box[0];
//    const auto& selection_inds = system->indicesOfAtom(selection);
//
//    const auto Rg = gyrationRadius();
//    real density = real(selection_inds.size()) / std::pow(box_length, 3);
//}
//
//real AnalysisGyration::gyrationRadius(std::vector<real> xs,
//                                        std::vector<real> ys,
//                                        std::vector<real> zs) {
//
//
//    real Rg = 0;
//
//    for( int i = 0; i < xs.size(); i++ ) {
//
//    }
//        for( int i = 0; i < xs.size(); i++ ){
//        Rg += 1;
//    }
//
//}

void Pressure::analyseAll() {

    TrajectoryFrame frame(systemx->world.getAllActiveParticleIndices().size());

    auto start = std::chrono::high_resolution_clock::now();

    while (reader->read(frame)) {
        if (frame.number % traj_data.dt == 0) {
            n_samples++;
            analyse(frame);
        }
    }

    writeOutput(output_filename_);
}

Pressure::Pressure(const nlohmann::json& input, const Systemx& systemx, Energy* energizer)
        :
        TrajectoryAnalysis(input.value("trajfile", ""), &systemx),
        sampling_period_(input.value("samplingperiod", 1000)),
        energizer_(energizer),
        output_filename_(input.value("output", "pressure.dat"))
{

}


void Pressure::analyse(const TrajectoryFrame &frame) {

/*    if(frame.step%sampling_period_!=0){
        return;
    }*/

    const real n_atoms = systemx->getAllActiveAtomIndices().size();
    const auto cutoff = systemx->cutoff;
    const gmx::RVec box_size = {frame.bounding_box[XX], frame.bounding_box[YY], frame.bounding_box[ZZ]};
    const auto volume = frame.bounding_box[XX]*frame.bounding_box[YY]*frame.bounding_box[ZZ];
    const auto beta = systemx->beta;
    const auto virial = energizer_->calcVirial(*systemx, box_size, cutoff*cutoff);
    const auto pressure_correction = pressureCorrectionScaledCutoff(*systemx, box_size, cutoff);
    const auto pressure = (virial /volume + ((real)n_atoms) / (beta*volume))*1.0e30/6.022e23*1e3*1e-5 + pressure_correction;

    std::cout << "correction = " << pressure_correction << "\n";
    std::cout << "virial = " << virial /volume * 1.0e30/6.022e23*1e3*1e-5 <<"\n";
    pressures_data_.add(pressure);
}

real Pressure::pressureCorrection(const Systemx& system,
                                  const gmx::RVec box_size, const real cutoff)const{

    const auto volume = box_size[XX] * box_size[YY] * box_size[ZZ];

    const size_t n_atom_types = system.atom_type_to_lj_params.size();
    real pressure_correction = 0.0;
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


real Pressure::pressureCorrectionScaledCutoff(const Systemx &system, const gmx::RVec box_size,
                                              const real cutoff) const {
    const auto volume = box_size[XX] * box_size[YY] * box_size[ZZ];
    return energyCorrection(system, *energizer_, volume, cutoff)/volume * (1e30*1e3/6.022e23*1e-5);
}

void Pressure::writeOutput(){
    writeOutput(output_filename_);
}

void Pressure::writeOutput(std::string out_file_name) const{
    std::ofstream out_file(out_file_name);
    const auto& pressures = pressures_data_.getData();
    const auto avg_pressures = pressures_data_.calcRunningAverage();
    for( int i = 0; i < pressures.size(); ++i ){
        out_file << i << "\t" <<  pressures.at(i) <<  "\t" << avg_pressures.at(i) << "\n";
    }
    out_file.close();

    std::ofstream acf_file(output_filename_ + "_acf.dat");
    const auto acf = pressures_data_.calcAutoCorrelationFunction(pressures.size()/3);
    for( int i = 0; i < acf.size(); ++i ){
        acf_file << i << "\t" <<  acf.at(i) <<"\n";
    }
    acf_file.close();
}

