#pragma once
#include <vector>
#include <iostream>
#include <memory>
#include <functional>
#include <tuple>
#include <nlohmann/json.hpp>


class Reader;
class Coordinates;
class Box;
class TrajectoryFrame;

struct TrajAnalysisData{
    size_t first_frame = 0;
    size_t dt = 1;
    std::string  traj_file_name;
    std::string output_filename;
};


template <class Type>
class TimeSeries{

    std::vector<Type> data_;

public:
    void add(Type datum){
        data_.push_back(datum);
    }

    [[nodiscard]] const std::vector<Type>& getData() const{
        return data_;
    }

    [[nodiscard]] Type calcAverage() const{
        Type average_data = 0.;
        int n = 0;
        for(const auto& datum : data_){
            average_data = (average_data*n + datum)/(n+1);
            n++;
        }
        return average_data;
    }

    [[nodiscard]] std::vector<Type> calcRunningAverage(size_t n_begin = 0) const{
        int n_samples = data_.size();
        std::vector<Type> average_data(n_samples - n_begin, 0);
        for(int n = 0; n<n_samples-n_begin - 1; ++n){

            auto sample_index = n + n_begin;
            average_data.at(n+1) = (average_data.at(n)*n + data_.at(sample_index))/(n+1);
        }
        return average_data;
    }

    Type calcSig2() const{
        Type avg_data = calcAverage();
        Type avg_datasq =0;// calcAverage();
        real n = 0;
        for(const auto& datum : data_){
            avg_datasq = (avg_datasq*n + datum*datum)/(n+1);
            n++;
        }
        return avg_datasq - avg_data*avg_data;
    }


    [[nodiscard]] std::vector<Type> calcAutoCorrelationFunction(const size_t n_max) const{
        const auto avg_datum = calcAverage();
        const auto sig2_datum = calcSig2();
        std::vector<Type> autocorrelation(n_max);
        assert(n_max < data_.size());
        for( int dt = 0; dt < n_max; ++dt){
            for(int t = 0; t<n_max-dt; ++t){
                autocorrelation[dt] += (data_.at(t) - avg_datum)
                                       * (data_.at(t+dt) - avg_datum)/((real)(n_max - dt)*sig2_datum);
            }
        }
        return autocorrelation;
    }

    TimeSeries<Type> operator *(const TimeSeries<Type>& ts1) const{
        TimeSeries<Type> ts_result;
        int i = 0;
        for(const auto& datum1 : ts1.getData()){
            ts_result.add(datum1 * this->data_.at(i));
            i++;
        }
        return ts_result;
    }

    TimeSeries<Type> operator /(const TimeSeries<Type>& ts1) const{
        TimeSeries<Type> ts_result;
        int i = 0;
        for(const auto& datum1 : ts1.getData()){
            ts_result.add(datum1 / this->data_.at(i));
            i++;
        }
        return ts_result;
    }
};


class System;
class Systemx;
class Energy;

class TrajectoryAnalysis{

protected:
    TrajAnalysisData traj_data;
    size_t n_samples = 0;
    std::unique_ptr<Reader> reader;
    std::shared_ptr<System> system;
    const Systemx* systemx;

public:
    TrajectoryAnalysis() = default;
    TrajectoryAnalysis(const std::string& trajectory_name, const Systemx* systemx);
    virtual ~TrajectoryAnalysis() = default;

    virtual void analyse(const TrajectoryFrame& frame) = 0;
    virtual void analyseAll() = 0;

    virtual void writeOutput() = 0 ;

};


class NumberFluctsAnalysisOld : public TrajectoryAnalysis{

    std::vector<real> avg_n1_particles_in_boxes;
    std::vector<real> avg_n2_particles_in_boxes;
    std::vector<real> avg_n1n2_particles_in_boxes;
    std::vector<real> avg_n1sq_particles_in_boxes;
    std::vector<real> avg_n2sq_particles_in_boxes;
    std::vector<real> test_box_lengths;
    std::vector<std::vector<real>> test_box_llcs; // coordinates of lower left corners
    std::vector<real> probe_dr;
    real sphere_radius;
    std::string selection1;
    std::string selection2;
    std::string output_filename;
    std::vector<size_t> n1_particles_in_sphere;
    size_t first_step;

public:
    NumberFluctsAnalysisOld() : TrajectoryAnalysis(){

    }

    NumberFluctsAnalysisOld(const size_t dt, size_t first_step, const size_t n_boxes, const std::vector<real> probe_dr,
                            const real min_box_length, const real max_box_length,
                            const std::string& traj_name, std::string output_name,
                            const std::string& selection1, const std::string& selection2,
                            const real sphere_radius, std::shared_ptr<System> system);
    NumberFluctsAnalysisOld(const nlohmann::json& input, std::shared_ptr<System> system);

    bool is_in_region(const std::vector<real>& point, const real test_box_length, std::vector<real> test_box_llc ){
        const auto is_in_x = (point[0] - test_box_llc[0])<test_box_length and (point[0] - test_box_llc[0])>0;
        const auto is_in_y = (point[1] - test_box_llc[1])<test_box_length and (point[1] - test_box_llc[1])>0;
        const auto is_in_z = (point[2] - test_box_llc[2])<test_box_length and (point[2] - test_box_llc[2])>0;
        return is_in_x and is_in_y and is_in_z;
    }

    bool is_in_region_sphere(const std::vector<real>& point, const std::vector<real>& ion_center,
                             const real box_length, const std::vector<real> sphere_distance ){

        const real half_box_length = box_length / 2.0;

        real probe_sphere_center_x = ion_center[0] + sphere_distance[0];
        real probe_sphere_center_y = ion_center[1] + sphere_distance[1];
        real probe_sphere_center_z = ion_center[2] + sphere_distance[2];
        if(probe_sphere_center_x > box_length){
            probe_sphere_center_x -= box_length;
        } else if(probe_sphere_center_x < 0.0){
            probe_sphere_center_x += box_length;
        }

        if(probe_sphere_center_y > box_length){
            probe_sphere_center_y -= box_length;
        } else if(probe_sphere_center_x < 0.0){
            probe_sphere_center_y += box_length;
        }

        if(probe_sphere_center_z > box_length){
            probe_sphere_center_z -= box_length;
        } else if(probe_sphere_center_z < 0.0){
            probe_sphere_center_z += box_length;
        }
        const auto dx = distancePBC(point[0], probe_sphere_center_x, half_box_length, box_length);
        const auto dy = distancePBC(point[1], probe_sphere_center_y, half_box_length, box_length);
        const auto dz = distancePBC(point[2], probe_sphere_center_z, half_box_length, box_length);

        return (dx*dx + dy*dy + dz*dz)<sphere_radius*sphere_radius;
    }

    void analyse(const TrajectoryFrame& frame) override;
    void analyseAll() override;

    void writeOutput() override{

        std::ofstream out_file(output_filename);
        for( auto j = 0; j<test_box_lengths.size(); ++j){
            const auto V = test_box_lengths[j]*test_box_lengths[j]*test_box_lengths[j];
            const auto n1 = avg_n1_particles_in_boxes[j];
            const auto n2 = avg_n2_particles_in_boxes[j];
            const auto n1sq = avg_n1sq_particles_in_boxes[j];
            const auto n2sq = avg_n2sq_particles_in_boxes[j];
            const auto n1n2 = avg_n1n2_particles_in_boxes[j];
//            out_file << test_box_lengths[j] << "\t" <<  n1 << "\t" << n2 << "\t";
//            out_file << V *(n1sq - n1*n1 - n1)/(n1*n1) << "\t" <<  V*(n2sq - n2*n2 - n2)/(n2*n2) << "\t" << V*(n1n2 - n1*n2)/(n1*n2) << "\n";

            out_file << test_box_lengths[j] << "\t" <<  n1 << "\t";
            out_file << (n1sq - n1*n1) << "\n";
        }

    }

        private:
        void updateNParticlesInBoxes(const TrajectoryFrame& frame, size_t n_atoms);
        void updateNParticlesInSphere(const TrajectoryFrame& frame, size_t n_atoms);

};

class Volume{

public:
    virtual bool isIn(gmx::RVec r) const{
        return true;
    };
};

class Box : public Volume{

    gmx::RVec lower_left;
    gmx::RVec upper_right;

public:
    Box(gmx::RVec r_center, gmx::RVec l_sides){
        lower_left  = r_center - l_sides/2.0;
        upper_right = r_center + l_sides/2.0;
    }

    bool isIn(gmx::RVec r) const override{
        return  r >= lower_left and r < upper_right ;
    }
};

class Sphere : public Volume{

    gmx::RVec r_center;
    real radius_sq;

public:
    Sphere(gmx::RVec r_center, real radius_sq) :
        r_center(r_center), radius_sq(radius_sq) {}

    bool isIn(gmx::RVec r) const override{
        return  (r_center - r).norm2() <= radius_sq;
    }
};

class NumberFluctsAnalysis : public TrajectoryAnalysis{


    std::unique_ptr<Volume> probe_volume;
    TimeSeries<real> n_in_volume;
    TimeSeries<real> volume_size;
    std::string selection;
    std::string output_filename;

public:

    NumberFluctsAnalysis(const std::string& traj_name, std::string out_name, const std::string& selection, const Systemx& systemx);
    NumberFluctsAnalysis(const nlohmann::json& input, const Systemx& system);

    void analyse(const TrajectoryFrame& frame) override;
    void analyseAll() override;

    void writeOutput() override;

private:
};



class AnalysisRdf : public TrajectoryAnalysis {

    AnalysisRdf(AnalysisRdf &rdf);

    real dr=0.08;
    std::vector<real> rdf;
    std::vector<uint64_t> hist;
    std::string output_filename;
    std::string reference;
    std::string selection;
public:

    real avg_density = 0;
    std::vector<size_t> reference_inds;
    std::vector<size_t> selection_inds;


public:
    AnalysisRdf(size_t dt, real dr, real rmax, const std::string& file_name, std::string out_name,
                const std::string& reference, const std::string& selection, std::shared_ptr<System> system);
    AnalysisRdf(size_t dt, real dr, real rmax, const std::string& file_name, std::string out_name,
                const std::string& reference, const std::string& selection, const Systemx& systemx);
    AnalysisRdf( real dr = 0.02, real rmax = 40, std::string out_name = "", std::string  reference = "", std::string  selection = "");
    AnalysisRdf(const nlohmann::json& input, std::shared_ptr<System> system);
    AnalysisRdf(const nlohmann::json& input, const Systemx& system);

    void analyse(const TrajectoryFrame& frame) override;
    void analyseAll() override;

    void writeOutput() override;
    void writeOutput(std::string);

    [[nodiscard]] std::vector<real> correctgrVdV(real half_box_length, real avg_density);
private:
    void updateRdf(std::vector<real>& histogram, real density,
                        real half_box_length,  size_t n_atoms);

    void histogramToRdf(std::vector<u_int64_t>& histogram, real density,
                        real half_box_length,  size_t n_atoms);
    void updateHistogramGRO( std::vector<real>& histogram,  const std::vector<real>& values,
                            size_t n_bins) const;

    void updateHistogramGRO( std::vector<u_int64_t>& histogram,  const std::vector<real>& values,
                             size_t n_bins) const;

    std::vector<real> integrateRDF(const std::vector<real>& rdf, real half_box_length) const;

    std::vector<real> correctgrKruger(real half_box_length) const;

    std::vector<real> calculateDistances(const size_t reference_particle_index,
                       const real* __restrict__ x, const real* __restrict__ y,
                       const real* __restrict__ z,
                       const real box_length, const size_t n_selected_atoms,
                       const size_t* __restrict__ selection_inds) const;

};


class Pressure : public TrajectoryAnalysis {

    size_t sampling_period_ = 0;
    TimeSeries<real> pressures_data_;
    std::string output_filename_;
    Energy* energizer_;

public:
    Pressure(const nlohmann::json& input, const Systemx& system, Energy* energizer);

    void analyse(const TrajectoryFrame& frame) override;
    void analyseAll() override;
    void writeOutput() override;
    void writeOutput(std::string out_file_name) const;

private:
    real pressureCorrection(const Systemx& system,
                            const gmx::RVec box_size, const real cutoff)const;
    real pressureCorrectionScaledCutoff(const Systemx& system,
                            const gmx::RVec box_size, const real cutoff)const;

};


class AnalysisWidom : public TrajectoryAnalysis {


    std::vector<uint64_t> hist;
    std::string output_filename;
    std::vector<Particle> inserted_molecule_coordinates;

    size_t n_insertions;
    const Energy* energizer_;

public:
    size_t inserted_molecule_type;

    TimeSeries<real> exp_mbU;
    TimeSeries<real> V;


public:

    AnalysisWidom(const nlohmann::json& input, const Systemx& system, Energy* energizer);

    void analyse(const TrajectoryFrame& frame) override;
    void analyseAll() override;

    void writeOutput() override;
    void writeOutput(std::string);
    real widomEnergy() const;



private:

};


//class AnalysisGyration : public TrajectoryAnalysis {
//
//    std::vector<real> Rg;
//    std::vector<uint64_t> hist;
//    std::string output_filename;
//    real avg_density = 0;
//    std::string selection;
//    std::vector<u_int> selection_inds;
//
//public:
//    AnalysisGyration(size_t dt, real dr, real rmax, const std::string& file_name, std::string out_name,
//                const std::string& reference, const std::string& selection, std::shared_ptr<System> system);
//    AnalysisGyration(const nlohmann::json& input, std::shared_ptr<System> system);
//
//    void analyse(const TrajectoryFrame& frame) override;
//    void analyseAll() override;
//
//    void writeOutput() override;
//    void writeOutput(std::string);
//
//private:
//
//    real gyrationRadius(const TrajectoryFrame& frame);
//
//};

