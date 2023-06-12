#include "core.h"
#include "Readers.h"
#include <sstream>


void CrdReader::initialize(){

    std::string line;
    file.open(file_name);

    if (file.is_open()) {

        timestep.n_atoms = total_atoms;
        timestep.trajectory.x.resize(total_atoms);
        timestep.trajectory.y.resize(total_atoms);
        timestep.trajectory.z.resize(total_atoms);
        std::getline(file, line); // this line is a comment with TITLE
        size_t offset = 1;

        std::string coordinate;
        size_t i = 0;
        offsets.push_back(offset);
        while (std::getline(file, line) ) {
            std::stringstream ss(line);

            while(!ss.eof()){
                i++;
                ss >> coordinate;
            }

            offset += 1;
            
            if(i == 3*total_atoms){ // we encountered a line with box lengths
                std::getline(file, line);
                offset += 1 ;
                n_atoms_in_frame.push_back(i/3);
                offsets.push_back(offset);
                n_frames += 1;
                i = 0;
                continue;
            } 
        }
        file.close();
        first_frame_read = false;
    }
}

bool CrdReader::read(TrajectoryFrame& frame){

    if(!first_frame_read){
        file.open(file_name);
        if(file.is_open()){
            updateBoxAndCoordinates(frame, offsets[frame.step]);
            timestep.step = frame.step;
            first_frame_read = true;
            return true;
        } else{
            throw std::runtime_error("file: " + file_name + " not found or inaccessible!");
        }
    }
    if(file.is_open()) {
        frame.step++;
        timestep.step = frame.step;
        return updateBoxAndCoordinates(frame, 0);
    }else{
        return false;
    }
}

//const TrajectoryFrame& CrdReader::next(){
//
//    if(file.is_open()){
//        updateBoxAndCoordinates(0);
//    }
//    timestep.step += 1;
//    return timestep;
//}

bool CrdReader::updateBoxAndCoordinates(TrajectoryFrame& frame, size_t offset) {

    std::string line;

    for( std::size_t i = 0; i < offset; ++i ){
        std::getline(file, line);
    }

    size_t i = 0;
    while (std::getline(file, line) ) {
        std::stringstream ss(line);
        while(!ss.eof()){
            auto coord_index = i%3;
            auto atom_index = i/3;
            switch (coord_index) {
                case 0:
                    ss >> frame.trajectory.x.at(atom_index);
                    break;
                case 1:
                    ss >> frame.trajectory.y.at(atom_index);
                    break;
                case 2:
                    ss >> frame.trajectory.z.at(atom_index);
                    break;
            }
            i++;
        }

        if(i == 3*frame.n_atoms){ // we encountered a line with box lengths
            std::getline(file, line);
            std::stringstream ss2(line);
            ss2 >> frame.bounding_box[0] ;
            ss2 >> frame.bounding_box[1] ;
            ss2 >> frame.bounding_box[2] ;
            return true;
        }
        else if( i > 3*frame.n_atoms){
            throw std::runtime_error("reader failed due to wrong number of atoms");
        }
    }
    file.close();
    first_frame_read = false;
    return false;
}
//
//bool CrdReader::isThereNext() const{
//    if (timestep.step + 1 == n_frames){
//        return false;
//    }
//    return true;
//}



bool XYZReader::read(TrajectoryFrame& frame){
/*
    file.open(file_name);
    if (file.is_open()) {
        first_frame_read = true;
        timestep.step = frame.step;
        updateBoxAndCoordinates(offsets.at(frame.step));
    }
    else{
        throw std::runtime_error("file: " + file_name + " not found or inaccessible");
    }
*/
    if (first_frame_read) {
        timestep.step += 1;
        updateBoxAndCoordinates(0);
        return true;
    } else {
        throw std::runtime_error(file_name + " could not be opened");
    }
}

bool XYZReader::updateBoxAndCoordinates(const size_t offset) {
    std::string line;
    for( std::size_t i = 0; i < offset; ++i ){
        std::getline(file, line);
    }

    std::getline(file, line); // n_atoms
    timestep.reset(std::stoi(line));

    std::getline(file, line); // time info

    for(std::size_t i = 0; i < n_atoms_in_frame[timestep.step]; ++i )
    {
        std::getline(file, line);
        const auto separated_line = separateLine(line, " ");
        timestep.trajectory.x.push_back(std::stof(separated_line[1]));
        timestep.trajectory.y.push_back(std::stof(separated_line[2]));
        timestep.trajectory.z.push_back(std::stof(separated_line[3]));
        if( i==n_atoms_in_frame[timestep.step]){
            return true;
        }
    }
    return false;
}

void XYZReader::initialize(){

    std::string line;
    file.open(file_name);
    if (file.is_open()) {

        size_t offset = 0;
        while (std::getline(file, line)) {
            if (line.find(' ') == std::string::npos) {
                n_frames += 1;
                n_atoms_in_frame.push_back(std::stoi(line));
                offsets.push_back(offset);
            }
            offset += 1;
        }
        file.close();
        first_frame_read = false;
    }
}

/*
const TrajectoryFrame& XYZReader::next(){
    if (first_frame_read) {
        timestep.step += 1;
        updateBoxAndCoordinates(0);
    } else {
        throw std::runtime_error(file_name + " could not be opened");
    }
    return timestep;
}
*/

/*
bool XYZReader::isThereNext() const{
    if (timestep.step + 1 == n_frames){
        return false;
    }
    return true;
}
*/

std::vector<std::string> XYZReader::separateLine(std::string& line, std::string&& delimiter) {
    size_t pos;
    size_t start_pos = line.find_first_not_of(' ', 0);
    std::vector<std::string> separated_line;

    while ((pos = line.find(' ', start_pos)) != std::string::npos) {
        separated_line.push_back(line.substr(start_pos, pos-start_pos));
        start_pos = line.find_first_not_of(' ', pos);
    }
    separated_line.push_back(line.substr(start_pos, line.size()));
    return separated_line;
}


// ========== XTCTrajectoryFrame ==========

XTCTrajectoryFrame::XTCTrajectoryFrame(int number_of_atoms) { initNumberOfAtoms(number_of_atoms); }

XTCTrajectoryFrame::XTCTrajectoryFrame(const TrajectoryFrame& frame) {
    initNumberOfAtoms(static_cast<int>(frame.trajectory.x.size()));
    importFrame(frame);
}

XTCTrajectoryFrame& XTCTrajectoryFrame::operator=(const TrajectoryFrame& frame) {
    if (frame.trajectory.x.size() != number_of_atoms) {
        throw std::runtime_error("wrong number of particles to be assign into the XTC step");
    }
    importFrame(frame);
    return *this;
}

void XTCTrajectoryFrame::importFrame(const TrajectoryFrame& frame) {
    importTimestamp(frame.step, frame.timestamp);
    importBox(frame.box);
    importCoordinates(frame.trajectory);
}

void XTCTrajectoryFrame::importTimestamp(const int step, const float time) {
    xtc_step = step;
    xtc_time = time / static_cast<float>(1.0);
}

void XTCTrajectoryFrame::importBox(const std::vector<real>& box) {
    // empty box tensor
    float xtc_box_matrix [DIM][DIM] ;
    // only XYZ dimensions in nanometers on diagonal, as floats
    xtc_box_matrix[0][0] = static_cast<float>(box[0]);
    xtc_box_matrix[1][1] = static_cast<float>(box[1]);
    xtc_box_matrix[2][2] = static_cast<float>(box[2]);

    xtc_box_matrix[0][1] = 0.0;  xtc_box_matrix[0][2] = 0.0;
    xtc_box_matrix[1][0] = 0.0;  xtc_box_matrix[1][2] = 0.0;
    xtc_box_matrix[2][0] = 0.0;  xtc_box_matrix[2][1] = 0.0;

    // copy underlaying eigen structure (1D array, row-major) to the C-style 2D array
    std::copy(&xtc_box_matrix[0][0], &(xtc_box_matrix[0][0])+DIM*DIM, &(xtc_box[0][0]));
}

void XTCTrajectoryFrame::importCoordinates(const Coordinates& coordinates) {
    // setNumberOfAtoms(coordinates.size());
    if (coordinates.x.size() != number_of_atoms) {
        // to avoid mistakes, the number_of_atoms is immutable
        throw std::runtime_error("wrong number of particles to be saved in the XTC step");
    }
    for (int i = 0; i < number_of_atoms; ++i) {
        // coordinates in nanometers, as floats
        const std::vector<XTCFloat> xtc_pos = {static_cast<XTCFloat>(coordinates.x[i]),
                                               static_cast<XTCFloat>(coordinates.y[i]),
                                               static_cast<XTCFloat>(coordinates.z[i])};        // copy underlaying eigen structure (1D array) to the correct place in C-style 2D array
        std::copy(xtc_pos.data(), xtc_pos.data() + DIM, xtc_coordinates.get()[i]);
    }
}

void XTCTrajectoryFrame::exportFrame(TrajectoryFrame& frame) const {
    exportTimestamp(frame.step, frame.timestamp);
    exportBox(frame.bounding_box);
    exportCoordinates(frame.trajectory);
}

void XTCTrajectoryFrame::exportTimestamp(size_t& step, float& time) const {
    step = xtc_step;
    time = xtc_time * static_cast<float>(1.0);
}

void XTCTrajectoryFrame::exportBox(real box[3] ) const {
    box[0] = static_cast<real>(xtc_box[0][0]) * 10.;
    box[1] = static_cast<real>(xtc_box[1][1]) * 10.;
    box[2] = static_cast<real>(xtc_box[2][2]) * 10.;
}

void XTCTrajectoryFrame::exportCoordinates(Coordinates& coordinates) const {
    if (coordinates.x.size() != number_of_atoms) {
        throw std::runtime_error("wrong number of particles in the loaded XTC step");
    }
    for (int i = 0; i < number_of_atoms; ++i) {
        coordinates.x[i] = static_cast<real>(xtc_coordinates.get()[i][0]) * 10.;
        coordinates.y[i] = static_cast<real>(xtc_coordinates.get()[i][1]) * 10.;
        coordinates.z[i] = static_cast<real>(xtc_coordinates.get()[i][2]) * 10.;
    }
}

void XTCTrajectoryFrame::initNumberOfAtoms(int new_number_of_atoms) {
//    assert(new_number_of_atoms >= 0);
    if (number_of_atoms != new_number_of_atoms) {
        number_of_atoms = new_number_of_atoms;
        xtc_coordinates = std::unique_ptr<XDRfile::rvec[]>(new XDRfile::rvec[number_of_atoms]);
    }
}

// ========== XTCReader ==========

namespace XDRfile {
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"
} // namespace XDRfile

XTCReader::XTCReader(const std::string& filename) : Reader(filename) {
    int number_of_atoms;
    if (XDRfile::read_xtc_natoms(const_cast<char *>(filename.c_str()), &number_of_atoms) == XDRfile::exdrOK) {
        xtc_frame = std::make_unique<XTCTrajectoryFrame>(number_of_atoms);
        xdrfile = XDRfile::xdrfile_open(filename.c_str(), "r");
    }
    if (!xtc_frame || (xdrfile == nullptr)) {
        throw std::runtime_error("xtc file " + filename + " could not be opened");
    }
}

XTCReader::~XTCReader() { XDRfile::xdrfile_close(xdrfile); }

int XTCReader::getNumberOfCoordinates() { return xtc_frame->number_of_atoms; }

bool XTCReader::readFrame() {
    return_code = XDRfile::read_xtc(xdrfile, xtc_frame->number_of_atoms, &xtc_frame->xtc_step, &xtc_frame->xtc_time,
                                    xtc_frame->xtc_box, xtc_frame->xtc_coordinates.get(), &xtc_frame->precision);
    if (return_code != XDRfile::exdrENDOFFILE && return_code != XDRfile::exdrOK) {
        throw std::runtime_error("xtc file " + filename + " could not be read (error code )");
    }
    return return_code == XDRfile::exdrOK;
}

bool XTCReader::read(TrajectoryFrame& frame) {
    if (readFrame()) {
        xtc_frame->exportFrame(frame);
        frame.number++;
        return true;
    }
    return false;
}

std::unique_ptr<Reader> createReader(const std::string& traj_name, const size_t n_atoms){
    auto extension = traj_name.substr(traj_name.find_last_of('.'));
    if(extension == ".crdbox"){
        return std::make_unique<CrdReader>(traj_name, n_atoms);
    }
    else if (extension == ".xyz"){
        return std::make_unique<XYZReader>(traj_name);
    }
    else if (extension == ".xtc"){
        return std::make_unique<XTCReader>(traj_name);
    }
    else{
        throw std::runtime_error("extension: " + extension + " not implemented yet");
    }
}

// ========== XTCWriter ==========

XTCWriter::XTCWriter(const std::string& filename)
        : xdrfile(XDRfile::xdrfile_open(filename.c_str(), "w"))
        , filename(filename) {
    if (!xdrfile) {
        throw std::runtime_error("file " + filename + " not found");
    }
}

void XTCWriter::writeFrameAt(int step, float time) {
    return_code = XDRfile::write_xtc(xdrfile, xtc_frame->number_of_atoms, step, time, xtc_frame->xtc_box,
                                     xtc_frame->xtc_coordinates.get(), xtc_frame->precision);
    if (return_code != XDRfile::exdrOK) {
        throw std::runtime_error(
                "xtc file " + filename + " could not be written!" );
    }
}

void XTCWriter::writeFrame() {
    return_code =
            XDRfile::write_xtc(xdrfile, xtc_frame->number_of_atoms, xtc_frame->xtc_step, xtc_frame->xtc_time,
                               xtc_frame->xtc_box, xtc_frame->xtc_coordinates.get(), xtc_frame->precision);
    if (return_code != XDRfile::exdrOK) {
        throw std::runtime_error(
                "xtc file " + filename + " could not be written!" );
    }
}

void XTCWriter::write(const TrajectoryFrame& frame) {
    if (!xtc_frame) {
        xtc_frame = std::make_unique<XTCTrajectoryFrame>(frame);
    }
    *xtc_frame = frame;
    writeFrame();
    step_counter = frame.step + 1;
}

void XTCWriter::writeNext(const TrajectoryFrame& frame) {
    if (!xtc_frame) {
        xtc_frame = std::make_unique<XTCTrajectoryFrame>(frame);
    }
    *xtc_frame = frame;
    writeFrameAt(step_counter, step_counter * time_delta);
    ++step_counter;
}