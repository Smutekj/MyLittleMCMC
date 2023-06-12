#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>


class TrajectoryFrame;

class Reader{
public:
    std::string file_name;
    TrajectoryFrame timestep;
    std::vector<size_t> n_atoms_in_frame;
    size_t n_frames = 0;

protected:
    std::ifstream file;
    bool first_frame_read = false;
public:
    explicit Reader(const std::string& file_name) : file_name(file_name) {}

    virtual ~Reader(){
        if(file.is_open()){
            file.close();
        }
    }
    virtual bool read(TrajectoryFrame& frame) = 0;

//    virtual const TrajectoryFrame& next() = 0;

//    virtual bool isThereNext() const = 0;

};


class CrdReader : public Reader{

public:
   std::vector<size_t> offsets;
   size_t total_atoms;

    CrdReader(const std::string& file_name, size_t n_atoms) : Reader(file_name), total_atoms(n_atoms) {
        auto extension_pos = file_name.find_last_of('.');
        if(file_name.substr(extension_pos) != ".crdbox"){
            throw std::runtime_error( file_name + std::string(" file does not have .crdbox extension!"));
        }
        initialize();
    }

    ~CrdReader(){
        if(first_frame_read) {
            file.close();
        }
    }

    bool read(TrajectoryFrame& frame) override;

//    bool isThereNext() const override;

//    const TrajectoryFrame& next() override;

private:
    bool updateBoxAndCoordinates(TrajectoryFrame& frame, size_t offset) ;
    void initialize();
};

class XYZReader : public Reader{

public:
    std::vector<size_t> offsets;

    XYZReader(const std::string& file_name) : Reader(file_name) {
        initialize();
    }

    ~XYZReader(){
        if(first_frame_read) {
            file.close();
        }
    }

    bool read(TrajectoryFrame& frame) override;

//    const TrajectoryFrame& next() override;

//    bool isThereNext() const override;

private:
    std::vector<std::string> separateLine(std::string& line, std::string&& delimiter) ;
    void initialize();
    bool updateBoxAndCoordinates(size_t offset);
};

namespace XDRfile {
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"
} // namespace XDRfile


/**
 * @brief Base data structure for native XTC format as used in Gromacs and supported by the C library. Import methods
 * do the data conversion from the Faunus native format to the XTC format, and export methods do the oposite.
 *
 * By convention, XTC has coordinates' origin in a corner (main box's coordinates are always positive), while Faunus
 * has coordinates' origin in the center of the simulation box. During conversion the corresponding offset is
 * subtracted, or added, respectively.
 *
 * XTC format uses floats (Faunus reals) and dimensions are in nanometers (Faunus ångströms). XTC library requires
 * raw 2D C-style arrays for coordinates in row-major format. XTC tensor of the simulation box is converted
 * to an XYZ point pressuming orthogonal geometry bacause of current limitation of Faunus. XTC format does not support
 * variable number of coordinates (atoms) between frames.
 */
struct XTCTrajectoryFrame {
    XDRfile::matrix xtc_box;                          //!< box tensor; only diagonal elements are used
    std::unique_ptr<XDRfile::rvec[]> xtc_coordinates; //!< C-style array of particle coordinates
    int xtc_step = 0;                                 //!< current step number
    float xtc_time = 0.0;                             //!< current time (unit?)
    int number_of_atoms = 0;                          //!< number of coordinates (atoms) in each step
    float precision = 1000.0;                         //!< output precision
    /**
     * @brief Creates an empty XTC trajectory step for given number of coordinates (atoms).
     * @param number_of_atoms  number of coordinates (atoms)
     */
    explicit XTCTrajectoryFrame(int number_of_atoms);
    /**
     * @brief Creates an XTC trajectory step from TrajectoryFrame and converts data accordingly.
     * @param frame  source trajectory step
     */
    explicit XTCTrajectoryFrame(const TrajectoryFrame& frame);
    /**
     * @brief Creates an XTC trajectory step from scalar parameters and input iterator and converts data accordingly.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[in] step  step step
     * @param[in] time  timestamp in picoseconds
     * @param[in] box  box dimensions (xyz) in nanometers
     * @param[in] coordinates_begin  input iterator with coordinates in nanometers
     * @param[in] coordinates_end  input iterator's end
     */
    template <class begin_iterator, class end_iterator>
    XTCTrajectoryFrame(int step, float time, const std::vector<real>& box, begin_iterator coordinates_begin,
                       end_iterator coordinates_end) {
        initNumberOfAtoms(std::distance(coordinates_begin, coordinates_end));
        importFrame(step, time, box, coordinates_begin, coordinates_end);
    }
    /**
     * @brief Copies TrajectoryFrame and converts data. Calls importFrame.
     *
     * The number of coordinates is immutable to prevent mistakes.
     * @param frame  source step
     * @throw std::runtime_error  when the number of coordinates does not match
     */
    XTCTrajectoryFrame& operator=(const TrajectoryFrame& frame);
    /**
     * @brief Imports data from a TrajectoryFrame.
     * @param frame  source step
     * @throw std::runtime_error  when the number of coordinates does not match
     */
    void importFrame(const TrajectoryFrame& frame);
    /**
     * @brief Imports data from scalar parameters and an input iterator over coordinates.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[in] step  step step
     * @param[in] time  step timestamp
     * @param[in] box  box dimensions (xyz) in nanometers
     * @param[in] coordinates_begin  input iterator with coordinates in nanometers
     * @param[in] coordinates_end  input iterator's end
     * @throw std::runtime_error when the number of coordinates does not match
     */
    template <class begin_iterator, class end_iterator>
    void importFrame(const int step, const float time, const std::vector<real>& box, begin_iterator coordinates_begin,
                     end_iterator coordinates_end) {
        importTimestamp(step, time);
        importBox(box);
        importCoordinates(coordinates_begin, coordinates_end);
    }
    /**
     * @brief Exports data from a TrajectoryFrame.
     * @param frame  target step
     * @throw std::runtime_error  when the number of coordinates does not match
     */
    void exportFrame(TrajectoryFrame& frame) const;
    /**
     * @brief Exports data to scalar paramers and output iterator over atomic coordinates.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[out] step  step step
     * @param[out] time  step timestamp
     * @param[out] box  box dimensions (xyz) in nanometers
     * @param[out] coordinates_begin  output iterator with coordinates in nanometers
     * @param[out] coordinates_end  output iterator's end
     * @throw std::runtime_error  when the number of coordinates does not match
     */
    template <class begin_iterator, class end_iterator>
    void exportFrame(size_t& step, float& time, real box[3], begin_iterator coordinates_begin,
                     end_iterator coordinates_end) const {
        exportTimestamp(step, time);
        exportBox(box);
        exportCoordinates(coordinates_begin, coordinates_end);
    }

protected:
    using XTCFloat = float;
    /**
     * @brief Imports and converts step and timestamp.
     * @param[in] step  step step
     * @param[in] time  step timestamp
     */
    void importTimestamp(int step, float time);
    /**
     * @brief Imports and converts simulation box dimensions.
     * @param[in] box  simulation box dimensions in nanometers (xyz)
     */
    void importBox(const std::vector<real>& box);
    /**
     * @brief Imports and converts atomic coordinates. Offset is added to all coordinates to account different
     * coordinates' origin.
     * @param[in] coordinates  atomic coordinates in nanometers
     * @param[in] offset  offset in nanometers to add to all coordinates upon conversion
     */
    void importCoordinates(const Coordinates& coordinates);
    /**
     * @brief Exports and converts step and timestamp.
     * @param[out] step  step step
     * @param[out] time  step timestamp
     */
    void exportTimestamp(size_t& step, float& time) const;
    /**
     * @brief Exports and converts simulation box dimensions.
     * @param[out] box  simulation box dimensions in nanometers (xyz)
     */
    void exportBox(real box[3]) const;
    /**
     * @brief Exports and converts atomic coordinates. Offset is subtracted from all coordinates to account different
     * coordinates' origin.
     * @param[out] coordinates  atomic coordinates in nanometers
     * @param[in] offset  offset in nanometers to subtract from all coordinates upon conversion
     * @throw std::runtime_error  when the source box is not orthogonal
     */
    void exportCoordinates(Coordinates& coordinates) const;

private:
    /**
     * @brief Set number of coordinates (atoms) and allocate memory for them.
     */
    void initNumberOfAtoms(int);
};

/**
 * @brief Reads frames from an XTC file (GROMACS compressed trajectory file format). It is a wrapper around
 * C function calls.
 *
 * Frames are stored into a TrajectoryFrame structure or as a list of positions in an output iterator. The class
 * is responsible for I/O operations, not data conversion. For details about data conversion XTCTrajectoryFrame.
 */
class XTCReader : public Reader {
    int return_code = XDRfile::exdrOK;   //!< last return code of a C function
    XDRfile::XDRFILE* xdrfile = nullptr; //!< file handle
    //! data structure for C functions; the number of coordinates is immutable
    std::unique_ptr<XTCTrajectoryFrame> xtc_frame;

public:
    std::string filename; //!< name of the trajectory file, mainly for error reporting
    /**
     * @param filename  a name of the XTC file to open
     */
    explicit XTCReader(const std::string& filename);
    ~XTCReader() override;
    /**
     * @brief Returns number of coordinates (atoms) in each step. Immutable during object lifetime.
     * @return number of coordinates (atoms)
     */
    int getNumberOfCoordinates();
    /**
     * @brief Reads the next step in the trajectory from a file
     * @param[out] frame   target step
     * @return true on success, false at the end of file
     * @throw std::runtime_error  when other I/O error occures
     */
    bool read(TrajectoryFrame& frame) override;
    /**
     * @brief Reads the next step in the trajectory.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[out] step  step step
     * @param[out] time  step timestamp in picoseconds
     * @param[out] box  box dimensions (xyz) in nanometers
     * @param[out] coordinates_begin  output iterator to store coordinates
     * @param[out] coordinates_end  output iterator's end
     * @return true on success, false at the end of file
     * @throw std::runtime_error  when other I/O error occures
     */
    template <class begin_iterator, class end_iterator>
    bool read(int& step, float& time, std::vector<real>& box, begin_iterator coordinates_begin, end_iterator coordinates_end) {
        if (readFrame()) {
            xtc_frame->exportFrame(step, time, box, coordinates_begin, coordinates_end);
            return true;
        }
        return false;
    }

protected:
    /**
     * @brief Actual wrapper around C function that reads the next step into xtc_frame.
     * @return true on success, false at the end of file
     * @throw std::runtime_error  when other I/O error occures
     */
    bool readFrame();
};

/**
 * @brief Writes frames into an XTC file (GROMACS compressed trajectory file format). It is a wrapper around
 * C function calls.
 *
 * The frames can be provided as a TrajectoryFrame structure or as a list of positions in an input iterator. The class
 * is responsible for I/O operations, not data conversion.
 */
class XTCWriter {
    int return_code = XDRfile::exdrOK;             //!< last return code of a C function
    XDRfile::XDRFILE* xdrfile = nullptr; //!< file handle
    std::unique_ptr<XTCTrajectoryFrame> xtc_frame; //!< data structure for C functions;
    //!< the number of coordinates is immutable
    int step_counter = 0;                          //!< frame counter for automatic increments
    float time_delta = 1.0;                     //!< timestamp of a frame is computed as step * time_delta

public:
    const std::string filename; //!< name of the trajectory file, mainly for error reporting
    /**
     * @param filename  a name of the XTC file to open
     */
    explicit XTCWriter(const std::string& filename);
    /**
     * @brief Writes a frame into the file.
     * @param[in] frame  frame to be written
     * @throw std::runtime_error  when other I/O error occures
     */
    void write(const TrajectoryFrame& frame);
    /**
     * @brief Writes a next frame into the file using own automatic counter for step and timestamp.
     * The corresponding values in the frame are ignored.
     * @param[in] frame  frame to be written
     * @throw std::runtime_error  when other I/O error occures
     */
    void writeNext(const TrajectoryFrame& frame);
    /**
     * @brief Writes a next frame into the file using own automatic counter for step and timestamp.
     * @tparam begin_iterator
     * @tparam end_iterator
     * @param[in] box  dimensions of the cubic box (xyz)
     * @param[in] coordinates_begin  input iterator with coordinates (not particles)
     * @param[in] coordinates_end  input iterator's end
     */
    template <class begin_iterator, typename end_iterator>
    void writeNext(const std::vector<real>& box, const size_t n_atoms, begin_iterator coordinates_begin, end_iterator coordinates_end) {
        if (!xtc_frame) {
            xtc_frame = std::make_unique<XTCTrajectoryFrame>(n_atoms);
        }
        xtc_frame->importFrame(step_counter, static_cast<float>(step_counter) * time_delta, box, coordinates_begin,
                               coordinates_end);
        writeFrame();
        ++step_counter;
    }

protected:
    /**
     * @brief Actual wrapper around C function that writes the current frame into xtc_frame.
     * @throw std::runtime_error  when other I/O error occures
     */
    void writeFrame();
    /**
     * @brief Actual wrapper around C function that writes the current frame into xtc_frame overriding step
     * and timestamp.
     * @throw std::runtime_error  when other I/O error occures
     */
    void writeFrameAt(int step, float time);
};

std::unique_ptr<Reader> createReader(const std::string& traj_name, size_t n_atoms);
