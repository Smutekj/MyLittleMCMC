#include "System.h"
#include "core.h"
#include "vectypes.h"
#ifndef GR_GCMC_CLION_INTEGRATOR_H
#define GR_GCMC_CLION_INTEGRATOR_H

struct Point{
    double x;
    double y;
    double z;
};

class Grid {
    using cell_index = size_t;

    int n_cellsx;
    int n_cellsy;
    int n_cellsxy;
    int n_cellsz;
    int n_cells;

    double dx;
    double dy;
    double dz;

    double Lx;
    double Ly;
    double Lz;

    std::vector<std::vector<cell_index>> neighbouring_cell_indices;

    struct cell_point{
        int ix;
        int iy;
        int iz;
    };


private:
    size_t modulo (const int i, const int n);

public:
    Grid(int n_cellsx, int n_cellsy, int n_cellsz, double dx, double dy, double dz) :
            n_cellsx(n_cellsx), n_cellsy(n_cellsy), n_cellsxy(n_cellsx*n_cellsy), n_cellsz(n_cellsz),
            n_cells(n_cellsx * n_cellsy * n_cellsz),
            dx(dx), dy(dy), dz(dz), Lx(dx*n_cellsx), Ly(dy*n_cellsy), Lz(dz*n_cellsz),
            neighbouring_cell_indices(n_cells)
    {
        for( cell_index index=0; index<n_cells; index++) {
            const auto c_point = gridCellCoordsOf(index);

            for (auto i = -1; i <= 1; ++i) {
                for (auto j = -1; j <= 1; ++j) {
                    for (auto k = -1; k <= 1; ++k) {
                        const auto new_cell_index = indexOfCellPoint(modulo(c_point.ix + i, n_cellsx),
                                                                     modulo(c_point.iy + j, n_cellsy),
                                                                     modulo(c_point.iz+ k, n_cellsz));
                        neighbouring_cell_indices[index].push_back(new_cell_index);
                    }
                }
            }
        }
    }
    Grid(int n_cellsx, double dx) :
            Grid(n_cellsx, n_cellsx, n_cellsx, dx, dx, dx)
    {}

    cell_point gridCellCoordsAt(const Point& point){
        cell_point res;
        res.ix =std::floor(point.x/dx);
        res.iy =std::floor(point.y/dy);
        res.iz =std::floor(point.z/dz);
        return res;
    }

    cell_point gridCellCoordsOf(const cell_index index){
        cell_point res;
        res.iz = index/(n_cellsxy);
        res.iy = (index/n_cellsx)%n_cellsy;
        res.ix = (index%n_cellsx);
        return res;
    }

    cell_index indexOfCellPoint(const size_t nx, const size_t ny, const size_t nz){
        return  nx + ny*n_cellsx + nz*n_cellsxy;
    }
    cell_index cellIndexAt(const Point& point){
        const int nx = std::floor(point.x/dx);
        const int ny = std::floor(point.y/dy);
        const int nz = std::floor(point.z/dz);
        return nx + ny*n_cellsx + nz*n_cellsxy;
    }

};


class cList {
    using cell_index = size_t;
    using index_type = size_t;

public:

    Grid grid;

    Point simulation_box;
    Point grid_box;

    std::vector<std::vector<index_type>> cell2particles;


public:
    cList(int n_cellsx, int n_cellsy, int n_cellsz, double dx, double dy, double dz)
       : grid(n_cellsx, n_cellsy, n_cellsz, dx, dy, dz)
    {}
    cList(int n_cellsx, double dx) :
            cList(n_cellsx, n_cellsx, n_cellsx, dx, dx, dx)
    {}

    void insertMember(index_type particle_index, const Point& coord);

    const std::vector<size_t>& particlesInCell(cell_index index) const{
        return cell2particles.at(index);
    }
};


class PairList{
    std::vector<size_t> neighbours;
    double buffer_size;
    double r_cut;

public:
    PairList(double r_cut, double buffer_size) :
    neighbours(0), buffer_size(buffer_size), r_cut(r_cut) {}

    const std::vector<size_t>& getNeighbours() const{
        return neighbours;
    }

    void addToList(const size_t index){
        neighbours.push_back(index);
    }
};


struct force{
    double fx;
    double fy;
    double fz;
};


double LJ_force(double r2, double A, double B) {

    const auto r4 = (r2*r2);
    const double r8 = r4*r4;
    const double r14 = r4*r4*r4*r2;
    return - 12.* A/r14 + 6.*B/r8;
}



force LJ_force(double* xs, double* ys, double* zs, double x, double y, double z, double* sigs, double* eps,
               int n_particles, double box_length){
    force f;
    double half_box_length = box_length/2.;

#pragma clang loop vectorize(enable)
    for (auto i = 0; i < n_particles; ++i) {
        const auto dx = distancePBC(x, xs[i], half_box_length, box_length);
        const auto dy = distancePBC(y, ys[i], half_box_length, box_length);
        const auto dz = distancePBC(z, zs[i], half_box_length, box_length);
        const auto r2 = dx * dx + dy * dy + dz * dz;

        double f_size = LJ_force(r2, 1., 1.);
        f.fx += dx * f_size;
        f.fy += dy * f_size;
        f.fz += dz * f_size;
    }
    return f;
}

class Integrator {

    double dt; //! ps
    System& system;
    std::string name;

    void step(double* xs, double* ys, double* zs,
              double* vxs, double* vys, double* vzs,
              double* axs, double* ays, double* azs){
        for(int i = 0; i< system.n_atoms; ++i){
            xs[i] += vxs[i] * dt + 0.5 * dt * dt * axs[i];
            ys[i] += vys[i] * dt + 0.5 * dt * dt * ays[i];
            zs[i] += vzs[i] * dt + 0.5 * dt * dt * azs[i];
        }

        for(int i = 0; i< system.n_atoms; ++i){
            vxs[i] += vxs[i] + dt * axs[i];
            vys[i] += vys[i] + dt * ays[i];
            vzs[i] += vzs[i] + dt * azs[i];
        }
    }
};


#endif //GR_GCMC_CLION_INTEGRATOR_H
