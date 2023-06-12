//
// Created by smutekj on 23.04.23.
//

#ifndef GR_GCMC_CLION_GRID_H
#define GR_GCMC_CLION_GRID_H


#include <stdexcept>
#include "core.h"


struct Grid{
    real cell_size=10.0;
    real cell_sizexy=10.0;
    gmx::RVec box_size;
    size_t nx, ny, nz;

    std::vector<std::vector<size_t>> cell2nearest_neighbours;
    std::vector<std::vector<size_t>> cell2forward_nearest_neighbours;

    Grid(size_t n_cells_x, real box_size);

    [[nodiscard]] size_t coordToCell( real x,  real y,  real z) const;
    [[nodiscard]] size_t coordToCell( gmx::RVec r_coord) const;
    [[nodiscard]] size_t cellIndex(size_t nx,  size_t ny,  size_t nz) const;

    [[nodiscard]] size_t cellCoordX( size_t cell_index)const;
    [[nodiscard]] size_t cellCoordY( size_t cell_index)const;
    [[nodiscard]] size_t cellCoordZ( size_t cell_index)const;

    [[nodiscard]] size_t cellCoordX(gmx::RVec r_coord)const;
    [[nodiscard]] size_t cellCoordY(gmx::RVec r_coord)const;
    [[nodiscard]] size_t cellCoordZ(gmx::RVec r_coord)const;

    [[nodiscard]] gmx::IVec cellCoords(gmx::RVec r_coord) const;

};


class SearchGrid{

public:

    gmx::RVec cell_size_;
    gmx::RVec box_size_;
    gmx::IVec n_cells_;
    std::vector<std::vector<size_t>> cell2nearest_neighbours;

    SearchGrid(gmx::RVec box_size, gmx::RVec cell_size);

    [[nodiscard]] size_t coordToCell( gmx::RVec r_coord) const;
    [[nodiscard]] size_t cellIndex(size_t nx,  size_t ny,  size_t nz) const;


    [[nodiscard]] size_t cellCoordX( size_t cell_index)const;
    [[nodiscard]] size_t cellCoordY( size_t cell_index)const;
    [[nodiscard]] size_t cellCoordZ( size_t cell_index)const;
    [[nodiscard]] gmx::IVec cellCoord( size_t cell_index) const;

    [[nodiscard]] size_t cellCoordX(gmx::RVec r_coord)const;
    [[nodiscard]] size_t cellCoordY(gmx::RVec r_coord)const;
    [[nodiscard]] size_t cellCoordZ(gmx::RVec r_coord)const;

    [[nodiscard]] gmx::IVec cellCoords(gmx::RVec r_coord) const;


    size_t coordToCell(const real x, const real y, const real z) const;
};




#endif //GR_GCMC_CLION_GRID_H
