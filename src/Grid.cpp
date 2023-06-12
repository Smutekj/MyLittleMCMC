//
// Created by smutekj on 23.04.23.
//

#include "Grid.h"
#include <set>
#include <vector>
#include <unordered_map>
#include <map>
#include <nlohmann/json.hpp>



Grid::Grid(size_t n_cellx_x, real box_size) : nx(n_cellx_x), ny(n_cellx_x), nz(n_cellx_x),
                cell_size(box_size/(real)n_cellx_x), box_size({box_size, box_size, box_size}){
    cell2nearest_neighbours.resize(nx*ny*nz);
    cell2forward_nearest_neighbours.resize(nx*ny*nz);

    for( int ix = 0; ix < nx; ++ix) {
        for( int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                std::set<int> uff;
                const int cell_index = ix + iy*nx + iz*ny*nx;
                for (int i = 0; i < 27; ++i) {
                    int dnx0 = i % 3 - 1;
                    int dny0 = ((i / 3) % 3) - 1;
                    int dnz0 = ((i / 9) % 9) - 1;

                    int dnx = ix + dnx0;
                    int dny = iy + dny0;
                    int dnz = iz + dnz0;

                    const int nx2 = dnx + ((dnx<0) - (dnx >= (int)nx))*nx;
                    const int ny2 = dny + ((dny<0) - (dny >= (int)ny))*ny;
                    const int nz2 = dnz + ((dnz<0) - (dnz >= (int)nz))*nz;
                    assert(nx2>=0 and ny2 >= 0 and nz2 >= 0);
                    auto pico = nx2 + ny2*nx + nz2*nx*ny;
                    uff.insert(pico);

                    cell2nearest_neighbours[cell_index].push_back(nx2 + ny2*nx + nz2*nx*ny);
                    if(dnx0 >=0 and dny0>=0 and dnz0>=0){
                        cell2forward_nearest_neighbours[cell_index].push_back(nx2 + ny2*nx + nz2*nx*ny);
                    }
                }
            }
        }
    }
}

size_t Grid::coordToCell(const real x, const real y, const real z) const{
    const size_t ix = static_cast<size_t>(std::floor(x/cell_size))%nx;
    const size_t iy = static_cast<size_t>(std::floor(y/cell_size))%ny;
    const size_t iz = static_cast<size_t>(std::floor(z/cell_size))%nz;

    assert(ix + iy * nx + iz * nx*ny < nx*ny*nz);
    return ix + iy * nx + iz * nx*ny;
}

size_t Grid::coordToCell(const gmx::RVec r_coord) const{
    const size_t ix = static_cast<size_t>(std::floor(r_coord[XX]/cell_size))%nx;
    const size_t iy = static_cast<size_t>(std::floor(r_coord[YY]/cell_size))%ny;
    const size_t iz = static_cast<size_t>(std::floor(r_coord[ZZ]/cell_size))%nz;

    assert(ix + iy * nx + iz * nx*ny < nx*ny*nz);
    return ix + iy * nx + iz * nx*ny;
}

[[nodiscard]] size_t Grid::cellIndex(const size_t ix, const size_t iy, const size_t iz) const{
    return ix + iy * nx + iz * nx*ny;
}


size_t Grid::cellCoordX(const size_t cell_index)const{
    return cell_index % nx ;
}

size_t Grid::cellCoordY(const size_t cell_index)const{
    return (cell_index/nx) % ny ;
}

size_t Grid::cellCoordZ(const size_t cell_index)const{
    return (cell_index/(nx*ny) % nz) ;
}


size_t Grid::cellCoordX(const gmx::RVec r_coord)const{
    return static_cast<size_t>(r_coord[XX]/cell_sizexy) ;
}

size_t Grid::cellCoordY(const gmx::RVec r_coord)const{
    return static_cast<size_t>(r_coord[YY]/cell_sizexy) ;
}
size_t Grid::cellCoordZ(const gmx::RVec r_coord)const{
    return static_cast<size_t>(r_coord[ZZ]/cell_size) ;
}

gmx::IVec Grid::cellCoords(const gmx::RVec r_coord) const{
    const auto cell_coord_x =  static_cast<int>(r_coord[XX]/cell_sizexy) ;
    const auto cell_coord_y = static_cast<int>(r_coord[YY]/cell_sizexy) ;
    const auto cell_coord_z = static_cast<int>(r_coord[ZZ]/cell_size) ;

    return {cell_coord_x, cell_coord_y, cell_coord_z} ;
}


size_t SearchGrid::coordToCell(const real x, const real y, const real z) const{
    const size_t ix = static_cast<size_t>(std::floor(x/cell_size_[XX]))%n_cells_[XX];
    const size_t iy = static_cast<size_t>(std::floor(y/cell_size_[YY]))%n_cells_[YY];
    const size_t iz = static_cast<size_t>(std::floor(z/cell_size_[ZZ]))%n_cells_[ZZ];

    assert(ix + iy * n_cells_[XX] + iz * n_cells_[XX]*n_cells_[YY] < n_cells_[XX]*n_cells_[YY]*n_cells_[ZZ]);
    return ix + iy * n_cells_[XX] + iz * n_cells_[XX]*n_cells_[YY];
}

size_t SearchGrid::coordToCell(const gmx::RVec r_coord) const{
    const size_t ix = static_cast<size_t>(std::floor(r_coord[XX]/cell_size_[XX]))%n_cells_[XX];
    const size_t iy = static_cast<size_t>(std::floor(r_coord[YY]/cell_size_[YY]))%n_cells_[YY];
    const size_t iz = static_cast<size_t>(std::floor(r_coord[ZZ]/cell_size_[ZZ]))%n_cells_[ZZ];

    assert(ix + iy * n_cells_[XX] + iz * n_cells_[XX]*n_cells_[YY] < n_cells_[XX]*n_cells_[YY]*n_cells_[ZZ]);
    return ix + iy * n_cells_[XX] + iz * n_cells_[XX]*n_cells_[YY];
}

[[nodiscard]] size_t SearchGrid::cellIndex(const size_t ix, const size_t iy, const size_t iz) const{
    return ix + iy * n_cells_[XX] + iz * n_cells_[XX]*n_cells_[YY];
}


size_t SearchGrid::cellCoordX(const size_t cell_index)const{
    return cell_index % n_cells_[XX] ;
}

size_t SearchGrid::cellCoordY(const size_t cell_index)const{
    return (cell_index/n_cells_[XX]) % n_cells_[YY] ;
}

size_t SearchGrid::cellCoordZ(const size_t cell_index)const{
    return (cell_index/(n_cells_[XX]*n_cells_[YY]) % n_cells_[ZZ]) ;
}

gmx::IVec SearchGrid::cellCoord(const size_t cell_index)const{
    return { static_cast<int>(cellCoordX(cell_index)),
             static_cast<int>(cellCoordY(cell_index)),
             static_cast<int>(cellCoordZ(cell_index))};
}

size_t SearchGrid::cellCoordX(const gmx::RVec r_coord)const{
    return static_cast<size_t>(r_coord[XX]/cell_size_[XX]) ;
}
size_t SearchGrid::cellCoordY(const gmx::RVec r_coord)const{
    return static_cast<size_t>(r_coord[YY]/cell_size_[YY]) ;
}
size_t SearchGrid::cellCoordZ(const gmx::RVec r_coord)const{
    return static_cast<size_t>(r_coord[ZZ]/cell_size_[ZZ]) ;
}

gmx::IVec SearchGrid::cellCoords(const gmx::RVec r_coord) const{
    const auto cell_coord_x = static_cast<int>(r_coord[XX]/cell_size_[XX]) ;
    const auto cell_coord_y = static_cast<int>(r_coord[YY]/cell_size_[YY]) ;
    const auto cell_coord_z = static_cast<int>(r_coord[ZZ]/cell_size_[ZZ]) ;

    return {cell_coord_x, cell_coord_y, cell_coord_z} ;
}

SearchGrid::SearchGrid(gmx::RVec box_size, gmx::RVec cell_size) :
    box_size_(box_size), cell_size_(cell_size)
{
    n_cells_[XX] = std::floor(box_size_[XX]/cell_size_[XX]);
    n_cells_[YY] = std::floor(box_size_[YY]/cell_size_[YY]);
    n_cells_[ZZ] = std::floor(box_size_[ZZ]/cell_size_[ZZ]);
    cell_size_[XX] = box_size_[XX]  / static_cast<real>(n_cells_[XX]);
    cell_size_[YY] = box_size_[YY]  / static_cast<real>(n_cells_[YY]);
    cell_size_[ZZ] = box_size_[ZZ]  / static_cast<real>(n_cells_[ZZ]);
    cell2nearest_neighbours.resize(n_cells_[XX]*n_cells_[YY]*n_cells_[ZZ]);
}

