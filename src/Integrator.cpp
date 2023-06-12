//
// Created by smutekj on 16.07.22.
//

#include "Integrator.h"

size_t cList::modulo(const int i, const int n){
    assert(n>0);
    return ((i%n)+n)%n; }


void cList::createCellList(const Geometry& spc){
    const auto active_particles = spc.activeParticles();
    const auto n = spc.particles.size();
    particle2cell_pos.resize(n);

    const auto grid_box = spc.geometry.getLength();
    std::for_each(active_particles.begin(), active_particles.end(),
                  [&, half_box = 0.5 * grid_box](const Particle& particle) {
                      insertMember(indexOf(particle), particle.pos + half_box);
                  });
    this->grid_box = grid_box;
}


void cList::insertMember(index_type particle_index, const Point& coord){
    const auto cell_index = cellIndexAt(coord);
    particle2cell_pos.at(particle_index).first = cell_index;
    particle2cell_pos.at(particle_index).second = cell2particles.at(cell_index).size();
    cell2particles.at(cell_index).push_back(particle_index);
}
