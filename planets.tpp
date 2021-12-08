#ifndef SYSTEMS_HPP
#define SYSTEMS_HPP

#include <libcomputational/particles.hpp>

template<typename scalar=double>
struct cbody_t : public particle<scalar> {
    cbody_t() { }
    scalar M;
    vec3<scalar> momentum()
    {
        return this->vel * this->M;
    }
};
typedef cbody_t<double> cbody; // another fix for celestialsystem


#endif