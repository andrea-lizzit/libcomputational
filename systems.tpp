#ifndef SYSTEMS
#define SYSTEMS

#include <cmath>
#include <numeric>
#include <libcomputational/planets.tpp>

// Coord is a coordinate in phase space
// System(coord, time) returns the first derivative of coord
template<typename Coord>
class System {
    public:
    virtual Coord operator()(Coord y, double t) = 0;
};


template<typename Coord>
class IVP {
    // Initial Value Problem
    public:
    IVP(System<Coord>& s_, Coord y0_) : s(s_), y0(y0_), t(0) { }
    IVP(System<Coord>& s_, Coord y0_, double t_) : s(s_), y0(y0_), t(t_) { }
    IVP(System<Coord>&& s_, Coord y0_) = delete;
    IVP(System<Coord>&& s_, Coord y0_, double t_) = delete;
    System<Coord>& s;
    Coord y0;
    double t;
};


typedef struct coord_t {
    double x, v;
    struct coord_t operator+(struct coord_t b)
    {
        struct coord_t res = *this;
        res.x += b.x;
        res.v += b.v;
        return res;
    }
    struct coord_t operator+=(struct coord_t b)
    {
        x += b.x;
        v += b.v;
        return *this;
    }
    template<typename Scalar>
    struct coord_t operator*(Scalar a)
    {
        struct coord_t res = *this;
        res.x *= a;
        res.v *= a;
        return res;
    }
} coord;


class Pendulum : public System<coord> {
public:
    Pendulum(double g, double l) : g_(g), l_(l) { }
    coord operator()(coord y, double t)
    {
        coord res;
        res.x = y.v;
        res.v = -g_ / l_ * std::sin(y.x);
        return res;
    }
    coord operator()(coord y) // this system has a time-independent Hamiltonian
    {
        return operator()(y, 0);
    }
private:
    double g_, l_;
};


class HarmonicOsc : public System<coord> {
public:
    HarmonicOsc(double k, double m) : k_(k), m_(m) { }
    coord operator()(coord y, double t)
    {
        coord res;
        res.x = y.v;
        res.v = -k_/m_ * y.x;
        return res;
    }
    coord operator()(coord y) // this system has a time-independent Hamiltonian
    {
        return operator()(y, 0);
    }
private:
    double k_, m_;
};


class CelestialSystem : public System<cbody*> {
public:
    CelestialSystem(size_t size_, double G=6.674e-11) : G_(G), size(size_) {
        // heap allocation only at the construction of the object
        // avoids malloc overhead during evaluation of the system
        // should not slow down execution
        res = new cbody[size_];
    }
    cbody* operator()(cbody *bodies, double t)
    {
        for (size_t i = 0; i < size; ++i) {
            res[i].M = bodies[i].M;
            res[i].pos = bodies[i].vel;
            res[i].vel = {0, 0, 0};
            for (size_t j = 0; j < size; ++j) {
                if (i == j)
                    continue;
                res[i].vel += force(bodies[i], bodies[j]);
            }
            res[i].vel *= 1/bodies[i].M;
        }
        return res;
    }
    vec3<double> force(cbody b1, cbody b2) // force of b2 on b1
    {
        double distance = (b1.pos - b2.pos).length();
        double magnitude = G_ * b1.M * b2.M / (distance * distance);
        vec3<double> f = (b2.pos - b1.pos) * magnitude / distance;
        return f;
    }
    double energy(cbody *bodies)
    {
        double e = potential_energy(bodies) + kinetic_energy(bodies);
        return e;
    }
    double kinetic_energy(cbody* bodies)
    {
        double e = 0;
        for (size_t i = 0; i < size; ++i) {
            e += 0.5 * bodies[i].M * std::pow(bodies[i].vel.length(), 2);
        }
        return e;
    }
    double potential_energy(cbody *bodies)
    {
        double e = 0;
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < i; ++j) {
                // calculate reciprocal potential energy for every couple of bodies
                double distance = (bodies[i].pos - bodies[j].pos).length();
                e -= G_*bodies[i].M*bodies[j].M / distance;
            }
        }
        return e;
    }
    vec3<double> momentum(cbody *bodies)
    {
        vec3<double> m = {0.0, 0.0, 0.0};
        for (size_t i = 0; i < size; ++i) {
            m += bodies[i].momentum();
        }
        return m;
    }
    double mass(cbody *bodies)
    {
        double M;
        for (size_t i = 0; i < size; ++i) {
            M += bodies[i].M;
        }
        return M;
    }

    double G_;
    cbody* res;
    size_t size;
};

#endif