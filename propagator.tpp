#ifndef PROPAGATOR
#define PROPAGATOR

#include <cstddef>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <libcomputational/planets.tpp>
#include <libcomputational/systems.tpp>
using std::size_t;

template<typename Coord>
class Propagator {
public:
    Propagator(IVP<Coord>& problem_, double timestep_) : problem(problem_), timestep(timestep_), t(problem_.t), yi(problem_.y0) { }
    // declare public current state yi because in general it might be different from subclass-specific variable "state"
    Coord yi;
    double t; // keep track of time for time-dependent Hamiltonians
    virtual Coord step() = 0; // return next value of the system coordinate and update internal state
    Coord propagate(int steps)
    {
        // initialize the variable in case steps == 0
        Coord yi = problem.y0;
        for (int i = 0; i < steps; ++i) {
            yi = step();
        }
        return yi;
    }
protected:
    IVP<Coord>& problem;
    double timestep;
};

template<typename Coord>
class Prop_EEuler : public Propagator<Coord> {
    public:
    Prop_EEuler(IVP<Coord>& problem_, double timestep_) : Propagator<Coord>(problem_, timestep_) {
        state = this->problem.y0;
    }
    Coord step()
    {
        Coord dy = this->problem.s(state, this->t) * this->timestep;
        state += dy;
        this->yi = state;
        this->t += this->timestep;
        return state;
    }
    private:
    Coord state;
};


// warning: this implementation is done ad hoc for the CelestialSystem
// it will not work according to the general interface
class Prop_VelVer : public Propagator<cbody*> {
public:
    Prop_VelVer(IVP<cbody*>& problem_, double timestep_, size_t size_) : 
                Propagator<cbody*>(problem_, timestep_), size(size_) {
        state = this->problem.y0;
        // heap allocation only at the construction of the object
        // avoids malloc overhead during evaluation of the system
        // should not slow down execution
        res = new cbody[size_];
    }
    cbody* step()
    {
        // calculate dw_n which contains f_n
        // f_n is in dw_n[i].vel
        cbody dw_n[size];
        cbody* src = this->problem.s(state, this->t);
        std::copy(src, src+size, dw_n);
        for (size_t i = 0; i < size; ++i) {
            // calculate r_i,n+1
            state[i].pos += dw_n[i].pos * this->timestep + 
                dw_n[i].vel * 1/(2) * std::pow(this->timestep, 2);
        }
        
        // calculate dw_n1 which contains f_n+1
        cbody dw_n1[size];
        src = this->problem.s(state, this->t);
        std::copy(src, src+size, dw_n1);
        for (size_t i = 0; i < size; ++i) {
            // calculate v_i,n+1
            state[i].vel += (dw_n[i].vel + dw_n1[i].vel) * 1 / (2)* this->timestep;
        }
        this->t += this->timestep;
        return state;
    }
private:
    size_t size;
    cbody* res;
    cbody* state;
};

template<typename Coord>
class Prop_RK : public Propagator<Coord> {
    public:
    Prop_RK(IVP<Coord>& problem_, double timestep_) : Propagator<Coord>(problem_, timestep_) {
        state = this->problem.y0;
    }
    Coord step()
    {
        Coord y1 = state;
        Coord y2 = state + this->problem.s(y1, this->t) * (this->timestep / 2);
        Coord y3 = state + this->problem.s(y2, this->t + this->timestep / 2) * (this->timestep / 2);
        Coord y4 = state + this->problem.s(y3, this->t + this->timestep / 2) * this->timestep;
        Coord y = state + (
                    this->problem.s(y1, this->t) + 
                    this->problem.s(y2, this->t + this->timestep / 2) * 2 +
                    this->problem.s(y3, this->t + this->timestep / 2) * 2 +
                    this->problem.s(y4, this->t)
                    ) * (this->timestep / 6);
        state = y;
        this->yi = state;
        this->t += this->timestep;
        return state;
    }
    private:
    Coord state;
};

template<typename Coord>
class Prop_Numerov : public Propagator<Coord> {
    public:
    Prop_Numerov(NumerovIVP<Coord>& problem_, double dx_) : Propagator<Coord>(problem_, dx_) {
        NumerovIVP<Coord>& problem = static_cast<NumerovIVP<Coord>&>(this->problem);
        state = problem.y0;
        yim1 = problem.ym1;
    }
    Coord step()
    {
        NumerovIVP<Coord>& problem = static_cast<NumerovIVP<Coord>&>(this->problem);
        NumerovSystem<Coord>& s = static_cast<NumerovSystem<Coord>&>(this->problem.s);
        Coord factor = 1 / (12 - (s.E - s.c(this->t + this->timestep)) * std::pow(this->timestep, 2));
        Coord tmp = (10 * std::pow(this->timestep, 2) * (s.E - s.c(this->t)) + 24) * this->yi +
                ((s.E - s.c(this->t - this->timestep)) * std::pow(this->timestep, 2) - 12) * yim1;
        state = factor * tmp;
        this->yim1 = this->yi;
        this->yi = state;
        this->t += this->timestep;
        return state;
    }
    private:
    Coord state;
    Coord yim1;
};


#endif