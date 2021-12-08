#include <COperator.h>

template <typename Tvector = float>
class CFunction : COperator<Tvector, Tvector>{
    virtual Tres operator() (Tres o);
};