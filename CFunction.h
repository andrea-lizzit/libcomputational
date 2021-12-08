#ifndef CFUNCTION
#define CFUNCTION

#include <COperator.h>


template <typename Tvector = float>
class CFunction : COperator<Tvector, Tvector>{
    virtual Tres operator() (Tres o);
};

// template <typename Tx, typename Ty>
// class Function {
//     virtual Ty fun(Tx);
//     Ty* operator();
// };

// template <typename Tx, typename Ty>
// Ty Function::fun(Tx x)
// {

// }
#endif