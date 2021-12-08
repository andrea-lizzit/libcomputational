#ifndef COPERATOR
#define COPERATOR

#include <functional>

template <template<class> typename Fin=std::function, template<class> typename Fout=std::function>
class COperator {
    virtual Fout operator() (Fin o);
};

std::function<int(int)> f = [](int a)->int {return 2*a;};

class OperatorH {
    OperatorH(std::function<int(int)> f) : tf(f) {};
    std::function<int(int)> tf;
    std::function<int(int)> operator()(std::function<int(int)> fi)
    {
        std::function<int(int)> ret = [](int x)->int {return };
    }

}


template<class F>
struct coperator_t {
    F f;
    template<class... Args>
    constexpr auto operator()(Args&&... args) &
        noexcept(noexcept(!std::invoke(f, std::forward<Args>(args)...)))
        -> decltype(!std::invoke(f, std::forward<Args>(args)...))
    {
        return !std::invoke(f, std::forward<Args>(args)...);
    }
};

template<class F>
constexpr coperator_t<std::decay_t<F>> COperator(F&& f)
{
    return { std::forward<F>(f) };
}

#endif