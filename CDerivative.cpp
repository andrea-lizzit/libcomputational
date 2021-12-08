#include <COperator.h>

template<class>
class CDerivative : COperator<std::function<class> {};

template<class>
std::function<class C> CDerivative<C>::operator() (std::function<class C> fun) {};

template<template<class> typename Fout>
Fout derivative()

// you need to rethink more clearly what needs to be achieved.
// For the time being, implement the easy and quick solution