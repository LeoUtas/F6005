#include <TMB.hpp>
// #include <math.h>

template <class Type>
Type power(Type x, int y) // x is defined as Type for later use of unspecified x (line 34)
{
  Type result = 1;
  for (int i = 0; i < y; i++)
  {
    result = result * x;
  }
  return result;
}

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(ssb);
  DATA_VECTOR(rec);
  DATA_VECTOR(lo