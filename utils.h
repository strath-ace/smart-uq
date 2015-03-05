#ifndef UTILS_H
#define UTILS_H

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

#endif // UTILS_H
