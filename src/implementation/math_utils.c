#include "../../headers/math_utils.h"

// Float implementations
float add_float(float a, float b) { return a + b; }
float sub_float(float a, float b) { return a - b; }
float mul_float(float a, float b) { return a * b; }
float div_float(float a, float b) { return a / b; }

float mod_float(float a, float b) { return a - (int)(a / b) * b; }

float power_float(float a, float n) {
  int exp = (int)n;
  if (exp == 0)
    return 1;
  float temp = power_float(a, exp / 2);
  if (exp % 2 == 0)
    return temp * temp;
  else
    return a * temp * temp;
}

float root_float(float a, float n) {
  float x = a / 2.0; // Initial guess
  for (int i = 0; i < 20; i++) {
    x = (1.0 / n) * ((n - 1) * x + a / power_float(x, n - 1));
  }
  return x;
}

// Int implementations
int add_int(int a, int b) { return a + b; }
int sub_int(int a, int b) { return a - b; }
int mul_int(int a, int b) { return a * b; }
int div_int(int a, int b) { return a / b; }

int mod_int(int a, int b) { return a % b; }

int power_int(int a, int n) {
  if (n == 0)
    return 1;
  int temp = power_int(a, n / 2);
  if (n % 2 == 0)
    return temp * temp;
  else
    return a * temp * temp;
}

float fact(float a) {
  if (a == 0)
    return 1;
  return a * fact(a - 1);
}
