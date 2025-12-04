#include "../include/math_utils.h"
#include <stdio.h>

int main() {
  int i1 = 10, i2 = 3;
  float f1 = 10.5f, f2 = 3.2f;

  printf("--- Integer Operations ---\n");
  printf("Add: %d + %d = %d\n", i1, i2, add(i1, i2));
  printf("Sub: %d - %d = %d\n", i1, i2, sub(i1, i2));
  printf("Mul: %d * %d = %d\n", i1, i2, mul(i1, i2));
  printf("Div: %d / %d = %d\n", i1, i2, divide(i1, i2));

  printf("\n--- Float Operations ---\n");
  printf("Add: %.2f + %.2f = %.2f\n", f1, f2, add(f1, f2));
  printf("Sub: %.2f - %.2f = %.2f\n", f1, f2, sub(f1, f2));
  printf("Mul: %.2f * %.2f = %.2f\n", f1, f2, mul(f1, f2));
  printf("Div: %.2f / %.2f = %.2f\n", f1, f2, divide(f1, f2));

  return 0;
}
