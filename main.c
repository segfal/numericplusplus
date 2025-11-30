#include "headers/math_utils.h"
#include <stdio.h>

int main() {
  int i1 = 10, i2 = 3;
  float f1 = 10.5f, f2 = 3.2f;

  printf("--- Integer Operations ---\n");
  printf("Add: %d + %d = %d\n", i1, i2, add(i1, i2));
  printf("Sub: %d - %d = %d\n", i1, i2, sub(i1, i2));
  printf("Mul: %d * %d = %d\n", i1, i2, mul(i1, i2));
  printf("Div: %d / %d = %d\n", i1, i2, div(i1, i2));
  printf("Mod: %d %% %d = %d\n", i1, i2, mod(i1, i2));
  printf("Pow: %d ^ %d = %d\n", i1, i2, power(i1, i2));

  printf("\n--- Float Operations ---\n");
  printf("Add: %.2f + %.2f = %.2f\n", f1, f2, add(f1, f2));
  printf("Sub: %.2f - %.2f = %.2f\n", f1, f2, sub(f1, f2));
  printf("Mul: %.2f * %.2f = %.2f\n", f1, f2, mul(f1, f2));
  printf("Div: %.2f / %.2f = %.2f\n", f1, f2, div(f1, f2));
  printf("Mod: %.2f %% %.2f = %.2f\n", f1, f2, mod(f1, f2));
  printf("Pow: %.2f ^ 2 = %.2f\n", f1, power(f1, 2.0f));
  printf("Root: Sqrt(%.2f) = %.2f\n", f1, root(f1, 2.0f));

  printf("\n--- Factorial ---\n");
  printf("Fact(5.0) = %.2f\n", fact(5.0f));

  return 0;
}
