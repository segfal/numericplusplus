#include "../include/numbers/BigInt.hpp"
#include <cassert>
#include <iostream>

int main() {
  BigInt a(123456789);
  BigInt b("987654321");
  BigInt c = a + b;

  std::cout << "a: " << a << "\n";
  std::cout << "b: " << b << "\n";
  std::cout << "a + b: " << c << "\n";

  BigInt d("10000000000000000000"); // 10^19
  BigInt e("1");
  std::cout << "d: " << d << "\n";
  std::cout << "d + e: " << d + e << "\n";
  std::cout << "d - e: " << d - e << "\n";

  // Subtraction tests
  BigInt sub1(100), sub2(50);
  std::cout << "100 - 50 = " << sub1 - sub2 << "\n";
  std::cout << "50 - 100 = " << sub2 - sub1 << "\n";

  // Multiplication tests
  BigInt mul1(12), mul2(12);
  std::cout << "12 * 12 = " << mul1 * mul2 << "\n";

  return 0;
}
