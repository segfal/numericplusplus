#pragma once
#include <iostream>
#include <string>
#include <vector>

class BigInt {
private:
  std::vector<int> digits;
  bool isNegative;

  void removeLeadingZeros();

public:
  BigInt();
  BigInt(long long n);
  BigInt(std::string s);

  BigInt operator+(const BigInt &other) const;
  BigInt operator-(const BigInt &other) const;
  BigInt operator*(const BigInt &other) const;

  bool operator==(const BigInt &other) const;
  bool operator!=(const BigInt &other) const;
  bool operator<(const BigInt &other) const;
  bool operator>(const BigInt &other) const;
  bool operator<=(const BigInt &other) const;
  bool operator>=(const BigInt &other) const;

  friend std::ostream &operator<<(std::ostream &os, const BigInt &bi);
};
