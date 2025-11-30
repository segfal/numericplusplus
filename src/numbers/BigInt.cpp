#include "../../include/numbers/BigInt.hpp"
#include <algorithm>
#include <iostream>

// Private Helpers
void BigInt::removeLeadingZeros() {
  while (digits.size() > 1 && digits.back() == 0) {
    digits.pop_back();
  }
  if (digits.size() == 1 && digits[0] == 0)
    isNegative = false;
}

// Constructors
BigInt::BigInt() {
  digits.push_back(0);
  isNegative = false;
}

BigInt::BigInt(long long n) {
  if (n < 0) {
    isNegative = true;
    n = -n;
  } else {
    isNegative = false;
  }
  if (n == 0)
    digits.push_back(0);
  while (n > 0) {
    digits.push_back(n % 10);
    n /= 10;
  }
}

BigInt::BigInt(std::string s) {
  if (s.length() == 0) {
    digits.push_back(0);
    isNegative = false;
    return;
  }

  isNegative = (s[0] == '-');
  int start = isNegative ? 1 : 0;

  for (int i = s.length() - 1; i >= start; i--) {
    if (!isdigit(s[i]))
      throw std::invalid_argument("Invalid character in BigInt string");
    digits.push_back(s[i] - '0');
  }
  removeLeadingZeros();
}

// Comparison
bool BigInt::operator==(const BigInt &other) const {
  return isNegative == other.isNegative && digits == other.digits;
}

bool BigInt::operator!=(const BigInt &other) const { return !(*this == other); }

bool BigInt::operator<(const BigInt &other) const {
  if (isNegative != other.isNegative)
    return isNegative;
  if (digits.size() != other.digits.size()) {
    return isNegative ? digits.size() > other.digits.size()
                      : digits.size() < other.digits.size();
  }
  for (int i = digits.size() - 1; i >= 0; i--) {
    if (digits[i] != other.digits[i]) {
      return isNegative ? digits[i] > other.digits[i]
                        : digits[i] < other.digits[i];
    }
  }
  return false;
}

bool BigInt::operator>(const BigInt &other) const { return other < *this; }
bool BigInt::operator<=(const BigInt &other) const { return !(*this > other); }
bool BigInt::operator>=(const BigInt &other) const { return !(*this < other); }

// Arithmetic
BigInt BigInt::operator+(const BigInt &other) const {
  if (isNegative == other.isNegative) {
    BigInt res;
    res.isNegative = isNegative;
    res.digits.clear();
    int carry = 0;
    int n1 = digits.size(), n2 = other.digits.size();
    for (int i = 0; i < std::max(n1, n2) || carry; ++i) {
      int sum =
          carry + (i < n1 ? digits[i] : 0) + (i < n2 ? other.digits[i] : 0);
      res.digits.push_back(sum % 10);
      carry = sum / 10;
    }
    return res;
  }
  // Handle signs: a + (-b) -> a - b
  if (isNegative) {
    BigInt temp = *this;
    temp.isNegative = false;
    return other - temp;
  }
  BigInt temp = other;
  temp.isNegative = false;
  return *this - temp;
}

BigInt BigInt::operator-(const BigInt &other) const {
  if (isNegative != other.isNegative) {
    // a - (-b) -> a + b
    BigInt temp = other;
    temp.isNegative = !other.isNegative;
    return *this + temp;
  }

  // Compare magnitudes
  bool thisIsSmaller = false;
  if (digits.size() != other.digits.size()) {
    thisIsSmaller = digits.size() < other.digits.size();
  } else {
    for (int i = digits.size() - 1; i >= 0; i--) {
      if (digits[i] != other.digits[i]) {
        thisIsSmaller = digits[i] < other.digits[i];
        break;
      }
    }
  }

  if (thisIsSmaller) {
    BigInt res = other - *this;
    res.isNegative = !isNegative;
    return res;
  }

  // Assume |this| >= |other|
  BigInt res;
  res.isNegative = isNegative;
  res.digits.clear();
  int borrow = 0;
  int n1 = digits.size(), n2 = other.digits.size();

  for (int i = 0; i < n1; ++i) {
    int sub = digits[i] - borrow - (i < n2 ? other.digits[i] : 0);
    if (sub < 0) {
      sub += 10;
      borrow = 1;
    } else {
      borrow = 0;
    }
    res.digits.push_back(sub);
  }
  res.removeLeadingZeros();
  return res;
}

// Naive Multiplication
BigInt BigInt::operator*(const BigInt &other) const {
  BigInt res;
  res.digits.resize(digits.size() + other.digits.size(), 0);
  res.isNegative = isNegative != other.isNegative;

  for (size_t i = 0; i < digits.size(); ++i) {
    int carry = 0;
    for (size_t j = 0; j < other.digits.size() || carry; ++j) {
      long long cur =
          res.digits[i + j] +
          digits[i] * 1LL * (j < other.digits.size() ? other.digits[j] : 0) +
          carry;
      res.digits[i + j] = cur % 10;
      carry = cur / 10;
    }
  }
  res.removeLeadingZeros();
  return res;
}

// Output
std::ostream &operator<<(std::ostream &os, const BigInt &bi) {
  if (bi.digits.empty()) {
    os << "0";
    return os;
  }
  if (bi.isNegative)
    os << "-";
  for (int i = bi.digits.size() - 1; i >= 0; i--) {
    os << bi.digits[i];
  }
  return os;
}
