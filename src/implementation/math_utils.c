
float add_float(float a, float b) { return a + b; }
float sub_float(float a, float b) { return a - b; }
float mul_float(float a, float b) { return a * b; }
float div_float(float a, float b) { return a / b; }
float mod_float(float a, float b) { return a - (int)(a / b) * b; }

float power_float(float a, float n) {
  float res = 1.0f;
  int exp = (int)n;
  for (int i = 0; i < exp; i++) {
    res *= a;
  }
  return res;
}

float root_float(float a, float n) {
  float x = a / 2.0;
  for (int i = 0; i < 20; i++) {
    x = (1.0 / n) * ((n - 1) * x + a / power_float(x, n - 1));
  }
  return x;
}

int add_int(int a, int b) { return a + b; }
int sub_int(int a, int b) { return a - b; }
int mul_int(int a, int b) { return a * b; }
int div_int(int a, int b) { return a / b; }

int mod_int(int a, int b) { return a % b; }

int power_int(int a, int n) {
  int res = 1;
  for (int i = 0; i < n; i++) {
    res *= a;
  }
  return res;
}

float fact(float a) {
  float res = 1.0f;
  for (int i = 1; i <= (int)a; i++) {
    res *= i;
  }
  return res;
}
