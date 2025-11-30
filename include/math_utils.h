#pragma once

#ifdef __cplusplus
extern "C" {
#endif

float add_float(float a, float b);
float sub_float(float a, float b);
float mul_float(float a, float b);
float div_float(float a, float b);
float mod_float(float a, float b);
float power_float(float a, float n);
float root_float(float a, float n);

int add_int(int a, int b);
int sub_int(int a, int b);
int mul_int(int a, int b);
int div_int(int a, int b);
int mod_int(int a, int b);
int power_int(int a, int n);

float fact(float a);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
inline int add(int a, int b) { return add_int(a, b); }
inline float add(float a, float b) { return add_float(a, b); }

inline int sub(int a, int b) { return sub_int(a, b); }
inline float sub(float a, float b) { return sub_float(a, b); }

inline int mul(int a, int b) { return mul_int(a, b); }
inline float mul(float a, float b) { return mul_float(a, b); }

inline int divide(int a, int b) { return div_int(a, b); }
inline float divide(float a, float b) { return div_float(a, b); }

inline int mod(int a, int b) { return mod_int(a, b); }
inline float mod(float a, float b) { return mod_float(a, b); }

inline int power(int a, int n) { return power_int(a, n); }
inline float power(float a, float n) { return power_float(a, n); }

inline float root(float a, float n) { return root_float(a, n); }

#else
#define add(a, b) _Generic((a), int: add_int, float: add_float)(a, b)
#define sub(a, b) _Generic((a), int: sub_int, float: sub_float)(a, b)
#define mul(a, b) _Generic((a), int: mul_int, float: mul_float)(a, b)
#define divide(a, b) _Generic((a), int: div_int, float: div_float)(a, b)
#define mod(a, b) _Generic((a), int: mod_int, float: mod_float)(a, b)
#define power(a, b) _Generic((a), int: power_int, float: power_float)(a, b)
#define root(a, b) root_float(a, b)
#endif
