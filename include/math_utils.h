#pragma once

#ifdef __cplusplus
extern "C" {
#endif

float add_float(float a, float b);
float sub_float(float a, float b);
float mul_float(float a, float b);
float div_float(float a, float b);

int add_int(int a, int b);
int sub_int(int a, int b);
int mul_int(int a, int b);
int div_int(int a, int b);

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

#else
#define add(a, b) _Generic((a), int: add_int, float: add_float)(a, b)
#define sub(a, b) _Generic((a), int: sub_int, float: sub_float)(a, b)
#define mul(a, b) _Generic((a), int: mul_int, float: mul_float)(a, b)
#define divide(a, b) _Generic((a), int: div_int, float: div_float)(a, b)
#endif
