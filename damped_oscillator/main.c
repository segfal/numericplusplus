#include "oscillator.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
  double dt = 0.01;
  double t_max = 10.0;
  double t = 0.0;

  FILE *fp;
  State s;
  Parameters p;
  State states[5] = {
      {1.0, 0.0}, {2.0, 0.5}, {3.0, 1.0}, {4.0, 1.5}, {5.0, 2.0}};
  Parameters params[5] = {{1.0, 0.1, 1.0},
                          {1.0, 0.2, 1.0},
                          {1.0, 0.3, 1.0},
                          {1.0, 0.4, 1.0},
                          {1.0, 0.5, 1.0}};

  for (int i = 0; i < 5; i++) {
    char fp_string[64];
    snprintf(fp_string, sizeof(fp_string), "oscillator_data_%d.csv", i);
    s = states[i];
    p = params[i];

    fp = fopen(fp_string, "w");
    if (fp == NULL) {
      perror("Error opening file");
      return EXIT_FAILURE;
    }
    for (t = 0; t <= t_max; t += dt) {
      fprintf(fp, "%f,%f,%f\n", t, s.position, s.velocity);
      s = euler_step(s, p, dt);
    }

    fclose(fp);
  }

  return EXIT_SUCCESS;
}
