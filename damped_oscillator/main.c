#include "oscillator.h"
#include <stdio.h>
#include <stdlib.h> 

int main() {
    double dt = 0.01;
    double t_max = 10.0;
    double t = 0.0;
    State s;
    FILE *fp;   
    Parameters p;
    p.mass = 1.0,p.damping = 0.1, p.stiffness = 1.0;
    s.position = 1.0, s.velocity = 0.0;



    fp = fopen("oscillator_data.csv", "w");
    if (fp == NULL) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }
    for (t = 0; t <= t_max; t += dt) {
        fprintf(fp, "%f,%f,%f\n", t, s.position, s.velocity);
        s = euler_step(s, p, dt);
    }

    fclose(fp);

    return EXIT_SUCCESS;
}
