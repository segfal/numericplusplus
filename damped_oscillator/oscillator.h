// damped_oscillator/oscillator.h
#ifndef OSCILLATOR_H
#define OSCILLATOR_H

#include <stdio.h>

typedef struct {
  double position;
  double velocity;
} State;

typedef struct {
  double mass;
  double damping;
  double stiffness;
} Parameters;

// Function Prototypes for numerical solving
State calculate_derivatives(State current, Parameters params);
State euler_step(State current, Parameters params, double dt);

#endif // OSCILLATOR_H
