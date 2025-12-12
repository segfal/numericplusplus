#include "oscillator.h"


State calculate_derivatives(State current, Parameters params) {
    double mass,damping,stiffness,position,velocity;
    
    mass = params.mass;
    damping = params.damping;
    stiffness = params.stiffness;

    position = current.velocity;
    velocity = -(stiffness/mass)*current.position - (damping/mass)*current.velocity;
    
    State derivatives;
    derivatives.position = position;
    derivatives.velocity = velocity;
    return derivatives;

}


State euler_step(State current, Parameters params, double dt) {
    State derrivatives = calculate_derivatives(current, params);
    State transitionState;
    transitionState.position = current.position + derrivatives.position * dt;
    transitionState.velocity = current.velocity + derrivatives.velocity * dt;
    return transitionState;
}
