# Writeup for P3 of Self-Driving Car Nanodegree Term2: Kidnapped Vehicle
[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

### Matt, MinGyu, Kim
---

[image1]: ./figures/particle_filter_flow.PNG "PARTICLE_FILTER_PROCESS"
[image2]: ./figures/tracking.png "SCREENSHOT"

### Brief Description

In this project, I built the localization program using [particle filter](https://en.wikipedia.org/wiki/Particle_filter). The following figure is the processing flow of a particle filter.

![alt_text][image1]

Each steps are implemented in separate functions in `particle_filter.cpp`

---

### Simulation

The compiled program which is in `build` received the observations of landmarks (the positions of landmarks in the vehicle's coordinate system) and processed them to update the weights of particles. Finally, the program returns the estimated location of the vehicles which is presented as blue circle by the simulator.


The program was estimating the position of the vehicle well during the entire simulation as you can see in the following figure.

![alt_text][image1]
