# CarND-Unscented-Kalman-Filter-Project

> Unscented Kalman Filter Project for Self-Driving Car ND

[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

![UKF results](https://user-images.githubusercontent.com/4352286/37217276-1f18d5f0-238b-11e8-9d1a-d686384e28e7.png)

The goal of this project is to use an Unscented Kalman Filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. The system will output the estimated x, y position from the Kalman filter's state vector, along with the RMSE, to the simulator provided by Udacity. Passing the project requires obtaining RMSE values that are lower than the tolerance outlined in the project rubric.

## Overview
Starting to work on this project consists of the following steps:


1. Install uWebSocketIO and all the required [dependencies](#installation-and-dependencies)
2. Clone this repository
3. Build the main program 
    - `mkdir build`
    - `cd build`
    - `cmake ..`
    - `make`
4. Launch `./UnscentedKF`
5. Launch the Udacity Term 2 simulator
6. Enjoy!

---

## Installation and Dependencies

This project involves the Udacity Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

This repository includes two files that can be used to set up and install [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems. For windows you can use either Docker, VMware, or even [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) to install uWebSocketIO. Please see [this concept in the classroom](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/16cf4a78-4fc7-49e1-8621-3450ca938b77) for the required version and installation scripts.

### Other Important Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
  
Once all the dependencies has been installed **clone** the project:

```sh
git clone https://github.com/gdangelo/CarND-Unscented-Kalman-Filter-Project/
```

and follow the steps 3 to 6 of the [Overview section](#overview) in order to build and run the main program.
  
---

## Questions or Feedback

> Contact me anytime for anything about my projects or machine learning in general. I'd be happy to help you :wink:

* Twitter: [@gdangel0](https://twitter.com/gdangel0)
* Linkedin: [Grégory D'Angelo](https://www.linkedin.com/in/gregorydangelo)
* Email: [gregory@gdangelo.fr](mailto:gregory@gdangelo.fr)



