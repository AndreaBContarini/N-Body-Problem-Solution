# N-Body Problem (and Satellite Landing Problem)

This repository contains the implementation of various astrophysical simulations in **C**, developed as part of a project for the Theoretical Astrophysics Course in the MSc Physics program. The project involves studying the dynamics of celestial bodies and satellites using numerical methods.

---

## Files and Their Purpose

### 1. `2body.c`
- **Description**: Solves the 2-body problem using the leapfrog algorithm. It predicts the motion of two celestial bodies interacting through gravitational forces.
- **Key Features**:
  - Implements gravitational acceleration computation.
  - Calculates kinetic and potential energy for analysis.
  - Outputs the motion trajectory to a file for visualization.
- **Numerical Method**: Leapfrog algorithm.

---

### 2. `3body.c`
- **Description**: Extends the 2-body problem to a 3-body system, using similar methods to study gravitational interactions among three bodies.
- **Key Features**:
  - Sets up initial conditions for three bodies positioned at the vertices of a Pythagorean triangle.
  - Analyzes motion trajectories with zero initial velocities.
  - Outputs results for 2D visualization.
- **Numerical Method**: Leapfrog algorithm.

---

### 3. `satellite.c`
- **Description**: Simulates the motion of a satellite falling toward Earth under the influence of gravity, atmospheric drag, and mass loss.
- **Key Features**:
  - Supports circular and radial motion simulations.
  - Solves equations using the Runge-Kutta 4th-order method.
  - Incorporates effects of atmospheric drag and satellite mass loss.
  - Outputs energy, radius, and trajectory for detailed analysis.
- **Applications**: Studies the impact of initial conditions and parameter changes on satellite motion.

---

### 4. `Nbody.c`
- **Description**: Generalizes the simulation to an N-body system, solving for the motion of multiple celestial bodies under mutual gravitational forces.
- **Key Features**:
  - Supports random initial position generation within a sphere.
  - Computes gravitational acceleration for N bodies.
  - Outputs positions and energies for analysis.
  - Uses the leapfrog algorithm for efficient computation.
- **Applications**: Analyzes large-scale gravitational systems, such as star clusters.

---

### 5. `Article_Astrophysics_Project.pdf`
- **Description**: A comprehensive report detailing the theoretical background, numerical methods, and results of the simulations.
- **Contents**:
  - Problem descriptions for satellite and N-body simulations.
  - Derivations of equations and their numerical solutions.
  - Detailed analysis and comparisons of theoretical vs. numerical results.

---

## Requirements

- **C Compiler**: GCC or any compatible compiler.
- **Libraries**: Standard C libraries (e.g., `math.h`, `stdio.h`, `stdlib.h`).

---

## Usage

1. Compile the desired file using:
   ```bash
   gcc -o <output> <file.c> -lm
---

## Authors
- Aurora Abbondanza
- Alessandro Agapito
- Andrea Belli Contarini
- Riccardo Caleno
