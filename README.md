# Channel Splicing for Wireless Sensing

This project explores the development and evaluation of channel splicing methods for wireless sensing applications. Channel splicing combines Channel State Information (CSI) from multiple frequency bands to create a comprehensive view of the channel response over a larger bandwidth. This approach enhances resolution and accuracy, which is essential for precise applications like indoor localization.

---

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Methodology](#methodology)
- [Algorithms Implemented](#algorithms-implemented)
- [File Structure](#file-structure)
- [Installation and Usage](#installation-and-usage)
- [References](#references)
- [License](#license)

---

## Introduction

Channel splicing is a technique that combines CSI measurements from multiple frequency bands, resulting in a high-resolution composite channel response. This method is particularly useful for applications like:
- High-precision indoor localization
- Environmental monitoring
- Motion detection and gesture recognition

The project was conducted at the **Telecommunications Circuits Laboratory, EPFL**, under the supervision of **Prof. Andreas Peter Burg**.

---

## Features

- **Channel Splicing Technique:** Combines CSI from multiple carriers to enhance sensing accuracy.
- **Multiple Reconstruction Algorithms:** Implemented and compared the following methods:
  - Chronos (Non-uniform Discrete Fourier Transform with sparsity constraint)
  - MUSIC (Multiple Signal Classification)
  - MVDR (Minimum Variance Distortionless Response)
- **Simulated and Real Data Analysis:** Tested with simulated Rayleigh fading data and real-world measurements.

---

## Methodology

1. **Data Generation:**
   - Simulated data using Rayleigh fading models.
   - Real data collected through controlled experiments with two distinct signal paths.

2. **Channel Splicing:**
   - Combines CSI data from different frequency bands to form a spliced CSI.
   - Computes Channel Impulse Response (CIR) using Inverse Discrete Fourier Transform (IDFT).

3. **Reconstruction Algorithms:**
   - **Chronos:** Utilizes Non-uniform Discrete Fourier Transform (NDFT) with sparsity constraints.
   - **MUSIC:** Separates signal and noise subspaces using eigenvalue decomposition.
   - **MVDR:** Enhances the desired signal path by minimizing interference.

4. **Performance Evaluation:**
   - Compared the accuracy and resolution of Chronos, MUSIC, and MVDR.
   - Combined MUSIC and MVDR for improved peak detection and distance accuracy.

---

## Algorithms Implemented

### 1. Chronos
- Uses Non-uniform Discrete Fourier Transform (NDFT) with sparsity constraints.
- Optimizes the solution using a gradient descent approach with sparsification steps.

### 2. MUSIC (Multiple Signal Classification)
- Identifies signal and noise subspaces using eigenvalue decomposition.
- Computes a pseudospectrum to detect time-of-arrival by identifying peaks.

### 3. MVDR (Minimum Variance Distortionless Response)
- Uses beamforming techniques to enhance the desired signal while minimizing interference.
- Computes a pseudospectrum for accurate peak detection.

---

## Installation and Usage

### Prerequisites
- MATLAB with Communication System Toolbox
- Clone this repository:
    ```sh
    git clone https://github.com/Portgas37/Channel-Splicing-For-Wireless-Sensing.git
    cd Channel-Splicing-For-Wireless-Sensing
    ```

### Running the Code
1. Open MATLAB and navigate to the project directory.
2. Run the desired script. For example, to test multi-band CSI splicing:
    ```matlab
    run('Scripts/scriptCsiMultiBand.m');
    ```
3. Modify parameters directly in the script files to test different scenarios.

---

## References
- [Bachelor Project Report](./Report/Bachelor_Project___Report.pdf)
- [Project Presentation](./Presentation/PrésentationVersionDIm.pdf)
- MATLAB Documentation:
  - [Rayleigh Channel](https://www.mathworks.com/help/comm/ref/comm.rayleighchannel-system-object.html)
  - [Rician Channel](https://www.mathworks.com/help/comm/ref/comm.ricianchannel-system-object.html)

---


## Acknowledgments
- **Prof. Andreas Peter Burg** - Project Supervisor, EPFL
- **Andreas T. Kristensen** - Project Mentor, EPFL
- **Telecommunications Circuits Laboratory (TCL), EPFL**

---

[Ecole Polytechnique Fédérale de Lausanne (EPFL)](https://www.epfl.ch)

