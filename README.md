# TISE Simulator (1D Version)

## Overview

This project is a **1D TISE (Time-Independent Schrödinger Equation) Simulator** built using Python, **Tkinter**, and **Matplotlib**. It allows you to simulate quantum mechanical wavefunctions in a 1-dimensional potential well. The program computes and visualizes the eigenstates and corresponding energy levels of a quantum particle within different types of potential wells (Infinite Well, Finite Well, Harmonic Oscillator), and provides options to visualize the probability density of the eigenfunctions.

## Features

* **1D Quantum Well Simulation**: Simulates and visualizes the wavefunctions of a quantum particle in 1D.
* **Customizable Potentials**: Users can choose from different types of potentials (Infinite Well, Finite Well, Harmonic Oscillator) or define custom potentials using Python expressions.
* **Eigenstate Computation**: Computes energy eigenvalues and eigenfunctions for the selected potential.
* **Visualizations**: Displays the probability densities for the selected eigenstates.

## Requirements

* Python 3.x
* Libraries:

  * `numpy`
  * `matplotlib`
  * `tkinter` (for GUI)

You can install the required libraries using **pip**:

```bash
pip install numpy matplotlib
```

`tkinter` is typically included in Python installations. If it's not available, you can install it using:

```bash
sudo apt-get install python3-tk
```

## Installation

1. Clone or download the project:

   ```bash
   git clone <repository-url>
   ```

2. Install the required dependencies using pip:

   ```bash
   pip install -r requirements.txt
   ```

3. Run the Python script:

   ```bash
   python tise_simulator.py
   ```

## Usage

### GUI Interface:

1. **Mass (m)**: Set the mass of the particle.

2. **ħ (hbar)**: Set the value of Planck's constant divided by 2π.

3. **Length (L)**: Set the length of the potential well.

4. **Grid points (N)**: Set the number of points used in the grid for numerical calculations.

5. **Potential Type**:

   * **Infinite Well**: A simple square potential well with infinitely high walls at the boundaries.
   * **Finite Well**: A finite square potential well, where you can specify the depth of the well.
   * **Harmonic Oscillator**: A potential of the form (V(x) = \frac{1}{2} m \omega^2 x^2), representing a harmonic oscillator.
   * **Custom**: Allow users to input custom Python expressions to define the potential function.

6. **Number of Eigenstates**: Choose how many eigenstates to compute and display.

7. **Select State**: Choose which eigenstate to visualize.

8. **Show State**: Click to visualize the selected eigenstate's probability density.

9. **Energy Eigenvalues**: Displays a list of computed energy eigenvalues.

### Visualization:

* The program will plot the **probability density** (|ψ(x)|²) of the selected eigenstate.
* For an infinite or finite well, this will show a graph of the standing wave.
* For a harmonic oscillator, the plot will show the Gaussian-shaped wavefunctions.
* You can save the plot as a PNG file for further analysis or documentation.

## Example Usage

1. **Setting Parameters**:

   * Set the **mass** of the particle (e.g., 1.0 kg).
   * Set **ħ** (Planck's constant over 2π, default is 1.0).
   * Set the **length** of the potential well (e.g., 10.0 m).
   * Set **grid points** (e.g., 800).
   * Choose the **potential type**: e.g., Infinite Well, Finite Well, or Harmonic Oscillator.
2. **Computing Eigenstates**:

   * Click **Compute** to calculate the eigenstates and eigenvalues for the selected potential.
3. **Displaying Results**:

   * Choose the **state** you want to visualize (e.g., first eigenstate).
   * Click **Show State** to plot the probability density of that eigenstate.
4. **Saving the Plot**:

   * Click **Save plot as PNG** to save the current plot as an image file.

## License

This project is licensed under the MIT License.

## Future Enhancements

* Support for **higher-dimensional simulations** (2D and 3D).
* Ability to export results to files (e.g., CSV, JSON).
* Implement **time-dependent simulations** to animate wavefunction evolution.

---

