# cem_peec

## Computational Electromagnetics – PEEC Method for Rectangular Conductors

This repository contains MATLAB codes developed as part of a Computational Electromagnetics project. The project applies the PEEC (Partial Element Equivalent Circuit) method to simulate electromagnetic field behavior in conductors. Using the PEEC formulation, the simulation allows for either terminal current or terminal voltage stimulation at the conductor extremities.

## What is PEEC (Partial Element Equivalent Circuit)?

PEEC is a numerical technique that discretizes complex conductors into smaller filamentary elements, each modeled by an equivalent circuit. This method enables accurate simulation of current diffusion, the skin effect, and electromagnetic coupling in conductors.

## Repository Overview

```
cem_peec/
├── method1/                 # Green's Function Approach Implementation
├── method2/                 # Energy Storage Approach Implementation
├── comparison/              # Comparative Analysis & Validation
└── presentation.pdf         # Documentation & Presentations
```
I recommend starting with the `presentation.pdf` file to better understand the different scripts.

### Core Components

#### Method 1: Green's Function Approach
```
method1/
├── CircularConductor/                    # Circular geometry implementation
│   ├── generateCircularConductorMesh.m   # Non-uniform mesh generation
│   └── SingleConductorCircularShape.m    # Single conductor simulation
├── RectangularConductor/                 # Rectangular geometry implementation
│   ├── generateNonUniformMesh.m          # Non-uniform mesh generation
│   ├── SingleConductor.m                 # Main simulation code
│   └── TestConvergenceLosses.m           # Convergence analysis
└── COMSOL/                               # Reference solutions
    └── SingleConductor.mph               # COMSOL model file
```

#### Method 2: Energy Storage Approach
```
method2/
├── single_conductor/                           # Single conductor analysis
│   ├── SingleConductor_InductanceComparison.m  # Inductance validation
│   └── SingleConductor_ReferenceLosses.xlsx    # Reference data
└── multi_conductor/                            # Multiple conductor analysis
    └── MultiConductor_UniformMesh.m            # Multi-conductor simulation
```

#### Validation & Comparison
```
comparison/
├── C_1C_Original_Dim_Current_Distribution_Fixed_Mesh/   # Current stimulation
│   ├── TestConvergenceLosses_CURRENT.m                  # Analysis code
│   └── Images_Method1_and_Method2/                      # Results
└── V_1C_Original_Dim_Current_Distribution_Fixed_Mesh/   # Voltage stimulation
    ├── TestConvergenceLosses_VOLTAGE.m                  # Analysis code
    └── Images_Method1_and_Method2/                      # Results
```

## Detailed Description

- **method1/**  
  Contains MATLAB codes implementing Method 1, based on the Green's Function Approach as described in [1].  
  - **CircularConductor/**  
    Includes scripts for simulating circular conductor geometries. The function `generateCircularConductorMesh` uses a cosine-based, non-uniform meshing approach to discretize the circular cross-section into rectangular sub-elements. This yields a finer mesh near the curved boundary and a coarser one near the center, ensuring an accurate representation of the conductor's geometry. The mesh data (element centers, areas, and vertices) are used to compute the current distribution under a specified terminal current, and the total losses are compared with COMSOL reference data.
  
  - **RectangularConductor/**  
    Contains scripts for simulating rectangular conductor geometries. These scripts generate a non-uniform mesh over the conductor's cross-section using a sinusoidal distribution, which provides finer resolution near the edges and coarser resolution near the center. The resulting mesh is then used to calculate the resistance and inductance matrices (with mutual inductance computed via a Green's function and self-inductance derived from the local geometry). The combined impedance matrix allows for the computation of filament currents under terminal excitation, and the total losses are evaluated against COMSOL reference data. Convergence tests confirm that results improve with mesh refinement.

  - **COMSOL/**  
    Contains COMSOL Multiphysics project files used to generate reference solutions for losses validation and comparison purposes.

- **method2/**  
  Contains MATLAB codes implementing Method 2, based on the energy storage approach detailed in [2] and [3]. 
  - **single_conductor/**  
    Contains the script `SingleConductor_InductanceComparison.m` and the reference data file `SingleConductor_ReferenceLosses.xlsx`. This folder focuses on single-conductor simulations under terminal voltage stimulation. It constructs the partial inductance and resistance matrices, computes the filament currents, and evaluates the total losses, which are then compared against reference data.
  
  - **multi_conductor/**  
    Extends the energy storage approach to multi-conductor configurations. The implementation handles conductor interactions and proximity effects through mutual inductance calculations. Results are validated against external reference data.

- **comparison/**  
  Contains validation scripts comparing both methods under different conditions. The analysis includes a systematic study of both methods' performance across:
  - Different mesh densities to evaluate convergence properties
  - A wide frequency range (from quasi DC to high frequencies) to assess the methods' ability to capture skin and proximity effects
  
  - **Current Stimulation/** (`C_1C_Original_Dim_Current_Distribution_Fixed_Mesh/`)  
    Analyzes mesh convergence and accuracy under terminal current excitation. Multiple mesh configurations are tested to evaluate how the discretization affects the solution accuracy. The frequency sweep helps identify the operational limits of each method, particularly in capturing high-frequency phenomena like the skin effect.
    
  - **Voltage Stimulation/** (`V_1C_Original_Dim_Current_Distribution_Fixed_Mesh/`)  
    Analyzes mesh convergence and accuracy under terminal voltage excitation. Both test cases include comprehensive result visualization and error analysis, stored in their respective `Images_Method1_and_Method2/` directories.

- **presentation.pdf**  
  Contains comprehensive documentation of the theoretical background and mathematical derivations of both methods.

## References

1. D. P. Morisco, S. Kurz, H. Rapp, and A. Möckel, “A hybrid modeling approach for current diffusion in rectangular conductors,” *IEEE Transactions on Magnetics*, vol. 55, no. 9, pp. 1–9, 2019.
2. C. R. Paul, *Inductance: Loop and Partial*, Hoboken, NJ: John Wiley & Sons, 2010.
3. A. E. Ruehli, “Inductance calculations in a complex integrated circuit environment,” *IBM Journal of Research and Development*, vol. 16, no. 5, pp. 470–481, 1972.
