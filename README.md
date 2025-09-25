# MEJ-ESS: An Enhanced Monte Carlo Simulator for Electron Transport in Gaseous Detectors


MEJ-ESS is an open-source Monte Carlo simulator designed to compute electron transport parameters in arbitrary gas mixtures. It builds on the existing METHES code. The tool is intended for researchers in high-energy physics and plasma modeling, particularly in the context of emerging eco-friendly gas mixtures.

---

### **Key Features**

* **GNU Octave Compatibility**: The code has been refactored for compatibility with GNU Octave. This eliminates the need for a MATLAB license, which restricted the use of METHES.
* **Expanded Output Parameters**: New routines were developed to calculate additional transport parameters, such as Townsend ionization ($\alpha$) and attachment ($\eta$) coefficients, and longitudinal ($D_L$) and transverse ($D_T$) diffusion coefficients.
* **User-Configurable Penning Transfer**: MEJ-ESS introduces user-configurable Penning transfer modeling. This allows simulations to include secondary ionization processes critical in gas mixtures with noble gases and quenchers.
* **Parallel Processing**: MEJ-ESS supports parallel execution using GNU Octave's built-in `parcellfun` function. This allows simulations to efficiently utilize multiple CPU cores and significantly reduce runtime.
* **LXCat Integration**: The tool uses LXCat cross-section data, providing a flexible platform for simulating novel or eco-friendly gas mixtures.

---

### **Getting Started**

#### **Prerequisites**

To use MEJ-ESS, you will need to have the following installed:

* **GNU Octave**: Version 9.3.0 or later is recommended.
* **Octave's Parallel Package**: You can install it from within Octave by running:
    ```octave
    pkg install -forge parallel
    ```

#### **Installation**

1.  Clone this repository to your local machine:
    ```bash
    git clone https://github.com/Asser-MohamedR/MEJ-ESS.git
    ```
2.  Navigate to the repository directory:
    ```bash
    cd MEJ-ESS
    ```

### **How to Run a Simulation**

1.  Open GNU Octave.

2.  Add the necessary directories to your Octave path to ensure the program can find all the required functions and data files.

    ```octave
    addpath('src');
    addpath('_functions');
    addpath('_Xsection');
    ```

3.  Run the main Monte Carlo simulation script, which is located in the `_mixture_run` directory.

    ```octave
    run _mixture_run/MonteCarlo_run.m
    ```

### **Acknowledgments**

MEJ-ESS is a derivative work based on the **METHES** code. We gratefully acknowledge the authors of METHES, M. Rabie and C. M. Franck, for their foundational work.

This project was made possible with financial support from Sultan Qaboos University and the computational resources provided by the Center for High Energy Physics, Fayoum University (CHEP-FU).

---

---

### **License**

This project is licensed under the **GNU General Public License v3.0**. As a derivative of the METHES code, which is also GPL-licensed, MEJ-ESS must also be licensed under the GPL. See the [LICENSE](LICENSE) file for details.
