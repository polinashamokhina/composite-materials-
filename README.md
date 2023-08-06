# composite-materials-
Several possible stacks for a carbon/epoxy composite material, in order to apply them to an aircraft fuselage
Preliminary Design of Fuselage Skin for Commercial Airliner

Introduction
The aim of this project is to perform a preliminary design of the fuselage skin for a commercial airliner. The fuselage is represented as a cylindrical tube with a diameter of 5.75 meters, taking into account the effects of the tail and nose. The fuselage skin is treated locally like a plate under biaxial tension with axial stress (σ1) and circumferential stress (σ2).

General Assumptions
The fuselage skin is represented by a closed cylindrical tube of 5.75m diameter under internal pressure.
The ends of the cylinder are considered flat with the same stresses as the cylinder.
The fuselage skin is locally treated like a plate under biaxial tension with σaxial = σ1 and σcircumferential = σ2.
Calculation of Stiffness and Compliance Matrices
We calculate the stiffness and compliance matrices for a Carbon/Epoxy ply in material coordinates. These matrices are further reduced under plane stress hypothesis for a Transverse isentropic material.

Classical Laminate Theory
Using the Classical Laminate Theory, we calculate:
Membrane strain
Curvatures
Bending stiffness matrices (A, B, and D)

Study with [45/0/-45/90]s Laminate
In this section, we analyze different laminate stacking sequences for a carbon/epoxy composite material and determine the damage sequence until final failure, along with the pressure at which it will occur.

From the manufacturing perspective, it is more suitable to use six plies instead of eight in the laminate structure, as the [45/0/90/-45]s configuration has fixed angles that cannot be changed without compromising the quasi-isotropic property. Additionally, using fewer plies reduces manufacturing costs, making it more advantageous for the company.

The optimal ply angle to achieve the maximum ply breaking pressure is found to be 42.6°, resulting in a first break pressure of 238.9 kPa.

Ultimately, the [90/43.5/-43.5]s laminate structure is identified as the strongest and least expensive option for the aircraft fuselage.

Requirements
The code requires MATLAB to be installed on your computer or using online version.

Usage
To run the MATLAB code, simply open the script and execute it in your MATLAB environment. The code is organized into sections, allowing you to run each section independently for specific calculations.

Notes
The code is provided for educational purposes and may require further refinement or customization for specific applications.
The information provided in this README is based on the analysis and results of the MATLAB code.

Author
The MATLAB code and README are authored by Polina Shamokhina, p.shamohina@mail.ru

Feedback
If you have any feedback, questions, or suggestions, feel free to contact me at p.shamohina@mail.ru
