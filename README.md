# CAPE-OPEN Reactor Model
This is a chemical process modelling component developed as a Dynamic-Link Library (DLL) in C++ using CAPE-OPEN Binary Interop Architecture (COBIA) middleware. This makes it compatible with any CAPE-OPEN Standard compliant simulator. The unit simulates the behaviour of conventional plug flow reactor with a custom model and operation parameters dedicated to a single process using its built-in reaction package.

#### Video Demo:  https://youtu.be/YHaJ11XPH9Q# CAPE-OPEN Reactor Model

## Definitions

### Chemical Process
A chemical process is the process of changing one or more chemical compounds in a way that involves chemical reactions. Engineers use reaction kinetics and reactors' governing equations to predict, design and optimise such processes. The process is usually a part of a larger manufacturing proccess forming an industrial plant.
An equation-oriented model of the chemical process usually involves non-linear differential equations and require robust tools such as numerical algorithms and computing environments. It also requirea large amounts of data such as material constants and state-dependent properties. This has been made easier with Process Simulation software.

### [Process Simulation](https://en.wikipedia.org/wiki/Process_simulation)
> Process simulation is used for the design, development, analysis, and optimization of technical processes such as: chemical plants, chemical processes, environmental systems, power stations, complex manufacturing operations, biological processes, and similar technical functions.

### [Computer-Aided Process Engineering Open Standard (CAPE-OPEN)](https://www.colan.org/general-information-on-co-lan/)
> CAPE-OPEN consists of a series of specifications to expand the range of application of process simulation technologies. The CAPE-OPEN specifications specify a set of software interfaces that allow plug and play inter-operability between a given process modelling environment (PME) and a third-party process modelling component (PMC).
> CAPE-OPEN is an EU funded project supported by the non-profit organization [CO-LaN](https://www.colan.org/).

### [CAPE-OPEN Binary Interop Architecture (COBIA)](https://www.colan.org/experiences-projects/cape-open-binary-interop-architecture-cobia/)
> A new middleware, the CAPE-OPEN Binary Interop Architecture (COBIA), is the next step in the evolution of CAPE-OPEN. COBIA includes registration components, binary interoperability standards, and middleware that acts as a bridge between software components. Development of COBIA involves a number of tasks, grouped in phases, which are performed incrementally.
> COBIA serves as a propritary replacement to Microsoft's Component Object Model (COM), upon which all earlier developments have relied.


## Dependencies
1. [COBIA-Development SDK](https://colan.repositoryhosting.com/trac/colan_cobia/downloads) v1.2.0.8 
2. [Boost](https://www.boost.org/users/history/version_1_77_0.html) v1.77
3. [Windows 10 SDK](https://developer.microsoft.com/en-us/windows/downloads/sdk-archive/) v10.0.19041.0

## Reactor Model
The reactor is a simplified plug flow (tubular) reactor.
### Assumptions
1. 1D (No radial gradient)
2. No axial dispersion
4. Homogenous
5. Empty bed
6. Adiabatic
7. Constant heat capacity
8. Constant heat of Reaction
9. No pressure drop
### ODE System
<img width="173" alt="image" src="https://user-images.githubusercontent.com/80135041/145170402-73cb978d-22f1-413f-84b6-5e90b3f2fa3c.png">

### Reaction
<img width="505" alt="image" src="https://user-images.githubusercontent.com/80135041/145170505-10ceb6c5-3620-481d-ae53-7f1c9bb99d39.png">

## Implementation
Within a CAPE-OPEN compliant Process Modelling Environemnt (PME), the reactor acts as a black box. It calls the PME via a subset of the avialble standards according to its needs and answers the PME's calls as defined in the specifications.
### Input
1. feed (material stream)
2. reaactorLength (floating point parameter)
3. reactorVolume (floating point parameter)
### Output
1. product (material stream)
2. energy (energy stream)
3. conversion (floating point parameter)

## Process Simulation Using CAPE-OPEN Reactor Model
1. Simulation Executive: COFE 3.5.0.1.4 x64
2. Property Package Manager: TEA 3.5.0.2 x64
3. Property Package: ChemSep PCD 
4. Model Set: Equation of State (Peng-Robinson)
5. Reaction Package Manager: None
6. Reaction Package: built-in
7. feed data:
   - Temperature = 880K
   - Pressure    = 1.378 bar
   - flowrate    = 152.2 gmol/s
   - Composition = 100% C8H10
8. Reactor input: 
   - Length      = 6 m
   - Volume      = 0.77 m3

<img width="347" alt="image" src="https://user-images.githubusercontent.com/80135041/147827559-ad6924e3-a7b3-4013-9cf0-3487c5567769.png">

## Benchmarking Against Commercial Simulator (ASPEN HYSYS v9.0)
Values are very close to results from commerical simulator with a much more complex models.

<img width="674" alt="HYSYS_PFR_Adiabatic" src="https://user-images.githubusercontent.com/80135041/145166322-36f82c31-f9ed-4963-acda-7f8e3fa74da7.png">
