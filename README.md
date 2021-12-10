# COBIA Unit Operation
CAPE-OPEN Unit Operation developed for COBIA middleware. The unit emulates the behaviour of conventional unit operation. The goal is to ultimately develop a non-conventional unit once the development of all modules is matured.

## Unit version
0.4 (PFR with integrated reaction package)

## Dependencies
1. COBIA SDK 1.2.0.8
2. Boost 1.77

## PFR Model
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
### Streams
1. feed
2. product
3. energy
### Unit Parameters
1. reaactorLength
2. reactorVolume
### Energy Parameters
1. work
2. temperatureLow  // Not utilised
3. temperatureHigh // Not utilised

## COCO simulation
1. Simulation Executive: COFE 3.5.0.1.4 x64
2. Property Package Manager: TEA 3.5.0.2 x64
3. Property Package: ChemSep PCD 
4. Model Set: Equation of State (Peng-Robinson)
5. Reaction Package Manager: None
6. Reaction Package: hard-coded
7. Inlet data:
   - Temperature = 880K
   - Pressure    = 1.378 bar
   - flowrate    = 152.2 gmol/s
   - Composition = 100% C8H10
8. Reactor data: 
   - Length      = 6 m
   - Volume      = 0.77 m3

<img width="479" alt="COBIA_PFR_Adiabatic" src="https://user-images.githubusercontent.com/80135041/145171803-a2a03a67-6bea-42db-8299-a779e9c6cfd1.png">

## HYSYS Benchmark
Values vary around 2% mainly due to model assuming constant heat of reaction and heat capacity throughout 210°C range

<img width="674" alt="HYSYS_PFR_Adiabatic" src="https://user-images.githubusercontent.com/80135041/145166322-36f82c31-f9ed-4963-acda-7f8e3fa74da7.png">
