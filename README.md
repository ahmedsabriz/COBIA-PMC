# COBIA
CAPE-OPEN Unit Operation using COBIA SDK

SDK version 1.2.0.7

Unit version 0.3 (heater)

The unit emulates the behaviour of conventional unit operations. The goal is to ultimately develop a non-conventional unit once the development of all modules is matured.

The current state of the unit emulates a simple heater with two parameters. _Outlet temperature_ is an input that determines the temparture of _product1_ stream. _heat duty_ is an output parameter. An option parameter will be introduced to reverse the input/output parameter operation.

The results have been verified against conventional heater in COCO simulator except for the displayed value of work parameter in the energy stream:


![image](https://user-images.githubusercontent.com/80135041/128592023-c75ba1d5-86fd-4869-a546-19fa4753e40d.png)

