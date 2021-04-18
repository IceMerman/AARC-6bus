# A New Affinely Adjustable Robust Model for Security Constrained Unit Commitment UnderUncertainty

Input data for the 6 bus test case and a modified version of the IEEE 24 bus system.  

This repository contains the information to optimize an affinely adjustable robust model for the unit commitment with "User Cuts".  

- [Description](#description)
- [Contact](#contact)
- [Files](#files)
- [Citation](#citation)
- [Clarification](#clarification)
- [License](#license)

## Description

Currently, optimization models for the safe and reliable operation of power systems deal with two major challenges: the first one is the reduction of the computational load when considering N-1 contingencies; the second one is the adequate modeling of the uncertainty of intermittent generation and demand.
This paper proposes a new affinely adjustable robust model to solve the security constrained unit commitment problem considering these sources of uncertainty. Linear decision rules, which take into account the forecasts and forecast errors of the different sources of uncertainty are used for the affine formulation of the dispatch variables, thus allowing the tractability of the model. Another major novelty is that the evaluation of the N-1 security constraints is performed by incorporating a novel method, proposed in the literature, based on the user-cuts concept. This method efficiently and dynamically adds only the binding N-1 security constraints, increasing the computational efficiency of the model when transmission line contingencies are considered. Finally, Monte Carlo simulations on the post-optimization results were run to demonstrate the effectiveness, feasibility and robustness of the solutions provided by the proposed model.

## Contact 

Juan Esteban Sierra-Aguilar, Universidad de Antioquia, juane.sierra@udea.edu.co  
Cristian Camilo Marín-Cano, Universidad de Antioquia, cristian1013@gmail.com      
Jesús M. López-Lezama, Universidad de Antioquia, jmaria.lopez@udea.edu.co   
Álvaro Jaramillo-Duque, Universidad de Antioquia, alvaro.jaramillod@udea.edu.co   
Juan G. Villegar, Universidad de Antioquia, juan.villegas@udea.edu.co  

## Files

# 6 Bus system

* Note for this system the bus index begging on 1  

6BUS_PTDF.inc: Power Transfer Distribution Factor for 6 bus system.  
6BUS_LODF.inc: Line Outage Distribution Factor for 6 bus system.  
6BUS_kptdf.inc: Outafe Transfer Distribution Factro for 6 bus system.  
6BUS_OTDF.inc: Outafe Transfer Distribution Factro for 6 bus system.  
wOut.inc: Wind power output in MW.  
wData.inc: Wind generator ubication.  
sOut.inc: PV power output in MW.  
sData.inc: PV generator ubication.  
lineData.inc: data for lines.  
genData.inc: data for hidrotermic generation units.  
d.inc: load profile for evey bus and hour.  

# 24 Bus system

* Note for this system the bus index begging on 0  

IEEE24-modified/24BUS_PTDF.inc: Power Transfer Distribution Factor for custom IEEE 24 bus system.  
IEEE24-modified/24BUS_LODF.inc: Line Outage Distribution Factor for custom IEEE 24 bus system.  
IEEE24-modified/24BUS_OTDF.inc: Outage Transfer Distribution Factor for custom IEEE 24 bus system.  
IEEE24-modified/24BUS_wOut.inc: Wind power output in MW.  
IEEE24-modified/24BUS_wData.inc: Wind generator ubication.  
IEEE24-modified/24BUS_sOut.inc: PV power output in MW.  
IEEE24-modified/24BUS_sData.inc: PV generator ubication.  
IEEE24-modified/24BUS_lineData.inc: data for lines.  
IEEE24-modified/24BUS_genData.inc: data for hidrotermic generation units  
IEEE24-modified/24BUS_d.inc: load profile for evey bus and hour.   

# images

img/6bus_FO_Tau.svg: objective function and execution time for different past information considerations, 6 bus system.  
img/6busMC_fullScale.svg: Montecarlo simulation, full scale, 6 bus system.  
img/6busMC_Zoom.svg: Zoomed montecarlo simulation, full scale, 6 bus system.  
img/24bus_FO_Tau.svg: objective function and execution time for different past information considerations, 24 bus system.   
img/MC_24b_fullScale.svg: Montecarlo simulation, full scale, 24 bus system.  
img/MC_24b_Zoom.svg: Zoomed montecarlo simulation, full scale, 24 bus system.  

# scripts

scripts/iter.sh: script to run multiple AARC optimization.  
scripts/iterMC.sh: script to run multiple montecarlo simulation.  
scripts/past_loop.bat: script to run multiple AARC formulation (windows).  
scripts/plotMC.py: script to plot values from montecarlo extracted information.  
scripts/plottingAARO.py: script to plot extracted values from AARC model.  
scripts/readGDX.py: script to extract data from a set of gdx files.   

## Citation

If you use this material, please refer to:  

[Autors, “A New Affinely Adjustable Robust Model forSecurity Constrained Unit Commitment UnderUncertainty,”Energies, vol. xx, no. x, 2021.](https://www.mdpi.com/xxxx-xxxx/xx/x/xxxx)

## Clarification

This work is a continuation of the  previous publications:  
* [C. C. Marín-Cano, J. E. Sierra-Aguilar, J. M. López-Lezama,  Jaramillo-Duque, andW. M. Villa-Acevedo, “Implementation of user cuts and linear sensitivity factors to im-prove the computational performance of the security-constrained unit commitment pro-blem,”Energies, vol. 12, no. 7, 2019.](https://www.mdpi.com/1996-1073/12/7/1399)
* [Marín-Cano CC, Sierra-Aguilar JE, López-Lezama JM, Jaramillo-Duque Á, Villegas JG, "A Novel Strategy to Reduce Computational Burden of the Stochastic Security Constrained Unit Commitment Problem," Energies, vol 13, no. 15, 2020;.](https://www.mdpi.com/1996-1073/13/15/3777)

## License

This code is provided under a BSD license as part of DEMIERI project.  

[License file](../master/LICENSE)

## Gratefulness

[Colciencias](https://colciencias.gov.co)  

![alt tag](https://minciencias.gov.co/sites/default/files/logo-minciencias_1.png)  

![alt tag](https://github.com/IceMerman/TransformerSoltion/blob/master/logoUDEA.png)  

![alt tag](https://github.com/IceMerman/TransformerSoltion/blob/master/gimel.png)  
