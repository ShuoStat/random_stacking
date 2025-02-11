# random_stacking

We propose Random Stacking, an extension of basic stacking that incorporates additional randomness through bootstrapped sampling and random feature selection. Random Stacking was applied to integrate clinical and omics data for risk stratification, and its predictive performance was compared against basic stacking and other methods using publicly available datasets. Additionally, a simulation study was conducted to assess its predictive performance across different scenarios and to explore its properties, with a particular focus on the sources of prediction error. 

### Project Structure

```
.
├── codes.realdata/
├── codes.simulation/
├── data/
├── output/
├── results/
└── README.md
```
`codes.realdata` codes for reproducing the results of real-world data application.  
`codes.simulation` codes for reproducing the results of simulation study.  
`data` real-world datasets used for model building.  
`output` intermediate outputs.    
`results` final results.

### Files in codes.realdata

`mainFunction.R` main functions to generate   
`fit.fun.R` functions for training models  
`random.stk_v4.R` random stacking functions (version 4)  
`stk.cox_v6.R` basic stacking for survival predictions (version 6)  
`stk.glm_v5.R` basic stacking for generalized models (version 5)  
`helper.R` helper function  
`ipflasso.R` IPFLasso  

### Files in codes.simulation

`simulation.R` main script for the simulation and result generation  
`fit.fun.simulation.R` training models   
`random.stk_v6.R` random stacking functions (version 6)  
`stk.glm_v6.R` basic stacking for generalized models (version 6)  
`baglasso.R` bagged lasso  

### Reproducing the Results

Run *mainFunction.R* step by step to generate the results of real-world application  
```
R /codes.realdata/mainFunction.R
```
Run *simulation.R* step by step to generate the results of simulations  
```
R /codes.simulation/simulation.R
```

### Contact

Maintainer: Shuo Wang  
GitHub: https://github.com/ShuoStat/  
Email: wangsures@foxmail.com  

### License

MIT License © Shuo Wang 2025

