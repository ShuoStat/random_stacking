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
`codes.realdata` Scripts for reproducing the results of the real-world data application.    
`codes.simulation` Scripts for reproducing the results of the simulation study.   
`data` Real-world datasets used for model development.    
`output` Intermediate results and temporary outputs.      
`results` Final processed results. 

### Files in codes.realdata

`mainFunction.R` Core functions for generating results  
`fit.fun.R` Functions for training models    
`random.stk_v4.R` Random stacking functions (version 4)   
`stk.cox_v6.R` Basic stacking implementation for survival prediction (version 6)   
`stk.glm_v5.R` Basic stacking implementation for generalized models (version 5)  
`helper.R` Utility and helper functions  
`ipflasso.R` Implementation of IPFLasso  

### Files in codes.simulation

`simulation.R` Main script for running simulations and generating results     
`fit.fun.simulation.R` Functions for training models in simulations    
`random.stk_v6.R` Random stacking functions (version 6)  
`stk.glm_v6.R` Basic stacking for generalized models (version 6)  
`baglasso.R` Bagged lasso  

### Reproducing the Results

Run *mainFunction.R* step by step to generate the results of real-world application  
```
R /codes.realdata/mainFunction.R
```
Run *simulation.R* step by step to generate the results of simulations  
```
R /codes.simulation/simulation.R
```

Of note, run on R version 4.3.3. Results may be slightly different with different R version.      

### Contact

Maintainer: Shuo Wang  
GitHub: https://github.com/ShuoStat/  
Email: wangsures@foxmail.com  

### License

CC BY-NC License © Shuo Wang 2025

