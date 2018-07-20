# GLM_and_Izhikevich

Associated with publication ([Weber & Pillow
2017](http://www.mitpressjournals.org/doi/abs/10.1162/neco_a_01021)), this repository contains Matlab code for:

1. generating simulated neural responses that reproduce a variety of dynamic spiking
behaviors using the Izhikevich neuron model [2,3].

2. fitting a Poisson generalized linear model (GLM) to neural responses.




To get started, open `demo_glm_fitting_and_simulation.m`.

`fit_glm.m` is the main function used to fit a Poisson GLM to data.




### References

1. Weber, A.I. & Pillow, J.W.  **Capturing the dynamical repertoire
   of single neurons with generalized linear
   models**. *Neural Computation*, 2017 29(12):3260-3289 [link](http://www.mitpressjournals.org/doi/abs/10.1162/neco_a_01021).

2. Izhikevich, E. M.  **Simple model of spiking neurons**. *IEEE Trans Neural Netw*, 2003 (14): 1569-1572.

3. Izhikevich, E. M.  **Which model to use for cortical spiking
   neurons?**  *IEEE Trans Neural
   Netw* 2004 (15): 1063-1070. 
