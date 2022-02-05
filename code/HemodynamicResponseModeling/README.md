# HemodynamicResponseModeling
modeling hemodynamic response functions of the BOLD signal using Windkessel-Balloon model

Friston KJ, Harrison L, Penny W (2003) Dynamic causal modelling. Neuroimage 19:1273–1302.



### Hemodynamic state equations

The remaining state variables of each region are biophysical states engendering the BOLD signal and mediate the
translation of neuronal activity into hemodynamic responses. Hemodynamic states are a function of, and only of,
the neuronal state of each region. These equations have been described elsewhere (Friston et al., 2000) and constitute a
hemodynamic model that embeds the _Balloon–Windkessel model_ (Buxton et al., 1998; Mandeville et al., 1999). 

In brief, for the ith region, neuronal activity $z_i$ causes an increase in a vasodilatory signal $s_i$ that is subject to autoregulatory feedback. 
Inflow $f_i$ responds in proportion to this signal with concomitant changes in blood volume $v_i$ and deoxyhemoglobin content $q_i$:

$$
\frac{s}{t}_i = z_i - 7kappa_I s_i - \gamma_i (f_i -1) 
\frac{f}{t}_i = s_i
$$
