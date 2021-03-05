# Monte-Carlo-WF-Mathematica

![](https://img.shields.io/static/v1?label=Version&message=0.1&color=red)

A Mathematica library to deal with Monte Carlo wave functions and Stochastic Master Equations. The package can be loaded with


```(Mathematica)
Package["Monte-Carlo-WF-0.1.wl"]
``` 

You can look up the syntax for the package functions using

```
? Function_Name
```

such as 

```(Mathematica)
?PhotoDetection
```

which returns

```
PhotoDetection[\[Psi]0_,H_, c_, nsteps_,\[CapitalDelta]t_] 
Returns the stochastic evolution given an initial wave function, the Hamiltonian, the jump operator and the integration parameters. 
The implementation is based on a jump-type Stochastic Schrödinger Equation
```
