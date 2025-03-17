### Loading Packages
using OrdinaryDiffEq,ModelingToolkit,MethodOfLines,DomainSets


## Variables, Parameters and Derivatives
@parameters t x
@variables u(..)
Dt=Differential(t) ## d/dt
Dxx=Differential(x)^2  ## d/dx^2


#### 1D PDE and boundary conditions
V(x)=0.0
eq= [[im*Dt(u(t,x))] ~ [Dxx(u(t,x))+V(x)*u(t,x)]]

bcs=[u(0,x)~sin(2*pi*x),
     u(t,0)~0.0,
     u(t,1)~0.0]


##### Space and time domains (like time span in ODEs)
domain=[t ∈ Interval(0.0,1.0),
        x ∈ Interval(0.0,1.0)]


### PDE System
@named pdesys=PDESystem(eq, bcs, domain, [t,x],[u(t,x)])


#### Method of Line discretization
dx=100
order=2
discretization=MOLFiniteDifference([x=>dx],t)

### COnvert the PDE problem to ODEProblem
prob=discretize(pdesys, discretization)

using OrdinaryDiffEq
sol=solve(prob, Tsit5(), saveat=0.2)


### Plot results and compare with exact Solution
discrete_x=sol[x]
discrete_t=sol[t]
solu=sol[u(t,x)]



### Plots
using Plots