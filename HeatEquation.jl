#### Loading packages
using OrdinaryDiffEq,ModelingToolkit,MethodOfLines,DomainSets

u_exact=(x,t)->exp.(-t)*cos.(x)

## Variables, Parameters and Derivatives@parameters t x   ## Dependent variables
@variables u(..)   ## InDependent variables
Dt=Differential(t)   # Dt=du/dt
Dxx=Differential(x)^2 # Dxx=d/dx^2


#### 1 PDE and boundary conditions

eq=Dt(u(t,x)) ~ Dxx(u(t,x))
bcs=[u(0,x)~cos(x),
     u(t,0)~exp(-t),
     u(t,1)~exp(-t)*cos(1)]


##### Space and time domains (like time span in ODEs)
domain=[t ∈ Interval(0.0,1.0),
        x ∈ Interval(0.0,1.0)]

### PDE System
@named pdesys=PDESystem(eq, bcs, domain, [t,x],[u(t,x)])

#### Method of Line discretization
dx=0.1
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

plt=plot()
for i in eachindex(discrete_t)
    plot!(discrete_x,solu[i,:],label="Numerical, t=$(discrete_t[i])")
    scatter!(discrete_x,u_exact(discrete_x,discrete_t[i]),label="Exact, t=$(discrete_t[i])")
end
plt


##### Neumann Boundary Conditions