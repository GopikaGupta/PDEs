using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomainSets, Plots

# Define Parameters and Variables
@parameters t x
@variables u(..)  

Dt = Differential(t)
Dxx = Differential(x)^2  # Second derivative in space

# Defining the Schrödinger Equation 
eq = [im*Dt(u(t,x))~Dxx(u(t,x))]  ### consider V(x)=0

# Initial and Boundary Conditions
bcs=[u(0,x)~sin(2*pi*x),
       u(t,0)~ 0,
       u(t,1)~0]

# Define the domain
domain = [t ∈ Interval(0.0, 1.0),
          x ∈ Interval(0.0, 1.0)]

# Define the PDE System
@named pdesys = PDESystem(eq, bcs, domain, [t, x], [u(t, x)])

# Discretization (Method of Lines)
discretization = MOLFiniteDifference([x => 100], t)

# Convert PDE to ODE Problem
prob = discretize(pdesys, discretization)

# Solve the ODE System
sol = solve(prob, TRBDF2(), saveat=0.01)

# Extract Discrete Space and Time Points
discx=sol[x]
disct=sol[t]
discu=sol[u(t,x)]

# Animate Real and Imaginary Parts of the Wave Function
anim = @animate for i in 1:length(disct)
       u = discu[i, :]
       plot(discx, [real.(u), imag.(u)], ylim = (-1.5, 1.5), title = "t = $(disct[i])", xlabel = "x", ylabel = "u(t,x)", label = ["Real (u)" "Imaginary (u)"])
   end
gif(anim, "schroedinger.gif", fps = 10)
