## Modelling  2D spots pattern

using DifferentialEquations, LinearAlgebra, Plots

## Generate the constants
p = (1,1,30,0.1,6.7,1,40,40) # DA, DI, DC, kA, kI, kS, kminus, kplus
N = 300
#Discretization of Laplacian with von Neumann boundary conditions
Ax = Array(Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1]))
Ay = copy(Ax)
Ax[2,1] = 2.0
Ax[end-1,end] = 2.0
Ay[1,2] = 2.0
Ay[end,end-1] = 2.0

function basic_version!(dr,r,p,t)
    DA, DI, DC, kA, kI, kS, kminus, kplus = p
  #Here defining A, I, S, C
  a = r[:,:,1]
  i = r[:,:,2]
  s = r[:,:,3]
  c = r[:,:,4]
  Da = DA*(Ay*a + a*Ax)
  Di = DI*(Ay*i + i*Ax)
  Dc = DC*(Ay*c + c*Ax)

  dr[:,:,1] = Da .+ kminus.*c.-kplus.*a.*i.-kA.*a .+ (1 ./(1 .+ (s ./0.01).^8))
  dr[:,:,2] = Di .+ kminus.*c.-kplus.*a.*i.-kI.*i .+ (10 ./(1 .+ (s ./0.01).^2))
  dr[:,:,3] = a .-kS.*s
  dr[:,:,4] = Dc .+ kplus.*a.*i - kminus.*c
end

DA, DI, DC, kA, kI, kS, kminus, kplus = p

r0 = zeros(N,N,4)
r0[:,:,1] .= 0.02 .+ 0.001.*rand(N,N)
r0[:,:,2] .= 0.26 .+ 0.001.*rand(N,N)
r0[:,:,3] .= 0.02 .+ 0.001.*rand(N,N)
r0[:,:,4] .= 0.006 .+ 0.0001.*rand(N,N)

prob = ODEProblem(basic_version!,r0,(0.0,3000.0),p)
sol = solve(prob,ROCK2())

p3 = heatmap(sol[end][:,:,1],clim=(0.0,0.04),title = "Activator concentration - end")


## Making the figure nicer, set tresholds
sol_transformed = deepcopy(sol)  # Make a copy of the solution to modify

# Iterate over the solution and apply the transformation
for i in eachindex(sol.u)
    sol_transformed.u[i] .= ifelse.(sol.u[i] .< 0.03, 0.0, 1.0)
end
sol_transformed

## Define custom color map and plot
custom_cmap = cgrad([RGB(1.0, 1.0, 1.0), RGB(0.0, 114/255, 178/255)], [0.0, 1.0])

p3 = heatmap(
    sol_transformed[end][:,:,1], 
    clim=(0.0, 1.0),  # Data range
    color=custom_cmap, 
    title="Activator concentration - end"
)
