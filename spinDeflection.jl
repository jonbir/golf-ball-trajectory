##Differential equation

#Starting parameters: velocity in m/s, tee angle in degrees, spin vector in rad/s
v0 = 76;
ϕ = 12;


# Conversion of starting parameters into Cartesian coordinates of starting point and velocity
ϕ_arc = ϕ*π/180
u0 = [0.0;0.0;0.02;cos(ϕ_arc)*v0;0.0;sin(ϕ_arc)*v0]


function movementParaSpin!(dudt,u,p,t)
    v⃗ = u[4:6] - wind(u[1],u[2],u[3])
    v = norm(v⃗);
    v⃗N = v⃗/v;

    ω⃗ = [0,-200*cos(p*π/180),200*sin(p*π/180)]*ℯ^(-t/τ)
    ω = norm(ω⃗)
    ω⃗N = ω⃗/ω;

    vM = v*sqrt(1-(ω⃗N⋅v⃗N)^2)

    ωRPM = ω*60/(2*π);

    CD = cDrag(ωRPM,v);
    CM = cMagnus(ωRPM,vM);

    aG = [0;0;-g] 
    aD = CD*α*v^2*(-v⃗N)
    aM = CM*α*vM^2*cross(ω⃗N,v⃗N)

    a = aG + aD + aM ;

    dudt[1] = u[4]
    dudt[2] = u[5]
    dudt[3] = u[6]
    dudt[4] = a[1]
    dudt[5] = a[2]
    dudt[6] = a[3]
end


##Solving the ODE and plotting

tspan = (0.0,10.0);
p=-30;
prob = ODEProblem(movementParaSpin!,u0,tspan,p);

condition(u,t,integrator) = u[3]
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)

sol = solve(prob,Tsit5(),callback=cb)
p1 = plot(sol,vars=(1,2,3); camera=(15,30), aspect_ratio=:equal,leg=false)

for i in 1:20
    p=-30+3*i;
    prob = ODEProblem(movementeq!,u0,tspan,p)

    sol = solve(prob,Tsit5(),callback=cb)
    plot!(sol,vars=(1,2,3); lcolor=:"Dodger Blue", camera=(15,30), aspect_ratio=:equal,leg=false);
end
plot(p1, xlims=(0,250), zlims=(0,30))