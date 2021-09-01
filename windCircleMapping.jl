##Wind
function logWind(z,p)
    r = 3;
    c = [-1,5,0];
    return (c + [r*cos(p),r*sin(p),0]).*(log(abs(z+1))/log(11))
end


##Differential equation

#Starting parameters: velocity in m/s, tee angle in degrees, spin vector in rad/s
v0 = 76;
ϕ = 12;
ω⃗0 = [0;-200;10];

# Conversion of starting parameters into Cartesian coordinates of starting point and velocity
ϕ_arc = ϕ*π/180
u0 = [0.0;0.0;0.02;cos(ϕ_arc)*v0;0.0;sin(ϕ_arc)*v0]

function movementParaWind!(dudt,u,p,t)
    v⃗ = u[4:6] - logWind(u[3],p)
    v = norm(v⃗);
    v⃗N = v⃗/v;

    ω⃗ = ω⃗0*ℯ^(-t/τ);
    ω = norm(ω⃗);
    ω⃗N = ω⃗/ω;
    ωRPM = ω*60/(2*π);

    vM = v*sqrt(1-(ω⃗N⋅v⃗N)^2);

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


##Solving the ODE

tspan = (0.0,10.0);
p=0;
landings = zeros(36,2);

condition(u,t,integrator) = u[3]
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)

prob = ODEProblem(movementParaWind!,u0,tspan,p)
sol = solve(prob,Tsit5(),callback=cb)
p1 = plot(sol,vars=(1,2,3);lcolor=:"Dodger Blue", camera=(15,30))
landings[1,:] = sol[end][1:2];

for i in 1:35
    p=i*π/18;
    prob = ODEProblem(movementParaWind!,u0,tspan,p)

    sol = solve(prob,Tsit5(),callback=cb)
    if iseven(i)
        plot!(sol,vars=(1,2,3); linecolor=:"Dodger Blue", camera=(15,30), zlims = (0,maximum(sol[3,:])+2), aspect_ratio=:equal,leg=false);
    end
    landings[i+1,:] = sol[end][1:2]
end
landingsX = push!(landings[:,1],landings[1,1]);
landingsY = push!(landings[:,2],landings[1,2]);
p2 = plot(landingsX,landingsY; lcolor=:"Dodger Blue", aspect_ratio=:equal, leg=false, framestyle=:origin);
#plot(p1,p2; layout = (2,1), leg = false, wsize=(600,800))
plot(p2, wsize=(600,200))