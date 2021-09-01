##Initializing packages
using DifferentialEquations
using Plots
using LinearAlgebra
using Dierckx
using DelimitedFiles


##Defining constants
#Most variables in this script are given in SI units, except from a couple "convenience units"

g = 9.81;
d = 0.0427;
m = 0.0459; 
ρ = 1.225;

α = π*d^2*ρ/(8*m);
ω⃗0 = [0;-200;0];
τ = 24.5;


##Reading data
#The data is imported from CSV files, which contain x, y, c(x,y) in three columns for all measured points 

dragData = readdlm("cDragHex.csv");
magnusData = readdlm("cMagnusHex.csv");

dragX = dragData[:,1];
dragY = dragData[:,2];
dragZ = dragData[:,3];

magnusX = magnusData[:,1];
magnusY = magnusData[:,2];
magnusZ = magnusData[:,3];


##Splining
#Because of fitting problems when doing a spline on unstructured points in 2D, the spline is constructed by first interpolating in 1D along a every measured velocity

#Drag coefficient
xdEval = LinRange(0,maximum(dragX),200);
yd = unique(dragY);
zdVal = zeros(200,length(yd));

indDrag = findfirst.(isequal.(yd), [dragY])
push!(indDrag,length(dragY)+1)

for i in 0:(length(yd)-1)
    xd1D = dragX[indDrag[i+1]:(indDrag[i+2]-1)];
    zd1D = dragZ[indDrag[i+1]:(indDrag[i+2]-1)];
    cd1D = Spline1D(xd1D,zd1D);

    zdVal[:,i+1] = evaluate(cd1D,xdEval);
end

cDrag = Spline2D(xdEval,yd,zdVal);

#Surface plot of cDrag
ydEval = LinRange(minimum(yd),maximum(yd),200); 
zSpl = evalgrid(cDrag,xdEval,ydEval);
surface(xdEval, ydEval, zSpl'; leg=:false, camera=(30,65), xl="ω", yl="v", zl="cD")
scatter!(dragX, dragY, dragZ; ms=2, mc=:red, ma=0.4)

#Magnus coefficient
xmEval = LinRange(0,maximum(magnusX),200);
ym = unique(magnusY);
zmVal = zeros(200,length(ym));

indMagnus = findfirst.(isequal.(ym), [magnusY])
push!(indMagnus,length(magnusY)+1)

for i in 0:(length(ym)-1)

    xm1D = magnusX[indMagnus[i+1]:(indMagnus[i+2]-1)];
    zm1D = magnusZ[indMagnus[i+1]:(indMagnus[i+2]-1)];
    cm1D = Spline1D(xm1D,zm1D);

    zmVal[:,i+1] = evaluate(cm1D,xmEval);
end

cMagnus = Spline2D(xmEval,ym,zmVal);


##Wind

function v⃗W(x,y,z)
    return [0,0,0];
end


##Differential equation

function movementeq!(dudt,u,p,t)
    v⃗ = u[4:6] - v⃗W(u[1],u[2],u[3])
    v = norm(v⃗);
    v⃗N = v⃗/v;

    ω⃗ = ω⃗0*ℯ^(-t/τ)
    ω = norm(ω⃗)
    ω⃗N = ω⃗/ω;
    ωRPM = ω*60/(2*π);

    vM = v*sqrt(1-(ω⃗N⋅v⃗N)^2)

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

#Starting parameters: velocity in m/s, tee angle in degrees, spin vector in rad/s
v0 = 76;
ϕ = 12;

# Conversion of starting parameters into Cartesian coordinates of starting point and velocity
u0 = [0.0;0.0;0.02;cos(ϕ*π/180)*v0;0.0;sin(ϕ*π/180)*v0];
tspan = (0.0,10.0);
prob = ODEProblem(movementeq!,u0,tspan);

condition(u,t,integrator) = u[3];
affect!(integrator) = terminate!(integrator);
cb = ContinuousCallback(condition,affect!);
sol = solve(prob,Tsit5(),callback=cb);


##Plotting

plot(sol,vars=(1,2,3); camera=(15,30), zlims = (0,maximum(sol[3,:])+2), aspect_ratio=:equal,label="Trajectory")
landing = push!(sol[end][1:2],0);
scatter!([landing[1]],[landing[2]],[landing[3]]; ms=4, mc=:red, ma=0.4,label="Landing Point")

