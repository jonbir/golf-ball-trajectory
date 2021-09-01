##Ball model comparison
#store solutions from golfBallTrajectory.jl in solHex and solCirc to call this plot

plot(solHex,vars=(1,2,3); camera=(15,30), ylims=(-60,60) ,zlims = (0,maximum(solHex[3,:])+2), label="Hexagonal")
landingHex = solHex[end][1:3];
scatter!([landingHex[1]],[landingHex[2]],[landingHex[3]]; ms=4, mc=:red, ma=0.4,label=false)

plot!(solCirc,vars=(1,2,3); label="Circular", xlims=(0,maximum(solHex[1,:])+20))
landingCirc = solCirc[end][1:3];
scatter!([landingCirc[1]],[landingCirc[2]],[landingCirc[3]]; ms=4, mc=:red, ma=0.4,label=false)


##No wind, headwind, tailwind comparison
#store solutions golfBallTrajectory.jl in solNoWind, solHeadWind and solTailwind to call this plot

plot(solNoWind,vars=(1,2,3); camera=(15,30), ylims=(-60,60) ,zlims = (0,maximum(solHeadWind[3,:])+2), label="no wind")
landingNoWind = solNoWind[end][1:3];
scatter!([landingNoWind[1]],[landingNoWind[2]],[landingNoWind[3]]; ms=4, mc=:red, ma=0.4,label=false)

plot!(solHeadWind,vars=(1,2,3); label="headwind")
landingHeadWind = solHeadWind[end][1:3];
scatter!([landingHeadWind[1]],[landingHeadWind[2]],[landingHeadWind[3]]; ms=4, mc=:red, ma=0.4,label=false)

plot!(solTailWind,vars=(1,2,3); label="tailwind", xlims=(0,maximum(solTailWind[1,:])+20))
landingTailWind = solTailWind[end][1:3];
scatter!([landingTailWind[1]],[landingTailWind[2]],[landingTailWind[3]]; ms=4, mc=:red, ma=0.4,label=false)

##Different wind models
function constWind(z)
	return [-5,0,0]
end

function linWind(z)
	return [-0.5,0,0].*z
end

function expWind(z)
	return [-5*ℯ^(1/10),0,0].*ℯ^(-1/z)
end

function logWind(z)
    if z>0
        return [-5,0,0].*((log(z+1))/log(11))
    else
        return [0,0,0]
    end
end

##Wind model comparison
#store solutions in solConst, solLin, solExp and solLog to call this plot
plot(solConst,vars=(1,2,3); camera=(15,30), ylims=(-60,60) ,zlims = (0,maximum(solLin[3,:])+2), label="constant")
landingConst = solConst[end][1:3];
scatter!([landingConst[1]],[landingConst[2]],[landingConst[3]]; ms=4, mc=:red, ma=0.4,label=false)

plot!(solLin,vars=(1,2,3); label="linear")
landingLin = solLin[end][1:3];
scatter!([landingLin[1]],[landingLin[2]],[landingLin[3]]; ms=4, mc=:red, ma=0.4,label=false)

plot!(solExp,vars=(1,2,3); label="exponential")
landingExp = solExp[end][1:3];
scatter!([landingExp[1]],[landingExp[2]],[landingExp[3]]; ms=4, mc=:red, ma=0.4,label=false)

plot!(solLog,vars=(1,2,3); label="logarithmic",  xlims=(0,maximum(solExp[1,:])+20))
landingLog = solLog[end][1:3];
scatter!([landingLog[1]],[landingLog[2]],[landingLog[3]]; ms=4, mc=:red, ma=0.4,label=false)