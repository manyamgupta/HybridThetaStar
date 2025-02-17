using PyPlot
using LinearAlgebra
using Dubins

function GenRandomObstacles(n, xy_bounds)
    
    obsList = Vector{Vector{Tuple{Float64, Float64}}}(undef,0)
    # centers
    cx = rand(1:xy_bounds[1], n)
    cy = rand(1:xy_bounds[2], n)
    centers = [cx cy]
    # centers = rand(1:xy_bounds[1],n,2)

    println(centers)
    for k in 1:n
        cntr = centers[k,:]
        
        numVerts = rand(1:5)+3
        psiList = rand(numVerts)*2*pi
        psiList = sort(psiList)
        obs = Vector{Tuple{Float64, Float64}}(undef, 0)
        for psi in psiList
            radius = minimum(xy_bounds)*.05*(1+rand())
            coords = (cntr[1]+radius*cos(psi), cntr[2]+radius*sin(psi))
            push!(obs, coords)
        
        end
        push!(obsList, obs)
    end
    return obsList
end

function RandomRiskFunctions(n, origin, xy_bounds)
    
    if origin[1]<-xy_bnds[1] || origin[2]<-xy_bnds[2]
        println("origin and xy_bounds are not compatibel")
        return Vector{Vector{Int64}}(undef, 0), Vector{Matrix{Float64}}(undef,0)
    end
    cx = rand(origin[1]:xy_bounds[1], n)
    cy = rand(origin[2]:xy_bounds[2], n)
    # centers = [cx cy]   
    centers = [[cx[i],cy[i]] for i in eachindex(cx)]
    covMatList = Vector{Matrix{Float64}}(undef,0)
    for k in 1:n
        
        var1 = rand(1:50000)+50000
        covMat = [[var1 0];[0 .1*var1]]
        rotAngle = rand(0:180)*pi/180
        rotMat = [[cos(rotAngle) -sin(rotAngle)];[sin(rotAngle) cos(rotAngle)]]
        covMat = rotMat*covMat*transpose(rotMat)
        push!(covMatList, covMat)
    
    end
  
  return centers, covMatList
end

function PlotObstacles(obstaclesList)
    # p = plot(aspect_ratio=:equal,legend=false)
    for obs in obstaclesList
        ov = Matrix{Float64}(undef, 0,2)
        for p in obs ov = vcat(ov, [p[1] p[2]]) end

        ov = vcat(ov, [ov[1,1] ov[1,2]] )
        
        PyPlot.plot(ov[:,1], ov[:,2],color="k", linewidth=2)
    end
    PyPlot.axis("equal")
end

function PlotRiskContour(args...)
    if length(args) == 2
        println("function no defined for two args")
    end
    if length(args) == 3
        PlotRiskSumContour(args[1], args[2], args[3])
    end
end

function PlotRiskSumContour(riskParams, dimsXY, ax, cmp=PyPlot.get_cmap("Reds"))
    
        
    xVec = LinRange(0, dimsXY[1], 501)
    yVec = LinRange(0, dimsXY[2], 501)
    X, Y = meshgridpy(xVec, yVec)
        
    riskMat = [[SumRiskFunctions(riskParams, (xCoord, yCoord)) for xCoord in xVec] for yCoord in yVec]
    ax.contour(X, Y, riskMat, cmap=cmp)

    return
end

function meshgridpy(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return collect(X'), collect(Y')
end


function InflexionPoints(dbPath)
    
    dubPathLength = dubins_path_length(dbPath)
    segLengths = dbPath.params*dbPath.ρ
    if segLengths[1] <0. sl1 = 0.
    elseif segLengths[1] > dubPathLength sl1 = dubPathLength
    else sl1 = segLengths[1] end

    if segLengths[1]+segLengths[2] <0. sl2 = 0.
    elseif segLengths[1]+segLengths[2] > dubPathLength sl2 = dubPathLength
    else sl2 = segLengths[1]+segLengths[2] end

    _,conf1 = dubins_path_sample(dbPath, sl1)
    _,conf2 = dubins_path_sample(dbPath, sl2)
    if isnothing(conf1) || isnothing(conf2)
        println("dbPath: ", dbPath)
    end
    return (conf1[1], conf1[2], conf1[3]), (conf2[1], conf2[2], conf2[3]), segLengths
end 

function DubinsPathTypeString(dbPath)

    if dbPath.path_type == DubinsPathType(0)
        return ["L", "S", "L"]    
    elseif dbPath.path_type == DubinsPathType(1)
        return ["L", "S", "R"]
    elseif dbPath.path_type == DubinsPathType(2)
        return ["R", "S", "L"]
    elseif dbPath.path_type == DubinsPathType(3)
        return ["R", "S", "R"]
    elseif dbPath.path_type == DubinsPathType(4)
        return ["R", "L", "R"]
    elseif dbPath.path_type == DubinsPathType(5)
        return ["L", "R", "L"]        
    else
        return ["N"] # unknown type
    end
end

function WayPointsSegment(startConf, segLength, segType, rho, deltaL)

    pt1 = startConf[1:2]
    t1 = startConf[3]
    ndisc = max(2,Int(ceil(segLength/deltaL)))
    if segType == "S"
        pt2 = [startConf[1], startConf[2]] + segLength*[cos(t1), sin(t1)]
        
        wayPts_x = LinRange(pt1[1], pt2[1], ndisc)
        wayPts_y = LinRange(pt1[2], pt2[2], ndisc)
        
    elseif segType == "L" || segType == "R"

        if segType =="L"
            rotSense = 1 
        else 
            rotSense = -1 
        end
        center = [pt1[1], pt1[2]] + rho*[cos(t1+rotSense*pi/2), sin(t1+rotSense*pi/2)]

        alVec = LinRange(t1-rotSense*pi/2 , t1 -rotSense*pi/2 +rotSense*segLength/rho, ndisc)

        wayPts_x = center[1] .+ rho*cos.(alVec)
        wayPts_y = center[2] .+ rho*sin.(alVec)

    end

    return [[x, wayPts_y[i]] for (i,x) in enumerate(wayPts_x)]
end

function PlotSegment(startConf, segLength, segType, rho, fmt)

    pt1 = startConf[1:2]
    t1 = startConf[3]
    if segType == "S"
        pt2 = [startConf[1], startConf[2]] + segLength*[cos(t1), sin(t1)]
        PyPlot.plot([pt1[1],pt2[1]], [pt1[2],pt2[2]], color=fmt.color, linewidth=fmt.linewidth, linestyle=fmt.linestyle)

        finalConf = [pt2[1], pt2[2], t1]
        if fmt.endPoints
            PyPlot.scatter([pt1[1],pt2[1]], [pt1[2],pt2[2]], marker=fmt.marker, color="k")
        end
    elseif segType == "L" || segType == "R"

        if segType =="L"
            rotSense = 1 
        else 
            rotSense = -1 
        end
        center = [pt1[1], pt1[2]] + rho*[cos(t1+rotSense*pi/2), sin(t1+rotSense*pi/2)]

        alVec = LinRange(t1-rotSense*pi/2 , t1 -rotSense*pi/2 +rotSense*segLength/rho,100)

        tc_x = center[1] .+ rho*cos.(alVec)
        tc_y = center[2] .+ rho*sin.(alVec)
        PyPlot.plot(tc_x, tc_y, color=fmt.color, linewidth=fmt.linewidth, linestyle=fmt.linestyle) 
        if fmt.endPoints
            PyPlot.scatter([tc_x[1], tc_x[end]], [tc_y[1], tc_y[end]], marker=fmt.marker, color="k")
        end
        t2 = mod(t1  + rotSense*segLength/rho, 2*pi)
        finalConf = (tc_x[end], tc_y[end], t2)
    end

    return finalConf
end


function PlotLineSegment(lstart, lend)

    PyPlot.plot([lstart[1], lend[1]], [lstart[2], lend[2]])
    PyPlot.scatter([lstart[1], lend[1]], [lstart[2], lend[2]],marker="x")

end

function PlotLine(pt1, pt2, lambounds)

    lamVec = collect(range(lambounds[1],lambounds[2],length=11))
    line = pt1' .* (1 .- lamVec) + pt2' .* lamVec
    PyPlot.plot(line[:,1], line[:,2])
    
end

function PlotDubins(dubPath, fmt)
    _, ds = dubins_path_sample_many(dubPath,.01)
    dp = vcat(ds'...)
    PyPlot.plot(dp[:,1], dp[:,2], color=fmt.color, linewidth=fmt.linewidth)
    PyPlot.axis("equal")

end

function PlotDubPathSegments(iniConf, pathMode, segLengths, rho, fmt)

    startConf = iniConf
    for k in eachindex(pathMode)
        startConf = PlotSegment(startConf, segLengths[k], pathMode[k], rho, fmt)
    end

    return
end

function PlotEllipse(center, axes, phi, pltFmt, confAlpha = 4.605)
    
    numSamples = 1000
    tVec = LinRange(0,2*pi, numSamples)
    x = zeros(numSamples)
    y = zeros(numSamples)
    # confAlpha = 4.605 # for 90% confidence level
    # confAlpha = 5.991 # for 95% confidence level
    # confAlpha = 9.210 # for 99% confidence level

    sqrtAlpha = sqrt(confAlpha)
    # for (ind, t) in enumerate(tVec)
    #     x[ind] = sqrtAlpha*(axes[1]*cos(t))*cos(phi) - sqrtAlpha*(axes[2]*sin(t))*sin(phi)+center[1]
    #     y[ind] = sqrtAlpha*(axes[2]*sin(t))*cos(phi) + sqrtAlpha*(axes[1]*cos(t))*sin(phi)+center[2]
    # end
    
    x = [sqrtAlpha*(axes[1]*cos(t))*cos(phi) - sqrtAlpha*(axes[2]*sin(t))*sin(phi)+center[1] for t in tVec]
    y = [sqrtAlpha*(axes[2]*sin(t))*cos(phi) + sqrtAlpha*(axes[1]*cos(t))*sin(phi)+center[2] for t in tVec]
    
    PyPlot.plot(x,y,color=pltFmt.color, linewidth=pltFmt.linewidth, linestyle=pltFmt.linestyle)
    return x,y
end

function PtLiesInConfEllipse(pt, mean, invCovMat, confAlpha=3.219)
    # confAlpha = 3.219 # for 80% confidence level
    # confAlpha = 4.605 # for 90% confidence level
    # confAlpha = 5.991 # for 95% confidence level
    # confAlpha = 9.210 # for 99% confidence level
    ellipseEqnVal = [pt[1]-mean[1], pt[2]-mean[2]]'*invCovMat*[pt[1]-mean[1], pt[2]-mean[2]]
    if ellipseEqnVal < confAlpha
        return true
    else
        return false
    end
end