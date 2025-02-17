include("HRTsFocalSearchFuncs.jl")
include("HRTStar_FocalSelect.jl")
include("utils.jl")
include("HeurCost.jl")
using PyPlot
using FileIO, JLD2

pygui(true)

# rhoVec = [100,150,200,250]
# cellResVec = [100,150,200,250]
riskPercentVec = [0.1, 0.2, 0.3, 0.4]
w_vec = [1.1, 1.2, 1.3, 1.4]
# for j=3
for k=1:2
    # j=4
    # k=6
    fileNameStr = string("Scenarios/RandScenario_",string(k),".jld2")
    println(" ")
    println(fileNameStr)
    data = FileIO.load(fileNameStr)

    startConf = data["startConf"]
    goalConf = data["goalConf"]
    obsList = data["obsList"]
    obsList =[]
    # rho = data["rho"]
    rho=200.

    originCoords = data["origin"]
    xy_bnds = data["xy_bnds"]
    riskParams = data["riskParams"]
    riskOnSline = data["riskOnSline"]
    
    # riskPercent = riskPercentVec[j]
    riskPercent = 0.2
    maxRisk = riskPercent*riskOnSline
    riskMaxDist = 400.
    ######### Algorithm parameters
    # w = w_vec[j]
    w=1.2
    # cellRes = cellResVec[j]
    cellRes = 100.
    hdngRes = .5    
    arcLen = 1.0*cellRes
    pltExplFlag = false

    riskParams = (riskCntrs=riskParams.riskCntrs, invCovMats=riskParams.invCovMats, riskMaxDist=riskMaxDist, maxRisk=maxRisk)

    ###################### Run the scenario ######################
    stats = @timed HRTStarFocal(startConf, goalConf, rho, arcLen, obsList, riskParams, originCoords, xy_bnds, cellRes, hdngRes, w, pltExplFlag) 
    hrtsResult = stats.value
    compTime = stats.time
    pathCost = hrtsResult.cost
    status = hrtsResult.status
    println(hrtsResult.status)
    

    htState = hrtsResult.htState;
    htsPath = hrtsResult.path;
    
    ###################### Print the results ######################
    # for p in htsPath println(p) end
    println("Path Cost: ", hrtsResult.cost)
    println("Risk limit: ", maxRisk)
    println("Risk accumulated along the path: ", hrtsResult.intgrCons)
    println("compTime: ", compTime)

    ##################### Plot the path ######################
    deltaL = 20 #discretization parameter: length between two waypoints
    if hrtsResult.status == :success
        wayPoints = PathWayPoints(hrtsResult.path, rho, deltaL)
        wayPointsMat = reduce(vcat,transpose.(wayPoints))
        fig, ax = PyPlot.subplots(1,1)
        PlotResultScenario(hrtsResult, ax)
        PyPlot.plot(wayPointsMat[:,1], wayPointsMat[:,2], color="blue", linewidth=2, linestyle="-")
    end
end
# end