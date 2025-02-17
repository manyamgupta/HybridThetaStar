import Base
# using GeometricalPredicates
using Luxor
using DataStructures

"""Node of the state tree to explore"""
mutable struct HRTSNode
    id::Int64
    config::Tuple{Float64, Float64, Float64}
    parent::Union{HRTSNode, Nothing}
    edgType::String #this is the type of edge coming from the parent, "S" or "L" or "R"
    snapHash2D::UInt64
    snapHash3D::UInt64
    g::Float64 # cost from start to current state
    f::Float64 # total cost, g + h (heuristic)
    intgrCons::Float64
    intgrConsTotalLB::Float64

    function HRTSNode(id, state, parnt)
      # HRTSNode(id, state, parnt ) = new(id, state, parnt, "NS", 0,0, NaN, NaN, NaN, NaN)
      return new(id, state, parnt, "NS", 0,0, NaN, NaN, NaN, NaN)
    end
end

"order by f = g + h"
Base.isless(n1::HRTSNode, n2::HRTSNode) = Base.isless(n1.f, n2.f)

mutable struct HRTSSearchState
    openheap::Vector{HRTSNode}
    nodedict::Dict{Int64, HRTSNode}
    closedset::Set{UInt64} 
    startTime::Float64
    originCoordinates::Tuple{Float64, Float64}
    xy_bnds::Tuple{Float64, Float64}
    xy_inds_max::Tuple{Int64, Int64}
    cellRes::Float64
    hdngRes::Float64
    rho::Float64
    arcLen::Float64
    startConf::Tuple{Float64, Float64, Float64}
    goalConf::Tuple{Float64, Float64, Float64}
    obstaclesPolyList::Vector{Vector{Point}}    
    infeas_xyinds_hashlist::Set{UInt64}
    riskParams::NamedTuple{}
    # riskParams are named tuple of (riskCntrs, covMats, riskMaxDist, maxRisk)
end


# "Results structure"
struct HRTSResult
  status::Symbol
  path:: Vector{Tuple{Tuple{Float64, Float64, Float64}, String, Float64, Float64}} 
  # path::Vector{TupleTState}
  # edgTypes::Vector{String}
  cost::Float64
  cost_unSmooth::Float64
  intgrCons::Float64
#   htState::Union{HRTSSearchState{Tuple{Float64, Float64, Float64}, Float64}, Nothing}
  htState::Any
end

function InitializeSearchState(startConf, goalConf, obstaclesPolyList, riskGaussParams, infeas_xyinds_hashlist, originCoordinates, xy_bnds, xy_inds_max, cellRes, hdngRes, rho, arcLen, startHeur)
  
  
  strtNodeId = GetNodeId4mConf(startConf)       
  startNode = HRTSNode(strtNodeId, startConf, nothing)
  startNode.edgType= "Start"
  startNode.snapHash2D = hash(PointToIndex(startConf[1:2], originCoordinates, cellRes))
  startNode.snapHash3D = hash(ConfToIndex(startConf, originCoordinates, cellRes, hdngRes)) 
  startNode.g = 0.
  startNode.f = startHeur
  startNode.intgrCons = 0.
  startNode.intgrConsTotalLB = 0.


  openheap = [startNode]
  startHash = hash(startNode.config)
  nodedict = Dict(strtNodeId => startNode)
  # openSnapdNodeHashList = Vector{UInt64}(undef,0)
  closedset = Set{UInt64}()
  hrtsState = HRTSSearchState(openheap, nodedict, closedset, time(), 
  originCoordinates,xy_bnds, xy_inds_max, cellRes, hdngRes, rho, arcLen, startConf, goalConf,
  obstaclesPolyList, infeas_xyinds_hashlist,riskGaussParams )
  # hrtsState = 0
  return hrtsState
end

function CheckStartGoalFeasible(startConf, goalConf, originCoordinates,xy_bnds, cellRes, infeas_xyinds_hashlist)

  # if !CheckPtFeasible(startConf[1:2], hts_state.obstaclesPolyList, hts_state.originCoordinates, hts_state.xy_bnds)
  # println("Checking the feasibility of start and goal")
  if !CheckPtFeasibleDisc(startConf[1:2], originCoordinates,xy_bnds, cellRes, infeas_xyinds_hashlist )
    println("Start position infeasible")
    return false
  end

  # if !CheckPtFeasible(goalConf[1:2], hts_state.obstaclesPolyList, hts_state.originCoordinates, hts_state.xy_bnds)
  if !CheckPtFeasibleDisc(goalConf[1:2], originCoordinates,xy_bnds, cellRes, infeas_xyinds_hashlist )
    println("Goal position infeasible")
    return false
  end
  return true
end

function ObsPolygonList(obsList)
  
  obsPolyList = Vector{Vector{Point}}(undef,0)
  for obs in obsList
      polygn = [Point(p) for p in obs]
      push!(obsPolyList, polygn)
  end
  return obsPolyList
end

function CheckPtFeasible(coordsPt, obsPolyList, originCoordinates, xy_bnds )
# Feasibility is checked without discretizing
# @show coordsPt
  if coordsPt[1]<originCoordinates[1] || coordsPt[1]>originCoordinates[1]+xy_bnds[1] || coordsPt[2]<originCoordinates[2] || coordsPt[2]>originCoordinates[2]+xy_bnds[2]
    return false
  end

  for obs in obsPolyList
    if Point(coordsPt[1], coordsPt[2]) in obs
      return true
    end
    if isinside(Point(coordsPt[1], coordsPt[2]), obs, allowonedge=false)
        return false
    end
  end     
  return true
end

function CheckPtFeasibleDisc(coordsPt, hts_state )
  # Feasibility is checked with discretizing
  # returns true if point is feasible
  # @show coordsPt
  if coordsPt[1]<hts_state.originCoordinates[1] || coordsPt[1]>hts_state.originCoordinates[1]+hts_state.xy_bnds[1] || coordsPt[2]<hts_state.originCoordinates[2] || coordsPt[2]>hts_state.originCoordinates[2]+hts_state.xy_bnds[2]
    return false
  end

  pt_indices = PointToIndex(coordsPt, hts_state.originCoordinates, hts_state.cellRes)
  if hash(pt_indices) in hts_state.infeas_xyinds_hashlist
    return false
  end

  return true
end

function CheckPtFeasibleDisc(coordsPt, originCoordinates,xy_bnds, cellRes, infeas_xyinds_hashlist )
    # Feasibility is checked with discretizing
    # returns true if point is feasible
    # @show coordsPt
    if coordsPt[1]<originCoordinates[1] || coordsPt[1]>originCoordinates[1]+xy_bnds[1] || coordsPt[2]<originCoordinates[2] || coordsPt[2]>originCoordinates[2]+xy_bnds[2]
      return false
    end
  
    pt_indices = PointToIndex(coordsPt, originCoordinates, cellRes)
    if hash(pt_indices) in infeas_xyinds_hashlist
      return false
    end

    return true
end

function GenerateObstacleInds(xy_bnds,cellRes, obsPolyList, originCoordinates )
  # Discretizes the space into cells and checks four corners of each cell lies in any obstacle
  # Designates the cell to be infeasible if any corner lies inside obstacle
  # crreates the has of every pair of xy_indices that are infasible, and returns the list
  # To do: Can be optimized by precompiling list of all corner points, currently there is lot of repitition 

  xy_inds_max = (Int(floor((xy_bnds[1]/cellRes))), Int(floor(xy_bnds[2]/cellRes)) )
  x_inds_list = 0:xy_inds_max[1]
  y_inds_list = 0:xy_inds_max[2]

  obs_xyinds_hashlist = Set{UInt64}()
  obs_snap_coordsList = Matrix{Float64}(undef, 0,2)
  # obs_y_inds = []
  for x_ind in x_inds_list
    for y_ind in y_inds_list
      if !CheckPtFeasible((originCoordinates[1]+x_ind*cellRes, originCoordinates[2]+y_ind*cellRes),  obsPolyList, originCoordinates, xy_bnds)
        
        push!(obs_xyinds_hashlist, hash((x_ind, y_ind)) )
        push!(obs_xyinds_hashlist, hash((x_ind+1, y_ind)) )
        push!(obs_xyinds_hashlist, hash((x_ind-1, y_ind)) )
        push!(obs_xyinds_hashlist, hash((x_ind, y_ind+1)) )
        push!(obs_xyinds_hashlist, hash((x_ind, y_ind-1)) )
        # push!(obs_xyinds_hashlist, hash((x_ind+1, y_ind+1)) )
        # push!(obs_xyinds_hashlist, hash((x_ind-1, y_ind+1)) )
        # push!(obs_xyinds_hashlist, hash((x_ind+1, y_ind-1)) )
        # push!(obs_xyinds_hashlist, hash((x_ind-1, y_ind-1)) )
        obs_snap_coordsList = [obs_snap_coordsList; [originCoordinates[1]+x_ind*cellRes originCoordinates[2]+y_ind*cellRes]]
        obs_snap_coordsList = [obs_snap_coordsList; [originCoordinates[1]+(x_ind+1)*cellRes originCoordinates[2]+y_ind*cellRes]]
        obs_snap_coordsList = [obs_snap_coordsList; [originCoordinates[1]+(x_ind-1)*cellRes originCoordinates[2]+y_ind*cellRes]]
        obs_snap_coordsList = [obs_snap_coordsList; [originCoordinates[1]+x_ind*cellRes originCoordinates[2]+(y_ind+1)*cellRes]]
        obs_snap_coordsList = [obs_snap_coordsList; [originCoordinates[1]+x_ind*cellRes originCoordinates[2]+(y_ind-1)*cellRes]]

      end
    end
  end
  return obs_xyinds_hashlist, obs_snap_coordsList
end

function NeighborParams(rho, arcLen)
  arcAngle = arcLen/rho        
  cordLen = 2*rho*sin(arcAngle/2)
  return arcAngle, cordLen
end

function ConfToIndex(conf, hts_state)
  hdng = mod(conf[3], 2*pi)
  dH = hts_state.hdngRes-1e-6
  return ( Int(floor(( (conf[1]-hts_state.originCoordinates[1])/hts_state.cellRes))), Int(floor(((conf[2]-hts_state.originCoordinates[2])/hts_state.cellRes))), Int(floor((hdng/dH))) )
  
end

function ConfToIndex(conf, originCoordinates, cellRes, hdngRes)
    hdng = mod(conf[3], 2*pi)
    return ( Int(floor(( (conf[1]-originCoordinates[1])/cellRes))), Int(floor(((conf[2]-originCoordinates[2])/cellRes))), Int(floor((hdng/hdngRes))) )    
end

function PointToIndex(pt, originCoordinates, cellRes)
  
  return ( Int(floor(( (pt[1]-originCoordinates[1])/cellRes))), Int(floor(((pt[2]-originCoordinates[2])/cellRes))) )
  
end

function GetNodeId4mConf(conf)
  hdng = mod(conf[3], 2*pi)*180/pi
  return Int64( floor(conf[1]*1e10 + conf[2]*1e6 + hdng*1e2))
end

function Neighbors(nodeConfig, hts_state, arcAngle, cordLen)

    # return [(nodeConfig[1]+arcLen, nodeConfig[2], 0.), (nodeConfig[1]-arcLen, nodeConfig[2], 0.),  (nodeConfig[1], nodeConfig[2]+arcLen, 0.), (nodeConfig[1], nodeConfig[2]-arcLen, 0.)]
    hdng = nodeConfig[3]
    ndisc = 2
    edgTypes = ["S", "L", "R"]
    arcLenVec = LinRange(0., hts_state.arcLen, ndisc)
    coordsSline = (nodeConfig[1]+hts_state.arcLen*cos(hdng), nodeConfig[2]+hts_state.arcLen*sin(hdng), hdng)
   
    arcAngleVec = arcLenVec/hts_state.rho
    cordLenVec = 2*hts_state.rho*sin.(arcAngleVec/2)
    
    coordsLtTurn = (nodeConfig[1]+cordLen*cos(hdng+arcAngle/2), nodeConfig[2]+cordLen*sin(hdng+arcAngle/2), hdng+arcAngle)
    
    coordsRtTurn = (nodeConfig[1]+cordLen*cos(hdng-arcAngle/2), nodeConfig[2]+cordLen*sin(hdng-arcAngle/2), hdng-arcAngle)

    chldCoordsList = [coordsSline, coordsLtTurn, coordsRtTurn]
    feasNbrStatesList = Vector{Tuple{Float64, Float64, Float64}}(undef,0)
    edgTypesList = Vector{String}(undef, 0)
    for (k, nbr_states) in enumerate(chldCoordsList)
      # nbr_snap_hash = hash(ConfToIndex(nbr_states, hts_state)[1:2])
      if nbr_states[1] >= hts_state.originCoordinates[1] && nbr_states[1] <= hts_state.xy_bnds[1] && nbr_states[2] >= hts_state.originCoordinates[2] && nbr_states[2] <= hts_state.xy_bnds[2]
        nbr_snap_hash = hash(PointToIndex(nbr_states, hts_state.originCoordinates, hts_state.cellRes) ) 
        if !(nbr_snap_hash in hts_state.infeas_xyinds_hashlist)
          push!(feasNbrStatesList, nbr_states)
          push!(edgTypesList, edgTypes[k])
        end
      end
    end

  return feasNbrStatesList, edgTypesList
end

function Heuristic(state, heurCostDict, originCoordinates, cellRes)

    # return sqrt((neighbour[1]-goal[1])^2+(neighbour[2]-goal[2])^2)
    stateInds = PointToIndex(state, originCoordinates, cellRes)
    if hash(stateInds) in keys(heurCostDict)
      return 1.0*heurCostDict[hash(stateInds)].f      
    else
      println("key not found in the heuristic cost for the state: ", state, "; state indices: ", stateInds)
      return nothing
    end
end

function HeuristicDubins(curNodeConf, goalConf, rho)

  _, dbPath = dubins_shortest_path(collect(curNodeConf), collect(goalConf), rho)
  return dubins_path_length(dbPath)
end

function IsFeasibleAnalyticExpn(node, goal, hts_state)
  if abs(node.config[1]-goal[1])<50 && abs(node.config[2]-goal[2])<50 && abs(mod(node.config[3],2*pi)-mod(goal[3],2*pi)) < .5
    return true
  end
  _, dubPath = dubins_shortest_path(collect(node.config), collect(goal), hts_state.rho)
  if CheckDubinsPathFeas(dubPath, hts_state) 
    dbEdgIntgrCons, _ = IntgrConsDubPath(dubPath, hts_state)
    if node.intgrCons+dbEdgIntgrCons <= hts_state.riskParams.maxRisk  
      return true
    end
  end
  return false
end

function CheckDubinsPathFeas(dubPath, hts_state)
  if dubins_path_length(dubPath) < hts_state.arcLen
    _, ds = dubins_path_sample_many(dubPath, dubins_path_length(dubPath)/4)
  else
    _, ds = dubins_path_sample_many(dubPath, hts_state.arcLen/5.)
  end
  if isnothing(ds)
    println(dubPath)
    println(hts_state.arcLen)

  end
  dp = [(p[1],p[2]) for p in ds]
  # CheckFeas(p) = CheckPtFeasibleDisc(p, hts_state )
  # feasVec = CheckFeas.(dp)
  feasVec = [CheckPtFeasibleDisc(p, hts_state ) for p in dp]
  return prod(feasVec)
end

function CheckDubinsPathFeas2(dubPath, obsPolyList, originCoordinates, xy_bnds )
  # This function sampels the dubins path and checks the feasibility of each sample analytically
  # _, ds = dubins_path_sample_many(dubPath, 10.)
  _, ds = dubins_path_sample_many(dubPath, dubPath.ρ/4)
  if isnothing(ds)
    _, ds = dubins_path_sample_many(dubPath, dubins_path_length(dubPath)/4)
  end
  dp = [(p[1],p[2]) for p in ds]
  feasVec = [CheckPtFeasible(p, obsPolyList, originCoordinates, xy_bnds ) for p in dp]
  return prod(feasVec)
end

function Reconstructpath(node, hrts_state)
    pathConfigs = Vector{Tuple{Float64, Float64, Float64}}(undef,0)
    pathEdgeTypes = Vector{String}(undef,0)
    push!(pathConfigs, hrts_state.goalConf)
    push!(pathEdgeTypes, "D")
    while !isnothing(node.parent)
        push!(pathConfigs, node.config)
        push!(pathEdgeTypes, node.edgType)

        node=node.parent
    end
    push!(pathConfigs, node.config)
    push!(pathEdgeTypes, node.edgType)
    
    path, pathLength_unSmooth, intgrConsTotal = FormatPath(reverse(pathConfigs), reverse(pathEdgeTypes), hrts_state)
    # println("Path Length before smoothenign: ", pathLength_unSmooth )
    # println("Integer Cons before smoothenign: ", intgrConsTotal )

    # smoothPath = PathSmoothening(path, hrts_state)
    smoothPath = path
    println("Path smoothening skipped")
    # println("Path Length after smoothenign: ", ComputePathLength(smothPath) )
    # println("Integer Cons after smoothenign: ", ComputeIntgrCons(smothPath) )
    return smoothPath, ComputePathLength(smoothPath), pathLength_unSmooth, ComputeIntgrCons(smoothPath)
end

function FormatPath(pathConfigs, pathEdgeTypes, hrts_state)

  formatdPath = [(pathConfigs[1], pathEdgeTypes[1], 0., 0.)]
  pathLength = 0.
  intgrConsTotal = 0.
  for k in 2:length(pathConfigs)

    if pathEdgeTypes[k] == "D"
      _, dubPath = dubins_shortest_path(collect(pathConfigs[k-1]), collect(pathConfigs[k]), hrts_state.rho)
      # dubLength = dubins_path_length(dubPath)
      inf1, inf2, segLengths = InflexionPoints(dubPath)
      pathTypeStr = DubinsPathTypeString(dubPath)
      segEndPts = [inf1, inf2, pathConfigs[k]]
      segEndPts2 = [pathConfigs[k-1], inf1, inf2, pathConfigs[k]]
      for i in 1:3   
        if segLengths[i]>1e-6
          intgrConsSeg = EdgeIntgrCons(segEndPts2[i], segEndPts2[i+1], hrts_state, segLengths[i], pathTypeStr[i] )     
          push!(formatdPath, (segEndPts[i], pathTypeStr[i], segLengths[i], intgrConsSeg))
          intgrConsTotal += intgrConsSeg
        end
      
      end
      pathLength += sum(segLengths)
    else
      intgrConsSeg = EdgeIntgrCons(collect(pathConfigs[k-1]), collect(pathConfigs[k]), hrts_state, hrts_state.arcLen, pathEdgeTypes[k] )     

      push!(formatdPath, (pathConfigs[k], pathEdgeTypes[k], hrts_state.arcLen, intgrConsSeg))
      pathLength += hrts_state.arcLen
      intgrConsTotal += intgrConsSeg
    end

  end

  return formatdPath, pathLength, intgrConsTotal
end


function PlotHRTSPath(pathString, rho, pathFmt)
  # input is a vector of tuples, each tuple is a config
  # p = plot(aspect_ratio=:equal,legend=false)
  # arrLen = 50
  for j in 1:length(pathString)-1
      conf = pathString[j][1]
      if pathString[j+1][2] == "D"
        _, dubPath = dubins_shortest_path(collect(conf), collect(pathString[j+1]), rho)
        PlotDubins(dubPath, pathFmt)
      else
        PlotSegment(conf, pathString[j+1][3], pathString[j+1][2], rho, pathFmt)
      end      
  end
  # savefig(p,"path_plot.png")
end

function PathWayPoints(pathString, rho, deltaL)
  # input is a vector of tuples, each tuple is a config
  
  wyPtspath = Vector{Vector{Float64}}(undef, 0)
  push!(wyPtspath, collect(pathString[1][1][1:2]))
  for j in 1:length(pathString)-1
      conf = pathString[j][1]
      wpSeg = WayPointsSegment(conf, pathString[j+1][3], pathString[j+1][2], rho, deltaL)
      append!(wyPtspath, wpSeg[2:end])
  end
  return wyPtspath
end

function SegType(config1, config2, arcLen, rho)
        
    config1 = collect(config1)
    config2 = collect(config2)
    t2 = mod(config2[3],2*pi)
    dphi = arcLen/rho
    chordLen = 2*rho*sin(dphi/2)
    eucDisSeg = sqrt((config1[1]-config2[1])^2+(config1[2]-config2[2])^2)
    if abs(eucDisSeg-arcLen)<1e-6
        return "S"
    else
        c2_leftTurn = config1[1:2]+2*rho*sin(dphi/2)*[cos(dphi/2), sin(dphi/2)]
        c2_finalHead = mod(config1[3]+dphi, 2*pi)
        
        if norm(c2_leftTurn-config2[1:2])<1e-4 && abs(c2_finalHead-t2)<1e-4
            return "L"
        end
        c2_rightTurn = config1[1:2]+2*rho*sin(dphi/2)*[-cos(dphi/2), sin(dphi/2)]
        c2_finalHead = mod(config1[3]-dphi, 2*pi)       
        if norm(c2_rightTurn-config2[1:2])<1e-4 && abs(c2_finalHead-t2)<1e-4
            return "R"
        end
    end
    return "D"
end


function SumRiskFunctions(riskParams, posCoords, confAlpha=4.605 )
    
    riskVals = [RiskFunction(riskParams.riskCntrs[k], invCovMat, posCoords, confAlpha) for (k, invCovMat) in enumerate(riskParams.invCovMats)]
   
    if isempty(riskVals)
      return 0.
    else
      return sum(riskVals)
    end
end

function RiskFunction(center, invCovMat, posCoords, confAlpha = 4.605 )
    # confAlpha = 3.219 # for 80% confidence level
    # confAlpha = 4.605 # for 90% confidence level
    # confAlpha = 5.991 # for 95% confidence level
    # confAlpha = 9.210 # for 99% confidence level

  #    Risk as Gaussian distribution    
    # xmu = collect(posCoords)- collect(center)        
    # return exp(-0.5*transpose(xmu)*invCovMat*xmu)
    # if norm(xmu) < maxDist
    if PtLiesInConfEllipse(posCoords, center, invCovMat, confAlpha)
        xmu = collect(posCoords)- collect(center)
        return exp(-0.5*transpose(xmu)*invCovMat*xmu)
    else 
        return 0.0
    end
end


function EdgeIntgrCons(config1, config2, hrts_state, arcLen=nothing, segType=nothing)
        
    # config1 = [node1.coords[0], node1.coords[1], node1.heading]
    # config2 = [node2.coords[0], node2.coords[1], node2.heading]   
    if isnothing(segType)           
        segType = SegType(config1, config2, hrts_state.arcLen, hrts_state.rho)
    end
    if isnothing(arcLen)
        arcLen = hrts_state.arcLen
    end
    if arcLen <1e-6 return 0. end

    numSamples = max(3, Int(3*round(arcLen/hrts_state.arcLen)))
    arcLenVec = LinRange(0, arcLen, numSamples)                       
    arcAngleVec = arcLenVec./hrts_state.rho
    # println("arcAngleVec: ", arcAngleVec)
    cordLenVec = 2*hrts_state.rho*sin.(arcAngleVec./2)
    delL = hrts_state.rho*(arcAngleVec[2]-arcAngleVec[1])


    if segType == "D"
        _, dbPath  =dubins_shortest_path(collect(config1), collect(config2), hrts_state.rho)
        riskIntgrCons, _ = IntgrConsDubPath(dbPath, hrts_state)
    else
        if segType == "L"
            coordsSegVec = [(config1[1]+cordLenVec[i]*cos(config1[3]+arcAngleVec[i]/2), config1[2]+cordLenVec[i]*sin(config1[3]+arcAngleVec[i]/2)) for i in 1:numSamples]
        elseif segType == "R"
            coordsSegVec = [(config1[1]+cordLenVec[i]*cos(config1[3]-arcAngleVec[i]/2), config1[2]+cordLenVec[i]*sin(config1[3]-arcAngleVec[i]/2)) for i in 1:numSamples]
        elseif segType == "S"
            coordsSegVec = [(config1[1]+arcLenVec[i]*cos(config1[3]), config1[2]+arcLenVec[i]*sin(config1[3])) for i in 1:numSamples]
        end
        riskFnEvals = [SumRiskFunctions(hrts_state.riskParams, coordsSegVec[i]) for i in 1:numSamples]
        # println("coordsSegVec: ", coordsSegVec)        
        # println("riskFnEvals: ", riskFnEvals)
        riskIntgrCons = sum(riskFnEvals[2:end-1])*delL + (riskFnEvals[1]+riskFnEvals[end])*0.5*delL
        # PyPlot.scatter([e[1] for e in coordsSegVec], [e[2] for e in coordsSegVec])
    end
    return riskIntgrCons
end

function IntgrConsDubPath(dbPath, hrts_state)

    pathTypeStr = DubinsPathTypeString(dbPath)
    # segLengths = dbPath.params*dbPath.ρ
    startConf = dbPath.qi    
    endConf = dubins_path_endpoint(dbPath)[2]
    # infPt1 = dbPath.sample(dbPath.segment_length(0))
    # infPt2 = dbPath.sample(dbPath.segment_length(0)+dbPath.segment_length(1)-1e-8)
    infPt1, infPt2, segLengths = InflexionPoints(dbPath)
    configsVec = [startConf, infPt1, infPt2, endConf]
    riskDubSegs = Vector{Float64}(undef, 0)
    for k in 1:3
        if segLengths[k]>1e-6
          segRisk = EdgeIntgrCons(configsVec[k], configsVec[k+1], hrts_state, segLengths[k], pathTypeStr[k]) 
        # riskIntgrl += segRisk
          push!(riskDubSegs,segRisk)
        else
          push!(riskDubSegs,0.)
        end
    end
    return sum(riskDubSegs), riskDubSegs
end

function IntgConsToGoLB(pos, hrts_state)

    # distance from risk center to current position
    pos = collect(pos)[1:2]
    centers = hrts_state.riskParams.riskCntrs
    maxDist = hrts_state.riskParams.riskMaxDist        
    numDiscs = 10
    dVec = LinRange(0, maxDist, numDiscs)
    delD = maxDist/(numDiscs-1)
    riskLB = 0
    for (k, riskCntr) in enumerate(centers)
        dist_op = norm(pos-collect(riskCntr))
        if dist_op < maxDist
            invCovMat = hrts_state.riskParams.invCovMats[k]
            cov = inv(invCovMat)
            e, ev = eigen(cov)                
            ev2_theta = atan(ev[2,2], ev[1,2]) # second axis direction
            
            discPts1 = [(pos[1]+dist*cos(ev2_theta), pos[2]+dist*sin(ev2_theta)) for dist in dVec[1:end-1]]
            discPts2 = [(pos[1]+dist*cos(ev2_theta+pi), pos[2]+dist*sin(ev2_theta+pi)) for dist in dVec[1:end-1]]
            
            riskVec1 = [RiskFunction(riskCntr, hrts_state.riskParams.invCovMats[k], discPos) for discPos in discPts1]
            riskVec2 = [RiskFunction(riskCntr, hrts_state.riskParams.invCovMats[k], discPos) for discPos in discPts2]
            
            riskLB  += min(sum(riskVec1)*delD, sum(riskVec2)*delD)
        end
    end
            
    # riskLB=0
    return riskLB
end

function RemZeroSegs(path)
  # removing the zero segment lengths
  newPath = [path[1]]
  for p in path
    if p[3] > 1e-5
        push!(newPath, p)
    end
  end
  return newPath
end

function PathMergeSameSegs(path)
# Combining the consecutive same segment types
  newPath = [path[1]]
  for curSeg in path[2:end]
    prvSeg = newPath[end]
    if curSeg[2] == prvSeg[2]
        pop!(newPath)        
        push!(newPath, (curSeg[1], curSeg[2], curSeg[3]+prvSeg[3], curSeg[4]+prvSeg[4] ))
    else
        push!(newPath, curSeg)        
    end
  end
  return newPath
end

function ReverseSegment(segType)
    
    if segType == "L"
        return "R"
    elseif segType == "R"
        return "L"
    else
        return segType
    end
end
function ReversePath(path)
    
  revString = reverse(path)
  p = revString[1]
  revPath = [ ( (p[1][1], p[1][2], mod(p[1][3]+pi, 2*pi)), "Start", 0., 0.) ]
  
  for (ind, p) in enumerate(revString[2:end])
      p_prev = revString[ind]
      p_rev = ( (p[1][1], p[1][2], mod(p[1][3]+pi, 2*pi) ), ReverseSegment(p_prev[2]), p_prev[3], p_prev[4] )
      push!(revPath, p_rev)
  end
  return revPath
end

function PathSmoothening(hybPath, hrts_state, smoothResMax=6, nIter=1)

  # smoothRes; #number of segment considered at each iteration for smoothenign
  
  hybPath = RemZeroSegs(hybPath)
  hybPath = PathMergeSameSegs(hybPath)
          
  iter=1
  while iter<=nIter
      iter +=1
      for smoothRes in smoothResMax:-1:3

          for segNum in 1:smoothRes
              hybPath = PathSmoothSegFwd(hybPath, smoothRes, segNum, hrts_state)
              hybPathRev = ReversePath(hybPath)
              hybPathRev = PathSmoothSegFwd(hybPathRev, smoothRes, segNum, hrts_state)
              hybPath = ReversePath(hybPathRev)
          end
        end
  end
  return hybPath
end

function PathSmoothSegFwd(hybPath, smoothRes, segNum, hrts_state)

  pathLengthsVec = [p[3] for p in hybPath]
  pathIntgrVec = [p[4] for p in hybPath]
  pathIntgrCons_slack = hrts_state.riskParams.maxRisk-sum(pathIntgrVec)
                      
  numHybPathSegs = length(hybPath)-1 # the first entry in hybPath is start
  if numHybPathSegs < smoothRes
      return hybPath
  end
  smoothPath = [ hybPath[k] for k in 1:segNum ]
  while segNum+smoothRes <= numHybPathSegs+1
    frConf = hybPath[segNum][1]
    toConf = hybPath[segNum+smoothRes][1]
    
    _, dbPath = dubins_shortest_path(collect(frConf), collect(toConf), hrts_state.rho)
    dubPathLength = dubins_path_length(dbPath)
    if dubPathLength == 0
        print("err")
    end
    dbPathIntgrCons, riskDubSegs = IntgrConsDubPath(dbPath, hrts_state)
    segLength = sum(pathLengthsVec[segNum+1:segNum+smoothRes])
    segIntgrCons = sum(pathIntgrVec[segNum+1:segNum+smoothRes])
    pathIntgrCons_slack = pathIntgrCons_slack-dbPathIntgrCons+segIntgrCons

    if dubPathLength <= segLength && pathIntgrCons_slack>=0 && 
      CheckDubinsPathFeas2(dbPath, hrts_state.obstaclesPolyList, hrts_state.originCoordinates, hrts_state.xy_bnds )
        pathTypeStr = DubinsPathTypeString(dbPath) 
        _, _, dubSegLengths = InflexionPoints(dbPath)        
        dubSegLengthsSum = cumsum(dubSegLengths)

        for k in 1:3             
            if dubSegLengths[k] > 1e-10
                dubSegLenSum = dubSegLengthsSum[k]
                if abs(dubSegLenSum-dubPathLength)<1e-8
                    pt = toConf
                else
                    _, pt = dubins_path_sample(dbPath, dubSegLenSum)  
                end 
                
                edgIntgrCons = riskDubSegs[k]
                push!(smoothPath, ( (pt[1],pt[2],pt[3]), pathTypeStr[k], dubSegLengths[k], edgIntgrCons ) )
            end
        end
    else
        for k in 1:smoothRes                                                             
            push!(smoothPath,  hybPath[segNum+k] )
        end
    end
    segNum = segNum+smoothRes
    # segNum = segNum+1
  end
  if segNum < numHybPathSegs+1
      for ind in segNum+1:numHybPathSegs+1
          push!(smoothPath, hybPath[ind] )
      end
  end

  return smoothPath
end

function ComputePathLength(path)
  lengthsVec = [l[3] for l in path]
  return sum(lengthsVec)
end

function ComputeIntgrCons(path)
  intgrConsVec = [l[4] for l in path]
  return sum(intgrConsVec)
end

function PlotResultScenario(hrtsResult, ax, cmp=get_cmap("Reds"))

  riskGaussParams = hrtsResult.htState.riskParams
  xy_bnds = hrtsResult.htState.xy_bnds
  startConf = hrtsResult.htState.startConf
  goalConf = hrtsResult.htState.goalConf
  rho = hrtsResult.htState.rho
  obsList = hrtsResult.htState.obstaclesPolyList

  PlotObstacles(obsList)
  PlotRiskSumContour(riskGaussParams, xy_bnds, ax, cmp)
  # pathFmt = (color = "b", linewidth=2, linestyle = "-", marker = "x",  endPoints = true)    
  pathFmt = (color = rand(3), linewidth=2, linestyle = "-", marker = "x",  endPoints = true)    

  PyPlot.scatter([startConf[1]], [startConf[2]], color = "green", marker="^", label="Start")
  PyPlot.scatter([goalConf[1]], [goalConf[2]], color = "orange", marker="^", label="Goal")
  
  # if hrtsResult.status == :success    
  #     PlotHRTSPath(hrtsResult.path, rho, pathFmt)   
  #     PyPlot.axis([-100, xy_bnds[1], -100, xy_bnds[2]]);
  # end
end
