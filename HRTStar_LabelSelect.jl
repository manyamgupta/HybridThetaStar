function HRTStarLabelSelect(    
  start,
  goal,
  rho, 
  arcLen, 
  obsList,
  riskGaussParams,
  originCoordinates,
  xy_bnds,
  cellRes,
  hdngRes,
  pltExplFlag,
  ) 
  timeout = Inf
  
  arcAngle, cordLen = NeighborParams(rho, arcLen)
  obstaclesPolyList = ObsPolygonList(obsList)
  infeas_xyinds_hashlist, obs_snap_coordsList = GenerateObstacleInds(xy_bnds,cellRes, obstaclesPolyList, originCoordinates )
  xy_inds_max = (Int(floor(((xy_bnds[1]-originCoordinates[1])/cellRes))), Int(floor((xy_bnds[2]-originCoordinates[2])/cellRes)) )
  goalInds2D = PointToIndex(goal, originCoordinates, cellRes)

  heurCostDict = Calc_HeurCost(goalInds2D, infeas_xyinds_hashlist, xy_inds_max, cellRes)
#   PyPlot.figure()
#   PyPlot.scatter(obs_snap_coordsList[:,1], obs_snap_coordsList[:,2], marker=".", color="blue")
#   PlotObstacles(obsList)
#   return HRTSResult(:infeasible, [(start, "N", 0.)], NaN, NaN, NaN, nothing)
  
  if !CheckStartGoalFeasible(start, goal, originCoordinates,xy_bnds, cellRes, infeas_xyinds_hashlist )
    return HRTSResult(:infeasible, [(start, "N", 0., 0.)], NaN, NaN, NaN, nothing)
  end
  startHeur = Heuristic(start[1:2], heurCostDict, originCoordinates, cellRes)
  if isnothing(startHeur)
    return HRTSResult(:infeasible, [(start, "N", 0., 0.)], NaN, NaN, NaN, nothing)
  end
  hrts_state = InitializeSearchState(start, goal, obstaclesPolyList, riskGaussParams, infeas_xyinds_hashlist, originCoordinates, xy_bnds, xy_inds_max, cellRes, hdngRes, rho, arcLen, startHeur)
  
  # return HRTSResult(:infeasible, [(start, "N", 0.)], NaN, NaN, nothing)
  
  if pltExplFlag
    # PyPlot.figure()
    fig, ax = PyPlot.subplots(1,1)
    PlotObstacles(obsList)
    PlotRiskSumContour(riskGaussParams, xy_bnds, ax)
    PyPlot.scatter([start[1]], [start[2]], color = "green", marker="^")
    PyPlot.scatter([goal[1]], [goal[2]], color = "green", marker="^")
  end
  
  while !isempty(hrts_state.openheap)
    curNode = heappop!(hrts_state.openheap)
    if ConfToIndex(curNode.config, hrts_state) == (5, 14, 42)
      println("curNode: ", curNode.config, ", g: ", curNode.g, ", f: ", curNode.f, "intgrConf: $(curNode.intgrCons) ", "intgrConsTotalLB: $(curNode.intgrConsTotalLB)")    
    end
    # if norm(collect(curNode.config[1:2])-[1207., 3099.]) < 10
    #   println("curNode: ", curNode.config, ", g: ", curNode.g, ", f: ", curNode.f, "intgrConf: $(curNode.intgrCons) ", "intgrConsTotalLB: $(curNode.intgrConsTotalLB)")
    # end
    # println("curNode: ", curNode.config, ", g: ", curNode.g, ", f: ", curNode.f)
    if IsFeasibleAnalyticExpn(curNode, goal, hrts_state)
    # if feasAE
      # println("Solution found")
      # if pltExplFlag   savefig("explore.png") end
      htsPath, pathLength, pathLength_unSmooth, intgrConsTotal =  Reconstructpath(curNode, hrts_state)
      return HRTSResult(:success, htsPath, pathLength, pathLength_unSmooth, intgrConsTotal, hrts_state)
      
    end

    # nodehash = hash(curNode.config)
    nodeSnapHash = hash(ConfToIndex(curNode.config, hrts_state))
    # delete!(hrts_state.nodedict, nodeSnapHash)

    if timeout < Inf && time() - hrts_state.startTime >= timeout
      # if pltExplFlag  savefig("explore.png") end
      return HRTSResult(:timeout, [(start, "N", 0., 0.)], NaN, NaN, NaN, hrts_state)
    end
    push!(hrts_state.closedset, nodeSnapHash)

    neighbor_states, edgeTypesList = Neighbors(curNode.config, hrts_state, arcAngle, cordLen)
    # println("neighbour states: $neighbor_states")

    for (k, nbrState) in enumerate(neighbor_states)
      # println("neighbor state that is being treated: $nbrState")
      nbrSnapHash = hash(ConfToIndex(nbrState, hrts_state))
      # if ConfToIndex(nbrState, hrts_state) == (5, 14, 42)
      #   println("curNode: ", curNode.config, ", g: ", curNode.g, ", f: ", curNode.f, "intgrConf: $(curNode.intgrCons) ", "intgrConsTotalLB: $(curNode.intgrConsTotalLB)")
      # end
      nbrNodeId = GetNodeId4mConf(nbrState)
      #the conf of the node is snapped to the discretized grid, and if it is in closed list, then no further exploration of that node is done
      if nbrSnapHash in hrts_state.closedset 
        continue
      end
      # newSegIntgrCons = EdgeIntgrCons(curNode.config, nbrState, hrts_state)
      # nbrIntgrConsTotalLB = IntgConsToGoLB(nbrState[1:2], hrts_state)
      nbrNode = HRTSNode(nbrNodeId, nbrState, curNode)
      nbrNode.edgType = edgeTypesList[k]
      nbrNode.snapHash2D = hash(PointToIndex(nbrState[1:2], originCoordinates, cellRes))
      nbrNode.snapHash3D = nbrSnapHash
      push!(hrts_state.nodedict, nbrNodeId => nbrNode)

      # HRTSNode(nbrNodeId, nbrState, curNode, edgeTypesList[k], nbrSnapHash2d, nbrSnapHash3d, gfromCurNode,startHeur,0., 0.)
      nbrNode = CheckGpToNodeEdge(nbrNode, curNode, hrts_state, heurCostDict)
      nbrNode.intgrConsTotalLB = nbrNode.intgrCons + IntgConsToGoLB(nbrNode.config[1:2], hrts_state) 
      nbrNodeExistsInOL, exstngNodes = NodeExistsInList(nbrNode, hrts_state.openheap)

      newChildNodeReddnt = false
      if nbrNodeExistsInOL
        for nbrNodeInOL in exstngNodes
          
          if nbrNode.g < nbrNodeInOL.g && nbrNode.intgrCons < nbrNodeInOL.intgrCons 
            
            deleteat!(hrts_state.openheap, findall(x->x==nbrNodeInOL, hrts_state.openheap))        
            heapify!(hrts_state.openheap)
          elseif nbrNode.g >= nbrNodeInOL.g && nbrNode.intgrCons >= nbrNodeInOL.intgrCons      
          # elseif nbrNode.g >= nbrNodeInOL.g
            newChildNodeReddnt = true            
            # println("new child node $(nbrNode.config) found to be dominated")
            # println("*****************************************************")
            break
          end          
        end

      end      
      if !newChildNodeReddnt
        AddNodeToQue(nbrNode, hrts_state, heurCostDict)
      end

      if pltExplFlag   
        # curConf = curNode.config
        pathFmt = (color = "b", linewidth=1, linestyle = "-", marker = "x",  endPoints = true)
        # PlotSegment(curConf, arcLen, edgeTypesList[k], rho, pathFmt)
        PlotBranch(nbrNode.parent.config, nbrNode.config, nbrNode.edgType, arcLen, rho, pathFmt)
        
        PyPlot.axis("equal")
      end
    end
    # if pltExplFlag PyPlot.pause(0.001) end
  end
  return HRTSResult(:infeasible, [(start, "N", 0., 0.)], NaN, NaN, NaN, hrts_state)
end

function AddNodeToQue(newNode, hrts_state, heurCostDict)
  # This function adds the f cost to the newNode (g cost + heuristic), and the intgrConsTotalLB
  # newNode.intgrConsTotalLB = newNode.intgrCons + IntgConsToGoLB(newNode.config[1:2], hrts_state) 
  if newNode.intgrConsTotalLB < hrts_state.riskParams.maxRisk      
    # neighbourheuristic = Heuristic(newNode.config, heurCostDict, hrts_state.originCoordinates, hrts_state.cellRes)   
    neighbourheuristic = HeuristicDubins(newNode.config, hrts_state.goalConf, hrts_state.rho)        
    newNode.f = newNode.g + neighbourheuristic    
    heappush!(hrts_state.openheap, newNode)
    # hrts_state.nodedict[nbrSnapHash] =  newNode
  end
end
function NodeExistsInList(testNode, nodeList)
    # returns all the nodes in nodeList that matches with the test node snapHash3D,
    # there could be multiple nodes with same index id
    returnNodeList = [node for node in nodeList if node.snapHash3D == testNode.snapHash3D]

    return !isempty(returnNodeList), returnNodeList
end

function PlotBranch(startConf, endConf, edgType, arcLen, rho, pathFmt)

    if edgType =="D"
        _, dubPath = dubins_shortest_path(collect(startConf), collect(endConf), rho)
        PlotDubins(dubPath, pathFmt)
    else
        PlotSegment(startConf, arcLen, edgType, rho, pathFmt)
    end
end

function CheckGpToNodeEdge(newNode, currentNode, hrts_state, heurCostDict)
        
    grandParNode = currentNode.parent
    newNodeAttrAdded = false
    newSegIntgrCons = EdgeIntgrCons(currentNode.config, newNode.config, hrts_state)
    newNode.intgrCons = currentNode.intgrCons+newSegIntgrCons
    if !isnothing(grandParNode)           
        _, dubPath = dubins_shortest_path(collect(grandParNode.config), collect(newNode.config), hrts_state.rho)
        dubFeas = CheckDubinsPathFeas(dubPath, hrts_state)
        if dubFeas
        dubEdgIntgrCons, _ = IntgrConsDubPath(dubPath, hrts_state)        
    #             # if self.CheckDubPathFeas(dubPath) and grandParNode.intgrCons+dubEdgIntgrCons <= childNewIntgrCons:           
            if grandParNode.intgrCons+dubEdgIntgrCons <= newNode.intgrCons
                newNode.g = grandParNode.g + dubins_path_length(dubPath)
                newNode.parent = grandParNode
                newNode.edgType = "D"
                newNode.intgrCons = grandParNode.intgrCons+dubEdgIntgrCons        
                newNodeAttrAdded = true
            # else
            #     AddReplicaNodeFrmGP(newNode, grandParNode, grandParNode.g+dubins_path_length(dubPath), grandParNode.intgrCons+dubEdgIntgrCons, hrts_state, heurCostDict)
            end
        end
    end 
    if !newNodeAttrAdded      
        newNode.g = currentNode.g + hrts_state.arcLen #newNodeCost
        newNode.parent = currentNode        
    end
    # newNode.intgrConsTotalLB = newNode.intgrCons+IntgConsToGoLB(newNode.config[1:2], hrts_state)        
    # newNode.f = newNode.g + Heuristic(nbrNode.config, hrts_state.heurCostDict, hrts_state.originCoordinates, hrts_state.cellRes)
    return newNode
end

function AddReplicaNodeFrmGP(newNode, grandParNode, childCostFrmStart, childIntgrCons, hrts_state, heurCostDict)
        
  repNodeId = GetNodeIdWithExistngIndx(newNode.id, newNode.snapHash3D, hrts_state)

  repNode = HRTSNode(repNodeId, newNode.config, grandParNode)
  
  repNode.edgType = "D"
  repNode.snapHash2D = hash(PointToIndex(newNode.config[1:2], hrts_state.originCoordinates, hrts_state.cellRes))
  repNode.snapHash3D = hash(ConfToIndex(newNode.config, hrts_state))
  push!(hrts_state.nodedict, repNode.id => repNode)

  repNode.intgrCons = childIntgrCons
  repNode.g = childCostFrmStart
  repNode.intgrConsTotalLB = childIntgrCons + IntgConsToGoLB(newNode.config[1:2], hrts_state)
  
  AddNodeToQue(repNode, hrts_state, heurCostDict)

end

function GetNodeIdWithExistngIndx(nbrNodeId, newSnapHash3D, hrts_state)
        
  numExistingNodes = length([id for id in keys(hrts_state.nodedict) if hrts_state.nodedict[id].snapHash3D == newSnapHash3D])
  newNodeId = nbrNodeId+numExistingNodes

  return newNodeId
end

function FocalPop(opnheap, w)

  focalHeap = []
  minLabelCost = first(opnheap).f
  inds = findall(x->x.f< minLabelCost*w, opnheap )

  for nd in opnheap[inds]
    # push!(focalHeap, (nd.intgrConsTotalLB, nd.id))
    push!(focalHeap, (nd.intgrCons, nd.id))

  end
  returnNodeId = first(focalHeap)[2]
  popNode = opnheap[findall(x->x.id == returnNodeId, opnheap )]
  deleteat!(opnheap, findall(x->x.id == returnNodeId, opnheap ))

  return popNode[1]
end