import Base
using DataStructures

"""Node of the state tree to explore"""
mutable struct NodeG{}
  x::Int64
  y::Int64
  f::Float64
  parent::UInt64
end
"order by f = g + h"
Base.isless(n1::NodeG, n2::NodeG) = Base.isless(n1.f, n2.f)

function get_motion_model()
    # dx, dy, cost
    motion = [[1, 0, 1],
              [0, 1, 1],
              [-1, 0, 1],
              [0, -1, 1],
              [-1, -1, sqrt(2)],
              [-1, 1, sqrt(2)],
              [1, -1, sqrt(2)],
              [1, 1, sqrt(2)]]

    return motion
end

function GetNodeId(x_ind, y_ind, xIndsMax=10000)
    # return x_ind * xIndsMax + y_ind
    return hash((x_ind, y_ind))
end

    
function VerifyNode(node, obsMap, xyIndsMax)
    if node.x < 0
        return false
    elseif node.y < 0
        return false
    elseif node.x > xyIndsMax[1]
        return false
    elseif node.y > xyIndsMax[2]
        return false
    end

    # if obsMap[node.x, node.y]:
    if GetNodeId(node.x, node.y) in obsMap
        return false
    end
    return true
    
end

function PrintHeurCost(heurCostMap, goalInd, xyIndsMax, cellRes)
    PyPlot.axis("equal")    
    PyPlot.scatter([goalInd[1]*cellRes], [goalInd[2]*cellRes], marker= "^", color="r" )
    
    for i in 0:xyIndsMax[1]
        for j in 0:xyIndsMax[2]
            ind = GetNodeId(i, j)
            if ind in keys(heurCostMap)
                PyPlot.scatter([i*cellRes],[j*cellRes], marker=".", color="m")
                PyPlot.text(i*cellRes,j*cellRes, string(Int(floor(heurCostMap[ind].f))), fontsize=8)
                # println("Heuristic cost of inds: ",i,", ", j,": ", heurCostMap[ind].f)
            end
        end
    end
end
                
function Calc_HeurCost(goalInds, obsIndsHash, xyIndsMax, cellRes)
    
    # gx: goal x index
    # gy: goal y index
    # ox: x position list of Obstacles indices
    # oy: y position list of Obstacles indices
    # resolution: grid resolution [m]
    # The algorithm follows similar to dynamic program by incremental search starting from the goal node
    

    goal_node = NodeG(goalInds[1], goalInds[2], 0.0, 0)

    motion = get_motion_model()

    # open_set, closed_set = Dict(), Dict()
    # closed_set = Set{UInt64}()
    open_set = Dict(GetNodeId(goal_node.x,goal_node.y) => goal_node)
    closed_set = Dict(GetNodeId(goal_node.x,goal_node.y) => goal_node)
    
    # priority_queue = [(0, GetNodeId(goal_node.x,goal_node.y))]  
    openheap = [goal_node]  
    while !isempty(openheap)

        curNode = heappop!(openheap)
        # println("length of open heap: ", length(openheap))

        c_id =  GetNodeId(curNode.x,curNode.y)
        push!(closed_set, c_id =>curNode )
        delete!(open_set, c_id)

        for m in motion
            node = NodeG(curNode.x + m[1],
                        curNode.y + m[2],
                        curNode.f + cellRes*m[3], c_id)
            n_id = GetNodeId(node.x, node.y)

            if n_id in keys(closed_set)
                continue
            end
            if !VerifyNode(node, obsIndsHash, xyIndsMax)               
                continue
            end
            if !(n_id in keys(open_set))
                open_set[n_id] = node  # Discover a new node
                heappush!(openheap, node)
            else
                if open_set[n_id].f >= node.f
                    # This path is the best until now. record it!                    
                    deleteat!(openheap, findall(x->x==open_set[n_id], openheap))
                    open_set[n_id] = node
                    heappush!(openheap, node)
                end
            end
        end
    end
    # PrintHeurCost(closed_set, goalInds, xyIndsMax, cellRes) #prints the costs in plot
    # println("length of heuristic cost dict: ", length(closed_set))
    return closed_set
    
end