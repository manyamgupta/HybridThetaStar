# Hybrid Theta*
This is a motion planning algorithm for vehicles with curvature constraints, such as minimum turn radius, a primary cost such as path length, and  a secondary cost such as resource cost. 
Inputs to the algorithm are: 
map with obstacles, risk field (resource cost), start and end configurations, minimum turn radius, resource cost limit
Output: minimum cost path satisfying curvature constraints and integral resource constraint

For further details refer to the paper:
S. G. Manyam, D. W. Casbeer and C. Taylor, "Hybrid Theta∗: Motion Planning for Dubins Vehicles With Integral Constraints," in IEEE Robotics and Automation Letters, vol. 10, no. 2, pp. 1497-1504, Feb. 2025, doi: 10.1109/LRA.2024.3522786.

Distribution Statement A: Approved for Public Release, PA# AFRL-2024-6456
