# Import the class Position
from position import Position

pos = Position(2,3)
# Access data from class
print("x=%.6f"%(pos.x))
print("y=%.6f"%(pos.y))
# Compute distance by calling the function from class
distance = pos.get_dist_from_origin()
print("distance from origin(0,0)=",distance)
