# Indentation has to be consistent
# We are using two spaces as indentation.

# Creating a position class
class Position(object):
  # Special method __init__ to initialize your data 
  # attributes
  # Note as opposed to normal function, __init__ contains
  # self;
  # self: parameter to refer to an instance of a class
  # x,y: what you provide while creating this class/calling it
  def __init__(self, x, y):
    # self.x or self.y: Look for x/y that belong to this class
    self.x = x
    self.y = y

  # Methods that you can use to compute things
  # Following method computes distance of (x,y) from origin
  def get_dist_from_origin(self):
    # Note how the x/y values are called (with self.)
    dist = (self.x**2.0 + self.y**2.0)**0.5
    return dist

  # You can add as-many methods as you want
