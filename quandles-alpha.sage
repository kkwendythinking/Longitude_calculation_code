def spp(x,y):
    return Matrix([x,y]).det()

FIXED_VECTOR = vector([-1,0])

def quandleElementFromMatrix(mat, m):
    if mat[0,1].is_zero():
        alpha = 0
        z = mat[1,0]/(m-1/m)
        return quandleElement(z, alpha, m)
    else:
        alpha = mat[0,1]
        z = (mat[1,1] - m)/(-alpha)
        return quandleElement(z, alpha, m)

class quandleElement():

    """An element of the generic quandle"""

    def __init__(self, z, alpha, m):
        """Initialize with parameters

        :z: distinguished eigenvector is [z, 1]
        :alpha: other eigenvector is (proportional to) alpha v + (m - 1/m) FIXED_VECTOR
        :m: eigenvalue

        """
        self.z = z
        self.alpha = alpha
        self.m = m
        self.v = vector([z,1])

        #get matrix
        self.matrix = Matrix([
            [ alpha*z + 1/m , alpha ],
            [ -alpha*z^2 + m*z - z/m, -alpha*z + m],
        ])

    def vectorAction(self, y, left=False, inverse=False):
        """Return yX where X is the underlying matrix

        :y: vector
        :left: whether it's the right action on row vectors (default) or the left action on column vectors
        :inverse: whether it's the action (default) or inverse action
        :returns: yX

        """
        if left:
            raise NotImplementedError

        if inverse:
            return y*self.matrix
        else:
            return y*self.matrix^(-1)

    def quandleAction(self, other, inverse=False):
        """Act on another quandle element and return the result

        :other: the thing to act on
        :inverse: if true, return the inverse action
        :returns: a quandleElement

        """
        if not inverse:
            out = self.matrix^(-1)*other.matrix*self.matrix
        else:
            out = self.matrix*other.matrix*self.matrix^(-1)

        return quandleElementFromMatrix(out, other.m)

    def quandleUnder(self, other, inverse=False):
        return other.quandleAction(self, inverse)

    def getEquations(self, other):
        """Return polynomials whose zero sets assert that self and other are equal

        :other: another quandle elements
        :returns: a list of polynomials (in the entries of the quandle elements)

        Maybe just do rational functions for now
        """
        return [ self.z - other.z, self.alpha - other.alpha, self.m-other.m]
