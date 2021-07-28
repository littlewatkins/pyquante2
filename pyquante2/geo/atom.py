"""
Class to create an atom object

>>> h = atom(1,0,0,0)
>>> h
1 H     0.000000     0.000000     0.000000
>>> h.r
array([ 0.,  0.,  0.])

"""
import numpy as np
from pyquante2.geo.elements import floatcolor,radius,mass,symbol
from pyquante2.constants import ang2bohr
from pyquante2.utils import norm2

class atom(object):
    """creates an atom that can be used in calculations
    """
    def __init__(self,atno,x,y,z,**kwargs):
        """initializes the atom with the atomic number and its coordinates in 
        cartesian space

        Parameters
        ----------
        atno : int
            atomic number of the atom
        x : int, float
            x cartesian coordinate
        y : int, float
            y cartesian coordinate
        z : int, float
            z cartesian coordinate
        **kwargs : str
            used for passing in units??
        """
        self.atno = int(round(atno))
        self.Z = atno
        self.r = np.array([x,y,z],'d')
        self.units = kwargs.get('units','bohr').lower()
        assert self.units[:4] in ['bohr','angs']
        if not self.units == 'bohr':
            self.r *= ang2bohr
        return

    def atuple(self):
        """Creates a tuple out of the atomic number, and the cartesian coordinates

        Returns
        -------
        tuple
        """
        return (self.atno,self.r[0],self.r[1],self.r[2])

    def html_row(self,table,i=0):
        """Add a row in a table for an html representation of an atom

        Returns
        -------
        None
        """
        import xml.etree.ElementTree as ET
        tr = ET.SubElement(table,"tr")
        td = ET.SubElement(tr,"td")
        td.text = str(i+1)
        td = ET.SubElement(tr,"td")
        td.text = str(self.atno)
        td = ET.SubElement(tr,"td")
        td.text = symbol[self.atno]
        td = ET.SubElement(tr,"td")
        td.text = "%.5f" % self.r[0]
        td = ET.SubElement(tr,"td")
        td.text = "%.5f" % self.r[1]
        td = ET.SubElement(tr,"td")
        td.text = "%.5f" % self.r[2]
        return

    def __repr__(self):
        """returns the atomic number, symbol, and cartesian coordinates

        Returns
        -------
        str
        """
        return "%d %s %12.6f %12.6f %12.6f" % (
                    self.atno,symbol[self.atno],self.r[0],self.r[1],self.r[2])
    
    def __getitem__(self, i):
        """returns the specific cartesian coordinate for x=0, y=1, and z=2

        Returns
        -------
        int, float
        """
        return self.r[i]

    def xyz(self):
        """compiles a string of the atomic number and the cartesian coordinates

        Returns
        -------
        str
            atomic number and cartesian coordinates
        """
        return "%4s %12.6f %12.6f %12.6f" %\
               (symbol[self.atno],self.r[0],self.r[1],self.r[2])

    def distance(self,other):
        """computes the norm of the difference between the atom's specified
        location and another location in cartesian coordinates

        Returns
        -------
        float
        """
        return np.linalg.norm(self.r-other.r)

    def color(self):
        """returns the color associated with the atomic number in elements.py

        Returns
        -------
        str
        """
        return floatcolor[self.atno]

    def radius(self):
        """returns the radius associated with the atomic number in elements.py
        """
        return radius[self.atno]

    def mass(self):
        """returns the mass associated with the atomic number in elements.py
        
        Returns
        -------
        float
        """
        return mass[self.atno]

    def symbol(self):
        """returns the symbol associated with the atomic number in elements.py

        Returns
        -------
        str
        """
        return symbol[self.atno]

    def tag(self,index=1):
        """\
        >>> h = atom(1,0,0,0)
        >>> h.tag()
        'H1'
        >>> li = atom(3,0,0,0)
        >>> li.tag(4)
        'Li4'
        """
        return "%s%d" % (self.symbol(),index)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
