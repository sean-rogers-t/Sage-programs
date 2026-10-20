"""Generic-degree certificates for four normalized exceptional Springer maps.

The proof and its published intersection-theory inputs are documented in
EXCEPTIONAL_DEGREES.md. These degrees are for simply connected groups.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class SpringerDegreeCertificate:
    cartan_type: str
    representation: str
    highest_weight: tuple
    torus_lattice: str
    normalized_volume: int
    boundary_contribution: int
    central_cover: int
    method: str = 'toric intersection with explicit boundary Segre class'
    proof_reference: str = 'EXCEPTIONAL_DEGREES.md'

    @property
    def torus_degree(self):
        """Degree on the torus with the character lattice recorded above."""
        return self.normalized_volume-self.boundary_contribution

    @property
    def degree(self):
        """Degree of the Springer map on the simply connected group."""
        return self.central_cover*self.torus_degree


def exceptional_springer_degree(cartan_type, *, representation='minimal'):
    """Return the proved E6/E7 minimal or adjoint generic-degree calculation.

    For adjoint representations the volume and boundary calculation use Q;
    central_cover restores the simply connected group's degree. For minimal
    representations they use P. This is not a minimum Cayley-degree API.
    """
    if cartan_type not in ('E6','E7'):
        raise ValueError("degree certificates are available for 'E6' and 'E7'")
    if representation not in ('minimal','adjoint'):
        raise ValueError("representation must be 'minimal' or 'adjoint'")
    if cartan_type=='E6':
        if representation=='minimal':
            values=((1,0,0,0,0,0),'P',72*3+27*32,0,1)
        else:
            curves,triples,doubles=54*15,54*10,270*3
            boundary=7*curves-4*triples-2*doubles
            values=((0,1,0,0,0,0),'Q',54*52,boundary,3)
    elif representation=='minimal':
        values=((0,0,0,0,0,0,1),'P',576*3+126*64,126*20,1)
    else:
        values=((1,0,0,0,0,0,0),'Q',576*2+56*360,56*(240//3),2)
    return SpringerDegreeCertificate(cartan_type,representation,*values)
