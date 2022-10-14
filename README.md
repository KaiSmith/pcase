pcase is a SageMath module built for computing and studying Lagrangians in an orbifold known as the pillowcase. By work initiated by Hedden, Herald, and Kirk, Lagrangians in the pillowcase can be associated to 2-stranded tangles in the 3-ball and this can be used to study knots.
pcase can be used to compute upper bounds on the rank of Kronheimer and Mrowka's singular instanton homology. It is conjectured that there is a way to assign bounding cochains to these Lagrangians so that the resulting homology is isomorphic to singular instanton homology. The author plans to use this program to investigate possible assignments. pcaseLag is also being used to compute Chern-Simons invariants of knots in [Ronn] and investigate SU(2)-simplicity [Smith]. This program was written to perform computations for [Smith].


### What the program can do:

+ Generate the perturbed traceless SU(2) character varieties of tangles associated to Montesinos knots and their images in the pillowcase. Uses results from [Smith].

+ Given a character variety, it can modify it in the way described in [CHKK], which is associated to adding an earring to the tangle.

+ Computes the Lagrangian Floer chain groups of a pair of Lagrangians in the pillowcase. Implements the method in [HHK18, Section 5] to endow this group with a Z/4 grading.

+ Computes Lagrangian Floer homology in the pillowcase. Does this by implementing the algorithm from [Blank] to find immersed bigons.

### Future things:

+ Optimize by implementing a sweep-line algorithm to speed up finding intersections, which is the current bottle-neck.

+ Generate character varieties for tangles associated to arbirtary arborescent knots (in progress). Can be done by implementing more results from [Smith].

+ Make it easier to input PL approximations of other character varieties (for example, the character varieties associated to torus knots in [HHK18]).

+ Allow for computation of homology with bounding cochains.

### Usage

+ The file Example.ipynb is a file which has commented examples of the various functions available in pcase.

+ If there is something that is not clear, feel free to contact me: kjs7 (at) iu (dot) edu

### References:

[Blank] Samuel Blank: Extending immersions and regular homotopies in codimension 1. PhD Thesis, 1967.

[CHKK] Guillem Cazassus, Christopher Herald, Paul Kirk, Artem Kotelskiy: The correspondence induced on the pillowcase by the earring tangle. arXiv:2010.04320.

[HHK18] Matthew Hedden, Christopher Herald, Paul Kirk: The pillowcase and traceless representations of knot groups II: a Lagrangian-Floer theory in the pillowcase.  J. Symplectic Geom. 16 (2018), no. 3, 721–815. 

[Ronn] Mark Ronnenberg: The Chern-Simons bundle and traceless character varieties of punctured surfaces and tangles. In preparation.

[Smith] Kai Smith: The perturbed traceless SU(2) character varieties of tangle sums. In preparation.