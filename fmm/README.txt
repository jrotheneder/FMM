Brief overview of the content of the various source files: 

abstract_orthtree.hpp
    Abstract base class for quad-/octrees 
adaptive_orthtree.hpp
    Abstract base class for adaptive quad-/octrees
balanced_orthtree.hpp
    Abstract base class for balanced quad-/octrees

abstract_fmm_tree.hpp
    Abstract base class for the FMM tree data structure 
adaptive_fmm_tree.hpp
    Implementation of the adaptive FMM
balanced_fmm_tree.hpp
    Implementation of the balanced FMM

series_expansion.hpp
    Abstract base class for series expansions 
multipole_expansion.hpp
    Implementation of the multipole expansion
local_expansion.hpp
    Implementation of the local expansion

vector.hpp
    Vector and pointsource classes used in the FMM
fields.hpp
    Potential and force field functions, direct algorithm 
fmm_tables.hpp
    Classes that function as lookup tables which help reduce the number of 
    expensive function evaluations
fmm_general.hpp
    Several functions used in multiple locations but fall into none of the
    other categories
debugging.hpp
    Code that aided in debugging, printing of output, writing to file etc


test
    Folder containing some of the code used for testing
examples
    Folder containing some minimal examples
misc
    Folder containing things that belong nowhere else
mathematica
    Folder containing the code for the Mathematica interface
logs
    Target directory for logs
