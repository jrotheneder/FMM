namespace fmm {

template<typename Vector, typename Source, std::size_t d, typename = void>
struct LocalExpansion {
    static_assert(d==2 || d==3, 
        "This implementation supports only 2 or 3 dimensions.\n"
    ); 
};

// 2-D Implementation
template<typename Vector, typename Source, std::size_t d>
struct LocalExpansion<Vector, Source, d, typename std::enable_if<d==2>::type> {

    using ME = MultipoleExpansion<Vector, Source, d>;

    std::vector<double> coefficients; 

    LocalExpansion(ME me = {}) {} // TODO 
    std::vector<double> shift(Vector v);
};

// 3-D Implementation
template<typename Vector, typename Source, std::size_t d >
struct LocalExpansion<Vector, Source, d, typename std::enable_if<d==3>::type> {

    using ME = MultipoleExpansion<Vector, Source, d>;

    std::vector<double> coefficients; 

    LocalExpansion(ME me = {}) {} // TODO 
    std::vector<double> shift(Vector v);
};

} // namespace fmm
