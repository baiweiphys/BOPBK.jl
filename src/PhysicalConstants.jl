"""
Description: PA module containing fundamental physical constants in SI units..
Author: Bai Wei (baiweiphys@gmail.com)
Date: 2025-12-04
LastEditors: Bai Wei
Last Modified: 2026-02-03
"""

# PhysicalConstants.jl
module PhysicalConstants

export PlasmaConfig

"""
    PlasmaConfig
A thread-safe, type-stable container for physical constants.
"""
struct PlasmaConfig
    c2::Float64
    epsilon0::Float64
    mu0::Float64
    kB::Float64
    qe::Float64
    me::Float64
    mp::Float64

    # Internal constructor to ensure consistency
    function PlasmaConfig()
        # Fundamental constants
        _c2 = (2.99792458e8)^2 # Speed of light squared [m²/s²]
        # Electromagnetic constants
        _eps0 = 8.854187817e-12  # Vacuum permittivity [F/m]
        _mu0 = 1.0 / (_c2 * _eps0) # Vacuum permeability [N/A²]
        # Thermodynamic constants
        _kB = 1.38064852e-23 # Boltzmann constant [J/K]
        # Elementary particle constants
        _qe = 1.602176634e-19 # Elementary charge [C]
        _mp = 1.67262192369e-27 # Proton mass [kg]
        _me = 9.1093837015e-31 # Electron mass [kg] 
        return new(_c2, _eps0, _mu0, _kB, _qe, _me, _mp)
    end
end

# Define a constant global instance for easy access
const SI = PlasmaConfig()

end # module
