ADD_LIBRARY(Material Material.cpp LinearElasticMaterial.cpp 
            NeoHookean.cpp MooneyRivlin.cpp
            IsotropicLinearElasticJ2PlasticMaterial.cpp
            IsotropicLinearElasticJ2PlasticPlaneStressMaterial.cpp
            KorkolisKyriakidesPlaneStressMaterial.cpp
            KorkolisKyriakidesPlaneStressMaterialWithExperimentalYielding.cpp )
#add_library(Material Material.cpp LinearElasticMaterial.cpp 
#                   NeoHookean.cpp  )
