from conans import ConanFile, CMake

class HydroConan(ConanFile) :
    settings = "os", "compiler", "build_type", "arch"
    requires = "OpenMPI/3.0.0@demo/testing3"
    generators = "cmake"
    def imports(self):
        #self.copy("hydra_pmi_proxy*", dst="bin", src="bin")
        #self.copy("mpiexec.hydra", dst="bin", src="bin")
        self.copy("mpic++", dst="bin", src="bin")
        self.copy("mpicc", dst="bin", src="bin")
        self.copy("mpicxx", dst="bin", src="bin")
        self.copy("orterun", dst="bin", src="bin")
        self.copy("mpirun", dst="bin", src="bin")
        self.copy("mpiexec", dst="bin", src="bin")
        self.copy("opal_wrapper", dst="bin", src="bin")