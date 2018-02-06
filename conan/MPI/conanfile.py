from conans import ConanFile, CMake, tools
import os
from conans.tools import download
from conans.tools import unzip
from conans import AutoToolsBuildEnvironment

class OpenmpiConan(ConanFile):
    name = "OpenMPI"
    version = "3.0.0"
    license = "https://www.open-mpi.org/community/license.php"
    url = "https://www.open-mpi.org/software/ompi/v3.0/"
    source_url = "https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.bz2"
    description = "OpenMPI Message Passing library."
    unzipped_path = "openmpi-3.0.0/"
    zip_name = os.path.basename(source_url)
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False]}
    default_options = "shared=False"
    generators = "cmake"

    def source(self):
        if(not os.path.isfile(self.zip_name)):
            download(self.source_url, self.zip_name)
        else:
            print("Found the source file in house")
        unzip(self.zip_name)
        # os.unlink(self.zip_name)

    def build(self):
        self.target = {
            "Macos":"macosx",
            "Linux": "linux",
            "Windows":"mingw"
        }[str(self.settings.os)]

        env = AutoToolsBuildEnvironment(self)
        # env.configure(configure_dir="openmpi-3.0.0/")
        self.run("openmpi-3.0.0/configure CC=clang CXX=clang++")
        env.make()
        # self.run("pwd")
        # self.run("%s %s/%s/configure --enable-mpi-thread-multiple --prefix=%s CC=clang CXX=clang++" % (env.command_line,self.conanfile_directory, self.unzipped_path, self.package_folder))
        # self.run("%s/%s/configure --enable-mpi-cxx --prefix=%s CC=clang CXX=clang++" % (self.conanfile_directory, self.unzipped_path, self.package_folder))
        # self.run("%s make -j 8 install" % env.command_line)


    def package(self):
        self.copy("*.h", dst="include", src="install/include")
        self.copy("*.lib", dst="lib", src="install/lib")
        self.copy("*.a", dst="lib", src="install/lib")

    def package_info(self):
        self.cpp_info.libs = ["mpi", "mpi_cxx"]
#     def source(self):
#         self.run("git clone https://github.com/memsharded/hello.git")
#         self.run("cd hello && git checkout static_shared")
#         # This small hack might be useful to guarantee proper /MT /MD linkage in MSVC
#         # if the packaged project doesn't have variables to set it properly
#         tools.replace_in_file("hello/CMakeLists.txt", "PROJECT(MyHello)", '''PROJECT(MyHello)
# include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
# conan_basic_setup()''')
#
#     def build(self):
#         cmake = CMake(self)
#         cmake.configure(source_folder="hello")
#         cmake.build()
#
#         # Explicit way:
#         # self.run('cmake %s/hello %s' % (self.source_folder, cmake.command_line))
#         # self.run("cmake --build . %s" % cmake.build_config)
#
#     def package(self):
#         self.copy("*.h", dst="include", src="hello")
#         self.copy("*hello.lib", dst="lib", keep_path=False)
#         self.copy("*.dll", dst="bin", keep_path=False)
#         self.copy("*.so", dst="lib", keep_path=False)
#         self.copy("*.dylib", dst="lib", keep_path=False)
#         self.copy("*.a", dst="lib", keep_path=False)
#
#     def package_info(self):
#         self.cpp_info.libs = ["hello"]
