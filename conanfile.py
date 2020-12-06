from conans import ConanFile, CMake, tools

class UnderwaterRayTracerConan(ConanFile):
    name = "underwater-ray-tracer"
    version = "0.1.0"
    description = "A C++ header-only implementation of an underwater ray tracer"
    author = "Brede-J"
    license = "MIT"
    settings = "os", "arch", "compiler", "build_type"
    generators = "cmake"
    exports_sources = "include/*", "test_package/*", "CMakeLists.txt"
    no_copy_source = True

    def _is_env_true(self, env_key, default="False"):
        return tools.get_env(env_key, default).lower() in ["true", "on"]

    def _is_build_tests_enabled(self):
        return any(self._is_env_true(x) for x in ["BUILD_TESTING", "CONAN_RUN_TESTS"])

    def _is_run_tests_enabled(self):
        return self._is_env_true("CONAN_RUN_TESTS")

    def _is_build_examples_enabled(self):
        return self._is_env_true("BUILD_EXAMPLES")

    def build_requirements(self):
        if self._is_build_tests_enabled():
            self.build_requires("gtest/1.10.0")
            self.build_requires("benchmark/1.5.0")
#            self.build_requires("date/2.4.1")
#            self.build_requires_options['date'].header_only = True
        if self._is_build_examples_enabled():
            self.build_requires("gnuplot-iostream/1.1@qvandijk/stable")
            self.build_requires("boost/1.74.0")     # used by gnuplot-iostream

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        self.copy("*.hpp")

    def package_id(self):
        self.info.header_only()
