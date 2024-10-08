import os
import pathlib
from setuptools import setup
from setuptools import Extension
from setuptools.command.build_ext import build_ext
import distutils.sysconfig
import subprocess

pyinc = distutils.sysconfig.get_python_inc()

class CMakeExtension(Extension):
    def __init__(self, name):
        super().__init__(name, sources=[])


class BuildExt(build_ext):
    def run(self):
        for ext in self.extensions:
            if isinstance(ext, CMakeExtension):
                self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):

        cwd = pathlib.Path().absolute()

        build_temp = f"{pathlib.Path(self.build_temp)}/{ext.name}"
        os.makedirs(build_temp, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.mkdir(parents=True, exist_ok=True)

        config = "Debug" if self.debug else "Release"
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + str(extdir.parent.absolute())+"/pyrtklib5",
            "-DCMAKE_BUILD_TYPE=" + config,
            "-DPYTHON_INCLUDE_DIR="+pyinc
        ]
        if self.debug:
            cmake_args.append("-DDEBUG=ON")

        if os.sys.platform == "darwin":
                cmake_args.append("-DDARWIN=ON")
                print("set gcc compiler successfully")

        if os.name == "nt":
            cmake_args += [
                "-DCMAKE_GENERATOR=Visual Studio 16 2019",  # or your specific Visual Studio version
                "-A", "x64",  # or your specific architecture
                "-DCMAKE_C_FLAGS_RELEASE=/MT",
                "-DCMAKE_CXX_FLAGS_RELEASE=/MT",
                "-DWIN32=ON",
                "-DCMAKE_CXX_FLAGS=/bigobj /DWIN32",  # 添加 /bigobj 和 /DWIN32 选项
                "-DCMAKE_C_FLAGS=/bigobj /DWIN32",  # 添加 /bigobj 和 /DWIN32 选项
                "-DCMAKE_EXE_LINKER_FLAGS=/bigobj",  # 添加 /bigobj 选项
                "-DCMAKE_SHARED_LINKER_FLAGS=/bigobj",  # 添加 /bigobj 选项
                "-DADDITIONAL_LIBRARIES=winmm;ws2_32"
            ]
            print("Windows config successfully")

        build_args = ["--config", config]

        if os.name != "nt":
            build_args += ["--", "-j8"]

        os.chdir(build_temp)
        self.spawn(["cmake", f"{str(cwd)}/{ext.name}"] + cmake_args)
        if not self.dry_run:
            self.spawn(["cmake", "--build", "."] + build_args)
        os.chdir(str(cwd))


pyrtklib5 = CMakeExtension("pyrtklib5")

setup(name="pyrtklib5",
      version="0.2.7",
      description="This is a python binding for rtklib_demo5",
      author="Runzhi Hu",
      author_email = "run-zhi.hu@connect.polyu.hk",
      url = "https://github.com/IPNL-POLYU/pyrtklib_demo5",
      packages=["pyrtklib5"],
      package_data={
        'pyrtklib5':['pyrtklib5.pyi','__init__.py']
      },
      ext_modules=[pyrtklib5],  
      cmdclass={"build_ext": BuildExt},
      long_description=open('readme.md').read(),
      long_description_content_type='text/markdown',
)
