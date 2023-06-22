from cx_Freeze import setup, Executable

setup(name="spectrum",
      version="0",
      description="test",
      executables=[Executable("spectrum.py")])

