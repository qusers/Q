# --- commented out --- not done yet 🤭
# from pathlib import Path
# import subprocess
# from setuptools import setup, find_packages
# from setuptools.command.build_py import build_py

# class QCustomBuild(build_py):
#     # Numpy is used to compile Q, coming already installed with the package.
#     def run(self):
#         # Run the Makefile
#         if not Path('src/q6/makefile').exists():
#             raise FileNotFoundError("makefile not found in 'src/q6/'")
#         subprocess.check_call(['make', '-C', 'src/q6/'])
#         # Run the standard build process
#         super().run()

# setup(
#     name='qligfep',
#     version='0.1',
#     packages=find_packages('src/pypackage'),
#     package_dir={'': 'src/pypackage'},
#     cmdclass={
#         'build_py': QCustomBuild,
#     },
#     # Add additional metadata and parameters as needed
# )
