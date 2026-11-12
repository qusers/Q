@ECHO OFF
pushd %~dp0

if "%SPHINXBUILD%" == "" set SPHINXBUILD=sphinx-build
set SOURCEDIR=source
set BUILDDIR=build

%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
popd
