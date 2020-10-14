#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
		name         = 'divcleaner',
		version      = '0.1.0',
		description  = 'Divergency cleaning module using a Poisson solver',
		author       = 'Flavio Calvo',
		author_email = 'flavio.calvo@irsol.ch',
		ext_modules  = [Extension("divcleaner", ["divcleaner.c", "bicgstab.c"], libraries=["m"], extra_compile_args=["-ggdb"])],
		include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs()
)
